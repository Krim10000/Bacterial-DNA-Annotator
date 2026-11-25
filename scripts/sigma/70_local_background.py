#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Comparación de ventanas y fondos para la capa 1 (localización) sigma.

La capa 1 del evaluador actual mide si el score del sigma verdadero en el TSS
supera el 95 % de ventanas genómicas aleatorias del mismo largo. Esa definición
mezcla dos cosas: señal del motivo y composición del genoma. Este script barre
ventanas y fondos sin tocar el anotador:

  ventanas:  -60/+20, -80/+20, -100/+20, -120/+20
  fondos:
    genomic  ventanas aleatorias del genoma completo (definición actual)
    shuffle  shuffles mononucleotídicos de la ventana del TSS (100 por TSS)
    flank    subventanas locales de -500/-100 y +100/+500 (100 por TSS)

Para cada combinación (genoma, ventana, fondo) se mide, por sigma y global:

  n, fracción p<=0.05, fracción percentil>=0.95, z mediano y soporte
  compatible con la capa 1 actual (percentil>=0.95 y z>=1.0).

También se reutiliza el ranking z para top-1/top-3 de asignación.

Salidas:
  data/sigma/localization/comparison.json
  data/sigma/localization/comparison.tsv
  data/sigma/localization/assignment.tsv

Uso:
    python scripts/sigma/70_local_background.py --all
    python scripts/sigma/70_local_background.py --all --jobs 8
    python scripts/sigma/70_local_background.py --genome ecoli_k12 --limit 50
    python scripts/sigma/70_local_background.py --selftest
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import json
import statistics
import sys
import zlib
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma.sigma_annotator import (  # noqa: E402
    SIGMA_ORDER, SIGMA_SPACER, _encode, load_pwm_library, load_sequence,
    revcomp, score_sigma,
)

VERSION = "1.0.0"
DEFAULT_WINDOWS = ((60, 20), (80, 20), (100, 20), (120, 20))
DEFAULT_BACKGROUNDS = ("genomic", "shuffle", "flank")
DEFAULT_GENOMES = ("ecoli_k12", "subtilis")
DEFAULT_DECOYS = 200          # fondo genómico global
DEFAULT_SHUFFLES = 100        # shuffles / subventanas de flanco por TSS
DEFAULT_SEED = 17
DEFAULT_PMAX = 0.05
DEFAULT_PCT = 0.95
DEFAULT_Z = 1.0
MIN_FLANK_DECOYS = 20
FLANK_UP = (-500, -100)       # en orientación de transcripción
FLANK_DOWN = (100, 500)

# The evaluator is imported as a module so that gold grouping and sigma labels
# are exactly the same as in 50_eval_layers.py (no duplicated semantics).
def _load_eval_module():
    path = Path(__file__).with_name("50_eval_layers.py")
    spec = importlib.util.spec_from_file_location("sigma_eval_layers_local", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"no se pudo cargar {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_EVAL = _load_eval_module()
load_gold = _EVAL.load_gold


# --------------------------------------------------------------------------
# Secuencia y ventanas
# --------------------------------------------------------------------------
def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name: str | None = None
    buf: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def _fasta_candidates(data_root: Path, genome: str) -> list[Path]:
    if genome == "ecoli_k12":
        return [
            data_root / "allostery" / "work" / "ecoli_k12.fna",
            data_root / "annotations" / "ecoli_k12_genomic.fna",
            data_root / "raw" / "ecoli_k12_genomic.fna",
        ]
    if genome == "subtilis":
        return [
            data_root / "dbtbs" / "genome.fna",
            data_root / "annotations" / "subtilis_genomic.fna",
        ]
    return []


def load_genome(data_root: Path, genome: str,
                gold: list[dict[str, Any]] | None = None
                ) -> tuple[dict[str, str], str]:
    """FASTA explícito primero; `load_sequence` como respaldo."""
    wanted = {r["contig"] for r in (gold or []) if r.get("contig")}
    for path in _fasta_candidates(data_root, genome):
        if not path.exists():
            continue
        seqs = _read_fasta(path)
        if not seqs:
            continue
        if wanted and not (wanted & set(seqs)):
            continue
        return seqs, str(path)
    seqs = load_sequence(data_root, genome)
    if seqs:
        if wanted and not (wanted & set(seqs)):
            # El nombre del contig no coincide: es mejor fallar que evaluar
            # coordenadas contra otra secuencia.
            raise ValueError(f"contigs de {genome} no coinciden con el gold: "
                             f"{sorted(wanted)} vs {sorted(seqs)}")
        return seqs, "load_sequence"
    raise FileNotFoundError(f"sin secuencia para {genome}")


def oriented_window(seq: str, tss: int, strand: str,
                    start: int, end: int) -> str:
    """Intervalo [start,end) en orientación de transcripción, 0 = TSS.

    start/end son offsets: negativos aguas arriba, positivos aguas abajo.
    Para hebra '-' las coordenadas genómicas van al revés y se revcomp.
    """
    if not seq:
        return ""
    p = tss - 1
    if strand == "+":
        a, b = p + start, p + end
    else:
        a, b = p - end + 1, p - start + 1
    a = max(0, a)
    b = min(len(seq), b)
    if b <= a:
        return ""
    sub = seq[a:b]
    return sub if strand == "+" else revcomp(sub)


def tss_window(seq: str, tss: int, strand: str, up: int, down: int) -> str:
    return oriented_window(seq, tss, strand, -up, down + 1)


def _tss_seed(genome: str, tss: int, strand: str, up: int, down: int,
              salt: str) -> int:
    key = f"{genome}:{tss}:{strand}:{up}:{down}:{salt}".encode()
    return int(zlib.crc32(key) & 0xFFFFFFFF)


# --------------------------------------------------------------------------
# Fondos locales
# --------------------------------------------------------------------------
def local_shuffles(window: str, n: int, seed: int) -> list[str]:
    """n shuffles mononucleotídicos de `window`, mismos largo y composición."""
    if not window or n <= 0:
        return []
    arr = np.frombuffer(window.encode("ascii", errors="replace"), dtype=np.uint8)
    rng = np.random.default_rng(seed)
    order = rng.random((n, arr.size)).argsort(axis=1)
    return [arr[row].tobytes().decode("ascii") for row in order]


def flank_segments(seq: str, tss: int, strand: str) -> tuple[str, str]:
    """Flancos locales en orientación de transcripción."""
    return (
        oriented_window(seq, tss, strand, FLANK_UP[0], FLANK_UP[1]),
        oriented_window(seq, tss, strand, FLANK_DOWN[0], FLANK_DOWN[1]),
    )


def sample_flank_windows(seq: str, tss: int, strand: str, length: int,
                         n: int, seed: int) -> tuple[list[str], dict[str, Any]]:
    """Subventanas de largo `length` muestreadas de los dos flancos locales.

    Se reparten ~mitad y mitad; si un flanco no alcanza, se rellena con el otro.
    Cuando no hay material suficiente se devuelve lo disponible y el llamador
    decide si cae a shuffles.
    """
    up, down = flank_segments(seq, tss, strand)
    segs = [s for s in (up, down) if len(s) >= length]
    info: dict[str, Any] = {
        "up_len": len(up), "down_len": len(down),
        "n": 0, "fallback": "",
    }
    if not segs or n <= 0:
        return [], info
    rng = np.random.default_rng(seed)
    if len(segs) == 2:
        quota = [n // 2, n - n // 2]
    else:
        quota = [n]
    out: list[str] = []
    for seg, want in zip(segs, quota):
        max_start = len(seg) - length
        if max_start < 0 or want <= 0:
            continue
        take = min(want, max_start + 1)
        if max_start + 1 >= take:
            starts = rng.choice(max_start + 1, size=take, replace=False)
        else:
            starts = rng.integers(0, max_start + 1, size=take)
        out.extend(seg[int(a):int(a) + length] for a in starts)
    # Rellenar el déficit con el flanco que aún tenga capacidad.
    while len(out) < n:
        added = False
        for seg in segs:
            if len(out) >= n:
                break
            max_start = len(seg) - length
            if max_start < 0:
                continue
            start = int(rng.integers(0, max_start + 1))
            out.append(seg[start:start + length])
            added = True
        if not added:
            break
    info["n"] = len(out)
    return out, info


# --------------------------------------------------------------------------
# Scorer vectorizado para los fondos locales
# --------------------------------------------------------------------------
def _encode_batch(windows: list[str]) -> np.ndarray:
    length = len(windows[0])
    codes = np.empty((len(windows), length), dtype=np.int16)
    for i, window in enumerate(windows):
        if len(window) != length:
            raise ValueError("las ventanas del lote deben tener el mismo largo")
        codes[i] = _encode(window)
    return codes


def _pwm_scores_batch(codes: np.ndarray, pwm: np.ndarray) -> np.ndarray:
    n_seq, length = codes.shape
    width = pwm.shape[0]
    n = length - width + 1
    if n <= 0:
        return np.empty((n_seq, 0), dtype=float)
    out = np.zeros((n_seq, n), dtype=float)
    for k in range(width):
        block = codes[:, k:k + n]
        valid = block >= 0
        out += np.where(valid, pwm[k][np.where(valid, block, 0)], -6.0)
    return out


def _best_geometry_batch(up_scores: np.ndarray, down_scores: np.ndarray,
                         up_w: int, down_w: int,
                         spacers: tuple[int, int]) -> np.ndarray:
    best = np.full(up_scores.shape[0], -np.inf, dtype=float)
    for spacer in range(spacers[0], spacers[1] + 1):
        offset = up_w + spacer
        n = min(up_scores.shape[1] - offset, down_scores.shape[1] - offset)
        if n <= 0:
            continue
        totals = up_scores[:, :n] + down_scores[:, offset:offset + n]
        np.maximum(best, totals.max(axis=1), out=best)
    best[~np.isfinite(best)] = np.nan
    return best


def score_batch(windows: list[str],
                pwms: dict[str, dict[str, Any]]) -> dict[str, np.ndarray]:
    """Score de todos los sigma para un lote de ventanas del mismo largo.

    Es equivalente a `[score_sigma(w, rec, len(w))["score"] for w in windows]`
    pero codifica una sola vez; el selftest lo verifica contra `score_sigma`.
    """
    if not windows:
        return {}
    codes = _encode_batch(windows)
    out: dict[str, np.ndarray] = {}
    for sigma, rec in pwms.items():
        up, down = rec["up"], rec["down"]
        up_scores = _pwm_scores_batch(codes, up["pwm"])
        down_scores = _pwm_scores_batch(codes, down["pwm"])
        if up_scores.shape[1] == 0 or down_scores.shape[1] == 0:
            out[sigma] = np.full(len(windows), np.nan)
            continue
        out[sigma] = _best_geometry_batch(
            up_scores, down_scores, up["width"], down["width"],
            SIGMA_SPACER.get(sigma, (14, 19)))
    return out


# --------------------------------------------------------------------------
# Estadística contra un fondo
# --------------------------------------------------------------------------
def null_stats(score: float, null: Iterable[float] | np.ndarray) -> dict[str, Any]:
    arr = np.asarray(list(null) if not isinstance(null, np.ndarray) else null,
                     dtype=float)
    arr = arr[np.isfinite(arr)]
    n = int(arr.size)
    if n == 0:
        return {"n": 0, "pvalue": 1.0, "pct": 0.0, "z": 0.0,
                "mean": 0.0, "sd": 0.0}
    pvalue = float((1 + int(np.sum(arr >= score))) / (1 + n))
    pct = float(np.mean(arr < score))
    mean = float(arr.mean())
    sd = float(arr.std(ddof=1)) if n > 1 else 0.0
    z = float((score - mean) / (sd if sd > 0 else 1.0))
    return {"n": n, "pvalue": pvalue, "pct": pct, "z": z,
            "mean": mean, "sd": sd}


def _is_supported(stats: dict[str, Any], pmax: float, pct_min: float,
                  zmin: float) -> bool:
    return bool(stats["pct"] >= pct_min and stats["z"] >= zmin)


# --------------------------------------------------------------------------
# Fondos
# --------------------------------------------------------------------------
def genomic_null(seqs: dict[str, str], sigma: str, rec: dict[str, Any],
                 length: int, n: int, seed: int,
                 cache: dict[tuple[str, int], np.ndarray]) -> np.ndarray:
    key = (sigma, length)
    if key in cache:
        return cache[key]
    rng = np.random.default_rng(seed + (zlib.crc32(sigma.encode()) % 1000))
    contigs = [s for s in seqs.values() if len(s) >= length]
    values: list[float] = []
    if contigs:
        for _ in range(n):
            seq = contigs[int(rng.integers(0, len(contigs)))]
            start = int(rng.integers(0, len(seq) - length + 1))
            hit = score_sigma(seq[start:start + length], rec, length)
            if hit is not None:
                values.append(float(hit["score"]))
    arr = np.asarray(values, dtype=float)
    cache[key] = arr
    return arr


# --------------------------------------------------------------------------
# Evaluación de una combinación
# --------------------------------------------------------------------------
def _summarize_stats(rows: list[dict[str, Any]], pmax: float,
                     pct_min: float, zmin: float) -> dict[str, Any]:
    n = len(rows)
    if not n:
        return {"n": 0, "n_p05": 0, "frac_p05": None,
                "n_pct95": 0, "frac_pct95": None,
                "n_supported": 0, "frac_supported": None,
                "median_z": None, "median_pct": None}
    n_p05 = sum(1 for r in rows if r["pvalue"] <= pmax)
    n_pct = sum(1 for r in rows if r["pct"] >= pct_min)
    n_sup = sum(1 for r in rows if _is_supported(r, pmax, pct_min, zmin))
    return {
        "n": n,
        "n_p05": n_p05, "frac_p05": n_p05 / n,
        "n_pct95": n_pct, "frac_pct95": n_pct / n,
        "n_supported": n_sup, "frac_supported": n_sup / n,
        "median_z": statistics.median(r["z"] for r in rows),
        "median_pct": statistics.median(r["pct"] for r in rows),
    }


def _per_sigma_assignment(rows: list[dict[str, Any]]) -> dict[str, Any]:
    labels = sorted({s for r in rows for s in r["true_sigmas"]})
    out: dict[str, Any] = {}
    for sigma in labels:
        tp = fp = fn = 0
        for row in rows:
            true = set(row["true_sigmas"])
            pred = ({row["pred_top1"]} if row.get("n_true", 1) == 1
                    else set(row["pred_top3"]))
            if sigma in true and sigma in pred:
                tp += 1
            elif sigma in pred:
                fp += 1
            elif sigma in true:
                fn += 1
        precision = tp / (tp + fp) if tp + fp else None
        recall = tp / (tp + fn) if tp + fn else None
        f1 = (2 * precision * recall / (precision + recall)
              if precision and recall else None)
        out[sigma] = {"support": sum(1 for r in rows if sigma in r["true_sigmas"]),
                      "tp": tp, "fp": fp, "fn": fn,
                      "precision": precision, "recall": recall, "f1": f1}
    return out


def evaluate_config(data_root: Path, genome: str, up: int, down: int,
                    background: str, *, decoys: int = DEFAULT_DECOYS,
                    shuffles: int = DEFAULT_SHUFFLES,
                    seed: int = DEFAULT_SEED,
                    pmax: float = DEFAULT_PMAX,
                    pct_min: float = DEFAULT_PCT,
                    zmin: float = DEFAULT_Z,
                    limit: int = 0,
                    collect_rows: bool = False) -> dict[str, Any]:
    """Evalúa una combinación (genoma, ventana, fondo)."""
    gold = load_gold(data_root / "sigma" / "gold" / genome / "tss_sigma.tsv",
                     "gold")
    if limit:
        gold = gold[:limit]
    seqs, seq_source = load_genome(data_root, genome, gold)
    pwms = load_pwm_library(data_root)
    if not pwms:
        raise ValueError("sin librería PWM")

    sigma_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    tss_rows: list[dict[str, Any]] = []
    assign_rows: list[dict[str, Any]] = []
    null_sizes: list[int] = []
    fallbacks = 0
    missing_window = 0
    no_scores = 0
    genomic_cache: dict[tuple[str, int], np.ndarray] = {}

    for gold_rec in gold:
        contig = gold_rec["contig"] or next(iter(seqs), "")
        seq = seqs.get(contig, "")
        window = tss_window(seq, gold_rec["tss"], gold_rec["strand"], up, down)
        if not window:
            missing_window += 1
            continue
        length = len(window)
        observed: dict[str, dict[str, Any]] = {}
        for sigma, rec in pwms.items():
            hit = score_sigma(window, rec, length)
            if hit is not None:
                observed[sigma] = hit
        if not observed:
            no_scores += 1
            continue

        # Fondo y estadísticos por sigma.
        stats: dict[str, dict[str, Any]] = {}
        if background == "genomic":
            for sigma, hit in observed.items():
                null = genomic_null(seqs, sigma, pwms[sigma], length,
                                    decoys, seed, genomic_cache)
                stats[sigma] = null_stats(float(hit["score"]), null)
        else:
            salt = "shuffle" if background == "shuffle" else "flank"
            tss_seed = _tss_seed(genome, gold_rec["tss"], gold_rec["strand"],
                                 up, down, salt + f":{seed}")
            null_windows: list[str]
            if background == "shuffle":
                null_windows = local_shuffles(window, shuffles, tss_seed)
            else:
                null_windows, _info = sample_flank_windows(
                    seq, gold_rec["tss"], gold_rec["strand"], length,
                    shuffles, tss_seed)
                if len(null_windows) < MIN_FLANK_DECOYS:
                    null_windows = local_shuffles(window, shuffles, tss_seed)
                    fallbacks += 1
            if not null_windows:
                no_scores += 1
                continue
            null_scores = score_batch(null_windows, pwms)
            null_sizes.append(len(null_windows))
            for sigma, hit in observed.items():
                scores = null_scores.get(sigma)
                if scores is None:
                    stats[sigma] = null_stats(float(hit["score"]), np.array([]))
                else:
                    stats[sigma] = null_stats(float(hit["score"]), scores)

        true_sigmas = [s for s in gold_rec["sigmas"] if s in observed]
        if not true_sigmas:
            continue
        for sigma in true_sigmas:
            sigma_rows[sigma].append(stats[sigma])

        # Métrica global: "algún sigma verdadero soportado" por TSS.
        tss_rows.append({
            "pvalue": min(stats[s]["pvalue"] for s in true_sigmas),
            "pct": max(stats[s]["pct"] for s in true_sigmas),
            "z": max(stats[s]["z"] for s in true_sigmas),
            "supported": any(_is_supported(stats[s], pmax, pct_min, zmin)
                             for s in true_sigmas),
            "p05": any(stats[s]["pvalue"] <= pmax for s in true_sigmas),
            "pct95": any(stats[s]["pct"] >= pct_min for s in true_sigmas),
        })

        # Ranking reutilizado del evaluador: z, -pvalue, score.
        ranking = sorted(
            ((float(stats[s]["z"]), -float(stats[s]["pvalue"]),
              float(observed[s]["score"]), s) for s in observed),
            reverse=True)
        pred = [s for _, _, _, s in ranking]
        assign_rows.append({
            "true_sigmas": true_sigmas,
            "pred_top1": pred[0] if pred else "",
            "pred_top3": pred[:3],
            "n_true": len(true_sigmas),
        })

    per_sigma: dict[str, Any] = {}
    for sigma in SIGMA_ORDER:
        per_sigma[sigma] = _summarize_stats(sigma_rows.get(sigma, []),
                                            pmax, pct_min, zmin)
    n_tss = len(tss_rows)
    global_summary = {
        "n_tss": n_tss,
        "n_p05": sum(1 for r in tss_rows if r["p05"]),
        "frac_p05": (sum(1 for r in tss_rows if r["p05"]) / n_tss
                     if n_tss else None),
        "n_pct95": sum(1 for r in tss_rows if r["pct95"]),
        "frac_pct95": (sum(1 for r in tss_rows if r["pct95"]) / n_tss
                       if n_tss else None),
        "n_supported": sum(1 for r in tss_rows if r["supported"]),
        "frac_supported": (sum(1 for r in tss_rows if r["supported"]) / n_tss
                           if n_tss else None),
        "median_z": (statistics.median(r["z"] for r in tss_rows)
                     if n_tss else None),
        "median_pct": (statistics.median(r["pct"] for r in tss_rows)
                       if n_tss else None),
    }

    n_eval = len(assign_rows)
    top1 = sum(1 for r in assign_rows if r["pred_top1"] in r["true_sigmas"])
    top3_any = sum(1 for r in assign_rows
                   if set(r["pred_top3"]) & set(r["true_sigmas"]))
    top3_all = sum(1 for r in assign_rows
                   if set(r["true_sigmas"]) <= set(r["pred_top3"]))
    assignment = {
        "n_evaluable": n_eval,
        "top1_accuracy": top1 / n_eval if n_eval else None,
        "top3_any": top3_any / n_eval if n_eval else None,
        "top3_all": top3_all / n_eval if n_eval else None,
        "per_sigma": _per_sigma_assignment(assign_rows),
    }

    null_info = {
        "kind": background,
        "n_decoys_genomic": decoys if background == "genomic" else None,
        "n_shuffles_target": shuffles if background != "genomic" else None,
        "n_null_median": (statistics.median(null_sizes) if null_sizes
                          else (decoys if background == "genomic" else None)),
        "n_null_min": min(null_sizes) if null_sizes else None,
        "n_tss_low_null": sum(1 for x in null_sizes if x < MIN_FLANK_DECOYS),
        "n_flank_fallback_shuffle": fallbacks if background == "flank" else None,
    }

    result = {
        "genome": genome,
        "window": [-up, down],
        "window_up": up,
        "window_down": down,
        "background": background,
        "sequence_source": seq_source,
        "n_gold_tss": len(gold),
        "n_tss_evaluated": len(tss_rows),
        "n_missing_window": missing_window,
        "n_no_scores": no_scores,
        "null_info": null_info,
        "global_summary": global_summary,
        "per_sigma": per_sigma,
        "assignment": assignment,
    }
    if collect_rows:
        result["rows"] = assign_rows
    return result


# --------------------------------------------------------------------------
# Ejecución y salidas
# --------------------------------------------------------------------------
def _run_task(args: tuple[str, str, int, int, str, dict[str, Any]]) -> dict[str, Any]:
    genome, data_root_s, up, down, background, params = args
    return evaluate_config(Path(data_root_s), genome, up, down, background,
                           **params)


def _fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_outputs(results: list[dict[str, Any]], data_root: Path, out_dir: Path,
                  config: dict[str, Any], baseline: dict[str, Any],
                  recommendation: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "config": config,
        "gold": {
            g: {
                "path": str(data_root / "sigma" / "gold" / g / "tss_sigma.tsv"),
                "n_tss": next((r["n_gold_tss"] for r in results
                               if r["genome"] == g), None),
            }
            for g in sorted({r["genome"] for r in results})
        },
        "baseline_reproduction": baseline,
        "recommendation": recommendation,
        "results": results,
        "limits": [
            "La capa 1 mide señal de motivo contra un fondo; no prueba el TSS "
            "por sí mismo.",
            "score_sigma busca la mejor geometría en toda la ventana, sin "
            "restringir -35/-10 a una distancia biológica del TSS.",
            "Los fondos locales (shuffle/flank) son condicionales a cada TSS; "
            "sus porcentajes no son directamente intercambiables con el fondo "
            "genómico global.",
            "El gold de K-12 es RegulonDB y el de S. subtilis es DBTBS; no son "
            "colecciones exhaustivas ni homogéneas.",
        ],
    }
    (out_dir / "comparison.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    cols = ["genome", "window_up", "window_down", "background", "scope", "n",
            "n_p05", "frac_p05", "n_pct95", "frac_pct95", "median_z",
            "median_pct", "n_supported", "frac_supported"]
    lines = ["\t".join(cols)]
    for r in results:
        for scope in ["global"] + list(SIGMA_ORDER):
            m = (r["global_summary"] if scope == "global"
                 else r["per_sigma"][scope])
            lines.append("\t".join(_fmt(x) for x in [
                r["genome"], r["window_up"], r["window_down"], r["background"],
                scope, m.get("n_tss", m.get("n")), m.get("n_p05"),
                m.get("frac_p05"), m.get("n_pct95"), m.get("frac_pct95"),
                m.get("median_z"), m.get("median_pct"), m.get("n_supported"),
                m.get("frac_supported")]))
    (out_dir / "comparison.tsv").write_text("\n".join(lines) + "\n",
                                             encoding="utf-8")

    acols = ["genome", "window_up", "window_down", "background", "scope", "n",
             "top1_accuracy", "top3_any", "top3_all", "precision", "recall",
             "f1"]
    alines = ["\t".join(acols)]
    for r in results:
        a = r["assignment"]
        alines.append("\t".join(_fmt(x) for x in [
            r["genome"], r["window_up"], r["window_down"], r["background"],
            "global", a["n_evaluable"], a["top1_accuracy"], a["top3_any"],
            a["top3_all"], "", "", ""]))
        for sigma in SIGMA_ORDER:
            m = a["per_sigma"].get(sigma)
            if not m:
                continue
            alines.append("\t".join(_fmt(x) for x in [
                r["genome"], r["window_up"], r["window_down"], r["background"],
                sigma, m["support"], "", "", "", m["precision"], m["recall"],
                m["f1"]]))
    (out_dir / "assignment.tsv").write_text("\n".join(alines) + "\n",
                                             encoding="utf-8")


def _baseline_check(results: list[dict[str, Any]], data_root: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"window": [-80, 20], "background": "genomic",
                           "per_genome": {}}
    for r in results:
        if r["window"] != [-80, 20] or r["background"] != "genomic":
            continue
        path = data_root / "sigma" / "layers" / r["genome"] / "layer_metrics.json"
        existing = None
        if path.exists():
            try:
                existing = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                existing = None
        old = (existing or {}).get("layer1_localization", {}).get("frac_supported")
        new = r["global_summary"]["frac_supported"]
        out["per_genome"][r["genome"]] = {
            "n_tss": r["global_summary"]["n_tss"],
            "frac_supported_new": new,
            "frac_supported_existing_layers": old,
            "delta_vs_existing": (new - old) if (new is not None and old is not None)
            else None,
        }
    return out


def _recommend(results: list[dict[str, Any]], backgrounds: list[str],
               guardrail: float = 0.01) -> dict[str, Any]:
    """Regla explícita con guardarraíl de asignación.

    La capa 1 es el objetivo, pero no se acepta una combinación que empeore
    top-1/top-3 de asignación más de `guardrail` (1 punto porcentual) respecto
    de la referencia actual (-80/+20, fondo genómico). Entre las aceptadas se
    maximiza la media macro de soporte TSS; desempate por p<=0.05 y ventana
    más corta.
    """
    by_key: dict[tuple[int, int, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for r in results:
        by_key[(r["window_up"], r["window_down"], r["background"])][r["genome"]] = r
    genomes = sorted({r["genome"] for r in results})

    baseline = by_key.get((80, 20, "genomic"), {})
    baseline_top1 = baseline_top3 = None
    if len(baseline) == len(genomes):
        top1 = [baseline[g]["assignment"]["top1_accuracy"] for g in genomes]
        top3 = [baseline[g]["assignment"]["top3_any"] for g in genomes]
        if all(x is not None for x in top1 + top3):
            baseline_top1 = statistics.mean(top1)
            baseline_top3 = statistics.mean(top3)

    candidates = []
    for (up, down, bg), per_genome in by_key.items():
        if bg not in backgrounds or len(per_genome) != len(genomes):
            continue
        supports = [per_genome[g]["global_summary"]["frac_supported"]
                    for g in genomes]
        p05 = [per_genome[g]["global_summary"]["frac_p05"] for g in genomes]
        top1 = [per_genome[g]["assignment"]["top1_accuracy"] for g in genomes]
        top3 = [per_genome[g]["assignment"]["top3_any"] for g in genomes]
        if any(x is None for x in supports):
            continue
        macro_top1 = (statistics.mean(top1)
                      if all(x is not None for x in top1) else None)
        macro_top3 = (statistics.mean(top3)
                      if all(x is not None for x in top3) else None)
        allowed = True
        if baseline_top1 is not None and macro_top1 is not None:
            allowed = (macro_top1 >= baseline_top1 - guardrail and
                       macro_top3 is not None and
                       macro_top3 >= baseline_top3 - guardrail)
        candidates.append({
            "window": [-up, down], "window_up": up, "window_down": down,
            "background": bg,
            "macro_frac_supported": statistics.mean(supports),
            "macro_frac_p05": statistics.mean(p05),
            "macro_top1": macro_top1,
            "macro_top3": macro_top3,
            "allowed_by_guardrail": allowed,
            "per_genome": {
                g: {
                    "frac_supported": per_genome[g]["global_summary"]["frac_supported"],
                    "frac_p05": per_genome[g]["global_summary"]["frac_p05"],
                    "top1_accuracy": per_genome[g]["assignment"]["top1_accuracy"],
                    "top3_any": per_genome[g]["assignment"]["top3_any"],
                } for g in genomes},
        })
    if not candidates:
        return {"rule": "sin combinaciones completas", "chosen": None}

    accepted = [c for c in candidates if c["allowed_by_guardrail"]]
    pool = accepted or candidates
    pool.sort(key=lambda c: (-c["macro_frac_supported"],
                             -c["macro_frac_p05"], c["window_up"]))
    chosen = pool[0]
    runner = pool[1] if len(pool) > 1 else None
    return {
        "rule": ("maximizar la media macro entre genomas de frac_supported "
                 "(TSS con >=1 sigma verdadero con percentil>=0.95 y z>=1.0), "
                 "sin empeorar top-1 ni top-3 más de 1 punto porcentual "
                 "respecto de -80/+20 genómico; desempate por frac_p05 y "
                 "ventana más corta"),
        "guardrail": {
            "reference": "window [-80,20], background genomic",
            "baseline_macro_top1": baseline_top1,
            "baseline_macro_top3": baseline_top3,
            "max_drop": guardrail,
        },
        "n_candidates": len(candidates),
        "n_accepted": len(accepted),
        "chosen": chosen,
        "runner_up": runner,
    }


def selftest() -> str:
    seq = "ACGTACGTACGTACGTACGT"
    assert tss_window(seq, 10, "+", 4, 2) == seq[5:12]
    assert tss_window(seq, 10, "-", 4, 2) == revcomp(seq[7:14])
    # Un shuffle conserva composición y largo, y es reproducible.
    a = local_shuffles("AACCGGTT", 5, 1)
    b = local_shuffles("AACCGGTT", 5, 1)
    assert a == b and all(sorted(x) == sorted("AACCGGTT") for x in a)
    # Los flancos no se solapan con la ventana central -20/+20.
    long_seq = "ACGT" * 1000
    flanks = flank_segments(long_seq, 1000, "+")
    assert len(flanks[0]) == 400 and len(flanks[1]) == 400
    central = tss_window(long_seq, 1000, "+", 20, 20)
    assert flanks[0] == long_seq[499:899]
    assert flanks[1] == long_seq[1099:1499]
    assert central == long_seq[979:1020]
    # null_stats reproduce la definición del evaluador.
    st = null_stats(5.0, np.array([1.0, 2.0, 3.0, 4.0]))
    assert st["n"] == 4 and st["pvalue"] == 0.2 and st["pct"] == 1.0
    assert st["z"] > 0
    # El scorer vectorizado coincide exactamente con score_sigma.
    pwms = load_pwm_library(ROOT / "data")
    if pwms:
        window = seq * 8
        windows = [window[i:i + 60] for i in range(3)] + [revcomp(window)[:60]]
        batched = score_batch(windows, pwms)
        for sig, rec in pwms.items():
            expected = [score_sigma(w, rec, len(w))["score"] for w in windows]
            assert np.allclose(batched[sig], expected, atol=1e-9), sig
    return "ventanas, shuffles, flancos, null_stats y scorer batch OK"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", action="append", choices=list(DEFAULT_GENOMES),
                    help="Genoma a evaluar (repetible). Por defecto ambos.")
    ap.add_argument("--all", action="store_true",
                    help="Evalúa todos los genomas por defecto.")
    ap.add_argument("--windows",
                    default=",".join(str(up) for up, _ in DEFAULT_WINDOWS),
                    help="Ventanas aguas arriba separadas por coma.")
    ap.add_argument("--backgrounds", default="genomic,shuffle,flank",
                    help="Fondos separados por coma.")
    ap.add_argument("--down", type=int, default=20)
    ap.add_argument("--decoys", type=int, default=DEFAULT_DECOYS)
    ap.add_argument("--shuffles", type=int, default=DEFAULT_SHUFFLES)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--pmax", type=float, default=DEFAULT_PMAX)
    ap.add_argument("--pct", type=float, default=DEFAULT_PCT)
    ap.add_argument("--zmin", type=float, default=DEFAULT_Z)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0,
                    help="Limita TSS por genoma (smoke test).")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "sigma" / "localization")
    ap.add_argument("--from-json", type=Path, default=None,
                    help="Recomputa recomendación y salidas desde un "
                         "comparison.json ya generado (no re-evalúa TSS).")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0

    data_root = ROOT / "data"
    if args.from_json:
        payload = json.loads(args.from_json.read_text(encoding="utf-8"))
        results = payload.get("results") or []
        if not results:
            ap.error(f"{args.from_json} no contiene resultados")
        config = payload.get("config") or {}
        backgrounds = config.get("backgrounds") or list(DEFAULT_BACKGROUNDS)
        baseline = _baseline_check(results, data_root)
        recommendation = _recommend(results, backgrounds)
        write_outputs(results, data_root, args.out, config, baseline,
                      recommendation)
        print(f"escrito {args.out}/comparison.json (desde {args.from_json})")
        chosen = recommendation.get("chosen") or {}
        if chosen:
            print(f"recomendación: ventana {chosen['window']} fondo "
                  f"{chosen['background']} · macro soporte "
                  f"{chosen['macro_frac_supported']:.1%}")
        return 0

    genomes = args.genome or list(DEFAULT_GENOMES)
    windows = []
    for token in str(args.windows).split(","):
        token = token.strip()
        if not token:
            continue
        up = int(token)
        if up <= 0:
            ap.error("las ventanas aguas arriba deben ser positivas")
        windows.append((up, args.down))
    backgrounds = [b.strip() for b in str(args.backgrounds).split(",")
                   if b.strip()]
    unknown = set(backgrounds) - set(DEFAULT_BACKGROUNDS)
    if unknown:
        ap.error(f"fondos desconocidos: {sorted(unknown)}")
    params = {"decoys": args.decoys, "shuffles": args.shuffles,
              "seed": args.seed, "pmax": args.pmax, "pct_min": args.pct,
              "zmin": args.zmin, "limit": args.limit}
    tasks = [(g, str(data_root), up, down, bg, params)
             for g in genomes for up, down in windows for bg in backgrounds]
    results: list[dict[str, Any]] = []
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(_run_task, t): t for t in tasks}
            for fut in as_completed(futures):
                results.append(fut.result())
    else:
        for task in tasks:
            results.append(_run_task(task))
    order_bg = {b: i for i, b in enumerate(DEFAULT_BACKGROUNDS)}
    order_g = {g: i for i, g in enumerate(genomes)}
    results.sort(key=lambda r: (order_g.get(r["genome"], 99), -r["window_up"],
                                order_bg.get(r["background"], 99)))
    baseline = _baseline_check(results, data_root)
    recommendation = _recommend(results, backgrounds)
    config = {
        "windows": [[-u, d] for u, d in windows],
        "backgrounds": backgrounds,
        "decoys_genomic": args.decoys,
        "shuffles_local": args.shuffles,
        "seed": args.seed,
        "pmax": args.pmax,
        "pct_min": args.pct,
        "zmin": args.zmin,
        "limit": args.limit,
        "metric_supported": "pct>=pct_min and z>=zmin (definición capa 1 actual)",
    }
    write_outputs(results, data_root, args.out, config, baseline, recommendation)
    print(f"escrito {args.out}/comparison.json, comparison.tsv, assignment.tsv")
    chosen = recommendation.get("chosen") or {}
    if chosen:
        print(f"recomendación: ventana {chosen['window']} fondo "
              f"{chosen['background']} · macro soporte "
              f"{chosen['macro_frac_supported']:.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
