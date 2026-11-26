#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evaluación de la capa de contexto sigma con validación fuera de muestra.

Protocolo:
  * gold curado `data/sigma/gold/<genoma>/tss_sigma.tsv` (no se modifica);
  * PWM de la misma especie (K-12: RegulonDB; B. subtilis: DBTBS propio);
  * validación cruzada bloqueada por coordenada (5 bloques por defecto);
  * las tablas de log-odds y la temperatura de calibración se ajustan SOLO
    con las particiones de entrenamiento;
  * la escala global de contexto se elige en entrenamiento (0 = apagado);
  * el top-3 se calcula con el ranking completo, no con el modo acotado.

Salidas por genoma:
  data/sigma/context/<genoma>/
    context_metrics.json          métricas, ablación, calibración, sigma54
    context_sites.tsv             una fila por TSS
    context_predictions.tsv       una fila por TSS x sigma con LLR por evidencia
    context_evidence.tsv          tablas de log-odds (ajuste completo, descriptivo)
    context_evidence_coverage.tsv cobertura, fuente y peso declarado
    context_report.md             informe legible

Uso:
    python scripts/sigma/51_eval_context.py --genome ecoli_k12
    python scripts/sigma/51_eval_context.py --genome subtilis
    python scripts/sigma/51_eval_context.py --all
    python scripts/sigma/51_eval_context.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma import sigma_context as C  # noqa: E402
from scripts.sigma.sigma_annotator import (  # noqa: E402
    SIGMA_ORDER, load_pwm_library, load_sequence,
)

VERSION = "1.0.0"
DEFAULT_UP = 80
DEFAULT_DOWN = 20
DEFAULT_DECOYS = 200
DEFAULT_SEED = 17
DEFAULT_FOLDS = 5
DEFAULT_GAP = 1.0
DEFAULT_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0)
DEFAULT_MIN_COVERAGE = 0.05
DEFAULT_MIN_SUPPORT = 5


def _load_layers():
    spec = importlib.util.spec_from_file_location(
        "sigma_eval_layers_ctx", ROOT / "scripts" / "sigma" / "50_eval_layers.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


LAYERS = _load_layers()


# --------------------------------------------------------------------------
# Secuencia: fallback local sin Biopython
# --------------------------------------------------------------------------
def _read_gbff_origin(path: Path) -> dict[str, str]:
    """Extrae secuencias de un GBFF sin depender de Biopython."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records: dict[str, str] = {}
    name: str | None = None
    version: str | None = None
    buf: list[str] = []
    in_origin = False
    with opener(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                if name and buf:
                    records[version or name] = "".join(buf).upper()
                name = line.split()[1]
                version = None
                buf = []
                in_origin = False
            elif line.startswith("VERSION"):
                parts = line.split()
                if len(parts) > 1:
                    version = parts[1]
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif line.startswith("//"):
                if name and buf:
                    records[version or name] = "".join(buf).upper()
                name, version, buf, in_origin = None, None, [], False
            elif in_origin:
                buf.append("".join(ch for ch in line if ch.isalpha()))
        if name and buf:
            records[version or name] = "".join(buf).upper()
    return records


def load_sequences(data_root: Path, genome: str) -> tuple[dict[str, str], str]:
    """Carga secuencia con varias fuentes locales; declara la elegida."""
    seqs = load_sequence(data_root, genome)
    if seqs:
        return seqs, "sigma_annotator.load_sequence"
    seqs = LAYERS._fallback_sequence(data_root, genome)
    if seqs:
        return seqs, "50_eval_layers._fallback_sequence"
    candidates = [
        data_root / "allostery" / "work" / f"{genome}.fna",
        data_root / "terminators" / f"{genome}.fna",
        data_root / "dbtbs" / "genome.fna" if genome == "subtilis" else None,
        data_root / "raw" / f"{genome}_genomic.fna",
        data_root / "raw" / f"{genome}_genomic.fna.gz",
    ]
    for path in candidates:
        if path and path.exists():
            got = LAYERS._read_fasta(path)
            if got:
                return got, str(path)
    for path in (
        data_root / "annotations" / f"{genome}_genomic.gbff",
        data_root / "annotations" / f"{genome}_genomic.gbff.gz",
        data_root / "test_ecoli.gbff" if genome == "ecoli_k12" else None,
    ):
        if path and path.exists():
            got = _read_gbff_origin(path)
            if got:
                return got, str(path)
    return {}, ""


# --------------------------------------------------------------------------
# PWM por especie
# --------------------------------------------------------------------------
def resolve_pwm(data_root: Path, genome: str, pwm: Path | None) -> Path:
    if pwm is not None:
        return Path(pwm)
    if genome == "subtilis":
        return data_root / "dbtbs" / "promoter_pwms_subtilis_all.tsv"
    return data_root / "regulondb" / "gt" / "promoter_pwms.tsv"


def is_ecoli_transfer(genome: str, pwm: Path) -> bool:
    return genome == "subtilis" and "regulondb" in str(pwm).lower()


# --------------------------------------------------------------------------
# Puntuación del motivo
# --------------------------------------------------------------------------
def score_gold(rows: Sequence[dict[str, Any]], seqs: dict[str, str],
               pwms: dict[str, dict[str, Any]], up: int, down: int,
               decoys: int, seed: int) -> dict[tuple[str, int, str], dict[str, Any]]:
    cache: dict[tuple[str, int], dict[str, Any]] = {}
    out: dict[tuple[str, int, str], dict[str, Any]] = {}
    for gold in rows:
        key = (gold["contig"], int(gold["tss"]), gold["strand"])
        out[key] = LAYERS._score_tss(seqs, gold, pwms, up, down,
                                     decoys, seed, cache)
    return out


def blocked_folds(keys: Sequence[tuple[str, int, str]], k: int) -> dict[tuple[str, int, str], int]:
    """Bloques contiguos en orden genómico para no filtrar TSS vecinos."""
    order = sorted(keys, key=lambda x: (x[0], x[1], x[2]))
    fold: dict[tuple[str, int, str], int] = {}
    n = len(order)
    for i, key in enumerate(order):
        fold[key] = min(k - 1, i * k // n) if n else 0
    return fold


# --------------------------------------------------------------------------
# Métricas
# --------------------------------------------------------------------------
def _evaluable(site: dict[str, Any]) -> bool:
    return bool(site["true"]) and bool(site["z"])


def _classification(site: dict[str, Any], top1: str, top3: Sequence[str]) -> str:
    true = set(site["true"])
    if not true:
        return "unlabeled"
    if not site["z"]:
        return "no_pwm"
    if not top1:
        return "no_prediction"
    if top1 in true:
        return "correct" if len(true) == 1 else "hit_one_of_multi"
    if set(top3) & true:
        return "partial_top3"
    return "missed"


def _per_sigma(sites: Sequence[dict[str, Any]],
               pred: dict[tuple[str, int, str], dict[str, Any]]) -> dict[str, Any]:
    labels = sorted({s for site in sites for s in site["true"]})
    out: dict[str, Any] = {}
    for sigma in labels:
        tp1 = fp1 = fn1 = tp3 = fp3 = fn3 = 0
        for site in sites:
            if not _evaluable(site):
                continue
            key = site["key"]
            true = set(site["true"])
            top1 = pred[key]["top1"]
            top3 = pred[key]["top3"]
            if sigma in true and sigma == top1:
                tp1 += 1
            elif sigma == top1:
                fp1 += 1
            elif sigma in true:
                fn1 += 1
            if sigma in true and sigma in top3:
                tp3 += 1
            elif sigma in top3:
                fp3 += 1
            elif sigma in true:
                fn3 += 1

        def _f1(tp: int, fp: int, fn: int) -> tuple[float | None, float | None, float | None]:
            p = tp / (tp + fp) if tp + fp else None
            r = tp / (tp + fn) if tp + fn else None
            f = (2 * p * r / (p + r)) if p and r else None
            return p, r, f

        p1, r1, f1 = _f1(tp1, fp1, fn1)
        p3, r3, f3 = _f1(tp3, fp3, fn3)
        out[sigma] = {
            "support": sum(1 for s in sites if sigma in s["true"]),
            "top1": {"tp": tp1, "fp": fp1, "fn": fn1,
                     "precision": p1, "recall": r1, "f1": f1},
            "top3": {"tp": tp3, "fp": fp3, "fn": fn3,
                     "precision": p3, "recall": r3, "f1": f3},
        }
    return out


def _multilabel(sites: Sequence[dict[str, Any]],
                pred: dict[tuple[str, int, str], dict[str, Any]]) -> dict[str, Any]:
    ev = [s for s in sites if _evaluable(s)]
    multi = [s for s in ev if len(s["true"]) > 1]
    n = len(ev)
    top1_any = sum(1 for s in ev if pred[s["key"]]["top1"] in s["true"])
    top3_any = sum(1 for s in ev if set(pred[s["key"]]["top3"]) & s["true"])
    top3_all = sum(1 for s in ev if set(s["true"]) <= set(pred[s["key"]]["top3"]))
    exact_top1 = sum(1 for s in multi if set(s["true"]) == {pred[s["key"]]["top1"]})
    jaccard = [
        len(set(pred[s["key"]]["top3"]) & s["true"]) / len(set(pred[s["key"]]["top3"]) | s["true"])
        for s in ev if pred[s["key"]]["top3"] or s["true"]
    ]
    return {
        "n": n, "n_single": n - len(multi), "n_multi": len(multi),
        "top1_any": top1_any / n if n else None,
        "top3_any": top3_any / n if n else None,
        "top3_all": top3_all / n if n else None,
        "multi_top1_any": (sum(1 for s in multi if pred[s["key"]]["top1"] in s["true"]) / len(multi)
                           if multi else None),
        "multi_exact_top1": exact_top1 / len(multi) if multi else None,
        "jaccard_top3": round(float(np.mean(jaccard)), 6) if jaccard else None,
    }


def _majority(sites: Sequence[dict[str, Any]]) -> tuple[str, float | None]:
    ev = [s for s in sites if _evaluable(s)]
    counts = Counter(sigma for site in ev for sigma in site["true"])
    if not ev or not counts:
        return "", None
    sigma, n = counts.most_common(1)[0]
    return sigma, n / len(ev)


def _classes(sites: Sequence[dict[str, Any]],
             pred: dict[tuple[str, int, str], dict[str, Any]]) -> dict[str, int]:
    return dict(Counter(
        _classification(site, pred[site["key"]]["top1"], pred[site["key"]]["top3"])
        for site in sites if _evaluable(site)
    ).most_common())


def summarize(sites: Sequence[dict[str, Any]],
              pred: dict[tuple[str, int, str], dict[str, Any]]) -> dict[str, Any]:
    ml = _multilabel(sites, pred)
    majority_sigma, majority = _majority(sites)
    return {
        "top1": ml["top1_any"],
        "top3_any": ml["top3_any"],
        "top3_all": ml["top3_all"],
        "n_evaluable": ml["n"],
        "n_single": ml["n_single"],
        "n_multi": ml["n_multi"],
        "multi_top1_any": ml["multi_top1_any"],
        "multi_exact_top1": ml["multi_exact_top1"],
        "jaccard_top3": ml["jaccard_top3"],
        "majority_sigma": majority_sigma,
        "majority_baseline": majority,
        "per_sigma": _per_sigma(sites, pred),
        "classes": _classes(sites, pred),
    }


# --------------------------------------------------------------------------
# Predicción y calibración
# --------------------------------------------------------------------------
def precompute_llr(sites: Sequence[dict[str, Any]], ctx: C.ContextTable,
                   models: dict[int, C.LogOddsModel],
                   fold: dict[tuple[str, int, str], int],
                   evidence_names: Sequence[str]
                   ) -> dict[tuple[str, int, str], dict[str, dict[str, float]]]:
    """Cachea el LLR crudo (peso 1.0) por TSS x sigma x evidencia.

    Las ablaciones solo cambian los pesos; el LLR no depende del peso.
    """
    cache: dict[tuple[str, int, str], dict[str, dict[str, float]]] = {}
    for site in sites:
        key = site["key"]
        f = fold[key]
        row = ctx.row(key)
        cache[key] = {
            sigma: {name: models[f].llr(row, sigma, {name: 1.0})
                    for name in evidence_names}
            for sigma in site["z"]
        }
    return cache


def predict(sites: Sequence[dict[str, Any]], ctx: C.ContextTable,
            llr_cache: dict[tuple[str, int, str], dict[str, dict[str, float]]],
            fold: dict[tuple[str, int, str], int],
            weights: dict[str, float], scale_by_fold: dict[int, float],
            gap: float | None) -> dict[tuple[str, int, str], dict[str, Any]]:
    out: dict[tuple[str, int, str], dict[str, Any]] = {}
    for site in sites:
        key = site["key"]
        f = fold[key]
        z = site["z"]
        llr = llr_cache[key]
        scale = float(scale_by_fold[f])
        final = {
            sigma: float(zv) + scale * sum(
                float(weights.get(name, 0.0)) * float(llr.get(sigma, {}).get(name, 0.0))
                for name in weights if weights.get(name, 0.0))
            for sigma, zv in z.items()
        }
        motif_top1 = max(z, key=z.get) if z else ""
        motif_top3 = C.rank_topk(z, 3)
        out[key] = {
            "final": final,
            "top1": C.select_top1(final, z, gap),
            "top3": C.rank_topk(final, 3),
            "motif_top1": motif_top1,
            "motif_top3": motif_top3,
            "scale": scale,
            "llr": llr,
        }
    return out


def _temperature_records(sites: Sequence[dict[str, Any]],
                         ctx: C.ContextTable,
                         llr_cache: dict[tuple[str, int, str], dict[str, dict[str, float]]],
                         fold: dict[tuple[str, int, str], int],
                         weights: dict[str, float],
                         scale_by_fold: dict[int, float],
                         which: str) -> dict[int, list[dict[str, Any]]]:
    records: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for site in sites:
        if not _evaluable(site):
            continue
        key = site["key"]
        f = fold[key]
        if which == "motif":
            scores = dict(site["z"])
        else:
            scale = float(scale_by_fold[f])
            llr = llr_cache[key]
            scores = {
                sigma: zv + scale * sum(
                    float(weights.get(name, 0.0)) * float(llr.get(sigma, {}).get(name, 0.0))
                    for name in weights if weights.get(name, 0.0))
                for sigma, zv in site["z"].items()
            }
        records[f].append({"scores": scores, "true": set(site["true"])})
    return records


def calibrated(sites: Sequence[dict[str, Any]], pred: dict[tuple[str, int, str], dict[str, Any]],
               train_records: dict[int, list[dict[str, Any]]],
               fold: dict[tuple[str, int, str], int],
               sigmas: Sequence[str]) -> dict[str, Any]:
    """Temperatura ajustada por fold en entrenamiento; métricas en test."""
    temperatures = {f: C.fit_temperature(recs) for f, recs in train_records.items()}
    test_records: list[dict[str, Any]] = []
    for site in sites:
        if not _evaluable(site):
            continue
        key = site["key"]
        f = fold[key]
        probs = C.softmax(pred[key]["final"], temperatures[f], sigmas)
        test_records.append({"probs": probs, "true": set(site["true"])})
    metrics = C.calibration_metrics(test_records, sigmas=sigmas)
    metrics["temperatures"] = {str(f): round(float(t), 3) for f, t in sorted(temperatures.items())}
    return metrics


def select_scales(sites: Sequence[dict[str, Any]], ctx: C.ContextTable,
                  llr_cache: dict[tuple[str, int, str], dict[str, dict[str, float]]],
                  fold: dict[tuple[str, int, str], int], folds: int,
                  weights: dict[str, float], scale_grid: Sequence[float],
                  gap: float | None, min_gain: float = 0.01) -> dict[int, float]:
    """Elige una escala por fold maximizando top-1 en ENTRENAMIENTO.

    Se queda con 0.0 salvo que otra escala mejore > `min_gain`. Así la capa
    se apaga sola donde el contexto no aporta y no se elige sobre test.
    """
    selected: dict[int, float] = {}
    for f in range(folds):
        train = [s for s in sites if fold[s["key"]] != f and _evaluable(s)]
        if not train:
            selected[f] = 0.0
            continue
        values: dict[float, float] = {}
        for scale in scale_grid:
            hits = 0
            for site in train:
                key = site["key"]
                z = site["z"]
                llr = llr_cache[key]
                final = {
                    sigma: zv + scale * sum(
                        float(weights.get(name, 0.0)) * float(llr.get(sigma, {}).get(name, 0.0))
                        for name in weights if weights.get(name, 0.0))
                    for sigma, zv in z.items()
                }
                hits += int(C.select_top1(final, z, gap) in site["true"])
            values[float(scale)] = hits / len(train)
        best_scale, best_value = 0.0, values.get(0.0, 0.0)
        for scale in scale_grid:
            scale = float(scale)
            if scale == 0.0:
                continue
            if values.get(scale, 0.0) > best_value + min_gain:
                best_scale, best_value = scale, values[scale]
        selected[f] = best_scale
    return selected


def fixed_scales(value: float, folds: int) -> dict[int, float]:
    return {f: float(value) for f in range(folds)}


# --------------------------------------------------------------------------
# Calibración sin cambiar top-1 (modo confidence-only)
# --------------------------------------------------------------------------
def _binary_nll(conf: np.ndarray, y: np.ndarray) -> float:
    conf = np.clip(conf, 1e-6, 1 - 1e-6)
    return float(-(y * np.log(conf) + (1 - y) * np.log(1 - conf)).mean())


def _fixed_label_calibration(conf: np.ndarray, y: np.ndarray, n_bins: int = 10) -> dict[str, Any]:
    if conf.size == 0:
        return {"n": 0, "ece": None, "brier": None, "nll": None,
                "mean_confidence": None, "accuracy": None}
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        mask = (conf >= lo) & ((conf < hi) if i < n_bins - 1 else (conf <= hi))
        if mask.any():
            ece += float(mask.sum()) / len(conf) * abs(float(conf[mask].mean())
                                                       - float(y[mask].mean()))
    return {
        "n": int(conf.size),
        "ece": round(float(ece), 6),
        "brier": round(float(((conf - y) ** 2).mean()), 6),
        "nll": round(_binary_nll(conf, y), 6),
        "mean_confidence": round(float(conf.mean()), 6),
        "accuracy": round(float(y.mean()), 6),
    }


def confidence_only(sites: Sequence[dict[str, Any]], ctx: C.ContextTable,
                    llr_cache: dict[tuple[str, int, str], dict[str, dict[str, float]]],
                    fold: dict[tuple[str, int, str], int], folds: int,
                    weights: dict[str, float], scales: Sequence[float],
                    sigmas: Sequence[str],
                    grid: Sequence[float] | None = None) -> dict[str, Any]:
    """Recalibra la confianza del top-1 del motivo SIN cambiar top-1/top-3.

    Por fold elige (escala, temperatura) minimizando NLL binaria del acierto
    del top-1 en entrenamiento. No usa etiquetas de test ni cambia la
    predicción: mide si el contexto mejora la confianza de lo que el motivo
    ya eligió.
    """
    grid = list(grid) if grid is not None else [i / 10 for i in range(2, 101)]
    ev = [s for s in sites if _evaluable(s)]
    result: dict[str, Any] = {"baseline": {}, "context": {},
                              "selected": {}, "temperatures": {}}
    base_conf: list[float] = []
    base_y: list[float] = []
    ctx_conf: list[float] = []
    ctx_y: list[float] = []
    for f in range(folds):
        train = [s for s in ev if fold[s["key"]] != f]
        test = [s for s in ev if fold[s["key"]] == f]
        if not train or not test:
            continue
        train_keys = [s["key"] for s in train]
        top1_train = [predict_top1(s, llr_cache, fold, weights, 0.0) for s in train]
        # Matriz de scores de entrenamiento para la escala 0 (baseline).
        z_train = np.asarray([[s["z"].get(sig, -1e9) for sig in sigmas] for s in train])
        y_train = np.asarray([1.0 if s["true"] and top1_train[i] in s["true"] else 0.0
                              for i, s in enumerate(train)])
        idx_train = np.asarray([max(0, sigmas.index(top1_train[i]))
                                if top1_train[i] in sigmas else 0
                                for i in range(len(train))])
        best_base = (float("inf"), 1.0)
        for temp in grid:
            probs = C._softmax_matrix(z_train, temp)
            conf = probs[np.arange(len(train)), idx_train]
            nll = _binary_nll(conf, y_train)
            if nll < best_base[0]:
                best_base = (nll, float(temp))
        # Baseline en test (misma temperatura); no se mezcla con entrenamiento.
        z_test = np.asarray([[s["z"].get(sig, -1e9) for sig in sigmas] for s in test])
        top1_test = [predict_top1(s, llr_cache, fold, weights, 0.0) for s in test]
        y_test = np.asarray([1.0 if s["true"] and top1_test[i] in s["true"] else 0.0
                             for i, s in enumerate(test)])
        idx_test = np.asarray([max(0, sigmas.index(top1_test[i]))
                               if top1_test[i] in sigmas else 0
                               for i in range(len(test))])
        base_conf.extend(float(v) for v in (
            C._softmax_matrix(z_test, best_base[1])[np.arange(len(test)), idx_test]))
        base_y.extend(float(v) for v in y_test)

        # Contexto: elegir escala y temperatura en entrenamiento.
        best_ctx = None
        for scale in scales:
            scores = np.asarray([
                [s["z"].get(sig, -1e9) + scale * sum(
                    float(weights.get(name, 0.0))
                    * float(llr_cache[s["key"]].get(sig, {}).get(name, 0.0))
                    for name in weights if weights.get(name, 0.0))
                 for sig in sigmas] for s in train])
            for temp in grid:
                probs = C._softmax_matrix(scores, temp)
                conf = probs[np.arange(len(train)), idx_train]
                nll = _binary_nll(conf, y_train)
                if best_ctx is None or nll < best_ctx[0]:
                    best_ctx = (nll, float(temp), float(scale))
        if best_ctx is None:
            best_ctx = (0.0, 1.0, 0.0)
        result["selected"][str(f)] = best_ctx[2]
        result["temperatures"][str(f)] = best_ctx[1]
        scores_test = np.asarray([
            [s["z"].get(sig, -1e9) + best_ctx[2] * sum(
                float(weights.get(name, 0.0))
                * float(llr_cache[s["key"]].get(sig, {}).get(name, 0.0))
                for name in weights if weights.get(name, 0.0))
             for sig in sigmas] for s in test])
        ctx_conf.extend(float(v) for v in (
            C._softmax_matrix(scores_test, best_ctx[1])[np.arange(len(test)), idx_test]))
        ctx_y.extend(float(v) for v in y_test)
    result["baseline"] = _fixed_label_calibration(np.asarray(base_conf),
                                                  np.asarray(base_y))
    result["context"] = _fixed_label_calibration(np.asarray(ctx_conf),
                                                 np.asarray(ctx_y))
    return result


def predict_top1(site: dict[str, Any],
                 llr_cache: dict[tuple[str, int, str], dict[str, dict[str, float]]],
                 fold: dict[tuple[str, int, str], int],
                 weights: dict[str, float], scale: float) -> str:
    final = {
        sigma: zv + scale * sum(
            float(weights.get(name, 0.0)) * float(llr_cache[site["key"]].get(sigma, {}).get(name, 0.0))
            for name in weights if weights.get(name, 0.0))
        for sigma, zv in site["z"].items()
    }
    return C.select_top1(final, site["z"], None) if final else ""


# --------------------------------------------------------------------------
# Salidas
# --------------------------------------------------------------------------
def _write_tsv(path: Path, rows: Iterable[dict[str, Any]], cols: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            vals = []
            for col in cols:
                value = row.get(col, "")
                if isinstance(value, (list, tuple, set)):
                    value = ",".join(str(x) for x in value)
                elif isinstance(value, dict):
                    value = json.dumps(value, ensure_ascii=False, sort_keys=True)
                vals.append(str(value).replace("\t", " ").replace("\n", " "))
            fh.write("\t".join(vals) + "\n")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8")


def _fmt(value: Any, digits: int = 3) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _pct(value: Any) -> str:
    return "—" if value is None else f"{100 * float(value):.1f} %"


def _report(result: dict[str, Any]) -> str:
    base = result["baseline"]
    ctx = result["context"]
    delta = result["delta"]
    lines = [
        f"# Contexto sigma — {result['genome']}", "",
        f"PWM: `{result['pwm']}` · ventana −{result['config']['up']}/+{result['config']['down']} · "
        f"{result['config']['decoys']} señuelos · {result['config']['folds']} folds · "
        f"semilla {result['config']['seed']}",
        f"Secuencia: `{result['data']['sequence_source']}`",
        "", "## Antes / después", "",
        "| capa | top-1 | top-3 cualquiera | top-3 todos | mayoría trivial |",
        "|---|---:|---:|---:|---:|",
        f"| Motivo (línea base) | {_pct(base['top1'])} | {_pct(base['top3_any'])} | "
        f"{_pct(base['top3_all'])} | {_pct(base['majority_baseline'])} ({base['majority_sigma']}) |",
        f"| Contexto (escala elegida) | {_pct(ctx['top1'])} | {_pct(ctx['top3_any'])} | "
        f"{_pct(ctx['top3_all'])} | {_pct(ctx['majority_baseline'])} ({ctx['majority_sigma']}) |",
        f"| Δ top-1 | {_fmt(delta['top1'] * 100, 2)} pp | {_fmt(delta['top3_any'] * 100, 2)} pp | "
        f"{_fmt(delta['top3_all'] * 100, 2)} pp | — |",
        "",
        f"Escalas elegidas por fold: {result['selected_scales']}. "
        "0.0 significa capa apagada (no mejoró en entrenamiento).",
        "",
        ("IC bootstrap del Δ top-1: "
         f"[{_fmt(delta['bootstrap']['ci_low'] * 100, 2)}, "
         f"{_fmt(delta['bootstrap']['ci_high'] * 100, 2)}] pp "
         f"(p bilateral ≈ {_fmt(delta['bootstrap']['p_two_sided'], 3)}).")
        if delta['bootstrap'].get('ci_low') is not None
        else "IC bootstrap: —",
        "", "## Calibración (top-1)", "",
        "| capa | ECE | Brier | NLL | confianza media |",
        "|---|---:|---:|---:|---:|",
        f"| Motivo | {_fmt(base['calibration']['ece'], 4)} | "
        f"{_fmt(base['calibration']['brier'], 4)} | {_fmt(base['calibration']['nll'], 4)} | "
        f"{_fmt(base['calibration']['mean_confidence'], 3)} |",
        f"| Contexto | {_fmt(ctx['calibration']['ece'], 4)} | "
        f"{_fmt(ctx['calibration']['brier'], 4)} | {_fmt(ctx['calibration']['nll'], 4)} | "
        f"{_fmt(ctx['calibration']['mean_confidence'], 3)} |",
        "", "## Ablación (escala fija = 1,0, fuera de muestra)", "",
        "| evidencia | top-1 | Δ top-1 | top-3 | cambios | aciertos | errores |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["ablation"]:
        lines.append(
            f"| {row['label']} | {_pct(row['top1'])} | "
            f"{_fmt(row['delta_top1'] * 100, 2)} pp | {_pct(row['top3_any'])} | "
            f"{row['changed']} | {row['good']} | {row['bad']} |")
    lines += ["", "## Barrido de escala (fuera de muestra)", "",
              "| escala | top-1 | Δ top-1 | IC 95% Δ | top-3 |",
              "|---:|---:|---:|---:|---:|"]
    for row in result["sweep"]:
        boot = row.get("bootstrap", {})
        ci = ("—" if boot.get("ci_low") is None
              else f"[{_fmt(boot['ci_low'] * 100, 2)}, {_fmt(boot['ci_high'] * 100, 2)}]")
        lines.append(f"| {_fmt(row['scale'], 2)} | {_pct(row['top1'])} | "
                     f"{_fmt(row['delta_top1'] * 100, 2)} pp | {ci} | "
                     f"{_pct(row['top3_any'])} |")
    lines += ["", "## Por sigma (línea base → contexto)", "",
              "| sigma | support | P top-1 | R top-1 | F1 top-1 | P top-3 | R top-3 | F1 top-3 |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for sigma in sorted(base["per_sigma"]):
        b = base["per_sigma"][sigma]
        c = ctx["per_sigma"].get(sigma, {})
        lines.append(
            f"| {sigma} | {b['support']} | "
            f"{_fmt(b['top1']['precision'])} → {_fmt(c.get('top1', {}).get('precision'))} | "
            f"{_fmt(b['top1']['recall'])} → {_fmt(c.get('top1', {}).get('recall'))} | "
            f"{_fmt(b['top1']['f1'])} → {_fmt(c.get('top1', {}).get('f1'))} | "
            f"{_fmt(b['top3']['precision'])} → {_fmt(c.get('top3', {}).get('precision'))} | "
            f"{_fmt(b['top3']['recall'])} → {_fmt(c.get('top3', {}).get('recall'))} | "
            f"{_fmt(b['top3']['f1'])} → {_fmt(c.get('top3', {}).get('f1'))} |")
    lines += ["", "## Calibración sin cambiar top-1 (exploratorio)", "",
              "La escala y la temperatura se eligen en entrenamiento para calibrar",
              "la confianza del top-1 del motivo; la predicción no cambia.", "",
              "| modo | ECE | Brier | NLL | confianza media | acierto |",
              "|---|---:|---:|---:|---:|---:|"]
    co = result.get("confidence_only", {})
    for mode in ("baseline", "context"):
        m = co.get(mode, {})
        lines.append(f"| {mode} | {_fmt(m.get('ece'), 4)} | {_fmt(m.get('brier'), 4)} | "
                     f"{_fmt(m.get('nll'), 4)} | {_fmt(m.get('mean_confidence'), 3)} | "
                     f"{_pct(m.get('accuracy'))} |")
    lines.append(f"\nEscalas elegidas: {co.get('selected')}; "
                 f"temperaturas: {co.get('temperatures')}.")
    lines += ["", "## Multietiqueta", "",
              f"- TSS evaluables: {base['n_evaluable']} "
              f"({base['n_single']} single, {base['n_multi']} multi)",
              f"- top-1 cualquiera: {_pct(base['top1'])} → {_pct(ctx['top1'])}",
              f"- top-3 cualquiera: {_pct(base['top3_any'])} → {_pct(ctx['top3_any'])}",
              f"- top-3 todos: {_pct(base['top3_all'])} → {_pct(ctx['top3_all'])}",
              f"- multi top-1 (alguno): {_pct(base['multi_top1_any'])} → {_pct(ctx['multi_top1_any'])}",
              f"- multi top-1 exacto: {_pct(base['multi_exact_top1'])} → {_pct(ctx['multi_exact_top1'])}",
              f"- Jaccard top-3: {_fmt(base['jaccard_top3'])} → {_fmt(ctx['jaccard_top3'])}",
              "", "## σ54 con geometría propia", ""]
    s54 = result["sigma54"]
    lines += [
        f"- σ54 verdaderos: {s54['n_true']} · predichos top-1: {s54['n_pred']}",
        f"- recall top-1 / top-3: {_pct(s54['recall_top1'])} / {_pct(s54['recall_top3'])} "
        f"(motivo: {_pct(s54['motif_recall_top1'])} / {_pct(s54['motif_recall_top3'])})",
        f"- precisión σ54: {_pct(s54['precision'])}",
        f"- geometría en verdaderos: {s54['geometry_true']}",
        f"- geometría en predichos: {s54['geometry_pred']}",
        f"- nota: {s54['note']}",
        "", "## Evidencias: fuente, cobertura y peso", "",
        "| evidencia | disponible | cobertura | peso activo | fuente |",
        "|---|---:|---:|---:|---|",
    ]
    for name, spec in result["evidence"]["specs"].items():
        lines.append(
            f"| {name} | {result['evidence']['available'].get(name)} | "
            f"{_pct(result['evidence']['coverage'].get(name))} | "
            f"{_fmt(result['evidence']['weights'].get(name), 2)} | `{spec['source']}` |")
    gate = result["evidence"].get("geometry_gate", {})
    lines += ["", f"Puerta de geometría σ54: PWM n_sites={gate.get('pwm_sites_sigma54')}, "
              f"mínimo {gate.get('min_sites')} → activa={gate.get('active')}. "
              "Si está inactiva se reporta como diagnóstico y no mueve el ranking.",
              "", "## Límites", "",
              "- El gold es DBTBS/BSGatlas (B. subtilis) y RegulonDB (K-12); no hay",
              "  negativos experimentales. Un TSS fuera del gold no es falso positivo.",
              "- La validación es cruzada dentro de cada especie; no es validación",
              "  experimental ni transferencia a un genoma nuevo.",
              "- En B. subtilis, la anotación de genes y el gold comparten procedencia",
              "  (BSGatlas/DBTBS). La mejora es válida dentro de esa anotación.",
              "- La expresión génica no está disponible; `expression` queda en 0.",
              "- La estructura/dominio no cubre estos dos genomas; `domain` queda en 0.",
              "- La escala de contexto se elige en entrenamiento y puede apagarse.",
              "- σ54 de B. subtilis usa una PWM de una caja: la geometría se aproxima",
              "  por posición, no por dos cajas −24/−12 reales.",
              ""]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Evaluación de un genoma
# --------------------------------------------------------------------------
def evaluate_genome(data_root: Path, genome: str, out_root: Path,
                    pwm_path: Path | None = None, *, up: int = DEFAULT_UP,
                    down: int = DEFAULT_DOWN, decoys: int = DEFAULT_DECOYS,
                    seed: int = DEFAULT_SEED, folds: int = DEFAULT_FOLDS,
                    gap: float = DEFAULT_GAP,
                    scales: Sequence[float] = DEFAULT_SCALES,
                    min_coverage: float = DEFAULT_MIN_COVERAGE,
                    min_support: int = DEFAULT_MIN_SUPPORT,
                    min_gain: float = 0.01,
                    curated_context: bool = False) -> dict[str, Any]:
    gold_path = data_root / "sigma" / "gold" / genome / "tss_sigma.tsv"
    if not gold_path.exists():
        raise FileNotFoundError(f"falta {gold_path}")
    pwm = resolve_pwm(data_root, genome, pwm_path)
    if not pwm.exists():
        raise FileNotFoundError(f"falta PWM {pwm}")
    seqs, seq_source = load_sequences(data_root, genome)
    if not seqs:
        raise ValueError(f"sin secuencia local para {genome}")
    pwms = load_pwm_library(data_root, pwm)
    if not pwms:
        raise ValueError(f"sin PWMs válidas en {pwm}")

    gold_rows = LAYERS.load_gold(gold_path, "gold")
    scored = score_gold(gold_rows, seqs, pwms, up, down, decoys, seed)

    sites: list[dict[str, Any]] = []
    for row in gold_rows:
        key = (row["contig"], int(row["tss"]), row["strand"])
        sc = scored[key]
        sites.append({
            "key": key,
            "contig": key[0], "tss": key[1], "strand": key[2],
            "true": set(row["sigmas"]),
            "z": {sigma: float(h.get("z", 0.0)) for sigma, h in sc["scores"].items()},
            "hits": sc["scores"],
            "sequence_available": sc["sequence_available"],
        })
    keys = [s["key"] for s in sites]
    fold = blocked_folds(keys, folds)

    ctx = C.build_context(data_root, genome, keys, use_curated=curated_context)
    geom = 0
    for site in sites:
        hit = site["hits"].get("sigma54")
        category = C.sigma54_geometry_category(hit, genome) if hit else "none"
        ctx.sites[site["key"]]["sigma54_geometry"] = category
        if category != "none":
            geom += 1
    ctx.coverage["sigma54_geometry"] = geom / len(sites) if sites else 0.0
    ctx.available["sigma54_geometry"] = "sigma54" in pwms and geom > 0

    weights = C.active_weights(ctx, min_coverage=min_coverage)
    # Regla a priori: la geometría σ54 solo se usa como prior si la PWM tiene
    # >= 20 sitios de entrenamiento. B. subtilis tiene 11 (matriz de una caja,
    # inestable); se reporta como diagnóstico, no mueve el ranking.
    pwm_sites54 = int(pwms.get("sigma54", {}).get("n_sites", 0) or 0)
    geometry_active = bool(ctx.available.get("sigma54_geometry", False)
                           and pwm_sites54 >= 20)
    if not geometry_active:
        weights["sigma54_geometry"] = 0.0
    active_names = [name for name, weight in weights.items() if weight > 0.0]

    # Modelos fuera de muestra: uno por fold, ajustado solo con entrenamiento.
    models: dict[int, C.LogOddsModel] = {}
    for f in range(folds):
        train = {s["key"]: set(s["true"]) for s in sites
                 if fold[s["key"]] != f and s["true"]}
        models[f] = C.fit_logodds(ctx, train, evidence_names=active_names,
                                  min_support=min_support)

    # El LLR crudo se calcula una vez por TSS x sigma x evidencia.
    llr_cache = precompute_llr(sites, ctx, models, fold, active_names)

    # Escala del contexto elegida en entrenamiento (puede quedar en 0).
    selected = select_scales(sites, ctx, llr_cache, fold, folds, weights, scales,
                             gap, min_gain=min_gain)

    # Configuraciones.
    configs: list[dict[str, Any]] = []
    predictions: dict[str, dict[tuple[str, int, str], dict[str, Any]]] = {}
    predictions["baseline"] = predict(sites, ctx, llr_cache, fold, {},
                                      fixed_scales(0.0, folds), None)
    predictions["context"] = predict(sites, ctx, llr_cache, fold, weights, selected, gap)
    for name, weight in weights.items():
        if weight <= 0.0:
            continue
        predictions[f"only:{name}"] = predict(
            sites, ctx, llr_cache, fold, {name: weight}, fixed_scales(1.0, folds), gap)
        rest = {k: v for k, v in weights.items() if k != name}
        predictions[f"without:{name}"] = predict(
            sites, ctx, llr_cache, fold, rest, fixed_scales(1.0, folds), gap)
    predictions["all_fixed_1.0"] = predict(
        sites, ctx, llr_cache, fold, weights, fixed_scales(1.0, folds), gap)
    for scale in scales:
        predictions[f"sweep:{float(scale):.2f}"] = predict(
            sites, ctx, llr_cache, fold, weights, fixed_scales(float(scale), folds), gap)

    # Calibración.
    base_train = _temperature_records(sites, ctx, llr_cache, fold, {}, fixed_scales(0.0, folds), "motif")
    ctx_train = _temperature_records(sites, ctx, llr_cache, fold, weights, selected, "context")
    calibration: dict[str, dict[str, Any]] = {
        "baseline": calibrated(sites, predictions["baseline"], base_train, fold, SIGMA_ORDER),
        "context": calibrated(sites, predictions["context"], ctx_train, fold, SIGMA_ORDER),
    }
    for name, pred in predictions.items():
        if name in calibration or not name.startswith("only:"):
            continue
        evidence = name.split(":", 1)[1]
        train_records = _temperature_records(
            sites, ctx, llr_cache, fold, {evidence: weights[evidence]},
            fixed_scales(1.0, folds), "context")
        calibration[name] = calibrated(sites, pred, train_records, fold, SIGMA_ORDER)

    baseline = summarize(sites, predictions["baseline"])
    baseline["calibration"] = calibration["baseline"]
    context_summary = summarize(sites, predictions["context"])
    context_summary["calibration"] = calibration["context"]
    # Corrección de nombres: summarize usa top1/top3_any; se conservan ambos.
    delta = {
        "top1": (context_summary["top1"] or 0.0) - (baseline["top1"] or 0.0),
        "top3_any": (context_summary["top3_any"] or 0.0) - (baseline["top3_any"] or 0.0),
        "top3_all": (context_summary["top3_all"] or 0.0) - (baseline["top3_all"] or 0.0),
    }
    evaluable = [s for s in sites if _evaluable(s)]
    correct_base = [1 if predictions["baseline"][s["key"]]["top1"] in s["true"] else 0
                    for s in evaluable]
    correct_ctx = [1 if predictions["context"][s["key"]]["top1"] in s["true"] else 0
                   for s in evaluable]
    delta["bootstrap"] = C.bootstrap_delta(correct_base, correct_ctx, n_boot=1000, seed=seed)

    changed = good = bad = neutral = promoted_outside_top3 = 0
    for site in evaluable:
        key = site["key"]
        b = predictions["baseline"][key]["top1"]
        c = predictions["context"][key]["top1"]
        if b == c:
            continue
        changed += 1
        if c in site["true"]:
            good += 1
        elif b in site["true"]:
            bad += 1
        else:
            neutral += 1
        if c not in predictions["baseline"][key]["motif_top3"]:
            promoted_outside_top3 += 1

    ablation: list[dict[str, Any]] = []
    for name, pred in predictions.items():
        if not (name.startswith("only:") or name.startswith("without:") or name == "all_fixed_1.0"):
            continue
        summary = summarize(sites, pred)
        row_changed = sum(1 for site in evaluable
                          if pred[site["key"]]["top1"] != predictions["baseline"][site["key"]]["top1"])
        row_good = sum(1 for site in evaluable
                       if pred[site["key"]]["top1"] in site["true"]
                       and predictions["baseline"][site["key"]]["top1"] not in site["true"])
        row_bad = sum(1 for site in evaluable
                      if pred[site["key"]]["top1"] not in site["true"]
                      and predictions["baseline"][site["key"]]["top1"] in site["true"])
        label = ("todas (1,0)" if name == "all_fixed_1.0"
                 else name.replace("only:", "solo ").replace("without:", "sin "))
        ablation.append({
            "name": name, "label": label,
            "top1": summary["top1"], "top3_any": summary["top3_any"],
            "delta_top1": (summary["top1"] or 0.0) - (baseline["top1"] or 0.0),
            "changed": row_changed, "good": row_good, "bad": row_bad,
            "per_sigma": summary["per_sigma"],
        })

    sweep: list[dict[str, Any]] = []
    for scale in scales:
        name = f"sweep:{float(scale):.2f}"
        summary = summarize(sites, predictions[name])
        correct_sw = [1 if predictions[name][s["key"]]["top1"] in s["true"] else 0
                      for s in evaluable]
        sweep.append({
            "scale": float(scale),
            "top1": summary["top1"], "top3_any": summary["top3_any"],
            "delta_top1": (summary["top1"] or 0.0) - (baseline["top1"] or 0.0),
            "bootstrap": C.bootstrap_delta(correct_base, correct_sw,
                                           n_boot=1000, seed=seed),
        })

    # σ54.
    true54 = [s for s in sites if "sigma54" in s["true"]]
    pred54 = [s for s in evaluable if predictions["context"][s["key"]]["top1"] == "sigma54"]
    motif_pred54 = [s for s in evaluable if predictions["baseline"][s["key"]]["top1"] == "sigma54"]
    geom_true = Counter(ctx.row(s["key"]).get("sigma54_geometry", "none") for s in true54)
    geom_pred = Counter(ctx.row(s["key"]).get("sigma54_geometry", "none") for s in pred54)
    motif_recall1 = (sum(1 for s in true54
                         if predictions["baseline"][s["key"]]["top1"] == "sigma54")
                     / len(true54)) if true54 else None
    motif_recall3 = (sum(1 for s in true54
                         if "sigma54" in predictions["baseline"][s["key"]]["top3"])
                     / len(true54)) if true54 else None
    recall1 = (sum(1 for s in true54
                   if predictions["context"][s["key"]]["top1"] == "sigma54")
               / len(true54)) if true54 else None
    recall3 = (sum(1 for s in true54
                   if "sigma54" in predictions["context"][s["key"]]["top3"])
               / len(true54)) if true54 else None
    precision = (sum(1 for s in pred54 if "sigma54" in s["true"]) / len(pred54)
                 if pred54 else None)
    if genome == "subtilis":
        note = ("PWM de una caja; geometría aproximada por posición del extremo "
                "(9–32 nt desde el TSS). No hay dos cajas −24/−12 reales.")
    else:
        note = ("PWM de dos cajas; el escáner exige espaciador 10–11 nt (−24/−12). "
                "La geometría canónica es parte del score.")
    sigma54 = {
        "n_true": len(true54), "n_pred": len(pred54),
        "n_motif_pred": len(motif_pred54),
        "precision": precision, "recall_top1": recall1, "recall_top3": recall3,
        "motif_recall_top1": motif_recall1, "motif_recall_top3": motif_recall3,
        "geometry_true": dict(geom_true.most_common()),
        "geometry_pred": dict(geom_pred.most_common()),
        "note": note,
    }

    # Tablas descriptivas de evidencia (ajuste completo; no se usan en test).
    full_labels = {s["key"]: set(s["true"]) for s in sites if s["true"]}
    full_model = C.fit_logodds(ctx, full_labels, evidence_names=active_names,
                               min_support=min_support)
    evidence_rows = [dict(row, genome=genome, fit="all")
                     for row in full_model.to_rows()]

    conf_only = confidence_only(sites, ctx, llr_cache, fold, folds, weights,
                                scales, SIGMA_ORDER)

    out_dir = out_root / genome
    out_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "version": VERSION,
        "genome": genome,
        "pwm": str(pwm),
        "ecoli_transfer": is_ecoli_transfer(genome, pwm),
        "config": {
            "up": up, "down": down, "decoys": decoys, "seed": seed,
            "folds": folds, "gap": gap, "scales": list(map(float, scales)),
            "min_coverage": min_coverage, "min_support": min_support,
            "min_gain": min_gain,
            "curated_context": curated_context,
        },
        "data": {
            "n_gold_tss": len(sites),
            "n_evaluable": len(evaluable),
            "sequence_source": seq_source,
            "gold": str(gold_path),
        },
        "evidence": {
            "specs": {name: {"source": spec.source, "description": spec.description,
                             "weight": spec.weight, "applies_to": list(spec.applies_to)}
                      for name, spec in C.EVIDENCE_SPECS.items()},
            "available": dict(ctx.available),
            "coverage": {k: round(float(v), 6) for k, v in ctx.coverage.items()},
            "weights": dict(weights),
            "mode": ctx.sources.get("_mode", ""),
            "active": active_names,
            "geometry_gate": {
                "pwm_sites_sigma54": pwm_sites54,
                "min_sites": 20,
                "active": geometry_active,
                "note": "Geometría σ54 como prior solo si la PWM tiene >= 20 sitios.",
            },
        },
        "baseline": baseline,
        "context": context_summary,
        "delta": delta,
        "selected_scales": {str(f): float(s) for f, s in sorted(selected.items())},
        "ablation": ablation,
        "sweep": sweep,
        "changed": {"n_changed": changed, "good": good, "bad": bad,
                    "neutral": neutral,
                    "promoted_outside_motif_top3": promoted_outside_top3},
        "sigma54": sigma54,
        "classes_baseline": baseline["classes"],
        "classes_context": context_summary["classes"],
        "calibration": calibration,
        "confidence_only": conf_only,
    }

    # TSV de sitios.
    site_rows = []
    for site in sites:
        key = site["key"]
        b = predictions["baseline"][key]
        c = predictions["context"][key]
        site_rows.append({
            "genome": genome,
            "contig": site["contig"], "tss": site["tss"], "strand": site["strand"],
            "true_sigmas": sorted(site["true"]), "n_true": len(site["true"]),
            "motif_top1": b["top1"], "motif_top3": b["top3"],
            "context_top1": c["top1"], "context_top3": c["top3"],
            "scale_selected": c["scale"],
            "class_baseline": _classification(site, b["top1"], b["top3"]),
            "class_context": _classification(site, c["top1"], c["top3"]),
            "changed": "yes" if b["top1"] != c["top1"] else "no",
            "delta_correct": (int(c["top1"] in site["true"])
                              - int(b["top1"] in site["true"])),
            "calibration_prob_baseline": None,
            "calibration_prob_context": None,
        })
    # Probabilidades calibradas por sitio (con la temperatura de su fold).
    temp_b = calibration["baseline"].get("temperatures", {})
    temp_c = calibration["context"].get("temperatures", {})
    row_by_key = {(r["contig"], r["tss"], r["strand"]): r for r in site_rows}
    for site in sites:
        key = site["key"]
        f = str(fold[key])
        if not _evaluable(site):
            continue
        pb = C.softmax(predictions["baseline"][key]["final"], float(temp_b.get(f, 1.0)),
                       SIGMA_ORDER)
        pc = C.softmax(predictions["context"][key]["final"], float(temp_c.get(f, 1.0)),
                       SIGMA_ORDER)
        row = row_by_key[key]
        row["calibration_prob_baseline"] = round(max(pb.values()), 6)
        row["calibration_prob_context"] = round(max(pc.values()), 6)

    _write_tsv(out_dir / "context_sites.tsv", site_rows,
               ["genome", "contig", "tss", "strand", "true_sigmas", "n_true",
                "motif_top1", "motif_top3", "context_top1", "context_top3",
                "scale_selected", "class_baseline", "class_context", "changed",
                "delta_correct", "calibration_prob_baseline",
                "calibration_prob_context"])

    # TSV por evidencia (una fila por TSS x sigma).
    pred_rows = []
    for site in sites:
        key = site["key"]
        c = predictions["context"][key]
        for sigma, z in sorted(site["z"].items()):
            llr_raw = sum(c["llr"].get(sigma, {}).values())
            llr_weighted = sum(float(weights.get(name, 0.0))
                               * float(c["llr"].get(sigma, {}).get(name, 0.0))
                               for name in active_names)
            row = {
                "genome": genome, "contig": key[0], "tss": key[1], "strand": key[2],
                "sigma": sigma, "true": "yes" if sigma in site["true"] else "no",
                "motif_z": round(float(z), 6),
                "llr_total": round(float(llr_raw), 6),
                "llr_weighted": round(llr_weighted, 6),
                "scale_selected": c["scale"],
                "final": round(float(c["final"].get(sigma, 0.0)), 6),
                "is_motif_top1": "yes" if sigma == c["motif_top1"] else "no",
                "is_context_top1": "yes" if sigma == c["top1"] else "no",
                "in_context_top3": "yes" if sigma in c["top3"] else "no",
            }
            for name in active_names:
                row[f"llr_{name}"] = round(float(c["llr"].get(sigma, {}).get(name, 0.0)), 6)
            pred_rows.append(row)
    pred_cols = (["genome", "contig", "tss", "strand", "sigma", "true", "motif_z",
                  "llr_total", "llr_weighted", "scale_selected", "final",
                  "is_motif_top1", "is_context_top1", "in_context_top3"]
                 + [f"llr_{n}" for n in active_names])
    _write_tsv(out_dir / "context_predictions.tsv", pred_rows, pred_cols)

    _write_tsv(out_dir / "context_evidence.tsv", evidence_rows,
               ["genome", "evidence", "sigma", "value", "log_odds", "n_train",
                "min_support", "weight", "source", "fit"])
    coverage_rows = [
        {"genome": genome, "evidence": name, "available": ctx.available.get(name),
         "coverage": round(float(ctx.coverage.get(name, 0.0)), 6),
         "weight": weights.get(name, 0.0), "source": spec.source,
         "description": spec.description}
        for name, spec in C.EVIDENCE_SPECS.items()
    ]
    _write_tsv(out_dir / "context_evidence_coverage.tsv", coverage_rows,
               ["genome", "evidence", "available", "coverage", "weight",
                "source", "description"])

    _write_json(out_dir / "context_metrics.json", result)
    (out_dir / "context_report.md").write_text(_report(result), encoding="utf-8")
    return result


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def selftest() -> str:
    assert C.selftest()
    assert blocked_folds([("c", 10, "+"), ("c", 20, "+"), ("c", 30, "+")], 2) == {
        ("c", 10, "+"): 0, ("c", 20, "+"): 0, ("c", 30, "+"): 1}
    return "evaluador de contexto: folds y módulo OK"


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--pwm", type=Path, default=None,
                    help="PWM; por defecto la de la especie (B. subtilis: DBTBS).")
    ap.add_argument("--gold-root", type=Path, default=ROOT / "data" / "sigma" / "gold")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "sigma" / "context")
    ap.add_argument("--up", type=int, default=DEFAULT_UP)
    ap.add_argument("--down", type=int, default=DEFAULT_DOWN)
    ap.add_argument("--decoys", type=int, default=DEFAULT_DECOYS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--folds", type=int, default=DEFAULT_FOLDS)
    ap.add_argument("--gap", type=float, default=DEFAULT_GAP,
                    help="Margen z máximo que el contexto puede remontar (None = sin límite).")
    ap.add_argument("--scales", default="0,0.25,0.5,0.75,1.0",
                    help="Escalas de contexto a barrer/elegir.")
    ap.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE)
    ap.add_argument("--min-support", type=int, default=DEFAULT_MIN_SUPPORT)
    ap.add_argument("--min-gain", type=float, default=0.01,
                    help="Mejora mínima en entrenamiento para encender el contexto.")
    ap.add_argument("--curated-context", action="store_true",
                    help="B. subtilis: usar subtilis_tss_tts.tsv (sensibilidad).")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0
    data_root = ROOT / "data"
    if args.all and args.pwm is not None:
        ap.error("--all usa la PWM propia de cada especie; no combines --all con --pwm")
    scales = [float(x) for x in str(args.scales).split(",") if x.strip()]
    genomes = ["ecoli_k12", "subtilis"] if args.all else [args.genome]
    if not genomes or genomes == [None]:
        ap.error("indica --genome, --all o --selftest")
    ok = 0
    for genome in genomes:
        if genome is None:
            continue
        try:
            result = evaluate_genome(
                data_root, genome, args.out, args.pwm,
                up=args.up, down=args.down, decoys=args.decoys, seed=args.seed,
                folds=args.folds, gap=args.gap, scales=scales,
                min_coverage=args.min_coverage, min_support=args.min_support,
                min_gain=args.min_gain, curated_context=args.curated_context)
            b = result["baseline"]
            c = result["context"]
            print(f"{genome}: top1 {b['top1']:.4f} -> {c['top1']:.4f} · "
                  f"top3 {b['top3_any']:.4f} -> {c['top3_any']:.4f} · "
                  f"escalas {result['selected_scales']}")
            if result.get("ecoli_transfer"):
                print(f"  AVISO: {genome} evaluado con PWM de E. coli; "
                      "no usar como resultado principal.", file=sys.stderr)
            ok += 1
        except Exception as exc:  # noqa: BLE001
            print(f"FALLO {genome}: {exc}", file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
