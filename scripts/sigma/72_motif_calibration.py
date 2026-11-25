#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibración de la capa 2 (motivo) para sigma54 y sigma38.

El análisis responde tres preguntas, sin reentrenar el núcleo ni tocar la GUI:

1. ¿Ayuda recortar cada PWM a sus 6/7/8 posiciones de mayor información?
   Se enmascaran las columnas de baja información (la geometría de cajas y el
   espaciador no cambia) y se comparan AUC, precisión y recall del soporte en
   K-12 y en B. subtilis.
2. ¿Los umbrales de z, p y margen calibrados en K-12 transfieren a B. subtilis?
   Se elige un umbral por sigma y feature en K-12 y se aplica tal cual
   (sin reajustar) a B. subtilis.
3. ¿Qué regla simple de dos features funciona por sigma?
   Se buscan pares (z, p, margen, distancia al inicio) sobre K-12 y se evalúa
   la regla elegida en B. subtilis.

Convenciones:
  * Una fila del análisis es un TSS x sigma de
    ``data/sigma/layers/<genoma>/layer_sigma_features.tsv``.
  * ``true=yes`` es la etiqueta del gold. No se usan ``n_true``,
    ``gold_confidence`` ni ``source`` como señales.
  * ``top-1`` es la asignación actual del núcleo (``is_top1=yes``); se usa
    solo como baseline para calcular el lift.
  * El criterio de "soporte" de la capa 2 es ``pct_decoy >= 0.95 & z >= 1``,
    equivalente al que reporta ``50_eval_layers.py``.

Uso:
    .venv/bin/python scripts/sigma/72_motif_calibration.py
    .venv/bin/python scripts/sigma/72_motif_calibration.py --selftest
    .venv/bin/python scripts/sigma/72_motif_calibration.py --limit 100 --out /tmp/cal
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import json
import sys
import zlib
from itertools import combinations, product
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma.sigma_annotator import (  # noqa: E402
    SIGMA_SPACER, load_pwm_library, load_sequence, revcomp, score_sigma,
)

VERSION = "1.0.0"
DEFAULT_GENOMES = ("ecoli_k12", "subtilis")
DEFAULT_SIGMAS = ("sigma54", "sigma38")
DEFAULT_UP = 80
DEFAULT_DOWN = 20
DEFAULT_DECOYS = 200
DEFAULT_SEED = 17
K_VALUES = (6, 7, 8)
CORE_MODES = ("box", "joint")
SUPPORT_PCT = 0.95
SUPPORT_Z = 1.0
MIN_PRED = 5
MIN_TP = 3
MIN_RECALL = 0.20

#: Metadatos del gold que NUNCA deben entrar en una regla.
FORBIDDEN_COLUMNS = ("n_true", "gold_confidence", "source")

FEATURE_COLUMNS = {
    "z": "z",
    "p": "pvalue",
    "margin": "margin",
    "dist": "dist_to_start",
}
FEATURE_OPS = {"z": ">=", "p": "<=", "margin": ">=", "dist": "<="}
FEATURE_GRID: dict[str, tuple[float, ...]] = {
    "z": tuple(round(-1.0 + 0.25 * i, 3) for i in range(25)),          # -1..5
    "p": (0.20, 0.15, 0.10, 0.075, 0.06, 0.05, 0.04, 0.03, 0.025,
          0.02, 0.015, 0.01, 0.0075, 0.005),
    "margin": tuple(round(-1.0 + 0.25 * i, 3) for i in range(17)),     # -1..3
    "dist": (30, 40, 50, 60, 70, 80),
}
CORE_TSV_COLS = [
    "genome", "sigma", "variant", "mode", "k", "n_positions", "kept_up",
    "kept_down", "auc", "n_pred", "tp", "fp", "fn", "precision", "recall",
    "f1", "base_rate", "n_true", "n_tss",
]
CAL_TSV_COLS = [
    "sigma", "feature", "objective", "direction", "threshold",
    "k12_n_pred", "k12_tp", "k12_fp", "k12_precision", "k12_recall",
    "k12_f1", "k12_lift_vs_top1", "k12_base_rate",
    "sub_n_pred", "sub_tp", "sub_fp", "sub_precision", "sub_recall",
    "sub_f1", "sub_lift_vs_top1", "sub_base_rate",
]
TWO_TSV_COLS = [
    "sigma", "objective", "selected", "rule", "k12_n_pred", "k12_tp",
    "k12_fp", "k12_precision", "k12_recall", "k12_f1", "k12_lift_vs_top1",
    "sub_n_pred", "sub_tp", "sub_fp", "sub_precision", "sub_recall",
    "sub_f1", "sub_lift_vs_top1",
]


# --------------------------------------------------------------------------
# Utilidades
# --------------------------------------------------------------------------
def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


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


def load_genome_sequences(data_root: Path, genome: str) -> dict[str, str]:
    """Secuencia del genoma: primero el cargador del núcleo, luego FASTA."""
    seqs = load_sequence(data_root, genome)
    if seqs:
        return seqs
    candidates = [
        data_root / "reference_genomes" / genome / f"{genome}_genomic.fna.gz",
        data_root / "reference_genomes" / genome / f"{genome}_genomic.fna",
        data_root / "annotations" / f"{genome}_genomic.fna.gz",
        data_root / "annotations" / f"{genome}_genomic.fna",
        data_root / "raw" / f"{genome}_genomic.fna.gz",
        data_root / "raw" / f"{genome}_genomic.fna",
    ]
    if genome.lower() == "subtilis":
        candidates.append(data_root / "dbtbs" / "genome.fna")
    for path in candidates:
        if path.exists():
            got = _read_fasta(path)
            if got:
                return got
    return {}


def read_features(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Falta {path}; corré primero scripts/sigma/50_eval_layers.py")
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    required = {"contig", "tss", "strand", "sigma", "true", "is_top1",
                "z", "pvalue", "margin", "dist_to_start"}
    missing = required - set(rows[0] if rows else {})
    if missing:
        raise ValueError(f"{path} no tiene columnas requeridas: {sorted(missing)}")
    return rows


def _window(seq: str, tss: int, strand: str, up: int, down: int) -> str:
    """Ventana orientada al TSS, igual que 50_eval_layers._window."""
    if not seq:
        return ""
    p = tss - 1
    if strand == "+":
        a = max(0, p - up)
        b = min(len(seq), p + down + 1)
        return seq[a:b]
    a = max(0, p - down)
    b = min(len(seq), p + up + 1)
    return revcomp(seq[a:b])


def _limit_rows(rows: list[dict[str, str]], limit: int) -> list[dict[str, str]]:
    """Recorta el número de TSS únicos (solo para smoke tests)."""
    if limit <= 0:
        return rows
    seen: set[tuple[str, str, str]] = set()
    keep: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (row.get("contig", ""), row.get("tss", ""), row.get("strand", ""))
        if key not in seen:
            seen.add(key)
            keep.add(key)
            if len(keep) >= limit:
                break
    return [r for r in rows
            if (r.get("contig", ""), r.get("tss", ""), r.get("strand", "")) in keep]


# --------------------------------------------------------------------------
# PWM core
# --------------------------------------------------------------------------
def position_information(pwm: np.ndarray) -> np.ndarray:
    """Información por posición (bits) a partir del log-odds contra uniforme."""
    arr = np.asarray(pwm, dtype=float)
    probs = 0.25 * np.power(2.0, arr)
    return (probs * arr).sum(axis=1)


def build_core_pwm(rec: dict[str, Any], k: int,
                   mode: str = "box") -> tuple[dict[str, Any], dict[str, list[int]]]:
    """Devuelve una copia del modelo con solo las k posiciones de más información.

    Las columnas descartadas se ponen a log-odds 0 (todas las bases equiprobables),
    de modo que el ancho de caja y la geometría del espaciador no cambian.

    mode="box": top-k dentro de cada caja (si k >= ancho, la caja queda intacta).
    mode="joint": top-k entre las dos cajas juntas.
    """
    if k <= 0:
        raise ValueError("k debe ser positivo")
    if mode not in CORE_MODES:
        raise ValueError(f"mode inválido: {mode}")
    up_pwm = np.array(rec["up"]["pwm"], dtype=float, copy=True)
    down_pwm = np.array(rec["down"]["pwm"], dtype=float, copy=True)
    up_info = position_information(up_pwm)
    down_info = position_information(down_pwm)
    pivot = len(up_info)
    masks = {
        "up": np.zeros(len(up_info), dtype=bool),
        "down": np.zeros(len(down_info), dtype=bool),
    }
    if mode == "box":
        for box, info in (("up", up_info), ("down", down_info)):
            take = min(k, len(info))
            order = np.argsort(-info, kind="stable")[:take]
            masks[box][order] = True
    else:
        entries = ([("up", i, float(up_info[i])) for i in range(len(up_info))]
                   + [("down", i, float(down_info[i])) for i in range(len(down_info))])
        entries.sort(key=lambda e: (-e[2], e[0], e[1]))
        for box, i, _ in entries[:k]:
            masks[box][i] = True
    for box, pwm in (("up", up_pwm), ("down", down_pwm)):
        pwm[~masks[box]] = 0.0
    new = dict(rec)
    new["up"] = dict(rec["up"])
    new["down"] = dict(rec["down"])
    new["up"]["pwm"] = up_pwm
    new["down"]["pwm"] = down_pwm
    kept = {
        "up": [i + 1 for i in range(len(masks["up"])) if masks["up"][i]],
        "down": [i + 1 for i in range(len(masks["down"])) if masks["down"][i]],
    }
    _ = pivot  # documenta que el offset entre cajas no se toca
    return new, kept


def _fmt_positions(kept: dict[str, list[int]]) -> str:
    return ("up:" + ",".join(str(i) for i in kept["up"])
            + "|down:" + ",".join(str(i) for i in kept["down"]))


# --------------------------------------------------------------------------
# Null y métricas
# --------------------------------------------------------------------------
def _decoy_windows(seqs: dict[str, str], sigma: str, seed: int, n: int,
                   length: int) -> list[str]:
    """Ventanas señuelo con el mismo esquema que 50_eval_layers._decoy_scores."""
    rng = np.random.default_rng(seed + zlib.crc32(sigma.encode()) % 1000)
    contigs = [s for s in seqs.values() if len(s) >= length]
    if not contigs:
        return []
    out: list[str] = []
    for _ in range(n):
        seq = contigs[int(rng.integers(0, len(contigs)))]
        start = int(rng.integers(0, len(seq) - length + 1))
        out.append(seq[start:start + length])
    return out


def _null_distribution(windows: Sequence[str], model: dict[str, Any]) -> dict[str, Any]:
    vals = []
    for window in windows:
        hit = score_sigma(window, model, len(window))
        if hit is not None:
            vals.append(float(hit["score"]))
    arr = np.asarray(vals, dtype=float)
    if not arr.size:
        return {"sorted": arr, "mean": 0.0, "sd": 1.0, "n": 0}
    return {
        "sorted": np.sort(arr),
        "mean": float(arr.mean()),
        "sd": float(arr.std(ddof=1)) if arr.size > 1 else 1.0,
        "n": int(arr.size),
    }


def score_windows_with_null(windows: Sequence[str], model: dict[str, Any],
                            decoys_by_length: dict[int, list[str]]
                            ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Score, z y percentil señuelo para cada ventana."""
    n = len(windows)
    scores = np.full(n, -1e9, dtype=float)
    zs = np.zeros(n, dtype=float)
    pcts = np.zeros(n, dtype=float)
    by_len: dict[int, list[int]] = {}
    for i, window in enumerate(windows):
        by_len.setdefault(len(window), []).append(i)
    for length, idx in by_len.items():
        null = _null_distribution(decoys_by_length.get(length, []), model)
        for i in idx:
            hit = score_sigma(windows[i], model, length)
            if hit is None:
                continue
            score = float(hit["score"])
            scores[i] = score
            if null["n"]:
                zs[i] = (score - null["mean"]) / (null["sd"] or 1.0)
                pcts[i] = float(np.searchsorted(null["sorted"], score, side="left")
                                / null["n"])
    return scores, zs, pcts


def auc_score(labels: Sequence[bool], scores: Sequence[float]) -> float | None:
    """AUC ROC por rangos (Mann-Whitney), con corrección de empates."""
    y = np.asarray(labels, dtype=bool)
    s = np.asarray(scores, dtype=float)
    n_pos = int(y.sum())
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    order = np.argsort(s, kind="mergesort")
    sorted_scores = s[order]
    ranks_sorted = np.arange(1, len(s) + 1, dtype=float)
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_scores) != 0) + 1]
    ends = np.r_[starts[1:], len(s)]
    for a, b in zip(starts, ends):
        if b - a > 1:
            ranks_sorted[a:b] = (a + 1 + b) / 2.0
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = ranks_sorted
    return float((ranks[y].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def binary_metrics_from_mask(mask: Sequence[bool],
                             labels: Sequence[bool]) -> dict[str, Any]:
    m = np.asarray(mask, dtype=bool)
    y = np.asarray(labels, dtype=bool)
    tp = int(np.sum(m & y))
    fp = int(np.sum(m & ~y))
    fn = int(np.sum(~m & y))
    n_pred = tp + fp
    n_true = tp + fn
    precision = tp / n_pred if n_pred else None
    recall = tp / n_true if n_true else None
    f1 = (2.0 * precision * recall / (precision + recall)
          if precision and recall else None)
    return {
        "n_pred": n_pred, "tp": tp, "fp": fp, "fn": fn,
        "n_true": n_true,
        "precision": precision, "recall": recall, "f1": f1,
        "base_rate": n_true / len(y) if len(y) else None,
    }


def prepare_arrays(rows: list[dict[str, str]]) -> dict[str, np.ndarray]:
    labels = np.asarray([str(r.get("true", "")).strip().lower() == "yes"
                         for r in rows], dtype=bool)
    is_top1 = np.asarray([str(r.get("is_top1", "")).strip().lower() == "yes"
                          for r in rows], dtype=bool)
    arrays = {"labels": labels, "is_top1": is_top1}
    for feature, column in FEATURE_COLUMNS.items():
        arrays[feature] = np.asarray([_as_float(r.get(column), float("nan"))
                                      for r in rows], dtype=float)
    return arrays


def make_rule(conditions: Iterable[dict[str, Any]]) -> dict[str, Any]:
    conds = [dict(c) for c in conditions]
    label = " AND ".join(
        f"{c['feature']} {c['op']} {c['threshold']:g}" for c in conds)
    return {"conditions": conds, "label": label}


def rule_mask(rule: dict[str, Any], arrays: dict[str, np.ndarray]) -> np.ndarray:
    mask = np.ones(len(arrays["labels"]), dtype=bool)
    for cond in rule["conditions"]:
        values = arrays[cond["feature"]]
        threshold = cond["threshold"]
        if cond["op"] == ">=":
            mask &= values >= threshold
        elif cond["op"] == "<=":
            mask &= values <= threshold
        else:  # pragma: no cover
            raise ValueError(f"operador inválido: {cond['op']}")
    return mask


def evaluate_rule(rule: dict[str, Any], arrays: dict[str, np.ndarray]
                  ) -> dict[str, Any]:
    return binary_metrics_from_mask(rule_mask(rule, arrays), arrays["labels"])


def _with_lift(metrics: dict[str, Any], top1: dict[str, Any]) -> dict[str, Any]:
    out = dict(metrics)
    top1_precision = top1.get("precision")
    top1_f1 = top1.get("f1")
    out["lift_vs_top1"] = (
        metrics["precision"] / top1_precision
        if metrics.get("precision") is not None and top1_precision else None)
    out["f1_delta_vs_top1"] = (
        metrics["f1"] - top1_f1
        if metrics.get("f1") is not None and top1_f1 is not None else None)
    return out


# --------------------------------------------------------------------------
# 1) Core truncation
# --------------------------------------------------------------------------
def _variant_specs() -> list[tuple[str, int | None, str]]:
    specs: list[tuple[str, int | None, str]] = [("full", None, "none")]
    for mode in CORE_MODES:
        for k in K_VALUES:
            specs.append((f"{mode}_k{k}", k, mode))
    return specs


def evaluate_core(rows_by_genome: dict[str, list[dict[str, str]]],
                  seqs_by_genome: dict[str, dict[str, str]],
                  pwms: dict[str, dict[str, Any]],
                  sigmas: Sequence[str], up: int, down: int,
                  decoys: int, seed: int) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for genome, rows_all in rows_by_genome.items():
        seqs = seqs_by_genome[genome]
        default_contig = next(iter(seqs), "")
        for sigma in sigmas:
            rows = [r for r in rows_all if r.get("sigma") == sigma]
            if not rows:
                continue
            key_rows: dict[tuple[str, int, str], dict[str, str]] = {}
            for row in rows:
                key = (row.get("contig") or default_contig,
                       int(_as_float(row.get("tss"))), row.get("strand", ""))
                key_rows[key] = row
            keys = sorted(key_rows)
            windows: list[str] = []
            labels: list[bool] = []
            for key in keys:
                contig, tss, strand = key
                seq = seqs.get(contig)
                if seq is None:
                    raise KeyError(f"{genome}: contig {contig} no está en la secuencia")
                windows.append(_window(seq, tss, strand, up, down))
                labels.append(str(key_rows[key].get("true", "")).strip().lower() == "yes")
            label_arr = np.asarray(labels, dtype=bool)
            lengths = sorted({len(w) for w in windows})
            decoy_windows = {
                length: _decoy_windows(seqs, sigma, seed, decoys, length)
                for length in lengths
            }
            for variant, k, mode in _variant_specs():
                if variant == "full":
                    model = pwms[sigma]
                    kept = {
                        "up": list(range(1, len(pwms[sigma]["up"]["pwm"]) + 1)),
                        "down": list(range(1, len(pwms[sigma]["down"]["pwm"]) + 1)),
                    }
                else:
                    assert k is not None
                    model, kept = build_core_pwm(pwms[sigma], k, mode)
                decoys_by_length = {
                    length: decoy_windows[length] for length in lengths
                }
                scores, zs, pcts = score_windows_with_null(
                    windows, model, decoys_by_length)
                support = (pcts >= SUPPORT_PCT) & (zs >= SUPPORT_Z)
                metrics = binary_metrics_from_mask(support, label_arr)
                row = {
                    "genome": genome,
                    "sigma": sigma,
                    "variant": variant,
                    "mode": mode,
                    "k": k,
                    "n_positions": len(kept["up"]) + len(kept["down"]),
                    "kept_up": ",".join(str(i) for i in kept["up"]),
                    "kept_down": ",".join(str(i) for i in kept["down"]),
                    "auc": auc_score(label_arr, scores),
                    "n_pred": metrics["n_pred"],
                    "tp": metrics["tp"],
                    "fp": metrics["fp"],
                    "fn": metrics["fn"],
                    "precision": metrics["precision"],
                    "recall": metrics["recall"],
                    "f1": metrics["f1"],
                    "base_rate": metrics["base_rate"],
                    "n_true": metrics["n_true"],
                    "n_tss": len(keys),
                }
                out.append(row)
    return out


# --------------------------------------------------------------------------
# 2) Calibración de umbrales
# --------------------------------------------------------------------------
def calibrate_single(arrays: dict[str, np.ndarray], feature: str,
                     objective: str = "f1") -> dict[str, Any] | None:
    """Elige un umbral para una feature usando solo las etiquetas de `arrays`."""
    if feature not in FEATURE_GRID:
        raise ValueError(f"feature inválida: {feature}")
    op = FEATURE_OPS[feature]
    values = arrays[feature]
    candidates = []
    for threshold in FEATURE_GRID[feature]:
        mask = values >= threshold if op == ">=" else values <= threshold
        metrics = binary_metrics_from_mask(mask, arrays["labels"])
        if metrics["n_pred"] < MIN_PRED or metrics["tp"] < MIN_TP:
            continue
        rule = make_rule([{"feature": feature, "op": op,
                           "threshold": float(threshold)}])
        candidates.append({"threshold": float(threshold), "rule": rule,
                           "metrics": metrics})
    if not candidates:
        return None
    # Entre umbrales con las mismas métricas, preferir el más exigente:
    # z/margen más alto, p más bajo. Da margen de decisión y evita elegir
    # un umbral pegado al borde inferior de la grilla.
    strictness = lambda c: (c["threshold"] if op == ">=" else -c["threshold"])
    if objective == "f1":
        pool = [c for c in candidates if c["metrics"]["f1"] is not None]
        if not pool:
            return None
        return max(pool, key=lambda c: (c["metrics"]["f1"],
                                        c["metrics"]["precision"] or -1.0,
                                        -c["metrics"]["n_pred"],
                                        strictness(c)))
    if objective == "precision_at_recall":
        pool = [c for c in candidates
                if (c["metrics"]["recall"] or 0.0) >= MIN_RECALL
                and c["metrics"]["tp"] >= MIN_TP]
        if not pool:
            return None
        return max(pool, key=lambda c: (c["metrics"]["precision"] or -1.0,
                                        c["metrics"]["f1"] or -1.0,
                                        -c["metrics"]["n_pred"],
                                        strictness(c)))
    raise ValueError(f"objetivo inválido: {objective}")


def _calibration_row(sigma: str, feature: str, objective: str,
                     chosen: dict[str, Any], k12_top1: dict[str, Any],
                     sub_arrays: dict[str, np.ndarray],
                     sub_top1: dict[str, Any]) -> dict[str, Any]:
    k12 = chosen["metrics"]
    sub = evaluate_rule(chosen["rule"], sub_arrays)
    return {
        "sigma": sigma,
        "feature": feature,
        "objective": objective,
        "direction": FEATURE_OPS[feature],
        "threshold": chosen["threshold"],
        "rule": chosen["rule"]["label"],
        "k12_n_pred": k12["n_pred"], "k12_tp": k12["tp"],
        "k12_fp": k12["fp"], "k12_precision": k12["precision"],
        "k12_recall": k12["recall"], "k12_f1": k12["f1"],
        "k12_lift_vs_top1": _with_lift(k12, k12_top1)["lift_vs_top1"],
        "k12_base_rate": k12["base_rate"],
        "sub_n_pred": sub["n_pred"], "sub_tp": sub["tp"],
        "sub_fp": sub["fp"], "sub_precision": sub["precision"],
        "sub_recall": sub["recall"], "sub_f1": sub["f1"],
        "sub_lift_vs_top1": _with_lift(sub, sub_top1)["lift_vs_top1"],
        "sub_base_rate": sub["base_rate"],
    }


def calibrate_thresholds(rows_by_genome: dict[str, list[dict[str, str]]],
                         sigmas: Sequence[str],
                         calibration_genome: str = "ecoli_k12",
                         transfer_genome: str = "subtilis"
                         ) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for sigma in sigmas:
        k_rows = [r for r in rows_by_genome[calibration_genome]
                  if r.get("sigma") == sigma]
        s_rows = [r for r in rows_by_genome[transfer_genome]
                  if r.get("sigma") == sigma]
        if not k_rows or not s_rows:
            continue
        k_arr = prepare_arrays(k_rows)
        s_arr = prepare_arrays(s_rows)
        k_top1 = binary_metrics_from_mask(k_arr["is_top1"], k_arr["labels"])
        s_top1 = binary_metrics_from_mask(s_arr["is_top1"], s_arr["labels"])
        for feature in ("z", "p", "margin"):
            for objective in ("f1", "precision_at_recall"):
                chosen = calibrate_single(k_arr, feature, objective)
                if chosen is None:
                    continue
                out.append(_calibration_row(sigma, feature, objective, chosen,
                                            k_top1, s_arr, s_top1))
    return out


# --------------------------------------------------------------------------
# 3) Reglas de dos features
# --------------------------------------------------------------------------
def _condition_masks(arrays: dict[str, np.ndarray]
                     ) -> dict[str, dict[float, np.ndarray]]:
    masks: dict[str, dict[float, np.ndarray]] = {}
    for feature, grid in FEATURE_GRID.items():
        values = arrays[feature]
        op = FEATURE_OPS[feature]
        masks[feature] = {
            threshold: (values >= threshold if op == ">=" else values <= threshold)
            for threshold in grid
        }
    return masks


def two_feature_candidates(arrays: dict[str, np.ndarray],
                           max_per_pair: int = 40) -> list[dict[str, Any]]:
    """Evalúa pares de features con la grilla fija de umbrales."""
    cond_masks = _condition_masks(arrays)
    candidates: list[dict[str, Any]] = []
    features = ("z", "p", "margin", "dist")
    for f1, f2 in combinations(features, 2):
        pair_hits: list[dict[str, Any]] = []
        for t1 in FEATURE_GRID[f1]:
            left = cond_masks[f1][t1]
            for t2 in FEATURE_GRID[f2]:
                mask = left & cond_masks[f2][t2]
                metrics = binary_metrics_from_mask(mask, arrays["labels"])
                if metrics["n_pred"] < MIN_PRED or metrics["tp"] < MIN_TP:
                    continue
                pair_hits.append({
                    "rule": make_rule([
                        {"feature": f1, "op": FEATURE_OPS[f1],
                         "threshold": float(t1)},
                        {"feature": f2, "op": FEATURE_OPS[f2],
                         "threshold": float(t2)},
                    ]),
                    "features": [f1, f2],
                    "metrics": metrics,
                })
        pair_hits.sort(key=lambda c: (c["metrics"]["f1"] or -1.0,
                                      c["metrics"]["precision"] or -1.0),
                       reverse=True)
        candidates.extend(pair_hits[:max_per_pair])
    return candidates


def select_two_feature(candidates: list[dict[str, Any]],
                       objective: str = "f1") -> dict[str, Any] | None:
    if objective == "f1":
        pool = [c for c in candidates if c["metrics"]["f1"] is not None]
        if not pool:
            return None
        return max(pool, key=lambda c: (c["metrics"]["f1"],
                                        c["metrics"]["precision"] or -1.0,
                                        -c["metrics"]["n_pred"],
                                        c["rule"]["label"]))
    if objective == "precision_at_recall":
        pool = [c for c in candidates
                if (c["metrics"]["recall"] or 0.0) >= MIN_RECALL
                and c["metrics"]["tp"] >= MIN_TP]
        if not pool:
            return None
        return max(pool, key=lambda c: (c["metrics"]["precision"] or -1.0,
                                        c["metrics"]["f1"] or -1.0,
                                        -c["metrics"]["n_pred"],
                                        c["rule"]["label"]))
    raise ValueError(f"objetivo inválido: {objective}")


def _two_feature_row(sigma: str, objective: str, selected: bool,
                     rule: dict[str, Any], k12_metrics: dict[str, Any],
                     k12_top1: dict[str, Any],
                     sub_metrics: dict[str, Any],
                     sub_top1: dict[str, Any]) -> dict[str, Any]:
    return {
        "sigma": sigma,
        "objective": objective,
        "selected": bool(selected),
        "rule": rule["label"],
        "k12_n_pred": k12_metrics["n_pred"], "k12_tp": k12_metrics["tp"],
        "k12_fp": k12_metrics["fp"], "k12_precision": k12_metrics["precision"],
        "k12_recall": k12_metrics["recall"], "k12_f1": k12_metrics["f1"],
        "k12_lift_vs_top1": _with_lift(k12_metrics, k12_top1)["lift_vs_top1"],
        "sub_n_pred": sub_metrics["n_pred"], "sub_tp": sub_metrics["tp"],
        "sub_fp": sub_metrics["fp"], "sub_precision": sub_metrics["precision"],
        "sub_recall": sub_metrics["recall"], "sub_f1": sub_metrics["f1"],
        "sub_lift_vs_top1": _with_lift(sub_metrics, sub_top1)["lift_vs_top1"],
    }


def evaluate_two_feature_rules(rows_by_genome: dict[str, list[dict[str, str]]],
                               sigmas: Sequence[str],
                               calibration_genome: str = "ecoli_k12",
                               transfer_genome: str = "subtilis"
                               ) -> tuple[list[dict[str, Any]],
                                          list[dict[str, Any]]]:
    out: list[dict[str, Any]] = []
    top: list[dict[str, Any]] = []
    for sigma in sigmas:
        k_rows = [r for r in rows_by_genome[calibration_genome]
                  if r.get("sigma") == sigma]
        s_rows = [r for r in rows_by_genome[transfer_genome]
                  if r.get("sigma") == sigma]
        if not k_rows or not s_rows:
            continue
        k_arr = prepare_arrays(k_rows)
        s_arr = prepare_arrays(s_rows)
        k_top1 = binary_metrics_from_mask(k_arr["is_top1"], k_arr["labels"])
        s_top1 = binary_metrics_from_mask(s_arr["is_top1"], s_arr["labels"])
        candidates = two_feature_candidates(k_arr)
        for objective in ("f1", "precision_at_recall"):
            chosen = select_two_feature(candidates, objective)
            if chosen is None:
                continue
            sub_metrics = evaluate_rule(chosen["rule"], s_arr)
            out.append(_two_feature_row(
                sigma, objective, True, chosen["rule"], chosen["metrics"],
                k_top1, sub_metrics, s_top1))
        # Referencias: regla de dos features ya usada y regla de tres
        # features que menciona el contexto del proyecto.
        references = [
            ("reference_2f", make_rule([
                {"feature": "margin", "op": ">=", "threshold": 1.0},
                {"feature": "dist", "op": "<=", "threshold": 60.0},
            ])),
            ("reference_3f", make_rule([
                {"feature": "p", "op": "<=", "threshold": 0.05},
                {"feature": "margin", "op": ">=", "threshold": 1.0},
                {"feature": "dist", "op": "<=", "threshold": 60.0},
            ])),
        ]
        for label, rule in references:
            k_metrics = evaluate_rule(rule, k_arr)
            sub_metrics = evaluate_rule(rule, s_arr)
            out.append(_two_feature_row(
                sigma, label, False, rule, k_metrics, k_top1,
                sub_metrics, s_top1))
        candidates.sort(key=lambda c: (c["metrics"]["f1"] or -1.0,
                                       c["metrics"]["precision"] or -1.0),
                        reverse=True)
        for cand in candidates[:10]:
            top.append({
                "sigma": sigma,
                "rule": cand["rule"]["label"],
                "features": cand["features"],
                "k12": cand["metrics"],
            })
    return out, top


# --------------------------------------------------------------------------
# Salidas
# --------------------------------------------------------------------------
def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_jsonable(payload), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")


def _fmt_cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_tsv(path: Path, rows: Iterable[dict[str, Any]],
              columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(_fmt_cell(row.get(col)) for col in columns) + "\n")


def build_summary(core_rows: list[dict[str, Any]],
                  cal_rows: list[dict[str, Any]],
                  two_rows: list[dict[str, Any]],
                  baselines: dict[str, dict[str, Any]],
                  genomes: Sequence[str], sigmas: Sequence[str]) -> dict[str, Any]:
    core_findings: dict[str, Any] = {}
    for genome in genomes:
        for sigma in sigmas:
            sub = [r for r in core_rows
                   if r["genome"] == genome and r["sigma"] == sigma]
            if not sub:
                continue
            full = next((r for r in sub if r["variant"] == "full"), None)
            full_positions = full["n_positions"] if full else 0
            # "core" de verdad: variantes que descartan al menos una columna.
            cores = [r for r in sub
                     if r["variant"] != "full" and r["n_positions"] < full_positions]
            if not cores:
                cores = [r for r in sub if r["variant"] != "full"]
            best = max(cores, key=lambda r: (r["auc"] or -1.0)) if cores else None
            core_findings[f"{genome}/{sigma}"] = {
                "full_auc": full and full["auc"],
                "full_support_recall": full and full["recall"],
                "full_support_precision": full and full["precision"],
                "full_n_positions": full_positions,
                "best_core_variant_by_auc": best and best["variant"],
                "best_core_n_positions": best and best["n_positions"],
                "best_core_auc": best and best["auc"],
                "best_core_auc_delta": (
                    best["auc"] - full["auc"]
                    if best and full and best["auc"] is not None
                    and full["auc"] is not None else None),
                "best_core_support_recall": best and best["recall"],
            }
    # Transferencia de cada variante: mismo recorte en K-12 y en B. subtilis.
    by_key = {(r["genome"], r["sigma"], r["variant"]): r for r in core_rows}
    core_transfer: dict[str, Any] = {}
    for sigma in sigmas:
        full_k = by_key.get(("ecoli_k12", sigma, "full"))
        full_s = by_key.get(("subtilis", sigma, "full"))
        variants = sorted({r["variant"] for r in core_rows if r["sigma"] == sigma})
        core_transfer[sigma] = {}
        for variant in variants:
            k_row = by_key.get(("ecoli_k12", sigma, variant))
            s_row = by_key.get(("subtilis", sigma, variant))
            if not k_row or not s_row:
                continue
            core_transfer[sigma][variant] = {
                "n_positions": k_row["n_positions"],
                "k12_auc": k_row["auc"],
                "sub_auc": s_row["auc"],
                "k12_auc_delta": (
                    k_row["auc"] - full_k["auc"]
                    if full_k and k_row["auc"] is not None
                    and full_k["auc"] is not None else None),
                "sub_auc_delta": (
                    s_row["auc"] - full_s["auc"]
                    if full_s and s_row["auc"] is not None
                    and full_s["auc"] is not None else None),
                "k12_support_recall": k_row["recall"],
                "sub_support_recall": s_row["recall"],
            }
    cal_findings: dict[str, Any] = {}
    for sigma in sigmas:
        rows = [r for r in cal_rows if r["sigma"] == sigma
                and r["objective"] == "f1"]
        cal_findings[sigma] = {
            r["feature"]: {
                "threshold": r["threshold"],
                "k12_f1": r["k12_f1"],
                "sub_f1": r["sub_f1"],
                "sub_lift_vs_top1": r["sub_lift_vs_top1"],
            } for r in rows
        }
    two_findings: dict[str, Any] = {}
    for sigma in sigmas:
        rows = [r for r in two_rows if r["sigma"] == sigma
                and r["objective"] == "f1"]
        if rows:
            best = rows[0]
            two_findings[sigma] = {
                "rule": best["rule"],
                "k12_precision": best["k12_precision"],
                "k12_recall": best["k12_recall"],
                "k12_f1": best["k12_f1"],
                "k12_lift_vs_top1": best["k12_lift_vs_top1"],
                "sub_precision": best["sub_precision"],
                "sub_recall": best["sub_recall"],
                "sub_f1": best["sub_f1"],
                "sub_lift_vs_top1": best["sub_lift_vs_top1"],
            }
    return {
        "version": VERSION,
        "baselines_top1": baselines,
        "core": core_findings,
        "core_transfer": core_transfer,
        "calibration": cal_findings,
        "two_feature": two_findings,
        "files": [
            "core_truncation.tsv", "core_truncation.json",
            "threshold_calibration.tsv", "threshold_calibration.json",
            "two_feature_rules.tsv", "two_feature_rules.json",
            "summary.json",
        ],
    }


# --------------------------------------------------------------------------
# Selftest y CLI
# --------------------------------------------------------------------------
def selftest() -> str:
    pwm = np.asarray([
        [1.5, -5.0, -5.0, -5.0],
        [-5.0, 1.5, -5.0, -5.0],
        [-5.0, -5.0, 1.5, -5.0],
        [-5.0, -5.0, -5.0, 1.5],
    ], dtype=float)
    rec = {
        "sigma": "sigma54",
        "class": "alternative",
        "arch": "-24/-12",
        "up": {"width": 4, "pwm": pwm.copy()},
        "down": {"width": 4, "pwm": pwm.copy()},
    }
    core, kept = build_core_pwm(rec, 2, "joint")
    assert core["up"]["width"] == 4 and core["down"]["width"] == 4
    assert len(kept["up"]) + len(kept["down"]) == 2
    assert int(np.sum(~np.any(core["up"]["pwm"], axis=1))) == 4 - len(kept["up"])
    assert int(np.sum(~np.any(core["down"]["pwm"], axis=1))) == 4 - len(kept["down"])
    assert auc_score([True, False, True, False], [3, 1, 2, 0]) == 1.0
    assert auc_score([False, True], [3, 1]) == 0.0
    m = binary_metrics_from_mask([True, False, True, False],
                                 [True, True, False, False])
    assert (m["tp"], m["fp"], m["fn"]) == (1, 1, 1)
    rows = [
        {"true": "yes", "is_top1": "yes", "z": str(z), "pvalue": "0.001",
         "margin": "2", "dist_to_start": "30"}
        for z in (4, 3.5, 3, 2.5, 2)
    ] + [
        {"true": "no", "is_top1": "no", "z": str(z), "pvalue": "0.5",
         "margin": "-1", "dist_to_start": "70"}
        for z in (1, 0.5, 0, -0.5, -1)
    ]
    arrays = prepare_arrays(rows)
    chosen = calibrate_single(arrays, "z", "f1")
    assert chosen is not None and chosen["threshold"] >= 2.0
    assert "n_true" not in FEATURE_GRID and "gold_confidence" not in FEATURE_GRID
    return "calibración de motivo sigma OK"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", nargs="+", default=list(DEFAULT_GENOMES))
    ap.add_argument("--sigmas", nargs="+", default=list(DEFAULT_SIGMAS))
    ap.add_argument("--up", type=int, default=DEFAULT_UP)
    ap.add_argument("--down", type=int, default=DEFAULT_DOWN)
    ap.add_argument("--decoys", type=int, default=DEFAULT_DECOYS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "sigma" / "motif_calibration")
    ap.add_argument("--limit", type=int, default=0,
                    help="Limita los TSS por genoma (smoke test; no usar en serio).")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0
    if args.decoys < 20:
        ap.error("--decoys debe ser >= 20")
    data_root = ROOT / "data"
    pwms = load_pwm_library(data_root)
    missing = [s for s in args.sigmas if s not in pwms]
    if missing:
        print(f"Sin PWM para: {missing}", file=sys.stderr)
        return 2
    rows_by_genome: dict[str, list[dict[str, str]]] = {}
    seqs_by_genome: dict[str, dict[str, str]] = {}
    for genome in args.genomes:
        path = data_root / "sigma" / "layers" / genome / "layer_sigma_features.tsv"
        rows = _limit_rows(read_features(path), args.limit)
        rows_by_genome[genome] = rows
        seqs = load_genome_sequences(data_root, genome)
        if not seqs:
            print(f"Sin secuencia para {genome}", file=sys.stderr)
            return 2
        seqs_by_genome[genome] = seqs
    baselines: dict[str, dict[str, Any]] = {}
    for genome in args.genomes:
        baselines[genome] = {}
        for sigma in args.sigmas:
            rows = [r for r in rows_by_genome[genome] if r.get("sigma") == sigma]
            if rows:
                arrays = prepare_arrays(rows)
                baselines[genome][sigma] = binary_metrics_from_mask(
                    arrays["is_top1"], arrays["labels"])
    core_rows = evaluate_core(rows_by_genome, seqs_by_genome, pwms,
                              args.sigmas, args.up, args.down,
                              args.decoys, args.seed)
    cal_rows = calibrate_thresholds(rows_by_genome, args.sigmas)
    two_rows, top_two = evaluate_two_feature_rules(rows_by_genome, args.sigmas)

    out = args.out
    config = {
        "version": VERSION,
        "genomes": list(args.genomes),
        "sigmas": list(args.sigmas),
        "window": [-args.up, args.down],
        "decoys": args.decoys,
        "seed": args.seed,
        "support_criterion": {"pct_decoy": SUPPORT_PCT, "z": SUPPORT_Z},
        "min_pred": MIN_PRED,
        "min_tp": MIN_TP,
        "forbidden_columns": list(FORBIDDEN_COLUMNS),
        "limit": args.limit,
    }
    write_tsv(out / "core_truncation.tsv", core_rows, CORE_TSV_COLS)
    write_json(out / "core_truncation.json",
               {"config": config, "rows": core_rows})
    write_tsv(out / "threshold_calibration.tsv", cal_rows, CAL_TSV_COLS)
    write_json(out / "threshold_calibration.json",
               {"config": config, "rows": cal_rows})
    write_tsv(out / "two_feature_rules.tsv", two_rows, TWO_TSV_COLS)
    write_json(out / "two_feature_rules.json",
               {"config": config, "selected": two_rows, "top_k12": top_two})
    summary = build_summary(core_rows, cal_rows, two_rows, baselines,
                            args.genomes, args.sigmas)
    summary["config"] = config
    write_json(out / "summary.json", summary)
    print(f"salida: {out}")
    print(f"core: {len(core_rows)} filas · "
          f"calibración: {len(cal_rows)} filas · "
          f"reglas: {len(two_rows)} filas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
