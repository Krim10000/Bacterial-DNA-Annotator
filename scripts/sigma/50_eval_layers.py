#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evaluación por capas del pipeline sigma, con soporte para gold curado.

El objetivo no es un único accuracy: es medir dónde falla cada capa y separar
errores reales de procesos biológicos superpuestos.

Capas:
  0. Datos: TSS, secuencia, PWM, etiqueta.
  1. Localización: ¿el TSS correcto tiene mejor señal que ventanas señuelo?
  2. Motivo: ¿el PWM del sigma verdadero está soportado en ese TSS?
  3. Asignación: top-1/top-3, multietiqueta, matriz de confusión.
  4. Contexto: disponibilidad de inventario/ruta/expresión y su relación con
     el acierto (no se inventa contexto si no existe).
  5. Integración: enlace al regulón, si está disponible.

Gold esperado:
  data/sigma/gold/<id>/tss_sigma.tsv

Columnas aceptadas (sin orden fijo):
  contig, tss/pos/position, strand, sigma/sigma_factor/factor,
  operon/operon_id, gene/first_gene, evidence, confidence.

Una fila por TSS y sigma; si hay varios sigma para el mismo TSS, se agrupan.
Un TSS con varios sigma NO es un error: se evalúa como multietiqueta
(`hit_any`, `hit_all`, Jaccard), no como top-1 único.

Uso:
    python scripts/sigma/50_eval_layers.py --genome ecoli_k12 --source k12
    python scripts/sigma/50_eval_layers.py --genome subtilis --source subtilis
    python scripts/sigma/50_eval_layers.py --gold-root data/sigma/gold --all
    python scripts/sigma/50_eval_layers.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import statistics
import sys
import time
import zlib
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma.sigma_annotator import (  # noqa: E402
    SIGMA_ORDER, SIGMA_SPACER, _pvalue, _read_tsv, _zscore, build_null,
    load_pwm_library, load_sequence, revcomp, score_sigma,
)

VERSION = "1.0.0"
DEFAULT_WINDOW_UP = 80
DEFAULT_WINDOW_DOWN = 20
DEFAULT_DECOYS = 200
DEFAULT_SEED = 17
PWMS_PATH = None
DEFAULT_Z = 2.0
DEFAULT_P = 0.05

# Etiquetas de B. subtilis con equivalente funcional defendible.
BSU_MAP = {
    "A": "sigma70", "SIGA": "sigma70",
    "B": "sigma38", "SIGB": "sigma38",
    "D": "sigma28", "SIGD": "sigma28",
    "I": "sigma32", "SIGI": "sigma32",
    "L": "sigma54", "SIGL": "sigma54",
}

COL_ALIASES = {
    "contig": ("contig", "seqid", "chromosome", "replicon"),
    "tss": ("tss", "pos", "position", "start", "coordinate"),
    "strand": ("strand", "hebra"),
    "sigma": ("sigma", "sigma_factor", "factor", "sigmas", "sigma_type"),
    "operon": ("operon", "operon_id", "target_operon", "tu"),
    "gene": ("gene", "first_gene", "gene_name", "target_gene"),
    "evidence": ("evidence", "source", "evidencia"),
    "confidence": ("confidence", "confianza"),
}


def _pick(row: dict[str, str], key: str) -> str:
    for alias in COL_ALIASES.get(key, (key,)):
        value = row.get(alias)
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def _norm_strand(value: str) -> str:
    v = str(value or "").strip().lower()
    if v in {"+", "1", "forward", "fwd", "plus"}:
        return "+"
    if v in {"-", "-1", "reverse", "rev", "minus"}:
        return "-"
    return ""


def _split_sigmas(value: str) -> list[str]:
    raw = str(value or "").strip()
    if not raw or raw in {"?", "NA", "na", "unlabeled", "none", "-"}:
        return []
    parts = [p.strip() for p in re.split(r"[;,|/]+", raw) if p.strip()]
    return parts


def map_sigma(raw: str, source: str = "auto") -> tuple[str, str]:
    """Devuelve (sigma canónico, status). status: ok | no_pwm | unlabeled."""
    s = str(raw or "").strip()
    if not s or s in {"?", "NA", "na", "unlabeled", "none", "-"}:
        return "", "unlabeled"
    up = s.upper()
    if source == "subtilis" and up in BSU_MAP:
        return BSU_MAP[up], "ok"
    if source != "subtilis" and up in BSU_MAP:
        return BSU_MAP[up], "ok"
    for cand in (s, s.lower(), s.replace("σ", "sigma")):
        low = cand.lower()
        if low in SIGMA_ORDER:
            return low, "ok"
        if low.startswith("sigma") and low[5:] in {"70", "54", "32", "38", "28", "24"}:
            return "sigma" + low[5:], "ok"
        if low in {"sigmaa", "sigma a"}:
            return "sigma70", "ok"
    # B. subtilis letters only when source says so; no adivinar en K-12.
    if source == "subtilis":
        return "", "no_pwm"
    return "", "no_pwm"


def load_gold(path: Path, source: str = "auto") -> list[dict[str, Any]]:
    """Carga un TSV gold y agrupa por (contig, tss, strand)."""
    rows = _read_tsv(path)
    grouped: dict[tuple[str, int, str], dict[str, Any]] = {}
    for row in rows:
        contig = _pick(row, "contig")
        tss_raw = _pick(row, "tss")
        strand = _norm_strand(_pick(row, "strand"))
        try:
            tss = int(float(tss_raw))
        except (TypeError, ValueError):
            continue
        key = (contig, tss, strand)
        rec = grouped.setdefault(key, {
            "contig": contig, "tss": tss, "strand": strand,
            "sigmas": set(), "sigmas_raw": set(), "operons": set(),
            "genes": set(), "evidence": set(), "confidence": set(),
            "n_rows": 0,
        })
        rec["n_rows"] += 1
        for raw in _split_sigmas(_pick(row, "sigma")):
            rec["sigmas_raw"].add(raw)
            mapped, _ = map_sigma(raw, source)
            if mapped:
                rec["sigmas"].add(mapped)
        for field in ("operon", "gene"):
            value = _pick(row, field)
            if value:
                rec[f"{field}s"].add(value)
        ev = _pick(row, "evidence")
        if ev:
            rec["evidence"].add(ev)
        cf = _pick(row, "confidence")
        if cf:
            rec["confidence"].add(cf)
    out = []
    for rec in grouped.values():
        rec["sigmas"] = sorted(rec["sigmas"])
        rec["sigmas_raw"] = sorted(rec["sigmas_raw"])
        rec["operons"] = sorted(rec["operons"])
        rec["genes"] = sorted(rec["genes"])
        rec["evidence"] = sorted(rec["evidence"])
        rec["confidence"] = sorted(rec["confidence"])
        out.append(rec)
    out.sort(key=lambda r: (r["contig"], r["tss"], r["strand"]))
    return out


def load_k12(data_root: Path) -> list[dict[str, Any]]:
    path = data_root / "regulondb" / "gt" / "gt_promoters.tsv"
    rows = _read_tsv(path)
    out = []
    for row in rows:
        try:
            tss = int(float(row.get("tss") or ""))
        except (TypeError, ValueError):
            continue
        raw = row.get("sigma_factor") or ""
        sigma, status = map_sigma(raw, "k12")
        out.append({
            "contig": "", "tss": tss,
            "strand": _norm_strand(row.get("strand") or ""),
            "sigmas": [sigma] if sigma else [],
            "sigmas_raw": [raw] if raw else [],
            "operons": [], "genes": [row.get("first_gene") or ""],
            "evidence": [row.get("name") or ""],
            "confidence": [row.get("confidence") or ""],
            "n_rows": 1, "map_status": status,
        })
    return out


def load_subtilis(data_root: Path) -> list[dict[str, Any]]:
    path = data_root / "dbtbs" / "gt_tss.tsv"
    rows = _read_tsv(path)
    out = []
    for row in rows:
        try:
            tss = int(float(row.get("pos") or ""))
        except (TypeError, ValueError):
            continue
        raw = row.get("sigma") or ""
        sigma, status = map_sigma(raw, "subtilis")
        out.append({
            "contig": "", "tss": tss,
            "strand": _norm_strand(row.get("strand") or ""),
            "sigmas": [sigma] if sigma else [],
            "sigmas_raw": [raw] if raw else [],
            "operons": [], "genes": [],
            "evidence": [row.get("evidence") or ""],
            "confidence": [],
            "n_rows": 1, "map_status": status,
        })
    return out


def _window(seq: str, tss: int, strand: str, up: int, down: int) -> str:
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


def _decoy_scores(seqs: dict[str, str], sigma: str, rec: dict, length: int,
                  n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    contigs = [s for s in seqs.values() if len(s) >= length]
    if not contigs:
        return np.array([])
    vals = []
    for _ in range(n):
        s = contigs[int(rng.integers(0, len(contigs)))]
        a = int(rng.integers(0, len(s) - length + 1))
        w = s[a:a + length]
        h = score_sigma(w, rec, len(w))
        if h is not None:
            vals.append(h["score"])
    return np.asarray(vals, dtype=float)


def _decoy_stats(seqs: dict[str, str], sigma: str, rec: dict, length: int,
                 decoys: int, seed: int,
                 cache: dict[tuple[str, int], dict[str, Any]]) -> dict[str, Any]:
    key = (sigma, length)
    if key in cache:
        return cache[key]
    dec = _decoy_scores(seqs, sigma, rec, length, decoys, seed)
    if dec.size:
        cache[key] = {
            "sorted": np.sort(dec),
            "mean": float(dec.mean()),
            "sd": float(dec.std(ddof=1)) if dec.size > 1 else 1.0,
        }
    else:
        cache[key] = {"sorted": np.array([]), "mean": 0.0, "sd": 1.0}
    return cache[key]


def _score_tss(seqs: dict[str, str], gold: dict[str, Any],
               pwms: dict[str, dict], up: int, down: int,
               decoys: int, seed: int,
               decoy_cache: dict[tuple[str, int], dict[str, Any]]) -> dict[str, Any]:
    contig = gold["contig"] or next(iter(seqs), "")
    seq = seqs.get(contig, "")
    w = _window(seq, gold["tss"], gold["strand"], up, down)
    out = {
        "sequence_available": bool(w),
        "window": w,
        "window_len": len(w),
        "scores": {}, "true_motif": {},
    }
    if not w:
        return out
    for sigma, rec in pwms.items():
        h = score_sigma(w, rec, len(w))
        if h is None:
            continue
        stats = _decoy_stats(seqs, sigma, rec, len(w), decoys,
                             seed + (zlib.crc32(sigma.encode()) % 1000), decoy_cache)
        dec = stats["sorted"]
        if dec.size:
            z = float((h["score"] - stats["mean"]) / (stats["sd"] or 1.0))
            pct = float((dec < h["score"]).mean())
            pvalue = float((1 + int(np.sum(dec >= h["score"]))) / (1 + dec.size))
        else:
            z, pct, pvalue = 0.0, 0.0, 1.0
        h["z"] = round(z, 3)
        h["pvalue"] = round(pvalue, 5)
        h["pct_decoy"] = round(pct, 4)
        out["scores"][sigma] = h
        if sigma in gold["sigmas"]:
            out["true_motif"][sigma] = {
                "score": round(float(h["score"]), 3),
                "pct_decoy": round(pct, 4),
                "z_decoy": round(z, 3),
                "supported": bool(pct >= 0.95 and z >= 1.0),
            }
    return out


def _topk(scores: dict[str, dict], k: int = 3) -> list[str]:
    ranked = sorted(scores.items(), key=lambda kv: (kv[1].get("z", 0),
                                                    -kv[1].get("pvalue", 1),
                                                    kv[1].get("score", 0)),
                    reverse=True)
    return [s for s, _ in ranked[:k]]


def _per_sigma_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    labels = sorted({s for r in rows for s in r["true_sigmas"]})
    out = {}
    for sigma in labels:
        tp = fp = fn = 0
        for r in rows:
            true = set(r["true_sigmas"])
            # Single-label: se evalúa top-1. Multi-label: se permite top-3,
            # porque un solo top-1 no puede cubrir varios sigma verdaderos.
            pred = ({r["pred_top1"]} if r.get("n_true", 1) == 1
                    else set(r["pred_top3"]))
            if sigma in true and sigma in pred:
                tp += 1
            elif sigma in pred:
                fp += 1
            elif sigma in true:
                fn += 1
        p = tp / (tp + fp) if tp + fp else None
        rec = tp / (tp + fn) if tp + fn else None
        f1 = (2 * p * rec / (p + rec)) if p and rec else None
        out[sigma] = {"support": sum(1 for r in rows if sigma in r["true_sigmas"]),
                      "tp": tp, "fp": fp, "fn": fn,
                      "precision": p, "recall": rec, "f1": f1}
    return out


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = None
    buf: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        lines = fh.readlines()
    for line in lines:
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


def _fallback_sequence(data_root: Path, genome: str) -> dict[str, str]:
    candidates = [
        data_root / "reference_genomes" / genome / f"{genome}_genomic.fna.gz",
        data_root / "reference_genomes" / genome / f"{genome}_genomic.fna",
        data_root / "raw" / f"{genome}_genomic.fna.gz",
    ]
    if genome.lower() == "subtilis":
        candidates.append(data_root / "dbtbs" / "genome.fna")
    for path in candidates:
        if path.exists():
            got = _read_fasta(path)
            if got:
                return got
    return {}


def evaluate(rows: list[dict[str, Any]], data_root: Path, genome: str,
             source: str, up: int, down: int, decoys: int, seed: int,
             zmin: float, pmax: float) -> dict[str, Any]:
    seqs = load_sequence(data_root, genome)
    if not seqs:
        seqs = _fallback_sequence(data_root, genome)
    pwms = load_pwm_library(data_root, PWMS_PATH)
    if not seqs:
        raise ValueError(f"sin secuencia para {genome}")
    scored = []
    decoy_cache: dict[tuple[str, int], dict[str, Any]] = {}
    for gold in rows:
        sc = _score_tss(seqs, gold, pwms, up, down, decoys, seed, decoy_cache)
        # ranking global: usar z del score contra decoys por sigma
        ranking = []
        for sigma, h in sc["scores"].items():
            ranking.append((float(h.get("z", 0.0)),
                            -float(h.get("pvalue", 1.0)), sigma))
        ranking.sort(reverse=True)
        pred = [s for _, _, s in ranking]
        true = list(gold["sigmas"])
        rec = {
            **gold,
            "true_sigmas": true,
            "pred_top1": pred[0] if pred else "",
            "pred_top3": pred[:3],
            "pred_all": pred,
            "site": sc,
            "n_true": len(true),
            "has_pwm": bool(true) and all(s in pwms for s in true),
            "no_pwm": bool(true) and not all(s in pwms for s in true),
        }
        if not true and gold.get("map_status") == "no_pwm":
            rec["class"] = "no_pwm"
        elif not true:
            rec["class"] = "unlabeled"
        elif not sc["sequence_available"]:
            rec["class"] = "missing_sequence"
        elif not rec["has_pwm"]:
            rec["class"] = "no_pwm"
        elif not pred:
            rec["class"] = "no_prediction"
        elif rec["pred_top1"] in true:
            rec["class"] = "correct" if len(true) == 1 else "hit_one_of_multi"
        elif set(pred[:3]) & set(true):
            rec["class"] = "partial_top3"
        else:
            rec["class"] = "missed"
        rec["jaccard_top3"] = (
            len(set(pred[:3]) & set(true)) / len(set(pred[:3]) | set(true))
            if pred or true else None)
        scored.append(rec)

    # Tabla de features por (TSS, sigma) para minería de combinaciones.
    inventory_path = data_root / "sigma" / genome / "sigma_factors.tsv"
    inventory = {str(r.get("sigma_type") or "") for r in _read_tsv(inventory_path)}
    feature_rows: list[dict[str, Any]] = []
    for rec in scored:
        scores = rec.get("site", {}).get("scores", {})
        if not scores:
            continue
        zs = sorted((float(h.get("z", 0.0)), s) for s, h in scores.items())
        top_z = zs[-1][0] if zs else 0.0
        true_set = set(rec["true_sigmas"])
        for sigma, h in scores.items():
            others = [z for z, s in zs if s != sigma]
            margin = float(h.get("z", 0.0)) - (max(others) if others else 0.0)
            lo, hi = SIGMA_SPACER.get(sigma, (14, 19))
            spacer = int(h.get("spacer") or 0)
            spacer_dev = 0 if lo <= spacer <= hi else min(abs(spacer - lo), abs(spacer - hi))
            feature_rows.append({
                "genome": genome,
                "contig": rec.get("contig", ""),
                "tss": rec.get("tss", ""),
                "strand": rec.get("strand", ""),
                "sigma": sigma,
                "true": "yes" if sigma in true_set else "no",
                "is_top1": "yes" if sigma == rec.get("pred_top1") else "no",
                "z": round(float(h.get("z", 0.0)), 3),
                "pvalue": round(float(h.get("pvalue", 1.0)), 5),
                "score": round(float(h.get("score", 0.0)), 3),
                "spacer": spacer,
                "spacer_dev": spacer_dev,
                "dist_to_start": h.get("dist_to_start", ""),
                "margin": round(margin, 3),
                "top_z": round(top_z, 3),
                "inventory": "yes" if sigma in inventory else "no",
                "n_true": rec.get("n_true", 0),
                "gold_confidence": ",".join(rec.get("confidence", [])),
                "source": rec.get("source", ""),
            })

    evaluable = [r for r in scored if r["class"] not in
                 {"unlabeled", "missing_sequence", "no_pwm", "no_prediction"}]
    single = [r for r in evaluable if r["n_true"] == 1]
    multi = [r for r in evaluable if r["n_true"] > 1]
    labels = sorted({s for r in evaluable for s in r["true_sigmas"]})
    majority = Counter(s for r in evaluable for s in r["true_sigmas"]).most_common(1)
    n = len(evaluable)
    top1 = sum(1 for r in evaluable if r["pred_top1"] in r["true_sigmas"])
    top3_any = sum(1 for r in evaluable if set(r["pred_top3"]) & set(r["true_sigmas"]))
    top3_all = sum(1 for r in evaluable if set(r["true_sigmas"]) <= set(r["pred_top3"]))
    correct_single = sum(1 for r in single if r["pred_top1"] in r["true_sigmas"])
    correct_multi_any = sum(1 for r in multi if r["pred_top1"] in r["true_sigmas"])
    exact_multi = sum(1 for r in multi if set(r["true_sigmas"]) == {r["pred_top1"]})

    # Capa 1: localización
    site_rows = [r for r in scored if r["site"]["sequence_available"] and r["true_sigmas"]]
    site_supported = 0
    for r in site_rows:
        vals = [v for v in r["site"]["true_motif"].values() if v.get("supported")]
        if vals:
            site_supported += 1
    # Capa 2: motivo por sigma
    motif_by_sigma = {}
    for sigma in labels:
        vals = [r["site"]["true_motif"].get(sigma) for r in scored
                if sigma in r["true_sigmas"]]
        vals = [v for v in vals if v]
        motif_by_sigma[sigma] = {
            "n": len(vals),
            "supported": sum(1 for v in vals if v.get("supported")),
            "pct_supported": (sum(1 for v in vals if v.get("supported")) / len(vals)
                              if vals else None),
            "median_pct_decoy": statistics.median([v["pct_decoy"] for v in vals])
            if vals else None,
        }

    # Capa 4: contexto disponible (si existe sigma_promoters.tsv)
    context_path = data_root / "sigma" / genome / "sigma_promoters.tsv"
    context_rows = _read_tsv(context_path) if context_path.exists() else []
    context_genes = {(r.get("first_gene_name") or "").lower() for r in context_rows
                     if r.get("first_gene_name")}
    context_hits = 0
    for r in scored:
        genes = {g.lower() for g in r.get("genes", []) if g}
        if genes & context_genes:
            context_hits += 1

    result = {
        "version": VERSION,
        "genome": genome,
        "source": source,
        "config": {"window": [-up, down], "decoys": decoys, "seed": seed,
                   "zmin": zmin, "pmax": pmax},
        "layer0_data": {
            "n_gold_tss": len(rows),
            "n_with_sequence": sum(1 for r in scored if r["site"]["sequence_available"]),
            "n_unlabeled": sum(1 for r in scored if r["class"] == "unlabeled"),
            "n_no_pwm": sum(1 for r in scored if r["class"] == "no_pwm"),
            "n_missing_sequence": sum(1 for r in scored if r["class"] == "missing_sequence"),
            "coverage_sequence": (sum(1 for r in scored if r["site"]["sequence_available"])
                                  / len(rows) if rows else None),
        },
        "layer1_localization": {
            "n": len(site_rows),
            "n_supported": site_supported,
            "frac_supported": site_supported / len(site_rows) if site_rows else None,
            "note": "TSS con al menos un sigma verdadero cuyo score supera el 95% de señuelos",
        },
        "layer2_motif": motif_by_sigma,
        "layer3_assignment": {
            "n_evaluable": n,
            "n_single": len(single),
            "n_multi": len(multi),
            "top1_accuracy": top1 / n if n else None,
            "top3_any": top3_any / n if n else None,
            "top3_all": top3_all / n if n else None,
            "single_top1": correct_single / len(single) if single else None,
            "multi_top1_any": correct_multi_any / len(multi) if multi else None,
            "multi_exact_top1": exact_multi / len(multi) if multi else None,
            "majority_class": majority[0][0] if majority else "",
            "majority_baseline": majority[0][1] / n if n and majority else None,
            "per_sigma": _per_sigma_metrics(evaluable),
        },
        "layer4_context": {
            "sigma_promoters_available": bool(context_rows),
            "n_with_context_gene": context_hits,
            "frac_with_context_gene": context_hits / len(scored) if scored else None,
            "note": "Disponibilidad de contexto, no acierto: no se inventa contexto.",
        },
        "classes": dict(Counter(r["class"] for r in scored).most_common()),
        "rows": scored,
        "sigma_features": feature_rows,
    }
    return result


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8")


def _write_rows(path: Path, rows: Iterable[dict]) -> None:
    cols = ["contig", "tss", "strand", "true_sigmas", "pred_top1", "pred_top3",
            "n_true", "class", "jaccard_top3", "no_pwm", "has_pwm",
            "operons", "genes", "evidence", "confidence"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, (list, tuple, set)):
                    v = ",".join(str(x) for x in v)
                vals.append(str(v).replace("\t", " ").replace("\n", " "))
            fh.write("\t".join(vals) + "\n")


def _write_features(path: Path, rows: Iterable[dict]) -> None:
    cols = ["genome", "contig", "tss", "strand", "sigma", "true", "is_top1",
            "z", "pvalue", "score", "spacer", "spacer_dev", "dist_to_start",
            "margin", "top_z", "inventory", "n_true", "gold_confidence", "source"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ")
                                .replace("\n", " ") for c in cols) + "\n")


def _report(result: dict) -> str:
    l3 = result["layer3_assignment"]
    l1 = result["layer1_localization"]
    l0 = result["layer0_data"]
    lines = [
        f"# Evaluación por capas — {result['genome']}", "",
        f"Fuente: `{result['source']}` · ventana {result['config']['window']} · "
        f"{result['config']['decoys']} señuelos · semilla {result['config']['seed']}",
        "", "## Capa 0 — datos", "",
        f"- TSS gold: {l0['n_gold_tss']}",
        f"- Con secuencia: {l0['n_with_sequence']}",
        f"- Sin etiqueta: {l0['n_unlabeled']}",
        f"- Sin PWM: {l0['n_no_pwm']}",
        f"- Sin secuencia: {l0['n_missing_sequence']}",
        "", "## Capa 1 — localización", "",
        f"- TSS con soporte de motivo: {l1['n_supported']}/{l1['n']} "
        f"({l1['frac_supported']:.1%})" if l1['frac_supported'] is not None else "- sin datos",
        "", "## Capa 2 — motivo por sigma", "",
        "| sigma | n | soportados | % | mediana percentil señuelo |",
        "|---|---:|---:|---:|---:|",
    ]
    for sigma, m in sorted(result["layer2_motif"].items()):
        pct = f"{m['pct_supported']:.1%}" if m["pct_supported"] is not None else "—"
        med = f"{m['median_pct_decoy']:.3f}" if m["median_pct_decoy"] is not None else "—"
        lines.append(f"| {sigma} | {m['n']} | {m['supported']} | {pct} | {med} |")
    lines += [
        "", "## Capa 3 — asignación", "",
        f"- Evaluables: {l3['n_evaluable']} ({l3['n_single']} single, {l3['n_multi']} multi)",
        f"- Top-1: {l3['top1_accuracy']:.1%}" if l3['top1_accuracy'] is not None else "- Top-1: —",
        f"- Top-3 cualquiera: {l3['top3_any']:.1%}" if l3['top3_any'] is not None else "- Top-3: —",
        f"- Top-3 todos: {l3['top3_all']:.1%}" if l3['top3_all'] is not None else "- Top-3 todos: —",
        f"- Mayoría trivial: {l3['majority_class']} ({l3['majority_baseline']:.1%})"
        if l3['majority_baseline'] is not None else "- Mayoría: —",
        "", "| sigma | support | precision | recall | F1 |",
        "|---|---:|---:|---:|---:|",
    ]
    for sigma, m in sorted(l3["per_sigma"].items()):
        def f(x): return f"{x:.3f}" if x is not None else "—"
        lines.append(f"| {sigma} | {m['support']} | {f(m['precision'])} | "
                     f"{f(m['recall'])} | {f(m['f1'])} |")
    lines += [
        "", "## Capa 4 — contexto", "",
        f"- `sigma_promoters.tsv` disponible: {result['layer4_context']['sigma_promoters_available']}",
        f"- TSS con gen en contexto: {result['layer4_context']['n_with_context_gene']} "
        f"({result['layer4_context']['frac_with_context_gene']:.1%})"
        if result['layer4_context']['frac_with_context_gene'] is not None else "- contexto: —",
        "", "## Clases", "",
    ]
    for k, v in result["classes"].items():
        lines.append(f"- `{k}`: {v}")
    lines += ["", "## Lectura", "",
              "Un TSS con varios sigma no es un error: se evalúa como multietiqueta.",
              "Un `missed` no es automáticamente falso positivo: puede ser un TSS",
              "nuevo, un proceso superpuesto o una etiqueta incompleta del gold.",
              "La capa 4 solo mide disponibilidad de contexto; no lo inventa.", ""]
    return "\n".join(lines)


def selftest() -> str:
    assert map_sigma("sigma70", "k12") == ("sigma70", "ok")
    assert map_sigma("A", "subtilis") == ("sigma70", "ok")
    assert map_sigma("E", "subtilis") == ("", "no_pwm")
    assert map_sigma("?", "k12") == ("", "unlabeled")
    assert _norm_strand("-1") == "-"
    assert _split_sigmas("sigma70;sigma32") == ["sigma70", "sigma32"]
    return "mapeo, hebra y multietiqueta sigma OK"


def _complete_gold_ids(root: Path) -> set[str]:
    """IDs declarados completos en schema.json (campo genome_id.enum)."""
    schema = root / "schema.json"
    if not schema.exists():
        return set()
    try:
        data = json.loads(schema.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return set()
    for col in data.get("columns", []):
        if col.get("name") == "genome_id":
            return {str(x) for x in (col.get("enum") or [])}
    return set()


def _gold_is_partial(path: Path, complete_ids: set[str] | None = None) -> bool:
    genome = path.parent.name
    if complete_ids is not None and complete_ids and genome not in complete_ids:
        return True
    readme = path.parent / "README.md"
    if not readme.exists():
        return False
    text = readme.read_text(errors="replace")[:4000].lower()
    return any(word in text for word in ("parcial", "partial", "exploratorio", "exploratory"))


def _gold_paths(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted(p for p in root.glob("*/tss_sigma.tsv") if p.is_file())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", help="ID de genoma (ecoli_k12, subtilis, GCF...).")
    ap.add_argument("--source", choices=["auto", "k12", "subtilis", "gold"],
                    default="auto")
    ap.add_argument("--pwm", type=Path, default=None,
                    help="PWM específica; por defecto la de K-12")
    ap.add_argument("--gold-root", type=Path, default=ROOT / "data" / "sigma" / "gold")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--include-partial", action="store_true",
                    help="Incluye gold marcado como parcial en su README.")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "sigma" / "layers")
    ap.add_argument("--up", type=int, default=DEFAULT_WINDOW_UP)
    ap.add_argument("--down", type=int, default=DEFAULT_WINDOW_DOWN)
    ap.add_argument("--decoys", type=int, default=DEFAULT_DECOYS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--zmin", type=float, default=DEFAULT_Z)
    ap.add_argument("--pmax", type=float, default=DEFAULT_P)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0

    # ``--pwm`` debe llegar a ``evaluate``; antes se parseaba pero se ignoraba
    # (PWMS_PATH quedaba siempre en None y se usaba la matriz de K-12).
    global PWMS_PATH
    PWMS_PATH = args.pwm

    data_root = ROOT / "data"
    if args.all:
        paths = _gold_paths(args.gold_root)
        complete_ids = _complete_gold_ids(args.gold_root)
        if not paths:
            print(f"No hay gold en {args.gold_root}; nada que evaluar.", file=sys.stderr)
            return 2
        ok = 0
        for path in paths:
            genome = path.parent.name
            if _gold_is_partial(path, complete_ids) and not args.include_partial:
                print(f"skip parcial: {genome} (usar --include-partial)", file=sys.stderr)
                continue
            try:
                rows = load_gold(path, "gold")
                result = evaluate(rows, data_root, genome, "gold",
                                  args.up, args.down, args.decoys, args.seed,
                                  args.zmin, args.pmax)
                _write_json(args.out / genome / "layer_metrics.json", result)
                _write_rows(args.out / genome / "layer_predictions.tsv", result["rows"])
                _write_features(args.out / genome / "layer_sigma_features.tsv",
                                result.get("sigma_features", []))
                (args.out / genome / "layer_report.md").write_text(_report(result),
                                                                   encoding="utf-8")
                print(f"{genome}: {result['layer3_assignment']['n_evaluable']} evaluables")
                ok += 1
            except Exception as exc:  # noqa: BLE001
                print(f"FALLO {genome}: {exc}", file=sys.stderr)
        return 0 if ok else 1

    if not args.genome:
        ap.error("indica --genome, --all o --selftest")
    source = args.source
    if source == "auto":
        gold = args.gold_root / args.genome / "tss_sigma.tsv"
        if gold.exists():
            source = "gold"
        elif args.genome.lower() in {"ecoli_k12", "ecoli"}:
            source = "k12"
        elif args.genome.lower() == "subtilis":
            source = "subtilis"
        else:
            source = "gold"
    if source == "gold":
        path = args.gold_root / args.genome / "tss_sigma.tsv"
        if not path.exists():
            print(f"Falta {path}; esperando gold curado.", file=sys.stderr)
            return 2
        rows = load_gold(path, "gold")
    elif source == "k12":
        rows = load_k12(data_root)
    else:
        rows = load_subtilis(data_root)
    result = evaluate(rows, data_root, args.genome, source,
                      args.up, args.down, args.decoys, args.seed, args.zmin, args.pmax)
    _write_json(args.out / args.genome / "layer_metrics.json", result)
    _write_rows(args.out / args.genome / "layer_predictions.tsv", result["rows"])
    _write_features(args.out / args.genome / "layer_sigma_features.tsv",
                    result.get("sigma_features", []))
    (args.out / args.genome / "layer_report.md").write_text(_report(result),
                                                            encoding="utf-8")
    print(f"{args.genome}: {result['layer3_assignment']['n_evaluable']} evaluables · "
          f"top1={result['layer3_assignment']['top1_accuracy']:.3f}" if result['layer3_assignment']['top1_accuracy'] is not None else f"{args.genome}: sin evaluables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
