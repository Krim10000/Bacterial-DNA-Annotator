#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Benchmark TSS-level del anotador sigma para K-12 y *Bacillus subtilis*.

El benchmark NO reanota genomas: evalúa el scanner de PWMs sobre TSS curados.
Para cada TSS con etiqueta sigma mapeable:

1. extrae una ventana de 101 nt en orientación de transcripción
   (80 nt aguas arriba + TSS + 20 nt aguas abajo),
2. puntúa los seis sigma con ``score_sigma``,
3. calcula z/p contra un fondo nulo de ventanas aleatorias del mismo genoma,
4. elige top-1 por z-score y compara contra la etiqueta curada.

Uso:
    .venv/bin/python scripts/sigma/30_validate_sigma.py --k12
    .venv/bin/python scripts/sigma/30_validate_sigma.py --subtilis
    .venv/bin/python scripts/sigma/30_validate_sigma.py --all
    .venv/bin/python scripts/sigma/30_validate_sigma.py --selftest

Salidas por defecto en ``data/sigma/benchmark/``:
    k12_tss_predictions.tsv       una fila por TSS curado de K-12
    k12_confusion.json            matriz de confusión TSS-level
    k12_metrics.json              métricas, cobertura, calibración, comparación
    subtilis_tss_predictions.tsv  una fila por TSS curado de B. subtilis
    subtilis_confusion.json
    subtilis_metrics.json
    benchmark_summary.md          resumen legible

No descarga datos, no usa red y no escribe fuera de ``data/sigma/benchmark/``.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma.sigma_annotator import (  # noqa: E402
    DEFAULT_MAX_P,
    DEFAULT_MIN_MARGIN,
    DEFAULT_MIN_Z,
    DEFAULT_SEED,
    SIGMA_ORDER,
    _confidence,
    _load_fasta,
    _pvalue,
    _read_tsv,
    _write_tsv,
    _zscore,
    build_null,
    load_pwm_library,
    load_sequence,
    revcomp,
    score_sigma,
    selftest as scanner_selftest,
)

# --------------------------------------------------------------------------
# Configuración del benchmark
# --------------------------------------------------------------------------
TSS_UPSTREAM = 80
TSS_DOWNSTREAM = 20
TSS_WINDOW = TSS_UPSTREAM + TSS_DOWNSTREAM + 1  # 101 nt (-80..+20, TSS en 80)
BENCH_SEED = DEFAULT_SEED                         # 17
BENCH_NULL_WINDOWS = 500
SIGMA_LABELS = tuple(SIGMA_ORDER)

# Mapa conservador de etiquetas de B. subtilis a las seis PWMs de E. coli.
# No se mapean sigma de esporulación (E/F/G/K), competencia/transición (H)
# ni ECF (W/X/M/V/Y/Z): no hay una PWM equivalente fiable en el repositorio.
BSU_SIGMA_MAP: dict[str, tuple[str, str]] = {
    "A": ("sigma70", "SigA: sigma housekeeping, ortólogo funcional de RpoD"),
    "B": ("sigma38", "SigB: sigma de estrés general, análogo funcional de RpoS"),
    "D": ("sigma28", "SigD: sigma flagelar, ortólogo funcional de FliA"),
    "I": ("sigma32", "SigI: sigma de choque térmico, análogo funcional de RpoH"),
    "L": ("sigma54", "SigL: familia sigma54, relacionada con metabolismo de nitrógeno"),
}
BSU_NO_PWM_REASON = {
    "E": "SigE: esporulación; sin PWM equivalente",
    "F": "SigF: esporulación; sin PWM equivalente",
    "G": "SigG: esporulación; sin PWM equivalente",
    "K": "SigK: esporulación; sin PWM equivalente",
    "H": "SigH: competencia/transición; sin PWM equivalente fiable",
    "W": "SigW: ECF; la PWM ECF de E. coli no es transferible",
    "X": "SigX: ECF; la PWM ECF de E. coli no es transferible",
    "M": "SigM: ECF; la PWM ECF de E. coli no es transferible",
    "V": "SigV: ECF; sin PWM equivalente",
    "Y": "SigY: ECF; sin PWM equivalente",
    "Z": "SigZ: ECF; sin PWM equivalente",
}

PROMOTER_NAME_RE = re.compile(r"^(.*?)p\d*$")
PREDICTION_BASE_COLS = [
    "dataset", "tss_id", "source", "contig", "tss", "strand",
    "sigma_raw", "mapping_status", "mapping_reason", "label_status",
    "sigma_truth", "gene_key", "promoter_name", "operon",
    "pred_sigma", "pred_z", "pred_pvalue", "pred_score", "pred_margin",
    "pred_rank", "top3_sigmas", "correct_top1", "correct_top3",
    "confidence", "confidence_label", "confidence_score",
]
PREDICTION_COLS = PREDICTION_BASE_COLS + [
    f"{metric}_{sigma}" for sigma in SIGMA_LABELS for metric in ("z", "p")
]


# --------------------------------------------------------------------------
# Utilidades
# --------------------------------------------------------------------------
def _clean(text: Any) -> str:
    return str(text or "").strip()


def _to_int(value: Any) -> int | None:
    text = _clean(value).replace(",", "")
    if not text:
        return None
    try:
        number = int(text)
    except ValueError:
        try:
            number = int(float(text))
        except ValueError:
            return None
    return number


def _round(value: float | None, digits: int = 6) -> float | None:
    if value is None:
        return None
    return round(float(value), digits)


def _safe_div(num: float, den: float) -> float | None:
    if not den:
        return None
    return float(num) / float(den)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8")


def _sanitize_tsv_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    clean = []
    for row in rows:
        clean.append({key: ("" if value is None else value) for key, value in row.items()})
    return clean


# --------------------------------------------------------------------------
# Mapeo de etiquetas
# --------------------------------------------------------------------------
def map_k12_sigma(raw: Any, available: Iterable[str] = SIGMA_LABELS) -> dict[str, Any]:
    """Mapea una etiqueta de ``gt_promoters.tsv`` a una PWM disponible.

    Devuelve ``mapping_status``:
      evaluable  -> hay PWM y la etiqueta entra en el benchmark
      no_pwm     -> etiqueta sigma válida sin PWM equivalente
      unlabeled  -> la fila no tiene etiqueta sigma
    """
    label = _clean(raw)
    available_set = set(available)
    if not label:
        return {"mapping_status": "unlabeled", "sigma_truth": "",
                "mapping_reason": "sin etiqueta sigma"}
    if label in available_set:
        return {"mapping_status": "evaluable", "sigma_truth": label,
                "mapping_reason": ""}
    if label == "sigma19":
        return {"mapping_status": "no_pwm", "sigma_truth": "",
                "mapping_reason": "sigma19 no tiene PWM en promoter_pwms.tsv"}
    return {"mapping_status": "no_pwm", "sigma_truth": "",
            "mapping_reason": f"sin PWM para {label}"}


def map_bsu_sigma(raw: Any, available: Iterable[str] = SIGMA_LABELS) -> dict[str, Any]:
    """Mapea una etiqueta de B. subtilis (``A`` o ``SigA``) a PWM de E. coli."""
    label = _clean(raw)
    available_set = set(available)
    if not label or label in {"?", ".", "unknown", "NA", "None"}:
        return {"mapping_status": "unlabeled", "sigma_truth": "", "bsu_label": "",
                "mapping_reason": "sin etiqueta sigma"}
    bsu_label = label
    if bsu_label.lower().startswith("sig"):
        bsu_label = bsu_label[3:].strip()
    if bsu_label in BSU_SIGMA_MAP:
        mapped, note = BSU_SIGMA_MAP[bsu_label]
        if mapped in available_set:
            return {"mapping_status": "evaluable", "sigma_truth": mapped,
                    "bsu_label": bsu_label, "mapping_reason": note}
        return {"mapping_status": "no_pwm", "sigma_truth": "",
                "bsu_label": bsu_label,
                "mapping_reason": f"{note}; PWM {mapped} no disponible"}
    reason = BSU_NO_PWM_REASON.get(
        bsu_label, f"etiqueta {label!r} sin PWM equivalente")
    return {"mapping_status": "no_pwm", "sigma_truth": "",
            "bsu_label": bsu_label, "mapping_reason": reason}


def _k12_gene_key(row: dict[str, str]) -> str:
    """Nombre de gen para el cruce operón/gen (first_gene o nombre del promotor)."""
    first = _clean(row.get("first_gene"))
    if first:
        return first
    name = _clean(row.get("name"))
    match = PROMOTER_NAME_RE.match(name)
    if match:
        return match.group(1).strip()
    return name


# --------------------------------------------------------------------------
# Secuencia y ventanas
# --------------------------------------------------------------------------
def load_genome_sequence(data_root: Path, genome: str) -> tuple[dict[str, str], str]:
    """Carga la secuencia con ``load_sequence`` y fallbacks explícitos.

    ``load_sequence`` depende de ``data/regulation/<genoma>/genes_meta.json``;
    B. subtilis no tiene ese metadato, así que se usa ``data/dbtbs/genome.fna``.
    Devuelve ``(secuencias, fuente)``.
    """
    seqs = load_sequence(data_root, genome)
    if seqs:
        return seqs, "load_sequence"

    candidates: list[Path] = []
    if genome == "subtilis":
        candidates += [
            data_root / "dbtbs" / "genome.fna",
            data_root / "dbtbs" / "subtilis.gbff",
        ]
    if genome == "ecoli_k12":
        candidates += [
            data_root / "allostery" / "work" / "ecoli_k12.fna",
            data_root / "test_ecoli.gbff",
        ]
    candidates += [
        data_root / "annotations" / f"{genome}_genomic.fna",
        data_root / "annotations" / f"{genome}_genomic.fna.gz",
        data_root / "annotations" / f"{genome}_genomic.gbff",
        data_root / "annotations" / f"{genome}_genomic.gbff.gz",
    ]
    for path in candidates:
        if not path.exists():
            continue
        suffix = path.name.lower()
        loaded: dict[str, str] = {}
        if suffix.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz")):
            loaded = _load_fasta(path)
        elif suffix.endswith((".gbff", ".gbk", ".gbff.gz", ".gbk.gz")):
            try:
                from Bio import SeqIO
            except ImportError:
                continue
            opener = __import__("gzip").open if suffix.endswith(".gz") else open
            with opener(path, "rt") as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    loaded[record.id] = str(record.seq).upper()
        if loaded:
            return loaded, str(path)
    return {}, ""


def tss_window(seq: str, pos: int | None, strand: str,
               upstream: int = TSS_UPSTREAM,
               downstream: int = TSS_DOWNSTREAM) -> str:
    """Ventana orientada a transcripción: -upstream...+downstream del TSS.

    ``pos`` es 1-based. El TSS queda en el índice ``upstream`` de la ventana
    devuelta, tanto en hebra + como en hebra - (revcomp). Longitud esperada:
    ``upstream + downstream + 1``.
    """
    if not seq or pos is None or strand not in ("+", "-"):
        return ""
    index = int(pos) - 1
    if index < 0:
        return ""
    if strand == "+":
        start = index - upstream
        end = index + downstream + 1
        if start < 0 or end > len(seq):
            return ""
        return seq[start:end]
    start = index - downstream
    end = index + upstream + 1
    if start < 0 or end > len(seq):
        return ""
    return revcomp(seq[start:end])


# --------------------------------------------------------------------------
# Carga de TSS curados
# --------------------------------------------------------------------------
def load_k12_tss(data_root: Path,
                 pwms: dict[str, dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    path = data_root / "regulondb" / "gt" / "gt_promoters.tsv"
    rows = _read_tsv(path)
    available = tuple(pwms) if pwms else SIGMA_LABELS
    out: list[dict[str, Any]] = []
    for index, row in enumerate(rows, 1):
        raw = _clean(row.get("sigma_factor"))
        mapped = map_k12_sigma(raw, available)
        tss = _to_int(row.get("tss"))
        status = "missing_tss" if tss is None else mapped["mapping_status"]
        item: dict[str, Any] = {
            "dataset": "k12",
            "tss_id": _clean(row.get("promoter_id")) or f"k12:{index}",
            "source": "regulondb_gt_promoters",
            "contig": "",
            "tss": tss if tss is not None else "",
            "strand": _clean(row.get("strand")),
            "sigma_raw": raw,
            "mapping_status": mapped["mapping_status"],
            "mapping_reason": mapped["mapping_reason"],
            "label_status": status,
            "sigma_truth": mapped["sigma_truth"],
            "gene_key": _k12_gene_key(row),
            "promoter_name": _clean(row.get("name")),
            "operon": "",
            "confidence_gt": _clean(row.get("confidence")),
            "evidence": _clean(row.get("first_gene")),
        }
        out.append(item)
    return out


def load_bsu_tss(data_root: Path,
                 pwms: dict[str, dict[str, Any]] | None = None
                 ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Carga ``gt_tss.tsv`` y usa ``dbtbs_tss.tsv`` solo para enriquecer operón."""
    gt_path = data_root / "dbtbs" / "gt_tss.tsv"
    rows = _read_tsv(gt_path)
    available = tuple(pwms) if pwms else SIGMA_LABELS
    db_rows = _read_tsv(data_root / "dbtbs" / "dbtbs_tss.tsv")
    db_by_key: dict[tuple[int, str], dict[str, str]] = {}
    for row in db_rows:
        pos = _to_int(row.get("pos"))
        strand = _clean(row.get("strand"))
        if pos is not None and strand in ("+", "-") and (pos, strand) not in db_by_key:
            db_by_key[(pos, strand)] = row
    gt_keys: set[tuple[int, str]] = set()
    out: list[dict[str, Any]] = []
    for index, row in enumerate(rows, 1):
        raw = _clean(row.get("sigma"))
        mapped = map_bsu_sigma(raw, available)
        tss = _to_int(row.get("pos"))
        strand = _clean(row.get("strand"))
        status = "missing_tss" if tss is None else mapped["mapping_status"]
        evidence = _clean(row.get("evidence"))
        if tss is not None and strand in ("+", "-"):
            gt_keys.add((tss, strand))
        db_row = db_by_key.get((int(tss), strand)) if tss is not None else None
        source_bits = []
        if "DBTBS" in evidence:
            source_bits.append("DBTBS")
        if "BsubCyc" in evidence:
            source_bits.append("BsubCyc")
        if "Nicolas" in evidence:
            source_bits.append("Nicolas")
        item: dict[str, Any] = {
            "dataset": "subtilis",
            "tss_id": f"bsu:{tss}:{strand}" if tss is not None else f"bsu:{index}",
            "source": "+".join(source_bits) or "BSGatlas/gt_tss",
            "contig": "",
            "tss": tss if tss is not None else "",
            "strand": strand,
            "sigma_raw": raw,
            "mapping_status": mapped["mapping_status"],
            "mapping_reason": mapped["mapping_reason"],
            "label_status": status,
            "sigma_truth": mapped["sigma_truth"],
            "gene_key": "",
            "promoter_name": "",
            "operon": _clean(db_row.get("operon")) if db_row else "",
            "confidence_gt": "",
            "evidence": evidence,
            "operon_source": "dbtbs_tss" if db_row and _clean(db_row.get("operon")) else "",
        }
        out.append(item)
    db_keys = set(db_by_key)
    secondary = {
        "path": str(data_root / "dbtbs" / "dbtbs_tss.tsv"),
        "n_rows": len(db_rows),
        "n_unique_pos_strand": len(db_keys),
        "n_overlap_with_gt": len(db_keys & gt_keys),
        "n_only_in_dbtbs": len(db_keys - gt_keys),
        "used_for": "enriquecer columna operon; no se agregan TSS extra al ground truth",
    }
    return out, secondary


# --------------------------------------------------------------------------
# Scoring TSS-level
# --------------------------------------------------------------------------
def _empty_prediction_fields() -> dict[str, Any]:
    fields: dict[str, Any] = {
        "pred_sigma": "", "pred_z": "", "pred_pvalue": "", "pred_score": "",
        "pred_margin": "", "pred_rank": "", "top3_sigmas": "",
        "correct_top1": "", "correct_top3": "",
        "confidence": "", "confidence_label": "", "confidence_score": "",
    }
    for sigma in SIGMA_LABELS:
        fields[f"z_{sigma}"] = ""
        fields[f"p_{sigma}"] = ""
    return fields


def _inventory_for_confidence(data_root: Path, genome: str) -> set[str]:
    """Tipos sigma del inventario, con el alias sigmaE -> sigma24."""
    path = data_root / "sigma" / genome / "sigma_factors.tsv"
    out: set[str] = set()
    for row in _read_tsv(path):
        sigma = _clean(row.get("sigma_type"))
        if sigma == "sigmaE":
            sigma = "sigma24"
        if sigma in SIGMA_LABELS:
            out.add(sigma)
    return out


def score_tss_rows(rows: list[dict[str, Any]], seqs: dict[str, str],
                   pwms: dict[str, dict[str, Any]],
                   null: dict[str, dict[str, Any]], *,
                   inventory: set[str] | None = None,
                   upstream: int = TSS_UPSTREAM,
                   downstream: int = TSS_DOWNSTREAM,
                   max_p: float = DEFAULT_MAX_P,
                   min_z: float = DEFAULT_MIN_Z,
                   min_margin: float = DEFAULT_MIN_MARGIN) -> list[dict[str, Any]]:
    """Puntúa cada TSS evaluable y devuelve las filas con predicción/rank."""
    if not seqs:
        for row in rows:
            row.update(_empty_prediction_fields())
        return rows
    contig = max(seqs, key=lambda name: len(seqs[name]))
    seq = seqs[contig]
    window_len = upstream + downstream + 1
    inventory = set(inventory or ())

    for row in rows:
        row.update(_empty_prediction_fields())
        row["contig"] = contig
        if row.get("label_status") in {"missing_tss", "unlabeled", "no_pwm"}:
            continue
        if not seq:
            row["label_status"] = "no_sequence"
            continue
        tss = _to_int(row.get("tss"))
        window = tss_window(seq, tss, _clean(row.get("strand")),
                            upstream=upstream, downstream=downstream)
        if not window or len(window) != window_len:
            row["label_status"] = "no_window"
            continue

        hits: dict[str, dict[str, Any]] = {}
        for sigma in SIGMA_LABELS:
            rec = pwms.get(sigma)
            if not rec:
                continue
            hit = score_sigma(window, rec, len(window))
            if hit is None:
                continue
            hit["z"] = _zscore(hit["score"], null.get(sigma))
            hit["pvalue"] = _pvalue(hit["score"], null.get(sigma))
            hits[sigma] = hit
            row[f"z_{sigma}"] = _round(hit["z"], 4)
            row[f"p_{sigma}"] = _round(hit["pvalue"], 6)
        if not hits:
            row["label_status"] = "score_failed"
            continue

        ranked = sorted(
            hits.values(),
            key=lambda hit: (hit["z"], -hit["pvalue"], hit["score"]),
            reverse=True,
        )
        top = ranked[0]
        second = ranked[1] if len(ranked) > 1 else None
        margin = float(top["z"] - second["z"]) if second else 99.0
        truth = _clean(row.get("sigma_truth"))
        rank = next((i + 1 for i, hit in enumerate(ranked)
                     if hit["sigma"] == truth), None)
        confidence, label, cscore = _confidence(
            float(top["pvalue"]), float(top["z"]), margin,
            top["sigma"] in inventory, False, "not_available",
            max_p, min_z, min_margin)

        row.update({
            "label_status": "evaluable",
            "pred_sigma": top["sigma"],
            "pred_z": _round(top["z"], 4),
            "pred_pvalue": _round(top["pvalue"], 6),
            "pred_score": _round(top["score"], 4),
            "pred_margin": _round(margin, 4),
            "pred_rank": rank if rank is not None else "",
            "top3_sigmas": ";".join(hit["sigma"] for hit in ranked[:3]),
            "correct_top1": int(top["sigma"] == truth) if truth else "",
            "correct_top3": int(rank is not None and rank <= 3) if truth else "",
            "confidence": confidence,
            "confidence_label": label,
            "confidence_score": cscore,
        })
    return rows


# --------------------------------------------------------------------------
# Métricas
# --------------------------------------------------------------------------
def _confusion(evaluable: list[dict[str, Any]],
               labels: tuple[str, ...]) -> dict[str, Any]:
    index = {label: i for i, label in enumerate(labels)}
    matrix = [[0 for _ in labels] for _ in labels]
    for row in evaluable:
        truth = _clean(row.get("sigma_truth"))
        pred = _clean(row.get("pred_sigma"))
        if truth in index and pred in index:
            matrix[index[truth]][index[pred]] += 1
    row_totals = {label: sum(matrix[i]) for i, label in enumerate(labels)}
    col_totals = {label: sum(matrix[i][j] for i in range(len(labels)))
                  for j, label in enumerate(labels)}
    normalized = []
    for i in range(len(labels)):
        total = row_totals[labels[i]]
        normalized.append([_round(cell / total, 6) if total else None
                           for cell in matrix[i]])
    return {
        "labels": list(labels),
        "matrix": matrix,
        "row_totals": row_totals,
        "col_totals": col_totals,
        "normalized_by_row": normalized,
        "n": len(evaluable),
    }


def _calibration(evaluable: list[dict[str, Any]]) -> dict[str, Any]:
    if not evaluable:
        return {"status": "not_available", "reason": "sin predicciones evaluables"}

    def _stats(subset: list[dict[str, Any]]) -> dict[str, Any]:
        n = len(subset)
        correct = sum(1 for row in subset if int(row.get("correct_top1") or 0))
        scores = [float(row.get("confidence_score") or 0.0) for row in subset]
        mean_score = sum(scores) / n if n else None
        accuracy = correct / n if n else None
        return {
            "n": n,
            "correct": correct,
            "accuracy": _round(accuracy, 6),
            "mean_confidence_score": _round(mean_score, 6),
            "gap": _round(accuracy - mean_score, 6)
            if accuracy is not None and mean_score is not None else None,
        }

    by_label = {}
    for label in ("sospecha_fuerte", "sospecha_media", "sospecha_baja", "sin_sospecha"):
        subset = [row for row in evaluable if _clean(row.get("confidence_label")) == label]
        if subset:
            by_label[label] = _stats(subset)

    bins = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0001)]
    by_score = []
    for low, high in bins:
        subset = [row for row in evaluable
                  if low <= float(row.get("confidence_score") or 0.0) < high]
        if not subset:
            continue
        stats = _stats(subset)
        stats.update({"low": low, "high": min(high, 1.0)})
        by_score.append(stats)

    z_bins = [(-math.inf, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, math.inf)]
    by_z = []
    for low, high in z_bins:
        subset = [row for row in evaluable
                  if low <= float(row.get("pred_z") or -math.inf) < high]
        if not subset:
            continue
        stats = _stats(subset)
        stats.update({"low": None if math.isinf(low) else low,
                      "high": None if math.isinf(high) else high})
        by_z.append(stats)

    total = len(evaluable)
    ece = sum((item["n"] / total) * abs(item["gap"] or 0.0) for item in by_score)
    brier = sum((float(row.get("confidence_score") or 0.0)
                 - int(row.get("correct_top1") or 0)) ** 2
                for row in evaluable) / total
    return {
        "status": "available",
        "n": total,
        "by_confidence_label": by_label,
        "by_score_bin": by_score,
        "by_z_bin": by_z,
        "expected_calibration_error": _round(ece, 6),
        "brier_score": _round(brier, 6),
        "note": ("confidence_score NO es una probabilidad calibrada; "
                 "ECE/Brier se reportan como diagnóstico, no como probabilidad."),
    }


def compute_metrics(rows: list[dict[str, Any]], dataset: str,
                    labels: tuple[str, ...] = SIGMA_LABELS,
                    upstream: int = TSS_UPSTREAM,
                    downstream: int = TSS_DOWNSTREAM) -> dict[str, Any]:
    """Calcula cobertura, accuracy, P/R/F1, confusión y calibración."""
    status_counts = Counter(_clean(row.get("label_status")) for row in rows)
    mapping_counts = Counter(_clean(row.get("mapping_status")) for row in rows)
    n_rows = len(rows)
    n_missing = status_counts.get("missing_tss", 0)
    n_valid = n_rows - n_missing
    n_unlabeled = sum(1 for row in rows
                      if row.get("mapping_status") == "unlabeled"
                      and row.get("label_status") != "missing_tss")
    n_no_pwm = sum(1 for row in rows
                   if row.get("mapping_status") == "no_pwm"
                   and row.get("label_status") != "missing_tss")
    n_evaluable = sum(1 for row in rows if row.get("label_status") == "evaluable")
    n_no_window = status_counts.get("no_window", 0)
    n_score_failed = status_counts.get("score_failed", 0)
    n_no_sequence = status_counts.get("no_sequence", 0)
    n_labeled_valid = n_evaluable + n_no_pwm + n_no_window + n_score_failed + n_no_sequence
    evaluable = [row for row in rows if row.get("label_status") == "evaluable"]

    confusion = _confusion(evaluable, labels)
    correct1 = sum(1 for row in evaluable if int(row.get("correct_top1") or 0))
    correct3 = sum(1 for row in evaluable if int(row.get("correct_top3") or 0))
    truth_counts = Counter(_clean(row.get("sigma_truth")) for row in evaluable)
    pred_counts = Counter(_clean(row.get("pred_sigma")) for row in evaluable)
    majority = truth_counts.most_common(1)[0] if truth_counts else ("", 0)
    majority_baseline = _safe_div(majority[1], n_evaluable)
    accuracy1 = _safe_div(correct1, n_evaluable)
    accuracy3 = _safe_div(correct3, n_evaluable)

    per_sigma: dict[str, Any] = {}
    for sigma in labels:
        tp = sum(1 for row in evaluable
                 if row.get("sigma_truth") == sigma and row.get("pred_sigma") == sigma)
        fp = sum(1 for row in evaluable
                 if row.get("sigma_truth") != sigma and row.get("pred_sigma") == sigma)
        fn = sum(1 for row in evaluable
                 if row.get("sigma_truth") == sigma and row.get("pred_sigma") != sigma)
        support = tp + fn
        predicted = tp + fp
        precision = _safe_div(tp, predicted)
        recall = _safe_div(tp, support)
        f1 = (2 * precision * recall / (precision + recall)
              if precision is not None and recall is not None
              and (precision + recall) > 0 else None)
        per_sigma[sigma] = {
            "support": support,
            "predicted": predicted,
            "tp": tp, "fp": fp, "fn": fn,
            "precision": _round(precision, 6),
            "recall": _round(recall, 6),
            "f1": _round(f1, 6),
        }

    supported = [item for item in per_sigma.values() if item["support"] > 0]
    macro_supported = _safe_div(
        sum(item["f1"] for item in supported if item["f1"] is not None), len(supported)
    ) if supported else None
    # Macro conservador: las clases sin F1 definido cuentan como 0.
    macro_all = sum(item["f1"] or 0.0 for item in per_sigma.values()) / len(labels)
    weighted = _safe_div(
        sum((item["f1"] or 0.0) * item["support"] for item in per_sigma.values()),
        n_evaluable,
    )
    ranks = [int(row["pred_rank"]) for row in evaluable
             if str(row.get("pred_rank", "")).strip().isdigit()]
    no_pwm_counts = Counter(
        raw for row in rows
        if row.get("mapping_status") == "no_pwm"
        and row.get("label_status") != "missing_tss"
        for raw in [_clean(row.get("sigma_raw"))] if raw
    )

    metrics: dict[str, Any] = {
        "dataset": dataset,
        "status": "evaluated" if n_evaluable else "not_available",
        "label_space": list(labels),
        "window": {"upstream": upstream, "downstream": downstream,
                   "length": upstream + downstream + 1,
                   "orientation": "transcripción",
                   "tss_index_0based": upstream},
        "counts": {
            "n_rows": n_rows,
            "n_valid_tss": n_valid,
            "n_missing_tss": n_missing,
            "n_unlabeled": n_unlabeled,
            "n_no_pwm": n_no_pwm,
            "n_evaluable": n_evaluable,
            "n_no_window": n_no_window,
            "n_score_failed": n_score_failed,
            "n_no_sequence": n_no_sequence,
            "n_labeled_valid": n_labeled_valid,
            "n_distinct_truth_labels": len(truth_counts),
        },
        "status_counts": dict(status_counts),
        "mapping_status_counts": dict(mapping_counts),
        "coverage": {
            "evaluable_over_labeled_valid": _round(
                _safe_div(n_evaluable, n_labeled_valid), 6),
            "evaluable_over_valid_tss": _round(_safe_div(n_evaluable, n_valid), 6),
            "evaluable_over_all_rows": _round(_safe_div(n_evaluable, n_rows), 6),
            "no_pwm_over_labeled_valid": _round(
                _safe_div(n_no_pwm, n_labeled_valid), 6),
            "no_pwm_over_all_rows": _round(_safe_div(n_no_pwm, n_rows), 6),
            "unlabeled_over_all_rows": _round(_safe_div(n_unlabeled, n_rows), 6),
        },
        "accuracy": {
            "top1": _round(accuracy1, 6),
            "top3": _round(accuracy3, 6),
            "n_top1_correct": correct1,
            "n_top3_correct": correct3,
            "majority_class": majority[0],
            "majority_baseline": _round(majority_baseline, 6),
            "gain_over_majority": _round(
                (accuracy1 - majority_baseline)
                if accuracy1 is not None and majority_baseline is not None else None,
                6),
            "mean_truth_rank": _round(sum(ranks) / len(ranks), 6) if ranks else None,
        },
        "truth_counts": dict(truth_counts.most_common()),
        "prediction_counts": dict(pred_counts.most_common()),
        "no_pwm_counts": dict(no_pwm_counts.most_common()),
        "per_sigma": per_sigma,
        "macro_f1_supported": _round(macro_supported, 6),
        "macro_f1_all6": _round(macro_all, 6),
        "weighted_f1": _round(weighted, 6),
        "confusion_labels": list(labels),
        "confusion_matrix": confusion["matrix"],
        "calibration": _calibration(evaluable),
        "limitations": [
            "TSS-level: no equivale a validar promotores de operones completos.",
            "Solo seis tipos sigma con PWM; las etiquetas sin PWM se reportan como no_pwm.",
            "El fondo nulo es genómico y no un modelo de promotores reales.",
            "confidence_score es una heurística del módulo, no una probabilidad.",
            "B. subtilis usa PWMs de E. coli sobre etiquetas homologadas.",
        ],
    }
    return metrics


# --------------------------------------------------------------------------
# Comparación K-12 contra la salida existente por operón/gen
# --------------------------------------------------------------------------
def compare_k12_existing(data_root: Path,
                         tss_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Cruza ``data/sigma/ecoli_k12/sigma_promoters.tsv`` con el GT por gen.

    La salida existente asigna un único sigma por operón. El GT puede tener
    varios TSS (y varios sigma) para el mismo primer gen. Por eso se reportan
    dos coincidencias: ``exact`` (el sigma predicho está entre los curados para
    ese gen) y ``dominant`` (coincide con el sigma curado más frecuente).
    """
    existing_path = data_root / "sigma" / "ecoli_k12" / "sigma_promoters.tsv"
    if not existing_path.exists():
        return {"status": "not_available", "path": str(existing_path),
                "reason": "no existe la salida existente"}

    gt_by_gene: dict[str, dict[str, Any]] = {}
    for row in tss_rows:
        if row.get("label_status") == "missing_tss":
            continue
        key = _clean(row.get("gene_key"))
        if not key:
            continue
        gene = gt_by_gene.setdefault(key, {
            "n_tss": 0, "n_mappable": 0, "n_no_pwm": 0,
            "labels": Counter(), "raw_labels": Counter(),
        })
        gene["n_tss"] += 1
        raw = _clean(row.get("sigma_raw"))
        if raw:
            gene["raw_labels"][raw] += 1
        if row.get("mapping_status") == "evaluable":
            gene["n_mappable"] += 1
            gene["labels"][_clean(row.get("sigma_truth"))] += 1
        elif row.get("mapping_status") == "no_pwm":
            gene["n_no_pwm"] += 1

    existing_rows = _read_tsv(existing_path)
    matched: list[dict[str, Any]] = []
    unmatched: list[str] = []
    for row in existing_rows:
        key = _clean(row.get("first_gene_name"))
        if not key:
            label = _clean(row.get("operon_label"))
            key = label.split("-")[0].strip() if label else ""
        if key not in gt_by_gene or not gt_by_gene[key]["labels"]:
            unmatched.append(key or "<sin_gen>")
            continue
        gene = gt_by_gene[key]
        pred = _clean(row.get("sigma"))
        labels = set(gene["labels"])
        dominant = gene["labels"].most_common(1)[0][0]
        matched.append({
            "gene": key,
            "pred_sigma": pred,
            "gt_labels": sorted(labels),
            "gt_dominant": dominant,
            "n_tss": gene["n_tss"],
            "n_mappable_tss": gene["n_mappable"],
            "n_no_pwm_tss": gene["n_no_pwm"],
            "exact": pred in labels,
            "dominant_match": pred == dominant,
            "existing_z": _clean(row.get("zscore")),
            "existing_confidence": _clean(row.get("confidence_label")),
        })

    n_matched = len(matched)
    n_exact = sum(1 for item in matched if item["exact"])
    n_dominant = sum(1 for item in matched if item["dominant_match"])
    single = [item for item in matched if len(item["gt_labels"]) == 1]
    multi = [item for item in matched if len(item["gt_labels"]) > 1]
    single_exact = sum(1 for item in single if item["exact"])
    multi_any = sum(1 for item in multi if item["exact"])
    multi_dominant = sum(1 for item in multi if item["dominant_match"])
    multi_all = sum(1 for item in multi if item["exact"] and len(item["gt_labels"]) == 1)

    multi_sigma_genes = [gene for gene, data in gt_by_gene.items()
                         if len(data["labels"]) > 1]
    multi_tss_genes = [gene for gene, data in gt_by_gene.items()
                       if data["n_tss"] > 1]
    matched_genes = {item["gene"] for item in matched}
    gt_genes_with_mappable = {gene for gene, data in gt_by_gene.items() if data["labels"]}

    def _rate(num: int, den: int) -> float | None:
        return _round(_safe_div(num, den), 6)

    gene_confusion = _confusion(
        [{"sigma_truth": item["gt_dominant"], "pred_sigma": item["pred_sigma"]}
         for item in matched],
        SIGMA_LABELS,
    )
    by_sigma = {}
    for sigma in SIGMA_LABELS:
        subset = [item for item in matched if item["pred_sigma"] == sigma]
        if subset:
            by_sigma[sigma] = {
                "n": len(subset),
                "exact": sum(1 for item in subset if item["exact"]),
                "exact_rate": _rate(sum(1 for item in subset if item["exact"]), len(subset)),
            }

    examples = sorted(
        matched,
        key=lambda item: (not item["exact"], -len(item["gt_labels"]), item["gene"]),
    )[:12]
    return {
        "status": "evaluated",
        "path": str(existing_path),
        "mapping_method": ("first_gene_name de sigma_promoters.tsv contra "
                           "first_gene o nombre de promotor parseado de gt_promoters.tsv"),
        "n_existing_rows": len(existing_rows),
        "n_existing_matched": n_matched,
        "n_existing_unmatched": len(unmatched),
        "n_gt_genes_with_mappable": len(gt_genes_with_mappable),
        "n_gt_genes_mapped_to_existing": len(gt_genes_with_mappable & matched_genes),
        "n_gt_genes_without_existing": len(gt_genes_with_mappable - matched_genes),
        "n_gt_genes_with_multiple_tss": len(multi_tss_genes),
        "n_gt_genes_with_multiple_sigma": len(multi_sigma_genes),
        "pct_gt_genes_with_multiple_sigma": _rate(
            len(multi_sigma_genes), len(gt_genes_with_mappable)),
        "exact_matches": n_exact,
        "exact_rate": _rate(n_exact, n_matched),
        "dominant_matches": n_dominant,
        "dominant_rate": _rate(n_dominant, n_matched),
        "single_label": {
            "n_genes": len(single),
            "exact": single_exact,
            "exact_rate": _rate(single_exact, len(single)),
        },
        "multi_label": {
            "n_genes": len(multi),
            "any_exact": multi_any,
            "any_exact_rate": _rate(multi_any, len(multi)),
            "dominant_exact": multi_dominant,
            "dominant_exact_rate": _rate(multi_dominant, len(multi)),
            "all_labels_exact": multi_all,
            "mean_distinct_labels": _round(
                sum(len(item["gt_labels"]) for item in multi) / len(multi), 6)
            if multi else None,
        },
        "one_sigma_per_operon_limit": {
            "n_multi_sigma_genes": len(multi_sigma_genes),
            "n_multi_sigma_genes_matched": sum(
                1 for item in multi if item["gene"] in multi_sigma_genes),
            "n_multi_sigma_genes_any_match": multi_any,
            "n_multi_sigma_genes_all_labels_match": multi_all,
            "note": ("una fila por operón solo puede representar un sigma; "
                     "en genes con ≥2 sigma curados la coincidencia con todos "
                     "es imposible por construcción"),
        },
        "existing_sigma_counts": dict(Counter(_clean(row.get("sigma"))
                                              for row in existing_rows).most_common()),
        "gt_tss_sigma_counts": dict(Counter(
            _clean(row.get("sigma_truth")) for row in tss_rows
            if row.get("mapping_status") == "evaluable").most_common()),
        "by_predicted_sigma": by_sigma,
        "gene_level_confusion": gene_confusion,
        "examples": examples,
        "unmatched_examples": unmatched[:12],
    }


# --------------------------------------------------------------------------
# Orquestación y salidas
# --------------------------------------------------------------------------
def _mapping_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_raw: dict[str, dict[str, Any]] = {}
    for row in rows:
        raw = _clean(row.get("sigma_raw")) or "<vacío>"
        item = by_raw.setdefault(raw, {
            "n": 0, "mapping_status": _clean(row.get("mapping_status")),
            "sigma_truth": _clean(row.get("sigma_truth")),
            "reason": _clean(row.get("mapping_reason")),
        })
        item["n"] += 1
    return {
        "by_raw_label": dict(sorted(by_raw.items())),
        "n_mapped": sum(item["n"] for item in by_raw.values()
                        if item["mapping_status"] == "evaluable"),
        "n_no_pwm": sum(item["n"] for item in by_raw.values()
                        if item["mapping_status"] == "no_pwm"),
        "n_unlabeled": sum(item["n"] for item in by_raw.values()
                           if item["mapping_status"] == "unlabeled"),
    }


def _run_dataset(dataset: str, data_root: Path, out_dir: Path, pwms: dict[str, Any],
                 seed: int, null_windows: int,
                 upstream: int, downstream: int) -> dict[str, Any]:
    if dataset == "k12":
        rows = load_k12_tss(data_root, pwms)
        secondary = None
        genome = "ecoli_k12"
    else:
        rows, secondary = load_bsu_tss(data_root, pwms)
        genome = "subtilis"

    seqs, seq_source = load_genome_sequence(data_root, genome)
    if rows and seqs and pwms:
        null = build_null(seqs, pwms, upstream + downstream + 1,
                          null_windows, seed)
        inventory = (_inventory_for_confidence(data_root, genome)
                     if dataset == "k12" else set())
        scored = score_tss_rows(rows, seqs, pwms, null,
                                inventory=inventory,
                                upstream=upstream, downstream=downstream)
        metrics = compute_metrics(scored, dataset, upstream=upstream,
                                  downstream=downstream)
        metrics["sequence_source"] = seq_source
        metrics["contig"] = max(seqs, key=lambda name: len(seqs[name]))
        metrics["null"] = {
            "windows": null_windows, "seed": seed,
            "n_per_sigma": {sigma: int(null[sigma]["n"]) for sigma in null},
        }
        metrics["inventory_for_confidence"] = sorted(inventory)
        if secondary is not None:
            metrics["secondary_source"] = secondary
        if dataset == "subtilis":
            metrics["mapping"] = _mapping_summary(scored)
        if dataset == "k12":
            metrics["operon_comparison"] = compare_k12_existing(data_root, scored)
        confusion = {
            "dataset": dataset,
            "labels": list(SIGMA_LABELS),
            "matrix": metrics["confusion_matrix"],
            "row_totals": {label: sum(row) for label, row
                           in zip(SIGMA_LABELS, metrics["confusion_matrix"])},
            "col_totals": {label: sum(metrics["confusion_matrix"][i][j]
                                      for i in range(len(SIGMA_LABELS)))
                           for j, label in enumerate(SIGMA_LABELS)},
            "n": metrics["counts"]["n_evaluable"],
            "note": ("filas = etiqueta curada; columnas = predicción top-1; "
                     "solo TSS evaluables"),
        }
        _write_tsv(out_dir / f"{dataset}_tss_predictions.tsv",
                   _sanitize_tsv_rows(scored), PREDICTION_COLS)
        _write_json(out_dir / f"{dataset}_metrics.json", metrics)
        _write_json(out_dir / f"{dataset}_confusion.json", confusion)
        return {"status": "evaluated", "metrics": metrics,
                "confusion": confusion, "rows": scored}

    metrics = {
        "dataset": dataset,
        "status": "not_available",
        "counts": {"n_rows": len(rows), "n_evaluable": 0},
        "reason": ("faltan datos: filas, secuencia o PWMs"
                   if not rows or not seqs or not pwms else "sin TSS evaluables"),
        "sequence_source": seq_source,
        "confusion_matrix": [[0] * len(SIGMA_LABELS) for _ in SIGMA_LABELS],
        "confusion_labels": list(SIGMA_LABELS),
    }
    if secondary is not None:
        metrics["secondary_source"] = secondary
    _write_json(out_dir / f"{dataset}_metrics.json", metrics)
    _write_json(out_dir / f"{dataset}_confusion.json", {
        "dataset": dataset, "status": "not_available",
        "labels": list(SIGMA_LABELS),
        "matrix": metrics["confusion_matrix"], "n": 0,
    })
    return {"status": "not_available", "metrics": metrics, "rows": rows}


# --------------------------------------------------------------------------
# Resumen markdown
# --------------------------------------------------------------------------
def _f(value: Any, digits: int = 3) -> str:
    if value is None:
        return "not_available"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def _pct(value: Any) -> str:
    if value is None:
        return "not_available"
    return f"{100 * float(value):.2f}%"


def _pp(value: Any) -> str:
    if value is None:
        return "not_available"
    return f"{100 * float(value):+.2f} pp"


def _table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |",
             "|" + "|".join(["---"] * len(headers)) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(_f(cell) for cell in row) + " |")
    return "\n".join(lines)


def _dataset_section(name: str, result: dict[str, Any]) -> list[str]:
    metrics = result.get("metrics", {})
    if result.get("status") != "evaluated":
        return [f"## {name}", "",
                f"**not_available**: {metrics.get('reason', 'sin datos')}", ""]
    counts = metrics["counts"]
    coverage = metrics["coverage"]
    accuracy = metrics["accuracy"]
    lines = [f"## {name}", ""]
    lines += [
        _table(
            ["métrica", "valor"],
            [
                ["TSS en la tabla", counts["n_rows"]],
                ["TSS con coordenada válida", counts["n_valid_tss"]],
                ["TSS evaluables (sigma con PWM)", counts["n_evaluable"]],
                ["TSS sin PWM (no cuentan como fallo)", counts["n_no_pwm"]],
                ["TSS sin etiqueta", counts["n_unlabeled"]],
                ["TSS sin coordenada", counts["n_missing_tss"]],
                ["Cobertura sobre etiquetados con PWM", _pct(coverage["evaluable_over_labeled_valid"])],
                ["Cobertura sobre TSS válidos", _pct(coverage["evaluable_over_valid_tss"])],
                ["Tasa no_pwm sobre etiquetados", _pct(coverage["no_pwm_over_labeled_valid"])],
                ["Accuracy top-1", _pct(accuracy["top1"])],
                ["Accuracy top-3", _pct(accuracy["top3"])],
                ["Baseline mayoría", _pct(accuracy["majority_baseline"])],
                ["Ganancia sobre mayoría", _pp(accuracy["gain_over_majority"])],
                ["Macro-F1 (clases con apoyo)", _f(metrics["macro_f1_supported"], 3)],
                ["Macro-F1 (6 clases)", _f(metrics["macro_f1_all6"], 3)],
                ["Weighted-F1", _f(metrics["weighted_f1"], 3)],
            ],
        ),
        "",
    ]
    labels = metrics["confusion_labels"]
    matrix = metrics["confusion_matrix"]
    lines += ["### Matriz de confusión (filas = verdad, columnas = top-1)", ""]
    header = ["verdad \\ pred"] + labels
    rows = [[labels[i]] + matrix[i] for i in range(len(labels))]
    lines += [_table(header, rows), ""]
    lines += ["### Por tipo de sigma", ""]
    header = ["sigma", "support", "predichos", "TP", "FP", "FN", "precision", "recall", "F1"]
    rows = []
    for sigma in labels:
        item = metrics["per_sigma"][sigma]
        rows.append([sigma, item["support"], item["predicted"], item["tp"],
                     item["fp"], item["fn"], _f(item["precision"]),
                     _f(item["recall"]), _f(item["f1"])])
    lines += [_table(header, rows), ""]
    cal = metrics.get("calibration", {})
    lines += ["### Calibración por confianza", ""]
    if cal.get("status") == "available":
        lines += [
            _table(
                ["etiqueta", "n", "correctos", "accuracy", "confianza media", "gap"],
                [[label, data["n"], data["correct"], _pct(data["accuracy"]),
                  _f(data["mean_confidence_score"]), _f(data["gap"])]
                 for label, data in cal.get("by_confidence_label", {}).items()],
            ),
            "",
            _table(
                ["bin de score", "n", "accuracy", "confianza media", "gap"],
                [[f"{item['low']:.1f}-{item['high']:.1f}", item["n"],
                  _pct(item["accuracy"]), _f(item["mean_confidence_score"]),
                  _f(item["gap"])] for item in cal.get("by_score_bin", [])],
            ),
            "",
            f"- ECE (diagnóstico): **{_f(cal.get('expected_calibration_error'))}**",
            f"- Brier (diagnóstico): **{_f(cal.get('brier_score'))}**",
            f"- {cal.get('note', '')}",
            "",
        ]
    else:
        lines += [f"**not_available**: {cal.get('reason', 'sin datos')}", ""]
    return lines


def _operon_section(comparison: dict[str, Any]) -> list[str]:
    if comparison.get("status") != "evaluated":
        return ["### Comparación con la salida existente por operón/gen", "",
                f"**not_available**: {comparison.get('reason', 'sin datos')}", ""]
    lines = ["### Comparación con la salida existente por operón/gen", ""]
    lines += [
        _table(
            ["métrica", "valor"],
            [
                ["Filas existentes", comparison["n_existing_rows"]],
                ["Filas existentes mapeadas a un gen del GT", comparison["n_existing_matched"]],
                ["Genes del GT con sigma mapeable", comparison["n_gt_genes_with_mappable"]],
                ["Genes del GT mapeados a la salida existente", comparison["n_gt_genes_mapped_to_existing"]],
                ["Coincidencia exacta (predicho ∈ sigma del gen)", _pct(comparison["exact_rate"])],
                ["Coincidencia con el sigma dominante", _pct(comparison["dominant_rate"])],
                ["Genes con ≥2 TSS", comparison["n_gt_genes_with_multiple_tss"]],
                ["Genes con ≥2 sigma distintos", comparison["n_gt_genes_with_multiple_sigma"]],
                ["% genes con ≥2 sigma", _pct(comparison["pct_gt_genes_with_multiple_sigma"])],
            ],
        ),
        "",
        "**Limitación de un sigma por operón:** " + comparison["one_sigma_per_operon_limit"]["note"] + ".",
        "",
        _table(
            ["subconjunto", "n genes", "coincidencia exacta", "tasa"],
            [
                ["un solo sigma curado", comparison["single_label"]["n_genes"],
                 comparison["single_label"]["exact"],
                 _pct(comparison["single_label"]["exact_rate"])],
                ["≥2 sigma curados (coincide con alguno)", comparison["multi_label"]["n_genes"],
                 comparison["multi_label"]["any_exact"],
                 _pct(comparison["multi_label"]["any_exact_rate"])],
                ["≥2 sigma curados (coincide con el dominante)", comparison["multi_label"]["n_genes"],
                 comparison["multi_label"]["dominant_exact"],
                 _pct(comparison["multi_label"]["dominant_exact_rate"])],
            ],
        ),
        "",
        f"- Genes con ≥2 sigma donde la predicción coincide con **todos**: "
        f"**{comparison['one_sigma_per_operon_limit']['n_multi_sigma_genes_all_labels_match']}** "
        f"(imposible por construcción: una fila = un sigma).",
        f"- Método de mapeo: {comparison['mapping_method']}.",
        "",
    ]
    return lines


def _bsu_mapping_section(metrics: dict[str, Any]) -> list[str]:
    mapping = metrics.get("mapping")
    if not mapping:
        return []
    lines = ["### Mapeo de etiquetas de *B. subtilis*", ""]
    rows = []
    for raw, item in mapping["by_raw_label"].items():
        rows.append([raw or "<vacío>", item["n"], item["mapping_status"],
                     item["sigma_truth"] or "-", item["reason"] or "-"])
    lines += [_table(["etiqueta", "n", "estado", "sigma usado", "motivo"], rows), ""]
    return lines


def render_summary(results: dict[str, Any], config: dict[str, Any]) -> str:
    lines = [
        "# Benchmark sigma TSS-level — K-12 y *Bacillus subtilis*",
        "",
        "Generado por `scripts/sigma/30_validate_sigma.py`.",
        "",
        "## Método",
        "",
        f"- Ventana: **-{config['upstream']} / +{config['downstream']}** alrededor del TSS "
        f"({config['upstream'] + config['downstream'] + 1} nt), en orientación de transcripción.",
        "- Se puntúan las seis PWMs (`sigma70`, `sigma54`, `sigma32`, `sigma38`, `sigma28`, `sigma24`) "
        "con `score_sigma`; top-1 por z-score.",
        f"- Fondo nulo: **{config['null_windows']}** ventanas genómicas aleatorias, "
        f"semilla **{config['seed']}**, construido con `build_null`.",
        "- `top-3`: la etiqueta curada está entre los tres sigma con mayor z.",
        "- `no_pwm`: etiqueta sigma sin PWM equivalente; se excluye del denominador y "
        "**no se cuenta como fallo**.",
        "",
    ]
    lines += _dataset_section("K-12 (RegulonDB/GT)", results.get("k12", {}))
    if results.get("k12"):
        lines += _operon_section(results["k12"]["metrics"].get("operon_comparison", {}))
    lines += _dataset_section("*B. subtilis* (BSGatlas/DBTBS)", results.get("subtilis", {}))
    if results.get("subtilis"):
        lines += _bsu_mapping_section(results["subtilis"]["metrics"])
    lines += [
        "## Límites",
        "",
        "- **TSS vs operón:** este benchmark evalúa TSS aislados; la salida del módulo "
        "asigna un sigma por operón. Son unidades distintas y no deben compararse 1:1.",
        "- **Etiquetas sin PWM:** sigma19 en K-12 y E/F/G/K/H/W/X/M de *B. subtilis* no "
        "tienen PWM equivalente; se reportan como `no_pwm` y no como error.",
        "- **Transferencia entre especies:** *B. subtilis* se puntúa con PWMs de *E. coli* "
        "sobre etiquetas homologadas (A→sigma70, B→sigma38, D→sigma28, I→sigma32, L→sigma54). "
        "Una coincidencia no demuestra que la PWM sea la correcta para esa especie.",
        "- **Múltiples TSS:** un gen puede tener varios TSS y varios sigma; una predicción "
        "por operón no puede representarlos todos.",
        "- **Calibración:** `confidence_score` no es una probabilidad; ECE/Brier son "
        "diagnósticos de ordenamiento, no probabilidades calibradas.",
        "- **Fondo nulo:** las ventanas aleatorias incluyen promotores reales y regiones "
        "no promotoras; el z-score es una referencia interna, no un p-valor biológico.",
        "",
        "## Conclusión",
        "",
    ]
    k12 = results.get("k12", {})
    bsu = results.get("subtilis", {})
    if k12.get("status") == "evaluated":
        m = k12["metrics"]
        lines.append(
            f"- K-12: top-1 **{_pct(m['accuracy']['top1'])}**, top-3 "
            f"**{_pct(m['accuracy']['top3'])}** sobre **{m['counts']['n_evaluable']}** TSS "
            f"evaluables; baseline de mayoría {_pct(m['accuracy']['majority_baseline'])}.")
    else:
        lines.append("- K-12: **not_available**.")
    if bsu.get("status") == "evaluated":
        m = bsu["metrics"]
        lines.append(
            f"- *B. subtilis*: top-1 **{_pct(m['accuracy']['top1'])}**, top-3 "
            f"**{_pct(m['accuracy']['top3'])}** sobre **{m['counts']['n_evaluable']}** TSS "
            f"evaluables; baseline de mayoría {_pct(m['accuracy']['majority_baseline'])}.")
    else:
        lines.append("- *B. subtilis*: **not_available**.")
    lines += [
        "- La comparación por operón muestra la limitación estructural, no solo un "
        "desacuerdo de etiquetas: la salida existente no puede expresar varios sigma "
        "para el mismo gen.",
        "",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Autoprueba
# --------------------------------------------------------------------------
def selftest() -> str:
    """Prueba sin datos pesados: mapeo, ventana, métricas y calibración."""
    assert map_k12_sigma("sigma70")["mapping_status"] == "evaluable"
    assert map_k12_sigma("sigma19")["mapping_status"] == "no_pwm"
    assert map_k12_sigma("")["mapping_status"] == "unlabeled"
    assert map_bsu_sigma("A")["sigma_truth"] == "sigma70"
    assert map_bsu_sigma("SigB")["sigma_truth"] == "sigma38"
    assert map_bsu_sigma("D")["sigma_truth"] == "sigma28"
    assert map_bsu_sigma("E")["mapping_status"] == "no_pwm"
    assert map_bsu_sigma("?")["mapping_status"] == "unlabeled"

    seq = "ACGT" * 60
    plus = tss_window(seq, 121, "+")
    minus = tss_window(seq, 121, "-")
    assert len(plus) == TSS_WINDOW and len(minus) == TSS_WINDOW
    assert plus[80] == seq[120]
    assert minus == revcomp(seq[100:201])

    rows = [
        {"label_status": "evaluable", "sigma_truth": "sigma70", "pred_sigma": "sigma70",
         "correct_top1": 1, "correct_top3": 1, "pred_rank": 1, "pred_z": 4.0,
         "confidence_label": "sospecha_fuerte", "confidence_score": 0.9},
        {"label_status": "evaluable", "sigma_truth": "sigma70", "pred_sigma": "sigma38",
         "correct_top1": 0, "correct_top3": 0, "pred_rank": 5, "pred_z": 1.0,
         "confidence_label": "sin_sospecha", "confidence_score": 0.2},
        {"label_status": "evaluable", "sigma_truth": "sigma38", "pred_sigma": "sigma70",
         "correct_top1": 0, "correct_top3": 1, "pred_rank": 2, "pred_z": 3.0,
         "confidence_label": "sospecha_media", "confidence_score": 0.5},
        {"label_status": "no_pwm", "sigma_truth": "", "pred_sigma": "",
         "mapping_status": "no_pwm"},
    ]
    metrics = compute_metrics(rows, "selftest")
    assert metrics["counts"]["n_evaluable"] == 3
    assert metrics["counts"]["n_no_pwm"] == 1
    assert abs(metrics["accuracy"]["top1"] - 1 / 3) < 1e-6
    assert abs(metrics["accuracy"]["top3"] - 2 / 3) < 1e-6
    labels = list(SIGMA_LABELS)
    assert metrics["confusion_matrix"][0][0] == 1
    assert metrics["confusion_matrix"][0][labels.index("sigma38")] == 1
    assert metrics["per_sigma"]["sigma70"]["tp"] == 1
    cal = metrics["calibration"]
    assert cal["status"] == "available" and cal["n"] == 3
    assert "OK" in scanner_selftest()
    return "mapeo, ventana, métricas y calibración sigma OK"


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark TSS-level del módulo sigma (K-12 y B. subtilis).")
    parser.add_argument("--k12", action="store_true", help="benchmark K-12")
    parser.add_argument("--subtilis", action="store_true", help="benchmark B. subtilis")
    parser.add_argument("--all", action="store_true", help="ambos genomas")
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "sigma" / "benchmark",
                        help="directorio de salida (default: data/sigma/benchmark)")
    parser.add_argument("--selftest", action="store_true", help="prueba rápida sin datos")
    parser.add_argument("--seed", type=int, default=BENCH_SEED)
    parser.add_argument("--null-windows", type=int, default=BENCH_NULL_WINDOWS)
    parser.add_argument("--upstream", type=int, default=TSS_UPSTREAM)
    parser.add_argument("--downstream", type=int, default=TSS_DOWNSTREAM)
    args = parser.parse_args(argv)

    if args.selftest:
        print(selftest())
        return 0
    if not (args.k12 or args.subtilis or args.all):
        parser.error("indicá --k12, --subtilis, --all o --selftest")

    data_root = ROOT / "data"
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    pwms = load_pwm_library(data_root)
    if not pwms:
        print("ERROR: no hay PWMs en data/regulondb/gt/promoter_pwms.tsv",
              file=sys.stderr)
        return 1

    config = {
        "seed": args.seed,
        "null_windows": args.null_windows,
        "upstream": args.upstream,
        "downstream": args.downstream,
    }
    started = time.time()
    results: dict[str, Any] = {}
    if args.all or args.k12:
        results["k12"] = _run_dataset("k12", data_root, out_dir, pwms, **config)
    if args.all or args.subtilis:
        results["subtilis"] = _run_dataset("subtilis", data_root, out_dir, pwms, **config)

    summary = render_summary(results, config)
    (out_dir / "benchmark_summary.md").write_text(summary + "\n", encoding="utf-8")

    for name, result in results.items():
        if result["status"] != "evaluated":
            print(f"{name}: not_available")
            continue
        metrics = result["metrics"]
        print(f"{name}: n_evaluable={metrics['counts']['n_evaluable']} "
              f"top1={_pct(metrics['accuracy']['top1'])} "
              f"top3={_pct(metrics['accuracy']['top3'])} "
              f"no_pwm={metrics['counts']['n_no_pwm']}")
    print(f"salida: {out_dir} ({time.time() - started:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
