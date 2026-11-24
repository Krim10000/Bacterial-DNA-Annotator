#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Benchmark reproducible por capas usando SOLO datos locales.

Capas
-----
- ``abc``: conteos de las rutas A/B/C de los 400 genomas y cobertura
  descriptiva. No hay gold local de dominios/función, así que la precisión de
  función se declara ``not_evaluable``.
- ``membrane``: composición y validación estructural de los TSV/GFF de
  membrana. La validación local comprueba formato/coordenadas/1:1, no verdad
  biológica; la precisión de localización se declara ``not_evaluable``.
- ``complexes``: cobertura/concordancia de los clusters K-12 contra el gold
  ComplexPortal K-12. Es un benchmark **interno y no ciego**: el modelo de
  pares se ajustó con el mismo gold. No es precisión predictiva.
- ``oric``: estados de oriC y distancia a DoriC. ``validated`` usa DoriC como
  ancla, por lo que se separa de ``consensus`` y no se presenta como exactitud
  independiente.

Uso
---
    python tools/benchmark_layers.py
    python tools/benchmark_layers.py abc
    python tools/benchmark_layers.py --md reports/layer_benchmark.md \
        --json reports/layer_benchmark.json

No hay red, no se reentrena nada y no se modifican pipeline, GUI ni docs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MD = ROOT / "reports" / "layer_benchmark.md"
DEFAULT_JSON = ROOT / "reports" / "layer_benchmark.json"

ALL_LAYERS = ("abc", "membrane", "complexes", "oric")
ABC_LABELS = (
    "A_named_function",
    "B_uncharacterized_homolog",
    "C_no_swissprot_hit",
)
ABC_FRIENDLY = {
    "A_named_function": "A · función nombrada por homología",
    "B_uncharacterized_homolog": "B · homólogo no caracterizado / DUF",
    "C_no_swissprot_hit": "C · sin hit Swiss-Prot",
}
STRENGTH_ORDER = ("strong", "moderate", "weak", "none", "unknown")
THRESHOLDS_BP = (1000, 10000, 50000)
STATUS_LABEL = {
    "partial": "parcial",
    "internal_no_blind": "interno (no ciego)",
    "not_evaluable": "no evaluable",
}


# ---------------------------------------------------------------------------
# Utilidades de E/S y métricas puras
# ---------------------------------------------------------------------------
def iter_tsv(path: Path) -> Iterable[dict[str, str]]:
    """Itera un TSV con cabecera, sin cargarlo entero en memoria."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def read_tsv(path: Path) -> list[dict[str, str]]:
    return list(iter_tsv(path))


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _rel(root: Path, path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return Path(path).as_posix()


def _int(value: Any) -> int | None:
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        try:
            return int(float(text))
        except ValueError:
            return None


def _pct(numerator: int | float, denominator: int | float) -> float | None:
    if not denominator:
        return None
    return round(100.0 * float(numerator) / float(denominator), 2)


def _ratio(numerator: int | float, denominator: int | float) -> float | None:
    if not denominator:
        return None
    return round(float(numerator) / float(denominator), 4)


def _round(value: float | None, ndigits: int = 4) -> float | None:
    return None if value is None else round(float(value), ndigits)


def quantile(values: Sequence[float], q: float) -> float | None:
    """Cuantil lineal (interpolación entre posiciones) sin dependencias."""
    vals = sorted(float(v) for v in values)
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return vals[lo]
    frac = pos - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac


def jaccard(a: set[str], b: set[str]) -> float:
    """Jaccard de dos conjuntos; 0.0 si ambos están vacíos."""
    if not a and not b:
        return 0.0
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def set_coverage(predicted: set[str], gold: set[str]) -> dict[str, Any]:
    """Cobertura de gold y fracción de la predicción dentro de gold.

    Se evita llamar ``precision`` a la fracción de la predicción: en un mundo
    abierto (complejos, función) un no solapado no es necesariamente falso.
    """
    shared = predicted & gold
    return {
        "n_pred": len(predicted),
        "n_gold": len(gold),
        "n_shared": len(shared),
        "coverage_of_gold": _ratio(len(shared), len(gold)),
        "fraction_of_pred_in_gold": _ratio(len(shared), len(predicted)),
        "jaccard": _round(jaccard(predicted, gold), 4),
    }


def best_match(query: set[str], references: Mapping[str, set[str]]) -> dict[str, Any]:
    """Mejor referencia por Jaccard; desempate lexicográfico por id."""
    best_id: str | None = None
    best_j = 0.0
    best_shared = 0
    for ref_id in sorted(references):
        ref = references[ref_id]
        score = jaccard(query, ref)
        if score > best_j or (score == best_j and best_id is None):
            best_id = ref_id
            best_j = score
            best_shared = len(query & ref)
    return {
        "id": best_id,
        "jaccard": _round(best_j, 4),
        "n_shared": best_shared,
    }


def greedy_one_to_one(
    preds: Mapping[str, set[str]],
    golds: Mapping[str, set[str]],
    threshold: float = 0.5,
) -> list[dict[str, Any]]:
    """Emparejamiento 1:1 voraz por Jaccard descendente.

    Es determinista: ante empate ordena por ``(pred_id, gold_id)``. No pretende
    ser el óptimo global; se documenta para que las cifras sean reproducibles.
    """
    pairs: list[tuple[float, str, str]] = []
    for pred_id, pred_set in preds.items():
        for gold_id, gold_set in golds.items():
            score = jaccard(pred_set, gold_set)
            if score >= threshold and score > 0.0:
                pairs.append((score, pred_id, gold_id))
    pairs.sort(key=lambda item: (-item[0], item[1], item[2]))
    used_pred: set[str] = set()
    used_gold: set[str] = set()
    matches: list[dict[str, Any]] = []
    for score, pred_id, gold_id in pairs:
        if pred_id in used_pred or gold_id in used_gold:
            continue
        used_pred.add(pred_id)
        used_gold.add(gold_id)
        matches.append({
            "pred_id": pred_id,
            "gold_id": gold_id,
            "jaccard": _round(score, 4),
        })
    return matches


def circular_distance(position: int, reference: int, genome_length: int) -> int:
    """Distancia circular en pb entre dos coordenadas 1-based.

    Es la misma convención de ``scripts/oric/03_orca_predict.py``:
    ``min(|pos-ref|, L-|pos-ref|)``. DoriC puede tener ``start > end`` cuando
    el oriC cruza el origen; se usa ``start`` como ancla, igual que el script.
    """
    length = int(genome_length)
    if length <= 0:
        raise ValueError("genome_length debe ser > 0")
    a = int(position) % length
    b = int(reference) % length
    delta = abs(a - b)
    return min(delta, length - delta)


def distance_summary(
    distances: Iterable[int | None],
    thresholds: Sequence[int] = THRESHOLDS_BP,
) -> dict[str, Any]:
    """Resumen determinista de distancias (pb). Ignora ``None``."""
    vals = sorted(int(d) for d in distances if d is not None)
    if not vals:
        return {
            "n": 0,
            "min_bp": None,
            "median_bp": None,
            "mean_bp": None,
            "p90_bp": None,
            "max_bp": None,
            **{f"within_{t}bp": {"n": 0, "pct": None} for t in thresholds},
            "beyond_max_threshold": {"n": 0, "pct": None},
        }
    out: dict[str, Any] = {
        "n": len(vals),
        "min_bp": vals[0],
        "median_bp": _round(statistics.median(vals), 1),
        "mean_bp": _round(statistics.mean(vals), 1),
        "p90_bp": _round(quantile(vals, 0.9), 1),
        "max_bp": vals[-1],
    }
    for t in thresholds:
        n = sum(1 for d in vals if d <= t)
        out[f"within_{t}bp"] = {"n": n, "pct": _pct(n, len(vals))}
    n_beyond = sum(1 for d in vals if d > max(thresholds))
    out["beyond_max_threshold"] = {"n": n_beyond, "pct": _pct(n_beyond, len(vals))}
    return out


def parse_sequence_regions(gff_path: Path) -> dict[str, int]:
    """Lee ``##sequence-region <contig> 1 <len>`` de un GFF3."""
    regions: dict[str, int] = {}
    with Path(gff_path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("##sequence-region"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                regions[parts[1]] = int(parts[3])
            except ValueError:
                continue
    return regions


def _missing_section(layer: str, reason: str, sources: Sequence[str] | None = None) -> dict:
    return {
        "layer": layer,
        "status": "not_evaluable",
        "metrics": {},
        "sources": list(sources or []),
        "not_evaluable": [{"metric": "toda la capa", "reason": reason}],
        "limitations": [reason],
    }


def _union_sets(values: Iterable[set[str]]) -> set[str]:
    out: set[str] = set()
    for value in values:
        out |= value
    return out


# ---------------------------------------------------------------------------
# Capa A/B/C
# ---------------------------------------------------------------------------
def summarize_abc(root: Path = ROOT) -> dict:
    root = Path(root)
    annot = root / "data" / "annotations_ab_initio"
    files = sorted(annot.glob("*/protein_paths.tsv"))
    if not files:
        return _missing_section(
            "A/B/C",
            f"No hay protein_paths.tsv en {annot}",
            [_rel(root, annot)],
        )

    label_counts: Counter[str] = Counter()
    strength_counts: Counter[str] = Counter()
    crosstab: dict[str, Counter[str]] = defaultdict(Counter)
    total = 0
    with_uniprot = 0
    with_gene = 0
    subunit_named = 0
    per_genome: list[dict[str, Any]] = []

    for path in files:
        n_genome = 0
        local: Counter[str] = Counter()
        for row in iter_tsv(path):
            n_genome += 1
            label = (row.get("label") or "").strip()
            label_counts[label] += 1
            local[label] += 1
            strength = (row.get("evidence_strength") or "").strip() or "unknown"
            strength_counts[strength] += 1
            crosstab[label][strength] += 1
            if (row.get("uniprot") or "").strip():
                with_uniprot += 1
            if (row.get("gene") or "").strip():
                with_gene += 1
            if (row.get("subunit_hint") or "").strip() == "named_subunit":
                subunit_named += 1
        total += n_genome
        per_genome.append({
            "genome": path.parent.name,
            "proteins": n_genome,
            **{label: local.get(label, 0) for label in ABC_LABELS},
            "other": n_genome - sum(local.get(label, 0) for label in ABC_LABELS),
        })

    order = {label: i for i, label in enumerate(ABC_LABELS)}
    label_order = sorted(label_counts, key=lambda x: (order.get(x, 99), x))
    labels_payload = {
        label: {"n": label_counts[label], "pct": _pct(label_counts[label], total)}
        for label in label_order
    }
    strength_payload = {
        s: {"n": strength_counts[s], "pct": _pct(strength_counts[s], total)}
        for s in sorted(strength_counts, key=lambda x: (STRENGTH_ORDER.index(x) if x in STRENGTH_ORDER else 99, x))
    }
    label_x_evidence = {
        label: {s: crosstab[label].get(s, 0) for s in STRENGTH_ORDER}
        for label in label_order
    }

    proteins = [g["proteins"] for g in per_genome if g["proteins"]]
    a_pcts = [
        100.0 * g.get("A_named_function", 0) / g["proteins"]
        for g in per_genome if g["proteins"]
    ]

    batch_path = annot / "batch_summary.json"
    expected_cds: int | None = None
    if batch_path.exists():
        try:
            payload = read_json(batch_path)
            raw_expected = (payload.get("totals") or {}).get("cds")
            expected_cds = int(raw_expected) if raw_expected is not None else None
        except (ValueError, TypeError, json.JSONDecodeError):
            expected_cds = None

    metrics = {
        "genomes": len(files),
        "proteins": total,
        "labels": labels_payload,
        "evidence_strength": strength_payload,
        "label_x_evidence": label_x_evidence,
        "with_uniprot": {"n": with_uniprot, "pct": _pct(with_uniprot, total)},
        "with_gene": {"n": with_gene, "pct": _pct(with_gene, total)},
        "subunit_hint_named": {"n": subunit_named, "pct": _pct(subunit_named, total)},
        "per_genome_spread": {
            "proteins_min": min(proteins) if proteins else None,
            "proteins_median": _round(statistics.median(proteins), 1) if proteins else None,
            "proteins_max": max(proteins) if proteins else None,
            "a_pct_min": _round(min(a_pcts), 2) if a_pcts else None,
            "a_pct_median": _round(statistics.median(a_pcts), 2) if a_pcts else None,
            "a_pct_max": _round(max(a_pcts), 2) if a_pcts else None,
        },
        "batch_summary": {
            "path": _rel(root, batch_path) if batch_path.exists() else None,
            "expected_cds": expected_cds,
            "matches_protein_rows": (total == expected_cds) if expected_cds is not None else None,
        },
        "gold_domains_or_function": {
            "found": False,
            "checked": [
                _rel(root, root / "data" / "complexes" / "gold_k12.tsv"),
                _rel(root, root / "data" / "doric_ground_truth.json"),
                _rel(root, root / "data" / "benchmark" / "benchmark_summary.json"),
                _rel(root, root / "data" / "annotations"),
            ],
            "note": (
                "Los archivos encontrados son gold de complejos, de oriC, de "
                "llamada de genes o anotación NCBI de estructura. Ninguno se usa "
                "como gold de dominios/función para A/B/C: los productos NCBI no "
                "son una verdad de dominios y compararlos exigiría un mapeo y "
                "criterios que este benchmark no inventa."
            ),
        },
    }
    return {
        "layer": "A/B/C",
        "status": "partial",
        "metrics": metrics,
        "sources": [
            _rel(root, annot),
            _rel(root, batch_path) if batch_path.exists() else "",
        ],
        "not_evaluable": [
            {
                "metric": "precisión/recall de función nombrada (A) contra gold",
                "reason": (
                    "No existe una tabla local de gold de función para los 400 "
                    "genomas. 'A' significa hit Swiss-Prot con nombre, no función "
                    "confirmada experimentalmente."
                ),
            },
            {
                "metric": "precisión/recall de dominios (B/C) contra gold",
                "reason": (
                    "No hay gold local de Pfam/InterPro para estos genomas. El "
                    "panel Pfam de membrana es un panel de búsqueda, no verdad de "
                    "referencia."
                ),
            },
        ],
        "limitations": [
            "A/B/C es una ruta por homología, no una llamada de función confirmada.",
            "'C' significa sin hit Swiss-Prot; no significa que la proteína no tenga función.",
            "Los conteos dependen de la versión de Swiss-Prot/DIAMOND usada por el pipeline.",
            "No se evalúa la identidad de la función (solo si hay un nombre asignado).",
            "Existe anotación NCBI local (data/annotations) con product, pero se usa para estructura/llamada de genes; no se trata como gold de dominios/función.",
        ],
    }


# ---------------------------------------------------------------------------
# Capa membrana
# ---------------------------------------------------------------------------
TYPE_ORDER = ("tm", "signal_peptide", "lipoprotein", "transporter", "secretion")
CONF_ORDER = ("high", "medium", "low")


def summarize_membrane(root: Path = ROOT) -> dict:
    root = Path(root)
    full = root / "data" / "membrane" / "full"
    files = sorted(p for p in full.glob("*.tsv") if p.name.startswith("GCF_"))
    if not files:
        return _missing_section(
            "Membrana",
            f"No hay TSV de membrana en {full}",
            [_rel(root, full)],
        )

    types: Counter[str] = Counter()
    confidence: Counter[str] = Counter()
    rows = 0
    genomes: set[str] = set()
    proteins: set[str] = set()
    transporters = 0
    transporters_with_substrate = 0

    for path in files:
        for row in iter_tsv(path):
            rows += 1
            genome = (row.get("genome") or path.stem).strip()
            genomes.add(genome)
            proteins.add(f"{genome}|{(row.get('protein_id') or '').strip()}")
            rtype = (row.get("type") or "").strip()
            types[rtype] += 1
            confidence[(row.get("confidence") or "").strip()] += 1
            if rtype == "transporter":
                transporters += 1
                if (row.get("substrate") or "").strip():
                    transporters_with_substrate += 1

    validation_path = root / "data" / "membrane" / "validation_full.json"
    validation: dict[str, Any] | None = None
    if validation_path.exists():
        try:
            payload = read_json(validation_path)
            totals = payload.get("totals") or {}
            validation = {
                "source": _rel(root, validation_path),
                "ok": payload.get("ok"),
                "processed_ok": payload.get("processed_ok"),
                "checked_at": payload.get("checked_at"),
                "genomes": totals.get("genomes"),
                "rows": totals.get("rows"),
                "gff_features": totals.get("gff_features"),
                "errors": totals.get("errors"),
                "pending": totals.get("pending"),
                "genomes_ok": totals.get("genomes_ok"),
                "consistent_with_tables": bool(
                    totals.get("rows") == rows
                    and totals.get("gff_features") == rows
                    and totals.get("genomes") == len(genomes)
                ),
            }
        except (ValueError, TypeError, json.JSONDecodeError):
            validation = None

    type_payload = {
        t: {"n": types[t], "pct": _pct(types[t], rows)}
        for t in sorted(types, key=lambda x: (TYPE_ORDER.index(x) if x in TYPE_ORDER else 99, x))
    }
    conf_payload = {
        c: {"n": confidence[c], "pct": _pct(confidence[c], rows)}
        for c in sorted(confidence, key=lambda x: (CONF_ORDER.index(x) if x in CONF_ORDER else 99, x))
    }
    high_medium = confidence.get("high", 0) + confidence.get("medium", 0)

    not_evaluable = [
        {
            "metric": "precisión de localización (TM/SP/lipoproteína/secreción)",
            "reason": (
                "La validación local (05_validate_outputs.py) comprueba cabecera, "
                "coordenadas, protein_id y correspondencia TSV-GFF; no hay gold "
                "biológico local (TMHMM/SignalP/LipoP/UniProt) para medir aciertos."
            ),
        }
    ]
    if validation is None:
        not_evaluable.append({
            "metric": "validación estructural",
            "reason": "No se encontró data/membrane/validation_full.json.",
        })

    return {
        "layer": "Membrana",
        "status": "partial",
        "metrics": {
            "genomes": len(genomes),
            "rows": rows,
            "unique_proteins": len(proteins),
            "types": type_payload,
            "confidence": conf_payload,
            "high_medium": {"n": high_medium, "pct": _pct(high_medium, rows)},
            "transporters": {
                "n": transporters,
                "with_substrate": transporters_with_substrate,
                "pct_with_substrate": _pct(transporters_with_substrate, transporters),
            },
            "validation": validation,
        },
        "sources": [
            _rel(root, full),
            _rel(root, root / "data" / "membrane" / "gff" / "full"),
            _rel(root, root / "scripts" / "membrane" / "05_validate_outputs.py"),
            _rel(root, validation_path) if validation_path.exists() else "",
        ],
        "not_evaluable": not_evaluable,
        "limitations": [
            "TM/SP/lipobox son heurísticos; no se usaron TMHMM/SignalP/LipoP.",
            "El panel Pfam es de 255 modelos; algunas familias están subrepresentadas.",
            f"{_fmt_pct(_pct(confidence.get('low', 0), rows))} de las filas son confianza "
            "'low'; 'low' es hipótesis, no localización confirmada.",
            "La validación estructural no demuestra que la asignación sea correcta.",
        ],
    }


# ---------------------------------------------------------------------------
# Capa complejos
# ---------------------------------------------------------------------------
def _load_gold(path: Path) -> dict[str, dict[str, Any]]:
    gold: dict[str, dict[str, Any]] = {}
    for row in read_tsv(path):
        cid = (row.get("complex_id") or "").strip()
        genes = {g.strip() for g in (row.get("locus_tags") or "").split(",") if g.strip()}
        if not cid or len(genes) < 2:
            continue
        if cid in gold:
            gold[cid]["genes"] |= genes
        else:
            gold[cid] = {
                "id": cid,
                "name": (row.get("complex_name") or "").strip(),
                "genes": genes,
                "evidence": (row.get("evidence_class") or "otra").strip(),
            }
    for item in gold.values():
        item["n_genes"] = len(item["genes"])
    return gold


def _load_clusters(path: Path) -> dict[str, dict[str, Any]]:
    clusters: dict[str, dict[str, Any]] = {}
    for row in read_tsv(path):
        cid = (row.get("cluster_id") or "").strip()
        genes = {g.strip() for g in (row.get("genes") or "").split(",") if g.strip()}
        if not cid or len(genes) < 2:
            continue
        clusters[cid] = {
            "id": cid,
            "genes": genes,
            "n_genes": len(genes),
            "confidence": (row.get("confidence") or "?").strip(),
            "evidencia": (row.get("evidencia") or "?").strip(),
        }
    return clusters


def _size_bin(n: int) -> str:
    if n == 2:
        return "2"
    if n <= 4:
        return "3-4"
    return "5+"


def summarize_complexes(root: Path = ROOT) -> dict:
    root = Path(root)
    gold_path = root / "data" / "complexes" / "gold_k12.tsv"
    clusters_path = root / "data" / "complexes" / "clusters_ecoli_k12_evidencia.tsv"
    if not gold_path.exists() or not clusters_path.exists():
        missing = [str(p) for p in (gold_path, clusters_path) if not p.exists()]
        return _missing_section(
            "Complejos",
            f"Faltan archivos de complejos: {missing}",
            [_rel(root, gold_path), _rel(root, clusters_path)],
        )

    gold = _load_gold(gold_path)
    clusters = _load_clusters(clusters_path)
    if not gold or not clusters:
        return _missing_section(
            "Complejos",
            "El gold o los clusters quedaron vacíos tras filtrar >=2 genes.",
            [_rel(root, gold_path), _rel(root, clusters_path)],
        )

    gold_sets = {cid: item["genes"] for cid, item in gold.items()}
    cluster_sets = {cid: item["genes"] for cid, item in clusters.items()}
    gold_genes = _union_sets(gold_sets.values())
    cluster_genes = _union_sets(cluster_sets.values())
    best_for_gold = {cid: best_match(genes, cluster_sets) for cid, genes in gold_sets.items()}
    best_for_cluster = {cid: best_match(genes, gold_sets) for cid, genes in cluster_sets.items()}
    high_total = sum(1 for item in clusters.values() if item["confidence"] == "alta")

    def threshold_block(threshold: float) -> dict[str, Any]:
        matches = greedy_one_to_one(cluster_sets, gold_sets, threshold=threshold)
        matched_gold = {m["gold_id"] for m in matches}
        matched_clusters = {m["pred_id"] for m in matches}
        covered_genes = _union_sets(gold[gid]["genes"] for gid in matched_gold)
        matched_cluster_genes = _union_sets(
            clusters[cid]["genes"] for cid in matched_clusters
        )
        any_covered = sum(
            1 for cid in gold if best_for_gold[cid]["jaccard"] >= threshold
        )
        clusters_ge = sum(
            1 for cid in clusters if best_for_cluster[cid]["jaccard"] >= threshold
        )
        high_ge = sum(
            1 for cid, item in clusters.items()
            if item["confidence"] == "alta"
            and best_for_cluster[cid]["jaccard"] >= threshold
        )

        by_evidence: dict[str, dict[str, Any]] = {}
        for evidence in sorted({item["evidence"] for item in gold.values()}):
            ids = [cid for cid, item in gold.items() if item["evidence"] == evidence]
            covered = sum(1 for cid in ids if best_for_gold[cid]["jaccard"] >= threshold)
            by_evidence[evidence] = {
                "n": len(ids),
                "covered_any": covered,
                "pct": _pct(covered, len(ids)),
            }

        by_size: dict[str, dict[str, Any]] = {}
        for bin_name in ("2", "3-4", "5+"):
            ids = [
                cid for cid, item in gold.items()
                if _size_bin(item["n_genes"]) == bin_name
            ]
            if not ids:
                continue
            covered = sum(1 for cid in ids if best_for_gold[cid]["jaccard"] >= threshold)
            by_size[bin_name] = {
                "n": len(ids),
                "covered_any": covered,
                "pct": _pct(covered, len(ids)),
            }

        return {
            "any_cluster_covers_gold": {
                "n": any_covered,
                "pct": _pct(any_covered, len(gold)),
            },
            "one_to_one": {
                "pairs": len(matches),
                "gold_recall": _ratio(len(matched_gold), len(gold)),
                "cluster_concordance": _ratio(len(matched_clusters), len(clusters)),
                "gene_recall": _ratio(len(covered_genes), len(gold_genes)),
                "gene_purity": _ratio(
                    len(covered_genes & matched_cluster_genes),
                    len(matched_cluster_genes),
                ),
            },
            "clusters_with_best_jaccard_ge": {
                "n": clusters_ge,
                "pct": _pct(clusters_ge, len(clusters)),
            },
            "high_confidence_with_best_jaccard_ge": {
                "n": high_ge,
                "pct_of_high": _pct(high_ge, high_total),
            },
            "gold_by_evidence": by_evidence,
            "gold_by_size": by_size,
        }

    existing_evidence = Counter(item["evidencia"] for item in clusters.values())
    existing_confidence = Counter(item["confidence"] for item in clusters.values())
    transfer_files = len(list((root / "data" / "complexes").glob("transfer_*.tsv")))
    model_files = len(list((root / "data" / "complexes").glob("model_*.complexes.tsv")))

    metrics = {
        "benchmark_type": "internal_reconstruction_no_blind",
        "blind": False,
        "gold": {
            "n": len(gold),
            "genes": len(gold_genes),
            "evidence_class": dict(Counter(item["evidence"] for item in gold.values())),
            "size_bins": dict(Counter(_size_bin(item["n_genes"]) for item in gold.values())),
        },
        "clusters": {
            "n": len(clusters),
            "genes": len(cluster_genes),
            "confidence": dict(existing_confidence),
            "existing_evidence_labels": dict(existing_evidence),
        },
        "thresholds": {
            "0.8": threshold_block(0.8),
            "0.5": threshold_block(0.5),
        },
        "gold_with_any_cluster": {
            "n": sum(1 for cid in gold if best_for_gold[cid]["jaccard"] > 0.0),
            "pct": _pct(
                sum(1 for cid in gold if best_for_gold[cid]["jaccard"] > 0.0),
                len(gold),
            ),
        },
        "artifacts": {
            "transfer_files": transfer_files,
            "model_complex_files": model_files,
            "note": (
                "Inventario de artefactos para los 400. No hay gold de "
                "transferencia; no se calcula precisión sobre ellos."
            ),
        },
        "circularity": {
            "risk": "high",
            "reason": (
                "El modelo de pares se ajustó con gold_k12.tsv y los clusters K-12 "
                "se generaron aplicándolo al mismo genoma. El solapamiento con el "
                "gold no es una validación ciega."
            ),
            "what_it_does_not_measure": (
                "Precisión predictiva en genomas no vistos ni complejos ausentes "
                "del catálogo ComplexPortal."
            ),
        },
    }
    return {
        "layer": "Complejos",
        "status": "internal_no_blind",
        "metrics": metrics,
        "sources": [
            _rel(root, gold_path),
            _rel(root, clusters_path),
            _rel(root, root / "data" / "complexes" / "clusters_ecoli_k12.tsv"),
            _rel(root, root / "scripts" / "complexes" / "annotate_evidence.py"),
            "data/complexes/transfer_*.tsv",
            "data/complexes/model_*.complexes.tsv",
        ],
        "not_evaluable": [
            {
                "metric": "precisión predictiva ciega / externa",
                "reason": (
                    "El gold K-12 se usó para ajustar el modelo y para etiquetar "
                    "los clusters; no hay un conjunto K-12 held-out ni otro genoma "
                    "con gold de complejos."
                ),
            },
            {
                "metric": "falsos positivos en clusters 'nuevo'",
                "reason": (
                    "No hay verdad para complejos fuera de ComplexPortal; un cluster "
                    "no catalogado no puede contarse como falso positivo."
                ),
            },
        ],
        "limitations": [
            "Benchmark interno y no ciego: el solapamiento con gold está favorecido por construcción.",
            "Gold solo de E. coli K-12 (ComplexPortal), 273 complejos con >=2 genes.",
            "Complejos grandes (50S, 30S, complejo I) quedan casi todos fuera de Jaccard >=0,8.",
            "No se evalúa estequiometría ni interacción física; 'documentado' es coincidencia de subunidades.",
            "Los clusters 'nuevo'/'parcial' no son necesariamente falsos: el mundo está abierto.",
        ],
    }


# ---------------------------------------------------------------------------
# Capa oriC
# ---------------------------------------------------------------------------
def _genome_lengths(root: Path, genome: str, cache: dict[str, dict[str, int]]) -> dict[str, int]:
    if genome not in cache:
        path = root / "data" / "annotations_ab_initio" / genome / "annotation.gff3"
        cache[genome] = parse_sequence_regions(path) if path.exists() else {}
    return cache[genome]


def summarize_oric(root: Path = ROOT) -> dict:
    root = Path(root)
    status_path = root / "data" / "orca_predictions_fix" / "oriC_status.tsv"
    if not status_path.exists():
        return _missing_section(
            "oriC",
            f"No hay tabla de estados en {status_path}",
            [_rel(root, status_path)],
        )

    rows = read_tsv(status_path)
    lengths_cache: dict[str, dict[str, int]] = {}
    records: list[dict[str, Any]] = []
    without_length = 0
    for row in rows:
        genome = (row.get("genome") or "").strip()
        regions = _genome_lengths(root, genome, lengths_cache)
        length = max(regions.values()) if regions else None
        if length is None:
            without_length += 1
        position = _int(row.get("position"))
        doric = _int(row.get("doric"))
        distance = (
            circular_distance(position, doric, length)
            if position is not None and doric is not None and length is not None
            else None
        )
        records.append({
            "genome": genome,
            "status": (row.get("status") or "?").strip(),
            "validated_by": (row.get("validated_by") or "?").strip(),
            "position": position,
            "doric": doric,
            "length": length,
            "distance": distance,
        })

    status_counts = Counter(r["status"] for r in records)
    validated_by_counts = Counter(r["validated_by"] for r in records)
    doric_available = sum(1 for r in records if r["doric"] is not None)
    validated = [r for r in records if r["status"] == "validated"]
    anomalous = [r for r in records if r["status"] == "anomalous"]
    consensus = [r for r in records if r["status"] == "consensus"]
    lineage = [r for r in records if r["status"] == "lineage_specific"]
    validated_with_doric = [r for r in validated if r["distance"] is not None]
    anomalous_with_doric = [r for r in anomalous if r["distance"] is not None]
    all_with_doric = [r for r in records if r["distance"] is not None]
    validated_using_doric = sum(1 for r in validated if "DoriC" in r["validated_by"])
    validated_outside = [
        r for r in validated_with_doric
        if r["distance"] is not None and r["distance"] > max(THRESHOLDS_BP)
    ]

    # ORCA-only vs DoriC: otra lectura, independiente de la selección final.
    pred_dir = root / "data" / "orca_predictions"
    orca_distances: list[int] = []
    delta_checked = 0
    delta_mismatches = 0
    orca_without_length = 0
    if pred_dir.exists():
        for path in sorted(pred_dir.glob("*.json")):
            try:
                payload = read_json(path)
            except (ValueError, json.JSONDecodeError):
                continue
            if "orca_position" not in payload:
                continue
            orics = (payload.get("doric_truth") or {}).get("orics") or []
            if not orics:
                continue
            genome = (payload.get("genome") or path.stem).strip()
            regions = _genome_lengths(root, genome, lengths_cache)
            length = regions.get(payload.get("replicon")) or (
                max(regions.values()) if regions else None
            )
            if length is None:
                orca_without_length += 1
                continue
            distance = circular_distance(
                int(payload["orca_position"]),
                int(orics[0]["start"]),
                length,
            )
            orca_distances.append(distance)
            delta = payload.get("delta_bp")
            if delta is not None:
                delta_checked += 1
                if int(delta) != distance:
                    delta_mismatches += 1

    doric_dataset: dict[str, Any] | None = None
    manifest_path = root / "data" / "doric" / "DoriC14_release_manifest.json"
    if manifest_path.exists():
        try:
            manifest = read_json(manifest_path)
            doric_dataset = {
                "release": manifest.get("release"),
                "csv_file": manifest.get("csv_file"),
                "csv_rows": manifest.get("csv_rows"),
                "csv_sha256": manifest.get("csv_sha256"),
                "note": "DoriC no se redistribuye; se usa solo como referencia local.",
            }
        except (ValueError, json.JSONDecodeError):
            doric_dataset = None

    metrics = {
        "genomes_with_status": len(records),
        "doric_dataset": doric_dataset,
        "status_counts": {
            s: {"n": n, "pct": _pct(n, len(records))}
            for s, n in sorted(status_counts.items())
        },
        "validated_by_counts": {
            s: {"n": n, "pct": _pct(n, len(records))}
            for s, n in sorted(validated_by_counts.items())
        },
        "doric_available": {
            "n": doric_available,
            "pct": _pct(doric_available, len(records)),
        },
        "without_doric": {
            "n": len(records) - doric_available,
            "pct": _pct(len(records) - doric_available, len(records)),
            "status_counts": {
                s: n for s, n in sorted(Counter(
                    r["status"] for r in records if r["doric"] is None
                ).items())
            },
        },
        "validated_anchored": {
            **distance_summary(r["distance"] for r in validated_with_doric),
            "n_status": len(validated),
            "n_with_doric": len(validated_with_doric),
            "n_using_doric": validated_using_doric,
            "outside_50kb": {
                "n": len(validated_outside),
                "examples": [
                    {"genome": r["genome"], "distance_bp": r["distance"]}
                    for r in sorted(
                        validated_outside,
                        key=lambda x: x["distance"] or 0,
                        reverse=True,
                    )[:5]
                ],
            },
        },
        "anomalous_with_doric": {
            **distance_summary(r["distance"] for r in anomalous_with_doric),
            "n_status": len(anomalous),
            "note": (
                "Por definición el estado anomalous agrupa discrepancias; la "
                "distancia no es una medida de acierto del pipeline."
            ),
        },
        "all_with_doric": distance_summary(r["distance"] for r in all_with_doric),
        "orca_only_vs_doric": {
            **distance_summary(orca_distances),
            "delta_bp_checked": delta_checked,
            "delta_bp_mismatches": delta_mismatches,
            "delta_bp_recomputed_matches": (
                delta_checked > 0 and delta_mismatches == 0
            ) if delta_checked else None,
            "without_length": orca_without_length,
        },
        "consensus": {
            "n": len(consensus),
            "with_doric": sum(1 for r in consensus if r["doric"] is not None),
            "distance": "not_evaluable",
        },
        "lineage_specific": {
            "n": len(lineage),
            "with_doric": sum(1 for r in lineage if r["doric"] is not None),
            "distance": "not_evaluable",
        },
        "without_length": without_length,
    }

    not_evaluable = [
        {
            "metric": "exactitud independiente del pipeline oriC",
            "reason": (
                f"Los {len(validated)} estados validated usan DoriC "
                f"({validated_using_doric} con DoriC explícito en validated_by); "
                "no hay ningún validated sin DoriC. La distancia mide consistencia "
                "con el ancla, no acierto independiente."
            ),
        },
        {
            "metric": "distancia a verdad externa para consensus",
            "reason": (
                f"Los {len(consensus)} consensus no tienen DoriC "
                f"({sum(1 for r in consensus if r['doric'] is not None)} con DoriC). "
                "El estado combina ORCA y proximidad a dnaA, no una verdad externa."
            ),
        },
        {
            "metric": "distancia para lineage_specific y anomalous sin DoriC",
            "reason": (
                f"lineage_specific={len(lineage)}, anomalous sin DoriC="
                f"{sum(1 for r in anomalous if r['doric'] is None)}; no hay verdad "
                "externa para medir distancia."
            ),
        },
    ]
    return {
        "layer": "oriC",
        "status": "partial",
        "metrics": metrics,
        "sources": [
            _rel(root, status_path),
            _rel(root, pred_dir),
            _rel(root, root / "data" / "annotations_ab_initio"),
            _rel(root, root / "data" / "doric"),
            _rel(root, root / "scripts" / "oric" / "03_orca_predict.py"),
        ],
        "not_evaluable": not_evaluable,
        "limitations": [
            f"DoriC es una referencia curada, no una verdad universal; "
            f"{len(records) - doric_available}/{len(records)} genomas no tienen DoriC.",
            "DoriC puede cruzar el origen (start > end); se usa start como ancla, igual que 03_orca_predict.py.",
            "ORCA asume linealidad; no se rota el contig en este benchmark.",
            "La tabla de estados no trae replicón: para la distancia se usa el contig más largo como cromosoma.",
            f"El estado 'validated' no es accuracy: "
            f"{validated_using_doric}/{len(validated)} usan DoriC como ancla.",
            "El estado 'anomalous' se define por discrepancia; su distancia no mide desempeño global.",
        ],
    }


# ---------------------------------------------------------------------------
# Construcción del informe
# ---------------------------------------------------------------------------
def build_data(root: Path = ROOT, layers: Sequence[str] = ALL_LAYERS) -> dict[str, Any]:
    funcs = {
        "abc": summarize_abc,
        "membrane": summarize_membrane,
        "complexes": summarize_complexes,
        "oric": summarize_oric,
    }
    sections: dict[str, Any] = {}
    for layer in layers:
        if layer not in funcs:
            raise ValueError(f"capa desconocida: {layer}")
        sections[layer] = funcs[layer](root)
    return {
        "schema": "layer-benchmark/v1",
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "project_root": str(Path(root).resolve()),
        "scope": "Solo datos locales; sin red, sin reentrenar y sin modificar pipeline/GUI/docs.",
        "layers": list(layers),
        "sections": sections,
    }


def _fmt_int(value: Any) -> str:
    if value is None:
        return "—"
    try:
        return f"{int(value):,}".replace(",", ".")
    except (TypeError, ValueError):
        return str(value)


def _fmt_num(value: Any, ndigits: int = 1) -> str:
    if value is None:
        return "—"
    try:
        text = f"{float(value):,.{ndigits}f}"
    except (TypeError, ValueError):
        return str(value)
    return text.replace(",", "X").replace(".", ",").replace("X", ".")


def _fmt_pct(value: Any) -> str:
    if value is None:
        return "—"
    return f"{_fmt_num(value, 1)} %"


def _fmt_counts(values: Mapping[str, int] | None) -> str:
    if not values:
        return "—"
    return ", ".join(f"{key} {_fmt_int(value)}" for key, value in sorted(values.items()))


def _table(headers: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(str(h) for h in headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(cell) for cell in row) + " |")
    return "\n".join(lines)


def _render_not_evaluable(section: Mapping[str, Any]) -> list[str]:
    items = section.get("not_evaluable") or []
    if not items:
        return ["- Ninguna métrica declarada no evaluable.", ""]
    out = []
    for item in items:
        out.append(f"- **{item.get('metric', '?')}**: {item.get('reason', '')}")
    out.append("")
    return out


def _render_limitations(section: Mapping[str, Any]) -> list[str]:
    items = section.get("limitations") or []
    return [f"- {item}" for item in items] + ([""] if items else [])


def _render_abc(section: Mapping[str, Any]) -> list[str]:
    m = section.get("metrics") or {}
    lines = ["## 1. A/B/C", ""]
    lines += [
        "**Estado:** parcial. Se miden conteos y cobertura descriptiva; no hay "
        "gold local de dominios/función.",
        "",
        "### Conteos",
        "",
    ]
    labels = m.get("labels") or {}
    rows = [
        [ABC_FRIENDLY.get(label, label), _fmt_int(payload.get("n")), _fmt_pct(payload.get("pct"))]
        for label, payload in labels.items()
    ]
    lines += [
        _table(["etiqueta", "proteínas", "%"], rows),
        "",
        f"- Genomas: **{_fmt_int(m.get('genomes'))}**.",
        f"- Proteínas: **{_fmt_int(m.get('proteins'))}**.",
        f"- `batch_summary.totals.cds`: **{_fmt_int((m.get('batch_summary') or {}).get('expected_cds'))}**; "
        f"coincide con las filas: **{(m.get('batch_summary') or {}).get('matches_protein_rows')}**.",
        "",
        "### Etiqueta × fuerza de evidencia",
        "",
    ]
    by_label = m.get("label_x_evidence") or {}
    rows = []
    for label, counts in by_label.items():
        rows.append(
            [ABC_FRIENDLY.get(label, label)]
            + [_fmt_int(counts.get(s, 0)) for s in STRENGTH_ORDER]
            + [_fmt_int((labels.get(label) or {}).get("n"))]
        )
    lines += [
        _table(["etiqueta", "strong", "moderate", "weak", "none", "unknown", "total"], rows),
        "",
        "### Cobertura descriptiva",
        "",
    ]
    spread = m.get("per_genome_spread") or {}
    lines += [
        f"- Proteínas por genoma: mínimo **{_fmt_int(spread.get('proteins_min'))}**, "
        f"mediana **{_fmt_num(spread.get('proteins_median'))}**, máximo **{_fmt_int(spread.get('proteins_max'))}**.",
        f"- Porcentaje A por genoma: mínimo **{_fmt_pct(spread.get('a_pct_min'))}**, "
        f"mediana **{_fmt_pct(spread.get('a_pct_median'))}**, máximo **{_fmt_pct(spread.get('a_pct_max'))}**.",
        f"- Con UniProt: **{_fmt_int((m.get('with_uniprot') or {}).get('n'))}** "
        f"({_fmt_pct((m.get('with_uniprot') or {}).get('pct'))}).",
        f"- Con gen: **{_fmt_int((m.get('with_gene') or {}).get('n'))}** "
        f"({_fmt_pct((m.get('with_gene') or {}).get('pct'))}).",
        f"- `subunit_hint=named_subunit`: **{_fmt_int((m.get('subunit_hint_named') or {}).get('n'))}** "
        f"({_fmt_pct((m.get('subunit_hint_named') or {}).get('pct'))}).",
        "",
        "### No evaluable",
        "",
    ]
    lines += _render_not_evaluable(section)
    lines += ["### Limitaciones", ""]
    lines += _render_limitations(section)
    return lines


def _render_membrane(section: Mapping[str, Any]) -> list[str]:
    m = section.get("metrics") or {}
    lines = ["## 2. Membrana", ""]
    lines += [
        "**Estado:** parcial. La validación local es estructural; no hay gold "
        "biológico de localización.",
        "",
        "### Composición",
        "",
    ]
    rows = [
        [t, _fmt_int(payload.get("n")), _fmt_pct(payload.get("pct"))]
        for t, payload in (m.get("types") or {}).items()
    ]
    lines += [
        _table(["tipo", "filas", "%"], rows),
        "",
        f"- Genomas: **{_fmt_int(m.get('genomes'))}**; filas: **{_fmt_int(m.get('rows'))}**; "
        f"proteínas únicas: **{_fmt_int(m.get('unique_proteins'))}**.",
        "",
        "### Confianza y transportadores",
        "",
    ]
    rows = [
        [c, _fmt_int(payload.get("n")), _fmt_pct(payload.get("pct"))]
        for c, payload in (m.get("confidence") or {}).items()
    ]
    lines += [
        _table(["confianza", "filas", "%"], rows),
        "",
        f"- high+medium: **{_fmt_int((m.get('high_medium') or {}).get('n'))}** "
        f"({_fmt_pct((m.get('high_medium') or {}).get('pct'))}).",
    ]
    transporters = m.get("transporters") or {}
    lines += [
        f"- Transportadores: **{_fmt_int(transporters.get('n'))}**; con sustrato explícito "
        f"**{_fmt_int(transporters.get('with_substrate'))}** "
        f"({_fmt_pct(transporters.get('pct_with_substrate'))}).",
        "",
        "### Validación estructural local",
        "",
    ]
    validation = m.get("validation")
    if validation:
        lines += [
            f"- Fuente: `{validation.get('source')}` (checked_at `{validation.get('checked_at')}`).",
            f"- ok: **{validation.get('ok')}**; processed_ok: **{validation.get('processed_ok')}**.",
            f"- genomas: **{_fmt_int(validation.get('genomes'))}**; filas: "
            f"**{_fmt_int(validation.get('rows'))}**; features GFF: "
            f"**{_fmt_int(validation.get('gff_features'))}**.",
            f"- errores: **{_fmt_int(validation.get('errors'))}**; pendientes: "
            f"**{_fmt_int(validation.get('pending'))}**; genomas_ok: "
            f"**{_fmt_int(validation.get('genomes_ok'))}**.",
            f"- Consistente con los TSV leídos: **{validation.get('consistent_with_tables')}**.",
        ]
    else:
        lines.append("- No se encontró `validation_full.json`.")
    lines += ["", "### No evaluable", ""]
    lines += _render_not_evaluable(section)
    lines += ["### Limitaciones", ""]
    lines += _render_limitations(section)
    return lines


def _render_complexes(section: Mapping[str, Any]) -> list[str]:
    m = section.get("metrics") or {}
    gold = m.get("gold") or {}
    clusters = m.get("clusters") or {}
    lines = ["## 3. Complejos", ""]
    lines += [
        "**Estado:** benchmark interno y **no ciego**. El modelo se ajustó con el "
        "mismo gold K-12; el solapamiento no es precisión predictiva.",
        "",
        "### Conjuntos",
        "",
        _table(
            ["conjunto", "n", "genes", "detalle"],
            [
                ["gold K-12", _fmt_int(gold.get("n")), _fmt_int(gold.get("genes")),
                 f"evidencia: {_fmt_counts(gold.get('evidence_class'))}"],
                ["clusters K-12", _fmt_int(clusters.get("n")), _fmt_int(clusters.get("genes")),
                 f"confianza: {_fmt_counts(clusters.get('confidence'))}"],
            ],
        ),
        "",
        "### Concordancia con gold por umbral de Jaccard",
        "",
    ]
    rows = []
    for threshold in ("0.8", "0.5"):
        block = (m.get("thresholds") or {}).get(threshold) or {}
        one = block.get("one_to_one") or {}
        any_block = block.get("any_cluster_covers_gold") or {}
        rows.append([
            threshold.replace(".", ","),
            f"{_fmt_int(any_block.get('n'))} ({_fmt_pct(any_block.get('pct'))})",
            _fmt_int(one.get("pairs")),
            _fmt_pct((one.get("gold_recall") or 0) * 100) if one.get("gold_recall") is not None else "—",
            _fmt_pct((one.get("cluster_concordance") or 0) * 100) if one.get("cluster_concordance") is not None else "—",
            _fmt_pct((one.get("gene_recall") or 0) * 100) if one.get("gene_recall") is not None else "—",
            _fmt_pct((one.get("gene_purity") or 0) * 100) if one.get("gene_purity") is not None else "—",
        ])
    lines += [
        _table(
            ["Jaccard", "gold cubierto (cualquier cluster)", "pares 1:1", "recall gold",
             "concordancia clusters", "recall genes", "pureza genes"],
            rows,
        ),
        "",
        "`concordancia clusters` = pares 1:1 / total de clusters. Se evita llamarla "
        "precisión porque los clusters no catalogados no son falsos positivos conocidos.",
        "",
        "### Cobertura del gold a Jaccard ≥ 0,8",
        "",
    ]
    block = (m.get("thresholds") or {}).get("0.8") or {}
    rows = []
    for evidence, payload in (block.get("gold_by_evidence") or {}).items():
        rows.append([f"evidencia: {evidence}", _fmt_int(payload.get("n")),
                     _fmt_int(payload.get("covered_any")), _fmt_pct(payload.get("pct"))])
    for size, payload in (block.get("gold_by_size") or {}).items():
        rows.append([f"tamaño: {size}", _fmt_int(payload.get("n")),
                     _fmt_int(payload.get("covered_any")), _fmt_pct(payload.get("pct"))])
    lines += [
        _table(["grupo", "gold", "cubiertos", "%"], rows),
        "",
        "`cubiertos` = al menos un cluster con Jaccard ≥ 0,8; no es el emparejamiento 1:1.",
        "",
        "### Etiquetas existentes y artefactos",
        "",
        f"- Etiquetas `evidencia` de los clusters: **{_fmt_counts(clusters.get('existing_evidence_labels'))}**.",
        f"- Artefactos: **{_fmt_int((m.get('artifacts') or {}).get('transfer_files'))}** "
        f"`transfer_*.tsv` y **{_fmt_int((m.get('artifacts') or {}).get('model_complex_files'))}** "
        "`model_*.complexes.tsv`; no hay gold para medirlos.",
        "",
        "### Advertencia de circularidad",
        "",
        f"- {(m.get('circularity') or {}).get('reason')}",
        f"- No mide: {(m.get('circularity') or {}).get('what_it_does_not_measure')}",
        "",
        "### No evaluable",
        "",
    ]
    lines += _render_not_evaluable(section)
    lines += ["### Limitaciones", ""]
    lines += _render_limitations(section)
    return lines


def _distance_row(label: str, stats: Mapping[str, Any]) -> list[str]:
    def within(key: str) -> str:
        item = stats.get(key) or {}
        return f"{_fmt_int(item.get('n'))} ({_fmt_pct(item.get('pct'))})"

    return [
        label,
        _fmt_int(stats.get("n")),
        _fmt_num(stats.get("median_bp")),
        within("within_1000bp"),
        within("within_10000bp"),
        within("within_50000bp"),
    ]


def _render_oric(section: Mapping[str, Any]) -> list[str]:
    m = section.get("metrics") or {}
    lines = ["## 4. oriC", ""]
    lines += [
        "**Estado:** parcial. DoriC se separa del consenso; `validated` usa DoriC "
        "como ancla y no se presenta como exactitud independiente.",
        "",
        "### Estados",
        "",
    ]
    rows = [
        [status, _fmt_int(payload.get("n")), _fmt_pct(payload.get("pct"))]
        for status, payload in (m.get("status_counts") or {}).items()
    ]
    lines += [
        _table(["estado", "genomas", "%"], rows),
        "",
        f"- DoriC disponible: **{_fmt_int((m.get('doric_available') or {}).get('n'))}** "
        f"({_fmt_pct((m.get('doric_available') or {}).get('pct'))}).",
        f"- DoriC local: release **{(m.get('doric_dataset') or {}).get('release', '—')}**, "
        f"CSV `{(m.get('doric_dataset') or {}).get('csv_file', '—')}`, "
        f"filas **{_fmt_int((m.get('doric_dataset') or {}).get('csv_rows'))}**.",
        f"- Sin DoriC: **{_fmt_int((m.get('without_doric') or {}).get('n'))}** "
        f"({_fmt_pct((m.get('without_doric') or {}).get('pct'))}); "
        f"estados: {_fmt_counts((m.get('without_doric') or {}).get('status_counts'))}.",
        "",
        "### `validated_by`",
        "",
    ]
    rows = [
        [key, _fmt_int(payload.get("n")), _fmt_pct(payload.get("pct"))]
        for key, payload in (m.get("validated_by_counts") or {}).items()
    ]
    lines += [_table(["validated_by", "genomas", "%"], rows), ""]
    lines += ["### Distancia a DoriC", ""]
    rows = [
        _distance_row("validated (anclado a DoriC)", m.get("validated_anchored") or {}),
        _distance_row("anomalous con DoriC", m.get("anomalous_with_doric") or {}),
        _distance_row("todos con DoriC", m.get("all_with_doric") or {}),
        _distance_row("ORCA-only vs DoriC", m.get("orca_only_vs_doric") or {}),
    ]
    lines += [
        _table(
            ["grupo", "n", "mediana pb", "≤1 kb", "≤10 kb", "≤50 kb"],
            rows,
        ),
        "",
    ]
    validated = m.get("validated_anchored") or {}
    outside = validated.get("outside_50kb") or {}
    lines += [
        f"- `validated` con DoriC: **{_fmt_int(validated.get('n_with_doric'))}**; "
        f"usan DoriC: **{_fmt_int(validated.get('n_using_doric'))}**; "
        f"fuera de 50 kb: **{_fmt_int(outside.get('n'))}**.",
        f"- ORCA-only vs DoriC: `delta_bp` recomputado coincide: "
        f"**{(m.get('orca_only_vs_doric') or {}).get('delta_bp_recomputed_matches')}** "
        f"(mismatches={(m.get('orca_only_vs_doric') or {}).get('delta_bp_mismatches')}).",
        f"- `consensus`: **{_fmt_int((m.get('consensus') or {}).get('n'))}**; "
        f"con DoriC: **{_fmt_int((m.get('consensus') or {}).get('with_doric'))}**; "
        "distancia a verdad externa: **not_evaluable**.",
        f"- `lineage_specific`: **{_fmt_int((m.get('lineage_specific') or {}).get('n'))}**; "
        f"con DoriC: **{_fmt_int((m.get('lineage_specific') or {}).get('with_doric'))}**; "
        "distancia a verdad externa: **not_evaluable**.",
        "",
    ]
    examples = outside.get("examples") or []
    if examples:
        lines += [
            "Ejemplos de `validated` fuera de 50 kb: "
            + ", ".join(f"`{e.get('genome')}` ({_fmt_int(e.get('distance_bp'))} pb)" for e in examples),
            "",
            f"Estos {len(examples)} casos sugieren un desajuste de coordenadas/replicón; "
            "no se corrigen en este benchmark y no deben presentarse como aciertos.",
            "",
        ]
    lines += ["### No evaluable", ""]
    lines += _render_not_evaluable(section)
    lines += ["### Limitaciones", ""]
    lines += _render_limitations(section)
    return lines


def render_markdown(data: Mapping[str, Any]) -> str:
    sections = data.get("sections") or {}
    lines: list[str] = []
    lines += [
        "# Benchmark reproducible por capas",
        "",
        f"**Generado:** {data.get('generated_at', '—')}",
        "",
        f"**Raíz:** `{data.get('project_root', '—')}`",
        "",
        f"**Alcance:** {data.get('scope', '—')}",
        "",
        f"**Comando:** `python tools/benchmark_layers.py --md reports/layer_benchmark.md --json reports/layer_benchmark.json`",
        "",
        "Este informe no cierra brechas ni inventa gold: cada capa declara qué se "
        "pudo medir y qué queda como `not_evaluable`.",
        "",
        "## Resumen",
        "",
    ]
    summary_rows = []
    if "abc" in sections:
        m = sections["abc"].get("metrics") or {}
        labels = m.get("labels") or {}
        summary_rows.append([
            "A/B/C",
            STATUS_LABEL.get(sections["abc"].get("status"), sections["abc"].get("status")),
            f"{_fmt_int(m.get('proteins'))} proteínas; "
            + ", ".join(f"{label.split('_')[0]} {_fmt_pct((labels.get(label) or {}).get('pct'))}" for label in ABC_LABELS),
            "precisión de función/dominios",
        ])
    if "membrane" in sections:
        m = sections["membrane"].get("metrics") or {}
        summary_rows.append([
            "Membrana",
            STATUS_LABEL.get(sections["membrane"].get("status"), sections["membrane"].get("status")),
            f"{_fmt_int(m.get('rows'))} filas; high+medium {_fmt_pct((m.get('high_medium') or {}).get('pct'))}",
            "precisión de localización",
        ])
    if "complexes" in sections:
        m = sections["complexes"].get("metrics") or {}
        block = (m.get("thresholds") or {}).get("0.8") or {}
        one = block.get("one_to_one") or {}
        summary_rows.append([
            "Complejos",
            STATUS_LABEL.get(sections["complexes"].get("status"), sections["complexes"].get("status")),
            f"{_fmt_int(one.get('pairs'))}/{_fmt_int((m.get('gold') or {}).get('n'))} gold 1:1 a J≥0,8; "
            f"{_fmt_int(block.get('clusters_with_best_jaccard_ge', {}).get('n'))}/{_fmt_int((m.get('clusters') or {}).get('n'))} clusters con mejor J≥0,8",
            "precisión predictiva ciega",
        ])
    if "oric" in sections:
        m = sections["oric"].get("metrics") or {}
        sc = m.get("status_counts") or {}
        summary_rows.append([
            "oriC",
            STATUS_LABEL.get(sections["oric"].get("status"), sections["oric"].get("status")),
            f"validated {_fmt_int((sc.get('validated') or {}).get('n'))}, "
            f"consensus {_fmt_int((sc.get('consensus') or {}).get('n'))}, "
            f"anomalous {_fmt_int((sc.get('anomalous') or {}).get('n'))}; "
            f"DoriC {_fmt_int((m.get('doric_available') or {}).get('n'))}",
            "exactitud independiente",
        ])
    if summary_rows:
        lines += [
            _table(["capa", "estado", "métrica principal", "no evaluable"], summary_rows),
            "",
        ]
    if "abc" in sections:
        lines += _render_abc(sections["abc"])
    if "membrane" in sections:
        lines += _render_membrane(sections["membrane"])
    if "complexes" in sections:
        lines += _render_complexes(sections["complexes"])
    if "oric" in sections:
        lines += _render_oric(sections["oric"])
    lines += ["## Fuentes", ""]
    source_rows = []
    for layer, section in sections.items():
        for source in section.get("sources") or []:
            if source:
                source_rows.append([layer, f"`{source}`"])
    lines += [_table(["capa", "archivo"], source_rows), ""]
    lines += [
        "## Limitaciones globales",
        "",
        "- Solo se usaron datos locales; no hay descargas, entrenamiento ni gold externo.",
        "- Las métricas describen estos 400 genomas y este snapshot de datos, no una población.",
        "- Un `not_evaluable` no es un cero: significa que falta la verdad de referencia.",
        "- No se modificó el pipeline, la GUI ni la documentación.",
        "",
    ]
    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark por capas (A/B/C, membrana, complejos, oriC) con datos locales.",
    )
    parser.add_argument(
        "layer",
        nargs="?",
        default="all",
        choices=("all", *ALL_LAYERS),
        help="capa a ejecutar (por defecto: all)",
    )
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    layers = ALL_LAYERS if args.layer == "all" else (args.layer,)
    data = build_data(Path(args.root), layers)
    md = render_markdown(data)

    args.md.parent.mkdir(parents=True, exist_ok=True)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.md.write_text(md, encoding="utf-8")
    args.json.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if not args.quiet:
        print(f"[benchmark] capas={','.join(layers)}")
        for layer in layers:
            section = data["sections"][layer]
            print(f"[benchmark] {layer}: {section.get('status')}")
        print(f"[benchmark] escrito {args.md}")
        print(f"[benchmark] escrito {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
