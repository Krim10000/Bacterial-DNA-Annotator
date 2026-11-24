#!/usr/bin/env python3
"""Holdout explícito de complejos contra catálogos ComplexPortal locales.

El objetivo es evaluar clusters/transferencias **sin** confundir entrenamiento
con generalización. El script:

1. Descubre catálogos locales ``data/complexes/gold_*.tsv`` y
   ``data/complexes/raw/complexportal_*.tsv``.
2. Lee la procedencia del modelo K-12 (``model_k12.genomic.json``) para saber
   qué genoma y qué complejos se usaron en entrenamiento.
3. Aplica dos estrategias de holdout:
   - **por genoma**: solo catálogos de genomas fuera del entrenamiento;
   - **por complejo**: complejos del catálogo crudo que no están en el gold de
     entrenamiento.
4. Si no hay catálogo independiente, devuelve ``not_evaluable`` y no inventa
   precisión. Los cálculos sobre K-12 se etiquetan ``in_sample`` y no se
   presentan como generalización.
5. Mide cobertura de transferencias K-12, que es descriptiva y derivada del
   mismo catálogo; no es precisión.

Sin red. No modifica pipeline, app, GUI ni datos. Solo escribe los reportes que
se le indiquen.

Uso
---
    .venv/bin/python tools/complex_holdout.py
    .venv/bin/python tools/complex_holdout.py --limit 2
    .venv/bin/python tools/complex_holdout.py \
        --out-json reports/complex_holdout.json \
        --out-md reports/complex_holdout.md
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_THRESHOLDS = (0.3, 0.5, 0.7, 0.9)
DEFAULT_MIN_CATALOG_GENES = 2
TAXID_TO_GENOME = {
    "83333": "ecoli_k12",
    "511145": "ecoli_k12",
}
GENOME_ALIASES = {
    "k12": "ecoli_k12",
    "ecoli_k12": "ecoli_k12",
    "83333": "ecoli_k12",
    "511145": "ecoli_k12",
}


# --------------------------------------------------------------------------- #
# Utilidades
# --------------------------------------------------------------------------- #
def _open_tsv(path: Path):
    return path.open("r", encoding="utf-8-sig", errors="replace", newline="")


def _int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def _float(value: object, default: float = 0.0) -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def _round(value: float, ndigits: int = 4) -> float:
    return round(float(value), ndigits)


def _series_stats(values: Sequence[float]) -> dict:
    if not values:
        return {"n": 0, "mean": 0.0, "min": 0.0, "p25": 0.0,
                "median": 0.0, "p75": 0.0, "p90": 0.0, "max": 0.0}
    ordered = sorted(float(v) for v in values)
    n = len(ordered)

    def pick(q: float) -> float:
        pos = (n - 1) * q
        lo = math.floor(pos)
        hi = math.ceil(pos)
        if lo == hi:
            return ordered[lo]
        return ordered[lo] + (ordered[hi] - ordered[lo]) * (pos - lo)

    return {
        "n": n,
        "mean": _round(statistics.fmean(ordered)),
        "min": _round(ordered[0]),
        "p25": _round(pick(0.25)),
        "median": _round(pick(0.5)),
        "p75": _round(pick(0.75)),
        "p90": _round(pick(0.90)),
        "max": _round(ordered[-1]),
    }


def parse_participants(text: object) -> list[tuple[str, int]]:
    """Extrae ``(uniprot, estequiometría)`` de una lista ComplexPortal.

    Acepta ``P76052(1)|P77357(1)`` y proteoformas ``P0AAJ8-PRO_... (0)``.
    La parte de isoforma/proteoforma se descarta para mapear por acceso base.
    """
    out: list[tuple[str, int]] = []
    for chunk in str(text or "").split("|"):
        chunk = chunk.strip()
        if not chunk or chunk == "-":
            continue
        match = re.match(r"^([A-Za-z0-9]+(?:-PRO_[0-9]+)?)\s*\((\d+)\)$", chunk)
        if not match:
            continue
        accession = match.group(1).split("-", 1)[0]
        out.append((accession, int(match.group(2))))
    return out


def jaccard(a: frozenset[str] | set[str], b: frozenset[str] | set[str]) -> float:
    if not a or not b:
        return 0.0
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def _genes_from_csv(text: object) -> frozenset[str]:
    return frozenset(
        item.strip() for item in str(text or "").replace(";", ",").split(",")
        if item.strip()
    )


def _normalize_genome(value: object) -> str:
    key = str(value or "").strip()
    return GENOME_ALIASES.get(key, key)


# --------------------------------------------------------------------------- #
# Estructuras
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CatalogComplex:
    catalog: str
    complex_id: str
    name: str
    genome: str
    genes: frozenset[str]
    declared_subunits: int
    mapped_genes: int
    evidence_code: str
    source: str
    in_gold: bool = False
    in_training_gold: bool = False


@dataclass(frozen=True)
class Prediction:
    source: str
    prediction_id: str
    genes: frozenset[str]
    confidence: str = ""
    evidence: str = ""
    evidence_class: str = ""
    genome: str = ""
    complex_id: str = ""


# --------------------------------------------------------------------------- #
# Parseo de catálogos
# --------------------------------------------------------------------------- #
def load_uniprot_to_locus(path: Path) -> dict[str, frozenset[str]]:
    mapping: dict[str, set[str]] = defaultdict(set)
    if not path.is_file():
        return {}
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            uniprot = (row.get("uniprot") or "").strip()
            locus = (row.get("locus_tag") or row.get("gene_id") or "").strip()
            if uniprot and locus:
                mapping[uniprot].add(locus)
    return {key: frozenset(value) for key, value in mapping.items()}


def genome_from_gold_name(path: Path) -> str:
    stem = path.stem
    suffix = stem[len("gold_"):] if stem.startswith("gold_") else stem
    return _normalize_genome(suffix)


def load_gold_catalog(path: Path) -> list[CatalogComplex]:
    catalog = path.stem
    genome = genome_from_gold_name(path)
    out: list[CatalogComplex] = []
    if not path.is_file():
        return out
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("complex_id") or "").strip()
            if not cid:
                continue
            genes = _genes_from_csv(row.get("locus_tags"))
            declared = _int(row.get("n_genes"), len(genes)) or len(genes)
            out.append(CatalogComplex(
                catalog=catalog,
                complex_id=cid,
                name=(row.get("complex_name") or "").strip(),
                genome=genome,
                genes=genes,
                declared_subunits=declared,
                mapped_genes=len(genes),
                evidence_code=(row.get("evidence_class") or "").strip(),
                source=(row.get("source") or "").strip(),
            ))
    return out


def load_raw_catalog(path: Path,
                     uniprot_to_locus: dict[str, frozenset[str]],
                     ) -> list[CatalogComplex]:
    catalog = path.stem
    out: list[CatalogComplex] = []
    if not path.is_file():
        return out
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("#Complex ac") or "").strip()
            if not cid:
                continue
            taxid = (row.get("Taxonomy identifier") or "").strip()
            genome = TAXID_TO_GENOME.get(taxid, f"taxid_{taxid}" if taxid else "")
            participants_text = (
                row.get("Expanded participant list")
                or row.get("Identifiers (and stoichiometry) of molecules in complex")
                or ""
            )
            participants = parse_participants(participants_text)
            genes: set[str] = set()
            for accession, _ in participants:
                genes.update(uniprot_to_locus.get(accession, ()))
            declared = sum(stoich for _, stoich in participants)
            evidence = (row.get("Evidence Code") or "").split("(", 1)[0].strip()
            out.append(CatalogComplex(
                catalog=catalog,
                complex_id=cid,
                name=(row.get("Recommended name") or "").strip(),
                genome=genome,
                genes=frozenset(genes),
                declared_subunits=declared,
                mapped_genes=len(genes),
                evidence_code=evidence,
                source="ComplexPortal",
            ))
    return out


def discover_catalogs(complexes_dir: Path, raw_dir: Path,
                      uniprot_to_locus: dict[str, frozenset[str]],
                      ) -> tuple[list[CatalogComplex], list[dict]]:
    gold_paths = sorted(Path(p) for p in glob.glob(
        str(complexes_dir / "gold_*.tsv")))
    raw_paths = sorted(Path(p) for p in glob.glob(
        str(raw_dir / "complexportal_*.tsv")))

    gold: list[CatalogComplex] = []
    for path in gold_paths:
        gold.extend(load_gold_catalog(path))
    raw: list[CatalogComplex] = []
    for path in raw_paths:
        raw.extend(load_raw_catalog(path, uniprot_to_locus))

    gold_ids = {c.complex_id for c in gold}
    gold = [replace(c, in_gold=True) for c in gold]
    raw = [replace(c, in_gold=(c.complex_id in gold_ids)) for c in raw]

    catalogs: list[CatalogComplex] = []
    seen: set[tuple[str, str]] = set()
    for catalog_complex in gold + raw:
        key = (catalog_complex.catalog, catalog_complex.complex_id)
        if key in seen:
            continue
        seen.add(key)
        catalogs.append(catalog_complex)

    inventory: list[dict] = []
    for catalog in sorted({c.catalog for c in catalogs}):
        items = [c for c in catalogs if c.catalog == catalog]
        genomes = sorted({c.genome for c in items if c.genome})
        inventory.append({
            "catalogo": catalog,
            "genomas": genomes,
            "complejos": len(items),
            "con_2_o_mas_genes_mapeados": sum(
                1 for c in items if len(c.genes) >= 2),
            "con_menos_de_2_genes_mapeados": sum(
                1 for c in items if len(c.genes) < 2),
            "en_gold": sum(1 for c in items if c.in_gold),
            "fuente": items[0].source if items else "",
            "codigos_evidencia": dict(Counter(c.evidence_code for c in items)),
        })
    return catalogs, inventory


# --------------------------------------------------------------------------- #
# Procedencia del modelo
# --------------------------------------------------------------------------- #
def _gold_ids_from_path(path: Path) -> set[str]:
    ids: set[str] = set()
    if not path.is_file():
        return ids
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("complex_id") or "").strip()
            if cid:
                ids.add(cid)
    return ids


def load_model_provenance(model_json: Path, complexes_dir: Path,
                          all_gold_ids: set[str]) -> dict:
    info = {
        "fuente": "default",
        "trained_genome": "",
        "genomas_entrenamiento": [],
        "gold_path": "",
        "gold_complexes_declarados": None,
        "gold_complexes_encontrados": 0,
        "n_pos": None,
        "n_neg": None,
        "training_gold_ids": sorted(all_gold_ids),
        "advierte_k12": False,
    }
    if not model_json.is_file():
        info["trained_genome"] = "ecoli_k12"
        info["genomas_entrenamiento"] = ["ecoli_k12"]
        info["advierte_k12"] = True
        return info

    data = json.loads(model_json.read_text(encoding="utf-8"))
    trained = _normalize_genome(data.get("trained_genome"))
    gold_path = str(data.get("gold") or "")
    gold_ids: set[str] = set()
    if gold_path:
        path = Path(gold_path)
        if not path.is_absolute():
            path = ROOT / path
        gold_ids = _gold_ids_from_path(path)
    if not gold_ids:
        fallback = complexes_dir / "gold_k12.tsv"
        gold_ids = _gold_ids_from_path(fallback)
    info.update({
        "fuente": "model_json",
        "trained_genome": trained,
        "genomas_entrenamiento": [trained] if trained else [],
        "gold_path": gold_path,
        "gold_complexes_declarados": data.get("gold_complexes"),
        "gold_complexes_encontrados": len(gold_ids),
        "n_pos": data.get("n_pos"),
        "n_neg": data.get("n_neg"),
        "training_gold_ids": sorted(gold_ids or all_gold_ids),
        "advierte_k12": trained == "ecoli_k12",
    })
    return info


# --------------------------------------------------------------------------- #
# Predicciones
# --------------------------------------------------------------------------- #
def load_k12_clusters(path: Path) -> list[Prediction]:
    out: list[Prediction] = []
    if not path.is_file():
        return out
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes = _genes_from_csv(row.get("genes"))
            cid = (row.get("cluster_id") or "").strip()
            if not cid:
                continue
            out.append(Prediction(
                source="clusters_ecoli_k12",
                prediction_id=cid,
                genes=genes,
                confidence=(row.get("confidence") or "").strip().lower(),
                evidence=(row.get("evidence") or "").strip(),
                evidence_class=(row.get("evidencia") or "").strip(),
                genome="ecoli_k12",
            ))
    return out


def load_genome_predictions(complexes_dir: Path, genome: str) -> list[Prediction]:
    candidates = [
        complexes_dir / f"clusters_{genome}.tsv",
        complexes_dir / f"model_{genome}.complexes.tsv",
        complexes_dir / f"{genome}.complexes.tsv",
    ]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        return []
    out: list[Prediction] = []
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if "genes" in (row or {}):
                genes = _genes_from_csv(row.get("genes"))
                pid = (row.get("cluster_id") or row.get("complex_id") or "").strip()
            else:
                raw = (row.get("subunits") or "")
                genes = frozenset(
                    chunk.rsplit(":", 1)[-1].strip()
                    for chunk in raw.replace(",", ";").split(";")
                    if chunk.strip()
                )
                pid = (row.get("complex_id") or "").strip()
            if not pid:
                continue
            out.append(Prediction(
                source=path.stem,
                prediction_id=pid,
                genes=genes,
                confidence=(row.get("confidence") or "").strip().lower(),
                evidence=(row.get("evidence") or "").strip(),
                evidence_class=(row.get("evidencia") or row.get("evidence_class") or "").strip(),
                genome=genome,
            ))
    return out


# --------------------------------------------------------------------------- #
# Emparejamiento y métricas
# --------------------------------------------------------------------------- #
def _eligible_catalog(catalog: Sequence[CatalogComplex],
                      min_genes: int = DEFAULT_MIN_CATALOG_GENES,
                      ) -> list[CatalogComplex]:
    return [c for c in catalog if len(c.genes) >= min_genes]


def best_match(pred_genes: frozenset[str],
               catalog: Sequence[CatalogComplex]) -> tuple[float, str]:
    best_j = 0.0
    best_id = ""
    for item in catalog:
        if len(item.genes) < 2:
            continue
        inter = len(pred_genes & item.genes)
        if not inter:
            continue
        jac = inter / len(pred_genes | item.genes)
        if jac > best_j:
            best_j = jac
            best_id = item.complex_id
    return _round(best_j), best_id


def evaluate_predictions(predictions: Sequence[Prediction],
                         catalog: Sequence[CatalogComplex],
                         thresholds: Sequence[float] = DEFAULT_THRESHOLDS,
                         main_threshold: float = 0.5) -> dict:
    """Empareja cada predicción con el complejo de mayor Jaccard.

    No usa ``curado_id`` ni metadatos de la predicción: solo conjuntos de genes.
    Una predicción es TP si su mejor Jaccard >= umbral; un complejo del catálogo
    se recupera si alguna predicción alcanza el umbral contra él.
    """
    eligible = _eligible_catalog(catalog)
    catalog_best: dict[str, tuple[float, str]] = {
        item.complex_id: (0.0, "") for item in eligible
    }
    pred_best: list[tuple[str, float, str]] = []

    for pred in predictions:
        best_j = 0.0
        best_id = ""
        for item in eligible:
            inter = len(pred.genes & item.genes)
            if not inter:
                continue
            jac = inter / len(pred.genes | item.genes)
            if jac > best_j:
                best_j = jac
                best_id = item.complex_id
            if jac > catalog_best[item.complex_id][0]:
                catalog_best[item.complex_id] = (jac, pred.prediction_id)
        pred_best.append((pred.prediction_id, best_j, best_id))

    n_pred = len(predictions)
    n_cat = len(eligible)
    by_threshold: dict[str, dict] = {}
    for threshold in thresholds:
        tp_pred = sum(1 for _, jac, _ in pred_best if jac >= threshold)
        recalled = sum(
            1 for jac, _ in catalog_best.values() if jac >= threshold)
        precision = tp_pred / n_pred if n_pred else 0.0
        recall = recalled / n_cat if n_cat else 0.0
        f1 = (2 * precision * recall / (precision + recall)
              if precision + recall else 0.0)
        by_threshold[f"{threshold:.1f}"] = {
            "umbral": threshold,
            "predicciones_tp": tp_pred,
            "catalogo_recuperado": recalled,
            "precision": _round(precision),
            "recall": _round(recall),
            "f1": _round(f1),
        }

    principal = by_threshold.get(f"{main_threshold:.1f}", {})
    return {
        "predicciones": n_pred,
        "catalogo_elegible": n_cat,
        "umbral_principal": main_threshold,
        "principal": principal,
        "umbrales": by_threshold,
        "jaccard_prediccion": _series_stats([j for _, j, _ in pred_best]),
        "jaccard_catalogo": _series_stats(
            [j for j, _ in catalog_best.values()]),
        "predicciones_con_solapamiento": sum(
            1 for _, j, _ in pred_best if j > 0),
        "predicciones_jaccard_1": sum(
            1 for _, j, _ in pred_best if j >= 0.9999),
        "catalogo_recuperado_alguna_vez": sum(
            1 for j, _ in catalog_best.values() if j > 0),
    }


def evaluate_by_group(predictions: Sequence[Prediction],
                      catalog: Sequence[CatalogComplex],
                      key: str,
                      main_threshold: float = 0.5) -> dict:
    groups: dict[str, list[Prediction]] = defaultdict(list)
    for pred in predictions:
        value = getattr(pred, key, "") or "sin_dato"
        groups[value].append(pred)
    out: dict[str, dict] = {}
    for value, subset in sorted(groups.items()):
        metrics = evaluate_predictions(
            subset, catalog, thresholds=(main_threshold,),
            main_threshold=main_threshold)
        out[value] = {
            "predicciones": metrics["predicciones"],
            "precision": metrics["principal"].get("precision", 0.0),
            "recall_catalogo_por_grupo": metrics["principal"].get("recall", 0.0),
            "f1": metrics["principal"].get("f1", 0.0),
            "predicciones_tp": metrics["principal"].get("predicciones_tp", 0),
        }
    return out


def holdout_complex_detail(predictions: Sequence[Prediction],
                           catalog: Sequence[CatalogComplex],
                           thresholds: Sequence[float],
                           main_threshold: float) -> dict:
    eligible = _eligible_catalog(catalog)
    detail: list[dict] = []
    for item in eligible:
        best_j = 0.0
        best_pred = ""
        for pred in predictions:
            inter = len(pred.genes & item.genes)
            if not inter:
                continue
            jac = inter / len(pred.genes | item.genes)
            if jac > best_j:
                best_j = jac
                best_pred = pred.prediction_id
        detail.append({
            "complex_id": item.complex_id,
            "complex_name": item.name,
            "genes": len(item.genes),
            "mejor_prediccion": best_pred,
            "jaccard": _round(best_j),
            "recuperado_umbral_principal": best_j >= main_threshold,
        })
    detail.sort(key=lambda item: (-item["jaccard"], item["complex_id"]))
    recalls = {}
    for threshold in thresholds:
        recalled = sum(1 for item in detail if item["jaccard"] >= threshold)
        recalls[f"{threshold:.1f}"] = {
            "umbral": threshold,
            "catalogo_recuperado": recalled,
            "catalogo": len(detail),
            "recall": _round(recalled / len(detail) if detail else 0.0),
        }
    return {
        "catalogo": len(detail),
        "detalle": detail,
        "recall_por_umbral": recalls,
    }


# --------------------------------------------------------------------------- #
# Cobertura de transferencias
# --------------------------------------------------------------------------- #
def load_transfer_coverage(complexes_dir: Path) -> dict:
    rows = 0
    genomes: set[str] = set()
    by_complex: dict[str, set[str]] = defaultdict(set)
    all_genes: set[str] = set()
    for path in sorted(Path(p) for p in glob.glob(
            str(complexes_dir / "transfer_*.tsv"))):
        genome = path.stem[len("transfer_"):]
        genomes.add(genome)
        with _open_tsv(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rows += 1
                cid = (row.get("complex_id") or "").strip()
                if not cid:
                    continue
                genes = _genes_from_csv(row.get("genes_k12"))
                by_complex[cid].update(genes)
                all_genes.update(genes)
    return {
        "filas": rows,
        "genomas": len(genomes),
        "complejos_unicos": len(by_complex),
        "genes_k12_unicos": len(all_genes),
        "por_complejo": {cid: frozenset(genes)
                         for cid, genes in by_complex.items()},
    }


def transfer_catalog_coverage(transfer: dict,
                              gold_catalog: Sequence[CatalogComplex],
                              raw_catalog: Sequence[CatalogComplex]) -> dict:
    gold_ids = {c.complex_id for c in gold_catalog}
    raw_ids = {c.complex_id for c in raw_catalog}
    transfer_ids = set(transfer["por_complejo"])
    gold_genes = set().union(*(c.genes for c in gold_catalog)) if gold_catalog else set()
    raw_genes = set().union(*(c.genes for c in raw_catalog)) if raw_catalog else set()
    transfer_genes = set().union(*transfer["por_complejo"].values()) \
        if transfer["por_complejo"] else set()

    details: list[dict] = []
    for item in gold_catalog:
        if len(item.genes) < 2:
            continue
        observed = transfer["por_complejo"].get(item.complex_id, frozenset())
        frac = (len(observed & item.genes) / len(item.genes)
                if item.genes else 0.0)
        details.append({
            "complex_id": item.complex_id,
            "genes": len(item.genes),
            "observados": len(observed & item.genes),
            "fraccion": _round(frac),
        })
    fractions = [item["fraccion"] for item in details]
    return {
        "estado": "cobertura_derivada_no_holdout",
        "filas_transferencia": transfer["filas"],
        "genomas_con_transferencia": transfer["genomas"],
        "complejos_transferidos": transfer["complejos_unicos"],
        "complejos_transferidos_en_gold": len(transfer_ids & gold_ids),
        "complejos_transferidos_en_raw": len(transfer_ids & raw_ids),
        "complejos_transferidos_fuera_de_catalogo": len(
            transfer_ids - raw_ids - gold_ids),
        "cobertura_complejos_gold": _round(
            len(transfer_ids & gold_ids) / len(gold_ids) if gold_ids else 0.0),
        "cobertura_complejos_raw": _round(
            len(transfer_ids & raw_ids) / len(raw_ids) if raw_ids else 0.0),
        "genes_k12_en_transferencias": len(transfer_genes),
        "cobertura_genes_gold": _round(
            len(transfer_genes & gold_genes) / len(gold_genes)
            if gold_genes else 0.0),
        "cobertura_genes_raw": _round(
            len(transfer_genes & raw_genes) / len(raw_genes)
            if raw_genes else 0.0),
        "complejos_gold_completos": sum(1 for f in fractions if f >= 0.9999),
        "complejos_gold_parciales": sum(
            1 for f in fractions if 0 < f < 0.9999),
        "complejos_gold_sin_cobertura": sum(1 for f in fractions if f == 0),
        "complejos_gold_sin_cobertura_ids": sorted(
            item["complex_id"] for item in details if item["fraccion"] == 0),
        "fraccion_media_por_complejo": _round(
            statistics.fmean(fractions) if fractions else 0.0),
        "fraccion_mediana_por_complejo": _round(
            statistics.median(fractions) if fractions else 0.0),
        "detalle_complejos": details,
        "advertencia": (
            "Las transferencias se derivan del mismo catálogo K-12; esto mide "
            "cobertura, no precisión ni generalización."
        ),
    }


# --------------------------------------------------------------------------- #
# Holdouts
# --------------------------------------------------------------------------- #
def evaluate_genome_holdout(catalogs: Sequence[CatalogComplex],
                            training_genomes: Sequence[str],
                            complexes_dir: Path,
                            thresholds: Sequence[float],
                            main_threshold: float) -> dict:
    training = {_normalize_genome(g) for g in training_genomes if g}
    independent = [
        c for c in catalogs
        if c.genome and _normalize_genome(c.genome) not in training
    ]
    genomes = sorted({c.genome for c in independent})
    if not genomes:
        return {
            "estado": "not_evaluable",
            "motivo": (
                "No hay catálogo ComplexPortal para un genoma fuera del "
                "entrenamiento; el único catálogo local es K-12."
            ),
            "genomas_entrenamiento": sorted(training),
            "catalogos_independientes": [],
            "evaluaciones": [],
        }

    evaluations: list[dict] = []
    for genome in genomes:
        catalog = [c for c in independent if c.genome == genome]
        if not _eligible_catalog(catalog):
            evaluations.append({
                "genoma": genome,
                "estado": "not_evaluable",
                "motivo": "catálogo independiente sin complejos con >=2 genes",
                "catalogo": 0,
            })
            continue
        predictions = load_genome_predictions(complexes_dir, genome)
        if not predictions:
            evaluations.append({
                "genoma": genome,
                "estado": "not_evaluable",
                "motivo": "catálogo independiente sin predicciones mapeadas",
                "catalogo": len(catalog),
            })
            continue
        catalog_genes = set().union(*(c.genes for c in catalog))
        prediction_genes = set().union(*(p.genes for p in predictions))
        if catalog_genes and prediction_genes and not (
                catalog_genes & prediction_genes):
            evaluations.append({
                "genoma": genome,
                "estado": "not_evaluable",
                "motivo": (
                    "El catálogo y las predicciones no comparten ningún gen; "
                    "probablemente usan namespaces distintos."
                ),
                "catalogo": len(_eligible_catalog(catalog)),
                "predicciones": len(predictions),
            })
            continue
        metrics = evaluate_predictions(
            predictions, catalog, thresholds, main_threshold)
        evaluations.append({
            "genoma": genome,
            "estado": "evaluable",
            "catalogo": len(_eligible_catalog(catalog)),
            "predicciones": len(predictions),
            "metricas": metrics,
        })

    evaluable = any(item["estado"] == "evaluable" for item in evaluations)
    return {
        "estado": "evaluable" if evaluable else "not_evaluable",
        "motivo": "" if evaluable else (
            "Hay catálogo independiente pero no predicciones mapeadas."),
        "genomas_entrenamiento": sorted(training),
        "catalogos_independientes": genomes,
        "evaluaciones": evaluations,
    }


def evaluate_complex_holdout(catalogs: Sequence[CatalogComplex],
                             training_gold_ids: set[str],
                             predictions: Sequence[Prediction],
                             thresholds: Sequence[float],
                             main_threshold: float) -> dict:
    raw = [c for c in catalogs if c.source == "ComplexPortal"
           and c.catalog.startswith("complexportal_")]
    holdout_all = [
        c for c in raw
        if c.complex_id not in training_gold_ids
    ]
    eligible = _eligible_catalog(holdout_all)
    excluded = [c for c in holdout_all if len(c.genes) < 2]
    detail = holdout_complex_detail(
        predictions, eligible, thresholds, main_threshold)
    return {
        "estado": "exploratorio_no_independiente",
        "definicion": (
            "Complejos presentes en ComplexPortal crudo pero ausentes del gold "
            "usado para entrenar. El genoma K-12 sigue siendo el genoma de "
            "entrenamiento; no es un holdout independiente."
        ),
        "catalogo_total": len(raw),
        "catalogo_holdout": len(holdout_all),
        "catalogo_holdout_evaluable": len(eligible),
        "excluidos_menos_de_2_genes": len(excluded),
        "ids_excluidos": sorted(c.complex_id for c in excluded),
        "recall": detail["recall_por_umbral"],
        "detalle": detail["detalle"],
        "precision": {
            "estado": "not_evaluable",
            "motivo": (
                "La precisión requeriría puntuar contra todo el catálogo, que "
                "incluye complejos de entrenamiento; se reporta solo recall "
                "del holdout."
            ),
        },
        "advertencia": (
            "La selección no es aleatoria: muchos complejos raw-only son "
            "homómeros/unidades sin mapeo. No usar estos recall como "
            "generalización."
        ),
    }


# --------------------------------------------------------------------------- #
# Análisis
# --------------------------------------------------------------------------- #
def analyze(complexes_dir: Path, raw_dir: Path, k12_genes_path: Path,
            model_json: Path, clusters_path: Path,
            out_json: Path, out_md: Path,
            thresholds: Sequence[float] = DEFAULT_THRESHOLDS,
            main_threshold: float = 0.5) -> dict:
    complexes_dir = Path(complexes_dir)
    raw_dir = Path(raw_dir)
    k12_genes_path = Path(k12_genes_path)
    model_json = Path(model_json)
    clusters_path = Path(clusters_path)
    thresholds = tuple(sorted({float(t) for t in thresholds}))
    main_threshold = float(main_threshold)

    uniprot_map = load_uniprot_to_locus(k12_genes_path)
    catalogs, inventory = discover_catalogs(complexes_dir, raw_dir, uniprot_map)
    all_gold_ids = {c.complex_id for c in catalogs
                    if c.catalog.startswith("gold_")}
    provenance = load_model_provenance(
        model_json, complexes_dir, all_gold_ids)

    training_genomes = provenance["genomas_entrenamiento"]
    training_gold_ids = set(provenance["training_gold_ids"])
    catalogs = [
        replace(c, in_training_gold=(c.complex_id in training_gold_ids))
        for c in catalogs
    ]

    k12_predictions = load_k12_clusters(clusters_path)
    gold_catalog = [
        c for c in catalogs if c.catalog.startswith("gold_")
        and c.complex_id in training_gold_ids
    ]
    raw_catalog = [
        c for c in catalogs if c.catalog.startswith("complexportal_")
    ]

    genome_holdout = evaluate_genome_holdout(
        catalogs, training_genomes, complexes_dir, thresholds, main_threshold)
    complex_holdout = evaluate_complex_holdout(
        catalogs, training_gold_ids, k12_predictions,
        thresholds, main_threshold)

    in_sample: dict
    if k12_predictions and gold_catalog:
        metrics = evaluate_predictions(
            k12_predictions, gold_catalog, thresholds, main_threshold)
        in_sample = {
            "estado": "in_sample_k12_no_holdout",
            "predicciones": len(k12_predictions),
            "catalogo": len(_eligible_catalog(gold_catalog)),
            "metricas": metrics,
            "por_confianza": evaluate_by_group(
                k12_predictions, gold_catalog, "confidence", main_threshold),
            "por_evidencia": evaluate_by_group(
                k12_predictions, gold_catalog, "evidence_class", main_threshold),
            "advertencia": (
                "K-12 es el genoma de entrenamiento y el gold de 273 complejos "
                "se usó para entrenar. Estas cifras son in-sample; no son "
                "precisión de generalización."
            ),
        }
    else:
        in_sample = {
            "estado": "not_evaluable",
            "motivo": "faltan clusters K-12 o catálogo gold de entrenamiento",
        }

    transfer = load_transfer_coverage(complexes_dir)
    transfer_block = transfer_catalog_coverage(
        transfer, gold_catalog, raw_catalog)

    independent_genomes = genome_holdout["catalogos_independientes"]
    holdout_state = "evaluable" if genome_holdout["estado"] == "evaluable" \
        else "not_evaluable"
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    report = {
        "meta": {
            "herramienta": "tools/complex_holdout.py",
            "generado_utc": generated,
            "sin_red": True,
            "no_modifica_datos": True,
            "advertencia": (
                "Holdout de complejos: si no hay catálogo independiente, el "
                "resultado es not_evaluable. No se inventa precisión."
            ),
        },
        "inputs": {
            "complexes_dir": str(complexes_dir),
            "raw_dir": str(raw_dir),
            "k12_genes": str(k12_genes_path),
            "model_json": str(model_json),
            "clusters_k12": str(clusters_path),
            "genes_k12_mapeados": len(uniprot_map),
        },
        "procedencia_modelo": provenance,
        "catalogos": inventory,
        "holdout": {
            "estado": holdout_state,
            "evaluable_independiente": bool(independent_genomes),
            "motivo": genome_holdout["motivo"],
            "genomas_entrenamiento": training_genomes,
            "catalogos_independientes": independent_genomes,
            "estrategias": {
                "genoma": genome_holdout,
                "complejo": complex_holdout,
            },
        },
        "evaluaciones": {
            "in_sample_k12": in_sample,
            "holdout_complejo_k12": complex_holdout,
            "transferencias_k12_cobertura": transfer_block,
        },
        "limitaciones": _limitations(
            catalogs=len(catalogs),
            gold=len(gold_catalog),
            raw=len(raw_catalog),
            in_sample=in_sample,
            complex_holdout=complex_holdout,
            transfer_block=transfer_block,
            provenance=provenance,
        ),
    }
    return report


def _limitations(catalogs: int, gold: int, raw: int, in_sample: dict,
                 complex_holdout: dict, transfer_block: dict,
                 provenance: dict) -> list[str]:
    return [
        (
            f"Hay {gold} complejos gold + {raw} crudos locales "
            f"({catalogs} filas de catálogo) y todos pertenecen a K-12 "
            "(ComplexPortal taxón 83333). No hay catálogo curado de otra "
            "especie; por eso el holdout independiente es `not_evaluable`."
        ),
        (
            f"El modelo declara `trained_genome={provenance.get('trained_genome')}` "
            f"y {provenance.get('gold_complexes_encontrados', 0)} complejos gold "
            "usados en entrenamiento. Cualquier evaluación sobre K-12 es "
            "in-sample y no mide generalización."
        ),
        (
            f"El holdout por complejo usa {complex_holdout['catalogo_holdout_evaluable']} "
            f"complejos raw-only con >=2 genes mapeados; excluye "
            f"{complex_holdout['excluidos_menos_de_2_genes']} por falta de "
            "genes/estequiometría. No es una muestra aleatoria: los homómeros y "
            "complejos sin mapeo quedan fuera, y el genoma sigue siendo K-12."
        ),
        (
            "El emparejamiento es por Jaccard de conjuntos de genes. No usa "
            "`curado_id`, no evalúa estequiometría, isoformas ni interfaces, y "
            "puede fallar con complejos que comparten subunidades."
        ),
        (
            f"Las transferencias cubren {transfer_block['complejos_transferidos']} "
            "complejos y se derivan del mismo catálogo K-12; su cobertura no es "
            "precisión. RBH no trae cobertura de alineamiento ni demuestra "
            "interacción física."
        ),
        (
            "Los clusters `model_<genoma>.complexes.tsv` de los 400 genomas no "
            "se evalúan aquí porque no hay catálogo ComplexPortal de sus "
            "especies; solo K-12 tiene gold local."
        ),
        (
            "No hay validación experimental en otras especies, ni gold "
            "independiente para los 400 genomas. El informe no cierra brechas "
            "ni cambia el pipeline, la app/GUI o los datos."
        ),
    ]


# --------------------------------------------------------------------------- #
# Markdown
# --------------------------------------------------------------------------- #
def _es_int(value: object) -> str:
    try:
        return f"{int(value):,}".replace(",", ".")
    except (TypeError, ValueError):
        return str(value)


def _es_dec(value: object, ndigits: int = 2) -> str:
    try:
        return f"{float(value):.{ndigits}f}".replace(".", ",")
    except (TypeError, ValueError):
        return str(value)


def _pct(value: float, ndigits: int = 1) -> str:
    return f"{100.0 * float(value):.{ndigits}f}".replace(".", ",") + " %"


def _md_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    lines = ["| " + " | ".join(str(h) for h in headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)


def _short(text: object, limit: int = 55) -> str:
    value = str(text or "").strip().replace("|", "/")
    return value if len(value) <= limit else value[: limit - 1] + "…"


def build_markdown(report: dict) -> str:
    provenance = report["procedencia_modelo"]
    holdout = report["holdout"]
    in_sample = report["evaluaciones"]["in_sample_k12"]
    complex_holdout = report["evaluaciones"]["holdout_complejo_k12"]
    transfer = report["evaluaciones"]["transferencias_k12_cobertura"]
    metrics = in_sample.get("metricas", {})
    principal = metrics.get("principal", {})

    lines: list[str] = []
    lines.append("# Holdout de complejos contra catálogos ComplexPortal locales")
    lines.append("")
    lines.append(
        "> **Veredicto: `not_evaluable` para un holdout independiente.** "
        "Los únicos catálogos locales son de *E. coli* K-12 y el modelo se "
        "entrenó en K-12 con 273 complejos gold. No hay catálogo de otra "
        "especie ni validación experimental. Lo que sigue separa lo "
        "**in-sample** de lo **exploratorio**; nada de esto es precisión de "
        "generalización.")
    lines.append("")
    lines.append(f"Generado: `{report['meta']['generado_utc']}` · "
                 f"herramienta: `{report['meta']['herramienta']}`")
    lines.append("")

    lines.append("## 1. Veredicto del holdout")
    lines.append("")
    lines.append(_md_table(
        ["Estrategia", "Estado", "Detalle"],
        [
            ["Por genoma", holdout["estado"],
             holdout["motivo"] or "—"],
            ["Por complejo", complex_holdout["estado"],
             "Complejos raw-only ausentes del gold de entrenamiento"],
        ]))
    lines.append("")
    lines.append(
        f"- Genomas de entrenamiento declarados: "
        f"`{', '.join(holdout['genomas_entrenamiento']) or 'desconocido'}`.")
    lines.append(
        f"- Catálogos independientes encontrados: "
        f"{_es_int(len(holdout['catalogos_independientes']))}.")
    lines.append(
        f"- El holdout por genoma es **{holdout['estado']}** porque no existe "
        "catálogo de un genoma no entrenado.")
    lines.append(
        "- El holdout por complejo es **exploratorio**: los complejos "
        "raw-only no estaban en el gold, pero el genoma K-12 sí es el genoma "
        "de entrenamiento.")
    lines.append(
        "- **No se puede concluir** precisión/recall de generalización. "
        "**Sí se puede** reportar cobertura de catálogo, recall exploratorio "
        "de los complejos raw-only y contexto in-sample.")
    lines.append("")

    lines.append("## 2. Procedencia del modelo y catálogos")
    lines.append("")
    lines.append(_md_table(
        ["Dato", "Valor"],
        [
            ["Genoma de entrenamiento", provenance.get("trained_genome") or "—"],
            ["Gold de entrenamiento (declarado)", _es_int(provenance.get("gold_complexes_declarados") or 0)],
            ["Gold de entrenamiento (encontrado)", _es_int(provenance.get("gold_complexes_encontrados") or 0)],
            ["Pares positivos (`n_pos`)", _es_int(provenance.get("n_pos") or 0)],
            ["Pares negativos (`n_neg`)", _es_int(provenance.get("n_neg") or 0)],
            ["Fuente de procedencia", provenance.get("fuente") or "—"],
        ]))
    lines.append("")
    lines.append(_md_table(
        ["Catálogo", "Genoma", "Complejos", "≥2 genes", "<2 genes", "En gold"],
        [[item["catalogo"], ", ".join(item["genomas"]) or "—",
          _es_int(item["complejos"]),
          _es_int(item["con_2_o_mas_genes_mapeados"]),
          _es_int(item["con_menos_de_2_genes_mapeados"]),
          _es_int(item["en_gold"])]
         for item in report["catalogos"]]))
    lines.append("")
    lines.append(
        f"El catálogo crudo ComplexPortal contiene "
        f"{_es_int(sum(item['complejos'] for item in report['catalogos'] if item['catalogo'].startswith('complexportal_')))} "
        f"complejos; el gold de entrenamiento tiene "
        f"{_es_int(provenance.get('gold_complexes_encontrados') or 0)}. Los "
        f"{_es_int(complex_holdout['catalogo_holdout'])} ausentes del gold son "
        f"la base del holdout por complejo, pero solo "
        f"{_es_int(complex_holdout['catalogo_holdout_evaluable'])} tienen ≥2 "
        "genes mapeados.")
    lines.append("")

    lines.append("## 3. Holdout por genoma")
    lines.append("")
    if holdout["catalogos_independientes"]:
        for evaluation in holdout["estrategias"]["genoma"]["evaluaciones"]:
            lines.append(f"### {evaluation['genoma']}")
            if evaluation.get("estado") != "evaluable":
                lines.append(f"No evaluable: {evaluation.get('motivo', '—')}")
                continue
            m = evaluation["metricas"]
            p = m["principal"]
            lines.append(_md_table(
                ["Métrica", "Valor"],
                [["Predicciones", _es_int(m["predicciones"])],
                 ["Catálogo elegible", _es_int(m["catalogo_elegible"])],
                 [f"Precisión @ {m['umbral_principal']}", _es_dec(p.get("precision", 0), 3)],
                 [f"Recall @ {m['umbral_principal']}", _es_dec(p.get("recall", 0), 3)],
                 ["F1", _es_dec(p.get("f1", 0), 3)]]))
            lines.append("")
    else:
        lines.append(
            "**No evaluable.** No hay `gold_<genoma>.tsv` ni "
            "`raw/complexportal_<taxid>.tsv` para un genoma fuera de "
            "`ecoli_k12`. No se puede estimar precisión/recall de "
            "generalización con los datos locales.")
    lines.append("")

    lines.append("## 4. Holdout por complejo (exploratorio)")
    lines.append("")
    lines.append(
        f"Se usan {_es_int(complex_holdout['catalogo_holdout_evaluable'])} "
        f"complejos de ComplexPortal crudo que no están en el gold de "
        f"entrenamiento y tienen ≥2 genes mapeados. Se excluyen "
        f"{_es_int(complex_holdout['excluidos_menos_de_2_genes'])} por tener "
        "<2 genes mapeados (homómeros o subunidades no mapeadas). El genoma "
        "sigue siendo K-12; esto **no** "
        "es un holdout independiente.")
    lines.append("")
    recall_rows = []
    for key, value in complex_holdout["recall"].items():
        recall_rows.append([value["umbral"], _es_int(value["catalogo_recuperado"]),
                            _es_int(value["catalogo"]), _es_dec(value["recall"], 3)])
    lines.append(_md_table(
        ["Umbral Jaccard", "Recuperados", "Catálogo", "Recall"], recall_rows))
    lines.append("")
    detail_rows = []
    for item in complex_holdout["detalle"][:10]:
        detail_rows.append([
            item["complex_id"], _es_int(item["genes"]),
            item["mejor_prediccion"] or "—",
            _es_dec(item["jaccard"], 3),
            "sí" if item["recuperado_umbral_principal"] else "no",
        ])
    lines.append(_md_table(
        ["Complejo", "Genes", "Mejor predicción", "Jaccard", "Recuperado @ 0,5"],
        detail_rows))
    lines.append("")
    lines.append(
        f"**Precisión: `not_evaluable`.** {complex_holdout['precision']['motivo']}")
    lines.append("")
    lines.append(f"_{complex_holdout['advertencia']}_")
    lines.append("")

    lines.append("## 5. Contexto in-sample K-12 (no es holdout)")
    lines.append("")
    if in_sample.get("estado") == "in_sample_k12_no_holdout":
        lines.append(
            f"Predicciones `clusters_ecoli_k12_evidencia.tsv` contra "
            f"{_es_int(in_sample['catalogo'])} complejos gold. El modelo vio "
            "estos complejos al entrenar. **No usar como precisión de "
            "generalización.**")
        lines.append("")
        lines.append("### 5.1 Barrido de umbral")
        lines.append("")
        rows = []
        for _, value in metrics["umbrales"].items():
            rows.append([value["umbral"], _es_int(value["predicciones_tp"]),
                         _es_int(value["catalogo_recuperado"]),
                         _es_dec(value["precision"], 3),
                         _es_dec(value["recall"], 3),
                         _es_dec(value["f1"], 3)])
        lines.append(_md_table(
            ["Umbral", "Pred. TP", "Catálogo recuperado", "Precisión",
             "Recall", "F1"], rows))
        lines.append("")
        lines.append(
            f"En el umbral principal {metrics['umbral_principal']}: "
            f"precisión **{_es_dec(principal.get('precision', 0), 3)}**, "
            f"recall **{_es_dec(principal.get('recall', 0), 3)}**, "
            f"F1 **{_es_dec(principal.get('f1', 0), 3)}**.")
        lines.append("")
        lines.append("### 5.2 Por confianza")
        lines.append("")
        lines.append(_md_table(
            ["Confianza", "Predicciones", "Precisión", "Recall del catálogo",
             "F1"],
            [[conf, _es_int(data["predicciones"]),
              _es_dec(data["precision"], 3),
              _es_dec(data["recall_catalogo_por_grupo"], 3),
              _es_dec(data["f1"], 3)]
             for conf, data in in_sample["por_confianza"].items()]))
        lines.append("")
        lines.append("### 5.3 Por evidencia")
        lines.append("")
        lines.append(_md_table(
            ["Evidencia (clase)", "Predicciones", "Precisión",
             "Recall del catálogo", "F1"],
            [[ev, _es_int(data["predicciones"]),
              _es_dec(data["precision"], 3),
              _es_dec(data["recall_catalogo_por_grupo"], 3),
              _es_dec(data["f1"], 3)]
             for ev, data in in_sample["por_evidencia"].items()]))
        lines.append("")
        lines.append(
            "`nuevo` significa no documentado en el gold; su precisión 0 es "
            "por construcción y no un error del modelo.")
        lines.append("")
        lines.append(f"_{in_sample['advertencia']}_")
    else:
        lines.append(f"**No evaluable:** {in_sample.get('motivo', '—')}")
    lines.append("")

    lines.append("## 6. Cobertura de transferencias K-12 (no es precisión)")
    lines.append("")
    lines.append(_md_table(
        ["Métrica", "Valor"],
        [
            ["Filas de transferencia", _es_int(transfer["filas_transferencia"])],
            ["Genomas con transferencia", _es_int(transfer["genomas_con_transferencia"])],
            ["Complejos transferidos", _es_int(transfer["complejos_transferidos"])],
            ["Transferidos en gold", _es_int(transfer["complejos_transferidos_en_gold"])],
            ["Cobertura de complejos gold", _pct(transfer["cobertura_complejos_gold"])],
            ["Genes K-12 en transferencias", _es_int(transfer["genes_k12_en_transferencias"])],
            ["Cobertura de genes gold", _pct(transfer["cobertura_genes_gold"])],
            ["Complejos gold completos", _es_int(transfer["complejos_gold_completos"])],
            ["Complejos gold parciales", _es_int(transfer["complejos_gold_parciales"])],
            ["Complejos gold sin cobertura", _es_int(transfer["complejos_gold_sin_cobertura"])],
            ["Fracción media por complejo", _es_dec(transfer["fraccion_media_por_complejo"], 3)],
        ]))
    lines.append("")
    lines.append(f"_{transfer['advertencia']}_")
    lines.append("")

    lines.append("## 7. Limitaciones")
    lines.append("")
    for limitation in report["limitaciones"]:
        lines.append(f"- {limitation}")
    lines.append("")

    lines.append("## 8. Reproducir")
    lines.append("")
    lines.append("```bash")
    lines.append(".venv/bin/python -m pytest tests/test_complex_holdout.py -q")
    lines.append(".venv/bin/python tools/complex_holdout.py \\")
    lines.append("    --out-json reports/complex_holdout.json \\")
    lines.append("    --out-md reports/complex_holdout.md")
    lines.append("```")
    lines.append("")
    lines.append("El script no usa red y solo escribe los reportes indicados.")
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Holdout de complejos contra catálogos ComplexPortal "
                    "locales; no inventa precisión si no hay gold independiente.")
    parser.add_argument("--complexes-dir", type=Path,
                        default=ROOT / "data" / "complexes")
    parser.add_argument("--raw-dir", type=Path,
                        default=ROOT / "data" / "complexes" / "raw")
    parser.add_argument("--k12-genes", type=Path,
                        default=ROOT / "data" / "complexes" / "k12_genes.tsv")
    parser.add_argument("--model-json", type=Path,
                        default=ROOT / "data" / "complexes" / "model_k12.genomic.json")
    parser.add_argument("--clusters", type=Path,
                        default=ROOT / "data" / "complexes" /
                        "clusters_ecoli_k12_evidencia.tsv")
    parser.add_argument("--out-json", type=Path,
                        default=ROOT / "reports" / "complex_holdout.json")
    parser.add_argument("--out-md", type=Path,
                        default=ROOT / "reports" / "complex_holdout.md")
    parser.add_argument("--thresholds", type=str, default="0.3,0.5,0.7,0.9")
    parser.add_argument("--main-threshold", type=float, default=0.5)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    thresholds = tuple(
        float(item) for item in args.thresholds.split(",") if item.strip())
    report = analyze(
        complexes_dir=args.complexes_dir,
        raw_dir=args.raw_dir,
        k12_genes_path=args.k12_genes,
        model_json=args.model_json,
        clusters_path=args.clusters,
        out_json=args.out_json,
        out_md=args.out_md,
        thresholds=thresholds,
        main_threshold=args.main_threshold,
    )
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")
    args.out_md.write_text(build_markdown(report), encoding="utf-8")

    holdout = report["holdout"]
    print(
        f"holdout={holdout['estado']} "
        f"catalogos_independientes={len(holdout['catalogos_independientes'])} "
        f"holdout_complejo={report['evaluaciones']['holdout_complejo_k12']['catalogo_holdout_evaluable']} "
        f"in_sample_precision={report['evaluaciones']['in_sample_k12'].get('metricas', {}).get('principal', {}).get('precision', 'NA')}"
    )
    print(f"JSON: {args.out_json}")
    print(f"MD:   {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
