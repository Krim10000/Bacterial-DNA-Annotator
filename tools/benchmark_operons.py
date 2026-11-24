#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Benchmark local de operones predichos por reglas contra un gold tipo RegulonDB.

Objetivo
--------
Medir la estructura de operones de ``<genoma>.operons_auto.tsv`` (o de un GFF
de operones) contra un gold local. **No hay red, no se descarga nada y no se
inventa gold**: si no se encuentra un gold aplicable, el resultado es
``not_evaluable`` con la razón y las rutas buscadas.

Dos vistas, porque un operón es a la vez un conjunto de genes y un intervalo
-----------------------------------------------------------------------------
- **Nivel gen**: pares de genes declarados co-transcritos (métrica estándar,
  robusta a diferencias de borde) y *exact-set* (el conjunto de genes coincide
  exactamente). Requiere IDs de gen compartidos entre predicción y gold.
- **Nivel intervalo**: solapamiento recíproco mínimo configurable entre el
  intervalo predicho y el de gold, con emparejamiento 1:1; y *exact-interval*
  (mismas coordenadas). No requiere IDs de gen, así que sigue siendo evaluable
  cuando ``operons_auto`` usa ``cds_N`` y el gold usa b-numbers.

Cada vista reporta precisión, recall, F1 y exact-match, global y separada por
hebra ``+`` / ``-``.

Entradas
--------
- ``--pred``: TSV de operones (``operons_auto``, reglas, etc.) o GFF/GFF3 con
  features de operón. Si se omite, busca rutas locales conocidas.
- ``--gold``: TSV/GFF de gold. Si se omite, busca gold específico del genoma;
  solo usa RegulonDB K-12 cuando ``--genome`` es K-12.
- ``--genome``: nombre del genoma para autodetectar ``operons_auto``, GFF y
  gold. Si se omite, se infiere del nombre de ``--pred``/``--gold``.

Formatos TSV aceptados
----------------------
Se detectan alias de columnas: ``genes``/``gene_ids``/``gene_list``/``members``,
``strand``/``sign``, ``start``/``begin``, ``end``/``stop``,
``contig``/``seqid`` y una columna de ID. Los genes se separan por coma,
punto y coma, barra o espacios. Si existe ``mapping_ok``, se excluyen las filas
distintas de ``yes`` (convención de ``data/regulondb/gt/gt_operons.tsv``).

Límites declarados
------------------
El benchmark mide consistencia contra un gold curado, no generalización: para
K-12 las reglas pueden haberse calibrado con el mismo gold. No demuestra
cotranscripción, no evalúa promotores/TSS y no corrige diferencias de
nomenclatura de genes. Un ``not_evaluable`` no es un cero de rendimiento.

Uso
---
    python tools/benchmark_operons.py --genome ecoli_k12
    python tools/benchmark_operons.py --gold data/regulondb/gt/gt_operons.tsv \
        --pred data/operon_finder/ecoli_k12_operons_rules.tsv --genome ecoli_k12 \
        --json reports/operon_benchmark.json --md reports/operon_benchmark.md
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = ROOT / "reports" / "operon_benchmark.json"
DEFAULT_MD = ROOT / "reports" / "operon_benchmark.md"
SCHEMA = "operon_benchmark/v1"

GENE_KEYS = ("genes", "gene_ids", "gene_list", "members", "genelist",
             "gene_names")
STRAND_KEYS = ("strand", "sign", "dir", "orientation")
START_KEYS = ("start", "begin", "left", "first_gene_pos")
END_KEYS = ("end", "stop", "right", "last_gene_pos")
CONTIG_KEYS = ("contig", "seqid", "chromosome", "chr", "replicon")
ID_KEYS = ("operon_id", "operon", "id", "name", "label")
MAP_KEYS = ("mapping_ok", "map_ok")
GFF_OPERON_TYPES = {
    "operon", "operon_region", "transcription_unit", "transcriptional_unit",
    "tu", "regulon", "transcription_unit_region",
}
GFF_GENE_KEYS = ("genes", "gene_ids", "gene_list", "members", "genelist",
                 "gene")
GFF_GROUP_KEYS = ("operon_id", "operon", "tu_id", "transcription_unit",
                  "transcription_unit_id")
K12_NAMES = {"ecoli_k12", "ecoli-k12", "e_coli_k12", "escherichia_coli_k12",
             "k12", "k-12"}


# ---------------------------------------------------------------------------
# Estructuras y utilidades
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Operon:
    operon_id: str
    contig: str
    strand: str
    start: int | None
    end: int | None
    genes: frozenset[str]

    @property
    def has_coords(self) -> bool:
        return self.start is not None and self.end is not None


@dataclass
class ParsedTable:
    path: Path
    fmt: str
    operons: list[Operon]
    all_genes: set[str]
    warnings: list[str]
    n_rows: int
    n_skipped_mapping: int = 0
    n_skipped_no_genes: int = 0
    n_without_coords: int = 0
    n_unknown_strand: int = 0


def _detect(fieldnames: Sequence[str], keys: Sequence[str]) -> str | None:
    lower = {str(f).lower().strip(): f for f in fieldnames}
    for key in keys:
        if key in lower:
            return lower[key]
    for key in keys:
        for low, original in lower.items():
            if key in low:
                return original
    return None


def _cell(row: dict[str, Any], column: str | None) -> str:
    if not column:
        return ""
    value = row.get(column, "")
    return "" if value is None else str(value).strip()


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


def split_genes(text: str) -> list[str]:
    """Separa IDs de gen por coma, punto y coma, barra o espacios."""
    return [g for g in re.split(r"[,;|\s]+", (text or "").strip()) if g]


def normalize_strand(value: Any) -> str:
    text = str(value or "").strip().lower()
    if text in {"+", "1", "plus", "fwd", "forward", "sentido"}:
        return "+"
    if text in {"-", "-1", "minus", "rev", "reverse", "antisentido"}:
        return "-"
    return "unknown"


def _is_k12(genome: str) -> bool:
    return (genome or "").strip().lower() in K12_NAMES


def _rel(root: Path, path: Path | None) -> str:
    if not path:
        return ""
    try:
        return Path(path).resolve().relative_to(Path(root).resolve()).as_posix()
    except ValueError:
        return Path(path).as_posix()


def infer_genome(path: Path | None) -> str:
    """Infiere el genoma desde el nombre/ruta de un artefacto local."""
    if not path:
        return ""
    name = Path(path).name
    if name in {"gt_operons.tsv", "gt_operons.tsv.gz"}:
        return "ecoli_k12"
    for suffix in (".operons_auto.tsv", ".promoters_auto.tsv",
                   "_gold_operons.tsv", "_operons_rules.tsv",
                   "_operons.tsv", ".gff3", ".gff", ".tsv"):
        if name.endswith(suffix):
            stem = name[:-len(suffix)]
            break
    else:
        stem = name
    if stem in {"operons", "gold_operons", "operon_benchmark"}:
        parent = Path(path).parent
        if parent.parent.name in {"regulation", "benchmarks", "gold"}:
            return parent.name
    return stem or "unknown"


# ---------------------------------------------------------------------------
# Lectura TSV
# ---------------------------------------------------------------------------
def load_tsv(path: Path) -> ParsedTable:
    path = Path(path)
    warnings: list[str] = []
    operons: list[Operon] = []
    all_genes: set[str] = set()
    n_rows = n_skip_map = n_no_genes = n_no_coords = n_unknown = 0
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields:
            raise ValueError(f"{path}: TSV sin cabecera")
        gcol = _detect(fields, GENE_KEYS)
        if not gcol:
            raise ValueError(
                f"{path}: no encuentro columna de genes entre {fields}")
        scol = _detect(fields, STRAND_KEYS)
        stcol = _detect(fields, START_KEYS)
        encol = _detect(fields, END_KEYS)
        ccol = _detect(fields, CONTIG_KEYS)
        icol = _detect(fields, ID_KEYS)
        mcol = _detect(fields, MAP_KEYS)
        for index, row in enumerate(reader, 1):
            n_rows += 1
            genes = set(split_genes(_cell(row, gcol)))
            if not genes:
                n_no_genes += 1
                continue
            # Igual que el benchmark histórico: el universo de genes incluye
            # todas las filas, pero solo se evalúan las de mapping_ok=yes.
            all_genes |= genes
            if mcol:
                map_value = _cell(row, mcol).lower()
                if map_value not in {"", "yes", "true", "1", "si", "sí"}:
                    n_skip_map += 1
                    continue
            start = _int(_cell(row, stcol))
            end = _int(_cell(row, encol))
            if start is None or end is None or start < 1 or end < 1 \
                    or end < start:
                start = end = None
                n_no_coords += 1
            strand = normalize_strand(_cell(row, scol))
            if strand == "unknown":
                n_unknown += 1
            operon_id = _cell(row, icol) or f"{path.stem}_{index}"
            operons.append(Operon(
                operon_id=operon_id,
                contig=_cell(row, ccol),
                strand=strand,
                start=start,
                end=end,
                genes=frozenset(genes),
            ))
    if n_skip_map:
        warnings.append(f"{n_skip_map} filas con mapping_ok distinto de yes "
                        "excluidas del gold evaluable.")
    if n_no_genes:
        warnings.append(f"{n_no_genes} filas sin lista de genes ignoradas.")
    if n_no_coords:
        warnings.append(f"{n_no_coords} filas sin coordenadas válidas: "
                        "solo entran en métricas por gen.")
    if n_unknown:
        warnings.append(f"{n_unknown} filas sin hebra +/− reconocible.")
    return ParsedTable(path=path, fmt="tsv", operons=operons,
                       all_genes=all_genes, warnings=warnings,
                       n_rows=n_rows, n_skipped_mapping=n_skip_map,
                       n_skipped_no_genes=n_no_genes,
                       n_without_coords=n_no_coords,
                       n_unknown_strand=n_unknown)


# ---------------------------------------------------------------------------
# Lectura GFF/GFF3
# ---------------------------------------------------------------------------
def parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in (text or "").split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[unquote(key.strip())] = unquote(value.strip())
    return attrs


def _attr(attrs: dict[str, str], keys: Iterable[str]) -> str:
    for key in keys:
        if attrs.get(key):
            return attrs[key]
    return ""


def _gff_gene_id(attrs: dict[str, str]) -> str:
    return _attr(attrs, ("gene_id", "gene", "ID", "Name", "locus_tag"))


def load_gff(path: Path) -> ParsedTable:
    path = Path(path)
    warnings: list[str] = []
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype, start_s, end_s, _score, strand_s, \
                _phase, attr_s = parts[:9]
            rows.append({
                "contig": seqid.strip(),
                "type": ftype.strip().lower(),
                "start": _int(start_s),
                "end": _int(end_s),
                "strand": normalize_strand(strand_s),
                "attrs": parse_attributes(attr_s),
            })

    # Genes y límites de hijos por Parent/operon_id, para GFF sin genes=.
    child_genes: dict[str, set[str]] = defaultdict(set)
    child_bounds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    parent_keys = ("Parent", "parent", "operon_id", "operon", "tu_id",
                   "transcription_unit")
    for row in rows:
        attrs = row["attrs"]
        for key in parent_keys:
            raw = attrs.get(key, "")
            if not raw:
                continue
            for parent in re.split(r"[,;]", raw):
                parent = parent.strip()
                if not parent:
                    continue
                gene = _gff_gene_id(attrs)
                if gene:
                    child_genes[parent].add(gene)
                if row["start"] is not None and row["end"] is not None:
                    child_bounds[parent].append((row["start"], row["end"]))

    op_rows = [r for r in rows if r["type"] in GFF_OPERON_TYPES]
    operons: list[Operon] = []
    n_no_coords = 0
    if op_rows:
        for index, row in enumerate(op_rows, 1):
            attrs = row["attrs"]
            genes: set[str] = set()
            for key in GFF_GENE_KEYS:
                if attrs.get(key):
                    genes |= set(split_genes(attrs[key]))
            operon_id = _attr(attrs, ("ID", "operon_id", "Name")) \
                or f"gff_operon_{index}"
            if not genes:
                genes = set(child_genes.get(operon_id, set()))
            if not genes:
                warnings.append(
                    f"feature de operón {operon_id} sin genes= ni hijos; "
                    "se ignora.")
                continue
            start, end = row["start"], row["end"]
            if (start is None or end is None) and child_bounds.get(operon_id):
                bounds = child_bounds[operon_id]
                start = min(b[0] for b in bounds)
                end = max(b[1] for b in bounds)
            if start is None or end is None:
                n_no_coords += 1
            strand = row["strand"]
            if strand == "unknown" and child_genes.get(operon_id):
                child_strands = {
                    r["strand"] for r in rows
                    if _gff_gene_id(r["attrs"]) in child_genes[operon_id]
                }
                child_strands.discard("unknown")
                if len(child_strands) == 1:
                    strand = child_strands.pop()
            operons.append(Operon(
                operon_id=operon_id, contig=row["contig"], strand=strand,
                start=start, end=end, genes=frozenset(genes)))
    else:
        # Sin features de operón: reconstruir desde CDS que declaren
        # explícitamente operon_id/tu_id. No se usa Parent genérico porque
        # suele apuntar al gen, no al operón.
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            op_id = _attr(row["attrs"], GFF_GROUP_KEYS)
            if not op_id:
                continue
            gene = _gff_gene_id(row["attrs"])
            if gene:
                groups[op_id].append({"row": row, "gene": gene})
        if not groups:
            raise ValueError(
                f"{path}: GFF sin features de operón ni atributos "
                "operon_id/tu_id en los CDS")
        for op_id, members in groups.items():
            genes = {m["gene"] for m in members}
            coords = [(m["row"]["start"], m["row"]["end"])
                      for m in members
                      if m["row"]["start"] is not None
                      and m["row"]["end"] is not None]
            if coords:
                start = min(c[0] for c in coords)
                end = max(c[1] for c in coords)
            else:
                start = end = None
                n_no_coords += 1
            strands = {m["row"]["strand"] for m in members}
            strands.discard("unknown")
            strand = strands.pop() if len(strands) == 1 else "unknown"
            contigs = {m["row"]["contig"] for m in members if m["row"]["contig"]}
            contig = contigs.pop() if len(contigs) == 1 else ""
            operons.append(Operon(
                operon_id=op_id, contig=contig, strand=strand,
                start=start, end=end, genes=frozenset(genes)))

    if not operons:
        raise ValueError(f"{path}: GFF sin operones utilizables")
    if n_no_coords:
        warnings.append(f"{n_no_coords} operones GFF sin coordenadas: solo "
                        "entran en métricas por gen.")
    n_unknown = sum(1 for o in operons if o.strand == "unknown")
    if n_unknown:
        warnings.append(f"{n_unknown} operones GFF sin hebra +/− reconocible.")
    all_genes = {g for o in operons for g in o.genes}
    return ParsedTable(path=path, fmt="gff", operons=operons,
                       all_genes=all_genes, warnings=warnings,
                       n_rows=len(operons), n_without_coords=n_no_coords,
                       n_unknown_strand=n_unknown)


def load_operons(path: Path) -> ParsedTable:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() in {".gff", ".gff3"}:
        return load_gff(path)
    return load_tsv(path)


# ---------------------------------------------------------------------------
# Alias de genes
# ---------------------------------------------------------------------------
def load_alias(path: Path | None) -> dict[str, str]:
    if not path or not Path(path).exists():
        return {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        ncol = _detect(fields, ("name", "gene_name", "alias"))
        bcol = _detect(fields, ("bnumber", "gene_id", "id"))
        if not ncol or not bcol:
            return {}
        return {row[ncol].strip(): row[bcol].strip()
                for row in reader
                if _cell(row, ncol) and _cell(row, bcol)}


def apply_alias(operons: list[Operon],
                alias: dict[str, str]) -> tuple[list[Operon], int]:
    if not alias:
        return operons, 0
    out: list[Operon] = []
    n_aliased = 0
    for op in operons:
        genes: set[str] = set()
        for gene in op.genes:
            mapped = alias.get(gene, gene)
            if mapped != gene:
                n_aliased += 1
            genes.add(mapped)
        out.append(replace(op, genes=frozenset(genes)))
    return out, n_aliased


# ---------------------------------------------------------------------------
# Métricas puras
# ---------------------------------------------------------------------------
def pairs(operons: Iterable[Iterable[str]]) -> set[frozenset[str]]:
    out: set[frozenset[str]] = set()
    for genes in operons:
        for a, b in combinations(sorted(set(genes)), 2):
            out.add(frozenset((a, b)))
    return out


def prf(tp: int, n_pred: int, n_gt: int) -> dict[str, Any]:
    precision = tp / n_pred if n_pred else 0.0
    recall = tp / n_gt if n_gt else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if precision + recall else 0.0)
    return {
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "tp": tp,
        "n_pred": n_pred,
        "n_gt": n_gt,
    }


def _evaluated(tp: int, n_pred: int, n_gt: int,
               **extra: Any) -> dict[str, Any]:
    metric = {"status": "evaluated", "reason": ""}
    metric.update(prf(tp, n_pred, n_gt))
    metric.update(extra)
    return metric


def _not_evaluable(reason: str, n_pred: int = 0,
                   n_gt: int = 0, **extra: Any) -> dict[str, Any]:
    metric = {
        "status": "not_evaluable", "reason": reason,
        "precision": None, "recall": None, "f1": None,
        "tp": 0, "n_pred": n_pred, "n_gt": n_gt,
    }
    metric.update(extra)
    return metric


def reciprocal_overlap(a: Operon, b: Operon) -> float:
    """Mínimo de las fracciones de a y b cubiertas por su intersección."""
    if not a.has_coords or not b.has_coords:
        return 0.0
    inter = min(a.end, b.end) - max(a.start, b.start) + 1
    if inter <= 0:
        return 0.0
    len_a = a.end - a.start + 1
    len_b = b.end - b.start + 1
    return min(inter / len_a, inter / len_b)


def _compatible(a: Operon, b: Operon) -> bool:
    if a.strand != b.strand:
        return False
    if a.contig and b.contig and a.contig != b.contig:
        return False
    return a.has_coords and b.has_coords


def match_intervals(pred: list[Operon], gt: list[Operon], *,
                    exact: bool, min_overlap: float) -> int:
    """Emparejamiento 1:1 determinista por solapamiento o igualdad exacta."""
    candidates: list[tuple[float, str, str, int, int]] = []
    for i, p in enumerate(pred):
        for j, g in enumerate(gt):
            if not _compatible(p, g):
                continue
            if exact:
                if p.start != g.start or p.end != g.end:
                    continue
                score = 1.0
            else:
                score = reciprocal_overlap(p, g)
                if score + 1e-12 < min_overlap:
                    continue
            candidates.append((score, p.operon_id, g.operon_id, i, j))
    candidates.sort(key=lambda x: (-x[0], x[1], x[2]))
    used_p: set[int] = set()
    used_g: set[int] = set()
    matched = 0
    for _score, _pid, _gid, i, j in candidates:
        if i in used_p or j in used_g:
            continue
        used_p.add(i)
        used_g.add(j)
        matched += 1
    return matched


def gene_pair_metric(gt: list[Operon],
                     pred: list[Operon]) -> dict[str, Any]:
    if not gt or not pred:
        return _not_evaluable("sin_operones_para_el_nivel_gen",
                              n_pred=len(pred), n_gt=len(gt))
    gt_pairs = pairs(o.genes for o in gt)
    pred_pairs = pairs(o.genes for o in pred)
    return _evaluated(len(gt_pairs & pred_pairs), len(pred_pairs),
                      len(gt_pairs))


def gene_exact_metric(gt: list[Operon],
                      pred: list[Operon]) -> dict[str, Any]:
    if not gt or not pred:
        return _not_evaluable("sin_operones_para_el_nivel_gen",
                              n_pred=len(pred), n_gt=len(gt))
    gt_sets = {frozenset(o.genes) for o in gt}
    pred_sets = {frozenset(o.genes) for o in pred}
    exact = len(gt_sets & pred_sets)
    metric = _evaluated(exact, len(pred_sets), len(gt_sets))
    metric["n_exact"] = exact
    metric["exact_rate"] = round(exact / len(gt_sets), 4) if gt_sets else None
    return metric


def interval_metric(gt: list[Operon], pred: list[Operon], *,
                    exact: bool, min_overlap: float) -> dict[str, Any]:
    gt_c = [o for o in gt if o.has_coords]
    pred_c = [o for o in pred if o.has_coords]
    if not gt_c or not pred_c:
        if not gt:
            reason = "sin_operones_gold_en_la_hebra"
        elif not pred:
            reason = "sin_operones_predichos_en_la_hebra"
        elif not gt_c:
            reason = "sin_intervalos_gold_con_coordenadas"
        else:
            reason = "sin_intervalos_predichos_con_coordenadas"
        return _not_evaluable(
            reason,
            n_pred=len(pred_c), n_gt=len(gt_c),
            n_pred_without_coords=len(pred) - len(pred_c),
            n_gt_without_coords=len(gt) - len(gt_c))
    matched = match_intervals(pred_c, gt_c, exact=exact,
                              min_overlap=min_overlap)
    metric = _evaluated(matched, len(pred_c), len(gt_c))
    metric["n_matched"] = matched
    metric["n_pred_without_coords"] = len(pred) - len(pred_c)
    metric["n_gt_without_coords"] = len(gt) - len(gt_c)
    metric["min_overlap"] = 1.0 if exact else round(min_overlap, 4)
    if exact:
        metric["n_exact"] = matched
        metric["exact_rate"] = round(matched / len(gt_c), 4)
    return metric


def _bucket(gt: list[Operon], pred_gene: list[Operon],
            pred_interval: list[Operon], min_overlap: float) -> dict[str, Any]:
    gt_genes = {g for o in gt for g in o.genes}
    pred_genes = {g for o in pred_gene for g in o.genes}
    shared = gt_genes & pred_genes
    if not gt:
        gene_reason = "sin_operones_gold_en_la_hebra"
    elif not pred_gene:
        gene_reason = "sin_operones_predichos_con_genes_mapeados"
    elif not shared:
        gene_reason = "sin_ids_de_gen_compartidos_en_la_hebra"
    else:
        gene_reason = ""
    gene_pair = (gene_pair_metric(gt, pred_gene) if not gene_reason
                 else _not_evaluable(gene_reason, n_pred=len(pred_gene),
                                     n_gt=len(gt)))
    gene_exact = (gene_exact_metric(gt, pred_gene) if not gene_reason
                  else _not_evaluable(gene_reason, n_pred=len(pred_gene),
                                      n_gt=len(gt)))
    interval_overlap = interval_metric(gt, pred_interval, exact=False,
                                       min_overlap=min_overlap)
    interval_exact = interval_metric(gt, pred_interval, exact=True,
                                     min_overlap=min_overlap)
    statuses = [gene_pair["status"], gene_exact["status"],
                interval_overlap["status"], interval_exact["status"]]
    if all(s == "evaluated" for s in statuses):
        status = "evaluated"
    elif any(s == "evaluated" for s in statuses):
        status = "partial"
    else:
        status = "not_evaluable"
    return {
        "status": status,
        "n_gt_operons": len(gt),
        "n_pred_operons": len(pred_interval),
        "n_shared_genes": len(shared),
        "gene_pair": gene_pair,
        "gene_exact": gene_exact,
        "interval_overlap": interval_overlap,
        "interval_exact": interval_exact,
    }


# ---------------------------------------------------------------------------
# Evaluación
# ---------------------------------------------------------------------------
def _base_result(genome: str, gold: Path | None, pred: Path | None,
                 *, min_genes: int, min_overlap: float,
                 alias: Path | None, root: Path) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "status": "not_evaluable",
        "reason_code": "",
        "reason": "",
        "evaluation_level": "none",
        "genome": genome or "unknown",
        "gold": {"path": _rel(root, gold), "fmt": "", "searched": []},
        "pred": {"path": _rel(root, pred), "fmt": "", "searched": []},
        "parameters": {
            "min_genes": int(min_genes),
            "min_overlap": float(min_overlap),
            "alias": _rel(root, alias) if alias else "",
        },
        "counts": {},
        "by_strand": {},
        "warnings": [],
        "promoters": {
            "status": "not_evaluable",
            "reason": ("Este benchmark mide estructura de operones. No se "
                       "compara promoters_auto con gt_promoters porque una "
                       "ventana upstream no es un TSS y no se define aquí un "
                       "emparejamiento TSS→ventana; no se inventa gold."),
        },
    }


def _not_evaluable_result(base: dict[str, Any], code: str,
                          reason: str) -> dict[str, Any]:
    base["status"] = "not_evaluable"
    base["reason_code"] = code
    base["reason"] = reason
    return base


def evaluate(gold_path: Path, pred_path: Path, *, genome: str = "",
             min_genes: int = 2, min_overlap: float = 0.5,
             alias_path: Path | None = None,
             root: Path = ROOT,
             searched_gold: list[str] | None = None,
             searched_pred: list[str] | None = None) -> dict[str, Any]:
    """Evalúa una predicción local contra un gold local. Sin red."""
    gold_path = Path(gold_path)
    pred_path = Path(pred_path)
    base = _base_result(genome, gold_path, pred_path, min_genes=min_genes,
                        min_overlap=min_overlap, alias=alias_path, root=root)
    base["gold"]["searched"] = searched_gold or []
    base["pred"]["searched"] = searched_pred or []
    try:
        gold_table = load_operons(gold_path)
    except Exception as exc:  # noqa: BLE001
        return _not_evaluable_result(base, "gold_unreadable",
                                     f"No se pudo leer el gold: {exc}")
    try:
        pred_table = load_operons(pred_path)
    except Exception as exc:  # noqa: BLE001
        return _not_evaluable_result(base, "prediction_unreadable",
                                     f"No se pudo leer la predicción: {exc}")
    base["gold"]["fmt"] = gold_table.fmt
    base["pred"]["fmt"] = pred_table.fmt
    base["warnings"] = list(gold_table.warnings) + list(pred_table.warnings)

    if alias_path is None:
        candidate = gold_path.parent / "gt_gene_aliases.tsv"
        alias_path = candidate if candidate.exists() else None
    if alias_path and not Path(alias_path).exists():
        base["warnings"].append(f"Alias indicado pero no encontrado: "
                                f"{alias_path}")
    alias = load_alias(alias_path)
    base["parameters"]["alias"] = _rel(root, alias_path) if alias_path else ""
    pred_ops, n_aliased = apply_alias(pred_table.operons, alias)

    gold_all_genes = set(gold_table.all_genes)
    pred_all_genes = {g for o in pred_ops for g in o.genes}
    shared_genes = gold_all_genes & pred_all_genes

    gold_ops = [o for o in gold_table.operons if len(o.genes) >= min_genes]
    pred_ops_filtered = [o for o in pred_ops if len(o.genes) >= min_genes]

    pred_gene_ops: list[Operon] = []
    unmapped: Counter[str] = Counter()
    for op in pred_ops_filtered:
        mapped = op.genes & gold_all_genes
        for gene in op.genes - gold_all_genes:
            unmapped[gene] += 1
        if len(mapped) >= min_genes:
            pred_gene_ops.append(replace(op, genes=frozenset(mapped)))

    if not gold_ops:
        reason = (f"El gold no tiene operones con >= {min_genes} genes "
                  f"tras aplicar mapping_ok (filas leídas: "
                  f"{gold_table.n_rows}).")
        return _not_evaluable_result(base, "gold_without_evaluable_operons",
                                     reason)
    if not pred_ops_filtered:
        reason = (f"La predicción no tiene operones con >= {min_genes} "
                  f"genes (filas leídas: {pred_table.n_rows}).")
        return _not_evaluable_result(base, "prediction_without_operons",
                                     reason)

    predicate = {
        "all": lambda o: True,
        "+": lambda o: o.strand == "+",
        "-": lambda o: o.strand == "-",
    }
    pred_int_genes = {g for o in pred_ops_filtered for g in o.genes}
    shared_filtered = pred_int_genes & gold_all_genes
    base["counts"] = {
        "gt_rows": gold_table.n_rows,
        "pred_rows": pred_table.n_rows,
        "gt_operons_raw": len(gold_table.operons),
        "pred_operons_raw": len(pred_table.operons),
        "gt_operons_filtered": len(gold_ops),
        "pred_operons_filtered": len(pred_ops_filtered),
        "pred_operons_with_mapped_genes": len(pred_gene_ops),
        "gt_genes": len(gold_all_genes),
        "pred_genes": len(pred_all_genes),
        "pred_genes_filtered": len(pred_int_genes),
        "shared_genes": len(shared_genes),
        "shared_genes_filtered": len(shared_filtered),
        "pred_gene_mapping_coverage": (
            round(len(shared_filtered) / len(pred_int_genes), 4)
            if pred_int_genes else None),
        "pred_gene_tokens_unmapped": len(unmapped),
        "pred_genes_unmapped": len(pred_int_genes - gold_all_genes),
        "aliased_pred_genes": n_aliased,
        "gt_without_coords": sum(1 for o in gold_ops if not o.has_coords),
        "pred_without_coords": sum(1 for o in pred_ops_filtered
                                   if not o.has_coords),
        "gt_skipped_mapping": gold_table.n_skipped_mapping,
        "pred_skipped_mapping": pred_table.n_skipped_mapping,
    }
    if not shared_filtered:
        base["warnings"].append(
            "Sin IDs de gen compartidos: las métricas por gen quedan "
            "not_evaluable; las de intervalo siguen siendo evaluables.")
    elif len(shared_filtered) < len(pred_int_genes):
        base["warnings"].append(
            f"{len(pred_int_genes) - len(shared_filtered)} genes de la "
            "predicción no están en el universo gold y se excluyen de las "
            "métricas por gen.")
    if gold_table.n_without_coords or pred_table.n_without_coords:
        base["warnings"].append(
            "Hay operones sin coordenadas; no entran en las métricas de "
            "intervalo.")
    if any(not g.contig for g in gold_ops):
        base["warnings"].append(
            "El gold no declara contig en algunas filas: el contig se trata "
            "como comodín en el emparejamiento de intervalos.")

    for name, fn in predicate.items():
        gt_b = [o for o in gold_ops if fn(o)]
        pred_gene_b = [o for o in pred_gene_ops if fn(o)]
        pred_int_b = [o for o in pred_ops_filtered if fn(o)]
        base["by_strand"][name] = _bucket(gt_b, pred_gene_b, pred_int_b,
                                          min_overlap)
    if any(o.strand == "unknown" for o in gold_ops + pred_ops_filtered):
        base["by_strand"]["unknown"] = _bucket(
            [o for o in gold_ops if o.strand == "unknown"],
            [o for o in pred_gene_ops if o.strand == "unknown"],
            [o for o in pred_ops_filtered if o.strand == "unknown"],
            min_overlap)

    all_bucket = base["by_strand"]["all"]
    gene_ok = all_bucket["gene_pair"]["status"] == "evaluated"
    interval_ok = all_bucket["interval_overlap"]["status"] == "evaluated"
    if gene_ok and interval_ok:
        base["evaluation_level"] = "full"
    elif interval_ok:
        base["evaluation_level"] = "interval_only"
    elif gene_ok:
        base["evaluation_level"] = "gene_only"
    else:
        base["evaluation_level"] = "none"
    base["status"] = "evaluated"
    base["reason_code"] = ""
    base["reason"] = ""
    base["counts"]["top_unmapped_genes"] = [
        {"gene": gene, "n": n}
        for gene, n in unmapped.most_common(10)]
    return base


# ---------------------------------------------------------------------------
# Autodetección local (sin red)
# ---------------------------------------------------------------------------
def gold_candidates(root: Path, genome: str) -> list[Path]:
    root = Path(root)
    out: list[Path] = []
    if genome:
        out.extend([
            root / "data" / "regulation" / genome / "gold_operons.tsv",
            root / "data" / "regulation" / genome / "operons_gold.tsv",
            root / "data" / "benchmarks" / f"{genome}_gold_operons.tsv",
            root / "data" / "benchmarks" / genome / "gold_operons.tsv",
            root / "data" / "gold" / f"{genome}_operons.tsv",
            root / "data" / "gold" / f"{genome}_gold_operons.tsv",
            root / f"{genome}_gold_operons.tsv",
        ])
        if _is_k12(genome):
            out.append(root / "data" / "regulondb" / "gt" / "gt_operons.tsv")
    return out


def find_gold(root: Path, genome: str) -> tuple[Path | None, list[str]]:
    candidates = gold_candidates(root, genome)
    searched = [_rel(root, p) for p in candidates]
    for path in candidates:
        if path.exists():
            return path, searched
    return None, searched


def pred_candidates(root: Path, genome: str) -> list[Path]:
    root = Path(root)
    if not genome:
        return []
    return [
        root / f"{genome}.operons_auto.tsv",
        root / "results" / f"{genome}.operons_auto.tsv",
        root / "out" / f"{genome}.operons_auto.tsv",
        root / "out_fl" / f"{genome}.operons_auto.tsv",
        root / f"{genome}.gff3",
        root / "results" / f"{genome}.gff3",
        root / "out" / f"{genome}.gff3",
        root / "data" / "operon_finder" / f"{genome}_operons_rules.tsv",
        root / "data" / "operon_finder" / f"{genome}_operons.tsv",
        root / "data" / "regulation" / genome / "operons.tsv",
    ]


def find_pred(root: Path, genome: str) -> tuple[Path | None, list[str]]:
    root = Path(root)
    candidates = pred_candidates(root, genome)
    searched = [_rel(root, p) for p in candidates]
    for path in candidates:
        if path.exists():
            return path, searched
    if not genome:
        found: list[Path] = []
        for directory in (root, root / "results", root / "out",
                          root / "out_fl"):
            if directory.exists():
                found.extend(sorted(directory.glob("*.operons_auto.tsv")))
        if len(found) == 1:
            return found[0], searched
    return None, searched


# ---------------------------------------------------------------------------
# Informe
# ---------------------------------------------------------------------------
def _metric_cell(metric: dict[str, Any], key: str) -> str:
    if metric.get("status") != "evaluated":
        return "—"
    value = metric.get(key)
    return "—" if value is None else f"{value:.3f}" if isinstance(value, float) \
        else str(value)


def _detail_table(bucket: dict[str, Any]) -> list[str]:
    rows = [
        ("Pares de genes", bucket["gene_pair"]),
        ("Conjuntos exactos", bucket["gene_exact"]),
        ("Intervalos solapados", bucket["interval_overlap"]),
        ("Intervalos exactos", bucket["interval_exact"]),
    ]
    lines = [
        "| métrica | P | R | F1 | TP | n_pred | n_gt | estado |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for label, metric in rows:
        if metric["status"] == "evaluated":
            lines.append(
                f"| {label} | {metric['precision']:.3f} | "
                f"{metric['recall']:.3f} | {metric['f1']:.3f} | "
                f"{metric['tp']} | {metric['n_pred']} | {metric['n_gt']} | "
                f"evaluated |")
        else:
            lines.append(
                f"| {label} | — | — | — | — | {metric.get('n_pred', 0)} | "
                f"{metric.get('n_gt', 0)} | not_evaluable: "
                f"{metric.get('reason', '')} |")
    if bucket["gene_exact"]["status"] == "evaluated":
        ge = bucket["gene_exact"]
        lines.append("")
        lines.append(f"- exact-set: {ge['n_exact']} de {ge['n_gt']} "
                     f"operones gold ({ge['exact_rate']:.1%}).")
    if bucket["interval_exact"]["status"] == "evaluated":
        ie = bucket["interval_exact"]
        lines.append(f"- exact-intervalo: {ie['n_exact']} de {ie['n_gt']} "
                     f"operones gold ({ie['exact_rate']:.1%}).")
    return lines


def render_markdown(result: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("# Benchmark de operones por reglas")
    lines.append("")
    lines.append(f"**Generado:** {result.get('generated_at', '')}")
    lines.append("")
    lines.append("**Alcance:** solo datos locales; sin red, sin reentrenar y "
                 "sin modificar pipeline/GUI. No se inventa gold.")
    lines.append("")
    if result.get("command"):
        lines.append(f"**Comando:** `{result['command']}`")
        lines.append("")
    lines.append(f"**Genoma:** `{result.get('genome', '')}`")
    lines.append("")
    if result["status"] == "not_evaluable":
        lines.append("## Estado: no evaluable")
        lines.append("")
        lines.append(f"- **Código:** `{result.get('reason_code', '')}`")
        lines.append(f"- **Razón:** {result.get('reason', '')}")
        searched = (result.get("gold", {}).get("searched") or [])
        if searched:
            lines.append("- **Rutas de gold buscadas:** "
                         + ", ".join(f"`{p}`" for p in searched))
        searched_pred = (result.get("pred", {}).get("searched") or [])
        if searched_pred:
            lines.append("- **Rutas de predicción buscadas:** "
                         + ", ".join(f"`{p}`" for p in searched_pred))
        lines.append("")
        lines.append("Un `not_evaluable` no es un cero de rendimiento: faltan "
                     "los datos para medir, no se supone un resultado.")
        return "\n".join(lines) + "\n"

    lines.append("## Entradas")
    lines.append("")
    lines.append("| artefacto | formato | ruta |")
    lines.append("|---|---|---|")
    lines.append(f"| gold | {result['gold'].get('fmt', '')} | "
                 f"`{result['gold'].get('path', '')}` |")
    lines.append(f"| predicción | {result['pred'].get('fmt', '')} | "
                 f"`{result['pred'].get('path', '')}` |")
    lines.append("")
    lines.append(f"**Estado:** evaluado (`evaluation_level="
                 f"{result.get('evaluation_level', '')}`).")
    lines.append("")
    counts = result.get("counts", {})
    lines.append("## Conteos")
    lines.append("")
    lines.append("| campo | valor |")
    lines.append("|---|---:|")
    for key in ("gt_rows", "pred_rows", "gt_operons_filtered",
                "pred_operons_filtered", "pred_operons_with_mapped_genes",
                "gt_genes", "pred_genes", "pred_genes_filtered",
                "shared_genes", "shared_genes_filtered",
                "pred_gene_mapping_coverage", "gt_without_coords",
                "pred_without_coords", "aliased_pred_genes"):
        if key in counts:
            lines.append(f"| `{key}` | {counts[key]} |")
    lines.append("")
    lines.append("## Métricas por hebra")
    lines.append("")
    lines.append("| hebra | nivel gen P | R | F1 | exact-set | "
                 "nivel intervalo P | R | F1 | exact-intervalo |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for strand in ("all", "+", "-"):
        bucket = result["by_strand"].get(strand)
        if not bucket:
            continue
        gp = bucket["gene_pair"]
        ge = bucket["gene_exact"]
        io = bucket["interval_overlap"]
        ie = bucket["interval_exact"]
        lines.append(
            f"| `{strand}` | {_metric_cell(gp, 'precision')} | "
            f"{_metric_cell(gp, 'recall')} | {_metric_cell(gp, 'f1')} | "
            f"{_metric_cell(ge, 'exact_rate')} | "
            f"{_metric_cell(io, 'precision')} | "
            f"{_metric_cell(io, 'recall')} | {_metric_cell(io, 'f1')} | "
            f"{_metric_cell(ie, 'exact_rate')} |")
    unknown = result["by_strand"].get("unknown")
    if unknown:
        lines.append(
            f"| `unknown` | — | — | "
            f"{_metric_cell(unknown['gene_pair'], 'f1')} | — | — | — | "
            f"{_metric_cell(unknown['interval_overlap'], 'f1')} | — |")
    lines.append("")
    for strand in ("all", "+", "-", "unknown"):
        bucket = result["by_strand"].get(strand)
        if not bucket:
            continue
        lines.append(f"### Hebra `{strand}`")
        lines.append("")
        lines.extend(_detail_table(bucket))
        lines.append("")
    if result.get("warnings"):
        lines.append("## Avisos")
        lines.append("")
        for warning in result["warnings"]:
            lines.append(f"- {warning}")
        lines.append("")
    lines.append("## Promotores")
    lines.append("")
    promoters = result.get("promoters") or {}
    if promoters.get("status") == "not_evaluable":
        lines.append("- **Estado:** `not_evaluable`.")
        lines.append(f"- **Razón:** {promoters.get('reason', '')}")
    else:
        lines.append("- Estado: " + str(promoters.get("status", "n/d")))
    lines.append("")
    lines.append("## Límites")
    lines.append("")
    lines.append("- No es un holdout ciego: las reglas de K-12 pueden haberse "
                 "calibrado con el mismo gold; mide consistencia, no "
                 "generalización.")
    lines.append("- Nivel gen: los pares solo existen dentro de operones "
                 f"con ≥ {result['parameters']['min_genes']} genes; IDs no "
                 "mapeados se excluyen y eso puede favorecer la precisión.")
    lines.append("- Nivel intervalo: solapamiento recíproco mínimo "
                 f"{result['parameters']['min_overlap']}; coordenadas exactas "
                 "son una métrica aparte y muy sensible a los bordes.")
    lines.append("- Si el gold no declara contig, el contig se trata como "
                 "comodín. Hebras `+` y `-` se evalúan por separado.")
    lines.append("- Promotores/TSS no se evalúan aquí: una ventana upstream no "
                 "es un TSS y no se define un emparejamiento TSS→ventana.")
    lines.append("- Un `not_evaluable` no es cero: indica que falta gold, "
                 "predicción o IDs/coordenadas comparables.")
    return "\n".join(lines) + "\n"


def write_reports(result: dict[str, Any], json_path: Path | None,
                  md_path: Path | None) -> None:
    if json_path:
        json_path = Path(json_path)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
    if md_path:
        md_path = Path(md_path)
        md_path.parent.mkdir(parents=True, exist_ok=True)
        md_path.write_text(render_markdown(result), encoding="utf-8")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _print_summary(result: dict[str, Any]) -> None:
    if result["status"] != "evaluated":
        print(f"NOT_EVALUABLE [{result.get('reason_code', '')}]: "
              f"{result.get('reason', '')}")
        return
    print(f"\n=== Operones vs gold: {result.get('genome', '')} ===")
    print(f"  gold: {result['gold'].get('path', '')} "
          f"({result['counts'].get('gt_operons_filtered', 0)} operones)")
    print(f"  pred: {result['pred'].get('path', '')} "
          f"({result['counts'].get('pred_operons_filtered', 0)} operones)")
    print(f"  evaluation_level={result.get('evaluation_level', '')}")
    for strand in ("all", "+", "-"):
        bucket = result["by_strand"].get(strand)
        if not bucket:
            continue
        gp = bucket["gene_pair"]
        io = bucket["interval_overlap"]
        ie = bucket["interval_exact"]
        if gp["status"] == "evaluated":
            gene = f"gen P={gp['precision']:.3f} R={gp['recall']:.3f} F1={gp['f1']:.3f}"
        else:
            gene = f"gen not_evaluable ({gp['reason']})"
        if io["status"] == "evaluated":
            interval = (f"intervalo P={io['precision']:.3f} "
                        f"R={io['recall']:.3f} F1={io['f1']:.3f} "
                        f"exact={ie.get('exact_rate')}")
        else:
            interval = f"intervalo not_evaluable ({io['reason']})"
        print(f"  {strand:3s}  {gene}  |  {interval}")
    if result.get("warnings"):
        print(f"  avisos: {len(result['warnings'])}")


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Benchmark local de operones (sin red, sin inventar gold).")
    ap.add_argument("--gold", type=Path, default=None,
                    help="TSV/GFF de gold. Si falta, se autodetecta.")
    ap.add_argument("--pred", type=Path, default=None,
                    help="TSV de operones predichos o GFF/GFF3. Si falta, se "
                         "autodetecta (<genoma>.operons_auto.tsv, GFF, reglas).")
    ap.add_argument("--genome", default="",
                    help="Nombre del genoma para autodetectar artefactos.")
    ap.add_argument("--root", type=Path, default=ROOT,
                    help="Raíz del proyecto para autodetección (default: repo).")
    ap.add_argument("--min-genes", type=int, default=2,
                    help="Mínimo de genes por operón evaluado (default 2).")
    ap.add_argument("--min-overlap", type=float, default=0.5,
                    help="Solapamiento recíproco mínimo 0-1 para nivel "
                         "intervalo (default 0.5).")
    ap.add_argument("--alias", type=Path, default=None,
                    help="TSV name→bnumber; si falta, usa gt_gene_aliases.tsv "
                         "junto al gold cuando exista.")
    ap.add_argument("--json", type=Path, default=DEFAULT_JSON,
                    help="Salida JSON (default reports/operon_benchmark.json).")
    ap.add_argument("--md", type=Path, default=DEFAULT_MD,
                    help="Salida Markdown (default reports/operon_benchmark.md).")
    ap.add_argument("--quiet", action="store_true", help="No imprimir resumen.")
    raw_argv = list(sys.argv[1:]) if argv is None else list(argv)
    args = ap.parse_args(raw_argv)

    if args.min_genes < 1:
        ap.error("--min-genes debe ser >= 1")
    if not 0 < args.min_overlap <= 1:
        ap.error("--min-overlap debe estar en (0, 1]")

    root = Path(args.root).resolve()
    genome = (args.genome or "").strip()
    gold = Path(args.gold) if args.gold else None
    pred = Path(args.pred) if args.pred else None
    if not genome and pred:
        genome = infer_genome(pred)
    if not genome and gold:
        genome = infer_genome(gold)

    searched_gold: list[str] = []
    searched_pred: list[str] = []
    if pred is None:
        pred, searched_pred = find_pred(root, genome)
        if pred and not genome:
            genome = infer_genome(pred)
    if gold is None:
        gold, searched_gold = find_gold(root, genome)

    base = _base_result(genome, gold, pred, min_genes=args.min_genes,
                        min_overlap=args.min_overlap, alias=args.alias,
                        root=root)
    base["gold"]["searched"] = searched_gold
    base["pred"]["searched"] = searched_pred

    if gold is None:
        result = _not_evaluable_result(
            base, "gold_not_found",
            "No se encontró un gold local aplicable. Se buscaron: "
            + (", ".join(searched_gold) or "(sin candidatos para este genoma)")
            + ". No se usa red ni se inventa gold.")
    elif not Path(gold).exists():
        result = _not_evaluable_result(
            base, "gold_not_found",
            f"No existe el gold indicado: {gold}. No se usa red ni se "
            "inventa gold.")
    elif pred is None:
        result = _not_evaluable_result(
            base, "prediction_not_found",
            "No se encontró una predicción local. Se buscaron: "
            + (", ".join(searched_pred) or "(sin candidatos para este genoma)")
            + ". Usa --pred para indicarla.")
    elif not Path(pred).exists():
        result = _not_evaluable_result(
            base, "prediction_not_found",
            f"No existe la predicción indicada: {pred}.")
    else:
        result = evaluate(gold, pred, genome=genome,
                          min_genes=args.min_genes,
                          min_overlap=args.min_overlap,
                          alias_path=args.alias, root=root,
                          searched_gold=searched_gold,
                          searched_pred=searched_pred)

    command = "python tools/benchmark_operons.py"
    if raw_argv:
        command += " " + " ".join(str(a) for a in raw_argv)
    result["command"] = command
    write_reports(result, args.json, args.md)
    if not args.quiet:
        _print_summary(result)
        if args.json:
            print(f"  json -> {args.json}")
        if args.md:
            print(f"  md   -> {args.md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
