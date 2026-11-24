#!/usr/bin/env python3
"""Calibra los clusters del modelo K-12 sobre los 400 genomas del proyecto.

El insumo principal son los ``data/complexes/model_<genoma>.complexes.tsv``
generados por ``tools/run_complex_model_400.py``. Cada fila es una *sospecha*
(``complex_candidate``), no una interacción demostrada. Este script **no**
cierra brechas, **no** modifica el pipeline, la GUI ni los datos y **no usa
red**. Solo agrega, describe solapamientos y prioriza revisiones.

Qué calcula
-----------
1. Agregados por genoma, confianza (alta/media/baja), tamaño de cluster y
   tokens de evidencia.
2. Solapamiento con las transferencias K-12 ``transfer_<genoma>.tsv`` en dos
   planos independientes:
   - **gen**: intersección/Jaccard/cobertura entre los genes del cluster y los
     genes del genoma listados en la transferencia;
   - **intervalo**: fracción del intervalo del cluster cubierta por el
     intervalo derivado de los genes de la transferencia (mismo contig).
3. Solapamiento con las firmas K-12 ``clusters_ecoli_k12_evidencia.tsv``:
   - se construye un **puente K-12** con las filas de transferencia que
     comparten al menos un gen con el cluster;
   - **gen**: Jaccard entre el conjunto de genes K-12 de la transferencia y el
     conjunto de genes de cada firma;
   - **intervalo**: fracción del intervalo K-12 de la transferencia cubierta
     por el intervalo K-12 de la firma (coordenadas de ``k12_genes.tsv``).
   El puente es indirecto (RBH/complejo), no ortología demostrada.
4. Una tabla de priorización heurística: confianza + solapamiento genético con
   transferencia + tamaño (+ señales secundarias de firma e intervalo).

Los denominadores y pesos son explícitos en el JSON (``definiciones`` y
``priorizacion.pesos``). Ningún número de este informe es precisión, recall ni
probabilidad calibrada.

Uso
---
    .venv/bin/python tools/calibrate_model_clusters.py
    .venv/bin/python tools/calibrate_model_clusters.py --limit-genomes 5
    .venv/bin/python tools/calibrate_model_clusters.py \
        --out-json reports/model_clusters_calibration.json \
        --out-md reports/model_clusters_calibration.md
"""
from __future__ import annotations

import argparse
import csv
import glob
import heapq
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Sequence

ROOT = Path(__file__).resolve().parents[1]
MODEL_PREFIX = "model_"
MODEL_SUFFIX = ".complexes.tsv"
CONF_ORDER = ("alta", "media", "baja")
CONF_WEIGHT = {"alta": 2.0, "media": 1.0, "baja": 0.0}
DEFAULT_MAX_COLOCAL_SPAN = 50_000
SIZE_BINS = (
    (2, 2, "2"),
    (3, 3, "3"),
    (4, 4, "4"),
    (5, 6, "5-6"),
    (7, 9, "7-9"),
    (10, None, "10+"),
)
FLAG_KEYS = (
    "has_transfer_gene",
    "has_transfer_iv",
    "has_transfer_iv_any",
    "strong_transfer_gene",
    "strong_transfer_iv",
    "strong_transfer_iv_any",
    "has_signature_gene",
    "has_signature_iv",
    "has_bridge",
    "has_bridge_coloc",
    "bridge_signature_backed",
    "sin_genes_resueltos",
)
METRIC_KEYS = (
    "t_gene_cov",
    "t_gene_jac",
    "t_iv_frac",
    "t_iv_frac_any",
    "s_gene_jac",
    "s_iv_frac",
    "priority_score",
    "n_resolved",
)


# --------------------------------------------------------------------------- #
# Utilidades de lectura y normalización
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


def parse_gene_ids(text: object) -> list[str]:
    """Extrae identificadores de genes de ``subunits``/``genes``.

    Acepta ``nombre:gen;nombre:gen`` y ``gen,gen``. Si hay ``:`` se conserva
    la parte derecha (el identificador), tolerando nombres con ``:``.
    """
    out: list[str] = []
    seen: set[str] = set()
    for chunk in re.split(r"[;,]", str(text or "")):
        chunk = chunk.strip()
        if not chunk:
            continue
        gid = chunk.rsplit(":", 1)[-1].strip() if ":" in chunk else chunk
        if gid and gid not in seen:
            seen.add(gid)
            out.append(gid)
    return out


def normalize_confidence(value: object) -> str:
    conf = str(value or "").strip().lower()
    return conf if conf in CONF_ORDER else (conf or "desconocida")


def size_bin_label(size: int) -> str:
    for lo, hi, label in SIZE_BINS:
        if size >= lo and (hi is None or size <= hi):
            return label
    return "otro"


def jaccard(a: set[str] | frozenset[str], b: set[str] | frozenset[str]) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def interval_overlap_fraction(a_start: int, a_end: int,
                              b_start: int, b_end: int) -> float:
    """Fracción del intervalo A cubierta por B.

    Se usa el intervalo semiabierto [start, end): dos genes adyacentes
    (200 y 201) no se consideran solapados, y una cobertura exacta da 1,0.
    """
    if a_start <= 0 or a_end <= a_start or b_start <= 0 or b_end <= b_start:
        return 0.0
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    length = a_end - a_start
    return overlap / length if length > 0 else 0.0


def evidence_tokens(text: object) -> list[str]:
    """Normaliza la evidencia del modelo en tokens contables.

    Los valores numéricos de ``score medio``/``score mínimo`` se colapsan para
    no fragmentar el conteo; el resto se conserva tal cual.
    """
    out: list[str] = []
    seen: set[str] = set()
    for raw in re.split(r"[;]", str(text or "")):
        token = raw.strip()
        if not token:
            continue
        low = token.lower()
        if low.startswith("score medio"):
            token = "score medio"
        elif low.startswith("score mínimo") or low.startswith("score minimo"):
            token = "score mínimo"
        elif low.startswith("modelo k-12"):
            token = "modelo K-12"
        if token not in seen:
            seen.add(token)
            out.append(token)
    return out


# --------------------------------------------------------------------------- #
# Estructuras
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GeneCoord:
    contig: str
    start: int
    end: int


@dataclass
class ModelCluster:
    genome: str
    complex_id: str
    name: str
    contig: str
    start: int
    end: int
    confidence: str
    evidence: str
    status: str
    genes: tuple[str, ...]
    resolved: tuple[str, ...]
    size: int


@dataclass
class TransferCluster:
    genome: str
    complex_id: str
    name: str
    genes: frozenset[str]
    genes_k12: frozenset[str]
    contig: str
    start: int
    end: int
    target_span: int
    target_coloc: bool
    k12_contig: str
    k12_start: int
    k12_end: int
    k12_span: int
    k12_coloc: bool
    n_target_genes: int
    frac_mapped: float
    identity_mean: float
    score: float
    evidence_class: str


@dataclass
class Signature:
    cluster_id: str
    genes: frozenset[str]
    confidence: str
    evidencia: str
    evidence: str
    curado_id: str
    contig: str
    start: int
    end: int
    span: int
    coloc: bool


def load_gene_coords(path: Path) -> dict[str, GeneCoord]:
    """Coordenadas por ``gene_id`` y por ``locus_tag`` (ambos alias)."""
    coords: dict[str, GeneCoord] = {}
    if not path.is_file():
        return coords
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            contig = (row.get("contig") or "").strip()
            start = _int(row.get("start"))
            end = _int(row.get("end"))
            if not contig or start <= 0 or end < start:
                continue
            coord = GeneCoord(contig, start, end)
            for key in ("gene_id", "locus_tag"):
                key_value = (row.get(key) or "").strip()
                if key_value and key_value not in coords:
                    coords[key_value] = coord
    return coords


def _interval_for_genes(genes: Iterable[str],
                        coords: dict[str, GeneCoord]) -> tuple[str, int, int] | None:
    """Intervalo del contig con más genes mapeados; ``None`` si no hay."""
    by_contig: dict[str, list[GeneCoord]] = defaultdict(list)
    for gene in genes:
        coord = coords.get(gene)
        if coord is not None:
            by_contig[coord.contig].append(coord)
    if not by_contig:
        return None
    contig = sorted(by_contig, key=lambda c: (-len(by_contig[c]), c))[0]
    values = by_contig[contig]
    return contig, min(c.start for c in values), max(c.end for c in values)


def load_model_clusters(path: Path, genome: str,
                        coords: dict[str, GeneCoord]) -> list[ModelCluster]:
    clusters: list[ModelCluster] = []
    if not path.is_file():
        return clusters
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes = tuple(parse_gene_ids(row.get("subunits")))
            if not genes:
                genes = tuple(parse_gene_ids(row.get("missing")))
            resolved = tuple(g for g in genes if g in coords)
            contig = (row.get("contig") or "").strip()
            start = _int(row.get("start"))
            end = _int(row.get("end"))
            if (not contig or start <= 0 or end < start) and resolved:
                coord = coords[resolved[0]]
                contig = contig or coord.contig
                start = start if start > 0 else coord.start
                end = end if end >= start else coord.end
            size = len(genes) if genes else len(resolved)
            clusters.append(ModelCluster(
                genome=genome,
                complex_id=(row.get("complex_id") or "").strip(),
                name=(row.get("name") or "").strip(),
                contig=contig,
                start=start,
                end=end,
                confidence=normalize_confidence(row.get("confidence")),
                evidence=(row.get("evidence") or "").strip(),
                status=(row.get("status") or "").strip(),
                genes=genes,
                resolved=resolved,
                size=size,
            ))
    return clusters


def load_transfers(path: Path, genome: str,
                   coords: dict[str, GeneCoord],
                   k12_coords: dict[str, GeneCoord] | None = None,
                   max_colocal_span: int = DEFAULT_MAX_COLOCAL_SPAN,
                   ) -> list[TransferCluster]:
    transfers: list[TransferCluster] = []
    if not path.is_file():
        return transfers
    k12_coords = k12_coords or {}
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            target = parse_gene_ids(row.get("genes_genoma"))
            genes = frozenset(g for g in target if g in coords)
            interval = _interval_for_genes(genes, coords)
            contig, start, end = interval if interval else ("", 0, 0)
            target_span = max(0, end - start)
            target_coloc = bool(
                contig and target_span > 0
                and target_span <= max_colocal_span)
            k12 = frozenset(parse_gene_ids(row.get("genes_k12")))
            k12_interval = _interval_for_genes(k12, k12_coords)
            k12_contig, k12_start, k12_end = (
                k12_interval if k12_interval else ("", 0, 0))
            k12_span = max(0, k12_end - k12_start)
            k12_coloc = bool(
                k12_contig and k12_span > 0
                and k12_span <= max_colocal_span)
            transfers.append(TransferCluster(
                genome=genome,
                complex_id=(row.get("complex_id") or "").strip(),
                name=(row.get("complex_name") or "").strip(),
                genes=genes,
                genes_k12=k12,
                contig=contig,
                start=start,
                end=end,
                target_span=target_span,
                target_coloc=target_coloc,
                k12_contig=k12_contig,
                k12_start=k12_start,
                k12_end=k12_end,
                k12_span=k12_span,
                k12_coloc=k12_coloc,
                n_target_genes=len(target),
                frac_mapped=_float(row.get("frac_mapped")),
                identity_mean=_float(row.get("identity_mean")),
                score=_float(row.get("score")),
                evidence_class=(row.get("evidence_class") or "").strip(),
            ))
    return transfers


def load_signatures(path: Path, coords: dict[str, GeneCoord],
                    max_colocal_span: int = DEFAULT_MAX_COLOCAL_SPAN,
                    ) -> list[Signature]:
    signatures: list[Signature] = []
    if not path.is_file():
        return signatures
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes = frozenset(parse_gene_ids(row.get("genes")))
            interval = _interval_for_genes(genes, coords)
            contig, start, end = interval if interval else ("", 0, 0)
            span = max(0, end - start)
            coloc = bool(contig and span > 0 and span <= max_colocal_span)
            signatures.append(Signature(
                cluster_id=(row.get("cluster_id") or "").strip(),
                genes=genes,
                confidence=normalize_confidence(row.get("confidence")),
                evidencia=(row.get("evidencia") or "").strip(),
                evidence=(row.get("evidence") or "").strip(),
                curado_id=(row.get("curado_id") or "").strip(),
                contig=contig,
                start=start,
                end=end,
                span=span,
                coloc=coloc,
            ))
    return signatures


def load_gold_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    if not path.is_file():
        return ids
    with _open_tsv(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("complex_id") or "").strip()
            if cid:
                ids.add(cid)
    return ids


# --------------------------------------------------------------------------- #
# Índice de firmas
# --------------------------------------------------------------------------- #
class SignatureIndex:
    """Índice por gen y por bins de intervalo para consultar 960 firmas."""

    def __init__(self, signatures: Sequence[Signature], bin_size: int = 10_000):
        self.signatures = signatures
        self.bin_size = max(1, int(bin_size))
        self.gene_to_ids: dict[str, set[int]] = defaultdict(set)
        self.bins: dict[tuple[str, int], list[int]] = defaultdict(list)
        for idx, sig in enumerate(signatures):
            for gene in sig.genes:
                self.gene_to_ids[gene].add(idx)
            if sig.coloc and sig.contig and sig.start > 0 and sig.end >= sig.start:
                first = sig.start // self.bin_size
                last = sig.end // self.bin_size
                for b in range(first, last + 1):
                    self.bins[(sig.contig, b)].append(idx)

    def gene_candidates(self, genes: Iterable[str]) -> set[int]:
        out: set[int] = set()
        for gene in genes:
            out.update(self.gene_to_ids.get(gene, ()))
        return out

    def interval_candidates(self, contig: str, start: int, end: int) -> set[int]:
        if not contig or start <= 0 or end < start:
            return set()
        out: set[int] = set()
        first = start // self.bin_size
        last = end // self.bin_size
        for b in range(first, last + 1):
            out.update(self.bins.get((contig, b), ()))
        return out


# --------------------------------------------------------------------------- #
# Evaluación de un cluster
# --------------------------------------------------------------------------- #
def _best_transfer(cluster: ModelCluster,
                   transfers: Sequence[TransferCluster],
                   max_colocal_span: int = DEFAULT_MAX_COLOCAL_SPAN,
                   ) -> tuple[
                       TransferCluster | None, float, float, int,
                       TransferCluster | None, float,
                       TransferCluster | None, float]:
    """Devuelve mejor transferencia por gen, por intervalo y por intervalo amplio.

    La mejor por gen maximiza ``jaccard + cobertura`` para que la tabla no
    premie solo a clusters pequeños. ``t_iv_frac`` solo considera intervalos
    **colocalizados** (span <= ``max_colocal_span``); ``t_iv_frac_any`` usa
    cualquier intervalo y se conserva únicamente como diagnóstico de lo
    permisivo que resulta cuando las subunidades están dispersas.
    """
    cluster_genes = set(cluster.resolved)
    best_gene: TransferCluster | None = None
    best_rank = (-1.0, -1.0, -1.0)
    best_cov = 0.0
    best_jac = 0.0
    best_inter = 0
    best_iv: TransferCluster | None = None
    best_iv_frac = 0.0
    best_iv_any: TransferCluster | None = None
    best_iv_any_frac = 0.0

    if not cluster_genes:
        return None, 0.0, 0.0, 0, None, 0.0, None, 0.0

    cluster_span = max(0, cluster.end - cluster.start)
    cluster_coloc = bool(
        cluster.contig and cluster.start > 0 and cluster_span > 0
        and cluster_span <= max_colocal_span)

    for transfer in transfers:
        inter = len(cluster_genes & transfer.genes)
        if inter:
            cov = inter / len(cluster_genes)
            jac = inter / len(cluster_genes | transfer.genes)
            rank = (jac + cov, jac, cov)
            if rank > best_rank:
                best_rank = rank
                best_gene = transfer
                best_cov = cov
                best_jac = jac
                best_inter = inter
        if (transfer.contig != cluster.contig or transfer.start <= 0
                or cluster.start <= 0):
            continue
        frac = interval_overlap_fraction(
            cluster.start, cluster.end, transfer.start, transfer.end)
        if frac > best_iv_any_frac:
            best_iv_any_frac = frac
            best_iv_any = transfer
        if cluster_coloc and transfer.target_coloc and frac > best_iv_frac:
            best_iv_frac = frac
            best_iv = transfer

    return (best_gene, best_cov, best_jac, best_inter,
            best_iv, best_iv_frac, best_iv_any, best_iv_any_frac)


def _signature_overlap(bridge: Sequence[TransferCluster],
                       index: SignatureIndex,
                       max_colocal_span: int = DEFAULT_MAX_COLOCAL_SPAN,
                       ) -> tuple[float, str, float, str]:
    """Mejor Jaccard por gen y mejor fracción de intervalo colocalizado.

    El Jaccard usa todas las filas-puente. El intervalo K-12 solo se evalúa
    si la fila-puente tiene sus genes K-12 colocalizados (span <= umbral) y
    la firma también está indexada como colocalizada.
    """
    best_gene = 0.0
    best_gene_id = ""
    best_iv = 0.0
    best_iv_id = ""

    for transfer in bridge:
        k12_genes = set(transfer.genes_k12)
        if not k12_genes:
            continue
        for idx in index.gene_candidates(k12_genes):
            score = jaccard(k12_genes, index.signatures[idx].genes)
            if score > best_gene:
                best_gene = score
                best_gene_id = index.signatures[idx].cluster_id

        if not transfer.k12_coloc:
            continue
        for idx in index.interval_candidates(
                transfer.k12_contig, transfer.k12_start, transfer.k12_end):
            sig = index.signatures[idx]
            frac = interval_overlap_fraction(
                transfer.k12_start, transfer.k12_end, sig.start, sig.end)
            if frac > best_iv:
                best_iv = frac
                best_iv_id = sig.cluster_id

    return best_gene, best_gene_id, best_iv, best_iv_id


def priority_score(confidence: str, size: int, t_cov: float, t_jac: float,
                   s_jac: float, t_iv: float) -> float:
    """Score heurístico 0-5,5. No es una probabilidad."""
    conf = CONF_WEIGHT.get(confidence, 0.0)
    transfer_gene = t_cov + t_jac
    signature_gene = 0.5 * s_jac
    transfer_interval = 0.5 * t_iv
    size_component = 0.5 * (min(max(size, 2), 12) / 12.0)
    return _round(conf + transfer_gene + signature_gene
                  + transfer_interval + size_component, 4)


def priority_tier(score: float) -> str:
    if score >= 4.0:
        return "P1"
    if score >= 3.0:
        return "P2"
    return "P3"


def evaluate_cluster(cluster: ModelCluster,
                     transfers: Sequence[TransferCluster],
                     index: SignatureIndex,
                     sig_curado_ids: set[str],
                     max_colocal_span: int = DEFAULT_MAX_COLOCAL_SPAN,
                     ) -> dict:
    genes = cluster.resolved
    (t_best, t_cov, t_jac, t_inter, t_iv_best, t_iv,
     t_iv_any_best, t_iv_any) = _best_transfer(
        cluster, transfers, max_colocal_span)

    bridge = [
        t for t in transfers
        if genes and (set(genes) & t.genes)
    ]
    s_jac, s_best, s_iv, s_iv_best = _signature_overlap(
        bridge, index, max_colocal_span)

    score = priority_score(cluster.confidence, cluster.size,
                           t_cov, t_jac, s_jac, t_iv)
    flags = {
        "has_transfer_gene": bool(t_best),
        "has_transfer_iv": t_iv > 0.0,
        "has_transfer_iv_any": t_iv_any > 0.0,
        "strong_transfer_gene": t_cov >= 0.5,
        "strong_transfer_iv": t_iv >= 0.5,
        "strong_transfer_iv_any": t_iv_any >= 0.5,
        "has_signature_gene": s_jac > 0.0,
        "has_signature_iv": s_iv > 0.0,
        "has_bridge": bool(bridge),
        "has_bridge_coloc": any(t.k12_coloc for t in bridge),
        "bridge_signature_backed": any(
            t.complex_id in sig_curado_ids for t in bridge),
        "sin_genes_resueltos": not genes,
    }
    metrics = {
        "t_gene_cov": _round(t_cov),
        "t_gene_jac": _round(t_jac),
        "t_iv_frac": _round(t_iv),
        "t_iv_frac_any": _round(t_iv_any),
        "s_gene_jac": _round(s_jac),
        "s_iv_frac": _round(s_iv),
        "priority_score": score,
        "n_resolved": len(genes),
    }
    return {
        "genome": cluster.genome,
        "complex_id": cluster.complex_id,
        "name": cluster.name,
        "contig": cluster.contig,
        "start": cluster.start,
        "end": cluster.end,
        "confidence": cluster.confidence,
        "status": cluster.status,
        "size": cluster.size,
        "n_resolved": len(genes),
        "evidence": cluster.evidence,
        "evidence_tokens": evidence_tokens(cluster.evidence),
        "t_best_id": t_best.complex_id if t_best else (t_iv_best.complex_id if t_iv_best else ""),
        "t_best_name": t_best.name if t_best else (t_iv_best.name if t_iv_best else ""),
        "t_best_genes": t_inter,
        "t_iv_best_id": t_iv_best.complex_id if t_iv_best else "",
        "t_iv_any_best_id": t_iv_any_best.complex_id if t_iv_any_best else "",
        "s_best_id": s_best,
        "s_iv_best_id": s_iv_best,
        "bridge_rows": len(bridge),
        "bridge_transfer_ids": sorted({t.complex_id for t in bridge}),
        "metrics": metrics,
        "flags": flags,
        "priority_tier": priority_tier(score),
    }


# --------------------------------------------------------------------------- #
# Acumuladores
# --------------------------------------------------------------------------- #
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


class Stats:
    """Acumula un grupo (global, confianza, tamaño o genoma)."""

    def __init__(self, keep_series: bool = True, evidence_limit: int = 15):
        self.keep_series = keep_series
        self.evidence_limit = evidence_limit
        self.n = 0
        self.sizes: list[int] = []
        self.confidence: Counter[str] = Counter()
        self.evidence: Counter[str] = Counter()
        self.flags: Counter[str] = Counter()
        self.sums: dict[str, float] = defaultdict(float)
        self.maxima: dict[str, float] = defaultdict(float)
        self.series: dict[str, list[float]] = (
            defaultdict(list) if keep_series else {})

    def add(self, record: dict) -> None:
        self.n += 1
        self.sizes.append(int(record["size"]))
        self.confidence[record["confidence"]] += 1
        for token in record["evidence_tokens"]:
            self.evidence[token] += 1
        for key in FLAG_KEYS:
            if record["flags"].get(key):
                self.flags[key] += 1
        for key, value in record["metrics"].items():
            self.sums[key] += float(value)
            if float(value) > self.maxima[key]:
                self.maxima[key] = float(value)
            if self.keep_series:
                self.series[key].append(float(value))

    def _pct(self, count: int) -> float:
        return _round(100.0 * count / self.n, 2) if self.n else 0.0

    def to_dict(self) -> dict:
        metrics: dict[str, dict] = {}
        for key in METRIC_KEYS:
            if self.keep_series:
                metrics[key] = _series_stats(self.series.get(key, []))
            else:
                metrics[key] = {
                    "n": self.n,
                    "mean": _round(self.sums[key] / self.n) if self.n else 0.0,
                    "max": _round(self.maxima[key]),
                }
        overlap = {key: self.flags.get(key, 0) for key in FLAG_KEYS}
        overlap_pct = {key: self._pct(overlap[key]) for key in FLAG_KEYS}
        return {
            "clusters": self.n,
            "confianza": dict(self.confidence),
            "tamano": _series_stats(self.sizes),
            "tamano_bins": {
                label: sum(1 for size in self.sizes
                           if size_bin_label(size) == label)
                for _, _, label in SIZE_BINS
            },
            "evidencia_top": [
                {"token": token, "clusters": count}
                for token, count in self.evidence.most_common(self.evidence_limit)
            ],
            "solapamiento": overlap,
            "solapamiento_pct": overlap_pct,
            "metricas": metrics,
        }


# --------------------------------------------------------------------------- #
# Análisis principal
# --------------------------------------------------------------------------- #
def _model_files(complexes_dir: Path) -> list[Path]:
    pattern = str(complexes_dir / f"{MODEL_PREFIX}*{MODEL_SUFFIX}")
    return sorted(Path(p) for p in glob.glob(pattern))


def genome_from_model_path(path: Path) -> str:
    name = path.name
    if name.startswith(MODEL_PREFIX) and name.endswith(MODEL_SUFFIX):
        return name[len(MODEL_PREFIX):-len(MODEL_SUFFIX)]
    return path.stem


def analyze(complexes_dir: Path, regulation_dir: Path,
            signatures_path: Path, gold_path: Path, k12_genes_path: Path,
            top_n: int = 200, limit_genomes: int = 0,
            bin_size: int = 10_000,
            max_colocal_span: int = DEFAULT_MAX_COLOCAL_SPAN) -> dict:
    complexes_dir = Path(complexes_dir)
    regulation_dir = Path(regulation_dir)
    signatures_path = Path(signatures_path)
    gold_path = Path(gold_path)
    k12_genes_path = Path(k12_genes_path)

    k12_coords = load_gene_coords(k12_genes_path)
    signatures = load_signatures(
        signatures_path, k12_coords, max_colocal_span)
    sig_index = SignatureIndex(signatures, bin_size)
    gold_ids = load_gold_ids(gold_path)
    sig_curado_ids = {s.curado_id for s in signatures if s.curado_id}

    model_files = _model_files(complexes_dir)
    if limit_genomes and limit_genomes > 0:
        model_files = model_files[:limit_genomes]

    global_stats = Stats(keep_series=True, evidence_limit=30)
    by_conf: dict[str, Stats] = {
        conf: Stats(keep_series=True, evidence_limit=15)
        for conf in CONF_ORDER
    }
    by_size: dict[str, Stats] = {
        label: Stats(keep_series=True, evidence_limit=15)
        for _, _, label in SIZE_BINS
    }
    by_genome: dict[str, dict] = {}
    evidence_by_conf: dict[str, Counter[str]] = {
        conf: Counter() for conf in CONF_ORDER}
    tier_counts: Counter[str] = Counter()
    tier_conf: dict[tuple[str, str], int] = defaultdict(int)
    anchor_conf: Counter[str] = Counter()

    status_counts: Counter[str] = Counter()
    transfer_rows = 0
    transfer_rows_coloc = 0
    transfer_rows_k12_coloc = 0
    transfer_spans: list[int] = []
    transfer_ids: set[str] = set()
    matched_transfer_ids: set[str] = set()
    matched_transfer_iv_ids: set[str] = set()
    matched_signature_ids: set[str] = set()
    matched_signature_iv_ids: set[str] = set()
    genomes_sin_genes: list[str] = []
    genomes_sin_transfer: list[str] = []
    top_heap: list[tuple[float, int, float, int, dict]] = []
    seq = 0

    for model_path in model_files:
        genome = genome_from_model_path(model_path)
        genes_path = regulation_dir / genome / "genes.tsv"
        if not genes_path.is_file():
            genomes_sin_genes.append(genome)
        coords = load_gene_coords(genes_path)
        clusters = load_model_clusters(model_path, genome, coords)

        transfer_path = complexes_dir / f"transfer_{genome}.tsv"
        if not transfer_path.is_file():
            genomes_sin_transfer.append(genome)
        transfers = load_transfers(
            transfer_path, genome, coords, k12_coords, max_colocal_span)
        transfer_rows += len(transfers)
        transfer_rows_coloc += sum(1 for t in transfers if t.target_coloc)
        transfer_rows_k12_coloc += sum(1 for t in transfers if t.k12_coloc)
        transfer_spans.extend(t.target_span for t in transfers
                              if t.target_span > 0)
        transfer_ids.update(t.complex_id for t in transfers)
        genome_stats = Stats(keep_series=False, evidence_limit=10)
        for cluster in clusters:
            record = evaluate_cluster(
                cluster,
                transfers,
                sig_index,
                sig_curado_ids,
                max_colocal_span,
            )
            global_stats.add(record)
            genome_stats.add(record)
            conf = record["confidence"]
            by_conf.setdefault(
                conf, Stats(keep_series=True, evidence_limit=15)
            ).add(record)
            by_size.setdefault(
                size_bin_label(record["size"]),
                Stats(keep_series=True, evidence_limit=15),
            ).add(record)
            evidence_by_conf.setdefault(conf, Counter()).update(
                record["evidence_tokens"])
            status_counts[record["status"]] += 1

            tier = record["priority_tier"]
            tier_counts[tier] += 1
            tier_conf[(conf, tier)] += 1
            if record["flags"]["has_transfer_gene"]:
                anchor_conf[conf] += 1

            if record["t_best_id"]:
                matched_transfer_ids.add(record["t_best_id"])
            if record["t_iv_best_id"]:
                matched_transfer_iv_ids.add(record["t_iv_best_id"])
            if record["s_best_id"]:
                matched_signature_ids.add(record["s_best_id"])
            if record["s_iv_best_id"]:
                matched_signature_iv_ids.add(record["s_iv_best_id"])

            seq += 1
            heapq.heappush(top_heap, (
                record["metrics"]["priority_score"],
                record["size"],
                record["metrics"]["t_gene_cov"],
                seq,
                record,
            ))
            if len(top_heap) > max(1, top_n):
                heapq.heappop(top_heap)

        by_genome[genome] = genome_stats.to_dict()

    ranked = sorted(
        top_heap,
        key=lambda item: (item[0], item[1], item[2]),
        reverse=True,
    )
    top: list[dict] = []
    for rank, (_, _, _, _, record) in enumerate(ranked, 1):
        record["rank"] = rank
        top.append(record)

    total_clusters = global_stats.n
    total_sizes = global_stats.sizes
    transfer_complexes = len(transfer_ids)
    signature_curado_rows = sum(1 for s in signatures if s.curado_id)
    signature_curado_unique = len(sig_curado_ids)
    signature_with_match = len(matched_signature_ids)
    signature_with_iv_match = len(matched_signature_iv_ids)
    signature_coloc = sum(1 for s in signatures if s.coloc)

    gold_block = {
        "complejos_gold": len(gold_ids),
        "complejos_con_transferencia": len(gold_ids & transfer_ids),
        "complejos_gold_sin_transferencia": sorted(gold_ids - transfer_ids),
        "complejos_transferidos": transfer_complexes,
        "complejos_transferidos_en_gold": len(transfer_ids & gold_ids),
        "firmas_con_curado_id": signature_curado_rows,
        "firmas_curado_id_unicos": signature_curado_unique,
        "firmas_curado_id_en_gold": len(sig_curado_ids & gold_ids),
        "firmas_sin_curado_id": len(signatures) - signature_curado_rows,
    }

    overlap_transfer = {
        "filas_transferencia": transfer_rows,
        "filas_con_intervalo_colocalizado": transfer_rows_coloc,
        "filas_sin_intervalo_colocalizado": transfer_rows - transfer_rows_coloc,
        "filas_con_k12_colocalizado": transfer_rows_k12_coloc,
        "mediana_span_intervalo": _series_stats(transfer_spans)["median"],
        "complejos_k12": transfer_complexes,
        "complejos_con_al_menos_un_cluster_por_gen": len(matched_transfer_ids),
        "complejos_con_al_menos_un_cluster_por_intervalo": len(matched_transfer_iv_ids),
    }
    overlap_signature = {
        "firmas": len(signatures),
        "firmas_con_intervalo_colocalizado": signature_coloc,
        "firmas_sin_intervalo_colocalizado": len(signatures) - signature_coloc,
        "mediana_span_intervalo": _series_stats(
            [s.span for s in signatures if s.span > 0])["median"],
        "firmas_con_match_por_gen": signature_with_match,
        "firmas_con_match_por_intervalo": signature_with_iv_match,
        "clusters_con_puente_k12": global_stats.flags.get("has_bridge", 0),
        "clusters_con_puente_k12_colocalizado": global_stats.flags.get(
            "has_bridge_coloc", 0),
        "clusters_sin_puente_k12": total_clusters - global_stats.flags.get("has_bridge", 0),
        "clusters_con_match_firma_por_gen": global_stats.flags.get("has_signature_gene", 0),
        "clusters_con_match_firma_por_intervalo": global_stats.flags.get("has_signature_iv", 0),
    }

    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    report: dict = {
        "meta": {
            "herramienta": "tools/calibrate_model_clusters.py",
            "generado_utc": generated,
            "sin_red": True,
            "no_modifica_datos": True,
            "advertencia": (
                "Los clusters son hipótesis. Los solapamientos son descriptivos; "
                "no son precisión, recall ni validación experimental."
            ),
        },
        "inputs": {
            "directorio_complejos": str(complexes_dir),
            "directorio_regulacion": str(regulation_dir),
            "firmas_k12": str(signatures_path),
            "gold_k12": str(gold_path),
            "k12_genes": str(k12_genes_path),
            "genomas_modelo": len(model_files),
            "filas_modelo": total_clusters,
            "genomas_con_genes_tsv": len(model_files) - len(genomes_sin_genes),
            "genomas_sin_genes_tsv": genomes_sin_genes,
            "genomas_con_transferencia": len(model_files) - len(genomes_sin_transfer),
            "genomas_sin_transferencia": genomes_sin_transfer,
            "filas_transferencia": transfer_rows,
            "complejos_transferencia_unicos": transfer_complexes,
            "firmas_k12_total": len(signatures),
            "complejos_gold_k12": len(gold_ids),
            "genes_k12_con_coordenadas": len(k12_coords),
            "max_span_colocalizado": max_colocal_span,
            "filas_transferencia_intervalo_colocalizado": transfer_rows_coloc,
            "firmas_intervalo_colocalizado": signature_coloc,
        },
        "definiciones": {
            "cluster": "Una fila de model_<genoma>.complexes.tsv (complex_candidate).",
            "tamano_cluster": "Número de subunidades listadas en `subunits`; no es estequiometría.",
            "confianza": "Etiqueta del modelo (alta/media/baja); no es probabilidad calibrada.",
            "solapamiento_gen_transferencia": (
                "Intersección de los genes resueltos del cluster con `genes_genoma` "
                "de una fila de transferencia del mismo genoma."
            ),
            "jaccard_gen": "|intersección| / |unión| entre genes del cluster y de la transferencia.",
            "cobertura_gen_modelo": "|intersección| / |genes del cluster|.",
            "intervalo_colocalizado": (
                f"Un intervalo se considera colocalizado si sus genes están en un "
                f"mismo contig y abarcan <= {max_colocal_span} pb. Los intervalos "
                "dispersos se excluyen de la métrica de intervalo porque un span "
                "grande hace que casi cualquier cluster solape."
            ),
            "solapamiento_intervalo_transferencia": (
                "Fracción del intervalo [start,end) del cluster cubierta por el "
                "intervalo de una transferencia **colocalizada** en el mismo contig."
            ),
            "solapamiento_intervalo_transferencia_amplio": (
                "Igual que el anterior pero sin exigir colocalización; se reporta "
                "solo como diagnóstico de permisividad."
            ),
            "firma_k12": (
                "Cluster de clusters_ecoli_k12_evidencia.tsv (960 firmas K-12); "
                "sus genes son locus tags b#### de K-12."
            ),
            "puente_k12": (
                "Filas de transferencia que comparten al menos un gen del genoma con "
                "el cluster; sus `genes_k12` se comparan contra las firmas. Es un "
                "puente RBH/complejo, no ortología demostrada."
            ),
            "solapamiento_gen_firma": (
                "Máximo Jaccard entre los genes K-12 de una fila-puente y los genes "
                "de una firma K-12."
            ),
            "solapamiento_intervalo_firma": (
                "Fracción del intervalo K-12 de una fila-puente **colocalizada** "
                "cubierta por el intervalo K-12 de una firma **colocalizada**."
            ),
            "prioridad": (
                "Heurística de revisión, no probabilidad. Ver `priorizacion.pesos`."
            ),
        },
        "totales": {
            "genomas": len(model_files),
            "clusters": total_clusters,
            "subunidades_listadas": sum(total_sizes),
            "genes_resueltos": int(global_stats.sums.get("n_resolved", 0)),
            "clusters_sin_genes_resueltos": global_stats.flags.get(
                "sin_genes_resueltos", 0),
            "estados": dict(status_counts),
            "tamano_cluster": _series_stats(total_sizes),
            "confianza": dict(global_stats.confidence),
        },
        "agregado_global": global_stats.to_dict(),
        "por_confianza": {
            conf: stats.to_dict() for conf, stats in by_conf.items()
            if stats.n
        },
        "por_tamano": {
            label: stats.to_dict() for label, stats in by_size.items()
            if stats.n
        },
        "por_genoma": by_genome,
        "evidencia": {
            "top_global": [
                {"token": token, "clusters": count}
                for token, count in global_stats.evidence.most_common(30)
            ],
            "por_confianza": {
                conf: [
                    {"token": token, "clusters": count}
                    for token, count in counter.most_common(15)
                ]
                for conf, counter in evidence_by_conf.items()
            },
        },
        "solapamiento": {
            "transferencias": {
                **overlap_transfer,
                "global": global_stats.to_dict()["solapamiento"],
                "global_pct": global_stats.to_dict()["solapamiento_pct"],
            },
            "firmas": {
                **overlap_signature,
                "global": global_stats.to_dict()["solapamiento"],
                "global_pct": global_stats.to_dict()["solapamiento_pct"],
            },
            "gold": gold_block,
        },
        "priorizacion": {
            "formula": (
                "score = peso_confianza + (cobertura_gen + jaccard_gen) "
                "+ 0,5*jaccard_firma + 0,5*fraccion_intervalo_colocalizado "
                "+ 0,5*min(tamano,12)/12"
            ),
            "pesos": {
                "confianza": dict(CONF_WEIGHT),
                "solapamiento_gen_transferencia": 1.0,
                "solapamiento_gen_firma": 0.5,
                "solapamiento_intervalo_transferencia": 0.5,
                "tamano": 0.5,
                "techo_tamano": 12,
            },
            "tiers": {
                "P1": "score >= 4,0",
                "P2": "3,0 <= score < 4,0",
                "P3": "score < 3,0",
            },
            "conteo_tiers": dict(tier_counts),
            "conteo_tiers_por_confianza": {
                conf: {tier: tier_conf[(conf, tier)] for tier in ("P1", "P2", "P3")}
                for conf in by_conf
            },
            "clusters_con_anclaje_transferencia_por_confianza": dict(anchor_conf),
            "top_n": len(top),
            "top": top,
        },
        "limitaciones": _limitations(
            total_clusters=total_clusters,
            transfer_rows=transfer_rows,
            transfer_complexes=transfer_complexes,
            signatures=len(signatures),
            signature_curado_rows=signature_curado_rows,
            gold_complexes=len(gold_ids),
            overlap_pct=global_stats.to_dict()["solapamiento_pct"],
            signature_with_match=signature_with_match,
            signature_with_iv_match=signature_with_iv_match,
            transfer_rows_coloc=transfer_rows_coloc,
            signature_coloc=signature_coloc,
            max_colocal_span=max_colocal_span,
            transfer_span_median=overlap_transfer["mediana_span_intervalo"],
            signature_span_median=overlap_signature["mediana_span_intervalo"],
        ),
    }
    return report


def _limitations(total_clusters: int, transfer_rows: int,
                 transfer_complexes: int, signatures: int,
                 signature_curado_rows: int, gold_complexes: int,
                 overlap_pct: dict, signature_with_match: int,
                 signature_with_iv_match: int, transfer_rows_coloc: int,
                 signature_coloc: int, max_colocal_span: int,
                 transfer_span_median: float,
                 signature_span_median: float) -> list[str]:
    broad_iv = _es_dec(overlap_pct.get("has_transfer_iv_any", 0), 1)
    coloc_iv = _es_dec(overlap_pct.get("has_transfer_iv", 0), 1)
    gene_pct = _es_dec(overlap_pct.get("has_transfer_gene", 0), 1)
    return [
        (
            f"El modelo es un prior entrenado/derivado de K-12: los "
            f"{_es_int(total_clusters)} clusters de los 400 genomas heredan esa "
            "distribución. No hay gold independiente para las otras especies; "
            f"gold_k12.tsv cubre solo {_es_int(gold_complexes)} complejos de K-12."
        ),
        (
            "Los clusters son candidatos (`complex_candidate`): no demuestran "
            "interacción física, estequiometría ni función. La confianza "
            "alta/media/baja no es una probabilidad calibrada."
        ),
        (
            f"Las {_es_int(transfer_rows)} filas de transferencia cubren "
            f"{_es_int(transfer_complexes)} complejos K-12 y usan RBH; no hay "
            "cobertura de alineamiento en los TSV actuales. La evidencia "
            "transferida es `documentado_inferido`/ECO:0005544, no experimental."
        ),
        (
            f"Las {_es_int(signatures)} firmas K-12 solo cubren "
            f"{_es_int(signature_curado_rows)} filas con `curado_id`; el resto "
            f"no tiene ancla curada. Solo {_es_int(signature_with_match)} firmas "
            f"reciben al menos un cluster por gen y "
            f"{_es_int(signature_with_iv_match)} por intervalo colocalizado. "
            "La comparación con firmas es indirecta (puente de transferencia) y "
            "puede inflarse o perderse por parálogos y mapeos muchos-a-muchos."
        ),
        (
            "El intervalo amplio es engañoso: los genes de una transferencia "
            f"tienen span mediano de {_es_int(transfer_span_median)} pb y los de "
            f"una firma {_es_int(signature_span_median)} pb; por eso "
            f"{broad_iv}% de los clusters tocan algún intervalo amplio frente a "
            f"{coloc_iv}% con intervalo colocalizado (<= "
            f"{_es_int(max_colocal_span)} pb) y {gene_pct}% con al menos un gen "
            f"compartido. Solo {_es_int(transfer_rows_coloc)}/"
            f"{_es_int(transfer_rows)} transferencias y "
            f"{_es_int(signature_coloc)}/{_es_int(signatures)} firmas tienen "
            "intervalo colocalizado. Un solapamiento de intervalo por sí solo "
            "no es evidencia."
        ),
        (
            "La priorización usa pesos elegidos a mano (confianza, gen, "
            "intervalo colocalizado y tamaño). Cambiar los pesos cambia el "
            "orden. Los tiers P1/P2/P3 son cortes operativos, no clases "
            "biológicas."
        ),
        (
            "Los intervalos se derivan del mínimo/máximo de genes mapeados; en "
            "transferencias multi-contig se elige el contig con más genes y se "
            "pierde el resto. Los alias `gene_id`/`locus_tag` pueden colapsar "
            "duplicados."
        ),
        (
            "El informe es descriptivo y sin red. No se editó el pipeline, la "
            "GUI ni los datos; no se hizo commit."
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


def _pct(count: int, total: int, ndigits: int = 1) -> str:
    if not total:
        return "—"
    return f"{100.0 * count / total:.{ndigits}f}".replace(".", ",") + " %"


def _md_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    lines = ["| " + " | ".join(str(h) for h in headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(lines)


def _short(text: object, limit: int = 60) -> str:
    value = str(text or "").strip().replace("|", "/")
    return value if len(value) <= limit else value[: limit - 1] + "…"


def build_markdown(report: dict, top_rows: int = 30) -> str:
    inputs = report["inputs"]
    totals = report["totales"]
    por_conf = report["por_confianza"]
    por_tam = report["por_tamano"]
    global_stats = report["agregado_global"]
    transfer = report["solapamiento"]["transferencias"]
    firmas = report["solapamiento"]["firmas"]
    gold = report["solapamiento"]["gold"]
    prior = report["priorizacion"]
    top = prior["top"]

    lines: list[str] = []
    lines.append("# Calibración de clusters del modelo K-12 en los 400 genomas")
    lines.append("")
    lines.append(
        "> **Hipótesis, no validación.** Los clusters son candidatos del modelo "
        "K-12. Este informe mide solapamiento descriptivo y prioriza revisiones; "
        "no estima precisión, recall ni interacción física. Sin red.")
    lines.append("")
    lines.append(f"Generado: `{report['meta']['generado_utc']}` · "
                 f"herramienta: `{report['meta']['herramienta']}`")
    lines.append("")

    lines.append("## 1. Insumos y alcance")
    lines.append("")
    lines.append(_md_table(
        ["Insumo", "Valor"],
        [
            ["Genomas con `model_<genoma>.complexes.tsv`", _es_int(inputs["genomas_modelo"])],
            ["Filas de clusters", _es_int(inputs["filas_modelo"])],
            ["Genomas con `genes.tsv`", _es_int(inputs["genomas_con_genes_tsv"])],
            ["Genomas con `transfer_<genoma>.tsv`", _es_int(inputs["genomas_con_transferencia"])],
            ["Filas de transferencia", _es_int(inputs["filas_transferencia"])],
            ["Complejos K-12 transferidos (únicos)", _es_int(inputs["complejos_transferencia_unicos"])],
            ["Firmas K-12", _es_int(inputs["firmas_k12_total"])],
            ["Complejos gold K-12", _es_int(inputs["complejos_gold_k12"])],
            ["Genes K-12 con coordenadas", _es_int(inputs["genes_k12_con_coordenadas"])],
            ["Umbral de intervalo colocalizado", _es_int(inputs["max_span_colocalizado"]) + " pb"],
            ["Transferencias con intervalo colocalizado", _es_int(inputs["filas_transferencia_intervalo_colocalizado"])],
            ["Firmas con intervalo colocalizado", _es_int(inputs["firmas_intervalo_colocalizado"])],
        ]))
    lines.append("")
    lines.append(
        f"Los {_es_int(totals['clusters'])} clusters contienen "
        f"{_es_int(totals['subunidades_listadas'])} subunidades listadas y "
        f"{_es_int(totals['genes_resueltos'])} genes resueltos con coordenadas. "
        f"Hay {_es_int(totals['clusters_sin_genes_resueltos'])} clusters sin "
        "ningún gen resuelto. El estado es siempre `complex_candidate`.")
    lines.append("")
    t_global = transfer["global"]
    t_pct = transfer["global_pct"]
    alta = por_conf.get("alta", {})
    baja = por_conf.get("baja", {})
    top_first = prior["top"][0] if prior["top"] else None
    lines.append("**Hallazgos principales**")
    lines.append("")
    lines.append(
        f"- Confianza: {_es_dec(100 * por_conf.get('alta', {}).get('clusters', 0) / totals['clusters'], 1)}% "
        f"alta, {_es_dec(100 * por_conf.get('media', {}).get('clusters', 0) / totals['clusters'], 1)}% "
        f"media y {_es_dec(100 * por_conf.get('baja', {}).get('clusters', 0) / totals['clusters'], 1)}% "
        f"baja. La confianza alta es mayoritariamente de tamaño 2 "
        f"({_es_dec(100 * alta.get('tamano_bins', {}).get('2', 0) / max(1, alta.get('clusters', 0)), 1)}% "
        "de los clusters alta)."
    )
    lines.append(
        f"- Transferencias: {_es_int(t_global['has_transfer_gene'])} clusters "
        f"({_es_dec(t_pct['has_transfer_gene'], 2)}%) comparten al menos un gen "
        f"y {_es_int(t_global['strong_transfer_gene'])} "
        f"({_es_dec(t_pct['strong_transfer_gene'], 2)}%) tienen cobertura ≥ 0,5."
    )
    lines.append(
        f"- Intervalo: {_es_dec(t_pct['has_transfer_iv_any'], 2)}% solapa algún "
        f"intervalo amplio, pero solo {_es_dec(t_pct['has_transfer_iv'], 2)}% "
        f"uno colocalizado (≤ {_es_int(inputs['max_span_colocalizado'])} pb). "
        "El intervalo amplio no es evidencia."
    )
    lines.append(
        f"- Firmas: {_es_int(firmas['clusters_con_match_firma_por_gen'])} "
        f"clusters ({_pct(firmas['clusters_con_match_firma_por_gen'], totals['clusters'], 2)}) "
        f"tienen match por gen; se alcanzan "
        f"{_es_int(firmas['firmas_con_match_por_gen'])}/{_es_int(firmas['firmas'])} firmas."
    )
    if top_first:
        lines.append(
            f"- Priorización: {_es_int(prior['conteo_tiers'].get('P1', 0))} P1 y "
            f"{_es_int(prior['conteo_tiers'].get('P2', 0))} P2; el primer cluster "
            f"es `{top_first['complex_id']}` ({top_first['genome']}, "
            f"score {_es_dec(top_first['metrics']['priority_score'], 2)})."
        )
    lines.append("")

    lines.append("## 2. Agregación")
    lines.append("")
    lines.append("### 2.1 Por confianza")
    lines.append("")
    rows = []
    for conf in CONF_ORDER:
        data = por_conf.get(conf)
        if not data:
            continue
        size = data["tamano"]
        rows.append([
            conf, _es_int(data["clusters"]),
            _pct(data["clusters"], totals["clusters"], 1),
            _es_dec(size["mean"], 2), _es_dec(size["median"], 1),
            _es_int(data["tamano_bins"].get("2", 0)),
            _pct(data["tamano_bins"].get("2", 0), data["clusters"], 1),
            _es_dec(data["metricas"]["n_resolved"]["mean"], 2),
        ])
    lines.append(_md_table(
        ["Confianza", "Clusters", "% del total", "Tamaño medio", "Tamaño mediana",
         "Clusters de 2", "% de 2", "Genes resueltos (media)"],
        rows))
    lines.append("")
    lines.append(
        "La confianza alta está dominada por clusters de dos subunidades: "
        "el tamaño no crece con la confianza. Esto es una propiedad del modelo "
        "y debe considerarse al leer cualquier solapamiento. La confianza alta "
        f"tampoco implica más anclaje: {_es_dec(alta.get('solapamiento_pct', {}).get('has_transfer_gene', 0), 1)}% "
        "de los clusters alta comparte un gen con transferencia frente a "
        f"{_es_dec(baja.get('solapamiento_pct', {}).get('has_transfer_gene', 0), 1)}% "
        "de los baja.")
    lines.append("")

    lines.append("### 2.2 Por tamaño del cluster")
    lines.append("")
    rows = []
    for _, _, label in SIZE_BINS:
        data = por_tam.get(label)
        if not data:
            continue
        size = data["tamano"]
        rows.append([
            label.replace("+", " o más"),
            _es_int(data["clusters"]),
            _pct(data["clusters"], totals["clusters"], 1),
            _es_dec(size["mean"], 2),
            _es_dec(size["median"], 1),
            _es_dec(data["metricas"]["t_gene_cov"]["mean"], 3),
            _es_dec(data["metricas"]["t_iv_frac"]["mean"], 3),
        ])
    lines.append(_md_table(
        ["Subunidades", "Clusters", "% del total", "Media", "Mediana",
         "Cobertura gen transf. (media)", "Intervalo coloc. (media)"],
        rows))
    lines.append("")
    lines.append(
        "El solapamiento genético y el intervalo colocalizado crecen con el "
        "tamaño del cluster (cobertura media "
        f"{_es_dec(por_tam.get('2', {}).get('metricas', {}).get('t_gene_cov', {}).get('mean', 0), 3)} "
        "en tamaño 2 frente a "
        f"{_es_dec(por_tam.get('10+', {}).get('metricas', {}).get('t_gene_cov', {}).get('mean', 0), 3)} "
        "en 10+). Es consistente con que un cluster más grande tiene más "
        "oportunidades de contener un gen transferido; no implica mejor "
        "calidad.")
    lines.append("")

    lines.append("### 2.3 Por genoma")
    lines.append("")
    genome_sizes = [g["tamano"]["mean"] for g in report["por_genoma"].values()]
    genome_items = list(report["por_genoma"].items())
    genome_items.sort(key=lambda item: item[1]["clusters"], reverse=True)
    counts = [data["clusters"] for _, data in genome_items]
    lines.append(
        f"Hay {_es_int(len(report['por_genoma']))} genomas en el JSON. "
        f"Clusters por genoma: mínimo {_es_int(min(counts))}, "
        f"mediana {_es_dec(statistics.median(counts), 1)}, "
        f"máximo {_es_int(max(counts))}. "
        f"Tamaño medio por genoma: {_es_dec(min(genome_sizes), 2)}–"
        f"{_es_dec(max(genome_sizes), 2)} subunidades.")
    lines.append("")
    lines.append("Genomas con más clusters:")
    lines.append("")
    lines.append(_md_table(
        ["Genoma", "Clusters", "Confianza alta", "Tamaño medio", "Con gen transf."],
        [[name, _es_int(data["clusters"]),
          _es_int(data["confianza"].get("alta", 0)),
          _es_dec(data["tamano"]["mean"], 2),
          _pct(data["solapamiento"]["has_transfer_gene"], data["clusters"], 1)]
         for name, data in genome_items[:8]]))
    lines.append("")
    lines.append(
        "La tabla completa por genoma (confianza, tamaño, evidencia y "
        "solapamiento) está en el JSON, campo `por_genoma`.")
    lines.append("")

    lines.append("### 2.4 Evidencia")
    lines.append("")
    ev_rows = []
    for item in report["evidencia"]["top_global"][:12]:
        ev_rows.append([_short(item["token"], 55),
                        _es_int(item["clusters"]),
                        _pct(item["clusters"], totals["clusters"], 1)])
    lines.append(_md_table(
        ["Token de evidencia", "Clusters", "% del total"], ev_rows))
    lines.append("")
    lines.append(
        "La evidencia del modelo contiene el prior K-12 y señales genómicas "
        "(`operón`, `adyacente`, scores). No contiene STRING ni sintenia de "
        "K-12 aplicada a otros genomas.")
    lines.append("")

    lines.append("## 3. Solapamiento")
    lines.append("")
    lines.append(
        "**Gen**: Jaccard e intersección de genes. **Intervalo**: fracción del "
        "intervalo del cluster cubierta por el intervalo de referencia "
        f"**colocalizado** (mismo contig y span <= "
        f"{_es_int(inputs['max_span_colocalizado'])} pb). Son planos "
        "independientes: un intervalo puede solaparse sin compartir genes. No "
        "se afirma precisión.")
    lines.append("")

    lines.append("### 3.1 Con transferencias K-12 (mismo genoma)")
    lines.append("")
    t_global = transfer["global"]
    t_pct = transfer["global_pct"]
    lines.append(_md_table(
        ["Métrica", "Clusters", "% del total"],
        [
            ["Con al menos un gen compartido", _es_int(t_global["has_transfer_gene"]),
             _es_dec(t_pct["has_transfer_gene"], 1) + " %"],
            ["Cobertura de gen ≥ 0,5", _es_int(t_global["strong_transfer_gene"]),
             _es_dec(t_pct["strong_transfer_gene"], 1) + " %"],
            ["Con intervalo amplio (diagnóstico)", _es_int(t_global["has_transfer_iv_any"]),
             _es_dec(t_pct["has_transfer_iv_any"], 1) + " %"],
            ["Con intervalo colocalizado", _es_int(t_global["has_transfer_iv"]),
             _es_dec(t_pct["has_transfer_iv"], 1) + " %"],
            ["Intervalo colocalizado ≥ 0,5", _es_int(t_global["strong_transfer_iv"]),
             _es_dec(t_pct["strong_transfer_iv"], 1) + " %"],
            ["Con puente K-12 (≥1 fila)", _es_int(t_global["has_bridge"]),
             _es_dec(t_pct["has_bridge"], 1) + " %"],
        ]))
    lines.append("")
    lines.append(_md_table(
        ["Confianza", "Clusters", "Con gen", "%", "Int. amplio", "%",
         "Int. coloc.", "%", "Cob. gen mediana", "Jaccard mediano"],
        [[conf,
          _es_int(por_conf[conf]["clusters"]),
          _es_int(por_conf[conf]["solapamiento"]["has_transfer_gene"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_transfer_gene"], 1) + " %",
          _es_int(por_conf[conf]["solapamiento"]["has_transfer_iv_any"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_transfer_iv_any"], 1) + " %",
          _es_int(por_conf[conf]["solapamiento"]["has_transfer_iv"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_transfer_iv"], 1) + " %",
          _es_dec(por_conf[conf]["metricas"]["t_gene_cov"]["median"], 3),
          _es_dec(por_conf[conf]["metricas"]["t_gene_jac"]["median"], 3)]
         for conf in CONF_ORDER if conf in por_conf]))
    lines.append("")
    lines.append(
        f"Las transferencias cubren {_es_int(transfer['filas_transferencia'])} "
        f"filas y {_es_int(transfer['complejos_k12'])} complejos K-12. "
        f"{_es_int(transfer['filas_con_intervalo_colocalizado'])} filas tienen "
        f"intervalo colocalizado (span <= {_es_int(inputs['max_span_colocalizado'])} pb); "
        f"el span mediano de todas las transferencias es "
        f"{_es_int(transfer['mediana_span_intervalo'])} pb. "
        f"{_es_int(transfer['complejos_con_al_menos_un_cluster_por_gen'])} "
        "complejos tienen al menos un cluster con gen compartido; "
        f"{_es_int(transfer['complejos_con_al_menos_un_cluster_por_intervalo'])} "
        "tienen al menos un cluster con intervalo colocalizado.")
    lines.append("")
    lines.append("**Lectura:** el intervalo amplio es casi universal porque las "
                 "subunidades transferidas suelen estar dispersas; por eso se "
                 "reporta solo como diagnóstico. El intervalo colocalizado es la "
                 "lectura conservadora y aun así más permisiva que el gen. "
                 "Cualquier conclusión debe apoyarse en el gen.")
    lines.append("")

    lines.append("### 3.2 Con firmas K-12")
    lines.append("")
    lines.append(_md_table(
        ["Métrica", "Valor"],
        [
            ["Firmas K-12", _es_int(firmas["firmas"])],
            ["Firmas con match por gen", _es_int(firmas["firmas_con_match_por_gen"])],
            ["Firmas con match por intervalo colocalizado", _es_int(firmas["firmas_con_match_por_intervalo"])],
            ["Firmas con intervalo colocalizado", _es_int(firmas["firmas_con_intervalo_colocalizado"])],
            ["Firmas sin intervalo colocalizado", _es_int(firmas["firmas_sin_intervalo_colocalizado"])],
            ["Clusters con puente K-12", _es_int(firmas["clusters_con_puente_k12"])],
            ["Clusters con puente colocalizado", _es_int(firmas["clusters_con_puente_k12_colocalizado"])],
            ["Clusters con puente a complejo con `curado_id`", _es_int(global_stats["solapamiento"]["bridge_signature_backed"])],
            ["Clusters sin puente K-12", _es_int(firmas["clusters_sin_puente_k12"])],
            ["Clusters con match de firma por gen", _es_int(firmas["clusters_con_match_firma_por_gen"])],
            ["Clusters con match de firma por intervalo", _es_int(firmas["clusters_con_match_firma_por_intervalo"])],
        ]))
    lines.append("")
    lines.append(_md_table(
        ["Confianza", "Clusters con puente", "%", "Puente coloc.", "%",
         "Con gen firma", "%", "Con int. firma", "%"],
        [[conf,
          _es_int(por_conf[conf]["solapamiento"]["has_bridge"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_bridge"], 1) + " %",
          _es_int(por_conf[conf]["solapamiento"]["has_bridge_coloc"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_bridge_coloc"], 1) + " %",
          _es_int(por_conf[conf]["solapamiento"]["has_signature_gene"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_signature_gene"], 1) + " %",
          _es_int(por_conf[conf]["solapamiento"]["has_signature_iv"]),
          _es_dec(por_conf[conf]["solapamiento_pct"]["has_signature_iv"], 1) + " %"]
         for conf in CONF_ORDER if conf in por_conf]))
    lines.append("")
    lines.append(
        "El match con firmas exige un puente de transferencia; no es una "
        "comparación directa de secuencias. El match por intervalo exige "
        "además que la firma y la fila-puente tengan intervalo colocalizado, "
        "por eso es menor que el match por gen. Las firmas no cubren todos los "
        "clusters ni todos los complejos.")
    lines.append("")

    lines.append("### 3.3 Cobertura de gold K-12")
    lines.append("")
    lines.append(_md_table(
        ["Referencia", "Valor"],
        [
            ["Complejos gold K-12", _es_int(gold["complejos_gold"])],
            ["Gold con al menos una transferencia", _es_int(gold["complejos_con_transferencia"])],
            ["Complejos transferidos", _es_int(gold["complejos_transferidos"])],
            ["Transferidos que están en gold", _es_int(gold["complejos_transferidos_en_gold"])],
            ["Firmas con `curado_id`", _es_int(gold["firmas_con_curado_id"])],
            ["`curado_id` únicos", _es_int(gold["firmas_curado_id_unicos"])],
            ["`curado_id` presentes en gold", _es_int(gold["firmas_curado_id_en_gold"])],
            ["Firmas sin `curado_id`", _es_int(gold["firmas_sin_curado_id"])],
        ]))
    lines.append("")
    lines.append(
        "Gold K-12 es la única referencia curada; no hay gold para las otras "
        "especies. La cobertura de gold mide alcance del catálogo, no acierto "
        "de los clusters.")
    lines.append("")

    lines.append("## 4. Priorización para revisión")
    lines.append("")
    lines.append(
        f"**Fórmula:** `{prior['formula']}`. Es una heurística: suma confianza "
        f"(alta {_es_dec(prior['pesos']['confianza']['alta'], 1)}, media "
        f"{_es_dec(prior['pesos']['confianza']['media'], 1)}, baja "
        f"{_es_dec(prior['pesos']['confianza']['baja'], 1)}), solapamiento "
        "genético con transferencia, y componentes menores de firma, intervalo "
        "y tamaño. **No es una probabilidad.**")
    lines.append("")
    tiers = prior["conteo_tiers"]
    lines.append(_md_table(
        ["Tier", "Criterio", "Clusters"],
        [["P1", prior["tiers"]["P1"], _es_int(tiers.get("P1", 0))],
         ["P2", prior["tiers"]["P2"], _es_int(tiers.get("P2", 0))],
         ["P3", prior["tiers"]["P3"], _es_int(tiers.get("P3", 0))]]))
    lines.append("")
    anchored = [row for row in top if row["metrics"]["t_gene_cov"] > 0]
    shown = (anchored or top)[:top_rows]
    lines.append(
        f"Top {_es_int(len(shown))} de {_es_int(len(top))} candidatos en el JSON "
        "(filtrado a clusters con cobertura de gen de transferencia > 0). "
        "La tabla completa está en `priorizacion.top`.")
    lines.append("")
    rows = []
    for row in shown:
        metrics = row["metrics"]
        cluster_label = row["complex_id"].split(":")[-1]
        rows.append([
            _es_int(row["rank"]),
            _short(row["genome"], 28),
            _short(cluster_label, 18),
            row["confidence"],
            _es_int(row["size"]),
            _es_dec(metrics["t_gene_cov"], 2),
            _es_dec(metrics["t_gene_jac"], 2),
            _es_dec(metrics["t_iv_frac"], 2),
            _es_dec(metrics["s_gene_jac"], 2),
            _short(row.get("t_best_id") or "—", 14),
            _short(row.get("s_best_id") or "—", 20),
            _es_dec(metrics["priority_score"], 2),
        ])
    lines.append(_md_table(
        ["#", "Genoma", "Cluster", "Conf.", "Tam.", "Cob. gen",
         "Jacc. gen", "Int. coloc.", "Jacc. firma", "Transferencia",
         "Firma", "Score"],
        rows))
    lines.append("")
    lines.append(
        "`Transferencia` es el complejo K-12 de la mejor fila de transferencia "
        "por gen; `Firma`, la firma K-12 con mayor Jaccard. Ordenar por "
        "`priority_score` prioriza clusters con anclaje de transferencia y "
        "confianza alta, pero no garantiza que el cluster sea correcto. Un "
        "cluster sin transferencia puede seguir siendo una hipótesis válida; "
        "simplemente no tiene anclaje en este análisis.")
    lines.append("")

    lines.append("## 5. Limitaciones")
    lines.append("")
    for limitation in report["limitaciones"]:
        lines.append(f"- {limitation}")
    lines.append("")

    lines.append("## 6. Reproducir")
    lines.append("")
    lines.append("```bash")
    lines.append(".venv/bin/python -m pytest tests/test_calibrate_model_clusters.py -q")
    lines.append(".venv/bin/python tools/calibrate_model_clusters.py \\")
    lines.append("    --out-json reports/model_clusters_calibration.json \\")
    lines.append("    --out-md reports/model_clusters_calibration.md")
    lines.append("```")
    lines.append("")
    lines.append("El script no usa red y solo escribe los dos reportes indicados.")
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calibra clusters del modelo K-12 en los 400 genomas "
                    "(sin red, sin modificar datos).")
    parser.add_argument(
        "--complexes-dir", type=Path,
        default=ROOT / "data" / "complexes")
    parser.add_argument(
        "--regulation-dir", type=Path,
        default=ROOT / "data" / "regulation")
    parser.add_argument(
        "--signatures", type=Path,
        default=ROOT / "data" / "complexes" / "clusters_ecoli_k12_evidencia.tsv")
    parser.add_argument(
        "--gold", type=Path,
        default=ROOT / "data" / "complexes" / "gold_k12.tsv")
    parser.add_argument(
        "--k12-genes", type=Path,
        default=ROOT / "data" / "complexes" / "k12_genes.tsv")
    parser.add_argument(
        "--out-json", type=Path,
        default=ROOT / "reports" / "model_clusters_calibration.json")
    parser.add_argument(
        "--out-md", type=Path,
        default=ROOT / "reports" / "model_clusters_calibration.md")
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument(
        "--limit-genomes", type=int, default=0,
        help="Limita a los primeros N genomas (0 = todos).")
    parser.add_argument(
        "--bin-size", type=int, default=10_000,
        help="Tamaño de bin para el índice de intervalos de firmas.")
    parser.add_argument(
        "--max-colocal-span", type=int, default=DEFAULT_MAX_COLOCAL_SPAN,
        help="Span máximo (pb) para considerar colocalizado un intervalo.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    report = analyze(
        complexes_dir=args.complexes_dir,
        regulation_dir=args.regulation_dir,
        signatures_path=args.signatures,
        gold_path=args.gold,
        k12_genes_path=args.k12_genes,
        top_n=args.top_n,
        limit_genomes=args.limit_genomes,
        bin_size=args.bin_size,
        max_colocal_span=args.max_colocal_span,
    )

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")
    args.out_md.write_text(build_markdown(report), encoding="utf-8")

    totals = report["totales"]
    print(
        f"clusters={totals['clusters']} genomas={totals['genomas']} "
        f"confianza={totals['confianza']} "
        f"con_gen_transf={report['agregado_global']['solapamiento']['has_transfer_gene']} "
        f"con_intervalo_transf={report['agregado_global']['solapamiento']['has_transfer_iv']} "
        f"con_intervalo_amplio={report['agregado_global']['solapamiento']['has_transfer_iv_any']} "
        f"top={report['priorizacion']['top_n']}")
    print(f"JSON: {args.out_json}")
    print(f"MD:   {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
