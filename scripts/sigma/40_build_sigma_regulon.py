#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Construye la capa sigma -> operón -> regulón y módulos.

Une la asignación sigma por operón (`sigma_promoters.tsv`) con la estructura de
operones (`operons.tsv`, `operon_genes.tsv`, `genes.tsv`) y, si existe, con las
aristas TF->operón de `regulon.json`. El resultado es una tabla por genoma:

    data/sigma/regulon/<genoma>/sigma_regulon.tsv
    data/sigma/regulon/<genoma>/sigma_regulon_summary.json

y, con `--all`, una tabla global:

    data/sigma/regulon/sigma_regulon_all.tsv
    data/sigma/regulon/sigma_regulon_all_summary.json

La tabla **no** agrega dirección, fuerza ni validación experimental: solo
integra evidencia ya calculada por los pasos anteriores. Si un genoma no tiene
`regulon.json`, `tf_overlap`/`tf_names` quedan como `not_available`; no se
inventa solapamiento.

Uso:
    .venv/bin/python scripts/sigma/40_build_sigma_regulon.py --genome ecoli_k12
    .venv/bin/python scripts/sigma/40_build_sigma_regulon.py --all --jobs 8
    .venv/bin/python scripts/sigma/40_build_sigma_regulon.py --selftest

No usa red, no instala paquetes y no escribe fuera de `--out`
(por defecto `data/sigma/regulon/`).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable, Iterable

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

VERSION = "1.0.0"
DEFAULT_OUT = ROOT / "data" / "sigma" / "regulon"

OUTPUT_COLS = [
    "genome",
    "operon_id",
    "operon_label",
    "sigma",
    "sigma_class",
    "confidence_label",
    "sigma_status",
    "target_genes",
    "target_functions",
    "tf_overlap",
    "tf_names",
    "pathway_context",
    "evidence",
    "detected_by",
    "provenance",
]

# Solo se usa si `sigma_promoters.tsv` no trae `sigma_class` para un valor
# conocido. No reemplaza el valor de entrada cuando está presente.
SIGMA_CLASS_FALLBACK = {
    "sigma70": "housekeeping",
    "sigma54": "alternative",
    "sigma32": "stress",
    "sigma38": "stationary",
    "sigma28": "flagellar",
    "sigma24": "ecf",
    "sigma19": "ecf",
    "sigmaE": "ecf",
    "sigmaB": "stress",
    "sigmaH": "sporulation",
    "sigmaK": "sporulation",
    "sigmaF": "sporulation",
    "sigma_unresolved": "unknown",
}

# Límites de tamaño para que la tabla global de 400 genomas sea manejable.
# El recorte se declara explícitamente con "(+N más)" en `target_functions`.
MAX_FUNCTION_ENTRIES = 8
MAX_RHEA_PER_GENE = 2

STATUS_RANK = {"assigned": 0, "ambiguous": 1, "unresolved": 2}
REQUIRED_REG_FILES = ("genes.tsv", "operons.tsv", "operon_genes.tsv")


# --------------------------------------------------------------------------
# Lectura / escritura
# --------------------------------------------------------------------------
def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _cell(value: Any) -> str:
    """Normaliza una celda para TSV: sin tabs ni saltos de línea."""
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def _atomic_write_text(path: Path, writer: Callable[[Any], None]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp",
                               dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp)
    try:
        with tmp_path.open("w", encoding="utf-8", newline="") as fh:
            writer(fh)
        os.replace(tmp_path, path)
        try:
            os.chmod(path, 0o644)
        except OSError:
            pass
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def _write_tsv(path: Path, rows: Iterable[dict[str, Any]], cols: list[str]) -> None:
    def _writer(fh: Any) -> None:
        # Mismo contrato que `sigma_annotator._write_tsv`: TSV crudo, sin
        # comillas CSV. Las celdas ya se normalizan con `_cell`.
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(_cell(row.get(col, "")) for col in cols) + "\n")

    _atomic_write_text(path, _writer)


def _write_json(path: Path, obj: Any) -> None:
    def _writer(fh: Any) -> None:
        json.dump(obj, fh, ensure_ascii=False, indent=2)
        fh.write("\n")

    _atomic_write_text(path, _writer)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _input_hashes(paths: dict[str, Path]) -> dict[str, str]:
    return {name: (_sha256(p) if p.exists() else "not_available")
            for name, p in paths.items()}


def _dedupe(items: Iterable[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for item in items:
        if item and item not in seen:
            seen.add(item)
            out.append(item)
    return out


# --------------------------------------------------------------------------
# Parseo de campos de sigma_promoters.tsv
# --------------------------------------------------------------------------
def parse_pathway_hits(raw: Any) -> list[str]:
    """Convierte `pathway_hits` (JSON de lista) en términos limpios.

    Acepta lista JSON, string suelto, JSON inválido (se preserva como término)
    y valores vacíos. Deduplica preservando orden.
    """
    if raw is None:
        return []
    text = str(raw).strip()
    if not text:
        return []
    try:
        value = json.loads(text)
    except (json.JSONDecodeError, TypeError, ValueError):
        return _dedupe([_cell(text)])
    if isinstance(value, list):
        terms = [_cell(v) for v in value if _cell(v)]
    elif isinstance(value, str):
        terms = [_cell(value)] if _cell(value) else []
    elif value is None:
        terms = []
    else:
        terms = [_cell(value)]
    return _dedupe(terms)


def pathway_context(sigma_row: dict[str, Any]) -> str:
    """Devuelve `yes: term; term`, `no` o `not_available`.

    `yes`/`no` reflejan `pathway_support`; los términos salen de
    `pathway_hits`. Si hay hits pero el flag falta, se reporta `yes`.
    """
    support = _cell(sigma_row.get("pathway_support", "")).lower()
    terms = parse_pathway_hits(sigma_row.get("pathway_hits", ""))
    if support in {"yes", "true", "1", "si", "sí"}:
        return "yes: " + "; ".join(terms) if terms else "yes"
    if support in {"no", "false", "0"}:
        return "no"
    if terms:
        return "yes: " + "; ".join(terms)
    return "not_available"


def _gene_label(row: dict[str, Any]) -> str:
    for key in ("gene_name", "locus_tag", "gene_id"):
        value = _cell(row.get(key, ""))
        if value:
            return value
    return ""


# --------------------------------------------------------------------------
# Índices de entrada (funciones puras sobre filas)
# --------------------------------------------------------------------------
def _position_key(value: Any) -> tuple[int, int]:
    text = _cell(value)
    try:
        return (0, int(float(text)))
    except (TypeError, ValueError):
        return (1, 0)


def build_gene_index(rows: Iterable[dict[str, Any]]) -> dict[str, dict[str, str]]:
    """gene_id -> {gene_name, locus_tag, product, ec}; completa huecos."""
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        gid = _cell(row.get("gene_id", ""))
        if not gid:
            continue
        meta = index.setdefault(
            gid, {"gene_name": "", "locus_tag": "", "product": "", "ec": ""}
        )
        for key in ("gene_name", "locus_tag", "product", "ec"):
            value = _cell(row.get(key, ""))
            if value and not meta[key]:
                meta[key] = value
    return index


def build_operon_gene_index(
    rows: Iterable[dict[str, Any]],
) -> dict[str, list[dict[str, str]]]:
    """operon_id -> genes en orden de `position`, sin duplicados."""
    grouped: dict[str, list[tuple[tuple[int, int], int, dict[str, Any]]]] = defaultdict(list)
    for order, row in enumerate(rows):
        oid = _cell(row.get("operon_id", ""))
        if not oid:
            continue
        grouped[oid].append((_position_key(row.get("position", "")), order, row))

    index: dict[str, list[dict[str, str]]] = {}
    for oid, items in grouped.items():
        items.sort(key=lambda item: (item[0], item[1]))
        genes: list[dict[str, str]] = []
        seen: set[str] = set()
        for _, _, row in items:
            gene = {
                "gene_id": _cell(row.get("gene_id", "")),
                "gene_name": _cell(row.get("gene_name", "")),
                "locus_tag": _cell(row.get("locus_tag", "")),
                "contig": _cell(row.get("contig", "")),
                "start": _cell(row.get("start", "")),
                "end": _cell(row.get("end", "")),
            }
            key = gene["gene_id"] or gene["gene_name"]
            if not key:
                key = f"__row_{len(genes)}"
            if key in seen:
                continue
            seen.add(key)
            genes.append(gene)
        index[oid] = genes
    return index


def build_operon_index(rows: Iterable[dict[str, Any]]) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        oid = _cell(row.get("operon_id", ""))
        if oid and oid not in index:
            index[oid] = row
    return index


def build_operon_aliases(
    operon_gene_index: dict[str, list[dict[str, str]]],
    gene_index: dict[str, dict[str, str]],
) -> dict[str, set[str]]:
    """Alias (gene_id / gene_name / locus_tag, exacto y minúscula) -> operones.

    Un alias puede apuntar a más de un operón; el consumidor solo lo usa si el
    mapeo es inequívoco.
    """
    exact: dict[str, set[str]] = defaultdict(set)
    for oid, genes in operon_gene_index.items():
        for gene in genes:
            gid = gene["gene_id"]
            meta = gene_index.get(gid, {}) if gid else {}
            for value in (gid, gene["gene_name"], gene["locus_tag"],
                          meta.get("gene_name", ""), meta.get("locus_tag", "")):
                value = _cell(value)
                if value:
                    exact[value].add(oid)
    aliases: dict[str, set[str]] = defaultdict(set)
    for value, oids in exact.items():
        aliases[value].update(oids)
        aliases[value.lower()].update(oids)
    return aliases


def build_rhea_index(rows: Iterable[dict[str, Any]]) -> dict[str, dict[Any, list[str]]]:
    """Índice opcional de Rhea: por gene_id y por (contig, start, end)."""
    by_gene: dict[str, list[str]] = defaultdict(list)
    by_coord: dict[tuple[str, str, str], list[str]] = defaultdict(list)
    for row in rows:
        rhea = _cell(row.get("rhea_id", ""))
        if not rhea:
            continue
        gid = _cell(row.get("gene_id", ""))
        if gid:
            by_gene[gid].append(rhea)
        coord = (_cell(row.get("contig", "")), _cell(row.get("start", "")),
                 _cell(row.get("end", "")))
        if any(coord):
            by_coord[coord].append(rhea)
    return {
        "by_gene": {k: _dedupe(v) for k, v in by_gene.items()},
        "by_coord": {k: _dedupe(v) for k, v in by_coord.items()},
    }


def build_tf_index(
    edges: Iterable[dict[str, Any]],
    aliases: dict[str, set[str]],
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, int]]:
    """Agrupa aristas TF->operón por operón.

    Usa `target_operon` cuando existe. Si falta, intenta mapear `target_genes`
    (IDs o nombres) y solo acepta alias inequívocos. No inventa solapamiento.
    """
    index: dict[str, list[dict[str, Any]]] = defaultdict(list)
    stats = {
        "n_tf_edges_input": 0,
        "n_tf_edges_mapped": 0,
        "n_tf_edges_unmapped": 0,
        "n_tf_edges_ambiguous": 0,
        "n_tf_edges_without_tf": 0,
    }
    for edge in edges or []:
        if not isinstance(edge, dict):
            continue
        stats["n_tf_edges_input"] += 1
        tf = _cell(edge.get("tf", ""))
        if not tf:
            stats["n_tf_edges_without_tf"] += 1
            continue

        target_oids: set[str] = set()
        ambiguous = False
        oid = _cell(edge.get("target_operon", ""))
        if oid:
            target_oids.add(oid)
        else:
            for token in _cell(edge.get("target_genes", "")).split(","):
                token = token.strip()
                if not token:
                    continue
                hits = aliases.get(token) or aliases.get(token.lower()) or set()
                if len(hits) == 1:
                    target_oids.update(hits)
                elif len(hits) > 1:
                    ambiguous = True

        if not target_oids:
            if ambiguous:
                stats["n_tf_edges_ambiguous"] += 1
            else:
                stats["n_tf_edges_unmapped"] += 1
            continue
        stats["n_tf_edges_mapped"] += 1
        for target_oid in target_oids:
            index[target_oid].append(edge)
    return dict(index), stats


def unique_tfs(edges: Iterable[dict[str, Any]]) -> list[str]:
    names = {_cell(e.get("tf", "")) for e in edges if isinstance(e, dict)}
    names.discard("")
    return sorted(names, key=lambda value: (value.lower(), value))


# --------------------------------------------------------------------------
# Construcción de filas
# --------------------------------------------------------------------------
def _rhea_for_gene(gene: dict[str, Any],
                   rhea_index: dict[str, dict[Any, list[str]]]) -> list[str]:
    if not rhea_index:
        return []
    gid = _cell(gene.get("gene_id", ""))
    by_gene = rhea_index.get("by_gene", {})
    if gid and gid in by_gene:
        return by_gene[gid]
    coord = (_cell(gene.get("contig", "")), _cell(gene.get("start", "")),
             _cell(gene.get("end", "")))
    if any(coord):
        return rhea_index.get("by_coord", {}).get(coord, [])
    return []


def build_target_functions(
    genes: Iterable[dict[str, Any]],
    gene_index: dict[str, dict[str, str]],
    rhea_index: dict[str, dict[Any, list[str]]] | None = None,
) -> str:
    """`gen:producto [EC x] [RHEA:...]` por gen, deduplicado y acotado.

    Si no hay producto ni EC ni Rhea, se conserva el nombre del gen.
    """
    rhea_index = rhea_index or {}
    entries: list[str] = []
    for gene in genes:
        gid = _cell(gene.get("gene_id", ""))
        meta = gene_index.get(gid, {}) if gid else {}
        name = (_cell(meta.get("gene_name", ""))
                or _cell(gene.get("gene_name", ""))
                or _cell(meta.get("locus_tag", ""))
                or _cell(gene.get("locus_tag", ""))
                or gid)
        product = _cell(meta.get("product", ""))
        ec = _cell(meta.get("ec", ""))
        if not name and not product:
            continue
        entry = f"{name}:{product}" if name and product else (product or name)
        if ec:
            entry = f"{entry} [EC {ec}]"
        rheas = _rhea_for_gene(gene, rhea_index)
        if rheas:
            labels = [r if r.upper().startswith("RHEA:") else f"RHEA:{r}"
                      for r in rheas]
            shown = labels[:MAX_RHEA_PER_GENE]
            tag = ",".join(shown)
            if len(labels) > MAX_RHEA_PER_GENE:
                tag += f",+{len(labels) - MAX_RHEA_PER_GENE}"
            entry = f"{entry} [{tag}]"
        entries.append(entry)
    entries = _dedupe(entries)
    if not entries:
        return "not_available"
    if len(entries) > MAX_FUNCTION_ENTRIES:
        hidden = len(entries) - MAX_FUNCTION_ENTRIES
        entries = entries[:MAX_FUNCTION_ENTRIES] + [f"(+{hidden} más)"]
    return "; ".join(entries)


def _dedupe_sigma_rows(
    rows: Iterable[dict[str, Any]],
) -> tuple[list[dict[str, Any]], int]:
    """Una asignación por operón, con selección determinista si hay duplicados.

    Prioriza `sigma_status=assigned`, luego `confidence_score` y `score`.
    Conserva el orden de aparición de las filas elegidas.
    """
    def rank(item: tuple[int, dict[str, Any]]) -> tuple[int, float, float, int]:
        order, row = item
        status = _cell(row.get("sigma_status", "")).lower()
        try:
            conf = float(_cell(row.get("confidence_score", "")) or 0)
        except ValueError:
            conf = 0.0
        try:
            score = float(_cell(row.get("score", "")) or 0)
        except ValueError:
            score = 0.0
        return (STATUS_RANK.get(status, 3), -conf, -score, order)

    items = list(rows)
    best: dict[str, tuple[int, dict[str, Any]]] = {}
    for item in enumerate(items):
        oid = _cell(item[1].get("operon_id", ""))
        if not oid:
            continue
        if oid not in best or rank(item) < rank(best[oid]):
            best[oid] = item
    chosen = sorted(best.values(), key=lambda item: item[0])
    kept = [row for _, row in chosen]
    return kept, len(items) - len(kept)


def build_rows(
    genome: str,
    sigma_rows: Iterable[dict[str, Any]],
    operon_index: dict[str, dict[str, str]],
    operon_gene_index: dict[str, list[dict[str, str]]],
    gene_index: dict[str, dict[str, str]],
    tf_index: dict[str, list[dict[str, Any]]],
    tf_available: bool,
    rhea_index: dict[str, dict[Any, list[str]]] | None = None,
) -> list[dict[str, Any]]:
    """Filas de `sigma_regulon.tsv` (una por operón con asignación sigma)."""
    rows: list[dict[str, Any]] = []
    for sigma_row in sigma_rows:
        oid = _cell(sigma_row.get("operon_id", ""))
        if not oid:
            continue
        operon = operon_index.get(oid, {})
        genes = operon_gene_index.get(oid, [])

        target_genes = ",".join(label for label in (_gene_label(g) for g in genes) if label)
        if not target_genes:
            raw_names = _cell(operon.get("gene_names", ""))
            if raw_names:
                target_genes = ",".join(
                    name.strip() for name in raw_names.split(",") if name.strip()
                )
        if not target_genes:
            target_genes = _cell(sigma_row.get("first_gene_name", "")) or _cell(
                sigma_row.get("first_gene", "")
            )

        if tf_available:
            edges = tf_index.get(oid, [])
            tfs = unique_tfs(edges)
            tf_overlap: int | str = len(tfs)
            tf_names = "; ".join(tfs)
        else:
            tf_overlap = "not_available"
            tf_names = "not_available"

        # Si falta `operon_genes.tsv` para este operón, conserva al menos los
        # nombres de gen como "función desconocida"; no inventa productos.
        function_genes = genes
        if not function_genes and target_genes:
            function_genes = [
                {"gene_name": name.strip()}
                for name in target_genes.split(",") if name.strip()
            ]

        sigma = _cell(sigma_row.get("sigma", ""))
        sigma_class = (_cell(sigma_row.get("sigma_class", ""))
                       or SIGMA_CLASS_FALLBACK.get(sigma, "unknown"))
        operon_label = (_cell(sigma_row.get("operon_label", ""))
                        or _cell(operon.get("label", ""))
                        or _cell(sigma_row.get("first_gene_name", ""))
                        or oid)

        rows.append({
            "genome": genome,
            "operon_id": oid,
            "operon_label": operon_label,
            "sigma": sigma,
            "sigma_class": sigma_class,
            "confidence_label": (_cell(sigma_row.get("confidence_label", ""))
                                 or _cell(sigma_row.get("confidence", ""))),
            "sigma_status": (_cell(sigma_row.get("sigma_status", ""))
                             or _cell(sigma_row.get("status", ""))),
            "target_genes": target_genes,
            "target_functions": build_target_functions(
                function_genes, gene_index, rhea_index
            ),
            "tf_overlap": tf_overlap,
            "tf_names": tf_names,
            "pathway_context": pathway_context(sigma_row),
            "evidence": _cell(sigma_row.get("evidence", "")),
            "detected_by": _cell(sigma_row.get("detected_by", "")),
            "provenance": _cell(sigma_row.get("provenance", "")),
        })
    return rows


def _count_dict(counter: Counter) -> dict[str, int]:
    return {key: counter[key] for key in sorted(counter)}


def summarize(
    genome: str,
    rows: list[dict[str, Any]],
    operons: list[dict[str, Any]],
    sigma_rows_input: int,
    duplicate_count: int,
    tf_available: bool,
    tf_stats: dict[str, int],
    inputs: dict[str, str],
) -> dict[str, Any]:
    """Resumen por genoma (conteos y trazabilidad de insumos)."""
    sigma_counts = Counter(_cell(r.get("sigma", "")) for r in rows)
    class_counts = Counter(_cell(r.get("sigma_class", "")) for r in rows)
    conf_counts = Counter(_cell(r.get("confidence_label", "")) for r in rows)
    status_counts = Counter(_cell(r.get("sigma_status", "")) for r in rows)

    path_counts: Counter = Counter()
    for row in rows:
        context = _cell(row.get("pathway_context", ""))
        if context.startswith("yes"):
            path_counts["yes"] += 1
        elif context == "no":
            path_counts["no"] += 1
        else:
            path_counts["not_available"] += 1

    with_tf = 0
    without_tf = 0
    tf_not_available = 0
    tf_names: set[str] = set()
    for row in rows:
        overlap = row.get("tf_overlap")
        if overlap == "not_available":
            tf_not_available += 1
            continue
        try:
            value = int(overlap)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            value = 0
        if value > 0:
            with_tf += 1
            for name in _cell(row.get("tf_names", "")).split(";"):
                name = name.strip()
                if name:
                    tf_names.add(name)
        else:
            without_tf += 1

    operon_ids = {_cell(o.get("operon_id", "")) for o in operons}
    operon_ids.discard("")
    n_operons_total = len(operon_ids) if operon_ids else len(operons)

    return {
        "genome": genome,
        "version": VERSION,
        "n_edges": len(rows),
        "n_operons_total": n_operons_total,
        "n_operons_with_sigma": len(rows),
        "n_operons_without_sigma": max(0, n_operons_total - len(rows)),
        "n_sigma_rows_input": sigma_rows_input,
        "n_duplicates_resolved": duplicate_count,
        "sigma_counts": _count_dict(sigma_counts),
        "sigma_class_counts": _count_dict(class_counts),
        "confidence_counts": _count_dict(conf_counts),
        "sigma_status_counts": _count_dict(status_counts),
        "pathway_context_counts": _count_dict(path_counts),
        "tf_available": bool(tf_available),
        "n_tfs": len(tf_names),
        "tfs": sorted(tf_names, key=lambda value: (value.lower(), value)),
        "n_operons_with_tf_overlap": with_tf,
        "n_operons_sin_tf": without_tf,
        "n_operons_tf_not_available": tf_not_available,
        "tf_stats": dict(tf_stats),
        "inputs": inputs,
    }


# --------------------------------------------------------------------------
# Un genoma
# --------------------------------------------------------------------------
def _load_regulon(path: Path) -> tuple[bool, list[dict[str, Any]]]:
    """Devuelve (disponible, aristas). Un JSON inválido es un error explícito."""
    if not path.exists():
        return False, []
    try:
        obj = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"regulon.json inválido: {exc}") from exc
    if isinstance(obj, dict):
        edges = obj.get("edges", [])
    elif isinstance(obj, list):
        edges = obj
    else:
        raise ValueError("regulon.json no es un objeto ni una lista")
    if edges is None:
        edges = []
    if not isinstance(edges, list):
        raise ValueError("regulon.json: 'edges' no es una lista")
    return True, [e for e in edges if isinstance(e, dict)]


def _require(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"falta {label}: {path}")


def build_genome(data_root: Path, genome: str, out_root: Path) -> dict[str, Any]:
    """Construye y escribe la capa sigma->regulón de un genoma."""
    data_root = Path(data_root)
    out_root = Path(out_root)
    sigma_path = data_root / "sigma" / genome / "sigma_promoters.tsv"
    reg_dir = data_root / "regulation" / genome
    _require(sigma_path, "sigma_promoters.tsv")
    for name in REQUIRED_REG_FILES:
        _require(reg_dir / name, name)

    sigma_rows = _read_tsv(sigma_path)
    operons = _read_tsv(reg_dir / "operons.tsv")
    operon_genes = _read_tsv(reg_dir / "operon_genes.tsv")
    genes = _read_tsv(reg_dir / "genes.tsv")

    regulon_path = reg_dir / "regulon.json"
    tf_available, edges = _load_regulon(regulon_path)

    gene_map_path = (data_root / "metabolism" / "genemaps"
                     / f"{genome}_gene_map.tsv")
    rhea_rows = _read_tsv(gene_map_path)

    operon_index = build_operon_index(operons)
    operon_gene_index = build_operon_gene_index(operon_genes)
    gene_index = build_gene_index(genes)
    rhea_index = build_rhea_index(rhea_rows)
    aliases = build_operon_aliases(operon_gene_index, gene_index)
    tf_index, tf_stats = build_tf_index(edges, aliases) if tf_available else ({}, {
        "n_tf_edges_input": 0,
        "n_tf_edges_mapped": 0,
        "n_tf_edges_unmapped": 0,
        "n_tf_edges_ambiguous": 0,
        "n_tf_edges_without_tf": 0,
    })

    deduped, duplicate_count = _dedupe_sigma_rows(sigma_rows)
    rows = build_rows(
        genome=genome,
        sigma_rows=deduped,
        operon_index=operon_index,
        operon_gene_index=operon_gene_index,
        gene_index=gene_index,
        tf_index=tf_index,
        tf_available=tf_available,
        rhea_index=rhea_index,
    )

    inputs = _input_hashes({
        "sigma_promoters.tsv": sigma_path,
        "operons.tsv": reg_dir / "operons.tsv",
        "operon_genes.tsv": reg_dir / "operon_genes.tsv",
        "genes.tsv": reg_dir / "genes.tsv",
        "regulon.json": regulon_path,
        "gene_map.tsv": gene_map_path,
    })
    summary = summarize(
        genome=genome,
        rows=rows,
        operons=operons,
        sigma_rows_input=len(sigma_rows),
        duplicate_count=duplicate_count,
        tf_available=tf_available,
        tf_stats=tf_stats,
        inputs=inputs,
    )

    out_dir = out_root / genome
    _write_tsv(out_dir / "sigma_regulon.tsv", rows, OUTPUT_COLS)
    _write_json(out_dir / "sigma_regulon_summary.json", summary)
    return summary


def _run_one(payload: tuple[Path, str, Path]) -> dict[str, Any]:
    """Worker picklable para `--all --jobs N`."""
    data_root, genome, out_root = payload
    return build_genome(Path(data_root), genome, Path(out_root))


# --------------------------------------------------------------------------
# Lote global
# --------------------------------------------------------------------------
def discover_genomes(data_root: Path) -> list[str]:
    """Genomas con sigma_promoters.tsv y las tres tablas de regulación."""
    sigma_root = Path(data_root) / "sigma"
    if not sigma_root.exists():
        return []
    genomes: list[str] = []
    for entry in sorted(sigma_root.iterdir()):
        if not entry.is_dir():
            continue
        if not (entry / "sigma_promoters.tsv").exists():
            continue
        reg_dir = Path(data_root) / "regulation" / entry.name
        if not all((reg_dir / name).exists() for name in REQUIRED_REG_FILES):
            continue
        genomes.append(entry.name)
    return genomes


def write_global(
    out_root: Path,
    genomes: list[str],
    summaries: list[dict[str, Any]],
    failures: list[dict[str, str]],
) -> tuple[Path, Path]:
    """Escribe la tabla global por streaming y su resumen agregado."""
    out_root = Path(out_root)
    global_path = out_root / "sigma_regulon_all.tsv"

    def _writer(fh: Any) -> None:
        fh.write("\t".join(OUTPUT_COLS) + "\n")
        for genome in sorted(genomes):
            path = out_root / genome / "sigma_regulon.tsv"
            if not path.exists():
                continue
            with path.open(encoding="utf-8", errors="replace", newline="") as src:
                next(src, None)  # cabecera
                shutil.copyfileobj(src, fh)

    _atomic_write_text(global_path, _writer)

    sigma_counts: Counter = Counter()
    conf_counts: Counter = Counter()
    status_counts: Counter = Counter()
    path_counts: Counter = Counter()
    tfs: set[str] = set()
    for summary in summaries:
        sigma_counts.update(summary.get("sigma_counts", {}))
        conf_counts.update(summary.get("confidence_counts", {}))
        status_counts.update(summary.get("sigma_status_counts", {}))
        path_counts.update(summary.get("pathway_context_counts", {}))
        for name in summary.get("tfs", []) or []:
            if name:
                tfs.add(str(name))

    aggregate = {
        "version": VERSION,
        "n_genomes": len(summaries),
        "n_edges": sum(int(s.get("n_edges", 0)) for s in summaries),
        "sigma_counts": _count_dict(sigma_counts),
        "confidence_counts": _count_dict(conf_counts),
        "sigma_status_counts": _count_dict(status_counts),
        "pathway_context_counts": _count_dict(path_counts),
        "n_operons_with_tf_overlap": sum(
            int(s.get("n_operons_with_tf_overlap", 0)) for s in summaries
        ),
        "n_operons_sin_tf": sum(int(s.get("n_operons_sin_tf", 0)) for s in summaries),
        "n_operons_tf_not_available": sum(
            int(s.get("n_operons_tf_not_available", 0)) for s in summaries
        ),
        "n_tfs": len(tfs),
        "failures": failures,
    }
    summary_path = out_root / "sigma_regulon_all_summary.json"
    _write_json(summary_path, aggregate)
    return global_path, summary_path


# --------------------------------------------------------------------------
# Autoprueba y CLI
# --------------------------------------------------------------------------
def selftest() -> str:
    """Prueba determinista de joins y formato, sin tocar datos reales."""
    sigma_rows = [{
        "genome": "demo",
        "operon_id": "demo:1:+",
        "operon_label": "abc",
        "sigma": "sigma32",
        "sigma_class": "stress",
        "confidence_label": "sospecha_media",
        "sigma_status": "assigned",
        "confidence_score": "0.5",
        "score": "9.0",
        "first_gene_name": "abcA",
        "pathway_support": "yes",
        "pathway_hits": '["abcA:chaperone", "abcB:protease"]',
        "evidence": "PWM sigma32 z=2.2",
        "detected_by": "promoter_pwm:sigma32",
        "provenance": '{"window": 180}',
    }]
    operons = [{"operon_id": "demo:1:+", "label": "abc",
                "gene_names": "abcA,abcB"}]
    operon_genes = [
        {"operon_id": "demo:1:+", "position": "1", "gene_id": "g1",
         "gene_name": "abcA", "contig": "c1", "start": "10", "end": "20"},
        {"operon_id": "demo:1:+", "position": "2", "gene_id": "g2",
         "gene_name": "abcB", "contig": "c1", "start": "30", "end": "40"},
    ]
    genes = [
        {"gene_id": "g1", "gene_name": "abcA", "product": "chaperone A",
         "ec": "1.1.1.1"},
        {"gene_id": "g2", "gene_name": "abcB", "product": "protease B"},
    ]
    edges = [{"tf": "LexA", "target_operon": "demo:1:+"}]

    operon_index = build_operon_index(operons)
    operon_gene_index = build_operon_gene_index(operon_genes)
    gene_index = build_gene_index(genes)
    aliases = build_operon_aliases(operon_gene_index, gene_index)
    tf_index, stats = build_tf_index(edges, aliases)
    rows = build_rows(
        "demo", sigma_rows, operon_index, operon_gene_index, gene_index,
        tf_index, tf_available=True, rhea_index={},
    )
    assert len(rows) == 1, rows
    row = rows[0]
    assert row["target_genes"] == "abcA,abcB", row
    assert "abcA:chaperone A [EC 1.1.1.1]" in row["target_functions"], row
    assert row["tf_overlap"] == 1 and row["tf_names"] == "LexA", row
    assert row["pathway_context"] == "yes: abcA:chaperone; abcB:protease", row
    assert stats["n_tf_edges_mapped"] == 1, stats

    no_tf = build_rows(
        "demo", sigma_rows, operon_index, operon_gene_index, gene_index,
        {}, tf_available=False, rhea_index={},
    )
    assert no_tf[0]["tf_overlap"] == "not_available", no_tf[0]
    assert no_tf[0]["tf_names"] == "not_available", no_tf[0]

    summary = summarize(
        "demo", rows, operons, sigma_rows_input=1, duplicate_count=0,
        tf_available=True, tf_stats=stats, inputs={},
    )
    assert summary["n_edges"] == 1, summary
    assert summary["n_operons_with_tf_overlap"] == 1, summary
    assert summary["sigma_counts"] == {"sigma32": 1}, summary
    return f"OK sigma_regulon selftest (v{VERSION})"


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Construye la capa sigma -> operón -> regulón por genoma.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--genome", help="Genoma a procesar (p. ej. ecoli_k12).")
    ap.add_argument("--all", action="store_true",
                    help="Procesa todos los genomas con datos completos.")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT,
                    help="Directorio de salida (default: data/sigma/regulon).")
    ap.add_argument("--data-root", type=Path, default=ROOT / "data",
                    help="Raíz de datos de entrada (default: data).")
    ap.add_argument("--jobs", type=int, default=1,
                    help="Procesos en paralelo para --all (default: 1).")
    ap.add_argument("--limit", type=int, default=0,
                    help="Máximo de genomas en --all (0 = todos).")
    ap.add_argument("--offset", type=int, default=0,
                    help="Desplazamiento sobre la lista de genomas en --all.")
    ap.add_argument("--selftest", action="store_true",
                    help="Ejecuta una autoprueba sin datos reales.")
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    if args.selftest:
        print(selftest())
        return 0

    data_root = args.data_root.resolve()
    out_root = args.out.resolve()

    if args.all and args.genome:
        print("error: usa --genome o --all, no ambos", file=sys.stderr)
        return 2
    if not args.all and not args.genome:
        print("error: indica --genome, --all o --selftest", file=sys.stderr)
        return 2

    if args.all:
        genomes = discover_genomes(data_root)
        if args.offset:
            genomes = genomes[args.offset:]
        if args.limit:
            genomes = genomes[:args.limit]
        if not genomes:
            print(f"error: no hay genomas procesables en {data_root}", file=sys.stderr)
            return 1

        t0 = time.time()
        summaries: list[dict[str, Any]] = []
        failures: list[dict[str, str]] = []
        jobs = max(1, int(args.jobs))
        if jobs == 1:
            for index, genome in enumerate(genomes, 1):
                try:
                    summary = build_genome(data_root, genome, out_root)
                    summaries.append(summary)
                    print(f"[{index}/{len(genomes)}] {genome}: "
                          f"{summary['n_edges']} edges, "
                          f"{summary['n_operons_with_tf_overlap']} con TF",
                          flush=True)
                except Exception as exc:  # noqa: BLE001
                    failures.append({"genome": genome, "error": str(exc)})
                    print(f"FALLO {genome}: {exc}", file=sys.stderr, flush=True)
        else:
            payloads = [(data_root, genome, out_root) for genome in genomes]
            with ProcessPoolExecutor(max_workers=jobs) as executor:
                futures = {executor.submit(_run_one, payload): payload[1]
                           for payload in payloads}
                for index, future in enumerate(as_completed(futures), 1):
                    genome = futures[future]
                    try:
                        summary = future.result()
                        summaries.append(summary)
                        print(f"[{index}/{len(genomes)}] {genome}: "
                              f"{summary['n_edges']} edges, "
                              f"{summary['n_operons_with_tf_overlap']} con TF",
                              flush=True)
                    except Exception as exc:  # noqa: BLE001
                        failures.append({"genome": genome, "error": str(exc)})
                        print(f"FALLO {genome}: {exc}", file=sys.stderr, flush=True)

        summaries.sort(key=lambda item: item["genome"])
        if summaries:
            global_path, summary_path = write_global(
                out_root, genomes, summaries, failures
            )
        else:
            global_path = summary_path = None
        elapsed = time.time() - t0
        print(f"lote: {len(summaries)}/{len(genomes)} genomas en {elapsed:.1f}s")
        if global_path:
            print(f"global: {global_path}")
            print(f"resumen global: {summary_path}")
        if failures:
            print(f"fallos: {len(failures)} (ver stderr)", file=sys.stderr)
        return 0 if summaries else 1

    genome = args.genome
    try:
        summary = build_genome(data_root, genome, out_root)
    except Exception as exc:  # noqa: BLE001
        print(f"FALLO {genome}: {exc}", file=sys.stderr)
        return 1
    print(f"{genome}: {summary['n_edges']} edges · "
          f"{summary['n_operons_with_tf_overlap']} con TF · "
          f"{summary['n_operons_sin_tf']} sin TF")
    print(f"salida: {out_root / genome}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
