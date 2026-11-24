#!/usr/bin/env python3
"""Auditoría de cobertura de la capa móvil / defensa / CRISPR / AMR / virulencia.

Lee lo que ya existe en disco y **no ejecuta** HMMER, DefenseFinder,
CRISPRCasFinder, AMRFinderPlus, ISEScan ni ninguna base externa. Tampoco crea,
completa ni recoloca anotaciones: si un genoma no tiene datos, lo marca como
``not_available`` con el motivo concreto.

Qué mide
--------
1. Inventario real de genomas: lista de espera (``full_genomes.tsv`` o
   ``pilot_genomes.tsv``) y genomas descubiertos por sus artefactos.
2. Cobertura por capa: TSV tipado, GFF3, espaciadores CRISPR, capa previa de
   ``data/mobile_defense/features.tsv`` y FASTA de entrada.
3. Tipos agregados: ``mobile_element``, ``prophage``, ``defense``, ``crispr``,
   ``amr`` y ``virulence``, mapeados desde los ``element_type`` reales.
4. Atribución al pipeline: qué script de ``scripts/mobile/`` produce cada
   artefacto (por convención del repo; los TSV no embeben el generador).
5. Brechas para el inventario: qué genomas no tienen TSV/GFF y en cuáles no se
   pudo evaluar CRISPR por falta de FASTA.

Uso
---
    python tools/mobile_coverage.py
    python tools/mobile_coverage.py --data data --json reports/mobile_coverage.json \
        --md reports/mobile_coverage.md

Los informes por defecto se escriben en ``reports/mobile_coverage.json`` y
``reports/mobile_coverage.md``. El resto del árbol se trata como solo lectura.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "data"
DEFAULT_JSON = ROOT / "reports" / "mobile_coverage.json"
DEFAULT_MD = ROOT / "reports" / "mobile_coverage.md"

SCHEMA = "mobile_coverage/v1"
CATEGORIES = ("mobile_element", "prophage", "defense", "crispr", "amr", "virulence")
# Normalización explícita de los `element_type` del TSV y de los `type` de la
# capa previa. No añade tipos: solo agrupa nombres ya presentes.
TYPE_TO_CATEGORY = {
    "mobile_element": "mobile_element",
    "transposon": "mobile_element",
    "prophage": "prophage",
    "defense": "defense",
    "crispr": "crispr",
    "cas": "crispr",
    "amr": "amr",
    "virulence": "virulence",
}
CATEGORY_TYPES = {
    "mobile_element": ("mobile_element", "transposon"),
    "prophage": ("prophage",),
    "defense": ("defense",),
    "crispr": ("crispr", "cas"),
    "amr": ("amr",),
    "virulence": ("virulence",),
}
# ---------------------------------------------------------------------------
# Utilidades de lectura (tolerantes: un archivo ausente no rompe la auditoría)
# ---------------------------------------------------------------------------
def _iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    if not path.is_file():
        return
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None:
                return
            for row in reader:
                if row:
                    yield row
    except OSError:
        return


def _read_tsv(path: Path) -> list[dict[str, str]]:
    return list(_iter_tsv(path))


def _header(path: Path) -> list[str]:
    if not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            return fh.readline().rstrip("\r\n").split("\t")
    except OSError:
        return []


def _has_columns(path: Path, *columns: str) -> bool:
    header = set(_header(path))
    return all(col in header for col in columns)


def _count_rows(path: Path) -> int:
    """Filas de datos de un TSV (cabecera excluida). No valida el contenido."""
    if not path.is_file():
        return 0
    n = 0
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            for i, _ in enumerate(fh):
                if i:
                    n += 1
    except OSError:
        return 0
    return n


def _rel(path: Path | None, data_root: Path) -> str | None:
    if path is None:
        return None
    try:
        return path.resolve().relative_to(data_root.resolve()).as_posix()
    except (OSError, ValueError):
        return path.as_posix()


def _counter_dict(counter: Counter) -> dict[str, int]:
    return {key: counter[key] for key in sorted(counter, key=lambda k: (-counter[k], k))}


def _category_counts(counter: Counter) -> dict[str, int]:
    return {category: int(counter.get(category, 0)) for category in CATEGORIES}


def _pct(part: int | None, total: int | None) -> float | None:
    if part is None or not total:
        return None
    return round(100.0 * part / total, 2)


def _dedupe(items: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(str(item) for item in items if item))


def _layer(status: str, reason: str | None = None, **extra) -> dict:
    out = {"status": status, "reason": reason}
    out.update(extra)
    return out


# ---------------------------------------------------------------------------
# Diseño de rutas (full vs piloto) e inventario
# ---------------------------------------------------------------------------
def _layout(data_root: Path) -> dict:
    base = data_root / "mobile_typing"
    if (base / "full").is_dir() or (base / "full_genomes.tsv").exists():
        mode = "full"
        typed_dir = base / "full"
        gff_dir = base / "gff" / "full"
        spacer_dir = base / "crispr_spacers" / "full"
        work_dir = base / "_work" / "full"
        inventory = base / "full_genomes.tsv"
        summary = base / "summary_full.json"
        validation = base / "validation_full.json"
    elif base.is_dir():
        mode = "pilot"
        typed_dir = base
        gff_dir = base / "gff"
        spacer_dir = base / "crispr_spacers"
        work_dir = base / "_work"
        inventory = base / "pilot_genomes.tsv"
        summary = base / "summary.json"
        validation = base / "validation.json"
    else:
        mode = "none"
        typed_dir = base
        gff_dir = base / "gff"
        spacer_dir = base / "crispr_spacers"
        work_dir = base / "_work"
        inventory = base / "full_genomes.tsv"
        summary = base / "summary_full.json"
        validation = base / "validation_full.json"
    return {
        "mode": mode,
        "base": base,
        "typed_dir": typed_dir,
        "gff_dir": gff_dir,
        "spacer_dir": spacer_dir,
        "work_dir": work_dir,
        "inventory": inventory,
        "summary": summary,
        "validation": validation,
    }


def _load_inventory(data_root: Path, layout: dict) -> dict:
    base = data_root / "mobile_typing"
    candidates = _dedupe(
        [
            layout["inventory"],
            base / "full_genomes.tsv",
            base / "pilot_genomes.tsv",
        ]
    )
    empty_candidates: list[str] = []
    for raw in candidates:
        path = Path(raw)
        if not path.is_file():
            continue
        rows = _read_tsv(path)
        genomes = _dedupe(row.get("genome", "") for row in rows)
        if not genomes:
            empty_candidates.append(_rel(path, data_root) or path.as_posix())
            continue
        return {
            "source": _rel(path, data_root),
            "expected_genomes": len(genomes),
            "genomes": genomes,
            "rows": len(rows),
            "columns": list(rows[0].keys()) if rows else [],
            "meta": {row["genome"]: row for row in rows if row.get("genome")},
            "empty_candidates": empty_candidates,
        }
    return {
        "source": None,
        "expected_genomes": None,
        "genomes": [],
        "rows": 0,
        "columns": [],
        "meta": {},
        "empty_candidates": empty_candidates,
    }


def _discover(data_root: Path, layout: dict) -> dict:
    base = data_root / "mobile_typing"
    typed_full = {
        path.stem
        for path in (base / "full").glob("*.tsv")
        if _has_columns(path, "element_type")
    }
    typed_pilot = {
        path.stem
        for path in base.glob("*.tsv")
        if _has_columns(path, "element_type")
    }
    gff_full = {path.stem for path in (base / "gff" / "full").glob("*.gff3")}
    gff_pilot = {path.stem for path in (base / "gff").glob("*.gff3")}
    spacer_full = {path.stem for path in (base / "crispr_spacers" / "full").glob("*.tsv")}
    spacer_pilot = {path.stem for path in (base / "crispr_spacers").glob("*.tsv")}

    prior_path = data_root / "mobile_defense" / "features.tsv"
    prior_genomes = _dedupe(
        row.get("genome", "") for row in _iter_tsv(prior_path) if row.get("genome")
    )

    ann_root = data_root / "annotations_ab_initio"
    annotation_genomes: list[str] = []
    if ann_root.is_dir():
        annotation_genomes = sorted(
            path.name
            for path in ann_root.iterdir()
            if path.is_dir() and (path / "annotation.gff3").is_file()
        )

    # FASTA en data/raw; los 90 genomas pequeños/endosimbiontes se movieron a
    # data/raw_excluded. Se buscan ambas ubicaciones sin copiar ni mover nada.
    suffix = "_genomic.fna.gz"
    fasta_map: dict[str, tuple[Path, str]] = {}
    for root_name in ("raw", "raw_excluded"):
        root = data_root / root_name
        if not root.is_dir():
            continue
        for path in sorted(root.glob(f"*{suffix}")):
            genome = path.name[: -len(suffix)]
            fasta_map.setdefault(genome, (path, root_name))
    fasta_genomes = sorted(fasta_map)

    all_genomes = _dedupe(
        sorted(
            typed_full
            | typed_pilot
            | gff_full
            | gff_pilot
            | spacer_full
            | spacer_pilot
            | set(prior_genomes)
            | set(annotation_genomes)
            | set(fasta_genomes)
        )
    )
    return {
        "typed_full": typed_full,
        "typed_pilot": typed_pilot,
        "gff_full": gff_full,
        "gff_pilot": gff_pilot,
        "spacer_full": spacer_full,
        "spacer_pilot": spacer_pilot,
        "prior_genomes": set(prior_genomes),
        "annotation_genomes": set(annotation_genomes),
        "fasta_genomes": set(fasta_genomes),
        "fasta_map": fasta_map,
        "all_genomes": all_genomes,
        "prior_path": prior_path,
        "ann_root": ann_root,
    }


def _load_prior(discovery: dict) -> tuple[dict[str, Counter], Counter]:
    by_genome: dict[str, Counter] = defaultdict(Counter)
    by_type: Counter = Counter()
    for row in _iter_tsv(discovery["prior_path"]):
        genome = (row.get("genome") or "").strip()
        etype = (row.get("type") or "unknown").strip() or "unknown"
        if not genome:
            continue
        by_genome[genome][etype] += 1
        by_type[etype] += 1
    return by_genome, by_type


def _pick_path(data_root: Path, layout: dict, genome: str, kind: str) -> tuple[Path | None, list[Path]]:
    """Primera ruta existente y alternativas. ``full`` tiene prioridad."""
    base = data_root / "mobile_typing"
    if kind == "typed":
        full = base / "full" / f"{genome}.tsv"
        pilot = base / f"{genome}.tsv"
    elif kind == "gff":
        full = base / "gff" / "full" / f"{genome}.gff3"
        pilot = base / "gff" / f"{genome}.gff3"
    elif kind == "spacers":
        full = base / "crispr_spacers" / "full" / f"{genome}.tsv"
        pilot = base / "crispr_spacers" / f"{genome}.tsv"
    else:  # pragma: no cover - defensa ante llamadas inválidas
        raise ValueError(kind)
    if layout["mode"] == "pilot":
        order = [pilot, full]
    else:
        order = [full, pilot]
    existing = [path for path in _dedupe(order) if Path(path).is_file()]
    primary = Path(existing[0]) if existing else None
    alternates = [Path(path) for path in existing[1:]]
    return primary, alternates


def _file_bytes(path: Path | None) -> int | None:
    if path is None or not path.is_file():
        return None
    try:
        return path.stat().st_size
    except OSError:
        return None


# ---------------------------------------------------------------------------
# Cálculo de cobertura
# ---------------------------------------------------------------------------
def _inspect_genome(
    genome: str,
    data_root: Path,
    layout: dict,
    inventory: dict,
    discovery: dict,
    prior_by_genome: dict[str, Counter],
) -> dict:
    typed_path, typed_alt = _pick_path(data_root, layout, genome, "typed")
    gff_path, gff_alt = _pick_path(data_root, layout, genome, "gff")
    spacer_path, spacer_alt = _pick_path(data_root, layout, genome, "spacers")

    by_type: Counter = Counter()
    by_category: Counter = Counter()
    unknown_types: Counter = Counter()
    typed_rows = 0
    if typed_path is not None:
        for row in _iter_tsv(typed_path):
            typed_rows += 1
            etype = (row.get("element_type") or "unknown").strip() or "unknown"
            by_type[etype] += 1
            category = TYPE_TO_CATEGORY.get(etype)
            if category is None:
                unknown_types[etype] += 1
            else:
                by_category[category] += 1

    if typed_rows > 0:
        typed_layer = _layer(
            "available",
            None,
            path=_rel(typed_path, data_root),
            alternate_paths=[_rel(path, data_root) for path in typed_alt],
            rows=typed_rows,
            by_type=_counter_dict(by_type),
            by_category=_category_counts(by_category),
            unknown_types=_counter_dict(unknown_types),
        )
        genome_status = "available"
    elif typed_path is not None:
        typed_layer = _layer(
            "not_available",
            "TSV presente pero sin filas de datos",
            path=_rel(typed_path, data_root),
            alternate_paths=[_rel(path, data_root) for path in typed_alt],
            rows=0,
            by_type={},
            by_category=_category_counts(Counter()),
            unknown_types={},
        )
        genome_status = "not_available"
    else:
        typed_layer = _layer(
            "not_available",
            "sin TSV tipado en data/mobile_typing/full ni data/mobile_typing",
            path=None,
            alternate_paths=[],
            rows=0,
            by_type={},
            by_category=_category_counts(Counter()),
            unknown_types={},
        )
        genome_status = "not_available"

    gff_bytes = _file_bytes(gff_path)
    if gff_path is not None and gff_bytes:
        gff_layer = _layer(
            "available",
            None,
            path=_rel(gff_path, data_root),
            alternate_paths=[_rel(path, data_root) for path in gff_alt],
            bytes=gff_bytes,
        )
    else:
        gff_layer = _layer(
            "not_available",
            "GFF3 ausente" if gff_path is None else "GFF3 vacío",
            path=_rel(gff_path, data_root),
            alternate_paths=[_rel(path, data_root) for path in gff_alt],
            bytes=gff_bytes,
        )

    fasta_entry = discovery["fasta_map"].get(genome)
    fasta_path = fasta_entry[0] if fasta_entry else None
    fasta_source = fasta_entry[1] if fasta_entry else ""
    fasta_bytes = _file_bytes(fasta_path) if fasta_path else 0
    if fasta_bytes:
        fasta_layer = _layer(
            "available",
            None,
            path=_rel(fasta_path, data_root),
            bytes=fasta_bytes,
            source=fasta_source,
        )
        fasta_ok = True
    else:
        fasta_layer = _layer(
            "not_available",
            "sin FASTA en data/raw ni data/raw_excluded",
            path=None,
            bytes=0,
        )
        fasta_ok = False

    spacer_rows = _count_rows(spacer_path) if spacer_path is not None else 0
    if not fasta_ok:
        crispr_layer = _layer(
            "not_available",
            "FASTA ausente: CRISPR de novo no evaluado; no implica ausencia de arrays",
            path=_rel(spacer_path, data_root),
            alternate_paths=[_rel(path, data_root) for path in spacer_alt],
            rows=spacer_rows,
            evaluated=False,
        )
    elif spacer_path is None:
        crispr_layer = _layer(
            "not_available",
            "FASTA disponible pero falta el TSV de espaciadores",
            path=None,
            alternate_paths=[],
            rows=0,
            evaluated=False,
        )
    else:
        crispr_layer = _layer(
            "available",
            None if spacer_rows else "evaluado: 0 espaciadores",
            path=_rel(spacer_path, data_root),
            alternate_paths=[_rel(path, data_root) for path in spacer_alt],
            rows=spacer_rows,
            evaluated=True,
        )

    prior_counter = prior_by_genome.get(genome, Counter())
    prior_rows = sum(prior_counter.values())
    prior_category = Counter()
    for etype, count in prior_counter.items():
        category = TYPE_TO_CATEGORY.get(etype)
        if category is not None:
            prior_category[category] += count
    if prior_rows:
        prior_layer = _layer(
            "available",
            None,
            path=_rel(discovery["prior_path"], data_root),
            rows=prior_rows,
            by_type=_counter_dict(prior_counter),
            by_category=_category_counts(prior_category),
        )
    else:
        prior_layer = _layer(
            "not_available",
            "sin filas para este genoma en data/mobile_defense/features.tsv",
            path=_rel(discovery["prior_path"], data_root),
            rows=0,
            by_type={},
            by_category=_category_counts(Counter()),
        )

    layers = {
        "typed_mobile": typed_layer,
        "gff": gff_layer,
        "crispr_scan": crispr_layer,
        "prior_features": prior_layer,
        "fasta": fasta_layer,
    }
    not_available = {
        name: layer["reason"]
        for name, layer in layers.items()
        if layer["status"] != "available"
    }

    # `prior_features` es una capa auxiliar (sin contig/hebra): su ausencia se
    # reporta en `not_available`, pero no convierte en parcial una cobertura
    # tipada que sí tiene TSV, GFF y CRISPR evaluable.
    required_unavailable = [
        name
        for name in ("typed_mobile", "gff", "crispr_scan")
        if layers[name]["status"] != "available"
    ]
    if genome_status == "not_available":
        coverage = "none"
        if prior_rows:
            reason = (
                "sin TSV tipado; solo capa previa "
                f"({prior_rows} filas, sin contig ni hebra)"
            )
        elif genome in inventory["meta"]:
            reason = "genoma del inventario sin salida de mobile_typing"
        else:
            reason = "sin datos tipados ni previos"
    elif required_unavailable:
        coverage = "partial"
        reason = None
    else:
        coverage = "full"
        reason = None

    meta = inventory["meta"].get(genome, {})
    inputs = {
        "annotation_gff": (discovery["ann_root"] / genome / "annotation.gff3").is_file()
        if discovery["ann_root"].is_dir()
        else False,
        "proteins": (discovery["ann_root"] / genome / "proteins.faa").is_file()
        if discovery["ann_root"].is_dir()
        else False,
        "ncbi_gff": (data_root / "annotations" / f"{genome}.gff.gz").is_file(),
        "fasta": fasta_ok,
    }

    return {
        "genome": genome,
        "in_inventory": genome in inventory["meta"],
        "phylum": meta.get("phylum") or "unknown",
        "organism": meta.get("organism") or "unknown",
        "has_fasta": fasta_ok,
        "has_fasta_declared": (meta.get("has_fasta") == "yes") if meta else None,
        "status": genome_status,
        "coverage": coverage,
        "reason": reason,
        "not_available": not_available,
        "typed_rows": typed_rows,
        "typed_by_category": _category_counts(by_category),
        "typed_by_type": _counter_dict(by_type),
        "unknown_types": _counter_dict(unknown_types),
        "prior_rows": prior_rows,
        "prior_by_category": _category_counts(prior_category),
        "inputs": inputs,
        "layers": layers,
    }


def _declared_artifacts(data_root: Path, layout: dict) -> dict:
    summary_path = layout["summary"]
    validation_path = layout["validation"]
    declared = {
        "summary": {
            "path": _rel(summary_path, data_root),
            "present": summary_path.is_file(),
            "genomes_processed": None,
            "rows_total": None,
            "genomes_with_fasta": None,
        },
        "validation": {
            "path": _rel(validation_path, data_root),
            "present": validation_path.is_file(),
            "status": None,
            "errors": None,
        },
    }
    if summary_path.is_file():
        try:
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            coverage = payload.get("coverage", {}) if isinstance(payload, dict) else {}
            counts = payload.get("counts", {}) if isinstance(payload, dict) else {}
            declared["summary"].update(
                {
                    "genomes_processed": coverage.get("genomes_processed"),
                    "rows_total": counts.get("rows_total"),
                    "genomes_with_fasta": coverage.get("genomes_with_fasta"),
                }
            )
        except (OSError, ValueError):
            declared["summary"]["error"] = "no se pudo leer el summary"
    if validation_path.is_file():
        try:
            payload = json.loads(validation_path.read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                errors = payload.get("errors")
                declared["validation"].update(
                    {
                        "status": payload.get("status"),
                        "errors": len(errors) if isinstance(errors, list) else errors,
                    }
                )
        except (OSError, ValueError):
            declared["validation"]["error"] = "no se pudo leer la validación"
    return declared


def _stage_specs(layout: dict) -> list[dict]:
    # Rutas relativas a `--data`; así el mismo glob funciona con data/ y con un
    # árbol sintético en tests.
    base_rel = "mobile_typing"
    full = layout["mode"] != "pilot"
    if full:
        return [
            {
                "stage": "prepare",
                "script": "scripts/mobile/19_prepare_full.py",
                "command": "python scripts/mobile/19_prepare_full.py",
                "outputs": [f"{base_rel}/full_genomes.tsv"],
            },
            {
                "stage": "extract",
                "script": "scripts/mobile/01_extract_candidates.py",
                "command": "python scripts/mobile/01_extract_candidates.py --full --shard 1..3 --shards 3",
                "outputs": [f"{base_rel}/_work/full/*/candidates.tsv"],
            },
            {
                "stage": "pfam",
                "script": "scripts/mobile/20_pfam_pool.py",
                "command": "python scripts/mobile/20_pfam_pool.py --workers 5",
                "outputs": [f"{base_rel}/_work/full/*/pfam_hits.tsv"],
            },
            {
                "stage": "crispr",
                "script": "scripts/mobile/21_crispr_pool.py",
                "command": "python scripts/mobile/21_crispr_pool.py --workers 5",
                "outputs": [f"{base_rel}/_work/full/*/crispr_arrays.tsv"],
                "needs_fasta": True,
            },
            {
                "stage": "classify",
                "script": "scripts/mobile/05_classify.py",
                "command": "python scripts/mobile/05_classify.py --full",
                "outputs": [f"{base_rel}/_work/full/*/elements.tsv"],
            },
            {
                "stage": "write",
                "script": "scripts/mobile/06_write_outputs.py",
                "command": "python scripts/mobile/06_write_outputs.py --full",
                "outputs": [
                    f"{base_rel}/full/*.tsv",
                    f"{base_rel}/gff/full/*.gff3",
                    f"{base_rel}/crispr_spacers/full/*.tsv",
                ],
            },
            {
                "stage": "validate",
                "script": "scripts/mobile/07_validate.py",
                "command": "python scripts/mobile/07_validate.py --full",
                "outputs": [f"{base_rel}/validation_full.json"],
            },
        ]
    return [
        {
            "stage": "select",
            "script": "scripts/mobile/00_select_genomes.py",
            "command": "python scripts/mobile/00_select_genomes.py",
            "outputs": [f"{base_rel}/pilot_genomes.tsv"],
        },
        {
            "stage": "extract",
            "script": "scripts/mobile/01_extract_candidates.py",
            "command": "python scripts/mobile/01_extract_candidates.py",
            "outputs": [f"{base_rel}/_work/*/candidates.tsv"],
        },
        {
            "stage": "pfam_panel",
            "script": "scripts/mobile/02_build_pfam_panel.py",
            "command": "python scripts/mobile/02_build_pfam_panel.py",
            "outputs": [f"{base_rel}/_work/pfam_panel/panel.hmm"],
        },
        {
            "stage": "pfam_scan",
            "script": "scripts/mobile/03_scan_pfam.py",
            "command": "python scripts/mobile/03_scan_pfam.py",
            "outputs": [f"{base_rel}/_work/*/pfam_hits.tsv"],
        },
        {
            "stage": "crispr",
            "script": "scripts/mobile/04_crispr_arrays.py",
            "command": "python scripts/mobile/04_crispr_arrays.py",
            "outputs": [f"{base_rel}/_work/*/crispr_arrays.tsv"],
            "needs_fasta": True,
        },
        {
            "stage": "classify",
            "script": "scripts/mobile/05_classify.py",
            "command": "python scripts/mobile/05_classify.py",
            "outputs": [f"{base_rel}/_work/*/elements.tsv"],
        },
        {
            "stage": "write",
            "script": "scripts/mobile/06_write_outputs.py",
            "command": "python scripts/mobile/06_write_outputs.py",
            "outputs": [
                f"{base_rel}/*.tsv",
                f"{base_rel}/gff/*.gff3",
                f"{base_rel}/crispr_spacers/*.tsv",
            ],
        },
        {
            "stage": "validate",
            "script": "scripts/mobile/07_validate.py",
            "command": "python scripts/mobile/07_validate.py",
            "outputs": [f"{base_rel}/validation.json"],
        },
    ]


def _glob_count(data_root: Path, pattern: str) -> int:
    return sum(1 for path in data_root.glob(pattern) if path.is_file())


def _stage_state(found: list[int], expected: int | None) -> str:
    if not found:
        return "missing"
    if expected is None:
        return "present" if any(found) else "missing"
    if expected == 0:
        return "not_applicable"
    if all(count >= expected for count in found):
        return "present"
    return "partial" if any(found) else "missing"


def _pipeline_report(
    data_root: Path,
    layout: dict,
    inventory: dict,
    discovery: dict,
    coverage: dict,
    declared: dict,
) -> dict:
    n_expected = inventory["expected_genomes"]
    n_fasta = coverage["genomes_with_fasta"]
    stages = []
    for spec in _stage_specs(layout):
        outputs = []
        for pattern in spec["outputs"]:
            outputs.append(
                {
                    "path": pattern,
                    "found": _glob_count(data_root, pattern),
                }
            )
        if spec["stage"] == "prepare" or spec["stage"] == "select":
            expected = 1
        elif spec.get("needs_fasta"):
            expected = n_fasta if n_expected is not None else None
        elif spec["stage"] == "validate":
            expected = 1
        else:
            expected = n_expected if n_expected is not None else len(discovery["all_genomes"])
        state = _stage_state([item["found"] for item in outputs], expected)
        stages.append(
            {
                "stage": spec["stage"],
                "script": spec["script"],
                "command": spec["command"],
                "outputs": outputs,
                "expected": expected,
                "status": state,
            }
        )

    summary_declared = declared["summary"]
    consistency = {
        "summary_present": summary_declared["present"],
        "genomes_processed_declared": summary_declared["genomes_processed"],
        "genomes_typed_observed": coverage["genomes_with_typed_mobile_in_inventory"]
        if coverage["genomes_with_typed_mobile_in_inventory"] is not None
        else coverage["genomes_with_typed_mobile"],
        "rows_declared": summary_declared["rows_total"],
        "rows_observed": coverage["rows_total"],
        "match": None,
    }
    if summary_declared["present"] and summary_declared["rows_total"] is not None:
        consistency["match"] = (
            summary_declared["genomes_processed"] == coverage["genomes_with_typed_mobile"]
            and summary_declared["rows_total"] == coverage["rows_total"]
        )

    return {
        "attribution": (
            "Atribución estática por scripts/mobile/22_run_full.sh, "
            "scripts/mobile/paths.py y los nombres de salida del repo. "
            "Los TSV no guardan metadatos del generador."
        ),
        "chains": [
            {
                "mode": "full",
                "orchestrator": "scripts/mobile/22_run_full.sh",
                "command": "bash scripts/mobile/22_run_full.sh all",
                "scripts": "19_prepare_full, 01_extract_candidates, 20_pfam_pool, "
                "21_crispr_pool, 05_classify, 06_write_outputs, 07_validate",
                "active": layout["mode"] == "full",
            },
            {
                "mode": "pilot",
                "orchestrator": "scripts/mobile/run_all.sh",
                "command": "bash scripts/mobile/run_all.sh",
                "scripts": "00_select_genomes, 01_extract_candidates, "
                "02_build_pfam_panel, 03_scan_pfam, 04_crispr_arrays, "
                "05_classify, 06_write_outputs, 07_validate",
                "active": layout["mode"] == "pilot",
            },
        ],
        "stages": stages,
        "prior_sources": [
            {
                "section": "mobile_element",
                "json": "data/mobile_elements.json",
                "script": "scripts/annotation/22_annotate_transposons.py",
                "flattened": "data/mobile_defense/features.tsv",
            },
            {
                "section": "defense / crispr / prophage",
                "json": "data/defense_systems.json",
                "script": "scripts/annotation/20_annotate_defense.py",
                "flattened": "data/mobile_defense/features.tsv",
            },
            {
                "section": "amr / virulence",
                "json": "data/defense_virulence.json",
                "script": "scripts/annotation/21_annotate_amr_virulence.py",
                "flattened": "data/mobile_defense/features.tsv",
            },
        ],
        "declared": declared,
        "consistency": consistency,
        "notes": [
            "La capa previa (features.tsv) no tiene contig ni hebra y es keyword-based; "
            "no sustituye al TSV tipado.",
            "No se encontró en el repo el script que aplana los tres JSON en "
            "data/mobile_defense/features.tsv; se audita el archivo existente.",
            "Los scripts 20/21/22 de scripts/annotation/ son los generadores de los "
            "JSON previos, no de la capa tipada.",
        ],
    }


def build_report(data_root: Path | str = DEFAULT_DATA) -> dict:
    """Construye el informe de cobertura sin escribir nada."""
    data_root = Path(data_root)
    warnings: list[str] = []
    if not data_root.is_dir():
        warnings.append(f"data root inexistente: {data_root}")

    layout = _layout(data_root)
    inventory = _load_inventory(data_root, layout)
    discovery = _discover(data_root, layout)
    prior_by_genome, prior_by_type = _load_prior(discovery)

    inventory_set = set(inventory["genomes"])
    universe = _dedupe(sorted(inventory_set | set(discovery["all_genomes"])))
    if inventory["source"] is None:
        warnings.append(
            "No hay full_genomes.tsv ni pilot_genomes.tsv con filas; "
            "no se puede medir contra un total esperado."
        )

    per_genome = [
        _inspect_genome(
            genome,
            data_root,
            layout,
            inventory,
            discovery,
            prior_by_genome,
        )
        for genome in universe
    ]
    by_id = {item["genome"]: item for item in per_genome}

    if inventory["expected_genomes"] is None:
        scope = universe
    else:
        scope = inventory["genomes"]

    typed_available = [genome for genome in scope if by_id[genome]["status"] == "available"]
    typed_missing = [genome for genome in scope if by_id[genome]["status"] != "available"]
    gff_missing = [
        genome for genome in scope if by_id[genome]["layers"]["gff"]["status"] != "available"
    ]
    spacer_missing = [
        genome
        for genome in scope
        if by_id[genome]["layers"]["crispr_scan"]["status"] != "available"
        and by_id[genome]["has_fasta"]
    ]
    fasta_missing = [genome for genome in scope if not by_id[genome]["has_fasta"]]
    crispr_not_evaluated = [
        genome
        for genome in scope
        if not by_id[genome]["layers"]["crispr_scan"]["evaluated"]
    ]
    prior_missing = [
        genome
        for genome in scope
        if by_id[genome]["layers"]["prior_features"]["status"] != "available"
    ]
    without_any_output = [
        genome
        for genome in scope
        if all(
            by_id[genome]["layers"][layer]["status"] != "available"
            for layer in ("typed_mobile", "gff", "crispr_scan", "prior_features")
        )
    ]
    not_in_inventory = sorted(set(discovery["all_genomes"]) - inventory_set)
    inventory_without_annotation = [
        genome for genome in scope if not by_id[genome]["inputs"]["annotation_gff"]
    ]

    typed_genomes = len(typed_available)
    rows_total = sum(item["typed_rows"] for item in per_genome)
    spacer_rows_total = sum(item["layers"]["crispr_scan"]["rows"] for item in per_genome)
    genomes_with_crispr_arrays = sum(
        1 for item in per_genome if item["typed_by_type"].get("crispr", 0) > 0
    )
    genomes_with_cas = sum(
        1 for item in per_genome if item["typed_by_type"].get("cas", 0) > 0
    )
    full_count = sum(1 for item in per_genome if item["coverage"] == "full")
    partial_count = sum(1 for item in per_genome if item["coverage"] == "partial")
    none_count = sum(1 for item in per_genome if item["coverage"] == "none")

    typed_by_type: Counter = Counter()
    prior_by_type_all: Counter = Counter(prior_by_type)
    typed_by_category: Counter = Counter()
    prior_by_category: Counter = Counter()
    category_genomes: dict[str, dict[str, set[str]]] = {
        category: {"typed": set(), "prior": set()} for category in CATEGORIES
    }
    for item in per_genome:
        for etype, count in item["typed_by_type"].items():
            typed_by_type[etype] += count
        for category, count in item["typed_by_category"].items():
            typed_by_category[category] += count
            if count:
                category_genomes[category]["typed"].add(item["genome"])
        for category, count in item["prior_by_category"].items():
            prior_by_category[category] += count
            if count:
                category_genomes[category]["prior"].add(item["genome"])

    category_report = {}
    for category in CATEGORIES:
        typed_set = category_genomes[category]["typed"]
        prior_set = category_genomes[category]["prior"]
        category_report[category] = {
            "types": list(CATEGORY_TYPES[category]),
            "typed_rows": int(typed_by_category.get(category, 0)),
            "typed_genomes": len(typed_set),
            "prior_rows": int(prior_by_category.get(category, 0)),
            "prior_genomes": len(prior_set),
            "genomes_with_any": len(typed_set | prior_set),
        }

    unknown_typed = Counter()
    for item in per_genome:
        for etype, count in item["unknown_types"].items():
            unknown_typed[etype] += count
    if unknown_typed:
        warnings.append(
            "Hay element_type sin mapear a categoría: "
            + ", ".join(sorted(unknown_typed))
        )

    for item in per_genome:
        if item["has_fasta_declared"] is not None and item["has_fasta_declared"] != item["has_fasta"]:
            warnings.append(
                f"{item['genome']}: inventario dice has_fasta="
                f"{'yes' if item['has_fasta_declared'] else 'no'} pero el archivo "
                f"{'existe' if item['has_fasta'] else 'no existe'}"
            )

    inventory_with_typed = len(typed_available) if inventory["expected_genomes"] is not None else None
    if inventory["expected_genomes"] is not None:
        genome_denominator = inventory["expected_genomes"]
    else:
        genome_denominator = len(universe)
    coverage = {
        "genomes_expected": inventory["expected_genomes"],
        "genomes_accounted": len(universe),
        "genomes_discovered": len(discovery["all_genomes"]),
        "genomes_with_typed_mobile": typed_genomes,
        "genomes_without_typed_mobile": len(typed_missing),
        "genomes_with_typed_mobile_in_inventory": inventory_with_typed,
        "typed_pct": _pct(
            inventory_with_typed
            if inventory["expected_genomes"] is not None
            else typed_genomes,
            genome_denominator,
        ),
        "genomes_coverage_full": full_count,
        "genomes_coverage_partial": partial_count,
        "genomes_coverage_none": none_count,
        "genomes_with_gff": sum(
            1 for item in per_genome if item["layers"]["gff"]["status"] == "available"
        ),
        "genomes_with_crispr_scan": sum(
            1 for item in per_genome if item["layers"]["crispr_scan"]["status"] == "available"
        ),
        "genomes_with_prior_features": sum(
            1 for item in per_genome if item["layers"]["prior_features"]["status"] == "available"
        ),
        "genomes_with_fasta": sum(1 for item in per_genome if item["has_fasta"]),
        "genomes_crispr_not_evaluated": len(crispr_not_evaluated),
        "rows_total": rows_total,
        "crispr_arrays_total": int(typed_by_type.get("crispr", 0)),
        "cas_rows_total": int(typed_by_type.get("cas", 0)),
        "spacer_rows_total": int(spacer_rows_total),
        "genomes_with_crispr_arrays": genomes_with_crispr_arrays,
        "genomes_with_cas": genomes_with_cas,
        "typed_by_type": _counter_dict(typed_by_type),
        "typed_by_category": _category_counts(typed_by_category),
        "prior_by_type": _counter_dict(prior_by_type_all),
        "prior_by_category": _category_counts(prior_by_category),
        "categories": category_report,
    }

    if not universe and inventory["expected_genomes"] is None:
        status = "no_data"
    elif inventory["expected_genomes"] is None:
        status = "no_inventory"
    elif typed_missing:
        status = "gaps"
    else:
        status = "covered"

    typed_for_summary = (
        inventory_with_typed
        if inventory_with_typed is not None
        else typed_genomes
    )
    declared = _declared_artifacts(data_root, layout)
    if declared["summary"]["present"] and declared["summary"].get("rows_total") is not None:
        if declared["summary"]["rows_total"] != rows_total:
            warnings.append(
                "summary declara rows_total="
                f"{declared['summary']['rows_total']} y el recuento real es {rows_total}"
            )
    if declared["summary"]["present"] and declared["summary"].get("genomes_processed") is not None:
        if declared["summary"]["genomes_processed"] != typed_for_summary:
            warnings.append(
                "summary declara genomes_processed="
                f"{declared['summary']['genomes_processed']} y hay {typed_for_summary} TSV con filas"
            )

    if status == "gaps":
        warnings.append(
            f"{len(typed_missing)} genoma(s) del inventario sin TSV tipado."
        )
    if crispr_not_evaluated:
        if len(fasta_missing) == len(crispr_not_evaluated):
            detail = "por FASTA ausente"
        else:
            detail = "por FASTA ausente o salida de espaciadores incompleta"
        warnings.append(
            f"{len(crispr_not_evaluated)} genoma(s) sin CRISPR evaluable "
            f"({detail}); no significa ausencia de arrays."
        )
    if not_in_inventory:
        warnings.append(
            f"{len(not_in_inventory)} genoma(s) con artefactos fuera del inventario."
        )

    sources = {
        "mobile_typing_full_tsv": len(discovery["typed_full"]),
        "mobile_typing_pilot_tsv": len(discovery["typed_pilot"]),
        "mobile_typing_full_gff": len(discovery["gff_full"]),
        "mobile_typing_full_spacers": len(discovery["spacer_full"]),
        "mobile_defense_features_genomes": len(discovery["prior_genomes"]),
        "annotations_ab_initio_genomes": len(discovery["annotation_genomes"]),
        "raw_fasta_genomes": len(discovery["fasta_genomes"]),
    }

    report = {
        "schema": SCHEMA,
        "generated_by": "tools/mobile_coverage.py",
        "read_only": True,
        "data_root": str(data_root),
        "mode": layout["mode"],
        "status": status,
        "inventory": {
            "source": inventory["source"],
            "expected_genomes": inventory["expected_genomes"],
            "rows": inventory["rows"],
            "columns": inventory["columns"],
            "empty_candidates": inventory["empty_candidates"],
        },
        "discovery": {
            "sources": sources,
            "discovered_genomes": len(discovery["all_genomes"]),
            "not_in_inventory": not_in_inventory,
            "inventory_without_annotation": inventory_without_annotation,
        },
        "coverage": coverage,
        "pipeline": _pipeline_report(
            data_root, layout, inventory, discovery, coverage, declared
        ),
        "gaps": {
            "scope": "inventory" if inventory["expected_genomes"] is not None else "universe",
            "typed_missing": typed_missing,
            "gff_missing": gff_missing,
            "spacer_missing": spacer_missing,
            "fasta_missing": fasta_missing,
            "crispr_not_evaluated": crispr_not_evaluated,
            "prior_missing": prior_missing,
            "without_any_output": without_any_output,
            "to_reach_inventory": len(typed_missing)
            if inventory["expected_genomes"] is not None
            else None,
        },
        "warnings": warnings,
        "per_genome": per_genome,
    }
    return report


# ---------------------------------------------------------------------------
# Render Markdown
# ---------------------------------------------------------------------------
def _fmt_int(value: int | None) -> str:
    if value is None:
        return "—"
    return f"{value:,}".replace(",", ".")


def _fmt_float(value: float | None) -> str:
    return "—" if value is None else f"{value:.2f}".replace(".", ",")


def _md(text: object) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def _code_list(items: Iterable[str], limit: int = 100) -> str:
    values = list(items)
    if not values:
        return "—"
    shown = ", ".join(f"`{value}`" for value in values[:limit])
    if len(values) > limit:
        shown += f" … (+{len(values) - limit} más; ver JSON)"
    return shown


def _status_reading(status: str, report: dict) -> str:
    coverage = report["coverage"]
    if status == "no_data":
        return "No se encontraron datos móviles/defensa en el árbol auditado."
    if status == "no_inventory":
        return (
            "No hay lista de genomas esperados; la cobertura se mide sobre los "
            f"{coverage['genomes_discovered']} genomas descubiertos."
        )
    if status == "gaps":
        return (
            "El inventario no está cubierto: faltan TSV tipados para "
            f"{len(report['gaps']['typed_missing'])} genoma(s)."
        )
    return (
        "El inventario está cubierto por TSV tipados. Las brechas restantes son "
        "de evaluabilidad (FASTA/CRISPR) o de métodos externos, no de filas."
    )


def render_markdown(report: dict) -> str:
    """Renderiza el informe en Markdown determinista (sin fecha de ejecución)."""
    lines: list[str] = []
    add = lines.append

    add("# Cobertura móvil, defensa y AMR — estado real en disco")
    add("")
    add(
        "Auditoría de **solo lectura** generada por `tools/mobile_coverage.py`. "
        "No ejecuta HMMER, DefenseFinder, CRISPRCasFinder, AMRFinderPlus, "
        "ISEScan ni bases externas, y no crea ni completa anotaciones."
    )
    add("")
    add("## Veredicto")
    add("")
    add(f"**Estado global:** `{report['status']}` — {_status_reading(report['status'], report)}")
    add("")
    cov = report["coverage"]
    inv = report["inventory"]
    add("| métrica | valor |")
    add("|---|---:|")
    inventory_label = f"`{inv['source']}`" if inv["source"] else "no encontrado"
    add(f"| inventario | {inventory_label} |")
    add(f"| genomas esperados (inventario) | {_fmt_int(inv['expected_genomes'])} |")
    add(f"| genomas descubiertos por artefactos | {_fmt_int(cov['genomes_discovered'])} |")
    add(
        "| genomas con TSV tipado | "
        f"{_fmt_int(cov['genomes_with_typed_mobile'])} / {_fmt_int(inv['expected_genomes'])} "
        f"({_fmt_float(cov['typed_pct'])} %) |"
    )
    add(f"| genomas con GFF3 | {_fmt_int(cov['genomes_with_gff'])} |")
    add(f"| genomas con capa previa | {_fmt_int(cov['genomes_with_prior_features'])} |")
    add(
        "| FASTA disponible | "
        f"{_fmt_int(cov['genomes_with_fasta'])} / {_fmt_int(cov['genomes_accounted'])} |"
    )
    add(
        "| CRISPR evaluable (con FASTA) | "
        f"{_fmt_int(cov['genomes_with_crispr_scan'])}; "
        f"no evaluados: {_fmt_int(cov['genomes_crispr_not_evaluated'])} |"
    )
    add(
        "| arrays CRISPR / espaciadores | "
        f"{_fmt_int(cov['crispr_arrays_total'])} arrays en "
        f"{_fmt_int(cov['genomes_with_crispr_arrays'])} genomas / "
        f"{_fmt_int(cov['spacer_rows_total'])} espaciadores |"
    )
    add(
        "| filas `cas` / genomas con `cas` | "
        f"{_fmt_int(cov['cas_rows_total'])} / {_fmt_int(cov['genomes_with_cas'])} |"
    )
    add(f"| filas tipadas | {_fmt_int(cov['rows_total'])} |")
    add("")
    add(
        "Cobertura por genoma: "
        f"completa **{_fmt_int(cov['genomes_coverage_full'])}**, "
        f"parcial **{_fmt_int(cov['genomes_coverage_partial'])}**, "
        f"sin datos tipados **{_fmt_int(cov['genomes_coverage_none'])}**."
    )
    add("")

    add("## Inventario y fuentes")
    add("")
    add(
        "El total esperado **no está fijado en el código**: se lee de "
        f"`{inv['source']}`. Si ese archivo no existe, la cobertura se mide "
        "sobre los genomas descubiertos."
    )
    add("")
    add("| fuente | genomas/archivos |")
    add("|---|---:|")
    for key, value in report["discovery"]["sources"].items():
        add(f"| `{key}` | {_fmt_int(value)} |")
    add("")
    if report["discovery"]["not_in_inventory"]:
        add(
            "**Descubiertos fuera del inventario:** "
            + _code_list(report["discovery"]["not_in_inventory"])
        )
        add("")

    add("## Tipos agregados")
    add("")
    add("Mapeo aplicado (no añade tipos, solo agrupa los existentes):")
    add("")
    add("| categoría | `element_type` incluidos | filas tipadas | genomas tipados | filas previas | genomas previos |")
    add("|---|---|---:|---:|---:|---:|")
    for category in CATEGORIES:
        row = cov["categories"][category]
        add(
            f"| `{category}` | {', '.join('`' + t + '`' for t in row['types'])} | "
            f"{_fmt_int(row['typed_rows'])} | {_fmt_int(row['typed_genomes'])} | "
            f"{_fmt_int(row['prior_rows'])} | {_fmt_int(row['prior_genomes'])} |"
        )
    add("")
    add(
        "Tipos crudos en el TSV tipado: "
        + ", ".join(
            f"`{key}`: {_fmt_int(value)}"
            for key, value in cov["typed_by_type"].items()
        )
    )
    add("")
    add(
        "La capa previa (features.tsv) es keyword-based y no tiene contig ni "
        "hebra; sus recuentos se muestran aparte y **no** sustituyen al TSV."
    )
    add("")

    add("## Pipeline que generó los datos")
    add("")
    add(report["pipeline"]["attribution"])
    add("")
    add("| cadena | orquestador | comando | activa |")
    add("|---|---|---|---|")
    for chain in report["pipeline"]["chains"]:
        add(
            f"| {chain['mode']} | `{chain['orchestrator']}` | "
            f"`{chain['command']}` | {'sí' if chain['active'] else 'no'} |"
        )
    add("")
    add("| etapa | script | salida | encontrado | esperado | estado |")
    add("|---|---|---|---:|---:|---|")
    for stage in report["pipeline"]["stages"]:
        output = "<br>".join(
            f"`{item['path']}` ({_fmt_int(item['found'])})" for item in stage["outputs"]
        )
        found = "/".join(_fmt_int(item["found"]) for item in stage["outputs"])
        add(
            f"| {stage['stage']} | `{stage['script']}` | {output} | "
            f"{found} | {_fmt_int(stage['expected'])} | `{stage['status']}` |"
        )
    add("")
    consistency = report["pipeline"]["consistency"]
    add("**Contraste con artefactos declarados**")
    add("")
    add("| artefacto | declarado | observado | coincide |")
    add("|---|---:|---:|---|")
    add(
        f"| `{report['pipeline']['declared']['summary']['path']}` (genomas) | "
        f"{_fmt_int(consistency['genomes_processed_declared'])} | "
        f"{_fmt_int(consistency['genomes_typed_observed'])} | "
        f"{'sí' if consistency['match'] else ('—' if consistency['match'] is None else 'no')} |"
    )
    add(
        f"| `{report['pipeline']['declared']['summary']['path']}` (filas) | "
        f"{_fmt_int(consistency['rows_declared'])} | "
        f"{_fmt_int(consistency['rows_observed'])} | "
        f"{'sí' if consistency['match'] else ('—' if consistency['match'] is None else 'no')} |"
    )
    validation = report["pipeline"]["declared"]["validation"]
    add(
        f"| `{validation['path']}` | status={validation['status']} | "
        f"errors={validation['errors']} | — |"
    )
    add("")
    add("Fuentes previas (solo lectura; no son la capa tipada):")
    add("")
    add("| sección | script | JSON | aplanado |")
    add("|---|---|---|---|")
    for source in report["pipeline"]["prior_sources"]:
        add(
            f"| {source['section']} | `{source['script']}` | `{source['json']}` | "
            f"`{source['flattened']}` |"
        )
    add("")
    add(
        "No se encontró el script que aplana esos JSON en "
        "`data/mobile_defense/features.tsv`; se audita el archivo tal como está."
    )
    add("")

    add("## Brechas")
    add("")
    gaps = report["gaps"]
    add(f"Alcance de las listas: **{gaps['scope']}**.")
    add("")
    add("| brecha | n | detalle |")
    add("|---|---:|---|")
    add(
        f"| TSV tipado ausente | {len(gaps['typed_missing'])} | "
        f"{_code_list(gaps['typed_missing'])} |"
    )
    add(
        f"| GFF3 ausente | {len(gaps['gff_missing'])} | "
        f"{_code_list(gaps['gff_missing'])} |"
    )
    add(
        f"| Espaciadores ausentes con FASTA presente | {len(gaps['spacer_missing'])} | "
        f"{_code_list(gaps['spacer_missing'])} |"
    )
    add(
        f"| FASTA ausente | {len(gaps['fasta_missing'])} | "
        f"{_code_list(gaps['fasta_missing'])} |"
    )
    add(
        f"| CRISPR no evaluado | {len(gaps['crispr_not_evaluated'])} | "
        f"{_code_list(gaps['crispr_not_evaluated'])} |"
    )
    add(
        f"| Sin capa previa | {len(gaps['prior_missing'])} | "
        f"{_code_list(gaps['prior_missing'])} |"
    )
    add(
        f"| Sin ninguna salida (TSV/GFF/CRISPR/previo) | {len(gaps['without_any_output'])} | "
        f"{_code_list(gaps['without_any_output'])} |"
    )
    add("")
    if gaps["to_reach_inventory"] is not None:
        add(
            "Para cubrir el inventario faltan "
            f"**{gaps['to_reach_inventory']}** TSV tipados."
        )
    else:
        add("No hay inventario contra el cual calcular un faltante.")
    add("")
    if gaps["fasta_missing"]:
        add(
            "> Los genomas sin FASTA quedan `not_available` para CRISPR de novo: "
            "**no** se interpretan como “sin arrays”. El resto de sus tipos sigue "
            "disponible si existe TSV."
        )
        add("")

    add("## Cobertura por filo")
    add("")
    add("| filo | genomas | con TSV | completa | parcial | sin datos | filas | FASTA |")
    add("|---|---:|---:|---:|---:|---:|---:|---:|")
    phylum_stats: dict[str, dict] = defaultdict(
        lambda: {
            "genomes": 0,
            "typed": 0,
            "full": 0,
            "partial": 0,
            "none": 0,
            "rows": 0,
            "fasta": 0,
        }
    )
    for item in report["per_genome"]:
        stats = phylum_stats[item["phylum"] or "unknown"]
        stats["genomes"] += 1
        stats["typed"] += 1 if item["status"] == "available" else 0
        stats[item["coverage"]] += 1
        stats["rows"] += item["typed_rows"]
        stats["fasta"] += 1 if item["has_fasta"] else 0
    for phylum in sorted(
        phylum_stats, key=lambda name: (-phylum_stats[name]["genomes"], name)
    ):
        stats = phylum_stats[phylum]
        add(
            f"| {_md(phylum)} | {_fmt_int(stats['genomes'])} | {_fmt_int(stats['typed'])} | "
            f"{_fmt_int(stats['full'])} | {_fmt_int(stats['partial'])} | "
            f"{_fmt_int(stats['none'])} | {_fmt_int(stats['rows'])} | "
            f"{_fmt_int(stats['fasta'])} |"
        )
    add("")

    add("## Detalle por genoma")
    add("")
    add(
        "`M`=mobile_element, `P`=prophage, `D`=defense, `C`=crispr (incluye cas), "
        "`A`=amr, `V`=virulence. `no disponible` lista las capas `not_available`."
    )
    add("")
    add("<details>")
    add(f"<summary>Ver los {len(report['per_genome'])} genomas</summary>")
    add("")
    add("| genoma | filo | estado | cobertura | filas | M | P | D | C | A | V | FASTA | no disponible |")
    add("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|")
    for item in report["per_genome"]:
        counts = item["typed_by_category"]
        if item["status"] == "available":
            values = [
                _fmt_int(counts.get(category, 0))
                for category in ("mobile_element", "prophage", "defense", "crispr", "amr", "virulence")
            ]
        else:
            values = ["—"] * 6
        missing = ", ".join(item["not_available"]) or "—"
        add(
            f"| `{item['genome']}` | {_md(item['phylum'])} | `{item['status']}` | "
            f"`{item['coverage']}` | {_fmt_int(item['typed_rows'])} | "
            + " | ".join(values)
            + f" | {'sí' if item['has_fasta'] else 'no'} | {_md(missing)} |"
        )
    add("")
    add("</details>")
    add("")

    if report["warnings"]:
        add("## Advertencias")
        add("")
        for warning in report["warnings"]:
            add(f"- {warning}")
        add("")

    add("## Método y límites")
    add("")
    add("- Solo lectura: el informe se construye leyendo TSV/GFF/JSON/FASTA existentes.")
    add(
        "- `not_available` significa que **no hay dato en disco**, no que la "
        "biología esté ausente. En particular, CRISPR sin FASTA es “no evaluado”."
    )
    add(
        "- La atribución de cada artefacto a un script es estática "
        "(`22_run_full.sh`, `paths.py`, nombres de salida); no se ejecutó nada."
    )
    add(
        "- Los recuentos de la capa previa y del TSV tipado no son equivalentes: "
        "la primera es keyword-based y sin coordenadas por contig."
    )
    add("- Para reproducir: `python tools/mobile_coverage.py`.")
    add("")
    return "\n".join(lines)


def write_reports(report: dict, json_path: Path | str, md_path: Path | str) -> None:
    """Escribe JSON y Markdown; crea los directorios si hacen falta."""
    json_path = Path(json_path)
    md_path = Path(md_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    md_path.write_text(render_markdown(report), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audita la cobertura real de móviles/defensa/CRISPR/AMR/virulencia "
            "sin ejecutar herramientas externas ni crear anotaciones."
        )
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=DEFAULT_DATA,
        help="raíz de datos (por defecto data/)",
    )
    parser.add_argument(
        "--json",
        type=Path,
        default=DEFAULT_JSON,
        help="ruta del informe JSON",
    )
    parser.add_argument(
        "--md",
        type=Path,
        default=DEFAULT_MD,
        help="ruta del informe Markdown",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="no imprimir el resumen por stdout",
    )
    args = parser.parse_args(argv)

    report = build_report(args.data)
    write_reports(report, args.json, args.md)
    if not args.quiet:
        cov = report["coverage"]
        print(f"estado={report['status']} modo={report['mode']}")
        print(
            "genomas="
            f"{cov['genomes_discovered']} esperados={cov['genomes_expected']} "
            f"con_tsv={cov['genomes_with_typed_mobile']} "
            f"sin_tsv={cov['genomes_without_typed_mobile']} "
            f"con_fasta={cov['genomes_with_fasta']} "
            f"crispr_no_evaluado={cov['genomes_crispr_not_evaluated']} "
            f"filas={cov['rows_total']}"
        )
        print(f"json={args.json}")
        print(f"md={args.md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
