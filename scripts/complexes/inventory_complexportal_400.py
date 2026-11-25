#!/usr/bin/env python3
"""Inventario local de catálogos ComplexPortal para los genomas del proyecto.

Cruza, sin red y sin descargar catálogos:

- ``data/regulation/*/genes.tsv``          (qué genomas existen)
- ``data/taxonomy/full/*.json``            (especie y taxid NCBI)
- ``data/taxonomy/*.json``                 (fallback de los 10 pilotos)
- ``data/complexes/gold_<taxid>.tsv``      (catálogo ya mapeado a genes)
- ``data/complexes/raw/complexportal_<taxid>.tsv`` (TSV crudo de ComplexPortal)

La referencia K-12 no es un ensamblado GCF y no tiene JSON en
``data/taxonomy``. Su único vínculo local con ComplexPortal es el TSV crudo
``raw/complexportal_83333.tsv``; el taxid se lee de la cabecera de ese TSV
(columna ``Taxonomy identifier``), nunca se escribe a mano. Si el TSV no está,
el taxid queda ``not_available``.

Salidas por defecto:

- ``reports/complexportal_inventory.tsv``
- ``reports/complexportal_inventory.md``

Columnas: genome, species, taxid, local_catalog, status, reason.

Estados
-------
ready              catálogo gold local (mapeado) y taxid disponible.
raw_only           solo TSV crudo; falta mapearlo a genes.
catalog_missing    taxid disponible, pero no hay catálogo local.
taxid_missing      hay especie/JSON, pero no taxid; no se puede resolver.
taxonomy_missing   no hay especie ni taxid ni JSON de taxonomía.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NOT_AVAILABLE = "not_available"

STATUS_READY = "ready"
STATUS_RAW_ONLY = "raw_only"
STATUS_CATALOG_MISSING = "catalog_missing"
STATUS_TAXID_MISSING = "taxid_missing"
STATUS_TAXONOMY_MISSING = "taxonomy_missing"

STATUS_ORDER = (
    STATUS_READY,
    STATUS_RAW_ONLY,
    STATUS_CATALOG_MISSING,
    STATUS_TAXID_MISSING,
    STATUS_TAXONOMY_MISSING,
)
STATUS_MEANING = {
    STATUS_READY: "catálogo gold local + taxid",
    STATUS_RAW_ONLY: "TSV crudo local, sin mapear a genes",
    STATUS_CATALOG_MISSING: "taxid disponible, sin catálogo local",
    STATUS_TAXID_MISSING: "sin taxid; no se puede resolver catálogo",
    STATUS_TAXONOMY_MISSING: "sin especie, taxid ni JSON de taxonomía",
}

# El genoma de referencia K-12 no es un GCF del panel y no tiene JSON de
# taxonomía. El alias es explícito para no confundirlo con una especie del
# panel; el taxid se obtiene del propio catálogo, no de una constante escrita
# a mano. Si el TSV desaparece, el taxid vuelve a not_available.
REFERENCE_CATALOGS: dict[str, tuple[Path, Path]] = {
    "ecoli_k12": (
        Path("data/complexes/raw/complexportal_83333.tsv"),
        Path("data/complexes/gold_k12.tsv"),
    ),
}

TAXONOMY_DIRS = ("data/taxonomy/full", "data/taxonomy")
TSV_COLUMNS = ("genome", "species", "taxid", "local_catalog", "status", "reason")
_EMPTY_VALUES = {"", "unknown", "none", "null", "not_available", "na", "n/a"}


@dataclass(frozen=True)
class Taxonomy:
    species: str = NOT_AVAILABLE
    taxid: int | None = None
    json_present: bool = False


def _positive_int(value: object) -> int | None:
    """Devuelve un taxid entero positivo o None. No convierte floats."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value if value > 0 else None
    if isinstance(value, str):
        text = value.strip()
        if text.isdigit():
            number = int(text)
            return number if number > 0 else None
    return None


def _text(value: object) -> str:
    if value is None:
        return NOT_AVAILABLE
    text = str(value).strip()
    if text.lower() in _EMPTY_VALUES:
        return NOT_AVAILABLE
    return text


def _species_from_json(data: dict) -> str:
    for key in ("species", "species_ncbi", "species_gtdb"):
        species = _text(data.get(key))
        if species != NOT_AVAILABLE:
            return species
    evidence = data.get("evidence")
    if isinstance(evidence, dict):
        gtdb = evidence.get("gtdb")
        if isinstance(gtdb, dict):
            species = _text(gtdb.get("species"))
            if species != NOT_AVAILABLE:
                return species
    lineage = _text(data.get("gtdb_lineage"))
    if lineage != NOT_AVAILABLE:
        last = lineage.split(";")[-1].strip()
        if last.startswith("s__") and len(last) > 3:
            return last[3:].replace("_", " ")
    return NOT_AVAILABLE


def _taxid_from_json(data: dict) -> int | None:
    """Taxid del ensamblado; nunca inventa ni usa species_taxid como sustituto."""
    candidates: list[object] = [data.get("ncbi_taxid"), data.get("taxid")]
    evidence = data.get("evidence")
    if isinstance(evidence, dict):
        ncbi = evidence.get("ncbi")
        if isinstance(ncbi, dict):
            candidates.append(ncbi.get("taxid"))
    for value in candidates:
        taxid = _positive_int(value)
        if taxid is not None:
            return taxid
    return None


def _species_from_genbank(root: Path, genome: str) -> str:
    """Fallback local para genomas sin JSON: lee ``ORGANISM`` del GBFF fuente."""
    meta_path = root / "data" / "regulation" / genome / "genes_meta.json"
    if not meta_path.is_file():
        return NOT_AVAILABLE
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return NOT_AVAILABLE
    source = meta.get("source")
    gbff: object = None
    if isinstance(source, dict):
        gbff = source.get("gbff") or source.get("path")
    if not gbff:
        seq_source = meta.get("seq_source")
        if isinstance(seq_source, dict):
            gbff = seq_source.get("path")
    if not gbff:
        return NOT_AVAILABLE
    path = Path(str(gbff))
    if not path.is_absolute():
        path = root / path
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            for _, line in zip(range(200), fh):
                match = re.match(r"\s+ORGANISM\s+(.+?)\s*$", line)
                if match:
                    return _text(match.group(1))
    except OSError:
        return NOT_AVAILABLE
    return NOT_AVAILABLE


def resolve_taxonomy(root: Path, genome: str) -> Taxonomy:
    """Especie/taxid desde JSON de taxonomía; especie desde GBFF si falta JSON."""
    for relative in TAXONOMY_DIRS:
        path = root / relative / f"{genome}.json"
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return Taxonomy(_species_from_genbank(root, genome), None, True)
        if not isinstance(data, dict):
            return Taxonomy(_species_from_genbank(root, genome), None, True)
        species = _species_from_json(data)
        if species == NOT_AVAILABLE:
            species = _species_from_genbank(root, genome)
        return Taxonomy(species, _taxid_from_json(data), True)
    return Taxonomy(_species_from_genbank(root, genome), None, False)


def _catalog_taxid(path: Path) -> int | None:
    """Lee ``Taxonomy identifier`` del TSV; fallback al número del nombre."""
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            try:
                index = header.index("Taxonomy identifier")
            except ValueError:
                index = 3
            row = fh.readline().rstrip("\r\n").split("\t")
            if index < len(row):
                taxid = _positive_int(row[index])
                if taxid is not None:
                    return taxid
    except OSError:
        return None
    match = re.search(r"complexportal_(\d+)\.tsv$", path.name)
    return _positive_int(match.group(1)) if match else None


def _reference_catalog_paths(root: Path, genome: str) -> tuple[Path, Path] | None:
    relative = REFERENCE_CATALOGS.get(genome)
    if relative is None:
        return None
    raw, gold = relative
    return root / raw, root / gold


def find_local_catalog(root: Path, genome: str, taxid: int | None) -> tuple[Path, str] | None:
    """Devuelve (ruta, ``gold``|``raw``) o None. Prefiere el catálogo mapeado."""
    paths = _reference_catalog_paths(root, genome)
    if paths is not None:
        raw, gold = paths
        reference_taxid = _catalog_taxid(raw) if raw.is_file() else None
        if reference_taxid is not None and (taxid is None or taxid == reference_taxid):
            if gold.is_file():
                return gold, "gold"
            if raw.is_file():
                return raw, "raw"
    if taxid is not None:
        gold = root / "data" / "complexes" / f"gold_{taxid}.tsv"
        raw = root / "data" / "complexes" / "raw" / f"complexportal_{taxid}.tsv"
        if gold.is_file():
            return gold, "gold"
        if raw.is_file():
            return raw, "raw"
    return None


def _mapping_source(root: Path, genome: str) -> str:
    """Dónde hay un mapeo UniProt→gen para este genoma, si existe."""
    functional = root / "data" / "annotations_ab_initio" / genome / "functional.tsv"
    if functional.is_file():
        try:
            with functional.open(encoding="utf-8", errors="replace") as fh:
                first = fh.readline()
                second = fh.readline()
            # functional.tsv no tiene cabecera en este proyecto: la primera
            # línea ya es un hit. Se acepta también una cabecera + datos.
            if second.strip() or "sp|" in first or "tr|" in first:
                return "functional.tsv"
        except OSError:
            pass
    genes = root / "data" / "regulation" / genome / "genes.tsv"
    if genes.is_file():
        try:
            with genes.open(encoding="utf-8", errors="replace", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                if "acc" in (reader.fieldnames or []):
                    for index, row in enumerate(reader):
                        if (row.get("acc") or "").strip():
                            return "genes.tsv acc"
                        if index >= 1000:
                            break
        except (OSError, csv.Error):
            pass
    return NOT_AVAILABLE


def _mapping_label(source: str) -> str:
    if source == "functional.tsv":
        return "functional.tsv (UniProt→gene_id)"
    if source == "genes.tsv acc":
        return "genes.tsv acc (UniProt→locus_tag/gene_id)"
    return "no disponible"


def _reason(
    *,
    status: str,
    taxid: int | None,
    catalog_kind: str | None,
    species: str,
    mapping: str,
    taxid_from_catalog: bool,
) -> str:
    if status == STATUS_READY:
        text = f"catálogo {catalog_kind} local para taxid {taxid}"
    elif status == STATUS_RAW_ONLY:
        text = f"catálogo crudo local para taxid {taxid} sin gold mapeado"
    elif status == STATUS_CATALOG_MISSING:
        text = f"sin catálogo ComplexPortal local para taxid {taxid}"
    elif status == STATUS_TAXID_MISSING:
        text = "sin taxid en los metadatos de taxonomía"
    else:
        text = "sin especie, taxid ni JSON de taxonomía"
    if taxid_from_catalog:
        text += "; taxid leído del catálogo ComplexPortal local"
    if species == NOT_AVAILABLE:
        text += "; especie not_available"
    text += f"; mapeo UniProt→locus: {_mapping_label(mapping)}"
    return text


def discover_genomes(root: Path) -> list[str]:
    """Genomas con ``genes.tsv``, ordenados; ignora directorios ocultos."""
    regulation = root / "data" / "regulation"
    if not regulation.is_dir():
        return []
    genomes = {
        path.parent.name
        for path in regulation.glob("*/genes.tsv")
        if path.is_file() and not path.parent.name.startswith(".")
    }
    return sorted(genomes)


def build_inventory(root: Path = ROOT) -> list[dict[str, str]]:
    """Construye las filas del inventario. No escribe ni descarga nada."""
    rows: list[dict[str, str]] = []
    for genome in discover_genomes(root):
        taxonomy = resolve_taxonomy(root, genome)
        taxid = taxonomy.taxid
        taxid_from_catalog = False
        if taxid is None:
            paths = _reference_catalog_paths(root, genome)
            if paths is not None and paths[0].is_file():
                taxid = _catalog_taxid(paths[0])
                taxid_from_catalog = taxid is not None
        catalog = find_local_catalog(root, genome, taxid)
        catalog_path, catalog_kind = catalog if catalog else (None, None)
        mapping = _mapping_source(root, genome)

        if taxid is None:
            status = (
                STATUS_TAXID_MISSING
                if taxonomy.json_present or taxonomy.species != NOT_AVAILABLE
                else STATUS_TAXONOMY_MISSING
            )
        elif catalog_path is None:
            status = STATUS_CATALOG_MISSING
        elif catalog_kind == "raw":
            status = STATUS_RAW_ONLY
        else:
            status = STATUS_READY

        rows.append({
            "genome": genome,
            "species": taxonomy.species,
            "taxid": str(taxid) if taxid is not None else NOT_AVAILABLE,
            "local_catalog": (
                catalog_path.relative_to(root).as_posix()
                if catalog_path is not None else NOT_AVAILABLE
            ),
            "status": status,
            "reason": _reason(
                status=status,
                taxid=taxid,
                catalog_kind=catalog_kind,
                species=taxonomy.species,
                mapping=mapping,
                taxid_from_catalog=taxid_from_catalog,
            ),
        })
    return rows


def _md(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def render_markdown(rows: list[dict[str, str]], generated_at: str) -> str:
    panel = [row for row in rows if row["genome"].startswith("GCF_")]
    references = [row for row in rows if not row["genome"].startswith("GCF_")]
    with_taxid = [row for row in rows if row["taxid"] != NOT_AVAILABLE]
    panel_with_taxid = [row for row in panel if row["taxid"] != NOT_AVAILABLE]
    with_catalog = [row for row in rows if row["local_catalog"] != NOT_AVAILABLE]
    panel_with_catalog = [
        row for row in panel if row["local_catalog"] != NOT_AVAILABLE
    ]
    counts = {
        status: sum(1 for row in rows if row["status"] == status)
        for status in STATUS_ORDER
    }

    lines = [
        "# Inventario de catálogos ComplexPortal",
        "",
        f"Generado: `{generated_at}`",
        "",
        "Alcance: genomas con `data/regulation/*/genes.tsv`; especie/taxid desde "
        "`data/taxonomy/full/*.json` (fallback `data/taxonomy/*.json` y GBFF) y "
        "catálogos locales en `data/complexes/`. Este informe no descarga nada "
        "ni consulta la red.",
        "",
        "## Resumen",
        "",
        "| métrica | valor |",
        "|---|---:|",
        f"| genomas | {len(rows)} |",
        f"| panel GCF | {len(panel)} |",
        f"| referencias no-GCF | {len(references)} |",
        f"| con especie | {sum(row['species'] != NOT_AVAILABLE for row in rows)} |",
        f"| con taxid | {len(with_taxid)} |",
        f"| panel con taxid | {len(panel_with_taxid)} |",
        f"| con catálogo local | {len(with_catalog)} |",
        f"| panel con catálogo local | {len(panel_with_catalog)} |",
        f"| catálogo gold listo | {counts[STATUS_READY]} |",
        f"| TSV crudo sin mapear | {counts[STATUS_RAW_ONLY]} |",
        f"| sin catálogo local | {counts[STATUS_CATALOG_MISSING]} |",
        f"| sin taxid | {counts[STATUS_TAXID_MISSING]} |",
        f"| sin taxonomía | {counts[STATUS_TAXONOMY_MISSING]} |",
        "",
        "`not_available` significa que el dato no está en las fuentes locales; "
        "no significa que el taxid o el catálogo no existan en ComplexPortal.",
        "",
        "## Estados",
        "",
        "| status | n | significado |",
        "|---|---:|---|",
    ]
    for status in STATUS_ORDER:
        lines.append(
            f"| `{status}` | {counts[status]} | {STATUS_MEANING[status]} |"
        )

    lines.extend([
        "",
        "## Catálogos locales",
        "",
    ])
    if with_catalog:
        lines.extend([
            "| genome | species | taxid | local_catalog | status |",
            "|---|---|---|---|---|",
        ])
        for row in with_catalog:
            lines.append(
                f"| `{_md(row['genome'])}` | {_md(row['species'])} | "
                f"{_md(row['taxid'])} | `{_md(row['local_catalog'])}` | "
                f"`{_md(row['status'])}` |"
            )
    else:
        lines.append("No hay ningún catálogo ComplexPortal local.")

    lines.extend([
        "",
        "## Inventario completo",
        "",
        "| genome | species | taxid | local_catalog | status | reason |",
        "|---|---|---|---|---|---|",
    ])
    for row in rows:
        lines.append(
            f"| `{_md(row['genome'])}` | {_md(row['species'])} | "
            f"{_md(row['taxid'])} | `{_md(row['local_catalog'])}` | "
            f"`{_md(row['status'])}` | {_md(row['reason'])} |"
        )

    lines.extend([
        "",
        "## Regenerar",
        "",
        "```bash",
        "python scripts/complexes/inventory_complexportal_400.py",
        "```",
        "",
        "## Lectura",
        "",
        "- Un `catalog_missing` solo dice que no hay catálogo local. No prueba "
        "que ComplexPortal tenga o deje de tener una entrada para ese taxid; "
        "eso se comprueba contra el listado remoto (ver "
        "`docs/COMPLEXPORTAL_400.md`).",
        "- Los `model_<genoma>.complexes.tsv` y `transfer_<genoma>.tsv` no son "
        "catálogos ComplexPortal y no se cuentan como tales.",
        "- `raw_only` todavía no sirve para anotar: falta el mapeo "
        "UniProt→gen y, en su caso, el filtro de complejos incompletos.",
        "",
    ])
    return "\n".join(lines)


def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Inventario local de catálogos ComplexPortal (sin red)."
    )
    parser.add_argument("--root", type=Path, default=ROOT,
                        help="raíz del proyecto (default: la del script)")
    parser.add_argument("--out-tsv", type=Path, default=None,
                        help="default: <root>/reports/complexportal_inventory.tsv")
    parser.add_argument("--out-md", type=Path, default=None,
                        help="default: <root>/reports/complexportal_inventory.md")
    parser.add_argument("--generated-at", default=None,
                        help="timestamp fijo (para reproducibilidad/tests)")
    args = parser.parse_args(argv)

    root = args.root.resolve()
    rows = build_inventory(root)
    out_tsv = args.out_tsv or (root / "reports" / "complexportal_inventory.tsv")
    out_md = args.out_md or (root / "reports" / "complexportal_inventory.md")
    generated_at = args.generated_at or datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )

    write_tsv(rows, out_tsv)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(render_markdown(rows, generated_at), encoding="utf-8")

    panel = [row for row in rows if row["genome"].startswith("GCF_")]
    references = [row for row in rows if not row["genome"].startswith("GCF_")]
    with_taxid = [row for row in rows if row["taxid"] != NOT_AVAILABLE]
    panel_with_taxid = [row for row in panel if row["taxid"] != NOT_AVAILABLE]
    with_catalog = [row for row in rows if row["local_catalog"] != NOT_AVAILABLE]
    panel_with_catalog = [
        row for row in panel if row["local_catalog"] != NOT_AVAILABLE
    ]
    counts = {
        status: sum(1 for row in rows if row["status"] == status)
        for status in STATUS_ORDER
    }
    print(
        f"genomas: {len(rows)} (panel GCF: {len(panel)}, "
        f"referencias: {len(references)})"
    )
    print(f"con especie: {sum(row['species'] != NOT_AVAILABLE for row in rows)}")
    print(f"con taxid: {len(with_taxid)} (panel: {len(panel_with_taxid)})")
    print(
        f"con catálogo local: {len(with_catalog)} "
        f"(panel: {len(panel_with_catalog)})"
    )
    print("estados: " + ", ".join(f"{s}={counts[s]}" for s in STATUS_ORDER))
    print(f"escrito: {out_tsv}")
    print(f"escrito: {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
