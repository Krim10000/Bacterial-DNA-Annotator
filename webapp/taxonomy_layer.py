"""Capa de taxonomía local para un FASTA nuevo.

Normaliza la taxonomía **ya precalculada** en
``<data_root>/taxonomy/full/<genome>.json`` (o, como respaldo,
``<data_root>/taxonomy/<genome>.json``) a una tabla larga con las columnas
``taxon``, ``rank``, ``confidence``, ``method`` y ``source``.

Reglas de esta capa:

- **Solo datos locales.** No consulta NCBI, GTDB, UniProt ni ninguna API, y
  tampoco ejecuta fastANI/GTDB-Tk. Solo lee JSON y mapas locales.
- **No inventa especie ni linaje.** Si un rango no está en el JSON, no se
  escribe; si no hay JSON para el genoma, el resultado es
  ``status="not_available"`` con una razón explícita.
- **No recalcula confianza.** Copia la confianza declarada por el JSON; si
  falta, usa ``unknown``.
- **No fusiona NCBI y GTDB.** Elige un linaje (GTDB r232 si está disponible;
  si no, NCBI) y conserva en ``method`` de qué linaje viene cada fila. Las
  discrepancias entre ambos linajes quedan en el JSON original y en
  ``conflicts`` del resultado, no se mezclan en una fila.

Salida: ``<outdir>/<genome>.taxonomy.tsv``.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

FIELDS = ("taxon", "rank", "confidence", "method", "source")

RANK_ORDER = (
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)

PREFIX_TO_RANK = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
    "t": "strain",
}

METHOD_GTDB = "gtdb"
METHOD_NCBI = "ncbi"
METHOD_SPECIES = "species_json"
METHODS = (METHOD_GTDB, METHOD_NCBI, METHOD_SPECIES)

CONFIDENCE_VALUES = ("high", "medium", "low", "unknown")

# Valores que significan "no hay dato"; nunca se convierten en un taxón.
_PLACEHOLDERS = {"", "unknown", "unclassified", "unranked", "na", "n/a",
                 "none", "null", "not available", "not_available"}

LIMITATION = (
    "Taxonomía precalculada: el TSV normaliza un JSON local y no valida la "
    "especie por ANI ni por red. Un rango ausente no es evidencia negativa y "
    "una especie solo aparece si el JSON la declara."
)

# Mapas locales que `inventory` reconoce (no son necesarios para normalizar
# un JSON, pero se informan para auditar qué datos hay disponibles).
_ANI_MAPS = (
    ("taxonomy/ani_hits.tsv", "ani_hits_tsv"),
    ("taxonomy/full/ani_hits_full.tsv", "ani_hits_tsv"),
    ("taxonomy/full/ani_hits_full.jsonl", "ani_hits_jsonl"),
    ("taxonomy/full/ani_groups.json", "ani_groups"),
    ("taxonomy/full/scope_plan.tsv", "ani_scope_plan"),
)

_GENOME_MAPS = (
    ("selected.tsv", "selected_genomes"),
    ("taxonomy/refpanel.tsv", "reference_panel"),
    ("taxonomy/full/genome_meta.tsv", "genome_metadata"),
    ("taxonomy/full/panel_full.tsv", "reference_panel"),
    ("mobile_typing/pilot_genomes.tsv", "pilot_genomes"),
)


def _resolve_data_dir(data_root: Path | str) -> Path:
    """Acepta tanto el directorio ``data`` como la raíz del proyecto."""
    root = Path(data_root).expanduser().resolve()
    if (root / "taxonomy").is_dir():
        return root
    nested = root / "data"
    if (nested / "taxonomy").is_dir():
        return nested.resolve()
    return root


def _clean(value: Any) -> str:
    """Devuelve texto limpio o ``""`` para valores que significan "sin dato"."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in _PLACEHOLDERS:
        return ""
    return text


def _is_taxonomy_record(path: Path) -> bool:
    """True si el JSON parece un registro taxonómico y no un resumen agregado."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return False
    if not isinstance(data, dict):
        return False
    # Los registros reales tienen `genome`; se aceptan además JSON mínimos
    # con linaje o especie para no depender de un campo concreto.
    return any(key in data for key in
               ("genome", "gtdb_lineage", "ncbi_lineage", "species"))


def _taxonomy_files(directory: Path) -> list[Path]:
    if not directory.is_dir():
        return []
    return [p for p in sorted(directory.glob("*.json"))
            if _is_taxonomy_record(p)]


def _relative(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(path)


def _count_tsv_rows(path: Path) -> int:
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except OSError:
        return 0


def inventory(data_root: Path | str) -> dict[str, Any]:
    """Inventario de mapas locales de genomas, taxonomía y ANI.

    No modifica ni crea nada. Devuelve un diccionario serializable con:

    - ``available``: hay al menos un JSON taxonómico local.
    - ``n_full`` / ``n_pilot`` / ``n_taxonomy``: conteos de JSON válidos.
    - ``ani_maps`` / ``genome_maps``: listas de rutas y metadatos.
    - ``reason``: por qué no está disponible (vacío si sí lo está).
    """
    data_dir = _resolve_data_dir(data_root)
    taxonomy_dir = data_dir / "taxonomy"
    full_dir = taxonomy_dir / "full"

    full_files = _taxonomy_files(full_dir)
    pilot_files = _taxonomy_files(taxonomy_dir)
    n_taxonomy = len(full_files) + len(pilot_files)

    ani_maps: list[dict[str, Any]] = []
    for rel, kind in _ANI_MAPS:
        path = data_dir / rel
        if path.is_file():
            ani_maps.append({
                "path": rel,
                "kind": kind,
                "bytes": path.stat().st_size,
            })

    genome_maps: list[dict[str, Any]] = []
    for rel, kind in _GENOME_MAPS:
        path = data_dir / rel
        if path.is_file():
            genome_maps.append({
                "path": rel,
                "kind": kind,
                "n_rows": _count_tsv_rows(path),
            })

    if not data_dir.exists():
        reason = f"data_root no existe: {data_dir}"
    elif not taxonomy_dir.is_dir():
        reason = f"falta el directorio de taxonomía local: {taxonomy_dir}"
    elif n_taxonomy == 0:
        reason = (f"sin JSON de taxonomía en {full_dir} ni en "
                  f"{taxonomy_dir}")
    else:
        reason = ""

    return {
        "available": n_taxonomy > 0,
        "reason": reason,
        "data_root": str(data_dir),
        "taxonomy_dir": str(taxonomy_dir),
        "full_dir": str(full_dir),
        "n_full": len(full_files),
        "n_pilot": len(pilot_files),
        "n_taxonomy": n_taxonomy,
        "ani_maps": ani_maps,
        "genome_maps": genome_maps,
    }


def available(data_root: Path | str) -> tuple[bool, str]:
    """Indica si hay taxonomía local normalizable.

    Devuelve ``(True, "")`` si existe al menos un JSON de genoma bajo
    ``<data_root>/taxonomy/``; en caso contrario ``(False, razón)``.
    """
    info = inventory(data_root)
    return bool(info["available"]), str(info["reason"])


def _parse_lineage(lineage: Any) -> list[tuple[str, str]]:
    """Convierte ``d__Bacteria;p__...`` en pares ``(rank, taxón)``.

    Solo acepta tokens con prefijo de rango (``d__``, ``p__``, ...). Un token
    vacío, un prefijo desconocido o un texto sin prefijo se ignora: no se
    adivina el rango por posición. Los rangos repetidos conservan el primero.
    """
    text = _clean(lineage)
    if not text:
        return []
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    for token in text.split(";"):
        token = token.strip()
        if len(token) < 4 or token[1:3] != "__":
            continue
        rank = PREFIX_TO_RANK.get(token[0].lower())
        if rank is None:
            continue
        taxon = _clean(token[3:])
        if not taxon or rank in seen:
            continue
        seen.add(rank)
        rows.append((rank, taxon))
    return rows


def _ranks_from_evidence(data: dict[str, Any]) -> list[tuple[str, str]]:
    """Respaldo para JSON antiguos con ``evidence.ncbi`` en vez de linaje."""
    evidence = data.get("evidence")
    if not isinstance(evidence, dict):
        return []
    ncbi = evidence.get("ncbi")
    if not isinstance(ncbi, dict):
        return []

    rows: list[tuple[str, str]] = []
    lineage_ranks = ncbi.get("lineage_ranks")
    if isinstance(lineage_ranks, list):
        for item in lineage_ranks:
            if not isinstance(item, dict):
                continue
            rank = _clean(item.get("rank")).lower()
            taxon = _clean(item.get("name"))
            if rank in RANK_ORDER and taxon:
                rows.append((rank, taxon))

    ranks = ncbi.get("ranks")
    if isinstance(ranks, dict):
        for rank in RANK_ORDER:
            key = "species_ncbi" if rank == "species" else rank
            taxon = _clean(ranks.get(key))
            if taxon:
                rows.append((rank, taxon))

    deduped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for rank, taxon in rows:
        if rank in seen:
            continue
        seen.add(rank)
        deduped.append((rank, taxon))
    return deduped


def _fallback_species(data: dict[str, Any], method: str) -> str:
    """Especie declarada en el JSON, solo si el linaje elegido no la trae."""
    if method == METHOD_GTDB:
        candidates = (data.get("species_gtdb"), data.get("species"))
    elif method == METHOD_NCBI:
        candidates = (data.get("species_ncbi"), data.get("species"))
    else:
        candidates = (data.get("species"), data.get("species_ncbi"),
                      data.get("species_gtdb"))
    for value in candidates:
        species = _clean(value)
        if species:
            return species
    return ""


def normalize_taxonomy(data: dict[str, Any], source: str = "") -> list[dict[str, str]]:
    """Normaliza un registro taxonómico local a filas de cinco campos.

    Preferencia de linaje: GTDB r232 -> NCBI -> ``evidence.ncbi``. Si el
    linaje elegido no incluye especie, se usa el campo ``species`` del JSON
    (o ``species_ncbi``/``species_gtdb`` según el método) y se marca con
    ``method="species_json"``. No se completa ningún rango ausente.
    """
    if not isinstance(data, dict):
        raise TypeError("el registro taxonómico debe ser un diccionario")

    confidence = _clean(data.get("confidence")) or "unknown"
    evidence = data.get("evidence")
    if not isinstance(evidence, dict):
        evidence = {}
    gtdb_evidence = evidence.get("gtdb")
    if not isinstance(gtdb_evidence, dict):
        gtdb_evidence = {}
    ncbi_evidence = evidence.get("ncbi")
    if not isinstance(ncbi_evidence, dict):
        ncbi_evidence = {}

    parsed: list[tuple[str, str]] = []
    method = ""

    # 1) GTDB: campo principal y, si falta, evidencia GTDB local.
    for candidate in (data.get("gtdb_lineage"), gtdb_evidence.get("lineage")):
        parsed = _parse_lineage(candidate)
        if parsed:
            method = METHOD_GTDB
            break

    # 2) NCBI: campo principal y, si falta, evidencia NCBI local.
    if not parsed:
        for candidate in (data.get("ncbi_lineage"), ncbi_evidence.get("lineage")):
            parsed = _parse_lineage(candidate)
            if parsed:
                method = METHOD_NCBI
                break

    # 3) Respaldo para JSON sin cadena de linaje.
    if not parsed:
        parsed = _ranks_from_evidence(data)
        if parsed:
            method = METHOD_NCBI

    rows = [
        {
            "taxon": taxon,
            "rank": rank,
            "confidence": confidence,
            "method": method,
            "source": source,
        }
        for rank, taxon in parsed
    ]

    # Especie: solo desde un campo local del JSON y solo si el linaje no la
    # trae. No se deduce del género ni se consulta ninguna base externa.
    if not any(row["rank"] == "species" for row in rows):
        species = _fallback_species(data, method)
        if species:
            rows.append({
                "taxon": species,
                "rank": "species",
                "confidence": confidence,
                "method": METHOD_SPECIES,
                "source": source,
            })

    order = {rank: i for i, rank in enumerate(RANK_ORDER)}
    rows.sort(key=lambda row: order.get(row["rank"], len(RANK_ORDER)))
    return rows


def _not_available(genome: str, reason: str) -> dict[str, Any]:
    return {
        "genome": genome,
        "status": "not_available",
        "ok": False,
        "reason": reason,
        "rows": [],
        "n_rows": 0,
        "tsv": "",
        "source": "",
    }


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        with tmp.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=list(FIELDS),
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(path)
    finally:
        if tmp.exists():
            tmp.unlink()


def run(genome: str, data_root: Path | str,
        outdir: Path | str) -> dict[str, Any]:
    """Normaliza la taxonomía local de ``genome`` y escribe el TSV.

    Busca ``taxonomy/full/<genome>.json`` y, si no existe,
    ``taxonomy/<genome>.json``. Si no hay JSON, no es un JSON válido o no
    contiene ningún rango/especie utilizable, devuelve
    ``status="not_available"`` con ``reason`` y **no** escribe archivo.

    En éxito devuelve ``status="ok"``, ``rows``, ``n_rows``, ``tsv``,
    ``taxonomy_status`` (estado del JSON: assigned/provisional/...), la
    confianza, la especie declarada, ANI/referencia y conflictos, sin
    recalcular ni inventar nada.
    """
    genome = str(genome).strip()
    if not genome or "/" in genome or "\\" in genome:
        return _not_available(genome, f"genoma inválido: {genome!r}")

    data_dir = _resolve_data_dir(data_root)
    taxonomy_dir = data_dir / "taxonomy"
    full_dir = taxonomy_dir / "full"

    candidates = (full_dir / f"{genome}.json",
                  taxonomy_dir / f"{genome}.json")
    source = next((path for path in candidates if path.is_file()), None)
    if source is None:
        return _not_available(
            genome,
            f"sin taxonomía local para {genome}: falta "
            f"{full_dir / (genome + '.json')} (y "
            f"{taxonomy_dir / (genome + '.json')})",
        )

    try:
        data = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return _not_available(
            genome,
            f"JSON de taxonomía inválido o ilegible en {source}: {exc}",
        )
    if not isinstance(data, dict):
        return _not_available(
            genome,
            f"el JSON de taxonomía no es un objeto: {source}",
        )

    source_rel = _relative(source, data_dir)
    rows = normalize_taxonomy(data, source=source_rel)
    if not rows:
        return _not_available(
            genome,
            f"el JSON de taxonomía no contiene linaje ni especie utilizable: "
            f"{source_rel}",
        )

    outdir = Path(outdir)
    tsv = outdir / f"{genome}.taxonomy.tsv"
    _write_tsv(tsv, rows)

    evidence = data.get("evidence")
    if not isinstance(evidence, dict):
        evidence = {}
    conflicts = evidence.get("conflicts")
    if not isinstance(conflicts, list):
        conflicts = []

    return {
        "genome": genome,
        "status": "ok",
        "ok": True,
        "reason": "",
        "rows": rows,
        "n_rows": len(rows),
        "tsv": str(tsv),
        "source": str(source),
        "source_rel": source_rel,
        "taxonomy_status": _clean(data.get("status")) or "unknown",
        "confidence": _clean(data.get("confidence")) or "unknown",
        "species": _clean(data.get("species")) or None,
        "closest_ref": data.get("closest_ref"),
        "ani": data.get("ani"),
        "method": rows[0].get("method", ""),
        "conflicts": conflicts,
        "limitation": LIMITATION,
    }
