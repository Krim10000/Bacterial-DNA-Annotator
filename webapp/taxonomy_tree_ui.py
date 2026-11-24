"""Árbol taxonómico jerárquico a partir de la taxonomía local precalculada.

Alcance
-------
- Lee **solo** los JSON locales ``<data_root>/taxonomy/full/*.json`` (y, como
  respaldo, ``<data_root>/taxonomy/*.json``). No consulta NCBI, GTDB, UniProt
  ni ninguna API, y no ejecuta fastANI/GTDB-Tk.
- Construye un árbol ``domain -> phylum -> class -> order -> family -> genus
  -> species`` usando únicamente los rangos presentes. No inventa rangos
  ausentes ni fusiona linajes contradictorios: elige GTDB si tiene algún
  rango utilizable y, si no, NCBI (misma preferencia que
  :mod:`webapp.taxonomy_layer`).
- La interfaz es Streamlit y dibuja un árbol textual/indentado. No agrega
  dependencias (no usa Plotly ni Graphviz).

Límite central: la taxonomía es **precalculada**. El árbol no es ANI en vivo,
no valida especie ni filo y no convierte una ausencia de rango en evidencia
negativa. Cada vista lo declara.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "LIMITATION",
    "NOT_ANI_WARNING",
    "PREFIXES",
    "RANKS",
    "RANK_LABELS",
    "ROOT_NAME",
    "build_tree",
    "default_data_root",
    "find_path",
    "inventory",
    "load_lineages",
    "parse_lineage",
    "render",
    "tree_lines",
    "tree_text",
]

#: Orden canónico de los rangos que se muestran en el árbol.
RANKS: tuple[str, ...] = (
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)

RANK_LABELS = {
    "domain": "Dominio",
    "phylum": "Filo",
    "class": "Clase",
    "order": "Orden",
    "family": "Familia",
    "genus": "Género",
    "species": "Especie",
    "root": "Raíz",
}

PREFIX_TO_RANK = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
PREFIXES = {rank: prefix for prefix, rank in PREFIX_TO_RANK.items()}

ROOT_NAME = "Taxonomía"
NEW_GENOME_ID = "__new_genome__"

#: Valores que significan "sin dato"; nunca se convierten en un taxón.
_PLACEHOLDERS = {
    "",
    "unknown",
    "unclassified",
    "unranked",
    "na",
    "n/a",
    "none",
    "null",
    "not available",
    "not_available",
    "sin dato",
    "sin datos",
}

#: Advertencia que debe acompañar siempre al árbol.
LIMITATION = (
    "Taxonomía **precalculada**: el árbol se arma con JSON locales "
    "(GTDB/NCBI) y no es ANI en vivo. No se recalcula similitud, no se "
    "consulta la red y un rango ausente no es evidencia negativa."
)

NOT_ANI_WARNING = (
    "El ANI que se muestra (si existe) también es precalculado; no se "
    "recalcula al abrir esta vista."
)


# --------------------------------------------------------------------- utils
def _clean(value: Any) -> str:
    """Devuelve texto limpio o ``""`` para valores que significan "sin dato"."""
    if value is None:
        return ""
    text = str(value).strip()
    if text.casefold() in _PLACEHOLDERS:
        return ""
    return text


def _to_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def default_data_root() -> Path:
    """``DNA_DATA_ROOT`` si está definido; si no, ``<repo>/data``."""
    env = os.environ.get("DNA_DATA_ROOT")
    if env:
        return Path(env).expanduser()
    return Path(__file__).resolve().parents[1] / "data"


def _resolve_data_dir(data_root: Path | str | None) -> Path:
    """Acepta el directorio ``data`` o la raíz del proyecto."""
    root = Path(data_root).expanduser() if data_root is not None else default_data_root()
    if not root.is_absolute():
        root = (Path.cwd() / root).resolve()
    else:
        root = root.resolve()
    if (root / "taxonomy").is_dir():
        return root
    nested = root / "data"
    if (nested / "taxonomy").is_dir():
        return nested.resolve()
    return root


def _relative(path: Path, base: Path) -> str:
    try:
        return path.resolve().relative_to(base.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _is_taxonomy_record(data: Any) -> bool:
    """True si el JSON parece un registro de genoma y no un resumen agregado."""
    if not isinstance(data, Mapping):
        return False
    return any(key in data for key in (
        "genome", "gtdb_lineage", "ncbi_lineage", "species",
        "species_ncbi", "species_gtdb",
    ))


# ------------------------------------------------------------------ linajes
def _pairs_from_rank_mapping(mapping: Mapping[str, Any]) -> list[tuple[str, str]]:
    """Extrae ``(rank, taxón)`` de un diccionario indexado por rango."""
    rows: list[tuple[str, str]] = []
    seen: set[str] = set()
    for rank in RANKS:
        value = mapping.get(rank)
        if rank == "species" and not _clean(value):
            value = mapping.get("species_ncbi")
        if rank == "species" and not _clean(value):
            value = mapping.get("species_gtdb")
        taxon = _clean(value)
        if taxon and rank not in seen:
            rows.append((rank, taxon))
            seen.add(rank)
    return rows


def _pairs_from_value(value: Any) -> list[tuple[str, str]]:
    """Normaliza un linaje desde string, dict, lista de nombres o de pares.

    Solo acepta tokens con prefijo de rango (``d__``, ``p__``, ...) o un
    diccionario indexado por rango. Un string plano separado por ``;``/``>``
    se interpreta como nombres en el orden canónico solo si todos los tramos
    no están vacíos.
    """
    if value is None:
        return []

    if isinstance(value, Mapping):
        return _pairs_from_rank_mapping(value)

    if isinstance(value, str):
        text = _clean(value)
        if not text:
            return []
        rows: list[tuple[str, str]] = []
        seen: set[str] = set()
        for token in text.split(";"):
            token = token.strip()
            if len(token) < 4 or token[1:3] != "__":
                continue
            rank = PREFIX_TO_RANK.get(token[0].lower())
            taxon = _clean(token[3:])
            if rank is None or not taxon or rank in seen:
                continue
            seen.add(rank)
            rows.append((rank, taxon))
        if rows:
            return rows

        # Respaldo: "Bacteria > Pseudomonadota > ..." sin prefijos.
        parts = [part.strip() for chunk in text.split(";")
                 for part in chunk.split(">")]
        parts = [_clean(part) for part in parts]
        if parts and all(parts) and len(parts) <= len(RANKS):
            return list(zip(RANKS, parts))
        return []

    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        items = list(value)
        if not items:
            return []
        if all(isinstance(item, Mapping) for item in items):
            rows = []
            for item in items:
                rank = _clean(item.get("rank")).casefold()
                taxon = _clean(item.get("name") or item.get("taxon"))
                if rank in RANKS and taxon:
                    rows.append((rank, taxon))
            return _dedupe_pairs(rows)
        if all(isinstance(item, (tuple, list)) and len(item) == 2
               for item in items):
            rows = []
            for rank, taxon in items:
                rank_text = _clean(rank).casefold()
                taxon_text = _clean(taxon)
                if rank_text in RANKS and taxon_text:
                    rows.append((rank_text, taxon_text))
            return _dedupe_pairs(rows)
        if all(isinstance(item, str) for item in items):
            joined = ";".join(str(item) for item in items)
            parsed = _pairs_from_value(joined)
            if parsed:
                return parsed
            names = [_clean(item) for item in items]
            if names and all(names) and len(names) <= len(RANKS):
                return list(zip(RANKS, names))
        return []

    return []


def _dedupe_pairs(rows: Iterable[tuple[str, str]]) -> list[tuple[str, str]]:
    result: list[tuple[str, str]] = []
    seen: set[str] = set()
    for rank, taxon in rows:
        if rank in RANKS and taxon and rank not in seen:
            seen.add(rank)
            result.append((rank, taxon))
    return result


def parse_lineage(value: Any) -> list[tuple[str, str]]:
    """Convierte un linaje local en pares ``(rank, taxón)``.

    Acepta el string canónico ``d__Bacteria;p__...``, un diccionario por
    rango o una lista de nombres/registros. No adivina el rango de un token
    sin prefijo salvo en el respaldo de nombres planos documentado.
    """
    return _pairs_from_value(value)


def _canonical_lineage(ranks: Mapping[str, str]) -> str:
    return ";".join(
        f"{PREFIXES[rank]}__{ranks[rank]}"
        for rank in RANKS
        if _clean(ranks.get(rank))
    )


def _lineage_names(value: Any) -> list[str]:
    pairs = _pairs_from_value(value)
    if pairs:
        return [taxon for _, taxon in pairs]
    if isinstance(value, str):
        text = _clean(value)
        if text:
            parts = [part.strip() for chunk in text.split(";")
                     for part in chunk.split(">")]
            parts = [_clean(part) for part in parts]
            if parts and all(parts):
                return parts
    return []


def _make_record(
    genome: str,
    ranks: Mapping[str, Any],
    data: Mapping[str, Any] | None = None,
    *,
    method: str = "",
    species_from_json: bool = False,
    source: str = "",
    layer: str = "unknown",
    external: bool = False,
) -> dict[str, Any] | None:
    data = data or {}
    ordered = {
        rank: _clean(ranks.get(rank))
        for rank in RANKS
        if _clean(ranks.get(rank))
    }
    if not ordered:
        return None
    species_from_json = bool(species_from_json)
    return {
        "genome": genome,
        "ranks": ordered,
        "lineage": _canonical_lineage(ordered),
        "species": ordered.get("species", ""),
        "method": method or "unknown",
        "species_from_json": species_from_json,
        "confidence": _clean(data.get("confidence")) or "unknown",
        "status": _clean(data.get("status")) or "unknown",
        "ani": _to_float(data.get("ani")),
        "closest_ref": _clean(data.get("closest_ref")) or None,
        "source": source,
        "layer": layer,
        "missing_ranks": tuple(rank for rank in RANKS if rank not in ordered),
        "n_ranks": len(ordered),
        "external": external,
    }


def _record_from_data(
    genome: str,
    data: Mapping[str, Any],
    source: str,
    *,
    layer: str,
    external: bool = False,
) -> dict[str, Any] | None:
    """Normaliza un JSON taxonómico local sin inventar rangos ni especie."""
    evidence = data.get("evidence")
    if not isinstance(evidence, Mapping):
        evidence = {}
    gtdb_evidence = evidence.get("gtdb")
    if not isinstance(gtdb_evidence, Mapping):
        gtdb_evidence = {}
    ncbi_evidence = evidence.get("ncbi")
    if not isinstance(ncbi_evidence, Mapping):
        ncbi_evidence = {}

    candidates: list[tuple[str, Any]] = [
        ("gtdb", data.get("gtdb_lineage")),
        ("gtdb", gtdb_evidence.get("lineage")),
        ("ncbi", data.get("ncbi_lineage")),
        ("ncbi", ncbi_evidence.get("lineage")),
        ("ncbi", ncbi_evidence.get("lineage_ranks")),
        ("ncbi", ncbi_evidence.get("ranks")),
        ("", data),
    ]

    method = ""
    pairs: list[tuple[str, str]] = []
    for candidate_method, candidate in candidates:
        pairs = _pairs_from_value(candidate)
        if pairs:
            method = candidate_method
            break
    if not pairs:
        return None

    ranks = dict(pairs)
    if not method:
        method = "species_json" if set(ranks) == {"species"} else "local_fields"
    species_from_json = method == "species_json"
    if "species" not in ranks:
        if method == "gtdb":
            species_candidates = (
                data.get("species_gtdb"),
                data.get("species"),
                data.get("species_ncbi"),
            )
        elif method == "ncbi":
            species_candidates = (
                data.get("species_ncbi"),
                data.get("species"),
                data.get("species_gtdb"),
            )
        else:
            species_candidates = (
                data.get("species"),
                data.get("species_ncbi"),
                data.get("species_gtdb"),
            )
        for value in species_candidates:
            species = _clean(value)
            if species:
                ranks["species"] = species
                species_from_json = True
                if not method:
                    method = "species_json"
                break

    return _make_record(
        genome,
        ranks,
        data,
        method=method,
        species_from_json=species_from_json,
        source=source,
        layer=layer,
        external=external,
    )


# --------------------------------------------------------------- load_lines
def _scan(
    data_root: Path | str | None,
    *,
    include_pilot: bool,
) -> tuple[dict[str, dict[str, Any]], Path, Path, Path]:
    root = _resolve_data_dir(data_root)
    taxonomy_dir = root / "taxonomy"
    full_dir = taxonomy_dir / "full"

    records: dict[str, dict[str, Any]] = {}
    sources = [(full_dir, "full")]
    if include_pilot:
        sources.append((taxonomy_dir, "pilot"))

    for directory, layer in sources:
        if not directory.is_dir():
            continue
        for path in sorted(directory.glob("*.json")):
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError, json.JSONDecodeError):
                continue
            if not _is_taxonomy_record(data):
                continue
            genome = _clean(data.get("genome")) or path.stem
            if not genome or genome in records:
                continue
            record = _record_from_data(
                genome,
                data,
                _relative(path, root),
                layer=layer,
            )
            if record is not None:
                records[genome] = record

    return dict(sorted(records.items())), root, taxonomy_dir, full_dir


def load_lineages(
    data_root: Path | str | None = None,
    *,
    include_pilot: bool = True,
) -> dict[str, dict[str, Any]]:
    """Carga los linajes locales disponibles, sin red.

    Devuelve ``{genoma: registro}`` ordenado por genoma. Los JSON de
    ``taxonomy/full`` tienen prioridad; ``taxonomy/*.json`` solo rellena
    genomas ausentes (los pilotos duplicados no se cuentan dos veces).

    Cada registro contiene ``genome``, ``ranks`` (solo rangos presentes),
    ``lineage`` canónico, ``species``, ``method`` (``gtdb``/``ncbi``/
    ``species_json``), ``confidence``, ``status``, ``ani``, ``closest_ref``,
    ``source``, ``layer``, ``missing_ranks`` y ``n_ranks``.

    Un JSON inválido, un resumen agregado o un registro sin ningún rango
    utilizable se omite; ``{}`` significa "sin taxonomía local".
    """
    return _scan(data_root, include_pilot=include_pilot)[0]


def inventory(data_root: Path | str | None = None) -> dict[str, Any]:
    """Inventario serializable de la taxonomía local (no crea nada)."""
    records, root, taxonomy_dir, full_dir = _scan(
        data_root, include_pilot=True
    )
    n_full = sum(1 for record in records.values() if record["layer"] == "full")
    n_pilot = len(records) - n_full

    if not root.exists():
        reason = f"data_root no existe: {root}"
    elif not taxonomy_dir.is_dir():
        reason = f"falta el directorio de taxonomía local: {taxonomy_dir}"
    elif not records:
        reason = (
            f"sin JSON taxonómicos válidos en {full_dir} ni en {taxonomy_dir}"
        )
    else:
        reason = ""

    return {
        "available": bool(records),
        "reason": reason,
        "data_root": str(root),
        "taxonomy_dir": str(taxonomy_dir),
        "full_dir": str(full_dir),
        "n_full": n_full,
        "n_pilot": n_pilot,
        "n_lineages": len(records),
    }


# ------------------------------------------------------------------ árbol
def _new_node(name: str, rank: str) -> dict[str, Any]:
    return {
        "name": name,
        "rank": rank,
        "children": [],
        "genomes": [],
        "external_genomes": [],
        "_index": {},
    }


def _coerce_record(genome: str, raw: Any) -> dict[str, Any] | None:
    if isinstance(raw, str):
        return _record_from_data(
            genome,
            {"gtdb_lineage": raw},
            source="",
            layer="unknown",
            external=False,
        )
    if not isinstance(raw, Mapping):
        return None

    external = bool(raw.get("external"))
    if isinstance(raw.get("ranks"), Mapping):
        return _make_record(
            genome,
            raw["ranks"],
            raw,
            method=_clean(raw.get("method")) or "unknown",
            species_from_json=bool(raw.get("species_from_json")),
            source=_clean(raw.get("source")),
            layer=_clean(raw.get("layer")) or "unknown",
            external=external,
        )

    direct = {rank: raw.get(rank) for rank in RANKS if _clean(raw.get(rank))}
    if direct:
        return _make_record(
            genome,
            direct,
            raw,
            method=_clean(raw.get("method")) or "unknown",
            species_from_json=bool(raw.get("species_from_json")),
            source=_clean(raw.get("source")),
            layer=_clean(raw.get("layer")) or "unknown",
            external=external,
        )

    return _record_from_data(
        genome,
        raw,
        _clean(raw.get("source")),
        layer=_clean(raw.get("layer")) or "unknown",
        external=external,
    )


def _coerce_records(lineages: Any) -> list[dict[str, Any]]:
    """Acepta ``{genoma: registro}`` o una lista de registros."""
    records: list[dict[str, Any]] = []
    if isinstance(lineages, Mapping):
        for genome, raw in lineages.items():
            record = _coerce_record(str(genome), raw)
            if record is not None:
                records.append(record)
    else:
        for index, raw in enumerate(lineages or []):
            if not isinstance(raw, Mapping):
                continue
            genome = _clean(raw.get("genome")) or f"genome_{index}"
            record = _coerce_record(genome, raw)
            if record is not None:
                records.append(record)
    records.sort(key=lambda record: record["genome"])
    return records


def _external_record(new_genome: str, new_lineage: Any) -> dict[str, Any] | None:
    pairs = _pairs_from_value(new_lineage)
    if not pairs:
        return None
    return _make_record(
        new_genome,
        dict(pairs),
        {"status": "new_lineage"},
        method="new_lineage",
        source="new_lineage",
        layer="external",
        external=True,
    )


def build_tree(
    lineages: Mapping[str, Any] | Iterable[Mapping[str, Any]],
    new_genome: str | None = None,
    new_lineage: Any = None,
) -> dict[str, Any]:
    """Construye el árbol jerárquico ``domain -> ... -> species``.

    ``lineages`` puede ser el resultado de :func:`load_lineages` o una lista
    de registros. Solo se agregan los rangos presentes: si falta ``order``,
    el hijo de ``class`` se enlaza directamente con ``family`` (no se crea un
    nodo "unknown").

    Si se pasa ``new_lineage``, se agrega una rama **externa** al árbol para
    mostrar su posición (salvo que ``new_genome`` ya exista en los datos
    locales). La raíz contiene ``n_genomes`` (locales) y ``n_external``.
    """
    records = _coerce_records(lineages)
    existing = {record["genome"]: record for record in records}

    new_genome_text = _clean(new_genome)
    new_pairs = _pairs_from_value(new_lineage)
    external_record: dict[str, Any] | None = None
    if new_pairs:
        external_id = new_genome_text or NEW_GENOME_ID
        if external_id not in existing:
            external_record = _make_record(
                external_id,
                dict(new_pairs),
                {"status": "new_lineage"},
                method="new_lineage",
                source="new_lineage",
                layer="external",
                external=True,
            )
            if external_record is not None:
                records.append(external_record)
                existing[external_id] = external_record

    root = _new_node(ROOT_NAME, "root")
    for record in records:
        node = root
        for rank in RANKS:
            taxon = _clean(record["ranks"].get(rank))
            if not taxon:
                continue
            key = taxon.casefold()
            child = node["_index"].get(key)
            if child is None:
                child = _new_node(taxon, rank)
                node["children"].append(child)
                node["_index"][key] = child
            node = child
        node["genomes"].append(record["genome"])
        if record.get("external"):
            node["external_genomes"].append(record["genome"])

    def finalize(node: dict[str, Any]) -> tuple[int, int]:
        real = len(node["genomes"]) - len(node["external_genomes"])
        external = len(node["external_genomes"])
        for child in node["children"]:
            child_real, child_external = finalize(child)
            real += child_real
            external += child_external
        node["children"].sort(key=lambda child: child["name"].casefold())
        node["n_genomes"] = real
        node["n_external"] = external
        node.pop("_index", None)
        return real, external

    finalize(root)
    root["new_genome"] = new_genome_text or None
    root["new_lineage"] = _canonical_lineage(dict(new_pairs)) if new_pairs else ""
    root["has_external"] = bool(external_record)
    return root


def _find_nodes_by_genome(
    tree: Mapping[str, Any], genome: str | None
) -> list[dict[str, Any]]:
    """Nodos reales del árbol (no copias) desde el dominio hasta el genoma."""
    target = _clean(genome)
    if not target:
        return []

    def walk(node: Mapping[str, Any], path: list[dict[str, Any]]):
        current = path + [node]
        if target in (node.get("genomes") or ()):
            return current
        for child in node.get("children") or ():
            found = walk(child, current)
            if found:
                return found
        return []

    # Se recorre desde los hijos para que el camino no incluya la raíz.
    for child in tree.get("children") or ():
        found = walk(child, [])
        if found:
            return found
    return []


def _find_named_path(
    tree: Mapping[str, Any], name: str
) -> list[dict[str, Any]]:
    """Primer camino (DFS) hasta un nodo con ese nombre, sin la raíz."""
    key = name.casefold()

    def walk(node: Mapping[str, Any], path: list[dict[str, Any]]):
        for child in node.get("children") or ():
            current = path + [child]
            if str(child.get("name", "")).casefold() == key:
                return current
            found = walk(child, current)
            if found:
                return found
        return []

    return walk(tree, [])


def _find_nodes_by_lineage(
    tree: Mapping[str, Any], lineage: Any
) -> list[dict[str, Any]]:
    """Nodos reales del árbol que coinciden con el prefijo del linaje."""
    names = _lineage_names(lineage)
    if not names:
        return []

    node = tree
    path: list[dict[str, Any]] = []
    start = 0

    # Si el linaje no empieza en el dominio (p. ej. solo ``species``), el
    # primer nombre se busca en cualquier profundidad y luego se desciende.
    first = names[0].casefold()
    child = next(
        (candidate for candidate in node.get("children") or ()
         if str(candidate.get("name", "")).casefold() == first),
        None,
    )
    if child is not None:
        path.append(child)
        node = child
        start = 1
    else:
        found = _find_named_path(tree, names[0])
        if not found:
            return []
        path.extend(found)
        node = found[-1]
        start = 1

    for name in names[start:]:
        key = name.casefold()
        child = next(
            (candidate for candidate in node.get("children") or ()
             if str(candidate.get("name", "")).casefold() == key),
            None,
        )
        if child is None:
            break
        path.append(child)
        node = child
    return path


def find_path(
    tree: Mapping[str, Any],
    genome: str | None = None,
    lineage: Any = None,
) -> list[str]:
    """Devuelve el camino de nombres desde el dominio hasta el nodo hallado.

    Busca primero por ``genome`` exacto y, si no aparece, por ``lineage``
    (string canónico, dict por rango o lista de nombres). Devuelve ``[]`` si
    no se encuentra nada; si el linaje coincide solo parcialmente, devuelve
    el prefijo que sí coincide. La raíz sintética no se incluye.
    """
    if genome:
        nodes = _find_nodes_by_genome(tree, genome)
        if nodes:
            return [str(node["name"]) for node in nodes]
    if lineage is not None:
        nodes = _find_nodes_by_lineage(tree, lineage)
        if nodes:
            return [str(node["name"]) for node in nodes]
    return []


def _node_label(node: Mapping[str, Any], hidden: int = 0) -> str:
    label = str(node.get("name", ""))
    n_genomes = int(node.get("n_genomes") or 0)
    n_external = int(node.get("n_external") or 0)
    if n_genomes:
        label += f" ({n_genomes})"
    if n_external:
        label += f" (+{n_external} nuevo)"
    if hidden > 0:
        label += f" · {hidden} rama{'s' if hidden != 1 else ''} más"
    return label


def tree_lines(
    tree: Mapping[str, Any],
    *,
    selected: str | None = None,
    new_genome: str | None = None,
    new_lineage: Any = None,
    only_paths: bool = False,
    max_lines: int | None = None,
) -> list[str]:
    """Renderiza el árbol como líneas de texto indentado.

    Es una función pura (no importa Streamlit). ``selected`` y
    ``new_genome``/``new_lineage`` agregan marcas ``◀ seleccionado`` y
    ``★ nuevo`` en el nodo terminal correspondiente. Con ``only_paths=True``
    se muestran solo las ramas de esos caminos y se anota cuántas ramas
    quedaron ocultas.
    """
    selected_nodes = _find_nodes_by_genome(tree, selected)
    new_nodes = (
        _find_nodes_by_genome(tree, new_genome)
        if new_genome
        else []
    )
    if not new_nodes and new_lineage is not None:
        new_nodes = _find_nodes_by_lineage(tree, new_lineage)

    selected_ids = {id(node) for node in selected_nodes}
    new_ids = {id(node) for node in new_nodes}
    allowed = selected_ids | new_ids
    selected_terminal = id(selected_nodes[-1]) if selected_nodes else None
    new_terminal = id(new_nodes[-1]) if new_nodes else None

    n_genomes = int(tree.get("n_genomes") or 0)
    n_external = int(tree.get("n_external") or 0)
    root_label = f"{tree.get('name', ROOT_NAME)} ({n_genomes} genomas"
    if n_external:
        root_label += f" + {n_external} nuevo"
    root_label += ")"
    lines = [root_label]

    def markers(node: Mapping[str, Any]) -> str:
        marks = []
        if id(node) == selected_terminal:
            marks.append("◀ seleccionado")
        if id(node) == new_terminal:
            marks.append("★ nuevo")
        return ("  " + " · ".join(marks)) if marks else ""

    def walk(node: Mapping[str, Any], prefix: str, is_last: bool) -> None:
        children = list(node.get("children") or ())
        if only_paths:
            children = [child for child in children
                        if id(child) in allowed]
        hidden = len(node.get("children") or ()) - len(children)

        branch = "└─ " if is_last else "├─ "
        lines.append(prefix + branch + _node_label(node, hidden) + markers(node))
        child_prefix = prefix + ("   " if is_last else "│  ")
        for index, child in enumerate(children):
            walk(child, child_prefix, index == len(children) - 1)

    children = list(tree.get("children") or ())
    if only_paths:
        children = [child for child in children if id(child) in allowed]
    for index, child in enumerate(children):
        walk(child, "", index == len(children) - 1)

    if max_lines is not None and max_lines > 0 and len(lines) > max_lines:
        lines = lines[:max_lines]
        lines.append("… (árbol truncado)")
    return lines


def tree_text(
    tree: Mapping[str, Any],
    *,
    selected: str | None = None,
    new_genome: str | None = None,
    new_lineage: Any = None,
    only_paths: bool = False,
    max_lines: int | None = None,
) -> str:
    """Igual que :func:`tree_lines`, pero devuelve un único string."""
    return "\n".join(tree_lines(
        tree,
        selected=selected,
        new_genome=new_genome,
        new_lineage=new_lineage,
        only_paths=only_paths,
        max_lines=max_lines,
    ))


# ------------------------------------------------------------------- render
def _common_prefix(first: Sequence[str], second: Sequence[str]) -> list[str]:
    common: list[str] = []
    for left, right in zip(first, second):
        if left.casefold() != right.casefold():
            break
        common.append(left)
    return common


def render(
    data_root: Path | str | None = None,
    new_genome: str | None = None,
    new_lineage: Any = None,
    *,
    key: str = "taxonomy_tree",
) -> dict[str, Any]:
    """Dibuja la vista Streamlit del árbol taxonómico.

    Devuelve un resumen serializable (``available``, ``n_genomes``,
    ``selected``, ``path``, ``new_path``) útil para pruebas o para que otro
    componente sepa qué se mostró. No escribe archivos ni consulta la red.

    ``new_genome``/``new_lineage`` permiten señalar la posición de un genoma
    nuevo. Si ``new_lineage`` no está en los datos locales, se agrega una
    rama externa solo para la vista; nunca se guarda en disco.
    """
    import streamlit as st

    root = _resolve_data_dir(data_root)
    lineages = load_lineages(root)
    st.subheader("🌳 Árbol taxonómico")
    st.warning(LIMITATION)

    if not lineages:
        info = inventory(root)
        st.info(
            "Sin taxonomía local: no hay JSON utilizables en "
            "`taxonomy/full/` ni en `taxonomy/`."
        )
        st.caption(info["reason"])
        st.caption(
            "Se esperaba, por ejemplo, "
            f"`{Path(info['full_dir']) / 'GCF_*.json'}`. No se consulta la red."
        )
        return {
            "available": False,
            "reason": info["reason"],
            "n_genomes": 0,
            "selected": None,
            "path": [],
            "new_path": [],
        }

    new_genome_text = _clean(new_genome)
    new_pairs = _pairs_from_value(new_lineage)
    options = sorted(lineages)
    if new_genome_text and new_genome_text not in lineages and new_pairs:
        options.append(new_genome_text)

    default = new_genome_text if (new_genome_text in options) else options[0]
    st.caption(
        f"{len(lineages)} genomas con taxonomía local precalculada. "
        "El selector incluye todo `taxonomy/full/`."
    )
    selected = st.selectbox(
        "Genoma",
        options,
        index=options.index(default),
        key=f"{key}_genome",
        help="Genoma cuyo linaje y posición se resaltan en el árbol.",
    )

    tree = build_tree(
        lineages,
        new_genome=new_genome_text,
        new_lineage=new_lineage,
    )
    path = find_path(tree, genome=str(selected))
    record = lineages.get(str(selected))
    external_selected = False
    if record is None and new_genome_text and str(selected) == new_genome_text:
        record = _external_record(new_genome_text, new_lineage)
        external_selected = record is not None

    if record is None:
        st.warning(
            f"No hay linaje local para `{selected}`. No se puede ubicar en "
            "el árbol sin un JSON o un `new_lineage`."
        )
    else:
        columns = st.columns(4)
        columns[0].metric("Genoma", str(selected))
        columns[1].metric("Especie", record["species"] or "n/d")
        columns[2].metric(
            "ANI precalculado",
            "sin ANI" if record["ani"] is None else f"{record['ani']:.2f}",
        )
        columns[3].metric("Confianza", record["confidence"])
        if path:
            st.markdown("**Linaje:** " + " → ".join(path))
        else:
            st.warning("El genoma no aparece en el árbol construido.")
        origin = "nuevo / externo" if external_selected else "local"
        st.caption(
            f"Fuente: {record['source'] or 'n/d'} · método: "
            f"{record['method']} · estado: {record['status']} · origen: {origin}"
        )
        st.caption(NOT_ANI_WARNING)

        rank_rows = [
            {
                "Rango": RANK_LABELS.get(rank, rank),
                "Taxón": record["ranks"][rank],
            }
            for rank in RANKS
            if rank in record["ranks"]
        ]
        st.dataframe(rank_rows, hide_index=True)
        st.code(record["lineage"], language=None)

    new_path: list[str] = []
    if new_genome_text or new_pairs:
        new_path = find_path(
            tree,
            genome=new_genome_text,
            lineage=new_lineage,
        )
        if new_path:
            if new_genome_text:
                kind = "nuevo" if new_genome_text not in lineages else "local"
                label = f"Genoma {kind} `{new_genome_text}`"
            else:
                label = "Linaje nuevo"
            st.success(f"{label}: " + " → ".join(new_path))
            if path and selected != new_genome_text:
                common = _common_prefix(path, new_path)
                if common:
                    st.caption(
                        "Coincide con el seleccionado hasta "
                        f"**{common[-1]}**; después diverge."
                    )
                else:
                    st.caption(
                        "No comparte ningún nodo visible con el genoma "
                        "seleccionado."
                    )
        else:
            st.warning(
                "No se pudo ubicar `new_lineage` en el árbol. Revisa que "
                "tenga rangos con prefijo (`d__`, `p__`, ...) o un linaje "
                "canónico."
            )

    view = st.radio(
        "Vista del árbol",
        ("Ruta seleccionada", "Árbol completo"),
        horizontal=True,
        key=f"{key}_view",
    )
    only_paths = str(view).startswith("Ruta")
    lines = tree_lines(
        tree,
        selected=str(selected),
        new_genome=new_genome_text,
        new_lineage=new_lineage,
        only_paths=only_paths,
        max_lines=None if only_paths else 1500,
    )
    st.code("\n".join(lines), language=None)
    st.caption(LIMITATION)

    return {
        "available": True,
        "reason": "",
        "n_genomes": len(lineages),
        "selected": str(selected),
        "path": path,
        "new_path": new_path,
        "external": external_selected,
        "tree": tree,
    }
