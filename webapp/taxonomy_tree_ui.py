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
- La interfaz es Streamlit y ofrece un árbol textual/indentado y una vista
  gráfica SVG generada en Python (sin Plotly, Graphviz ni otras
  dependencias): columnas por rango, resaltado del genoma seleccionado y de
  la rama nueva, y descarga SVG/HTML. No consulta la red.

Límite central: la taxonomía es **precalculada**. El árbol no es ANI en vivo,
no valida especie ni filo y no convierte una ausencia de rango en evidencia
negativa. Cada vista lo declara.
"""
from __future__ import annotations

import html
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "HTML_MIME",
    "LIMITATION",
    "NOT_ANI_WARNING",
    "PREFIXES",
    "RANKS",
    "RANK_LABELS",
    "ROOT_NAME",
    "SVG_MIME",
    "SVG_NS",
    "build_tree",
    "default_data_root",
    "find_path",
    "inventory",
    "load_lineages",
    "parse_lineage",
    "render",
    "tree_html",
    "tree_layout",
    "tree_lines",
    "tree_svg",
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
        declared_method = _clean(data.get("method"))
        if declared_method:
            method = declared_method
        elif set(ranks) == {"species"}:
            method = "species_json"
        else:
            method = "local_fields"
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

    # No se usa un atajo con ``raw["species"]``: un registro crudo puede
    # traer ``species`` y ``gtdb_lineage`` a la vez, y el linaje completo
    # debe ganar. ``_record_from_data`` ya entiende campos por rango.
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


def _graph_node_label(node: Mapping[str, Any]) -> str:
    """Etiqueta compacta para el SVG: nombre y conteo local, sin sufijos.

    El marcador de rama externa se dibuja aparte (insignia ``★``) para no
    ensanchar todas las columnas del gráfico.
    """
    label = str(node.get("name", ""))
    n_genomes = int(node.get("n_genomes") or 0)
    if n_genomes:
        label += f" ({n_genomes})"
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


# ------------------------------------------------------------------ gráfico
SVG_NS = "http://www.w3.org/2000/svg"
SVG_MIME = "image/svg+xml"
HTML_MIME = "text/html"

_SVG_CSS = """
.taxonomy-tree .background { fill: #ffffff; }
.taxonomy-tree .title { fill: #0f172a;
                        font: 700 16px 'Segoe UI', system-ui, sans-serif; }
.taxonomy-tree .subtitle { fill: #475569;
                           font: 11px 'Segoe UI', system-ui, sans-serif; }
.taxonomy-tree .rank-label { fill: #64748b;
                             font: 700 9px 'Segoe UI', system-ui, sans-serif;
                             letter-spacing: .06em;
                             text-transform: uppercase; }
.taxonomy-tree .rank-guide { stroke: #e2e8f0; stroke-width: 1;
                             stroke-dasharray: 2 6; }
.taxonomy-tree .edge { fill: none; stroke: #94a3b8; stroke-width: 1.2;
                       stroke-linecap: round; stroke-linejoin: round; }
.taxonomy-tree .edge.selected { stroke: #1d4ed8; stroke-width: 2.6; }
.taxonomy-tree .edge.new { stroke: #15803d; stroke-width: 2.6;
                           stroke-dasharray: 7 4; }
.taxonomy-tree .edge.selected.new { stroke: #7c3aed; }
.taxonomy-tree .node circle { fill: #ffffff; stroke: #475569;
                              stroke-width: 1.2; }
.taxonomy-tree .node.root circle { fill: #0f172a; stroke: #0f172a;
                                   stroke-width: 2.4; }
.taxonomy-tree .node.selected circle { fill: #dbeafe; stroke: #1d4ed8;
                                       stroke-width: 2.6; }
.taxonomy-tree .node.new circle { fill: #dcfce7; stroke: #15803d;
                                  stroke-width: 2.6; }
.taxonomy-tree .node.selected.new circle { fill: #ede9fe;
                                           stroke: #7c3aed; }
.taxonomy-tree .node-label { fill: #0f172a;
                             font: 11px 'Segoe UI', system-ui, sans-serif;
                             paint-order: stroke; stroke: #ffffff;
                             stroke-width: 3px; stroke-linejoin: round; }
.taxonomy-tree .node.selected .node-label { fill: #1e3a8a;
                                            font-weight: 700; }
.taxonomy-tree .node.new .node-label { fill: #14532d; font-weight: 700; }
.taxonomy-tree .node.selected.new .node-label { fill: #6d28d9;
                                               font-weight: 700; }
.taxonomy-tree .node.root .node-label { fill: #0f172a;
                                        font-weight: 700; }
.taxonomy-tree .badge { font: 700 10px 'Segoe UI', system-ui, sans-serif;
                        paint-order: stroke; stroke: #ffffff;
                        stroke-width: 2.5px; stroke-linejoin: round; }
.taxonomy-tree .badge.selected { fill: #1d4ed8; }
.taxonomy-tree .badge.new { fill: #15803d; }
.taxonomy-tree .legend text { fill: #475569;
                              font: 10px 'Segoe UI', system-ui, sans-serif; }
.taxonomy-tree .footer { fill: #64748b;
                         font: 10px 'Segoe UI', system-ui, sans-serif; }
.taxonomy-tree .empty { fill: #64748b;
                        font: italic 12px 'Segoe UI', system-ui, sans-serif; }
"""


def _xml_escape(value: Any) -> str:
    """Escapa texto/atributos para que el SVG siempre sea XML válido."""
    return html.escape(str(value), quote=True)


def _estimate_text_width(text: str, font_size: float) -> float:
    """Ancho aproximado de un texto sin depender de fuentes del sistema.

    Es una estimación deliberadamente holgada (0,62 em por carácter) para
    reservar espacio a nombres largos y acentos; no se usa para recortar.
    """
    return max(len(text), 1) * font_size * 0.62


def _layout_ids(
    raw_nodes: Sequence[Mapping[str, Any]],
    raw_to_layout: Mapping[int, int],
) -> list[int]:
    return [
        raw_to_layout[id(node)]
        for node in raw_nodes
        if id(node) in raw_to_layout
    ]


def tree_layout(
    tree: Mapping[str, Any],
    *,
    selected: str | None = None,
    new_genome: str | None = None,
    new_lineage: Any = None,
    width: float = 1200.0,
    row_height: float = 22.0,
    margin: float = 24.0,
    column_width: float = 140.0,
    font_size: float = 11.0,
) -> dict[str, Any]:
    """Calcula la geometría del árbol horizontal (sin Streamlit, pura).

    La raíz queda a la izquierda y cada rango ocupa una columna
    ``Dominio`` … ``Especie``. El ancho de cada columna se adapta al nombre
    más largo de ese rango (con ``column_width`` como mínimo), para que las
    etiquetas no se pisen. Los rangos ausentes **no** crean columnas ni
    nodos intermedios: el hijo simplemente salta a la columna de su rango.
    Cada hoja recibe una ``y`` consecutiva (sin solapamiento) y cada nodo
    interno se centra entre su primera y su última hoja. ``width`` es un
    mínimo: el resultado se ensancha si un nombre o una insignia lo exige.

    Devuelve un diccionario serializable con ``nodes``, ``edges``,
    ``rank_columns``, ``paths`` (IDs del genoma seleccionado y de la rama
    nueva), ``width``, ``height`` y contadores. No muta ``tree``.
    """
    root_raw: Mapping[str, Any]
    if isinstance(tree, Mapping):
        root_raw = tree
    else:  # pragma: no cover - defensa para entradas externas inválidas
        root_raw = {}
    if not root_raw.get("name"):
        root_raw = dict(root_raw)
        root_raw["name"] = ROOT_NAME
        root_raw.setdefault("rank", "root")

    width = float(width) if width and width > 0 else 1200.0
    row_height = (
        float(row_height) if row_height is not None and row_height >= 8
        else 8.0
    )
    margin = float(margin) if margin is not None and margin >= 0 else 24.0
    column_width = (
        float(column_width) if column_width is not None and column_width >= 40
        else 40.0
    )
    font_size = (
        float(font_size) if font_size is not None and font_size >= 6
        else 6.0
    )

    title_y = margin + 12.0
    subtitle_y = margin + 30.0
    rank_header_y = margin + 54.0
    legend_y = margin + 78.0
    header_height = margin + 96.0
    top = header_height + row_height / 2.0

    # Caminos de resaltado antes de medir columnas: las insignias del nodo
    # terminal también necesitan espacio para no pisar la columna siguiente.
    selected_raw = _find_nodes_by_genome(tree, selected) if selected else []
    if new_genome:
        new_raw = _find_nodes_by_genome(tree, new_genome)
    else:
        new_raw = []
    if not new_raw and new_lineage is not None:
        new_raw = _find_nodes_by_lineage(tree, new_lineage)

    rank_extent: dict[str, float] = {}
    raw_info: dict[int, tuple[str, float]] = {}

    def measure(raw: Mapping[str, Any], is_root: bool) -> None:
        rank = _clean(raw.get("rank")) or ("root" if is_root else "")
        label_width = _estimate_text_width(_graph_node_label(raw), font_size)
        rank_extent[rank] = max(
            rank_extent.get(rank, 0.0), label_width + 22.0
        )
        raw_info[id(raw)] = (rank, label_width)
        for child in raw.get("children") or ():
            if isinstance(child, Mapping):
                measure(child, False)

    measure(root_raw, True)

    badge_texts: dict[int, list[str]] = {}
    if selected_raw:
        badge_texts.setdefault(id(selected_raw[-1]), []).append(
            f"◀ {_clean(selected)}"
        )
    if new_raw:
        new_name = _clean(new_genome)
        badge_texts.setdefault(id(new_raw[-1]), []).append(
            f"★ {new_name}" if new_name else "★ nuevo"
        )
    for raw_terminal, texts in badge_texts.items():
        rank, label_width = raw_info.get(raw_terminal, ("", 0.0))
        extent = label_width + 22.0
        for text in texts:
            extent += 14.0 + _estimate_text_width(text, font_size * 0.9)
        rank_extent[rank] = max(rank_extent.get(rank, 0.0), extent)

    x_positions: dict[str, float] = {"root": margin}
    cursor = margin + max(column_width, rank_extent.get("root", 0.0))
    for rank in RANKS:
        x_positions[rank] = cursor
        cursor += max(column_width, rank_extent.get(rank, 0.0))
    columns_right = cursor

    nodes: list[dict[str, Any]] = []
    raw_to_layout: dict[int, int] = {}
    leaf_cursor = [0]

    def visit(
        raw: Mapping[str, Any],
        parent_id: int | None,
        depth: int,
    ) -> int:
        node_id = len(nodes)
        raw_to_layout[id(raw)] = node_id
        rank = _clean(raw.get("rank")) or (
            "root" if parent_id is None else ""
        )
        if rank == "root":
            x = margin
        elif rank in x_positions:
            x = x_positions[rank]
        else:
            x = margin + column_width * max(depth, 1)
        if parent_id is not None:
            x = max(x, nodes[parent_id]["x"] + column_width * 0.55)

        node: dict[str, Any] = {
            "id": node_id,
            "name": str(raw.get("name", "")),
            "rank": rank,
            "depth": depth,
            "parent": parent_id,
            "children": [],
            "genomes": tuple(raw.get("genomes") or ()),
            "n_genomes": int(raw.get("n_genomes") or 0),
            "n_external": int(raw.get("n_external") or 0),
            "x": float(x),
            "y": 0.0,
            "is_leaf": False,
            "label": _graph_node_label(raw),
            "label_x": 0.0,
            "label_y": 0.0,
            "label_width": 0.0,
            "selected": False,
            "new": False,
            "badges": [],
        }
        nodes.append(node)

        child_ids: list[int] = []
        for child in raw.get("children") or ():
            if isinstance(child, Mapping):
                child_ids.append(visit(child, node_id, depth + 1))
        node["children"] = child_ids

        if child_ids:
            first_y = nodes[child_ids[0]]["y"]
            last_y = nodes[child_ids[-1]]["y"]
            node["y"] = (first_y + last_y) / 2.0
        else:
            node["y"] = top + leaf_cursor[0] * row_height
            leaf_cursor[0] += 1
            node["is_leaf"] = True
        return node_id

    root_id = visit(root_raw, None, 0)

    selected_ids = _layout_ids(selected_raw, raw_to_layout)
    new_ids = _layout_ids(new_raw, raw_to_layout)
    selected_set = set(selected_ids)
    new_set = set(new_ids)

    badge_map: dict[int, list[dict[str, Any]]] = {}
    if selected_ids:
        badge_map.setdefault(selected_ids[-1], []).append({
            "kind": "selected",
            "text": f"◀ {_clean(selected)}",
        })
    if new_ids:
        new_name = _clean(new_genome)
        badge_map.setdefault(new_ids[-1], []).append({
            "kind": "new",
            "text": f"★ {new_name}" if new_name else "★ nuevo",
        })

    content_right = 0.0
    for node in nodes:
        node["selected"] = node["id"] in selected_set
        node["new"] = node["id"] in new_set
        node["label_width"] = _estimate_text_width(
            node["label"], font_size
        )
        node["label_x"] = node["x"] + 9.0
        node["label_y"] = node["y"] + font_size * 0.36
        right = node["label_x"] + node["label_width"]
        badge_x = right + 14.0
        node["badges"] = badge_map.get(node["id"], [])
        for badge in node["badges"]:
            badge["x"] = badge_x
            badge["y"] = node["label_y"]
            badge_width = _estimate_text_width(badge["text"], font_size * 0.9)
            badge["width"] = badge_width
            badge_x += badge_width + 12.0
            right = badge_x
        content_right = max(content_right, right)

    edges: list[dict[str, Any]] = []
    for node in nodes:
        parent_id = node["parent"]
        if parent_id is None:
            continue
        parent = nodes[parent_id]
        midpoint = (parent["x"] + node["x"]) / 2.0
        edges.append({
            "parent": parent_id,
            "child": node["id"],
            "x1": parent["x"],
            "y1": parent["y"],
            "x2": node["x"],
            "y2": node["y"],
            "d": (
                f"M {parent['x']:.2f} {parent['y']:.2f} "
                f"H {midpoint:.2f} V {node['y']:.2f} "
                f"H {node['x']:.2f}"
            ),
            "selected": node["selected"],
            "new": node["new"],
        })

    minimum_width = columns_right + margin
    layout_width = max(width, minimum_width, content_right + margin)
    n_leaves = leaf_cursor[0]
    height = header_height + max(1, n_leaves) * row_height + margin

    rank_columns = [{
        "rank": "root",
        "label": RANK_LABELS["root"],
        "x": margin,
    }]
    for rank in RANKS:
        rank_columns.append({
            "rank": rank,
            "label": RANK_LABELS.get(rank, rank),
            "x": x_positions[rank],
        })

    return {
        "nodes": nodes,
        "edges": edges,
        "rank_columns": rank_columns,
        "paths": {"selected": selected_ids, "new": new_ids},
        "path_names": {
            "selected": [nodes[node_id]["name"] for node_id in selected_ids],
            "new": [nodes[node_id]["name"] for node_id in new_ids],
        },
        "terminals": {
            "selected": selected_ids[-1] if selected_ids else None,
            "new": new_ids[-1] if new_ids else None,
        },
        "root": root_id,
        "n_nodes": len(nodes),
        "n_edges": len(edges),
        "n_leaves": n_leaves,
        "n_genomes": int(tree.get("n_genomes") or 0)
        if isinstance(tree, Mapping) else 0,
        "n_external": int(tree.get("n_external") or 0)
        if isinstance(tree, Mapping) else 0,
        "empty": not any(node["children"] for node in nodes),
        "width": float(layout_width),
        "height": float(height),
        "margin": margin,
        "row_height": row_height,
        "column_width": column_width,
        "font_size": font_size,
        "header_height": header_height,
        "title_y": title_y,
        "subtitle_y": subtitle_y,
        "rank_header_y": rank_header_y,
        "legend_y": legend_y,
        "top": top,
    }


def tree_svg(
    tree: Mapping[str, Any],
    *,
    selected: str | None = None,
    new_genome: str | None = None,
    new_lineage: Any = None,
    width: float = 1200.0,
    row_height: float = 22.0,
    title: str = "Árbol taxonómico",
    description: str | None = None,
) -> str:
    """Genera un SVG autónomo del árbol, sin dependencias ni red.

    Usa :func:`tree_layout` para la geometría. Solo dibuja los taxones que
    existen en ``tree``: no agrega nodos "unknown" ni rangos ausentes. El
    camino del genoma ``selected`` se pinta en azul y la rama
    ``new_genome``/``new_lineage`` en verde discontinuo; para ver la rama
    nueva, ``tree`` debe haberse construido con
    :func:`build_tree` y esos mismos argumentos. Si no hay datos, devuelve
    un SVG válido con un estado vacío explícito.
    """
    layout = tree_layout(
        tree,
        selected=selected,
        new_genome=new_genome,
        new_lineage=new_lineage,
        width=width,
        row_height=row_height,
    )
    width_value = layout["width"]
    height_value = layout["height"]
    margin = layout["margin"]
    rank_header_y = layout["rank_header_y"]
    footer_y = height_value - 10.0

    if description is None:
        description = LIMITATION.replace("**", "")

    parts: list[str] = [
        (
            f'<svg xmlns="{SVG_NS}" class="taxonomy-tree" '
            f'width="{width_value:.0f}" height="{height_value:.0f}" '
            f'viewBox="0 0 {width_value:.0f} {height_value:.0f}" '
            f'role="img" aria-labelledby="taxonomy-tree-title '
            f'taxonomy-tree-desc" data-nodes="{layout["n_nodes"]}" '
            f'data-leaves="{layout["n_leaves"]}">'
        ),
        f'<title id="taxonomy-tree-title">{_xml_escape(title)}</title>',
        (
            '<desc id="taxonomy-tree-desc">'
            f'{_xml_escape(description)}</desc>'
        ),
        f"<style>{_SVG_CSS}</style>",
        (
            f'<rect class="background" x="0" y="0" '
            f'width="{width_value:.0f}" height="{height_value:.0f}" '
            'fill="#ffffff" />'
        ),
        (
            f'<text class="title" x="{margin:.2f}" '
            f'y="{layout["title_y"]:.2f}" fill="#0f172a" '
            "font-family=\"'Segoe UI', system-ui, sans-serif\" "
            'font-size="16" font-weight="700">'
            f'{_xml_escape(title)}</text>'
        ),
    ]

    subtitle = f"{layout['n_genomes']} genomas locales"
    if layout["n_external"]:
        plural = "nueva" if layout["n_external"] == 1 else "nuevas"
        subtitle += f" + {layout['n_external']} rama {plural}"
    subtitle += " · taxonomía precalculada · sin red"
    parts.append(
        f'<text class="subtitle" x="{margin:.2f}" '
        f'y="{layout["subtitle_y"]:.2f}" fill="#475569" '
        "font-family=\"'Segoe UI', system-ui, sans-serif\" "
        f'font-size="11">{_xml_escape(subtitle)}</text>'
    )

    # Leyenda de resaltado (siempre visible, aunque no haya coincidencias).
    legend_y = layout["legend_y"]
    legend_x = margin
    parts.extend([
        (
            f'<g class="legend"><line x1="{legend_x:.2f}" '
            f'y1="{legend_y:.2f}" x2="{legend_x + 24:.2f}" '
            f'y2="{legend_y:.2f}" stroke="#1d4ed8" stroke-width="2.6" />'
            f'<text x="{legend_x + 30:.2f}" '
            f'y="{legend_y + 3:.2f}" fill="#475569" font-size="10" '
            "font-family=\"'Segoe UI', system-ui, sans-serif\">"
            "seleccionado</text></g>"
        ),
        (
            f'<g class="legend"><line x1="{legend_x + 125:.2f}" '
            f'y1="{legend_y:.2f}" x2="{legend_x + 149:.2f}" '
            f'y2="{legend_y:.2f}" stroke="#15803d" stroke-width="2.6" '
            f'stroke-dasharray="7 4" />'
            f'<text x="{legend_x + 155:.2f}" '
            f'y="{legend_y + 3:.2f}" fill="#475569" font-size="10" '
            "font-family=\"'Segoe UI', system-ui, sans-serif\">"
            "nuevo / externo</text></g>"
        ),
        (
            f'<g class="legend"><line x1="{legend_x + 268:.2f}" '
            f'y1="{legend_y:.2f}" x2="{legend_x + 292:.2f}" '
            f'y2="{legend_y:.2f}" stroke="#7c3aed" stroke-width="2.6" '
            f'stroke-dasharray="7 4" />'
            f'<text x="{legend_x + 298:.2f}" '
            f'y="{legend_y + 3:.2f}" fill="#475569" font-size="10" '
            "font-family=\"'Segoe UI', system-ui, sans-serif\">"
            "coincide</text></g>"
        ),
    ])

    if layout["empty"]:
        parts.append(
            f'<text class="empty" x="{width_value / 2.0:.2f}" '
            f'y="{layout["top"] + 12.0:.2f}" text-anchor="middle" '
            'fill="#64748b" font-size="12" font-style="italic" '
            "font-family=\"'Segoe UI', system-ui, sans-serif\">"
            "Sin genomas con taxonomía local. No se inventan taxones."
            "</text>"
        )
    else:
        for column in layout["rank_columns"]:
            x = column["x"]
            parts.append(
                f'<line class="rank-guide" x1="{x:.2f}" '
                f'y1="{layout["header_height"]:.2f}" x2="{x:.2f}" '
                f'y2="{height_value - margin:.2f}" stroke="#e2e8f0" '
                'stroke-width="1" stroke-dasharray="2 6" />'
            )
            parts.append(
                f'<text class="rank-label" x="{x:.2f}" '
                f'y="{rank_header_y:.2f}" fill="#64748b" font-size="9" '
                'font-weight="700" letter-spacing=".06em" '
                "font-family=\"'Segoe UI', system-ui, sans-serif\">"
                f'{_xml_escape(column["label"])}</text>'
            )

        for edge in layout["edges"]:
            classes = ["edge"]
            if edge["selected"] and edge["new"]:
                classes.extend(("selected", "new"))
                edge_stroke = "#7c3aed"
                edge_width = "2.6"
                edge_dash = ' stroke-dasharray="7 4"'
            elif edge["selected"]:
                classes.append("selected")
                edge_stroke = "#1d4ed8"
                edge_width = "2.6"
                edge_dash = ""
            elif edge["new"]:
                classes.append("new")
                edge_stroke = "#15803d"
                edge_width = "2.6"
                edge_dash = ' stroke-dasharray="7 4"'
            else:
                edge_stroke = "#94a3b8"
                edge_width = "1.2"
                edge_dash = ""
            parts.append(
                f'<path class="{" ".join(classes)}" '
                f'data-parent="{edge["parent"]}" '
                f'data-child="{edge["child"]}" fill="none" '
                f'stroke="{edge_stroke}" stroke-width="{edge_width}" '
                f'stroke-linecap="round" stroke-linejoin="round"'
                f'{edge_dash} d="{edge["d"]}" />'
            )

        for node in layout["nodes"]:
            classes = ["node"]
            if node["rank"]:
                classes.append(f"rank-{node['rank']}")
            if node["id"] == layout["root"]:
                classes.append("root")
            if node["selected"]:
                classes.append("selected")
            if node["new"]:
                classes.append("new")
            radius = 4.2
            if node["id"] in {
                layout["terminals"]["selected"],
                layout["terminals"]["new"],
            }:
                radius = 5.4
            is_root = node["id"] == layout["root"]
            if is_root:
                radius = 5.8
            if node["selected"] and node["new"]:
                fill, stroke = "#ede9fe", "#7c3aed"
                text_fill = "#6d28d9"
            elif node["selected"]:
                fill, stroke = "#dbeafe", "#1d4ed8"
                text_fill = "#1e3a8a"
            elif node["new"]:
                fill, stroke = "#dcfce7", "#15803d"
                text_fill = "#14532d"
            elif is_root:
                fill, stroke = "#0f172a", "#0f172a"
                text_fill = "#0f172a"
            else:
                fill, stroke = "#ffffff", "#475569"
                text_fill = "#0f172a"
            stroke_width = "2.4" if is_root or node["selected"] or (
                node["new"]
            ) else "1.2"
            label_weight = (
                ' font-weight="700"'
                if is_root or node["selected"] or node["new"]
                else ""
            )
            parts.append(
                f'<g class="{" ".join(classes)}" '
                f'data-node-id="{node["id"]}" '
                f'data-rank="{_xml_escape(node["rank"])}" '
                f'data-name="{_xml_escape(node["name"])}">'
                f'<circle cx="{node["x"]:.2f}" cy="{node["y"]:.2f}" '
                f'r="{radius:.1f}" fill="{fill}" stroke="{stroke}" '
                f'stroke-width="{stroke_width}" />'
                f'<text class="node-label" x="{node["label_x"]:.2f}" '
                f'y="{node["label_y"]:.2f}" fill="{text_fill}" '
                'font-size="11" '
                "font-family=\"'Segoe UI', system-ui, sans-serif\""
                f'{label_weight}>'
                f'{_xml_escape(node["label"])}</text>'
            )
            for badge in node["badges"]:
                badge_fill = "#1d4ed8" if badge["kind"] == "selected" else (
                    "#15803d"
                )
                parts.append(
                    f'<text class="badge {badge["kind"]}" '
                    f'x="{badge["x"]:.2f}" y="{badge["y"]:.2f}" '
                    f'fill="{badge_fill}" font-size="10" font-weight="700" '
                    "font-family=\"'Segoe UI', system-ui, sans-serif\">"
                    f'{_xml_escape(badge["text"])}</text>'
                )
            parts.append("</g>")

    parts.append(
        f'<text class="footer" x="{margin:.2f}" y="{footer_y:.2f}" '
        'fill="#64748b" font-size="10" '
        "font-family=\"'Segoe UI', system-ui, sans-serif\">"
        "Taxonomía precalculada · sin red · un rango ausente no es "
        "evidencia negativa.</text>"
    )
    parts.append("</svg>")
    return "".join(parts)


def tree_html(
    tree: Mapping[str, Any],
    *,
    selected: str | None = None,
    new_genome: str | None = None,
    new_lineage: Any = None,
    title: str = "Árbol taxonómico",
    text: str | None = None,
    width: float = 1200.0,
    row_height: float = 22.0,
) -> str:
    """Documento HTML autónomo con el SVG y el árbol textual de respaldo.

    No incluye JavaScript ni recursos externos. El SVG se genera con
    :func:`tree_svg` y el respaldo usa :func:`tree_text`, de modo que si el
    navegador no muestra el gráfico la jerarquía sigue legible.
    """
    svg = tree_svg(
        tree,
        selected=selected,
        new_genome=new_genome,
        new_lineage=new_lineage,
        width=width,
        row_height=row_height,
        title=title,
    )
    if text is None:
        text = tree_text(
            tree,
            selected=selected,
            new_genome=new_genome,
            new_lineage=new_lineage,
        )
    return f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>{_xml_escape(title)}</title>
<style>
body {{ margin: 0; padding: 24px; background: #f8fafc; color: #0f172a;
       font: 14px 'Segoe UI', system-ui, sans-serif; }}
main {{ max-width: 1400px; margin: 0 auto; background: #ffffff;
        border: 1px solid #e2e8f0; border-radius: 12px;
        padding: 20px 24px; }}
h1 {{ font-size: 20px; margin: 0 0 8px; }}
p.limit {{ color: #475569; font-size: 12px; margin: 8px 0 16px; }}
.scroll {{ overflow: auto; max-height: 75vh; border: 1px solid #e2e8f0;
           border-radius: 8px; background: #ffffff; }}
details {{ margin-top: 16px; }}
summary {{ cursor: pointer; color: #1d4ed8; font-weight: 600; }}
pre {{ overflow: auto; max-height: 60vh; padding: 12px;
       background: #f8fafc; border: 1px solid #e2e8f0;
       border-radius: 8px; font-size: 12px; line-height: 1.45; }}
footer {{ margin-top: 18px; color: #64748b; font-size: 11px; }}
</style>
</head>
<body>
<main>
<h1>{_xml_escape(title)}</h1>
<p class="limit">{_xml_escape(LIMITATION.replace("**", ""))}</p>
<div class="scroll">{svg}</div>
<details>
<summary>Respaldo textual (misma jerarquía, sin gráfico)</summary>
<pre>{_xml_escape(text)}</pre>
</details>
<footer>Generado localmente en Python; no consulta la red y no recalcula
ANI.</footer>
</main>
</body>
</html>
"""


# ------------------------------------------------------------------- render
def _common_prefix(first: Sequence[str], second: Sequence[str]) -> list[str]:
    common: list[str] = []
    for left, right in zip(first, second):
        if left.casefold() != right.casefold():
            break
        common.append(left)
    return common


def zoom_html(svg: str, *, title: str = "Árbol taxonómico", height: int = 620) -> str:
    """Envuelve el SVG con zoom, reset y pantalla completa (sin dependencias)."""
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{title}</title>
<style>
body {{ margin:0; font-family: system-ui, sans-serif; background:#fff; }}
#bar {{ display:flex; gap:6px; align-items:center; padding:6px 8px; background:#f8fafc; border-bottom:1px solid #e2e8f0; position:sticky; top:0; z-index:2; }}
#bar button {{ cursor:pointer; border:1px solid #cbd5e1; background:#fff; border-radius:6px; padding:4px 9px; }}
#bar span {{ font-size:12px; color:#475569; }}
#viewport {{ overflow:auto; height:{height - 44}px; border:0; background:#fff; }}
#inner {{ transform-origin: 0 0; width:max-content; }}
</style></head><body>
<div id="bar">
  <button onclick="zoom(1.2)">+</button>
  <button onclick="zoom(1/1.2)">−</button>
  <button onclick="reset()">Ajustar</button>
  <button onclick="full()">⛶ Pantalla completa</button>
  <span id="pct">100%</span>
</div>
<div id="viewport"><div id="inner">{svg}</div></div>
<script>
let scale=1; const inner=document.getElementById('inner'), pct=document.getElementById('pct');
function apply(){{ inner.style.transform='scale('+scale+')'; pct.textContent=Math.round(scale*100)+'%'; }}
function zoom(f){{ scale=Math.max(0.2, Math.min(4, scale*f)); apply(); }}
function reset(){{ scale=1; apply(); document.getElementById('viewport').scrollTo(0,0); }}
function full(){{ const el=document.documentElement; if(el.requestFullscreen) el.requestFullscreen(); }}
document.getElementById('viewport').addEventListener('wheel', function(e){{ if(e.ctrlKey){{ e.preventDefault(); zoom(e.deltaY<0?1.12:1/1.12); }} }}, {{passive:false}});
apply();
</script>
</body></html>"""


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
            "graph": {
                "nodes": 0,
                "edges": 0,
                "leaves": 0,
                "width": 0.0,
                "height": 0.0,
                "svg_bytes": 0,
                "html_bytes": 0,
            },
            "svg": "",
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

    layout = tree_layout(
        tree,
        selected=str(selected),
        new_genome=new_genome_text,
        new_lineage=new_lineage,
    )
    svg = tree_svg(
        tree,
        selected=str(selected),
        new_genome=new_genome_text,
        new_lineage=new_lineage,
    )
    text_fallback = tree_text(
        tree,
        selected=str(selected),
        new_genome=new_genome_text,
        new_lineage=new_lineage,
    )
    html_doc = tree_html(
        tree,
        selected=str(selected),
        new_genome=new_genome_text,
        new_lineage=new_lineage,
        text=text_fallback,
    )

    st.markdown("#### Vista gráfica del árbol")
    st.caption(
        "Columnas = rango (Dominio → Especie). Azul: genoma seleccionado; "
        "verde discontinuo: rama nueva/externa. SVG generado en Python, "
        "sin dependencias extra ni red."
    )
    try:
        import streamlit.components.v1 as components
    except Exception:  # tests con un Streamlit mínimo
        components = None
    if components is not None:
        components.html(
            zoom_html(svg, title="Árbol taxonómico", height=620),
            height=640, scrolling=True,
        )
    else:
        st.markdown(
            '<div style="overflow:auto; max-height:72vh;">' + svg + '</div>',
            unsafe_allow_html=True,
        )
    st.download_button(
        "Descargar SVG",
        svg,
        file_name="arbol_taxonomico.svg",
        mime=SVG_MIME,
        key=f"{key}_download_svg",
    )
    st.download_button(
        "Descargar HTML (SVG + texto)",
        html_doc,
        file_name="arbol_taxonomico.html",
        mime=HTML_MIME,
        key=f"{key}_download_html",
    )
    st.caption(
        f"Gráfico: {layout['n_nodes']} nodos, {layout['n_leaves']} hojas; "
        "SVG y HTML se generan en memoria."
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
    st.caption("Respaldo textual del mismo árbol (funciona sin gráfico).")
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
        "graph": {
            "nodes": layout["n_nodes"],
            "edges": layout["n_edges"],
            "leaves": layout["n_leaves"],
            "width": layout["width"],
            "height": layout["height"],
            "svg_bytes": len(svg.encode("utf-8")),
            "html_bytes": len(html_doc.encode("utf-8")),
        },
        "svg": svg,
    }
