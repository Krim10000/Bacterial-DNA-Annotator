"""Ficha de PAE de interfaz para complejos, sin descargas ni inferencias.

Esta capa conecta :mod:`webapp.pae_complex` con Streamlit para que la ficha de
un complejo pueda mostrar su PAE **si ya existe en disco**. No descarga archivos,
no ejecuta AlphaFold y no inventa un mapeo cadena-subunidad.

Reglas deliberadas
------------------
- **PAE local o nada.** Solo se usan archivos ``*.pae.json`` locales. Una ruta
  remota se rechaza con ``not_available``; nunca se descarga.
- **El mapeo cadena-subunidad es una declaración.** Se acepta solo desde
  ``complex_row``, un sidecar junto al PAE o claves explícitas dentro del JSON
  del PAE. No se infiere por orden, por nombre ni por coincidencia de IDs.
- **Sin cadenas no hay interfaz.** Si el PAE no declara cadenas (o el mapeo
  falta), el estado es ``not_available``; no se rellena con pLDDT, ipTM ni
  estimaciones.
- **El PAE no demuestra interacción ni estequiometría.** La interfaz lo dice
  siempre; un valor bajo no es una prueba experimental.

API pública
-----------
``render_complex_pae(complex_row, data_root, key_prefix="")``
    Dibuja el bloque en Streamlit y devuelve el dict de la vista.
``complex_pae_view(complex_row, data_root)``
    Construye la vista sin tocar Streamlit (útil para pruebas e integraciones).
``find_complex_pae(complex_row, data_root)``
    Resuelve la ruta local del PAE, sin descargar.
``table_rows(view)``
    Filas listas para ``st.dataframe`` con media/máximo ``A -> B`` y ``B -> A``.

Ejemplo de integración en una ficha ya existente::

    from complex_pae_ui import render_complex_pae

    render_complex_pae(r.to_dict(), data_root, key_prefix=f"cx_{genome}")

``complex_row`` puede ser un ``dict`` o una fila de pandas (``Series``). Los
campos leídos son ``complex_id``/``id``, ``subunits`` y, si existe, ``pae`` (o
``pae_path``/``pae_json``/``pae_file``). El PAE se busca por ID de complejo y
por subunidades bajo ``data_root/metabolism/targets``.
"""
from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import streamlit as st

try:  # importación dentro del paquete webapp
    from . import pae_complex, pae_utils
except ImportError:  # pragma: no cover - ejecución plana desde webapp/
    import pae_complex
    import pae_utils


STATUS_OK = pae_complex.STATUS_OK
STATUS_NOT_AVAILABLE = pae_complex.STATUS_NOT_AVAILABLE
STATUS_INVALID = pae_complex.STATUS_INVALID
STATUSES: tuple[str, ...] = pae_complex.STATUSES

TITLE = "PAE de interfaz (predicho)"

#: Aviso breve que acompaña siempre al bloque.
CAVEAT = (
    "El PAE es error relativo del modelo (Å), no una distancia ni evidencia "
    "experimental. Un PAE de interfaz bajo no demuestra interacción física ni "
    "estequiometría; la asignación de cadenas a subunidades es una declaración "
    "y no se infiere del PAE."
)

LIMITATION = pae_complex.LIMITATION

#: Claves explícitas de mapeo cadena -> subunidad.
CHAIN_TO_SUBUNIT_KEYS: tuple[str, ...] = (
    "chain_subunits",
    "chain_to_subunit",
    "pae_chain_subunits",
)

#: Claves explícitas de mapeo subunidad -> cadena(s).
SUBUNIT_TO_CHAIN_KEYS: tuple[str, ...] = (
    "subunit_chains",
    "subunit_to_chains",
    "pae_subunit_chains",
)

#: Claves de ``complex_row`` que declaran una ruta local de PAE.
DECLARED_PATH_KEYS: tuple[str, ...] = (
    "pae",
    "pae_path",
    "pae_json",
    "pae_file",
)

MAPPING_SOURCES = ("complex_row", "sidecar", "pae_json")

_PAE_SUFFIX = ".pae.json"
_REMOTE_PREFIXES = ("http://", "https://", "ftp://", "sftp://", "s3://", "gs://")


# ---------------------------------------------------------------------------
# Normalización de la fila del complejo
# ---------------------------------------------------------------------------
def _as_mapping(row: Any) -> dict[str, Any]:
    """Convierte un dict o una fila pandas a un dict plano, sin lanzar."""
    if row is None:
        return {}
    if isinstance(row, Mapping):
        return dict(row)
    to_dict = getattr(row, "to_dict", None)
    if callable(to_dict):
        try:
            data = to_dict()
        except Exception:
            return {}
        return dict(data) if isinstance(data, Mapping) else {}
    return {}


def _text(value: Any) -> str:
    """Devuelve texto limpio; ignora None, NaN, booleanos y contenedores."""
    if value is None or isinstance(value, bool):
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, float) and math.isnan(value):
        return ""
    if isinstance(value, (list, tuple, dict, set)):
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "none", "<na>", "nat"}:
        return ""
    return text


def _has_content(value: Any) -> bool:
    """Indica si un valor de fila tiene contenido útil."""
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, tuple, dict, set)):
        return bool(value)
    if isinstance(value, float) and math.isnan(value):
        return False
    return _text(value) != ""


def complex_id_of(complex_row: Any) -> str:
    """Devuelve ``complex_id`` (o ``id``) de la fila, sin inventarlo."""
    row = _as_mapping(complex_row)
    for key in ("complex_id", "id"):
        text = _text(row.get(key))
        if text:
            return text
    return ""


def _split_parts(text: str) -> list[str]:
    """Separa subunidades por ``;`` o ``,``; conserva los dos puntos."""
    return [part.strip() for part in re.split(r"[;,]", text) if part.strip()]


def _unique(items: Sequence[str]) -> list[str]:
    """Deduplica preservando el orden y sin distinguir mayúsculas."""
    out: list[str] = []
    seen: set[str] = set()
    for item in items:
        text = item.strip()
        key = text.casefold()
        if text and key not in seen:
            seen.add(key)
            out.append(text)
    return out


def _entry_from_text(text: str) -> dict[str, Any]:
    """Normaliza una entrada ``gen:locus`` (o un ID suelto) para la ficha."""
    raw = text.strip()
    label = raw
    tail = ""
    if ":" in raw:
        head, rest = raw.split(":", 1)
        if head.strip():
            label = head.strip()
        tail = rest.strip()
    return {
        "label": label or raw,
        "raw": raw,
        "tokens": _unique([raw, label, tail]),
    }


def _entry_from_mapping(item: Mapping[str, Any]) -> dict[str, Any] | None:
    """Extrae una subunidad de un dict estructurado, si declara un ID."""
    fields = (
        "subunit",
        "subunit_id",
        "id",
        "name",
        "gene",
        "locus_tag",
        "uniprot",
        "acc",
        "accession",
    )
    values = [_text(item.get(key)) for key in fields]
    values = [value for value in values if value]
    if not values:
        return None
    label = values[0]
    return {
        "label": label,
        "raw": label,
        "tokens": _unique(values),
    }


def subunit_entries(complex_row: Any) -> list[dict[str, Any]]:
    """Devuelve las subunidades declaradas por la fila, sin inferirlas.

    Acepta el formato del proyecto (``"cyda:locus;cydb:locus"``), listas de
    IDs o dicts con campos como ``subunit``/``gene``/``uniprot``. No usa los
    valores de un dict como IDs de cadena: solo las claves declaradas.
    """
    row = _as_mapping(complex_row)
    raw: Any = None
    for key in ("subunits", "subunit_ids", "complex_subunits"):
        value = row.get(key)
        if _has_content(value):
            raw = value
            break
    if raw is None:
        return []

    entries: list[dict[str, Any]] = []
    if isinstance(raw, str):
        entries = [_entry_from_text(part) for part in _split_parts(raw)]
    elif isinstance(raw, Mapping):
        for key in raw:
            label = _text(key)
            if label:
                entries.append(_entry_from_text(label))
    elif isinstance(raw, (list, tuple)):
        for item in raw:
            if isinstance(item, Mapping):
                entry = _entry_from_mapping(item)
                if entry is not None:
                    entries.append(entry)
            elif isinstance(item, str):
                entries.extend(
                    _entry_from_text(part) for part in _split_parts(item)
                )
            elif item is not None and _text(item):
                entries.append(_entry_from_text(_text(item)))
    return entries


def complex_subunits(complex_row: Any) -> list[str]:
    """Etiquetas de subunidad para mostrar; desambigua duplicadas con el raw."""
    labels: list[str] = []
    for entry in subunit_entries(complex_row):
        label = entry["label"]
        if label in labels:
            fallback = entry["raw"]
            label = fallback if fallback not in labels else f"{label}#{len(labels) + 1}"
        labels.append(label)
    return labels


def _search_groups(complex_row: Any) -> list[list[str]]:
    """Grupos de tokens por subunidad para casar nombres de archivo."""
    groups: list[list[str]] = []
    for entry in subunit_entries(complex_row):
        tokens = _unique([token for token in entry["tokens"] if len(token) >= 2])
        if tokens:
            groups.append(tokens)
    return groups


def _sanitize(text: str) -> str:
    """Normaliza un ID para compararlo con un nombre de archivo."""
    cleaned = re.sub(r"[^0-9A-Za-z_.-]+", "_", text)
    return re.sub(r"_+", "_", cleaned).strip("_")


def _id_variants(complex_row: Any) -> set[str]:
    """Variantes razonables del ID de complejo, sin prefijos genéricos."""
    row = _as_mapping(complex_row)
    variants: set[str] = set()
    for key in ("complex_id", "id"):
        text = _text(row.get(key))
        if not text:
            continue
        variants.add(text)
        variants.add(_sanitize(text))
        variants.add(_sanitize(text.replace(":", "__")))
        variants.add(text.replace(":", "__"))
    return {variant for variant in variants if variant}


def _is_remote(value: str) -> bool:
    """Detecta una ruta remota para rechazarla; nunca se descarga."""
    return value.strip().casefold().startswith(_REMOTE_PREFIXES)


# ---------------------------------------------------------------------------
# Búsqueda local del PAE
# ---------------------------------------------------------------------------
def _declared_paths(complex_row: Any) -> tuple[list[str], Mapping[str, Any] | None]:
    """Lee rutas declaradas y, si ``pae`` es un dict, su mapeo embebido."""
    row = _as_mapping(complex_row)
    paths: list[str] = []
    mapping_container: Mapping[str, Any] | None = None
    for key in DECLARED_PATH_KEYS:
        if key not in row:
            continue
        value = row.get(key)
        if not _has_content(value):
            continue
        if isinstance(value, Mapping):
            if mapping_container is None:
                mapping_container = value
            for path_key in ("path", "file", "pae", "pae_path"):
                text = _text(value.get(path_key))
                if text:
                    paths.append(text)
            continue
        if isinstance(value, (list, tuple)):
            for item in value:
                if isinstance(item, Mapping):
                    if mapping_container is None:
                        mapping_container = item
                    for path_key in ("path", "file", "pae", "pae_path"):
                        text = _text(item.get(path_key))
                        if text:
                            paths.append(text)
                elif _text(item):
                    paths.append(_text(item))
            continue
        text = _text(value)
        if text:
            paths.append(text)
    return _unique(paths), mapping_container


def resolve_declared_path(
    value: str, data_root: Path | str | None
) -> list[Path]:
    """Resuelve una ruta declarada contra ``data_root`` y la raíz del proyecto.

    Devuelve candidatos ordenados; no comprueba existencia. Una ruta absoluta
    se usa tal cual. Las relativas se prueban primero bajo ``data_root`` (que
    en la app es ``data/``), luego bajo su padre (la raíz del proyecto) y
    finalmente bajo el directorio actual.
    """
    path = Path(value).expanduser()
    if path.is_absolute():
        return [path]

    roots: list[Path] = []
    if data_root:
        root = Path(data_root)
        roots.extend([root, root.parent])
    roots.append(Path.cwd())

    candidates: list[Path] = []
    for root in roots:
        candidate = root / path
        if candidate not in candidates:
            candidates.append(candidate)
    return candidates


def _targets_root(data_root: Path | str | None) -> Path | None:
    """Devuelve ``data/metabolism/targets``; acepta también la raíz del repo."""
    if not data_root:
        return None
    root = Path(data_root)
    for candidate in (
        root / "metabolism" / "targets",
        root / "data" / "metabolism" / "targets",
    ):
        if candidate.is_dir():
            return candidate
    return root / "metabolism" / "targets"


def _collect_pae_files(root: Path) -> list[Path]:
    """Lista PAE locales bajo ``root``; los errores de disco no rompen la UI."""
    try:
        files = list(root.rglob(f"*{_PAE_SUFFIX}"))
    except (OSError, ValueError):
        return []
    unique: dict[str, Path] = {}
    for path in files:
        try:
            unique[str(path)] = path
        except OSError:
            continue
    return sorted(unique.values(), key=lambda item: str(item).casefold())


def _match_score(
    path: Path, id_variants: set[str], groups: list[list[str]]
) -> int:
    """Puntúa un archivo: ID exacto vale más que todas las subunidades."""
    name = path.name
    stem = name[: -len(_PAE_SUFFIX)] if name.endswith(_PAE_SUFFIX) else name
    folded = stem.casefold()
    if any(variant.casefold() == folded for variant in id_variants):
        return 100
    if len(groups) >= 2 and all(
        any(token.casefold() in folded for token in group) for group in groups
    ):
        return 50
    return 0


def find_complex_pae(
    complex_row: Any, data_root: Path | str | None
) -> dict[str, Any]:
    """Busca el PAE local del complejo sin descargar nada.

    Orden: ruta declarada en la fila -> ``data/metabolism/targets``. En la
    búsqueda se casa el nombre del archivo con el ID del complejo (prioridad)
    o con **todas** las subunidades declaradas. Si hay más de un candidato con
    la misma prioridad, se devuelve ``invalid``: no se elige uno al azar.
    """
    row = _as_mapping(complex_row)
    declared, declared_mapping = _declared_paths(row)
    declared_missing: list[str] = []

    local_declared = [value for value in declared if not _is_remote(value)]
    remote_declared = [value for value in declared if _is_remote(value)]

    for value in local_declared:
        try:
            candidates = resolve_declared_path(value, data_root)
        except (OSError, TypeError, ValueError, RuntimeError):
            continue
        for candidate in candidates:
            try:
                if candidate.is_file():
                    return {
                        "status": STATUS_OK,
                        "ok": True,
                        "path": candidate,
                        "reason": "ok",
                        "declared": True,
                        "declared_missing": declared_missing,
                        "declared_mapping": declared_mapping,
                        "search_root": None,
                        "candidates": [str(candidate)],
                    }
            except (OSError, ValueError, TypeError):
                continue
            declared_missing.append(str(candidate))

    root = _targets_root(data_root)
    if root is None or not root.is_dir():
        if remote_declared:
            reason = (
                "la ruta PAE declarada es remota y esta capa no descarga "
                f"archivos: {remote_declared[0]}"
            )
        elif root is None:
            reason = "no hay data_root para buscar un PAE local"
        else:
            reason = f"no existe el directorio de PAE local: {root}"
        return {
            "status": STATUS_NOT_AVAILABLE,
            "ok": False,
            "path": None,
            "reason": reason,
            "declared": bool(declared),
            "declared_missing": declared_missing,
            "declared_mapping": declared_mapping,
            "search_root": str(root) if root is not None else None,
            "candidates": [],
        }

    variants = _id_variants(row)
    groups = _search_groups(row)
    scored: list[tuple[int, Path]] = []
    for path in _collect_pae_files(root):
        score = _match_score(path, variants, groups)
        if score:
            scored.append((score, path))

    if not scored:
        if remote_declared:
            reason = (
                "la ruta PAE declarada es remota y esta capa no descarga "
                f"archivos: {remote_declared[0]}"
            )
        elif not variants and not groups:
            reason = (
                "el complejo no declara ID ni subunidades para buscar un PAE "
                f"local en {root}; no se descarga"
            )
        else:
            reason = (
                f"no hay PAE local para este complejo en {root} "
                "(buscado por ID y por el conjunto de subunidades); "
                "esta ficha no descarga archivos"
            )
        if declared_missing:
            reason = (
                "ruta PAE declarada no encontrada: "
                + ", ".join(declared_missing[:3])
                + "; "
                + reason
            )
        return {
            "status": STATUS_NOT_AVAILABLE,
            "ok": False,
            "path": None,
            "reason": reason,
            "declared": bool(declared),
            "declared_missing": declared_missing,
            "declared_mapping": declared_mapping,
            "search_root": str(root),
            "candidates": [],
        }

    best = max(score for score, _path in scored)
    top = [path for score, path in scored if score == best]
    if len(top) > 1:
        names = ", ".join(path.name for path in top[:4])
        extra = f" (+{len(top) - 4})" if len(top) > 4 else ""
        return {
            "status": STATUS_INVALID,
            "ok": False,
            "path": None,
            "reason": (
                "varios PAE locales coinciden con el complejo y no se elige "
                f"uno al azar: {names}{extra}"
            ),
            "declared": bool(declared),
            "declared_missing": declared_missing,
            "declared_mapping": declared_mapping,
            "search_root": str(root),
            "candidates": [str(path) for path in top],
        }

    path = top[0]
    return {
        "status": STATUS_OK,
        "ok": True,
        "path": path,
        "reason": "ok",
        "declared": bool(declared),
        "declared_missing": declared_missing,
        "declared_mapping": declared_mapping,
        "search_root": str(root),
        "candidates": [str(path)],
    }


# ---------------------------------------------------------------------------
# Mapeo cadena <-> subunidad: solo declaraciones explícitas
# ---------------------------------------------------------------------------
def _mapping_from(
    container: Any, source: str
) -> dict[str, Any]:
    """Extrae un mapeo explícito de ``container``; no inventa uno."""
    empty = {
        "mapping": None,
        "source": source,
        "reason": None,
        "present": False,
    }
    if not isinstance(container, Mapping):
        return empty

    for key in CHAIN_TO_SUBUNIT_KEYS:
        if key not in container:
            continue
        value = container.get(key)
        if not _has_content(value):
            return {
                "mapping": None,
                "source": source,
                "reason": f"el mapeo cadena-subunidad ({key}) está vacío",
                "present": True,
            }
        try:
            mapping = pae_complex.subunits_from_chain_subunits(value)
        except (TypeError, ValueError) as exc:
            return {
                "mapping": None,
                "source": source,
                "reason": f"mapeo cadena-subunidad inválido ({key}): {exc}",
                "present": True,
            }
        return {
            "mapping": mapping,
            "source": source,
            "reason": None,
            "present": True,
        }

    for key in SUBUNIT_TO_CHAIN_KEYS:
        if key not in container:
            continue
        value = container.get(key)
        if not isinstance(value, Mapping) or not value:
            return {
                "mapping": None,
                "source": source,
                "reason": (
                    f"el mapeo subunidad-cadena ({key}) debe ser un objeto no "
                    "vacío"
                ),
                "present": True,
            }
        return {
            "mapping": {str(label): chains for label, chains in value.items()},
            "source": source,
            "reason": None,
            "present": True,
        }

    return empty


def _read_json(path: Path) -> Any:
    """Lee un JSON local; lanza ``ValueError`` con contexto si falla."""
    try:
        with path.open(encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise ValueError(f"no se pudo leer {path.name}: {exc}") from exc


def _sidecar_candidates(pae_path: Path) -> list[Path]:
    """Nombres de sidecar de mapeo admitidos junto al PAE."""
    name = pae_path.name
    if name.endswith(_PAE_SUFFIX):
        base = name[: -len(_PAE_SUFFIX)]
    elif name.endswith(".json"):
        base = name[: -len(".json")]
    else:
        base = name
    names = [
        f"{base}.chains.json",
        f"{base}.chain_subunits.json",
        f"{base}.subunit_chains.json",
        f"{base}.map.json",
        f"{name}.chains.json",
        f"{name}.map.json",
    ]
    out: list[Path] = []
    for candidate_name in names:
        candidate = pae_path.with_name(candidate_name)
        if candidate not in out:
            out.append(candidate)
    return out


def _json_containers(data: Any) -> list[Mapping[str, Any]]:
    """Devuelve el dict del PAE (o el primer modelo si es una lista)."""
    if isinstance(data, Mapping):
        return [data]
    if isinstance(data, (list, tuple)):
        for item in data:
            if isinstance(item, Mapping):
                return [item]
    return []


def resolve_chain_mapping(
    complex_row: Any,
    declared_mapping: Mapping[str, Any] | None,
    pae_path: Path,
) -> dict[str, Any]:
    """Resuelve el mapeo declarado: fila -> sidecar -> JSON del PAE.

    Nunca deriva el mapeo del orden de las cadenas, del nombre de la subunidad
    ni de que un ID de cadena coincida con una etiqueta. Si una fuente declara
    un mapeo inválido, se detiene y lo informa; no prueba la siguiente.
    """
    row = _as_mapping(complex_row)
    for container, source in (
        (declared_mapping, "complex_row"),
        (row, "complex_row"),
    ):
        result = _mapping_from(container, source)
        if result["present"]:
            return result

    # ``subunits`` como {subunidad: cadena(s)} es un mapeo explícito, no una
    # inferencia. Se acepta solo si es un objeto; una cadena o lista sigue
    # siendo solo la lista de subunidades.
    for key in ("subunits", "complex_subunits"):
        value = row.get(key)
        if isinstance(value, Mapping) and value:
            result = _mapping_from(
                {"subunit_chains": value}, "complex_row"
            )
            if result["present"]:
                return result

    for sidecar in _sidecar_candidates(pae_path):
        try:
            if not sidecar.is_file():
                continue
        except (OSError, ValueError, TypeError):
            continue
        try:
            data = _read_json(sidecar)
        except ValueError as exc:
            return {
                "mapping": None,
                "source": "sidecar",
                "reason": f"sidecar de mapeo ilegible: {exc}",
                "present": True,
            }
        result = _mapping_from(data, "sidecar")
        if result["present"]:
            return result

    try:
        data = _read_json(pae_path)
    except ValueError:
        # El lector de PAE dará el estado invalid con su propio motivo.
        return {
            "mapping": None,
            "source": None,
            "reason": None,
            "present": False,
        }
    for container in _json_containers(data):
        result = _mapping_from(container, "pae_json")
        if result["present"]:
            return result

    return {
        "mapping": None,
        "source": None,
        "reason": None,
        "present": False,
    }


def _declared_chains(pae_path: Path) -> tuple[list[str], str, str | None]:
    """Lee las cadenas declaradas por el PAE para explicar qué falta."""
    try:
        info = pae_utils.load_pae(pae_path)
    except Exception as exc:  # una matriz rara no debe romper la ficha
        return [], "", f"PAE no utilizable: {exc}"
    chains = [str(chain) for chain in (info.get("chains") or {})]
    chain_reason = str(info.get("chain_reason") or "")
    return chains, chain_reason, None


def _missing_mapping_reason(chains: list[str], chain_reason: str) -> str:
    """Mensaje accionable cuando hay PAE pero no mapeo cadena-subunidad."""
    base = (
        "falta el mapeo cadena-subunidad: hay que declarar qué cadena del PAE "
        "corresponde a cada subunidad en 'chain_subunits' o 'subunit_chains' "
        "(o en un sidecar *.chains.json); no se infiere por orden, nombre ni "
        "igualdad de ID"
    )
    if chains:
        return (
            base
            + f". El PAE declara las cadenas {', '.join(chains)}, pero no a qué "
            "subunidad pertenece cada una"
        )
    if chain_reason:
        return base + f". El PAE no permite leer cadenas: {chain_reason}"
    return base + ". El PAE no declara cadenas"


# ---------------------------------------------------------------------------
# Vista y tabla
# ---------------------------------------------------------------------------
def _empty_summary() -> dict[str, Any]:
    """Resumen vacío con la misma forma que ``pae_complex``."""
    return {
        "n_pairs": 0,
        "n_ok": 0,
        "n_not_available": 0,
        "n_invalid": 0,
        "n_residue_pairs": 0,
        "mean_pae": None,
        "max_pae": None,
    }


def _empty_view(
    *,
    status: str,
    reason: str,
    complex_id: str,
    subunits: list[str],
    pae_path: Path | None = None,
    pae_declared: bool = False,
    search_root: str | None = None,
    declared_missing: Sequence[str] = (),
    declared_chains: Sequence[str] = (),
) -> dict[str, Any]:
    """Vista base coherente para todos los estados."""
    return {
        "status": status,
        "ok": status == STATUS_OK,
        "reason": reason,
        "complex_id": complex_id,
        "subunits": list(subunits),
        "pae_path": str(pae_path) if pae_path is not None else None,
        "pae_found": pae_path is not None,
        "pae_declared": bool(pae_declared),
        "search_root": search_root,
        "declared_missing": list(declared_missing),
        "declared_chains": list(declared_chains),
        "mapping": None,
        "mapping_source": None,
        "mapping_reason": None,
        "analysis": None,
        "rows": [],
        "summary": _empty_summary(),
        "limitation": LIMITATION,
        "caveat": CAVEAT,
    }


def _pair_rows(analysis: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Una fila por par de subunidades, con ambas direcciones."""
    rows: list[dict[str, Any]] = []
    if not analysis:
        return rows
    for pair in analysis.get("pairs") or []:
        a_to_b = pair.get("a_to_b") or {}
        b_to_a = pair.get("b_to_a") or {}
        rows.append(
            {
                "pair": f"{pair.get('subunit_a')} / {pair.get('subunit_b')}",
                "subunit_a": pair.get("subunit_a"),
                "subunit_b": pair.get("subunit_b"),
                "chains_a": list(a_to_b.get("chains_a") or []),
                "chains_b": list(a_to_b.get("chains_b") or []),
                "mean_a_to_b": a_to_b.get("mean_pae"),
                "max_a_to_b": a_to_b.get("max_pae"),
                "mean_b_to_a": b_to_a.get("mean_pae"),
                "max_b_to_a": b_to_a.get("max_pae"),
                "n_pairs": int(pair.get("n_pairs") or 0),
                "status": pair.get("status"),
                "reason": pair.get("reason"),
            }
        )
    return rows


def complex_pae_view(
    complex_row: Any, data_root: Path | str | None
) -> dict[str, Any]:
    """Construye la vista del PAE de un complejo sin llamar a Streamlit.

    Devuelve siempre un dict con ``status``, ``reason``, ``pae_path``,
    ``mapping``, ``analysis``, ``rows``, ``summary`` y ``limitation``. No lanza
    por datos ausentes o inválidos.
    """
    row = _as_mapping(complex_row)
    cid = complex_id_of(row)
    subunits = complex_subunits(row)
    found = find_complex_pae(row, data_root)

    base = _empty_view(
        status=found["status"],
        reason=found["reason"],
        complex_id=cid,
        subunits=subunits,
        pae_path=found.get("path"),
        pae_declared=bool(found.get("declared")),
        search_root=found.get("search_root"),
        declared_missing=found.get("declared_missing") or (),
    )
    if found["status"] != STATUS_OK:
        return base

    pae_path = found["path"]
    mapping_result = resolve_chain_mapping(
        row, found.get("declared_mapping"), pae_path
    )
    base["mapping"] = mapping_result.get("mapping")
    base["mapping_source"] = mapping_result.get("source")
    base["mapping_reason"] = mapping_result.get("reason")

    if mapping_result.get("present") and mapping_result.get("mapping") is None:
        base["status"] = STATUS_INVALID
        base["ok"] = False
        base["reason"] = (
            "mapeo cadena-subunidad inválido: "
            + str(mapping_result.get("reason") or "sin detalle")
        )
        return base

    mapping = mapping_result.get("mapping")
    if mapping is None:
        chains, chain_reason, load_error = _declared_chains(pae_path)
        base["declared_chains"] = chains
        if load_error:
            base["status"] = STATUS_INVALID
            base["ok"] = False
            base["reason"] = load_error
            return base
        base["status"] = STATUS_NOT_AVAILABLE
        base["ok"] = False
        base["reason"] = _missing_mapping_reason(chains, chain_reason)
        return base

    try:
        analysis = pae_complex.analyze_complex(str(pae_path), mapping)
    except Exception as exc:  # la ficha no debe romperse por una matriz rara
        base["status"] = STATUS_INVALID
        base["ok"] = False
        base["reason"] = (
            f"error al analizar el PAE del complejo: {type(exc).__name__}: {exc}"
        )
        return base

    base["analysis"] = analysis
    base["rows"] = _pair_rows(analysis)
    base["summary"] = analysis.get("summary") or _empty_summary()
    base["status"] = analysis.get("status", STATUS_INVALID)
    base["ok"] = base["status"] == STATUS_OK
    base["reason"] = str(analysis.get("reason") or base["status"])
    return base


def _fmt(value: Any) -> str:
    """Formatea Å con un decimal; ``None`` se muestra como raya."""
    if value is None:
        return "—"
    try:
        return f"{float(value):.1f}"
    except (TypeError, ValueError):
        return "—"


def table_rows(view: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Convierte ``view['rows']`` en filas legibles para ``st.dataframe``.

    Columnas: par, subunidades, cadenas, media/máximo ``A -> B`` y
    ``B -> A``, pares de residuos, estado y motivo.
    """
    rows: list[dict[str, Any]] = []
    for row in view.get("rows") or []:
        chains_a = ", ".join(row.get("chains_a") or []) or "—"
        chains_b = ", ".join(row.get("chains_b") or []) or "—"
        status = str(row.get("status") or "")
        rows.append(
            {
                "Par": row.get("pair"),
                "Subunidad A": row.get("subunit_a"),
                "Cadenas A": chains_a,
                "Subunidad B": row.get("subunit_b"),
                "Cadenas B": chains_b,
                "A→B media (Å)": _fmt(row.get("mean_a_to_b")),
                "A→B máx (Å)": _fmt(row.get("max_a_to_b")),
                "B→A media (Å)": _fmt(row.get("mean_b_to_a")),
                "B→A máx (Å)": _fmt(row.get("max_b_to_a")),
                "Pares (residuos)": row.get("n_pairs") or 0,
                "Estado": status,
                "Motivo": "" if status == STATUS_OK else row.get("reason"),
            }
        )
    return rows


def _has_numbers(view: Mapping[str, Any]) -> bool:
    """Indica si la vista tiene al menos un valor numérico que mostrar."""
    for row in view.get("rows") or []:
        if any(
            row.get(key) is not None
            for key in (
                "mean_a_to_b",
                "max_a_to_b",
                "mean_b_to_a",
                "max_b_to_a",
            )
        ):
            return True
    return False


def _summary_caption(view: Mapping[str, Any]) -> str:
    """Resumen en una línea, sin ocultar pares no evaluados."""
    summary = view.get("summary") or _empty_summary()
    n_pairs = int(summary.get("n_pairs") or 0)
    n_ok = int(summary.get("n_ok") or 0)
    if not n_ok:
        return "Sin pares de subunidades evaluables."
    text = (
        f"Pares evaluados: {n_ok}/{n_pairs} · media del complejo "
        f"{_fmt(summary.get('mean_pae'))} Å · máximo "
        f"{_fmt(summary.get('max_pae'))} Å · "
        f"{int(summary.get('n_residue_pairs') or 0)} pares de residuos."
    )
    others: list[str] = []
    if summary.get("n_not_available"):
        others.append(f"{summary['n_not_available']} no disponibles")
    if summary.get("n_invalid"):
        others.append(f"{summary['n_invalid']} inválidos")
    if others:
        text += " Sin evaluar: " + ", ".join(others) + "."
    return text


def _source_caption(view: Mapping[str, Any]) -> str:
    """Fuente local y origen del mapeo, para poder auditar la tabla."""
    path = view.get("pae_path") or ""
    name = Path(path).name if path else "—"
    source = view.get("mapping_source") or "—"
    return f"Fuente local: {name} · mapeo cadena-subunidad: {source}."


def _widget_key(key_prefix: str, view: Mapping[str, Any]) -> str:
    """Clave estable para el dataframe, namespaced por ``key_prefix``."""
    slug = view.get("complex_id") or view.get("pae_path") or "complex"
    slug = re.sub(r"[^0-9A-Za-z_-]+", "_", str(slug)).strip("_") or "complex"
    return f"{key_prefix}complex_pae_{slug}"


# ---------------------------------------------------------------------------
# Render Streamlit
# ---------------------------------------------------------------------------
def render_complex_pae(
    complex_row: Any,
    data_root: Path | str | None,
    key_prefix: str = "",
) -> dict[str, Any]:
    """Dibuja el PAE de un complejo en la ficha y devuelve la vista.

    Parámetros
    ----------
    complex_row:
        Fila del complejo (dict o ``pandas.Series``). Debe declarar
        ``complex_id``/``id`` o ``subunits`` para poder buscar el PAE.
    data_root:
        ``data/`` del proyecto. La búsqueda ocurre en
        ``data/metabolism/targets``; no se descarga nada.
    key_prefix:
        Prefijo para la clave del dataframe cuando la ficha se dibuja varias
        veces.

    Devuelve
    --------
    dict
        La misma vista de :func:`complex_pae_view`.
    """
    view = complex_pae_view(complex_row, data_root)
    status = view["status"]

    if status == STATUS_OK:
        st.markdown(f"**{TITLE}**")
        st.caption(_source_caption(view))
        st.caption(_summary_caption(view))
        rows = table_rows(view)
        if rows:
            st.dataframe(
                rows,
                hide_index=True,
                width="stretch",
                key=_widget_key(key_prefix, view),
            )
        st.caption(CAVEAT)
        return view

    if status == STATUS_NOT_AVAILABLE:
        st.info(f"**{TITLE}:** `not_available` — {view['reason']}")
    else:
        st.warning(f"**{TITLE}:** `invalid` — {view['reason']}")
        if _has_numbers(view):
            st.caption(
                "Se muestran solo los pares con valores calculables; el estado "
                f"global del complejo es `{status}`."
            )
            st.dataframe(
                table_rows(view),
                hide_index=True,
                width="stretch",
                key=_widget_key(key_prefix, view),
            )
    st.caption(CAVEAT)
    return view


__all__ = [
    "CAVEAT",
    "CHAIN_TO_SUBUNIT_KEYS",
    "DECLARED_PATH_KEYS",
    "LIMITATION",
    "MAPPING_SOURCES",
    "STATUSES",
    "STATUS_INVALID",
    "STATUS_NOT_AVAILABLE",
    "STATUS_OK",
    "SUBUNIT_TO_CHAIN_KEYS",
    "TITLE",
    "complex_id_of",
    "complex_pae_view",
    "complex_subunits",
    "find_complex_pae",
    "render_complex_pae",
    "resolve_chain_mapping",
    "resolve_declared_path",
    "subunit_entries",
    "table_rows",
]
