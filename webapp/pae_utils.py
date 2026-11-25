"""Lectura local de PAE (Predicted Aligned Error) de AlphaFold/AFDB.

El PAE es una matriz NxN de errores esperados (Å) entre pares de residuos.
Este módulo **solo lee** archivos locales ``*.pae.json``; no descarga nada,
no ejecuta modelos y no inventa valores cuando el archivo falta o es inválido.

Formatos aceptados
------------------
- AFDB: ``{"predicted_aligned_error": [[...]], "max_predicted_aligned_error": 31}``
- AlphaFold 3 / Server / ColabFold: ``{"pae": [[...]], "token_chain_ids": [...]}``
- Variantes con ``pae_matrix``, ``aligned_error`` o ``matrix``.
- Lista de modelos (ColabFold): se usa el primer modelo que traiga matriz PAE.

Convenciones
------------
- La matriz se devuelve como ``list[list[float]]`` con índices 0-based.
- Los resúmenes por residuo y por par usan numeración 1-based.
- ``summarize`` incluye la diagonal en ``mean``/``max``; el resumen por
  residuo excluye el residuo consigo mismo y ``pair_summary`` solo considera
  pares ordenados fuera de la diagonal.
- ``interface_pae`` es direccional: ``chain_a -> chain_b`` usa filas de A y
  columnas de B. Si el JSON no declara cadenas, no se infiere una interfaz.
"""
from __future__ import annotations

import json
import math
import os
import statistics
from pathlib import Path
from typing import Any, Iterable, Sequence

Matrix = list[list[float]]

# Claves de matriz vistas en AFDB, AlphaFold 3/Server, ColabFold y variantes.
PAE_KEYS: tuple[str, ...] = (
    "predicted_aligned_error",
    "pae",
    "pae_matrix",
    "aligned_error",
    "matrix",
)

# Claves de metadatos por residuo/token que permiten separar cadenas.
CHAIN_SERIES_KEYS: tuple[str, ...] = (
    "token_chain_ids",
    "chain_ids",
    "chain_index",
    "chain_id",
)

LIMITATION = (
    "El PAE es una expectativa de error relativo del modelo (Å), no una "
    "distancia ni evidencia experimental: no demuestra interacción, no "
    "describe desorden ni conformaciones alternativas y no sustituye pLDDT, "
    "ipTM ni validación de laboratorio."
)

NO_CHAIN_REASON = (
    "el archivo no declara cadenas (token_chain_ids/chain_ids/chains); "
    "un PAE de monómero AFDB no permite separar una interfaz"
)


# ---------------------------------------------------------------------------
# Validación básica
# ---------------------------------------------------------------------------
def _to_float(value: Any, *, where: str) -> float:
    """Convierte a float finito; traduce desbordamientos a ``ValueError``."""
    try:
        out = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{where}: número no representable ({value!r})") from exc
    if not math.isfinite(out):
        raise ValueError(f"{where}: PAE no finito ({value!r})")
    return out


def _number(value: Any, *, where: str) -> float:
    """Convierte un valor JSON a float, rechazando bool, NaN/Inf y negativos."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{where}: se esperaba un número, no {type(value).__name__}")
    out = _to_float(value, where=where)
    if out < 0:
        raise ValueError(f"{where}: PAE negativo ({out:g})")
    return out


def _as_matrix(value: Any, *, source: str = "PAE") -> Matrix:
    """Valida una matriz PAE NxN y la devuelve como listas de floats.

    Rechaza matrices vacías, filas de tamaño inconsistente, valores no
    numéricos, negativos o no finitos. Un PAE no es cuadrado si alguna fila
    no tiene exactamente N columnas.
    """
    if isinstance(value, dict):
        for key in PAE_KEYS:
            if key in value:
                return _as_matrix(value[key], source=key)
        raise ValueError(
            f"{source}: objeto sin matriz PAE; claves aceptadas: {', '.join(PAE_KEYS)}"
        )

    if not isinstance(value, (list, tuple)):
        raise ValueError(
            f"{source}: se esperaba una matriz (lista de listas), "
            f"no {type(value).__name__}"
        )
    n = len(value)
    if n == 0:
        raise ValueError(f"{source}: matriz PAE vacía")

    rows: Matrix = []
    for i, raw_row in enumerate(value, start=1):
        if not isinstance(raw_row, (list, tuple)):
            raise ValueError(f"{source}: fila {i}: se esperaba una lista de números")
        if len(raw_row) != n:
            raise ValueError(
                f"{source}: fila {i} tiene {len(raw_row)} columnas; se esperaban {n} "
                "(se exige una matriz NxN)"
            )
        row = [
            _number(cell, where=f"{source}: fila {i}, columna {j + 1}")
            for j, cell in enumerate(raw_row)
        ]
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Carga de archivos
# ---------------------------------------------------------------------------
def _parse_data(data: Any) -> tuple[Matrix, str, int | None, dict[str, Any] | None]:
    """Extrae la matriz y los metadatos de un JSON ya decodificado.

    Devuelve ``(matriz, clave_usada, indice_modelo, metadatos)``. Para una
    lista de modelos (ColabFold) se usa el primer modelo con matriz PAE y
    ``indice_modelo`` es su posición 0-based; en los demás casos es ``None``.
    """
    if isinstance(data, dict):
        for key in PAE_KEYS:
            if key in data:
                return _as_matrix(data[key], source=key), key, None, data
        raise ValueError(
            "JSON sin matriz PAE; claves aceptadas: " + ", ".join(PAE_KEYS)
        )

    if isinstance(data, (list, tuple)):
        # Lista de modelos: [{"pae": [...], "plddt": [...]}, ...]
        if data and isinstance(data[0], dict):
            for index, model in enumerate(data):
                if not isinstance(model, dict):
                    continue
                for key in PAE_KEYS:
                    if key in model:
                        return (
                            _as_matrix(model[key], source=key),
                            key,
                            index,
                            model,
                        )
            raise ValueError(
                "lista de modelos sin matriz PAE; claves aceptadas: "
                + ", ".join(PAE_KEYS)
            )
        # Matriz directa: [[...], ...]
        return _as_matrix(data, source="matriz"), "matriz", None, None

    raise ValueError(
        f"JSON PAE debe ser un objeto o una lista, no {type(data).__name__}"
    )


def _read_json(path: str | os.PathLike[str]) -> Any:
    """Lee y decodifica un JSON local con errores en español y sin red."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"PAE no encontrado: {p}")
    try:
        text = p.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError(f"PAE no es UTF-8: {exc}") from exc
    except OSError as exc:
        raise ValueError(f"PAE no legible: {exc}") from exc
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"JSON inválido: {exc.msg} (línea {exc.lineno}, columna {exc.colno})"
        ) from exc


def load_matrix(path: str | os.PathLike[str]) -> Matrix:
    """Carga y valida ``*.pae.json``; devuelve la matriz NxN de floats.

    Lanza ``FileNotFoundError`` si la ruta no existe y ``ValueError`` si el
    JSON es inválido, no contiene PAE o la matriz no es NxN.
    """
    data = _read_json(path)
    matrix, _source, _model_index, _meta = _parse_data(data)
    return matrix


def _optional_number(data: dict[str, Any] | None, key: str) -> float | None:
    if not isinstance(data, dict):
        return None
    value = data.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    try:
        value = float(value)
    except (OverflowError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _optional_series(data: dict[str, Any] | None, key: str, n: int) -> list[float] | None:
    if not isinstance(data, dict):
        return None
    value = data.get(key)
    if not isinstance(value, (list, tuple)) or len(value) != n:
        return None
    out: list[float] = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            return None
        try:
            item = float(item)
        except (OverflowError, ValueError):
            return None
        if not math.isfinite(item):
            return None
        out.append(item)
    return out


# ---------------------------------------------------------------------------
# Cadenas
# ---------------------------------------------------------------------------
def _chain_label(position: int) -> str:
    """Etiqueta tipo PDB para cadenas sin identificador: A..Z, AA..AZ, ..."""
    label = ""
    value = position
    while True:
        label = chr(ord("A") + value % 26) + label
        value = value // 26 - 1
        if value < 0:
            return label


def _chain_key(value: Any) -> str:
    if value is None:
        raise ValueError("identificador de cadena nulo")
    if isinstance(value, bool):
        raise ValueError("identificador de cadena booleano")
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not value.is_integer():
            raise ValueError(f"identificador de cadena no entero: {value!r}")
        return str(int(value))
    text = str(value).strip()
    if not text:
        raise ValueError("identificador de cadena vacío")
    return text


def _indices(spec: Any, n: int, *, chain: str) -> list[int]:
    """Normaliza índices 0-based de una cadena y valida rango/duplicados."""
    if isinstance(spec, dict):
        for key in ("indices", "residues", "positions"):
            if key in spec:
                return _indices(spec[key], n, chain=chain)
        if "start" in spec and "end" in spec:
            start, end = spec["start"], spec["end"]
            if (
                isinstance(start, bool)
                or isinstance(end, bool)
                or not isinstance(start, int)
                or not isinstance(end, int)
            ):
                raise ValueError(f"cadena {chain!r}: start/end deben ser enteros")
            return _indices(range(start, end + 1), n, chain=chain)
        raise ValueError(
            f"cadena {chain!r}: se esperaban indices/residues/positions o start/end"
        )
    if isinstance(spec, range):
        values: Iterable[Any] = spec
    elif isinstance(spec, (list, tuple)):
        values = spec
    elif isinstance(spec, int) and not isinstance(spec, bool):
        values = [spec]
    else:
        raise ValueError(
            f"cadena {chain!r}: se esperaba una lista de índices, "
            f"no {type(spec).__name__}"
        )

    out: list[int] = []
    for raw in values:
        if isinstance(raw, bool) or not isinstance(raw, int):
            raise ValueError(f"cadena {chain!r}: índice no entero ({raw!r})")
        if not 0 <= raw < n:
            raise ValueError(
                f"cadena {chain!r}: índice {raw} fuera de rango [0, {n - 1}]"
            )
        out.append(raw)
    if not out:
        raise ValueError(f"cadena {chain!r}: sin residuos")
    if len(set(out)) != len(out):
        raise ValueError(f"cadena {chain!r}: índices duplicados")
    return sorted(out)


def _chains_from_series(series: Any, n: int, name: str) -> tuple[dict[str, list[int]], list[str]]:
    """Agrupa una lista por residuo (p. ej. ``token_chain_ids``)."""
    if not isinstance(series, (list, tuple)):
        raise ValueError(f"{name} debe ser una lista de {n} identificadores")
    if len(series) != n:
        raise ValueError(
            f"{name} tiene {len(series)} valores; se esperaban {n} "
            "(uno por residuo/token)"
        )
    chains: dict[str, list[int]] = {}
    order: list[str] = []
    for index, raw in enumerate(series):
        key = _chain_key(raw)
        if key not in chains:
            chains[key] = []
            order.append(key)
        chains[key].append(index)
    return chains, order


def _chains_explicit(raw: Any, n: int) -> tuple[dict[str, list[int]], list[str]]:
    """Normaliza ``{"A": [0, 1], "B": [2, 3]}`` o una lista de grupos."""
    if isinstance(raw, dict):
        items: list[tuple[Any, Any]] = list(raw.items())
    elif isinstance(raw, (list, tuple)):
        items = []
        for position, item in enumerate(raw):
            if isinstance(item, dict):
                chain_id = item.get(
                    "chain_id", item.get("id", item.get("chain", _chain_label(position)))
                )
                spec = item.get("indices", item.get("residues", item.get("positions")))
                if spec is None and "start" in item and "end" in item:
                    spec = {"start": item["start"], "end": item["end"]}
                items.append((chain_id, spec))
            elif isinstance(item, (list, tuple, range)):
                items.append((_chain_label(position), item))
            else:
                raise ValueError(
                    "chains: cada grupo debe ser lista, rango u objeto, "
                    f"no {type(item).__name__}"
                )
    else:
        raise ValueError(
            f"chains debe ser un objeto o lista de grupos, no {type(raw).__name__}"
        )

    chains: dict[str, list[int]] = {}
    order: list[str] = []
    used: set[int] = set()
    for raw_id, spec in items:
        key = _chain_key(raw_id)
        if key in chains:
            raise ValueError(f"cadena duplicada: {key!r}")
        indices = _indices(spec, n, chain=key)
        overlap = used.intersection(indices)
        if overlap:
            raise ValueError(
                f"cadena {key!r}: índices ya asignados a otra cadena: "
                + ", ".join(str(i) for i in sorted(overlap)[:5])
            )
        used.update(indices)
        chains[key] = indices
        order.append(key)
    if not chains:
        raise ValueError("chains vacío")
    return chains, order


def _extract_chains(
    data: dict[str, Any] | None, n: int
) -> tuple[dict[str, list[int]], list[str], str]:
    """Devuelve ``(chains, orden, motivo)``; nunca lanza por metadatos."""
    if not isinstance(data, dict):
        return {}, [], NO_CHAIN_REASON

    raw = data.get("chains")
    if isinstance(raw, dict):
        try:
            chains, order = _chains_explicit(raw, n)
            return chains, order, ""
        except ValueError as exc:
            return {}, [], f"metadatos de cadena inválidos: {exc}"
    if raw is not None:
        if isinstance(raw, (list, tuple)) and raw:
            if any(isinstance(item, (list, tuple, dict, range)) for item in raw):
                try:
                    chains, order = _chains_explicit(raw, n)
                    return chains, order, ""
                except ValueError as exc:
                    return {}, [], f"metadatos de cadena inválidos: {exc}"
            try:
                chains, order = _chains_from_series(raw, n, "chains")
                return chains, order, ""
            except ValueError as exc:
                return {}, [], f"metadatos de cadena inválidos: {exc}"
        if not isinstance(raw, (list, tuple)):
            return {}, [], (
                "metadatos de cadena inválidos: 'chains' debe ser objeto o lista"
            )

    for key in CHAIN_SERIES_KEYS:
        if key not in data:
            continue
        value = data[key]
        if isinstance(value, dict):
            try:
                chains, order = _chains_explicit(value, n)
                return chains, order, ""
            except ValueError as exc:
                return {}, [], f"metadatos de cadena inválidos: {exc}"
        try:
            chains, order = _chains_from_series(value, n, key)
            return chains, order, ""
        except ValueError as exc:
            return {}, [], f"metadatos de cadena inválidos: {exc}"

    return {}, [], NO_CHAIN_REASON


def load_pae(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Carga un ``*.pae.json`` y devuelve matriz, estadísticas y metadatos.

    El resultado incluye ``matrix``, ``n``, ``mean``, ``max``, ``chains``,
    ``chain_order``, ``chain_reason`` y metadatos opcionales (``plddt``,
    ``ptm``, ``iptm``, ``max_predicted_aligned_error``). No lanza por
    metadatos de cadena ausentes o inválidos: devuelve ``chains={}`` y el
    motivo en ``chain_reason``.
    """
    data = _read_json(path)

    matrix, source_key, model_index, meta = _parse_data(data)
    n = len(matrix)
    chains, chain_order, chain_reason = _extract_chains(meta, n)
    plddt = _optional_series(meta, "plddt", n)
    if plddt is None:
        plddt = _optional_series(meta, "token_plddts", n)

    return {
        "matrix": matrix,
        "n": n,
        "mean": sum(sum(row) for row in matrix) / (n * n),
        "max": max(max(row) for row in matrix),
        "source_key": source_key,
        "model_index": model_index,
        "chains": chains,
        "chain_order": chain_order,
        "chain_reason": chain_reason,
        "plddt": plddt,
        "ptm": _optional_number(meta, "ptm"),
        "iptm": _optional_number(meta, "iptm"),
        "max_predicted_aligned_error": _optional_number(
            meta, "max_predicted_aligned_error"
        ),
        "limitation": LIMITATION,
    }


# ---------------------------------------------------------------------------
# Estadísticas
# ---------------------------------------------------------------------------
def mean_pae(matrix: Sequence[Sequence[float]]) -> float:
    """Media aritmética de toda la matriz (incluye la diagonal)."""
    m = _as_matrix(matrix)
    n = len(m)
    return sum(sum(row) for row in m) / (n * n)


def max_pae(matrix: Sequence[Sequence[float]]) -> float:
    """Máximo de toda la matriz (incluye la diagonal)."""
    m = _as_matrix(matrix)
    return max(max(row) for row in m)


def residue_summary(matrix: Sequence[Sequence[float]]) -> list[dict[str, Any]]:
    """Resumen por residuo excluyendo el par consigo mismo.

    Cada fila devuelve ``residue`` (1-based), ``mean_pae``, ``max_pae`` y
    ``min_pae`` sobre las columnas de los demás residuos. Para ``n == 1``
    los tres valores son ``None`` (no hay pares).
    """
    m = _as_matrix(matrix)
    n = len(m)
    out: list[dict[str, Any]] = []
    for i, row in enumerate(m):
        values = [row[j] for j in range(n) if j != i]
        out.append(
            {
                "residue": i + 1,
                "mean_pae": statistics.fmean(values) if values else None,
                "max_pae": max(values) if values else None,
                "min_pae": min(values) if values else None,
            }
        )
    return out


def pair_summary(matrix: Sequence[Sequence[float]]) -> dict[str, Any]:
    """Resumen de pares ordenados fuera de la diagonal.

    ``highest`` es el par con mayor PAE (peor posición relativa) y
    ``lowest`` el de menor PAE. Los índices ``i``/``j`` son 1-based y el
    orden importa: ``(i, j)`` significa "alinear en i, predecir j".
    """
    m = _as_matrix(matrix)
    n = len(m)
    pairs = [
        (i, j, m[i][j])
        for i in range(n)
        for j in range(n)
        if i != j
    ]
    if not pairs:
        return {
            "n": n,
            "n_pairs": 0,
            "mean_pae": None,
            "max_pae": None,
            "min_pae": None,
            "median_pae": None,
            "highest": None,
            "lowest": None,
        }
    values = [pae for _i, _j, pae in pairs]
    highest = max(pairs, key=lambda item: item[2])
    lowest = min(pairs, key=lambda item: item[2])
    return {
        "n": n,
        "n_pairs": len(pairs),
        "mean_pae": statistics.fmean(values),
        "max_pae": max(values),
        "min_pae": min(values),
        "median_pae": statistics.median(values),
        "highest": {"i": highest[0] + 1, "j": highest[1] + 1, "pae": highest[2]},
        "lowest": {"i": lowest[0] + 1, "j": lowest[1] + 1, "pae": lowest[2]},
    }


# ---------------------------------------------------------------------------
# Resumen tolerante a fallos
# ---------------------------------------------------------------------------
def _summary(
    *,
    available: bool,
    ok: bool,
    n: int = 0,
    mean: float | None = None,
    max_value: float | None = None,
    reason: str,
) -> dict[str, Any]:
    return {
        "available": available,
        "ok": ok,
        "n": n,
        "mean": mean,
        "max": max_value,
        "reason": reason,
        "limitation": LIMITATION,
    }


def summarize(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Resumen tolerante a fallos de un ``*.pae.json`` local.

    ``available`` indica que el archivo existe; ``ok`` indica que la matriz
    es válida y que ``n``/``mean``/``max`` son reales. Si falta, está
    corrupto o no es NxN devuelve ``n=0``, ``mean=None``, ``max=None`` y un
    ``reason`` explícito: nunca inventa PAE.
    """
    try:
        p = Path(path)
        exists = p.is_file()
    except (TypeError, OSError):
        p = None
        exists = False

    if not exists:
        where = f": {p}" if p is not None else ""
        return _summary(
            available=False,
            ok=False,
            reason=f"archivo PAE ausente (no encontrado){where}",
        )

    try:
        info = load_pae(p)  # type: ignore[arg-type]
    except (OSError, ValueError, TypeError) as exc:
        return _summary(
            available=True,
            ok=False,
            reason=f"PAE no utilizable: {exc}",
        )

    return _summary(
        available=True,
        ok=True,
        n=info["n"],
        mean=info["mean"],
        max_value=info["max"],
        reason="ok",
    )


# ---------------------------------------------------------------------------
# Interfaz entre cadenas
# ---------------------------------------------------------------------------
def _is_path(source: Any) -> bool:
    return isinstance(source, (str, os.PathLike))


def _resolve_chain(
    chains: dict[str, list[int]], chain: Any, order: list[str] | None = None
) -> str | None:
    """Resuelve un ID o un índice de orden a una clave real de ``chains``."""
    if chain is None or isinstance(chain, bool):
        return None
    if isinstance(chain, int):
        if str(chain) in chains:
            return str(chain)
        keys = order or list(chains)
        if 0 <= chain < len(keys):
            return keys[chain]
        return None
    if isinstance(chain, float):
        if not chain.is_integer():
            return None
        return _resolve_chain(chains, int(chain), order)
    if isinstance(chain, str):
        if chain in chains:
            return chain
        matches = [key for key in chains if key.lower() == chain.lower()]
        if len(matches) == 1:
            return matches[0]
        keys = order or list(chains)
        if chain.isdigit():
            index = int(chain)
            if 0 <= index < len(keys):
                return keys[index]
        return None
    return None


def _interface_base(chain_a: Any, chain_b: Any) -> dict[str, Any]:
    return {
        "available": False,
        "ok": False,
        "chain_a": chain_a,
        "chain_b": chain_b,
        "n_a": 0,
        "n_b": 0,
        "n_pairs": 0,
        "mean": None,
        "max": None,
        "median": None,
        "mean_pae": None,
        "max_pae": None,
        "reason": "",
        "limitation": LIMITATION,
    }


def interface_pae(
    source: str | os.PathLike[str] | Sequence[Sequence[float]],
    chain_a: Any,
    chain_b: Any,
    chains: dict[str, Any] | list[Any] | None = None,
) -> dict[str, Any]:
    """PAE medio/máximo de la interfaz ``chain_a -> chain_b``.

    ``source`` puede ser una ruta a ``*.pae.json`` o una matriz NxN. Las
    cadenas se toman de los metadatos del JSON (``token_chain_ids``,
    ``chain_ids`` o ``chains``) o del argumento ``chains``; si no hay
    metadatos, devuelve ``ok=False`` y explica el motivo en vez de suponer
    límites de cadena.

    Es direccional: usa filas de ``chain_a`` y columnas de ``chain_b``. Por
    eso se recomienda mirar ambas direcciones. ``chain_a``/``chain_b``
    aceptan el ID de cadena o su índice de orden 0-based.
    """
    result = _interface_base(chain_a, chain_b)
    available = True
    matrix: Matrix
    chain_map: dict[str, list[int]]
    chain_order: list[str]
    chain_reason: str

    if _is_path(source):
        p = Path(source)  # type: ignore[arg-type]
        if not p.is_file():
            result.update(
                available=False,
                reason=f"archivo PAE ausente (no encontrado): {p}",
            )
            return result
        try:
            info = load_pae(p)
        except (OSError, ValueError, TypeError) as exc:
            result.update(available=True, reason=f"PAE no utilizable: {exc}")
            return result
        matrix = info["matrix"]
        chain_map = info["chains"]
        chain_order = info["chain_order"]
        chain_reason = info["chain_reason"]
    else:
        try:
            matrix = _as_matrix(source)
        except (ValueError, TypeError) as exc:
            result.update(available=True, reason=f"PAE no utilizable: {exc}")
            return result
        chain_map, chain_order, chain_reason = {}, [], NO_CHAIN_REASON

    if chains is not None:
        try:
            if (
                isinstance(chains, (list, tuple))
                and chains
                and not any(
                    isinstance(item, (list, tuple, dict, range)) for item in chains
                )
            ):
                # Atajo: lista de IDs por residuo, p. ej. ["A", "A", "B"].
                chain_map, chain_order = _chains_from_series(
                    chains, len(matrix), "chains"
                )
            else:
                chain_map, chain_order = _chains_explicit(chains, len(matrix))
            chain_reason = ""
        except (ValueError, TypeError) as exc:
            result.update(
                available=available,
                reason=f"metadatos de cadena inválidos: {exc}",
            )
            return result

    if not chain_map:
        result.update(
            available=available,
            reason=chain_reason or NO_CHAIN_REASON,
        )
        return result

    resolved_a = _resolve_chain(chain_map, chain_a, chain_order)
    resolved_b = _resolve_chain(chain_map, chain_b, chain_order)
    if resolved_a is None:
        result.update(
            available=available,
            reason=(
                f"cadena {chain_a!r} no encontrada; disponibles: "
                + ", ".join(chain_order)
            ),
        )
        return result
    if resolved_b is None:
        result.update(
            available=available,
            reason=(
                f"cadena {chain_b!r} no encontrada; disponibles: "
                + ", ".join(chain_order)
            ),
        )
        return result
    if resolved_a == resolved_b:
        result.update(
            available=available,
            reason=(
                f"chain_a y chain_b son la misma cadena ({resolved_a!r}); "
                "no hay interfaz entre cadenas distintas"
            ),
        )
        return result

    rows = chain_map[resolved_a]
    cols = chain_map[resolved_b]
    values = [matrix[i][j] for i in rows for j in cols]
    if not values:
        result.update(
            available=available,
            reason="alguna cadena no tiene residuos para calcular la interfaz",
        )
        return result

    mean_value = statistics.fmean(values)
    max_value = max(values)
    result.update(
        available=available,
        ok=True,
        chain_a=resolved_a,
        chain_b=resolved_b,
        n_a=len(rows),
        n_b=len(cols),
        n_pairs=len(values),
        mean=mean_value,
        max=max_value,
        median=statistics.median(values),
        mean_pae=mean_value,
        max_pae=max_value,
        reason="ok",
    )
    return result


__all__ = [
    "LIMITATION",
    "NO_CHAIN_REASON",
    "PAE_KEYS",
    "interface_pae",
    "load_matrix",
    "load_pae",
    "max_pae",
    "mean_pae",
    "pair_summary",
    "residue_summary",
    "summarize",
]
