"""PAE de interfaz para complejos, sobre ``webapp.pae_utils``.

Este módulo **no lee ni inventa PAE por su cuenta**: delega la lectura y el
cálculo de cada bloque de interfaz en :mod:`webapp.pae_utils`. Añade una capa
para complejos:

- normaliza un mapeo ``subunidad -> cadenas`` (una subunidad puede tener varias
  cadenas, por ejemplo un homómero);
- calcula la interfaz **direccional** ``A -> B`` y ``B -> A`` para cada par de
  subunidades;
- agrega media y máximo por par y para todo el complejo;
- declara estados explícitos ``ok``, ``not_available`` e ``invalid``.

Semántica de estados
--------------------
``ok``
    Hay PAE, hay cadenas y todos los pares de cadenas necesarios pudieron
    evaluarse. ``mean_pae``/``max_pae`` son reales.
``not_available``
    Falta el PAE, faltan las cadenas (el JSON no las declara y no se aportó un
    mapeo de residuos) o falta el mapeo de subunidades. No se rellena ningún
    valor.
``invalid``
    El PAE existe pero no es utilizable, o el mapeo de subunidades/cadenas es
    contradictorio (cadena desconocida, repetida, índices solapados, etc.).

Precedencia: un mapeo de subunidades mal formado es ``invalid`` aunque el PAE
falte. Un PAE ausente con un mapeo válido es ``not_available``.

Reglas deliberadas
------------------
- **Nunca se inventan cadenas.** Un ID de cadena debe existir en los metadatos
  del PAE o en el argumento ``chains``. No se acepta la resolución por índice
  ordinal de ``pae_utils`` como sustituto de un ID: si una cadena numérica
  termina resolviendo a otra cadena, el par se marca ``invalid``.
- **Medias conservadoras.** Si una subunidad tiene varias cadenas, solo se
  publica la media/máximo de una dirección cuando *todos* los pares de cadenas
  de esa dirección son ``ok``. Los detalles parciales quedan en
  ``chain_pairs``, pero no se presenta una media incompleta como si fuera la
  interfaz completa.
- **La media combinada está ponderada** por el número de pares de residuos de
  cada dirección; el máximo es el máximo de las direcciones.
- **El PAE no demuestra interacción ni estequiometría.** El mapeo de cadenas a
  subunidades es una declaración del usuario, no una inferencia del PAE.

Ejemplo
-------
>>> from webapp.pae_complex import analyze_complex
>>> result = analyze_complex(
...     "complejo.pae.json",
...     {"RpoA": ["A"], "RpoB": ["B"], "RpoC": ["C"]},
... )
>>> result["status"] in {"ok", "not_available", "invalid"}
True
"""
from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from itertools import combinations
from pathlib import Path
from typing import Any

try:  # importación plana cuando la app corre desde webapp/
    from . import pae_utils
except ImportError:  # pragma: no cover - paquete plano
    import pae_utils


STATUS_OK = "ok"
STATUS_NOT_AVAILABLE = "not_available"
STATUS_INVALID = "invalid"
STATUSES: tuple[str, ...] = (STATUS_OK, STATUS_NOT_AVAILABLE, STATUS_INVALID)

LIMITATION = (
    pae_utils.LIMITATION
    + " En un complejo, un PAE de interfaz bajo no demuestra interacción "
      "física, afinidad ni función, y no prueba la estequiometría ni la "
      "composición del ensamblado. La asignación de cadenas a subunidades es "
      "una declaración del usuario, no una inferencia del PAE, y el PAE no "
      "identifica qué residuos forman la interfaz."
)

NO_MAPPING_REASON = (
    "falta el mapeo de subunidades (subunidad -> cadenas); el PAE no permite "
    "inferir qué cadena corresponde a cada subunidad ni la estequiometría"
)

NO_CHAIN_REASON = pae_utils.NO_CHAIN_REASON

SubunitMap = dict[str, list[str]]


# ---------------------------------------------------------------------------
# Normalización del mapeo de subunidades
# ---------------------------------------------------------------------------
def _subunit_label(value: Any, *, where: str) -> str:
    """Valida una etiqueta de subunidad no vacía y la devuelve sin espacios."""
    if not isinstance(value, str):
        raise ValueError(
            f"{where}: la subunidad debe ser texto, no {type(value).__name__}"
        )
    text = value.strip()
    if not text:
        raise ValueError(f"{where}: etiqueta de subunidad vacía")
    return text


def _chain_id(value: Any, *, where: str) -> str:
    """Normaliza un ID de cadena a texto, sin aceptar booleanos ni vacíos."""
    if isinstance(value, bool):
        raise ValueError(f"{where}: identificador de cadena booleano")
    if isinstance(value, int):
        value = str(value)
    if not isinstance(value, str):
        raise ValueError(
            f"{where}: se esperaba un ID de cadena, no {type(value).__name__}"
        )
    text = value.strip()
    if not text:
        raise ValueError(f"{where}: identificador de cadena vacío")
    return text


def _normalize_subunits(
    subunits: Any, *, minimum: int = 2
) -> SubunitMap:
    """Valida ``{subunidad: cadena | [cadenas]}`` y lo normaliza.

    Rechaza etiquetas vacías, listas vacías, cadenas repetidas (dentro de una
    subunidad o entre subunidades), tipos no admitidos y mapas con menos de
    ``minimum`` subunidades. No modifica el diccionario original.
    """
    if not isinstance(subunits, Mapping):
        raise ValueError(
            "el mapeo de subunidades debe ser un objeto {subunidad: cadenas}, "
            f"no {type(subunits).__name__}"
        )

    normalized: SubunitMap = {}
    owner: dict[str, str] = {}

    for raw_label, raw_chains in subunits.items():
        label = _subunit_label(raw_label, where="mapeo de subunidades")
        if label in normalized:
            raise ValueError(f"subunidad duplicada: {label!r}")

        if isinstance(raw_chains, str) or (
            isinstance(raw_chains, int) and not isinstance(raw_chains, bool)
        ):
            candidates: list[Any] = [raw_chains]
        elif isinstance(raw_chains, (list, tuple)):
            candidates = list(raw_chains)
        else:
            raise ValueError(
                f"subunidad {label!r}: se esperaba una cadena o una lista de "
                f"cadenas, no {type(raw_chains).__name__}"
            )
        if not candidates:
            raise ValueError(f"subunidad {label!r}: sin cadenas")

        chains: list[str] = []
        for position, raw_chain in enumerate(candidates, start=1):
            chain = _chain_id(
                raw_chain, where=f"subunidad {label!r}, cadena {position}"
            )
            if chain in chains:
                raise ValueError(
                    f"subunidad {label!r}: cadena duplicada {chain!r}"
                )
            if chain in owner:
                raise ValueError(
                    f"cadena {chain!r} asignada a dos subunidades: "
                    f"{owner[chain]!r} y {label!r}"
                )
            chains.append(chain)
            owner[chain] = label
        normalized[label] = chains

    if len(normalized) < minimum:
        raise ValueError(
            f"el mapeo debe declarar al menos {minimum} subunidades; "
            f"recibidas {len(normalized)}"
        )
    return normalized


def subunits_from_chain_subunits(chain_to_subunit: Mapping[Any, Any]) -> SubunitMap:
    """Invierte ``{cadena: subunidad}`` a ``{subunidad: [cadenas]}``.

    Es útil cuando el mapeo disponible viene por cadena. Varias cadenas pueden
    apuntar a la misma subunidad (homómero). No analiza PAE ni exige que el
    resultado tenga dos subunidades; eso lo decide :func:`analyze_complex`.

    Lanza ``ValueError`` si el argumento no es un objeto, si una etiqueta de
    subunidad está vacía, si un ID de cadena se repite tras normalizarlo o si
    el mapeo está vacío.
    """
    if not isinstance(chain_to_subunit, Mapping):
        raise ValueError(
            "el mapeo cadena -> subunidad debe ser un objeto, "
            f"no {type(chain_to_subunit).__name__}"
        )

    grouped: SubunitMap = {}
    seen: set[str] = set()
    for raw_chain, raw_subunit in chain_to_subunit.items():
        chain = _chain_id(raw_chain, where="mapeo cadena -> subunidad")
        subunit = _subunit_label(
            raw_subunit, where=f"cadena {chain!r}"
        )
        if chain in seen:
            raise ValueError(f"cadena duplicada en el mapeo: {chain!r}")
        seen.add(chain)
        grouped.setdefault(subunit, []).append(chain)

    if not grouped:
        raise ValueError("mapeo cadena -> subunidad vacío")
    return grouped


# ---------------------------------------------------------------------------
# Carga y validación de la fuente PAE
# ---------------------------------------------------------------------------
def _is_path(source: Any) -> bool:
    return isinstance(source, (str, os.PathLike))


def _source_label(source: Any) -> str:
    return str(source) if _is_path(source) else "matriz"


def _source_kind(source: Any) -> str:
    return "file" if _is_path(source) else "matrix"


def _source_available(source: Any) -> bool:
    """Indica si la fuente existe (archivo) o si se aportó una matriz."""
    if not _is_path(source):
        return True
    try:
        return Path(source).is_file()
    except (TypeError, OSError):
        return False


def _invalid_info(reason: str, *, source: Any, available: bool) -> dict[str, Any]:
    return {
        "available": available,
        "status": STATUS_INVALID,
        "reason": reason,
        "matrix": None,
        "chains_arg": None,
        "has_chains": False,
        "chain_reason": "",
        "source": _source_label(source),
        "source_kind": _source_kind(source),
    }


def _load_source(source: Any, chains: Any) -> dict[str, Any]:
    """Carga el PAE una sola vez y deja lista la matriz y el mapeo de cadenas.

    Para archivos usa ``pae_utils.load_pae`` (valida JSON y matriz). Para
    matrices directas valida con ``pae_utils.mean_pae``. ``chains`` explícito
    tiene prioridad sobre los metadatos del archivo, igual que en
    ``pae_utils.interface_pae``.
    """
    label = _source_label(source)
    kind = _source_kind(source)

    if kind == "file":
        path = Path(source)
        try:
            exists = path.is_file()
        except (TypeError, OSError):
            exists = False
        if not exists:
            return {
                "available": False,
                "status": STATUS_NOT_AVAILABLE,
                "reason": f"archivo PAE ausente (no encontrado): {path}",
                "matrix": None,
                "chains_arg": None,
                "has_chains": False,
                "chain_reason": "",
                "source": label,
                "source_kind": kind,
            }
        try:
            info = pae_utils.load_pae(path)
        except (OSError, ValueError, TypeError) as exc:
            return _invalid_info(f"PAE no utilizable: {exc}", source=source,
                                 available=True)

        if chains is not None:
            chains_arg = chains
            chain_reason = ""
        else:
            chains_arg = info["chains"] or None
            chain_reason = info["chain_reason"] or NO_CHAIN_REASON

        return {
            "available": True,
            "status": STATUS_OK,
            "reason": "ok",
            "matrix": info["matrix"],
            "chains_arg": chains_arg,
            "has_chains": chains_arg is not None,
            "chain_reason": chain_reason,
            "source": label,
            "source_kind": kind,
        }

    # Matriz directa: se valida una vez para no confundir matriz rota con
    # "faltan cadenas".
    try:
        pae_utils.mean_pae(source)
    except (ValueError, TypeError, OverflowError) as exc:
        return _invalid_info(f"PAE no utilizable: {exc}", source=source,
                             available=True)

    if chains is None:
        return {
            "available": True,
            "status": STATUS_OK,
            "reason": "ok",
            "matrix": source,
            "chains_arg": None,
            "has_chains": False,
            "chain_reason": NO_CHAIN_REASON,
            "source": label,
            "source_kind": kind,
        }
    return {
        "available": True,
        "status": STATUS_OK,
        "reason": "ok",
        "matrix": source,
        "chains_arg": chains,
        "has_chains": True,
        "chain_reason": "",
        "source": label,
        "source_kind": kind,
    }


# ---------------------------------------------------------------------------
# Clasificación de cada bloque de ``pae_utils.interface_pae``
# ---------------------------------------------------------------------------
_NOT_AVAILABLE_MARKERS: tuple[str, ...] = (
    "no declara cadenas",
    "sin cadenas",
    "archivo PAE ausente",
    "falta el mapeo",
)


def _classify_interface_result(raw: Mapping[str, Any]) -> str:
    """Traduce ``available``/``ok``/``reason`` de pae_utils a un estado."""
    if raw.get("ok"):
        return STATUS_OK
    if not raw.get("available", True):
        return STATUS_NOT_AVAILABLE
    reason = str(raw.get("reason") or "")
    if reason == NO_CHAIN_REASON or any(
        marker in reason for marker in _NOT_AVAILABLE_MARKERS
    ):
        return STATUS_NOT_AVAILABLE
    return STATUS_INVALID


def _same_chain(left: Any, right: Any) -> bool:
    if left is None or right is None:
        return False
    return str(left).casefold() == str(right).casefold()


def _empty_chain_result(
    chain_a: str,
    chain_b: str,
    *,
    status: str,
    available: bool,
    reason: str,
    resolved_a: str | None = None,
    resolved_b: str | None = None,
) -> dict[str, Any]:
    return {
        "status": status,
        "ok": status == STATUS_OK,
        "available": available,
        "requested_chain_a": chain_a,
        "requested_chain_b": chain_b,
        "resolved_chain_a": resolved_a,
        "resolved_chain_b": resolved_b,
        "n_a": 0,
        "n_b": 0,
        "n_pairs": 0,
        "mean_pae": None,
        "max_pae": None,
        "median_pae": None,
        "reason": reason,
    }


def _chain_pair(
    matrix: Any,
    chains_arg: Any,
    chain_a: str,
    chain_b: str,
) -> dict[str, Any]:
    """Calcula un bloque ``chain_a -> chain_b`` y lo devuelve compacto."""
    try:
        raw = pae_utils.interface_pae(
            matrix, chain_a, chain_b, chains=chains_arg
        )
    except Exception as exc:  # la capa de UI no debe romperse por una entrada rara
        return _empty_chain_result(
            chain_a,
            chain_b,
            status=STATUS_INVALID,
            available=True,
            reason=f"error al calcular la interfaz: {type(exc).__name__}: {exc}",
        )

    status = _classify_interface_result(raw)
    resolved_a = raw.get("chain_a")
    resolved_b = raw.get("chain_b")

    # pae_utils acepta IDs numéricos como índice de orden. En un mapeo de
    # subunidades eso podría asignar una cadena equivocada en silencio; aquí
    # se exige que el ID solicitado coincida con la cadena resuelta.
    if status == STATUS_OK and not (
        _same_chain(resolved_a, chain_a) and _same_chain(resolved_b, chain_b)
    ):
        return _empty_chain_result(
            chain_a,
            chain_b,
            status=STATUS_INVALID,
            available=bool(raw.get("available", True)),
            reason=(
                f"la cadena solicitada {chain_a!r} -> {chain_b!r} resolvió a "
                f"{resolved_a!r} -> {resolved_b!r}; se rechaza la resolución "
                "por orden para no inventar cadenas"
            ),
            resolved_a=resolved_a if isinstance(resolved_a, str) else None,
            resolved_b=resolved_b if isinstance(resolved_b, str) else None,
        )

    if status != STATUS_OK:
        return _empty_chain_result(
            chain_a,
            chain_b,
            status=status,
            available=bool(raw.get("available", True)),
            reason=str(raw.get("reason") or status),
        )

    return {
        "status": STATUS_OK,
        "ok": True,
        "available": True,
        "requested_chain_a": chain_a,
        "requested_chain_b": chain_b,
        "resolved_chain_a": resolved_a,
        "resolved_chain_b": resolved_b,
        "n_a": int(raw.get("n_a") or 0),
        "n_b": int(raw.get("n_b") or 0),
        "n_pairs": int(raw.get("n_pairs") or 0),
        "mean_pae": raw.get("mean"),
        "max_pae": raw.get("max"),
        "median_pae": raw.get("median"),
        "reason": "ok",
    }


# ---------------------------------------------------------------------------
# Agregación por dirección y por par de subunidades
# ---------------------------------------------------------------------------
def _worst_status(statuses: Sequence[str]) -> str:
    """Prioriza ``invalid`` > ``not_available`` > ``ok``."""
    if STATUS_INVALID in statuses:
        return STATUS_INVALID
    if STATUS_NOT_AVAILABLE in statuses:
        return STATUS_NOT_AVAILABLE
    return STATUS_OK


def _unique_reasons(items: Sequence[Mapping[str, Any]], limit: int = 3) -> str:
    reasons: list[str] = []
    for item in items:
        if item.get("status") == STATUS_OK:
            continue
        reason = str(item.get("reason") or item.get("status") or "")
        if reason and reason not in reasons:
            reasons.append(reason)
    if len(reasons) > limit:
        extra = len(reasons) - limit
        reasons = reasons[:limit] + [f"... (+{extra} motivos)"]
    return "; ".join(reasons)


def _empty_direction(
    subunit_a: str,
    subunit_b: str,
    chains_a: Sequence[str],
    chains_b: Sequence[str],
    *,
    status: str,
    available: bool,
    reason: str,
) -> dict[str, Any]:
    expected = len(chains_a) * len(chains_b)
    return {
        "subunit_a": subunit_a,
        "subunit_b": subunit_b,
        "chains_a": list(chains_a),
        "chains_b": list(chains_b),
        "status": status,
        "ok": status == STATUS_OK,
        "available": available,
        "expected_chain_pairs": expected,
        "n_chain_pairs": 0,
        "n_chain_pairs_ok": 0,
        "n_pairs": 0,
        "mean_pae": None,
        "max_pae": None,
        "chain_pairs": [],
        "reason": reason,
    }


def _direction(
    matrix: Any,
    chains_arg: Any,
    subunit_a: str,
    chains_a: Sequence[str],
    subunit_b: str,
    chains_b: Sequence[str],
) -> dict[str, Any]:
    chain_pairs = [
        _chain_pair(matrix, chains_arg, chain_a, chain_b)
        for chain_a in chains_a
        for chain_b in chains_b
    ]
    statuses = [pair["status"] for pair in chain_pairs]
    status = _worst_status(statuses)
    expected = len(chains_a) * len(chains_b)

    if status == STATUS_OK:
        n_pairs = sum(pair["n_pairs"] for pair in chain_pairs)
        if n_pairs:
            weighted = sum(
                float(pair["mean_pae"]) * pair["n_pairs"]  # type: ignore[arg-type]
                for pair in chain_pairs
            )
            mean_value: float | None = weighted / n_pairs
            max_value: float | None = max(
                float(pair["max_pae"]) for pair in chain_pairs  # type: ignore[arg-type]
            )
        else:  # pragma: no cover - una cadena vacía ya sería invalid
            mean_value = None
            max_value = None
        reason = "ok"
    else:
        n_pairs = 0
        mean_value = None
        max_value = None
        reason = _unique_reasons(chain_pairs) or status

    return {
        "subunit_a": subunit_a,
        "subunit_b": subunit_b,
        "chains_a": list(chains_a),
        "chains_b": list(chains_b),
        "status": status,
        "ok": status == STATUS_OK,
        "available": all(pair["available"] for pair in chain_pairs),
        "expected_chain_pairs": expected,
        "n_chain_pairs": len(chain_pairs),
        "n_chain_pairs_ok": sum(
            1 for pair in chain_pairs if pair["status"] == STATUS_OK
        ),
        "n_pairs": n_pairs,
        "mean_pae": mean_value,
        "max_pae": max_value,
        "chain_pairs": chain_pairs,
        "reason": reason,
    }


def _pair(
    source_info: Mapping[str, Any],
    subunit_a: str,
    chains_a: Sequence[str],
    subunit_b: str,
    chains_b: Sequence[str],
) -> dict[str, Any]:
    if source_info["status"] != STATUS_OK:
        a_to_b = _empty_direction(
            subunit_a,
            subunit_b,
            chains_a,
            chains_b,
            status=str(source_info["status"]),
            available=bool(source_info["available"]),
            reason=str(source_info["reason"]),
        )
        b_to_a = _empty_direction(
            subunit_b,
            subunit_a,
            chains_b,
            chains_a,
            status=str(source_info["status"]),
            available=bool(source_info["available"]),
            reason=str(source_info["reason"]),
        )
    elif not source_info["has_chains"]:
        reason = str(source_info["chain_reason"] or NO_CHAIN_REASON)
        a_to_b = _empty_direction(
            subunit_a,
            subunit_b,
            chains_a,
            chains_b,
            status=STATUS_NOT_AVAILABLE,
            available=True,
            reason=reason,
        )
        b_to_a = _empty_direction(
            subunit_b,
            subunit_a,
            chains_b,
            chains_a,
            status=STATUS_NOT_AVAILABLE,
            available=True,
            reason=reason,
        )
    else:
        a_to_b = _direction(
            source_info["matrix"],
            source_info["chains_arg"],
            subunit_a,
            chains_a,
            subunit_b,
            chains_b,
        )
        b_to_a = _direction(
            source_info["matrix"],
            source_info["chains_arg"],
            subunit_b,
            chains_b,
            subunit_a,
            chains_a,
        )

    status = _worst_status((a_to_b["status"], b_to_a["status"]))
    if status == STATUS_OK:
        n_pairs = a_to_b["n_pairs"] + b_to_a["n_pairs"]
        if n_pairs:
            weighted = (
                float(a_to_b["mean_pae"]) * a_to_b["n_pairs"]
                + float(b_to_a["mean_pae"]) * b_to_a["n_pairs"]
            )
            mean_value: float | None = weighted / n_pairs
            max_value: float | None = max(
                float(a_to_b["max_pae"]), float(b_to_a["max_pae"])
            )
        else:  # pragma: no cover
            mean_value = None
            max_value = None
        reason = "ok"
    else:
        n_pairs = 0
        mean_value = None
        max_value = None
        reason = _unique_reasons((a_to_b, b_to_a)) or status

    return {
        "subunit_a": subunit_a,
        "subunit_b": subunit_b,
        "status": status,
        "ok": status == STATUS_OK,
        "a_to_b": a_to_b,
        "b_to_a": b_to_a,
        "n_pairs": n_pairs,
        "mean_pae": mean_value,
        "max_pae": max_value,
        "reason": reason,
    }


# ---------------------------------------------------------------------------
# Resultados vacíos / resumen
# ---------------------------------------------------------------------------
def _empty_summary() -> dict[str, Any]:
    return {
        "n_pairs": 0,
        "n_ok": 0,
        "n_not_available": 0,
        "n_invalid": 0,
        "n_residue_pairs": 0,
        "mean_pae": None,
        "max_pae": None,
    }


def _summary(pairs: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    ok_pairs = [pair for pair in pairs if pair["status"] == STATUS_OK]
    n_values = sum(int(pair["n_pairs"]) for pair in ok_pairs)
    if n_values:
        weighted = sum(
            float(pair["mean_pae"]) * int(pair["n_pairs"])  # type: ignore[arg-type]
            for pair in ok_pairs
        )
        mean_value: float | None = weighted / n_values
        max_value: float | None = max(
            float(pair["max_pae"]) for pair in ok_pairs  # type: ignore[arg-type]
        )
    else:
        mean_value = None
        max_value = None
    return {
        "n_pairs": len(pairs),
        "n_ok": len(ok_pairs),
        "n_not_available": sum(
            1 for pair in pairs if pair["status"] == STATUS_NOT_AVAILABLE
        ),
        "n_invalid": sum(
            1 for pair in pairs if pair["status"] == STATUS_INVALID
        ),
        "n_residue_pairs": n_values,
        "mean_pae": mean_value,
        "max_pae": max_value,
    }


def _analysis_reason(
    status: str, pairs: Sequence[Mapping[str, Any]]
) -> str:
    if status == STATUS_OK:
        return "ok"
    if not pairs:
        return "sin pares de subunidades"
    bad = [pair for pair in pairs if pair["status"] != STATUS_OK]
    prefix = (
        "PAE de interfaz no disponible"
        if status == STATUS_NOT_AVAILABLE
        else "PAE de interfaz inválido"
    )
    reasons = _unique_reasons(bad)
    if len(bad) == 1:
        return f"{prefix}: {reasons}"
    return f"{prefix} en {len(bad)} de {len(pairs)} pares: {reasons}"


def _empty_analysis(
    *,
    source: Any,
    status: str,
    reason: str,
    subunits: SubunitMap | None = None,
) -> dict[str, Any]:
    return {
        "status": status,
        "ok": status == STATUS_OK,
        "available": _source_available(source),
        "reason": reason,
        "source": _source_label(source),
        "source_kind": _source_kind(source),
        "subunits": dict(subunits or {}),
        "n_subunits": len(subunits or {}),
        "pairs": [],
        "summary": _empty_summary(),
        "limitation": LIMITATION,
    }


# ---------------------------------------------------------------------------
# API pública
# ---------------------------------------------------------------------------
def analyze_complex(
    source: str | os.PathLike[str] | Sequence[Sequence[float]],
    subunits: Mapping[str, Any] | None,
    chains: dict[str, Any] | list[Any] | None = None,
) -> dict[str, Any]:
    """Analiza todos los pares de subunidades de un complejo.

    Parámetros
    ----------
    source:
        Ruta a un ``*.pae.json`` local o matriz PAE ``NxN``. No se descarga
        nada.
    subunits:
        Mapeo ``{subunidad: cadena | [cadenas]}``. Se exigen al menos dos
        subunidades. Cada cadena debe aparecer una sola vez. El orden de las
        subunidades define el orden de ``pairs`` (combinaciones sin repetir).
    chains:
        Mapeo opcional de cadenas a índices de residuo, con el mismo formato
        que acepta :func:`webapp.pae_utils.interface_pae` (por ejemplo
        ``{"A": [0, 1], "B": [2, 3]}`` o una serie por residuo). Si se aporta,
        tiene prioridad sobre los metadatos del archivo.

    Devuelve
    --------
    dict
        ``status`` (``ok``/``not_available``/``invalid``), ``ok``,
        ``available``, ``reason``, ``source``, ``source_kind``, ``subunits``,
        ``n_subunits``, ``pairs``, ``summary`` y ``limitation``.

        Cada entrada de ``pairs`` contiene ``a_to_b``/``b_to_a`` con
        ``status``, ``mean_pae``, ``max_pae``, ``n_pairs`` y el detalle de
        ``chain_pairs``. La media combinada está ponderada por número de pares
        de residuos.
    """
    if subunits is None:
        return _empty_analysis(
            source=source,
            status=STATUS_NOT_AVAILABLE,
            reason=NO_MAPPING_REASON,
        )

    try:
        mapping = _normalize_subunits(subunits, minimum=2)
    except (ValueError, TypeError) as exc:
        return _empty_analysis(
            source=source,
            status=STATUS_INVALID,
            reason=f"mapeo de subunidades inválido: {exc}",
        )

    source_info = _load_source(source, chains)
    labels = list(mapping)
    pairs = [
        _pair(source_info, label_a, mapping[label_a], label_b, mapping[label_b])
        for label_a, label_b in combinations(labels, 2)
    ]
    status = _worst_status([pair["status"] for pair in pairs])

    return {
        "status": status,
        "ok": status == STATUS_OK,
        "available": bool(source_info["available"]),
        "reason": _analysis_reason(status, pairs),
        "source": source_info["source"],
        "source_kind": source_info["source_kind"],
        "subunits": mapping,
        "n_subunits": len(mapping),
        "pairs": pairs,
        "summary": _summary(pairs),
        "limitation": LIMITATION,
    }


def interface_pair(
    source: str | os.PathLike[str] | Sequence[Sequence[float]],
    subunit_a: str,
    subunit_b: str,
    subunits: Mapping[str, Any] | None,
    chains: dict[str, Any] | list[Any] | None = None,
) -> dict[str, Any]:
    """Analiza un único par de subunidades, en ambas direcciones.

    Es una vista de :func:`analyze_complex` restringida a ``subunit_a`` y
    ``subunit_b``. Devuelve el par con ``a_to_b``/``b_to_a`` más los metadatos
    de la fuente. Si el mapeo o las etiquetas no son válidos, devuelve un par
    con ``status`` explícito y direcciones vacías; no lanza por datos.
    """
    try:
        label_a = _subunit_label(subunit_a, where="subunidad_a")
        label_b = _subunit_label(subunit_b, where="subunidad_b")
    except ValueError as exc:
        empty = _empty_analysis(
            source=source, status=STATUS_INVALID, reason=str(exc)
        )
        return _pair_view_from_analysis(
            empty, str(subunit_a), str(subunit_b)
        )

    if label_a == label_b:
        empty = _empty_analysis(
            source=source,
            status=STATUS_INVALID,
            reason=(
                f"subunidad_a y subunit_b son la misma subunidad ({label_a!r}); "
                "no hay interfaz entre subunidades distintas"
            ),
        )
        return _pair_view_from_analysis(empty, label_a, label_b)

    if subunits is None:
        empty = _empty_analysis(
            source=source, status=STATUS_NOT_AVAILABLE, reason=NO_MAPPING_REASON
        )
        return _pair_view_from_analysis(empty, label_a, label_b)

    try:
        mapping = _normalize_subunits(subunits, minimum=2)
    except (ValueError, TypeError) as exc:
        empty = _empty_analysis(
            source=source,
            status=STATUS_INVALID,
            reason=f"mapeo de subunidades inválido: {exc}",
        )
        return _pair_view_from_analysis(empty, label_a, label_b)

    missing = [
        label for label in (label_a, label_b) if label not in mapping
    ]
    if missing:
        empty = _empty_analysis(
            source=source,
            status=STATUS_INVALID,
            reason=(
                "subunidad no encontrada: "
                + ", ".join(repr(label) for label in missing)
                + "; disponibles: "
                + ", ".join(mapping)
            ),
            subunits=mapping,
        )
        return _pair_view_from_analysis(empty, label_a, label_b)

    source_info = _load_source(source, chains)
    pair = _pair(
        source_info,
        label_a,
        mapping[label_a],
        label_b,
        mapping[label_b],
    )
    return {
        **pair,
        "available": bool(source_info["available"]),
        "source": source_info["source"],
        "source_kind": source_info["source_kind"],
        "subunits": mapping,
        "limitation": LIMITATION,
    }


def _pair_view_from_analysis(
    analysis: Mapping[str, Any], subunit_a: str, subunit_b: str
) -> dict[str, Any]:
    """Construye una vista de par cuando no se pudo llegar a ``_pair``."""
    status = str(analysis["status"])
    reason = str(analysis["reason"])
    available = bool(analysis["available"])
    a_to_b = _empty_direction(
        subunit_a,
        subunit_b,
        [],
        [],
        status=status,
        available=available,
        reason=reason,
    )
    b_to_a = _empty_direction(
        subunit_b,
        subunit_a,
        [],
        [],
        status=status,
        available=available,
        reason=reason,
    )
    return {
        "subunit_a": subunit_a,
        "subunit_b": subunit_b,
        "status": status,
        "ok": status == STATUS_OK,
        "a_to_b": a_to_b,
        "b_to_a": b_to_a,
        "n_pairs": 0,
        "mean_pae": None,
        "max_pae": None,
        "reason": reason,
        "available": available,
        "source": analysis["source"],
        "source_kind": analysis["source_kind"],
        "subunits": dict(analysis.get("subunits") or {}),
        "limitation": LIMITATION,
    }


__all__ = [
    "LIMITATION",
    "NO_CHAIN_REASON",
    "NO_MAPPING_REASON",
    "STATUSES",
    "STATUS_INVALID",
    "STATUS_NOT_AVAILABLE",
    "STATUS_OK",
    "analyze_complex",
    "interface_pair",
    "subunits_from_chain_subunits",
]
