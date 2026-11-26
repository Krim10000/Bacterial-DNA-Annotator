"""Políticas de decisión para sigma, separadas del escáner de motivos.

Ninguna política inventa un sigma que no esté en el inventario de la especie.
``argmax`` es el comportamiento histórico; ``two_stage`` exige que una
alternativa supere a sigma70 por un margen; ``hybrid`` usa una probabilidad
externa (p. ej. k-mer) y cae a two_stage si no es concluyente.
"""
from __future__ import annotations

from typing import Mapping


def _order(scores: Mapping[str, float], key: str = "z") -> list[str]:
    return [s for _, s in sorted(
        ((float(v or 0), s) for s, v in scores.items()), reverse=True)]


def argmax_order(scores: Mapping[str, float], **_: object) -> list[str]:
    return _order(scores)


def two_stage_order(scores: Mapping[str, float], *,
                    delta: float = 1.0,
                    fallback: str = "sigma70") -> list[str]:
    """Predice alternativa solo si su z supera al fallback por ``delta``."""
    order = _order(scores)
    if not order:
        return []
    base = float(scores.get(fallback, 0) or 0)
    alt = [s for s in order if s != fallback]
    if alt and float(scores.get(alt[0], 0) or 0) - base > float(delta):
        return [alt[0]] + [s for s in order if s != alt[0]]
    return [fallback] + [s for s in order if s != fallback]


def hybrid_order(scores: Mapping[str, float],
                 external: Mapping[str, float], *,
                 threshold: float = 0.75,
                 delta: float = 2.0,
                 fallback: str = "sigma70") -> list[str]:
    """Usa la señal externa solo si es concluyente.

    ``external`` puede ser una probabilidad k-mer o una calibración; no se
    inventa: si no supera ``threshold``, se usa la política two-stage.
    """
    if external:
        best = max(external, key=lambda s: float(external[s] or 0))
        if float(external.get(best, 0) or 0) >= float(threshold):
            return [best] + [s for s in _order(scores) if s != best]
    return two_stage_order(scores, delta=delta, fallback=fallback)


__all__ = ["argmax_order", "two_stage_order", "hybrid_order"]
