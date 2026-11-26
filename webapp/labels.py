"""Etiquetas de proceso anidadas para genes/operones (prototipo)."""
from __future__ import annotations

import re
from typing import Iterable

# Léxico controlado: términos que indican procesos, no funciones aisladas.
LABELS: dict[str, tuple[str, ...]] = {
    "heat_shock": (
        "heat shock", "chaperon", "chaperonin", "dnak", "dnaj", "groel",
        "groes", "grpe", "htpg", "clpb", "clpx", "clpp", "hslu", "hslv",
        "ibpa", "ibpb", "lon protease", "protease la", "rpoH",
    ),
    "flagellum": (
        "flagell", "flagellar", "flg", "flh", "fli", "motility", "chemotaxis",
        "chea", "cheb", "chey", "chez", "chew", "flia",
    ),
    "nitrogen": (
        "nitrogen", "nif", "glutamine", "glutamate synthase", "gln", "amt",
        "nitrate", "nitrite", "nar", "nir", "nor", "urease", "ure",
    ),
    "stationary": (
        "stationary", "rpos", "catalase", "superoxide dismutase", "oxidative",
        "osmotic", "osm", "trehalose", "glycogen", "rpos",
    ),
    "envelope": (
        "envelope", "outer membrane", "periplasm", "cell wall", "murein",
        "lipopolysaccharide", "lps", "capsule", "peptidoglycan", "rpoe",
    ),
    "sporulation": (
        "sporulation", "spore", "sigf", "sige", "sigk", "spo", "germination",
    ),
}
PARENTS = {
    "heat_shock": "stress_response",
    "envelope": "stress_response",
    "stationary": "stress_response",
    "nitrogen": "metabolism",
    "flagellum": "motility",
    "sporulation": "differentiation",
}
_ANTI = re.compile(r"anti[- ]?sigma|regulator of sigma|sigma[- ]?dependent", re.I)


def labels_for_text(text: str) -> list[str]:
    low = (text or "").lower()
    out = []
    for label, terms in LABELS.items():
        if any(t.lower() in low for t in terms):
            out.append(label)
    return out


def is_regulator(name: str, product: str) -> bool:
    """Evita etiquetar el propio factor sigma como si fuera su diana."""
    text = f"{name} {product}".lower()
    return bool(_ANTI.search(text)) or bool(re.search(
        r"\bsigma factor\b|\brpo[densh]\b|\bsig[a-z]\b|\bflia\b", text))


def enrich(rows: Iterable[tuple[bool, bool]]) -> tuple[float, float]:
    """Fisher exacto (OR, p) para 2x2; sin dependencias externas.

    ``rows``: (tiene_label, es_diana). Usa scipy si está; si no, devuelve
    aproximación conservadora.
    """
    a = sum(1 for x, y in rows if x and y)
    b = sum(1 for x, y in rows if x and not y)
    c = sum(1 for x, y in rows if not x and y)
    d = sum(1 for x, y in rows if not x and not y)
    try:
        from scipy.stats import fisher_exact
        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        return float(odds), float(p)
    except Exception:
        # OR de Haldane-Anscombe.
        return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)), 1.0
