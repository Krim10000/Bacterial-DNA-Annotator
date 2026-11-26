"""Estructura secundaria de ARN para tRNA y horquillas (predicción, no 3D).

API pura y determinista::

    from webapp.rna_structure import predict_structure

    r = predict_structure("GGGGAAAACCCC")
    r.dot_bracket   # "((((....))))"
    r.energy_kcal_mol  # p. ej. -5.4 si ViennaRNA está disponible
    r.source        # "vienna" | "nussinov" | "not_available"

Fuentes
-------
``vienna``
    ``RNAfold`` de ViennaRNA con sus parámetros por defecto (37 °C, modelo
    Turner 2004, sin restricciones ni pseudoknots). La energía es la MFE en
    kcal/mol. Es la fuente preferida.
``nussinov``
    Respaldo puro Python: maximiza el número de pares canónicos
    (AU, UA, GC, CG, GU, UG), sin pseudoknots y con tallo-bucle mínimo de
    3 nt. El ``score`` es el **número de pares**, no una energía. Sirve para
    tRNA/horquillas cortas cuando falta ViennaRNA, no para comparar ΔG.
``not_available``
    No se pudo predecir. ``dot_bracket`` queda vacío: no se inventan pares.

Ninguna fuente devuelve coordenadas 3D, valida una estructura experimental ni
dibuja nada. La GUI no usa este módulo todavía.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

Source = Literal["vienna", "nussinov", "not_available"]

#: Pares de bases permitidos (incluye el par wobble G-U).
CANONICAL_PAIRS = frozenset({
    ("A", "U"), ("U", "A"),
    ("G", "C"), ("C", "G"),
    ("G", "U"), ("U", "G"),
})

#: Tamaño mínimo del bucle de una horquilla, igual que el valor por defecto
#: de RNAfold. Un par (i, j) exige j - i - 1 >= MIN_HAIRPIN.
MIN_HAIRPIN = 3

#: Límite del respaldo Nussinov (O(n^3) en Python puro). Por encima se declara
#: ``not_available`` en lugar de tardar minutos.
FALLBACK_MAX_LENGTH = 200

#: Límite duro de la API para cualquier backend. tRNA y horquillas están muy
#: por debajo; evita que una secuencia cromosómica complete se cuele.
MAX_SEQUENCE_LENGTH = 1000

_BACKENDS = ("auto", "vienna", "nussinov")


@dataclass(frozen=True)
class RnaStructure:
    """Resultado de plegar una secuencia de ARN.

    ``pairs`` usa índices 0-based ``(i, j)`` con ``i < j`` sobre
    ``sequence``. ``energy_kcal_mol`` solo existe para ``vienna``.
    ``score`` es el número de pares de bases en todas las fuentes.
    """

    sequence: str
    dot_bracket: str
    energy_kcal_mol: float | None
    score: int | None
    source: Source
    pairs: tuple[tuple[int, int], ...] = ()
    reason: str = ""

    @property
    def available(self) -> bool:
        """True si hay una estructura (aunque sea sin pares)."""
        return self.source != "not_available"

    def as_dict(self) -> dict[str, Any]:
        """Versión JSON-serializable para conectar la API a la GUI."""
        return {
            "sequence": self.sequence,
            "dot_bracket": self.dot_bracket,
            "energy_kcal_mol": self.energy_kcal_mol,
            "score": self.score,
            "source": self.source,
            "pairs": [list(p) for p in self.pairs],
            "reason": self.reason,
        }


def normalize_sequence(sequence: str) -> str:
    """Normaliza a ARN en mayúsculas: quita espacios y convierte T -> U.

    Solo acepta A, C, G y U. Una base ambigua (N, R, Y, ...) no se puede
    emparejar sin inventar: se rechaza con ``ValueError`` en vez de adivinar.
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence debe ser str")
    clean = "".join(sequence.split()).upper().replace("T", "U")
    for pos, base in enumerate(clean):
        if base not in "ACGU":
            raise ValueError(
                f"base no soportada {base!r} en posición {pos}: "
                "solo se aceptan A, C, G, U/T (sin ambigüedades)"
            )
    return clean


def _can_pair(a: str, b: str) -> bool:
    return (a, b) in CANONICAL_PAIRS


def pairs_from_dot_bracket(
    dot_bracket: str, sequence: str | None = None
) -> tuple[tuple[int, int], ...]:
    """Convierte dot-bracket en pares 0-based ``(i, j)``.

    Si se pasa ``sequence``, comprueba longitud y que cada par sea canónico.
    Solo admite ``.``, ``(`` y ``)``: los pseudoknots (``[]``, ``{}``, ...)
    no están soportados y se rechazan en lugar de reinterpretarse.
    """
    if not isinstance(dot_bracket, str):
        raise TypeError("dot_bracket debe ser str")
    seq = normalize_sequence(sequence) if sequence is not None else None
    if seq is not None and len(seq) != len(dot_bracket):
        raise ValueError(
            f"longitud distinta: secuencia={len(seq)} dot_bracket={len(dot_bracket)}"
        )
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for pos, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(pos)
        elif char == ")":
            if not stack:
                raise ValueError(f"')' sin '(' en posición {pos}")
            start = stack.pop()
            if seq is not None and not _can_pair(seq[start], seq[pos]):
                raise ValueError(
                    f"par no canónico {seq[start]}-{seq[pos]} en ({start}, {pos})"
                )
            pairs.append((start, pos))
        elif char != ".":
            raise ValueError(
                f"carácter {char!r} no soportado en posición {pos}: "
                "solo dot-bracket simple sin pseudoknots"
            )
    if stack:
        raise ValueError(f"hay '(' sin cerrar en posiciones {stack}")
    return tuple(sorted(pairs))


def dot_bracket_from_pairs(
    length: int, pairs: tuple[tuple[int, int], ...] | list[tuple[int, int]]
) -> str:
    """Construye dot-bracket a partir de pares 0-based no cruzados.

    No comprueba que los pares sean canónicos (no recibe la secuencia); el
    respaldo Nussinov solo le pasa pares canónicos. Sí rechaza índices fuera
    de rango, bases repetidas y pseudoknots.
    """
    if not isinstance(length, int) or length < 0:
        raise ValueError("length debe ser un entero >= 0")
    chars = ["."] * length
    used: set[int] = set()
    accepted: list[tuple[int, int]] = []
    for pair in pairs:
        if len(pair) != 2:
            raise ValueError(f"par inválido {pair!r}")
        i, j = int(pair[0]), int(pair[1])
        if not (0 <= i < j < length):
            raise ValueError(f"par fuera de rango o desordenado: {(i, j)}")
        if i in used or j in used:
            raise ValueError(f"base repetida en el par {(i, j)}")
        for a, b in accepted:
            if (i < a < j < b) or (a < i < b < j):
                raise ValueError(f"pseudoknot entre {(a, b)} y {(i, j)}")
        used.update((i, j))
        accepted.append((i, j))
        chars[i], chars[j] = "(", ")"
    return "".join(chars)


def _fold_vienna(sequence: str, module: Any) -> tuple[str, float]:
    """Pliega con RNAfold. Devuelve (dot_bracket, MFE kcal/mol)."""
    folded = module.fold(sequence)
    if not isinstance(folded, (tuple, list)) or len(folded) < 2:
        raise RuntimeError("RNA.fold devolvió un resultado inesperado")
    structure, energy = folded[0], folded[1]
    if isinstance(structure, bytes):
        structure = structure.decode("ascii")
    structure = str(structure).strip()
    if len(structure) != len(sequence):
        raise RuntimeError(
            f"RNA.fold devolvió {len(structure)} posiciones para {len(sequence)} nt"
        )
    # Valida el dot-bracket de Vienna (y los pares) antes de confiar en él.
    pairs_from_dot_bracket(structure, sequence)
    return structure, float(energy)


def _fold_nussinov(sequence: str) -> tuple[str, int]:
    """Máximo de pares canónicos sin pseudoknots (programación dinámica).

    ``dp[i][j]`` = máximo de pares en ``sequence[i:j+1]``. Es determinista:
    ante empates se prefiere dejar ``i`` sin emparejar, luego ``j``, luego el
    par ``(i, j)`` y, por último, la primera bifurcación ``k``. No usa
    parámetros termodinámicos, así que el resultado no es una MFE.
    """
    n = len(sequence)
    if n < MIN_HAIRPIN + 2:
        return "." * n, 0

    dp = [[0] * n for _ in range(n)]
    for span in range(MIN_HAIRPIN + 1, n):
        for i in range(n - span):
            j = i + span
            best = dp[i + 1][j]
            skip_right = dp[i][j - 1]
            if skip_right > best:
                best = skip_right
            if _can_pair(sequence[i], sequence[j]):
                pair_inside = 1 + dp[i + 1][j - 1]
                if pair_inside > best:
                    best = pair_inside
            for k in range(i, j):
                split = dp[i][k] + dp[k + 1][j]
                if split > best:
                    best = split
            dp[i][j] = best

    pairs: list[tuple[int, int]] = []
    pending: list[tuple[int, int]] = [(0, n - 1)]
    while pending:
        i, j = pending.pop()
        if j - i <= MIN_HAIRPIN:
            continue
        if dp[i][j] == dp[i + 1][j]:
            pending.append((i + 1, j))
        elif dp[i][j] == dp[i][j - 1]:
            pending.append((i, j - 1))
        elif _can_pair(sequence[i], sequence[j]) and (
            dp[i][j] == 1 + dp[i + 1][j - 1]
        ):
            pairs.append((i, j))
            pending.append((i + 1, j - 1))
        else:
            for k in range(i, j):
                if dp[i][k] + dp[k + 1][j] == dp[i][j]:
                    pending.append((i, k))
                    pending.append((k + 1, j))
                    break
            else:  # pragma: no cover - solo si la DP y la traza divergen
                raise RuntimeError("traza Nussinov inconsistente")

    pairs.sort()
    return dot_bracket_from_pairs(n, pairs), len(pairs)


def _import_vienna() -> Any | None:
    """Importa RNA (ViennaRNA) o devuelve None. Aislado para poder testear."""
    try:
        import RNA  # type: ignore[import-not-found]
    except Exception:
        return None
    if not hasattr(RNA, "fold"):
        return None
    return RNA


def _not_available(sequence: str, reason: str) -> RnaStructure:
    return RnaStructure(
        sequence=sequence,
        dot_bracket="",
        energy_kcal_mol=None,
        score=None,
        source="not_available",
        pairs=(),
        reason=reason,
    )


def predict_structure(
    sequence: str,
    *,
    backend: str = "auto",
    allow_fallback: bool = True,
) -> RnaStructure:
    """Predice la estructura secundaria de una secuencia de ARN.

    Parámetros
    ----------
    sequence:
        Secuencia en ARN o ADN (T se convierte a U). Solo A/C/G/U; se quitan
        espacios y saltos de línea. Las bases ambiguas se rechazan.
    backend:
        ``"auto"`` (Vienna si está, si no Nussinov), ``"vienna"`` o
        ``"nussinov"`` (forzar el respaldo puro Python).
    allow_fallback:
        Si es ``False`` y Vienna no está disponible o falla, devuelve
        ``not_available`` en lugar de usar Nussinov.

    Returns
    -------
    RnaStructure
        ``source`` indica el algoritmo; ``energy_kcal_mol`` solo está en
        ``vienna``; ``score`` es el número de pares. Nunca rellena
        ``dot_bracket`` con pares no calculados.

    Raises
    ------
    ValueError
        Secuencia con bases no soportadas, ``backend`` desconocido o
        secuencia más larga que ``MAX_SEQUENCE_LENGTH``.
    """
    if backend not in _BACKENDS:
        raise ValueError(
            f"backend desconocido {backend!r}; usa uno de {', '.join(_BACKENDS)}"
        )
    seq = normalize_sequence(sequence)
    if not seq:
        return _not_available("", "secuencia vacía")
    if len(seq) > MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"secuencia de {len(seq)} nt: el límite de la API es "
            f"{MAX_SEQUENCE_LENGTH} nt (tRNA/horquillas)"
        )

    fallback_reason = ""
    if backend in ("auto", "vienna"):
        module = _import_vienna()
        if module is not None:
            try:
                structure, energy = _fold_vienna(seq, module)
                pairs = pairs_from_dot_bracket(structure, seq)
                return RnaStructure(
                    sequence=seq,
                    dot_bracket=structure,
                    energy_kcal_mol=energy,
                    score=len(pairs),
                    source="vienna",
                    pairs=pairs,
                    reason="",
                )
            except Exception as exc:  # RNA.fold puede fallar con entradas raras
                if not allow_fallback:
                    return _not_available(
                        seq, f"ViennaRNA/RNAfold falló: {type(exc).__name__}: {exc}"
                    )
                fallback_reason = (
                    "ViennaRNA/RNAfold falló "
                    f"({type(exc).__name__}: {exc}); se usó Nussinov."
                )
        elif not allow_fallback:
            return _not_available(
                seq,
                "ViennaRNA/RNAfold no está instalado y el respaldo "
                "Nussinov está desactivado.",
            )
        else:
            fallback_reason = (
                "ViennaRNA/RNAfold no está instalado; se usó Nussinov."
            )

    # Nussinov: explícito o respaldo.
    if len(seq) > FALLBACK_MAX_LENGTH:
        extra = " y ViennaRNA/RNAfold no está disponible" if fallback_reason else ""
        return _not_available(
            seq,
            f"secuencia de {len(seq)} nt: el respaldo Nussinov solo cubre "
            f"hasta {FALLBACK_MAX_LENGTH} nt{extra}.",
        )
    structure, n_pairs = _fold_nussinov(seq)
    reason = fallback_reason or (
        "backend=nussinov: máximo de pares canónicos; el score es el número "
        "de pares, no una energía."
    )
    return RnaStructure(
        sequence=seq,
        dot_bracket=structure,
        energy_kcal_mol=None,
        score=n_pairs,
        source="nussinov",
        pairs=pairs_from_dot_bracket(structure, seq),
        reason=reason,
    )


__all__ = [
    "CANONICAL_PAIRS",
    "FALLBACK_MAX_LENGTH",
    "MAX_SEQUENCE_LENGTH",
    "MIN_HAIRPIN",
    "RnaStructure",
    "dot_bracket_from_pairs",
    "normalize_sequence",
    "pairs_from_dot_bracket",
    "arc_svg",
    "predict_structure",
]


def arc_svg(sequence: str, dot_bracket: str, *, width: int = 900,
            height: int = 190, max_nt: int = 160) -> str:
    """Diagrama 2D de arcos a partir de un dot-bracket.

    No es una estructura 3D ni un layout físico: coloca las bases en una fila
    y dibuja un arco por par. Devuelve ``""`` si la estructura no es válida o
    es demasiado larga para que el dibujo sea legible.
    """
    seq = normalize_sequence(sequence)
    pairs = pairs_from_dot_bracket(dot_bracket, seq)
    n = len(seq)
    if n == 0 or n > max_nt:
        return ""
    margin_x, base_y = 30, height - 42
    step = (width - 2 * margin_x) / max(1, n - 1)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
           f'height="{height}" viewBox="0 0 {width} {height}" '
           'font-family="ui-monospace,monospace">']
    # Bases en una fila.
    for i, base in enumerate(seq):
        x = margin_x + i * step
        out.append(f'<text x="{x:.1f}" y="{base_y}" font-size="9" '
                   f'text-anchor="middle" fill="#334155">{base}</text>')
    # Arcos; los pares largos se elevan más.
    for i, j in pairs:
        x1 = margin_x + i * step
        x2 = margin_x + j * step
        span = max(1, j - i)
        lift = min(120.0, 14.0 + span * 2.2)
        mid = (x1 + x2) / 2
        y = base_y - 8
        color = "#2563eb" if span <= 12 else "#7c3aed"
        out.append(f'<path d="M {x1:.1f} {y:.1f} Q {mid:.1f} '
                   f'{y-lift:.1f} {x2:.1f} {y:.1f}" fill="none" '
                   f'stroke="{color}" stroke-width="1.5" opacity="0.75"/>')
    out.append(f'<text x="{margin_x}" y="18" font-size="11" fill="#475569">'
               f'{n} nt · {len(pairs)} pares · 5′→3′</text>')
    out.append("</svg>")
    return "".join(out)
