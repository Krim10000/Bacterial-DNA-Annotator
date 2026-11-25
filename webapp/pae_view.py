"""Heatmap PAE sin dependencias nuevas (SVG generado en Python)."""
from __future__ import annotations

from typing import Sequence


def _color(value: float, vmax: float) -> str:
    x = max(0.0, min(1.0, value / vmax if vmax > 0 else 0.0))
    if x < 0.5:
        # verde -> amarillo
        r = int(40 + 215 * (x / 0.5))
        g = int(160 + 60 * (x / 0.5))
        b = 60
    else:
        # amarillo -> rojo
        y = (x - 0.5) / 0.5
        r = 255
        g = int(220 - 180 * y)
        b = int(60 - 40 * y)
    return f"rgb({r},{g},{b})"


def heatmap_svg(matrix: Sequence[Sequence[float]], *, max_cells: int = 64,
                width: int = 560, height: int = 560) -> str:
    """Dibuja una matriz PAE reducida a una cuadrícula legible.

    No interpola: toma índices equiespaciados. No representa distancias ni
    contactos; solo el error relativo predicho.
    """
    n = len(matrix or [])
    if n == 0:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="320" height="60"><text x="8" y="28" font-family="sans-serif" font-size="12">PAE no disponible</text></svg>'
    k = min(int(max_cells), n)
    idx = sorted({round(i * (n - 1) / max(1, k - 1)) for i in range(k)})
    m = len(idx)
    pad_l, pad_t, pad_r, pad_b = 38, 16, 10, 26
    cell = min((width - pad_l - pad_r) / m, (height - pad_t - pad_b) / m)
    vals = [float(matrix[i][j]) for i in idx for j in idx]
    vmax = max(vals) if vals else 1.0
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}" font-family="sans-serif">',
           f'<text x="8" y="12" font-size="10" fill="#475569">PAE (Å) · {n}×{n} residuos</text>']
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            out.append(f'<rect x="{pad_l + b * cell:.1f}" y="{pad_t + a * cell:.1f}" '
                       f'width="{cell + 0.4:.1f}" height="{cell + 0.4:.1f}" '
                       f'fill="{_color(float(matrix[i][j]), vmax)}"/>')
    out.append(f'<text x="{pad_l}" y="{height - 8}" font-size="9" fill="#475569">1</text>')
    out.append(f'<text x="{pad_l + (m - 1) * cell}" y="{height - 8}" font-size="9" '
               f'fill="#475569">{n}</text>')
    out.append(f'<text x="2" y="{pad_t + 8}" font-size="9" fill="#475569">1</text>')
    out.append(f'<text x="2" y="{pad_t + (m - 1) * cell + 8}" font-size="9" '
               f'fill="#475569">{n}</text>')
    out.append("</svg>")
    return "".join(out)
