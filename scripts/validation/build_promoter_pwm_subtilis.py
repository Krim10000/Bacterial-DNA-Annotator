#!/usr/bin/env python3
"""build_promoter_pwm_subtilis.py — PWM del promotor σA de Bacillus subtilis.

Por qué existe
--------------
La regla R3 del buscador de operones usa el PWM de σ70 de *E. coli*. Aplicada a
*B. subtilis* se dispara 1.779 veces y fragmenta los operones largos. La
hipótesis es que el problema es la **matriz**, no la idea de usar promotores.

Hallazgo al construirla
-----------------------
Con los 1.949 TSS de σA de BSGatlas, la caja **−10** aparece nítida
(`TATAAT` alrededor de −12), pero **no hay caja −35**: al escanear hasta 90 pb
aguas arriba, el mejor 6-mero en esa zona es poli-T (`TTTTTT`, 0,54 bits), no
`TTGACA`. Los promotores σA de este set son en su mayoría **extended −10**.

Por eso la matriz es de **una sola caja**: 12 pb que cubren la caja −10 y la
región de inicio (−12..−1). Es lo que los datos sostienen; inventar una caja
−35 con poli-T sería ruido con forma de motivo.

Método
------
- TSS de σA con coordenada válida (1.949 de 1.949).
- Ventana de 45 pb aguas arriba, corregida por hebra.
- Anchura y desplazamiento elegidos por información, restringidos a la región
  del promotor (−20..−1), no a toda la ventana.
- Log-odds con pseudoconteo 0.5 y la composición real del genoma como fondo:
      pwm[b] = log2((n_b + 0.5) / (N + 2)) - log2(bg_b)

Salida: `data/dbtbs/promoter_pwms_subtilis.tsv` (mismo esquema que el de E. coli).

Uso
---
    python scripts/validation/build_promoter_pwm_subtilis.py
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TSS = ROOT / "data" / "dbtbs" / "gt_tss.tsv"
FASTA = ROOT / "data" / "dbtbs" / "genome.fna"
ECOLI_PWM = ROOT / "data" / "regulondb" / "gt" / "promoter_pwms.tsv"
OUT = ROOT / "data" / "dbtbs" / "promoter_pwms_subtilis.tsv"
BASES = "ACGT"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta(path: Path) -> str:
    parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def read_tss(path: Path, sigma: str) -> list[tuple[int, str]]:
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r.get("sigma", "").strip() != sigma:
            continue
        try:
            pos = int(r["pos"])
        except (KeyError, ValueError):
            continue
        strand = r.get("strand", "+").strip()
        if strand in ("+", "-") and pos >= 1:
            out.append((pos, strand))
    return out


def upstream(genome: str, pos: int, strand: str, width: int) -> str | None:
    """Secuencia aguas arriba del TSS, en orientación 5'->3' del transcrito."""
    if strand == "+":
        a, b = pos - 1 - width, pos - 1
        return genome[a:b] if a >= 0 else None
    a, b = pos, pos + width
    return genome[a:b].translate(COMP)[::-1] if b <= len(genome) else None


def build(windows: list[str], lo: int, hi: int, bg: dict[str, float]) -> dict:
    """PWM log-odds de una caja, con ventanas ya alineadas."""
    width = hi - lo
    counts = [[0] * 4 for _ in range(width)]
    for w in windows:
        for i, ch in enumerate(w[lo:hi]):
            if ch in BASES:
                counts[i][BASES.index(ch)] += 1
    n = len(windows)
    pwm, probs = [], []
    for i in range(width):
        col = [(counts[i][j] + 0.5) / (n + 2.0) for j in range(4)]
        probs.append(col)
        pwm.append([round(math.log2(col[j] / bg[BASES[j]]), 4) for j in range(4)])
    cons = "".join(BASES[max(range(4), key=lambda j: counts[i][j])] for i in range(width))
    info = round(sum(2.0 + sum(p * math.log2(p) for p in col if p > 0) for col in probs), 3)
    return {"width": width, "n_sites": n, "consensus": cons,
            "info_bits": info, "pwm": pwm}


def corr(a: list[list[float]], b: list[list[float]]) -> float:
    xa = [v for col in a for v in col]
    xb = [v for col in b for v in col]
    n = min(len(xa), len(xb))
    xa, xb = xa[:n], xb[:n]
    ma, mb = sum(xa) / n, sum(xb) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(xa, xb))
    da = math.sqrt(sum((x - ma) ** 2 for x in xa))
    db = math.sqrt(sum((y - mb) ** 2 for y in xb))
    return cov / (da * db) if da and db else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sigma", default="A")
    ap.add_argument("--up", type=int, default=45)
    a = ap.parse_args()

    genome = read_fasta(FASTA)
    bg = {b: genome.count(b) / len(genome) for b in BASES}
    tss = read_tss(TSS, a.sigma)
    windows = [w for w in (upstream(genome, p, s, a.up) for p, s in tss) if w]
    print(f"TSS σ{a.sigma}: {len(tss)} | ventanas válidas: {len(windows)} | "
          f"GC {100*(bg['G']+bg['C']):.1f}%")

    # ¿hay caja -35? se informa, no se inventa
    mejor35 = max(((build(windows, lo, lo + 6, bg)["info_bits"], lo,
                    build(windows, lo, lo + 6, bg)["consensus"])
                   for lo in range(0, a.up - 18)), key=lambda x: x[0])
    print(f"caja -35 (buscada hasta -{a.up}): mejor {mejor35[2]} en "
          f"{mejor35[1]-a.up:+d} con {mejor35[0]:.2f} bits -> "
          f"{'débil, se omite' if mejor35[0] < 0.8 else 'presente'}")

    # una sola caja: mejor anchura en la región del promotor
    best = None
    for width in (6, 7, 8, 9, 10, 12):
        for lo in range(a.up - 20, a.up):
            if lo + width > a.up:
                continue
            box = build(windows, lo, lo + width, bg)
            if best is None or box["info_bits"] > best[0]["info_bits"]:
                best = (box, lo)
    box, lo = best
    print(f"caja única: {box['consensus']} en {lo-a.up:+d}..{lo-a.up+box['width']:+d} "
          f"({box['width']} pb, {box['info_bits']:.2f} bits, {box['n_sites']} TSS)")

    # comparación con el -10 de E. coli (6 pb alineados por el extremo 3')
    if ECOLI_PWM.exists():
        for r in csv.DictReader(open(ECOLI_PWM), delimiter="\t"):
            if r["sigma"] == "sigma70" and r["box"] == "-10":
                eco = json.loads(r["pwm"])
                mine = box["pwm"][:6]
                print(f"  correlación con el -10 de E. coli (6 pb): "
                      f"r = {corr(mine, eco):.3f}")

    row = {"sigma": "sigma70", "box": "-10", "n_sites": box["n_sites"],
           "width": box["width"], "consensus": box["consensus"],
           "info_bits": box["info_bits"], "pwm": json.dumps(box["pwm"])}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        w.writeheader()
        w.writerow(row)
    print(f"-> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
