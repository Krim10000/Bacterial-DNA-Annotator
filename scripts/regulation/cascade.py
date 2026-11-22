#!/usr/bin/env python3
"""cascade.py — descarte por capas: ¿hasta dónde llega la precisión alta?

Idea (la del usuario): en vez de un umbral único que mezcla todo, partir las
fronteras en capas por confianza y reportar **precisión y cobertura por capa**.

Las capas son dos, porque las decisiones son dos:

    capa UNIR   probabilidad baja  de corte  -> "estos genes van juntos"
    capa CORTAR probabilidad alta  de corte  -> "aquí empieza otro operón"
    capa REVISAR lo demás                    -> sin evidencia suficiente

La capa de unión es la que importa para reconstruir operones: su precisión es
la fracción de pares que de verdad están en el mismo operón. Se busca el umbral
más permisivo que mantenga precisión >= objetivo (0.95 por defecto), y se
reporta cuánta cobertura se gana o se pierde.

Uso
---
    python scripts/regulation/cascade.py
    python scripts/regulation/cascade.py --target 0.98
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROBS = ROOT / "data" / "operon_finder" / "ecoli_k12_boundary_probs.tsv"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--target", type=float, default=0.95)
    ap.add_argument("--probs", type=Path, default=PROBS)
    a = ap.parse_args()

    if not a.probs.exists():
        print(f"falta {a.probs}; corre find_operons.py --train primero",
              file=sys.stderr)
        return 1

    rows = []
    for r in csv.DictReader(open(a.probs), delimiter="\t"):
        if r["label"] == "":
            continue
        rows.append((float(r["prob"]), int(r["label"])))
    if not rows:
        print("sin fronteras etiquetadas")
        return 1

    n = len(rows)
    n_cut = sum(lab for _, lab in rows)
    n_join = n - n_cut
    print(f"fronteras etiquetadas: {n}  (cortes {n_cut}, uniones {n_join})")

    # curva de fiabilidad por tramos de probabilidad
    print("\nfiabilidad por tramo de P(corte):")
    print(f"  {'P(corte)':>12s} {'n':>6s} {'cortes':>7s} {'P(corte|tramo)':>15s}")
    bins = [0, .02, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .98, 1.01]
    for lo, hi in zip(bins, bins[1:]):
        sub = [lab for p, lab in rows if lo <= p < hi]
        if not sub:
            continue
        print(f"  [{lo:.2f},{hi:.2f}) {len(sub):6d} {sum(sub):7d} "
              f"{sum(sub)/len(sub):15.3f}")

    # capa UNIR: los de P más baja, buscando el mayor conjunto con precisión >= target
    ordered = sorted(rows, key=lambda x: x[0])
    view = sol = 0
    t_join = 0.0
    for p, lab in ordered:
        view += 1
        sol += (1 - lab)
        if sol / view >= a.target:
            t_join = p
    join = [(p, lab) for p, lab in ordered if p <= t_join]
    j_prec = sum(1 - lab for _, lab in join) / len(join) if join else 0
    j_cov = len(join) / n_join if n_join else 0

    # capa CORTAR: simétrico por arriba
    rev = sorted(rows, key=lambda x: -x[0])
    view = sol = 0
    t_cut = 1.0
    for p, lab in rev:
        view += 1
        sol += lab
        if sol / view >= a.target:
            t_cut = p
    cut = [(p, lab) for p, lab in rev if p >= t_cut]
    c_prec = sum(lab for _, lab in cut) / len(cut) if cut else 0
    c_cov = len(cut) / n_cut if n_cut else 0

    print(f"\ncascada con objetivo de precisión >= {a.target:.2f}:")
    print(f"  UNIR    P(corte) <= {t_join:.3f}  n={len(join):5d}  "
          f"precisión {j_prec:.3f}  cobertura de uniones {j_cov:.1%}")
    print(f"  CORTAR  P(corte) >= {t_cut:.3f}  n={len(cut):5d}  "
          f"precisión {c_prec:.3f}  cobertura de cortes {c_cov:.1%}")
    review = n - len(join) - len(cut)
    print(f"  REVISAR resto                        n={review:5d}  "
          f"({100 * review / n:.1f}% de las fronteras)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
