#!/usr/bin/env python3
"""compare_structure.py — ¿aporta ViennaRNA sobre el proxy de horquilla?

Compara dos formas de puntuar un terminador intrínseco sobre **regiones
intergénicas reales** de K-12 (no sobre secuencias ya sabidas terminadoras, que
sería un test trivial):

  A) proxy heurístico: tallo GC (≤1 error) + tracto poli-T
  B) ΔG de ViennaRNA (RNAfold) sobre la misma ventana, más poli-T

Positivos: los 375 terminadores curados (`TerminatorSet.sequence`).
Negativos: regiones intergénicas aleatorias de K-12 con la misma distribución
de longitudes.

Se reporta AUC de cada puntuación y la precisión/recall con la regla clásica de
ARNold (ΔG ≤ −10 kcal/mol y poli-U ≥ 4).

Uso
---
    python scripts/regulation/compare_structure.py
"""
from __future__ import annotations

import argparse
import csv
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import find_operons as F  # noqa: E402

TERMS = ROOT / "data" / "regulondb" / "TerminatorSet.tsv"


def auc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    wins = sum(1 for p in pos for n in neg if p > n)
    ties = sum(1 for p in pos for n in neg if p == n)
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-neg", type=int, default=400)
    a = ap.parse_args()

    try:
        import RNA
    except Exception:
        print("falta ViennaRNA: pip install ViennaRNA", file=sys.stderr)
        return 1

    rows = []
    with open(TERMS, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = [c.split(")")[-1] for c in next(rd)]
        for r in rd:
            d = dict(zip(header, r))
            s = d.get("sequence", "").strip().upper()
            if s:
                rows.append(s)
    print(f"terminadores curados: {len(rows)}")

    seq, genes = F.read_gbff(ROOT / "data" / "test_ecoli.gbff")
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    inter = []
    for x, y in zip(genes, genes[1:]):
        lo, hi = min(x["end"], y["end"]), max(x["start"], y["start"])
        if hi - lo > 15:
            inter.append(seq[lo:hi])
    print(f"regiones intergénicas: {len(inter)}")

    rnd = random.Random(42)
    lengths = [len(s) for s in rows]
    negs = []
    for _ in range(a.n_neg):
        w = rnd.choice(inter)
        L = rnd.choice(lengths)
        if len(w) <= L:
            negs.append(w)
        else:
            start = rnd.randint(0, len(w) - L)
            negs.append(w[start:start + L])

    def dg(s: str) -> float:
        return -RNA.fold(s)[1]          # más positivo = más estable

    heur_pos = [F.terminator_scan(s)[0] for s in rows]
    heur_neg = [F.terminator_scan(s)[0] for s in negs]
    g_pos = [dg(s) for s in rows]
    g_neg = [dg(s) for s in negs]
    comb_pos = [h + 2 * g for h, g in zip(heur_pos, g_pos)]
    comb_neg = [h + 2 * g for h, g in zip(heur_neg, g_neg)]

    print("\nAUC (positivos = terminadores curados; negativos = intergénicas):")
    print(f"  proxy heurístico      {auc(heur_pos, heur_neg):.3f}")
    print(f"  ViennaRNA (-ΔG)       {auc(g_pos, g_neg):.3f}")
    print(f"  heurístico + 2·(-ΔG)  {auc(comb_pos, comb_neg):.3f}")
    print(f"\n  ΔG mediana: curados {statistics.median([-g for g in g_pos]):.1f} "
          f"kcal/mol | intergénicas {statistics.median([-g for g in g_neg]):.1f}")

    # regla clásica ARNold: ΔG <= -10 y poli-U >= 4
    tp = sum(1 for s, g in zip(rows, g_pos) if -g <= -10 and F.poly_t_score(s) >= 4)
    fp = sum(1 for s, g in zip(negs, g_neg) if -g <= -10 and F.poly_t_score(s) >= 4)
    print(f"\nregla ARNold (ΔG<=-10 y poliU>=4):")
    print(f"  sensibilidad {tp/len(rows):.1%} | falsos positivos {fp}/{len(negs)} "
          f"({fp/len(negs):.1%})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
