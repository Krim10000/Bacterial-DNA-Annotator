#!/usr/bin/env python3
"""terminator_classes.py — desagregar el terminador por clase estructural.

El objetivo no es un AUC global (que mezcla casos fáciles y difíciles) sino
saber **qué identificamos bien, qué a medias y qué mal**, considerando solo la
estructura de terminación.

Componentes que se miden por secuencia
--------------------------------------
  stem      longitud del tallo (pares complementarios, <=1 error)
  gc        fracción GC del tallo
  loop      longitud del bucle
  polyU     tracto de T/U más largo
  dG        energía libre del MFE (ViennaRNA)
  diversidad  ambigüedad conformacional del ensemble (ViennaRNA)

Clases (reglas explícitas, no ajustadas al resultado)
----------------------------------------------------
  A  fuerte   stem>=6 y gc>=0.60 y polyU>=4  -> "lo identificamos bien"
  B  media    stem>=5 y polyU>=3             -> "más o menos"
  C  débil    el resto                        -> "mal / invisible"

El contraste es contra intergénicas reales de K-12, no contra secuencias
aleatorias, porque el genoma está lleno de horquillas incidentales.

Uso
---
    python scripts/regulation/terminator_classes.py
    python scripts/regulation/terminator_classes.py --neg 800
"""
from __future__ import annotations

import argparse
import csv
import random
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import find_operons as F  # noqa: E402

TERMS = ROOT / "data" / "regulondb" / "TerminatorSet.tsv"
BASES = "ACGT"


def parse_structure(seq: str, mismatches: int = 1) -> dict:
    """Devuelve stem, gc, loop, polyU y posición del mejor candidato."""
    seq = seq.upper()
    n = len(seq)
    best = {"stem": 0, "gc": 0.0, "loop": 0, "start": 0}
    for i in range(n):
        for stem in range(min(12, n - i), 3, -1):
            head = seq[i:i + stem]
            if any(b not in BASES for b in head):
                continue
            for loop in range(2, 11):
                j = i + stem + loop
                if j + stem > n:
                    continue
                tail = seq[j:j + stem]
                if any(b not in BASES for b in tail):
                    continue
                mm = sum(1 for x, y in zip(head, F.revcomp(tail)) if x != y)
                if mm > mismatches:
                    continue
                gc = sum(1 for b in head if b in "GC") / stem
                cand = {"stem": stem, "gc": gc, "loop": loop, "start": i}
                if (cand["stem"] + 3 * cand["gc"]) > (best["stem"] + 3 * best["gc"]):
                    best = cand
    run = 0
    best_run = 0
    for b in seq + "A":
        run = run + 1 if b == "T" else 0
        best_run = max(best_run, run)
    best["polyU"] = best_run
    return best


def classify(s: dict) -> str:
    if s["stem"] >= 6 and s["gc"] >= 0.60 and s["polyU"] >= 4:
        return "A_fuerte"
    if s["stem"] >= 5 and s["polyU"] >= 3:
        return "B_media"
    return "C_debil"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--neg", type=int, default=600)
    a = ap.parse_args()

    try:
        import RNA
    except Exception:
        print("falta ViennaRNA", file=sys.stderr)
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
    tipos = Counter()
    with open(TERMS, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = [c.split(")")[-1] for c in next(rd)]
        for r in rd:
            d = dict(zip(header, r))
            tipos[d.get("type", "")] += 1
    print("tipos en el catálogo:", tipos.most_common())

    seq, genes = F.read_gbff(ROOT / "data" / "test_ecoli.gbff")
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    inter = [seq[x["end"]:y["start"]] for x, y in F.pairs(genes)
             if y["start"] - x["end"] > 40]
    rnd = random.Random(11)
    L = int(statistics.median(len(s) for s in rows))
    negs = []
    for _ in range(a.neg):
        w = rnd.choice(inter)
        st = rnd.randint(0, max(0, len(w) - L))
        negs.append(w[st:st + L])
    print(f"negativos: {len(negs)} ventanas de {L} nt de intergénicas reales")

    def measure(seqs):
        out = []
        for s in seqs:
            st = parse_structure(s)
            mfe_struct, mfe = RNA.fold(s)
            st["dG"] = -mfe
            try:
                # ambigüedad = cuánta energía se gana al permitir el ensemble
                _struct, e_ens = RNA.pf_fold(s)
                st["div"] = round(mfe - e_ens, 3)   # >=0; mayor = más alternativas
            except Exception:
                st["div"] = float("nan")
            st["cls"] = classify(st)
            out.append(st)
        return out

    pos = measure(rows)
    neg = measure(negs)

    print(f"\n{'clase':10s} {'pos':>5s} {'neg':>5s} {'%pos':>6s} {'%neg':>6s} {'enriquec.':>10s}")
    pc = Counter(x["cls"] for x in pos)
    nc = Counter(x["cls"] for x in neg)
    for cls in ("A_fuerte", "B_media", "C_debil"):
        fp = pc[cls] / len(pos)
        fn = nc[cls] / len(neg)
        enr = (fp / fn) if fn else float("inf")
        print(f"{cls:10s} {pc[cls]:5d} {nc[cls]:5d} {100*fp:5.0f}% {100*fn:5.0f}% {enr:9.1f}x")

    # matriz de confusión práctica sobre una selección de clase A
    for cls in ("A_fuerte", "B_media"):
        tp = sum(1 for x in neg if x["cls"] == cls)
        print(f"\nsi aceptamos solo '{cls}': {pc[cls]} verdaderos de {len(pos)}, "
              f"y se cuelan {tp} de {len(neg)} intergénicas")

    print("\nambigüedad conformacional (diversidad del ensemble):")
    for cls in ("A_fuerte", "B_media", "C_debil"):
        sub = [x["div"] for x in pos if x["cls"] == cls and x["div"] == x["div"]]
        sn = [x["div"] for x in neg if x["cls"] == cls and x["div"] == x["div"]]
        if sub:
            print(f"  {cls:10s} curados mediana {statistics.median(sub):.2f} | "
                  f"intergénicas {statistics.median(sn) if sn else float('nan'):.2f}")

    print("\nseñales por clase (medianas):")
    for cls in ("A_fuerte", "B_media", "C_debil"):
        sub = [x for x in pos if x["cls"] == cls]
        if not sub:
            continue
        print(f"  {cls:10s} stem {statistics.median([x['stem'] for x in sub]):.0f} | "
              f"gc {statistics.median([x['gc'] for x in sub]):.2f} | "
              f"polyU {statistics.median([x['polyU'] for x in sub]):.0f} | "
              f"dG {statistics.median([x['dG'] for x in sub]):.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
