#!/usr/bin/env python3
"""robustness.py — ¿cuánto valen los AUC de comparar_structure.py?

La comparación previa (proxy 0.852 vs ViennaRNA 0.829) usó **un genoma**, **un
set curado** y **un esquema de negativos**. Antes de concluir nada hay que ver:

  1. **Incertidumbre**: bootstrap sobre positivos y negativos.
  2. **Sensibilidad al negativo**: intergénicas, codificantes, barajadas y de
     otro genoma.
  3. **¿Porta el umbral?**: la misma puntuación sobre intergénicas de 8 genomas
     de filos distintos. Si la tasa implícita cambia mucho, el corte fijado en
     E. coli no vale para el genoma 401.

Uso
---
    python scripts/regulation/robustness.py
"""
from __future__ import annotations

import csv
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import find_operons as F  # noqa: E402

TERMS = ROOT / "data" / "regulondb" / "TerminatorSet.tsv"
CUTOFF = 6.0


def auc(pos: list[float], neg: list[float]) -> float:
    if not pos or not neg:
        return float("nan")
    wins = sum(1 for p in pos for n in neg if p > n)
    ties = sum(1 for p in pos for n in neg if p == n)
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def auc_fast(pos: list[float], neg: list[float]) -> float:
    """AUC por rangos (evita el doble bucle en bootstrap)."""
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks = {}
    i = 0
    rank_sum_pos = 0.0
    n_pos, n_neg = len(pos), len(neg)
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        avg = (i + j + 1) / 2
        for k in range(i, j):
            if allv[k][1] == 1:
                rank_sum_pos += avg
        i = j
    return (rank_sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def load_terms() -> list[str]:
    rows = []
    with open(TERMS, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = [c.split(")")[-1] for c in next(rd)]
        for r in rd:
            d = dict(zip(header, r))
            s = d.get("sequence", "").strip().upper()
            if s:
                rows.append(s)
    return rows


def intergenic(genome_fasta: Path, gff: Path, limit: int = 400) -> list[str]:
    seqs = F.read_fasta(genome_fasta)
    genes = F.read_gff_genes(gff)
    out = []
    for a, b in F.pairs(F.ordered(genes)):
        s = seqs.get(a.get("seqid", ""), "")
        if not s:
            continue
        lo, hi = a["end"], b["start"]
        if hi - lo > 20:
            out.append(s[lo:hi])
        if len(out) >= limit:
            break
    return out


def main() -> int:
    rnd = random.Random(7)
    terms = load_terms()
    seq, genes = F.read_gbff(ROOT / "data" / "test_ecoli.gbff")
    genes = F.ordered([g for g in genes if g.get("locus_tag")])
    inter = [seq[a["end"]:b["start"]] for a, b in F.pairs(genes)
             if b["start"] - a["end"] > 20]
    coding = [seq[g["start"]:g["end"]] for g in genes[:600]
              if g["end"] - g["start"] in range(30, 400)]
    print(f"terminadores {len(terms)} | intergénicas K-12 {len(inter)} | codificantes {len(coding)}")

    heur_pos = [F.terminator_scan(s)[0] for s in terms]
    # puntuaciones cacheadas: el bootstrap remuestrea, no re-escanea
    inter_scores = [F.terminator_scan(s)[0] for s in inter]
    coding_scores = [F.terminator_scan(s)[0] for s in coding]

    # --- 1. bootstrap
    boots = []
    for _ in range(200):
        p = [rnd.choice(heur_pos) for _ in heur_pos]
        n = [rnd.choice(inter_scores) for _ in range(300)]
        boots.append(auc_fast(p, n))
    boots.sort()
    print(f"\n1) incertidumbre del proxy (200 bootstraps):")
    print(f"   AUC mediana {statistics.median(boots):.3f} "
          f"IC95% [{boots[5]:.3f}, {boots[194]:.3f}]")

    # --- 2. sensibilidad al conjunto negativo
    print("\n2) AUC del proxy según el negativo:")
    shuffled = ["".join(rnd.sample(list(s), len(s)))
                for s in (rnd.choice(terms) for _ in range(150))]
    sets = {
        "intergénicas K-12": inter_scores,
        "codificantes K-12": coding_scores,
        "barajadas (mono)": [F.terminator_scan(s)[0] for s in shuffled],
    }
    for name, neg in sets.items():
        samp = [rnd.choice(neg) for _ in range(min(300, len(neg)))]
        print(f"   {name:22s} AUC {auc_fast(heur_pos, samp):.3f} | "
              f"mediana neg {statistics.median(samp):.1f}")

    # --- 3. ¿porta el umbral a otros genomas?
    print(f"\n3) la misma puntuación en intergénicas de otros genomas "
          f"(corte E. coli = {CUTOFF}):")
    print(f"   {'genoma':32s} {'n':>4s} {'mediana':>8s} {'%>corte':>8s}")
    gffs = sorted((ROOT / "data" / "annotations_ab_initio").glob("*/annotation.gff3"))[:8]
    for gff in gffs:
        g = gff.parent.name
        fa = ROOT / "data" / "raw" / f"{g}_genomic.fna.gz"
        if not fa.exists():
            continue
        wins = intergenic(fa, gff, limit=150)
        if not wins:
            continue
        sc = [F.terminator_scan(w)[0] for w in wins]
        print(f"   {g[:32]:32s} {len(sc):4d} {statistics.median(sc):8.1f} "
              f"{100 * sum(1 for x in sc if x >= CUTOFF) / len(sc):7.0f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
