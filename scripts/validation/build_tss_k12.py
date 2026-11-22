#!/usr/bin/env python3
"""build_tss_k12.py — TSS curados de K-12 en un TSV simple.

RegulonDB ya los trae en `PromoterSet.tsv` (4.056 promotores), pero con nombres
de columna incómodos y hebra como forward/reverse. Esto lo normaliza a
`pos, strand, sigma, first_gene, evidence` para que `rule_operons.py --tss`
lo consuma directamente.

Salida: data/regulondb/tss_k12.tsv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "regulondb" / "PromoterSet.tsv"
OUT = ROOT / "data" / "regulondb" / "tss_k12.tsv"


def main() -> int:
    if not SRC.exists():
        print(f"falta {SRC}")
        return 1
    n = 0
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "strand", "sigma", "first_gene", "evidence"])
        for r in csv.DictReader(open(SRC), delimiter="\t"):
            pos = (r.get("4)posTSS") or "").strip()
            st = (r.get("3)strand") or "").strip()
            if not pos.isdigit():
                continue
            w.writerow([int(pos), "-" if st == "reverse" else "+",
                        (r.get("5)sigmaFactor") or "").strip(),
                        (r.get("7)firstGeneName") or "").strip(),
                        (r.get("13)promoterEvidence") or "").strip()])
            n += 1
    print(f"{n} TSS -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
