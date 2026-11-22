#!/usr/bin/env python3
"""Valida la lista piloto y escribe data/membrane/pilot_genomes.tsv con metadatos.

La lista curada vive en scripts/membrane/pilot_genomes.tsv (10 filoes distintos).
Este paso no anota: solo comprueba que existan GFF/proteina/DIAMOND y agrega
tamano, contigs y CDS para que el resumen sea trazable.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def read_pilot(path: Path) -> list[dict]:
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                sys.exit(f"[select] linea invalida: {line!r}")
            rows.append(dict(genome=f[0], phylum=f[1], organism=f[2], rationale=f[3]))
    if len(rows) != 10:
        sys.exit(f"[select] se esperaban 10 genomas, hay {len(rows)}")
    if len({r['genome'] for r in rows}) != 10:
        sys.exit("[select] genomas duplicados")
    return rows


def gff_stats(path: Path):
    contigs = {}
    cds = 0
    for line in open(path):
        if line.startswith("##sequence-region"):
            p = line.split()
            contigs[p[1]] = int(p[3])
        elif not line.startswith("#") and "\tCDS\t" in line:
            cds += 1
    return contigs, cds


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default=str(HERE / "pilot_genomes.tsv"))
    ap.add_argument("--annot-root", default=str(ROOT / "data/annotations_ab_initio"))
    ap.add_argument("--out", default=str(ROOT / "data/membrane/pilot_genomes.tsv"))
    args = ap.parse_args()

    rows = read_pilot(Path(args.pilot))
    out = []
    errors = []
    for r in rows:
        d = Path(args.annot_root) / r["genome"]
        need = ["annotation.gff3", "proteins.faa", "functional.tsv"]
        missing = [n for n in need if not (d / n).exists()]
        if missing:
            errors.append(f"{r['genome']}: faltan {missing}")
            continue
        contigs, cds = gff_stats(d / "annotation.gff3")
        nprot = sum(1 for line in open(d / "proteins.faa") if line.startswith(">"))
        out.append(dict(**r, contigs=len(contigs), bp=sum(contigs.values()),
                        cds=cds, proteins=nprot))
    if errors:
        sys.exit("[select] " + "; ".join(errors))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["genome", "phylum", "organism", "contigs", "bp",
                                       "cds", "proteins", "rationale"])
        w.writeheader()
        w.writerows(out)
    print(f"[select] {len(out)} genomas -> {args.out}")
    for r in out:
        print(f"  {r['genome']}\t{r['phylum']}\t{r['organism']}\t{r['bp']} bp\t{r['cds']} CDS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
