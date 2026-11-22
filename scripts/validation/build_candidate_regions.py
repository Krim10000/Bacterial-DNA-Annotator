#!/usr/bin/env python3
"""build_candidate_regions.py — regiones donde buscar TSS y TTS.

Un promotor o un terminador viven en el espacio intergénico, así que analizar
el genoma entero es desperdiciar 87% del trabajo. Pero recortar en seco también
falla: de los 375 terminadores curados de K-12, **48 (13%) están dentro de un
CDS** y 74% a menos de 50 pb de uno. Un recorte estricto los perdería.

Solución: partir del espacio intergénico, descartar las regiones diminutas y
**añadir un margen dentro de los genes vecinos**. Con 50 pb de margen el área de
búsqueda pasa de 12,5% a 19% del genoma —y cubre los terminadores solapados.

Salida: TSV `contig start end length kind`, listo para escanear o para trocear
y pasar a un detector externo.

Uso
---
    python scripts/validation/build_candidate_regions.py \
        --gbff data/test_ecoli.gbff --margin 50 --out data/regulondb/regions_k12.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FEATURES = ("CDS", "tRNA", "rRNA", "ncRNA")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gbff", required=True)
    ap.add_argument("--margin", type=int, default=50,
                    help="pb dentro de los genes vecinos (default 50)")
    ap.add_argument("--min-len", type=int, default=20,
                    help="descartar regiones intergénicas menores (default 20 pb)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    from Bio import SeqIO
    total = 0
    rows = []
    for rec in SeqIO.parse(a.gbff, "genbank"):
        L = len(rec.seq)
        total += L
        contig = rec.id
        feats = []
        for f in rec.features:
            if f.type in FEATURES:
                for p in f.location.parts:
                    feats.append((int(p.start), int(p.end)))
        feats.sort()
        merged: list[tuple[int, int]] = []
        for x, y in feats:
            if merged and x <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], y))
            else:
                merged.append((x, y))
        regs = []
        prev = 0
        for x, y in merged:
            if x > prev:
                regs.append((prev, x))
            prev = max(prev, y)
        if prev < L:
            regs.append((prev, L))
        for x, y in regs:
            if y - x < a.min_len:
                continue
            rows.append((contig, max(0, x - a.margin), min(L, y + a.margin),
                         y - x, "intergenic+margin"))

    rows.sort()
    OUT = Path(a.out)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    covered = 0
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig", "start", "end", "span", "kind"])
        cur_contig, cur_start, cur_end = None, 0, 0
        for c, s, e, span, kind in rows:
            if c == cur_contig and s <= cur_end:      # fusionar solapes del margen
                cur_end = max(cur_end, e)
                continue
            if cur_contig is not None:
                w.writerow([cur_contig, cur_start, cur_end, cur_end - cur_start, kind])
                covered += cur_end - cur_start
            cur_contig, cur_start, cur_end = c, s, e
        if cur_contig is not None:
            w.writerow([cur_contig, cur_start, cur_end, cur_end - cur_start, kind])
            covered += cur_end - cur_start

    print(f"genoma: {total:,} pb")
    print(f"regiones: {len(rows)} (>= {a.min_len} pb) con margen {a.margin} pb")
    print(f"cubierto: {covered:,} pb = {100*covered/total:.1f}% del genoma")
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
