#!/usr/bin/env python3
"""Selecciona los ORFs no anotados de STRUCTURAL_CANDIDATES para plegado."""
import csv
from pathlib import Path
D=Path(__file__).resolve().parents[2]; T=D/"data"/"metabolism"/"targets"
rows=list(csv.DictReader(open(T/"STRUCTURAL_CANDIDATES.tsv"),delimiter="\t"))
want={(r["genome"],r["protein_id"]) for r in rows}
seq={}
for g in {x[0] for x in want}:
    hdr=None; buf=[]
    for line in open(D/"data"/"annotations_ab_initio"/g/"proteins.faa"):
        if line.startswith(">"):
            if hdr and (g,hdr) in want: seq[(g,hdr)]="".join(buf).replace("*","")
            hdr=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if hdr and (g,hdr) in want: seq[(g,hdr)]="".join(buf).replace("*","")
with open(T/"structural_candidates.faa","w") as fh:
    for (g,p),s in seq.items(): fh.write(f">{g}|{p}\n{s}\n")
print("selected",len(seq))
