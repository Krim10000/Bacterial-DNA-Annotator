#!/usr/bin/env python3
"""Descarga AFDB de los homólogos de STRUCTURAL_CANDIDATES (sin GPU)."""
import csv, json, urllib.request
from pathlib import Path
D=Path(__file__).resolve().parents[2]; T=D/"data"/"metabolism"/"targets"; O=T/"structures"/"afdb"; O.mkdir(parents=True, exist_ok=True)
rows=list(csv.DictReader(open(T/"STRUCTURAL_CANDIDATES.tsv"),delimiter="\t"))
for r in rows:
    acc=r.get("best_acc","").strip()
    if not acc: continue
    out=O/f"{acc}.pdb"
    if out.exists(): continue
    try:
        meta=json.load(urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",timeout=20))
        url=meta[0]["pdbUrl"]
        out.write_bytes(urllib.request.urlopen(url,timeout=30).read())
        print("ok",acc)
    except Exception as e:
        print("fail",acc,str(e)[:80])
