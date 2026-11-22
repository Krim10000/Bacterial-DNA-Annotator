#!/usr/bin/env python3
"""Extrae los perfiles HMM de los targets esenciales para su búsqueda."""
import csv, os
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=Path(os.environ.get("TARGET_DIR", D/"data"/"metabolism"/"targets"))
want={r["domain"] for r in csv.DictReader(open(O/"domains.tsv"),delimiter="\t")}
out=open(O/"panel.hmm","w"); block=[]; name=None; n=0
for line in open(D/".cache/pfam/Pfam-A.hmm",errors="ignore"):
    if line.startswith("HMMER3"): block=[line]; name=None
    elif block:
        block.append(line)
        if line.startswith("NAME"): name=line.split(None,1)[1].strip()
        if line.startswith("//"):
            if name in want: out.writelines(block); n+=1
            block=[]
out.close(); print("HMMs",n)
