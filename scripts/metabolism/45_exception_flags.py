#!/usr/bin/env python3
"""Marca gaps con excepciones documentadas (aviso, no exención automática)."""
import csv, json
from collections import Counter
from pathlib import Path

D=Path(__file__).resolve().parents[2]
E=D/"data"/"metabolism"/"essential"
exc={}
for r in csv.DictReader(open(E/"essential_exceptions.tsv"),delimiter="\t"):
    exc.setdefault(r["map_id"],[]).append(r)
paths={ (r["genome"],r["gap"]):r.get("pathways","") for r in csv.DictReader(open(D/"data/metabolism/real_gaps.tsv"),delimiter="\t") }
rows=[]
for g in csv.DictReader(open(D/"data/metabolism/gap_interpretation.tsv"),delimiter="\t"):
    p=paths.get((g["genome"],g["gap"]),"")
    if not p: continue
    for mid,xs in exc.items():
        if mid in p:
            for x in xs:
                rows.append({"genome":g["genome"],"gap":g["gap"],"map_id":mid,
                             "exception":x["exception"],"example_taxon":x["example_taxon"],
                             "implication":x["implication"],"reference":x["reference"]})
with open(E/"gap_exception_flags.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
s={"rows":len(rows),"genomes":len({r['genome'] for r in rows}),"gaps":len({(r['genome'],r['gap']) for r in rows}),
   "by_map":dict(Counter(r["map_id"] for r in rows))}
json.dump(s,open(E/"exceptions_flags_summary.json","w"),indent=2)
print(json.dumps(s,indent=2))
