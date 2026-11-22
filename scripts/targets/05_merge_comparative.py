#!/usr/bin/env python3
"""Añade co-ocurrencia 400 a los candidatos top."""
import csv, json
from collections import defaultdict, Counter
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=D/"data"/"metabolism"/"targets"; C=O/"comparative"; E=D/"data"/"metabolism"/"essential"
# presencia por genoma
pres=defaultdict(set); genomes=[]
for f in sorted(C.glob("*.tsv")):
    g=f.stem; genomes.append(g)
    for line in open(f,errors="ignore"):
        p=line.split("\t")
        if len(p)<5: continue
        try: pid,qcov=float(p[2]),float(p[3])
        except: continue
        if pid>=30 and qcov>=50: pres[p[1].strip()].add(g)
allg=set(genomes)
# genomas con gap por target
gap=defaultdict(set)
for r in csv.DictReader(open(E/"essential_gap_queue.tsv"),delimiter="\t"):
    gap[(r["map_id"],r["gap"])].add(r["genome"])
rows=list(csv.DictReader(open(O/"target_candidates.tsv"),delimiter="\t"))
out=[]
for r in rows:
    key=f"{r['genome']}|{r['protein_id']}"; present=pres.get(key,set())
    gs=gap.get((r["map_id"],r["target"]),set())
    ng=allg-gs
    pg=len(present & gs); pn=len(present & ng)
    rg=pg/len(gs) if gs else 0; rn=pn/len(ng) if ng else 0
    spec=round(rn-rg,3)
    sup="high" if spec>=0.2 and rn>=0.3 else "medium" if spec>0 else "low"
    out.append({**r,"present_gap":pg,"present_nongap":pn,"gap_genomes":len(gs),"nongap_genomes":len(ng),
                "rate_gap":round(rg,3),"rate_nongap":round(rn,3),"comparative_specificity":spec,
                "comparative_support":sup})
out.sort(key=lambda x:(x["target"], -float(x["score"]) - max(0,float(x["comparative_specificity"]))*0.2))
with open(O/"target_candidates_comparative.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(out)
json.dump({"genomes":len(genomes),"candidates":len(out),"support":dict(Counter(r["comparative_support"] for r in out))},open(O/"comparative_summary.json","w"),indent=2)
print({"genomes":len(genomes),"candidates":len(out),"support":dict(Counter(r["comparative_support"] for r in out))})
