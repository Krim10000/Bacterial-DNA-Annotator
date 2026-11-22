#!/usr/bin/env python3
"""Prepara targets P1/P2, genomas con gap y referencias Rhea/UniProt."""
import csv, json, os
from collections import defaultdict
from pathlib import Path

D=Path(__file__).resolve().parents[2]; OUT=Path(os.environ.get("TARGET_DIR", D/"data"/"metabolism"/"targets")); OUT.mkdir(parents=True, exist_ok=True); E=D/"data"/"metabolism"/"essential"
PRI=set(os.environ.get("TARGET_PRIORITIES","P1,P2").split(","))
pri=[r for r in csv.DictReader(open(E/"essential_priority_targets.tsv"),delimiter="\t") if r["priority"] in PRI]
targets=[]; by_key={}
for r in pri:
    tid=f"{r['map_id']}::{r['target']}"
    if tid in by_key: continue
    x={"target_id":tid,"map_id":r["map_id"],"pathway":r["pathway"],"target":r["target"],
       "tier":r["tier"],"priority":r["priority"],"rhea":r["rhea"],"ec":r["ec"],
       "known_exception":r.get("known_exception",""),"example_taxon":r.get("example_taxon","")}
    targets.append(x); by_key[(r["map_id"],r["target"])]=x
with open(OUT/"targets.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(targets[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(targets)
# genomas con gap en targets P1/P2
q=list(csv.DictReader(open(E/"essential_gap_queue.tsv"),delimiter="\t"))
pairs=[]
for r in q:
    t=by_key.get((r["map_id"],r["gap"]))
    if t: pairs.append({"target_id":t["target_id"],"map_id":t["map_id"],"target":t["target"],"genome":r["genome"],"support":r["support"]})
# unique
seen=set(); ug=[]
for r in pairs:
    k=(r["target_id"],r["genome"])
    if k in seen: continue
    seen.add(k); ug.append(r)
with open(OUT/"target_genomes.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(ug[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(ug)
# refs
r2a=defaultdict(set)
for line in open(D/"data/metabolism/rhea2uniprot_sprot.tsv"):
    f=line.rstrip().split("\t")
    if len(f)>=4 and f[0]!="RHEA_ID": r2a["RHEA:"+f[0]].add(f[3])
pf=json.load(open(D/"data/metabolism/uniprot_pfam_cache.json"))
refs={}
for t in targets:
    for rid in t["rhea"].split(";"):
        for acc in r2a.get(rid,set()):
            key=(t["target_id"],acc)
            refs.setdefault(key,{"target_id":t["target_id"],"acc":acc,"rhea":rid,"domains":";".join(pf.get(acc,[]))})
with open(OUT/"refs.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=["target_id","acc","rhea","domains"],delimiter="\t"); w.writeheader(); w.writerows(refs.values())
# fasta refs
want={v["acc"] for v in refs.values()}; out=open(OUT/"refs.faa","w"); keep=False; n=0
for line in open(D/"data/db/uniprot_sprot.fasta"):
    if line.startswith(">"):
        acc=line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
        keep=acc in want
        if keep:
            out.write(f">{acc}\n"); n+=1
    elif keep: out.write(line)
out.close()
# dominios
doms=sorted({d for v in refs.values() for d in v["domains"].split(";") if d})
with open(OUT/"domains.tsv","w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["domain"]); [w.writerow([d]) for d in doms]
print({"targets":len(targets),"target_genomes":len(ug),"genomes":len({r['genome'] for r in ug}),"refs":len(refs),"ref_fasta":n,"domains":len(doms)})
