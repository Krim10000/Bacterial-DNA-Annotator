#!/usr/bin/env python3
"""Sitios activos de la mejor referencia bacteriana por target esencial."""
import csv, glob, json, time, urllib.request
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

D=Path(__file__).resolve().parents[2]; E=D/"data"/"metabolism"/"essential"
info={}
for f in glob.glob(str(D/"data/metabolism/evidence/work*/acc_info.tsv")):
    for r in csv.DictReader(open(f),delimiter="\t"): info[r["acc"]]=r.get("euk","")
# candidate best_acc
best={}
for p in [D/"data/metabolism/targets/target_candidates_comparative.tsv", D/"data/metabolism/targets_p3/target_candidates.tsv"]:
    if p.exists():
        for r in csv.DictReader(open(p),delimiter="\t"):
            if r.get("best_acc") and (r["map_id"],r["target"]) not in best:
                if info.get(r["best_acc"]) == "0": best[(r["map_id"],r["target"])]=r["best_acc"]
# fallback refs
refs=defaultdict(list)
for p in [D/"data/metabolism/targets/refs.tsv", D/"data/metabolism/targets_p3/refs.tsv"]:
    if p.exists():
        for r in csv.DictReader(open(p),delimiter="\t"):
            if info.get(r["acc"])=="0": refs[r["target_id"].split("::")[0]+"::"+r["target_id"].split("::")[1]].append(r["acc"])
pri=list(csv.DictReader(open(E/"essential_priority_targets.tsv"),delimiter="\t"))
r2a=defaultdict(set)
for line in open(D/"data/metabolism/rhea2uniprot_sprot.tsv"):
    f=line.rstrip().split("\t")
    if len(f)>=4 and f[0]!="RHEA_ID": r2a["RHEA:"+f[0]].add(f[3])
refg={}
for r in csv.DictReader(open(D/"out_fl/refined/refined_gaps.tsv"),delimiter="\t"): refg[(r["genome"],r["gap"])]=r
q={}
for r in csv.DictReader(open(E/"essential_gap_queue.tsv"),delimiter="\t"): q.setdefault((r["map_id"],r["gap"]),r["genome"])
def fetch(acc):
    try:
        req=urllib.request.Request(f"https://rest.uniprot.org/uniprotkb/{acc}.json",headers={"User-Agent":"dna-poc/1.0"})
        return json.load(urllib.request.urlopen(req,timeout=20))
    except Exception:
        return {}
accs_to={}
for r in pri:
    k=(r["map_id"],r["target"]); a=best.get(k)
    if not a:
        x=refs.get(r["map_id"]+"::"+r["target"],[])
        a=x[0] if x else ""
    if not a:
        g=q.get(k)
        if g:
            ps=[v for v in refg.get((g,r["target"]),{}).get("producers","").split(";") if v]
            cs=[]
            for rid in ps: cs += [z for z in r2a.get(rid,set()) if info.get(z)=="0"]
            a=cs[0] if cs else ""
    if a: accs_to[k]=a
uniq=sorted(set(accs_to.values()))
print("targets",len(pri),"with_acc",len(accs_to),"unique_acc",len(uniq),flush=True)
cache={}
with ThreadPoolExecutor(max_workers=4) as ex:
    for acc,d in zip(uniq, ex.map(fetch,uniq)): cache[acc]=d
rows=[]
for r in pri:
    k=(r["map_id"],r["target"]); acc=accs_to.get(k,""); d=cache.get(acc,{})
    name=((d.get("proteinDescription",{}) or {}).get("recommendedName",{}) or {}).get("fullName",{}).get("value","")
    ecs=((d.get("proteinDescription",{}) or {}).get("recommendedName",{}) or {}).get("ecNumbers",[]) or []
    ec=";".join(x.get("value","") for x in ecs)
    feats=[f for f in (d.get("features",[]) or []) if f.get("type") in ("Active site","Binding site","Metal binding","Site")]
    if not feats:
        rows.append({"map_id":r["map_id"],"pathway":r["pathway"],"target":r["target"],"priority":r["priority"],
                     "acc":acc,"protein":name,"ec":ec,"organism":d.get("organism",{}).get("scientificName",""),
                     "feature":"none","description":"sin anotación de sitio en UniProt","start":"","end":"",
                     "exception":r.get("known_exception","")})
    for f in feats:
        loc=f.get("location",{}) or {}
        rows.append({"map_id":r["map_id"],"pathway":r["pathway"],"target":r["target"],"priority":r["priority"],
                     "acc":acc,"protein":name,"ec":ec,"organism":d.get("organism",{}).get("scientificName",""),
                     "feature":f.get("type",""),"description":f.get("description",""),
                     "start":(loc.get("start",{}) or {}).get("value",""),"end":(loc.get("end",{}) or {}).get("value",""),
                     "exception":r.get("known_exception","")})
with open(E/"ACTIVE_SITES_ESSENTIAL.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(E/"ACTIVE_SITES_ESSENTIAL.md","w") as fh:
    fh.write("# Sitios activos de references esenciales\n\n")
    fh.write("Una mejor referencia bacteriana por target; residuos del homólogo, no del ORF.\n\n")
    fh.write("| target | referencia | proteína | EC | features |\n|---|---|---|---|---:|\n")
    by=defaultdict(list)
    for x in rows: by[(x["map_id"],x["target"])].append(x)
    for (m,t),rs in by.items():
        r=rs[0]; fh.write(f"| {t} | {r['acc'] or '—'} | {r['protein'] or '—'} | {r['ec'] or '—'} | {len([x for x in rs if x['feature']!='none'])} |\n")
    fh.write("\n## Límites\n\n- Una sola referencia por target; no todas las isoformas.\n- UniProt no siempre anota el sitio activo.\n- Transferir solo por alineación estructural y conservación.\n")
s={"targets":len(pri),"with_ref":len(accs_to),"with_features":len({(x['map_id'],x['target']) for x in rows if x['feature']!='none'}),"rows":len(rows)}
json.dump(s,open(E/"ACTIVE_SITES_ESSENTIAL.json","w"),indent=2)
print(s)
