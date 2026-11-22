#!/usr/bin/env python3
"""Construye el panel esencial curado a partir de la literatura."""
import csv, json
from collections import Counter
from pathlib import Path

D=Path(__file__).resolve().parents[2]
E=D/"data"/"metabolism"/"essential"; R=E/"research"
rec=list(csv.DictReader(open(R/"panel_recommendations.tsv"),delimiter="\t"))
ev={r["map_id"]:r for r in csv.DictReader(open(R/"pathway_evidence.tsv"),delimiter="\t")}
cur={r["map_id"]:r for r in csv.DictReader(open(E/"essential_pathways.tsv"),delimiter="\t")}
exc={}
for r in csv.DictReader(open(R/"exceptions.tsv"),delimiter="\t"):
    exc.setdefault(r["map_id"],[]).append(r["exception"])
rows=[]
for r in rec:
    if r["action"]=="remove": continue
    mid=r["map_id"]; e=ev.get(mid,{}); c=cur.get(mid,{})
    tier=r["recommended_tier"] if r["recommended_tier"] in ("core","conditional") else c.get("tier","conditional")
    medium="defined_mineral" if tier=="core" else "conditional"
    if mid=="map00190": medium="aerobic_medium"
    if mid=="map00630": medium="acetate/C2"
    if mid=="map00480": medium="redox/lineage"
    if mid=="map00330": bvd="proline/polyamine; catabolism"
    elif mid=="map00190": bvd="energy"
    elif mid=="map00630": bvd="anaplerosis"
    elif mid=="map00480": bvd="redox"
    elif mid in ("map00550","map00540","map00520","map00061","map00564"): bvd="envelope"
    else: bvd="biosynthesis"
    rows.append({
        "map_id":mid,"pathway":e.get("pathway") or c.get("pathway") or r["pathway"],
        "tier":tier,"rationale":e.get("notes") or c.get("rationale",""),
        "applicability":"vida libre, medio mínimo; condicional por filo/estilo de vida",
        "medium_dependence":medium,"lineage_scope":"bacteria",
        "biosynthesis_vs_degradation":bvd,
        "evidence":e.get("evidence",""),"reference":e.get("reference",""),
        "exception":e.get("exception") or "; ".join(exc.get(mid,[])),
        "confidence":e.get("confidence","medium"),"notes":e.get("notes",""),
    })
with open(E/"essential_pathways.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(E/"essential_exceptions.tsv","w",newline="") as fh:
    src=list(csv.DictReader(open(R/"exceptions.tsv"),delimiter="\t"))
    w=csv.DictWriter(fh,fieldnames=list(src[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(src)
# completitud
summ={r["pathway"].split(" ")[0]:r for r in csv.DictReader(open(D/"data/metabolism/pathway_completeness_summary.tsv"),delimiter="\t")}
with open(E/"essential_completeness.tsv","w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["map_id","pathway","tier","mean_coverage","mean_missing","gaps","essential_gaps"])
    for x in rows:
        s=summ.get(x["map_id"],{}); w.writerow([x["map_id"],x["pathway"],x["tier"],s.get("mean_coverage",""),s.get("mean_missing",""),s.get("gaps",""),s.get("essential_gaps","")])
# cola
review={(r["genome"],r["gap"]):r for r in csv.DictReader(open(D/"data/metabolism/review/gap_review_current.tsv"),delimiter="\t")}
q=[]
for g in csv.DictReader(open(D/"data/metabolism/real_gaps.tsv"),delimiter="\t"):
    for x in rows:
        if x["map_id"] in g.get("pathways",""):
            rv=review.get((g["genome"],g["gap"]),{})
            q.append({"genome":g["genome"],"gap":g["gap"],"map_id":x["map_id"],"pathway":x["pathway"],
                      "tier":x["tier"],"support":g.get("support_400",""),"transport":g.get("transport_status",""),
                      "classification":rv.get("classification","unknown"),"action":rv.get("action",""),
                      "cause":rv.get("evidence","")[:200]})
with open(E/"essential_gap_queue.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(q[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(q)
s={"pathways":len(rows),"core":sum(1 for x in rows if x["tier"]=="core"),"conditional":sum(1 for x in rows if x["tier"]=="conditional"),
   "gaps_in_essential":len(q),"by_tier":dict(Counter(x["tier"] for x in q)),"by_map":dict(Counter(x["map_id"] for x in q))}
json.dump(s,open(E/"summary.json","w"),indent=2)
print(json.dumps(s,indent=2))
