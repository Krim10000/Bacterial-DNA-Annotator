#!/usr/bin/env python3
"""Calibración data-driven del panel esencial (no sustituye literatura)."""
import csv, json, statistics
from collections import Counter, defaultdict
from pathlib import Path

D = Path(__file__).resolve().parents[2]
E = D/"data"/"metabolism"/"essential"
pc = list(csv.DictReader(open(D/"data/metabolism/pathway_completeness.tsv"), delimiter="\t"))
by = defaultdict(list)
for r in pc:
    by[r["pathway"].split(" ")[0]].append(r)
paths = list(csv.DictReader(open(E/"essential_pathways.tsv"), delimiter="\t"))
gaps = list(csv.DictReader(open(D/"data/metabolism/real_gaps.tsv"), delimiter="\t"))
review = {(r["genome"], r["gap"]): r for r in csv.DictReader(open(D/"data/metabolism/review/gap_review_current.tsv"), delimiter="\t")}
cal=[]; outliers=[]
for p in paths:
    mid=p["map_id"]; rs=by.get(mid,[])
    cov=sorted(float(r["coverage"]) for r in rs) if rs else [0]
    med=statistics.median(cov); p10=cov[max(0,int(len(cov)*0.1)-1)]
    gp=[g for g in gaps if mid in g.get("pathways","")]
    strong=[g for g in gp if int(g.get("support_400",0) or 0)>=100 and g.get("transport_status")=="no_transport"]
    w=1.0 if p["tier"]=="core" else 0.2
    score=w*(len(strong)+2*len(gp))/max(1,sum(int(r.get("expected",0) or 0) for r in rs)/max(1,len(rs)))
    cal.append({**p,"genomes":len(rs),"mean_coverage":round(sum(cov)/len(cov),4),
                "median_coverage":round(med,4),"p10_coverage":round(p10,4),
                "gaps":len(gp),"strong_gaps":len(strong),"priority_score":round(score,4)})
    for r in rs:
        c=float(r["coverage"])
        if c < p10:
            outliers.append({"genome":r["genome"],"map_id":mid,"pathway":p["pathway"],
                             "coverage":c,"p10":round(p10,4),"median":round(med,4),
                             "gaps":r.get("gaps",""),"tier":p["tier"]})
# ranks
cal.sort(key=lambda x:(-x["priority_score"], -x["gaps"]))
for i,x in enumerate(cal,1): x["priority_rank"]=i
with open(E/"essential_calibration.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(cal[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(cal)
with open(E/"essential_outliers.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(outliers[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(outliers)
summary={"pathways":len(cal),"core":sum(1 for x in cal if x["tier"]=="core"),
         "conditional":sum(1 for x in cal if x["tier"]=="conditional"),
         "outlier_rows":len(outliers),"top10":[x["pathway"] for x in cal[:10]],
         "by_tier_gaps":dict(Counter(x["tier"] for x in cal for _ in range(x["gaps"])))}
json.dump(summary,open(E/"calibration_summary.json","w"),indent=2)
print(json.dumps(summary,indent=2))
