#!/usr/bin/env python3
"""Contexto conservado por cluster para genes no_data (IGC baja)."""
from __future__ import annotations
import argparse, csv, json
from collections import Counter, defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/context")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    # protein_id -> (genome,gene_id)
    pmap={r["protein_id"]:(r["genome"],r["gene_id"]) for r in read_tsv(a.data_root/"labels/go_indirect/proteins_map.tsv")}
    # gene terms (only annotated)
    terms=defaultdict(Counter); nodata=set(); key_of={}
    with (a.data_root/"labels/gene_evidence.tsv").open(encoding="utf-8",errors="replace") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            key=(r["genome"],r["gene_id"]); key_of[key]=r.get("protein_id","")
            if r.get("source")=="no_data": nodata.add(key)
            elif r.get("term_id"): terms[key][f"{r['term_id']} {r.get('term_name','')}".strip()]+=1
    # cluster membership
    cluster=defaultdict(list)
    with (a.data_root/"labels/go_indirect/clusters.tsv").open(errors="replace") as fh:
        for line in fh:
            f=line.rstrip().split("\t")
            if len(f)>=2: cluster[f[0]].append(f[1])
    # build pid->cluster efficiently
    pid2rep={}
    for rep,members in cluster.items():
        for m in members: pid2rep[m]=rep
    donors=defaultdict(list)
    for key,t in terms.items():
        pid=f"{key[0]}|{key[1]}"; rep=pid2rep.get(pid)
        if rep: donors[rep].append((key,t))
    rows=[]; with_ctx=0; no_cluster=0
    for key in sorted(nodata):
        pid=f"{key[0]}|{key[1]}"; rep=pid2rep.get(pid)
        if not rep: no_cluster+=1; continue
        ds=donors.get(rep)
        if not ds: continue
        c=Counter()
        for _,t in ds: c.update(t)
        top=[x for x,_ in c.most_common(5)]
        rows.append({"genome":key[0],"gene_id":key[1],"protein_id":pid,"cluster":rep,
                     "n_donors":len(ds),"context_terms":";".join(top),
                     "source":"cluster_context","evidence":"IGC","confidence":"low"})
        with_ctx+=1
    out=a.out/"gene_cluster_context.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"no_data":len(nodata),"with_cluster":with_ctx,"no_cluster":no_cluster,
             "clusters":len(cluster),"annotated_genes":len(terms),"rows":len(rows)}
    (a.out/"cluster_context_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
