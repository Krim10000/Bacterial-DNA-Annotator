#!/usr/bin/env python3
"""Contexto IGC para no_data usando gene_evidence actual y posición en operón."""
from __future__ import annotations
import argparse, csv, json
from collections import Counter, defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/context")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    has=set(); terms=defaultdict(Counter)
    with (a.data_root/"labels/gene_evidence.tsv").open(encoding="utf-8",errors="replace") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            if r.get("source")=="no_data": continue
            key=(r["genome"],r["gene_id"]); has.add(key)
            if r.get("term_id"):
                terms[key][f"{r['term_id']} {r.get('term_name','')}".strip()]+=1
    has_by=defaultdict(set)
    for k in has: has_by[k[0]].add(k[1])
    rows=[]; total=0; ctx=0; internal=0; edge=0; no_op=0; no_ann=0
    for d in sorted((a.data_root/"regulation").iterdir()):
        if not d.is_dir() or not (d/"genes.tsv").exists(): continue
        g=d.name; genes=list(csv.DictReader((d/"genes.tsv").open(encoding="utf-8",errors="replace"),delimiter="\t"))
        has_g=has_by.get(g,set())
        no=[x for x in genes if (x.get("gene_id") or x.get("locus_tag")) not in has_g]
        total+=len(no)
        opath=d/"operons.tsv"; ops=list(csv.DictReader(opath.open(encoding="utf-8",errors="replace"),delimiter="\t")) if opath.exists() else []
        if not ops: no_op+=len(no); continue
        bygene={}
        for op in ops:
            ids=[x for x in (op.get("gene_ids") or "").split(",") if x]
            for gid in ids: bygene[gid]=op
        for gene in no:
            gid=gene.get("gene_id") or gene.get("locus_tag") or ""
            op=bygene.get(gid)
            if not op: continue
            ids=[x for x in (op.get("gene_ids") or "").split(",") if x]
            ann=[x for x in ids if (g,x) in has]
            if not ann: no_ann+=1; continue
            pos="edge"; 
            try:
                i=ids.index(gid)
                if 0<i<len(ids)-1: pos="internal"; internal+=1
                else: edge+=1
            except ValueError: pass
            c=Counter()
            for x in ann: c.update(terms.get((g,x),{}))
            top=[t for t,_ in c.most_common(5)]
            rows.append({"genome":g,"gene_id":gid,"operon_id":op.get("operon_id",""),
                         "position":pos,"n_members":len(ids),"n_annotated":len(ann),
                         "neighbor_genes":",".join(ann[:20]),
                         "context_terms":";".join(top),"source":"operon_context",
                         "evidence":"IGC","confidence":"low"})
            ctx+=1
    with (a.out/"gene_context.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"total_no_data":total,"with_context":ctx,"internal":internal,"edge":edge,
             "no_operon_file":no_op,"operon_without_annotated":no_ann,"context_rows":len(rows)}
    (a.out/"context_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
