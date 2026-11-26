#!/usr/bin/env python3
"""Evalúa no_data por contexto de operón/vecindad (IGC, baja confianza)."""
from __future__ import annotations
import argparse, csv, json
from collections import Counter, defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/context")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    rows=[]; summary={}; total_no=0; with_ctx=0; singleton=0; no_op=0
    for gdir in sorted((a.data_root/"regulation").iterdir()):
        if not gdir.is_dir() or not (gdir/"genes.tsv").exists(): continue
        g=gdir.name
        ev=read(a.data_root/"labels/gene_evidence"/f"{g}.tsv")
        has=set(); terms=defaultdict(Counter)
        for r in ev:
            if r.get("source")=="no_data": continue
            gid=r.get("gene_id",""); has.add(gid)
            if r.get("term_id"): terms[gid][f"{r.get('term_id')} {r.get('term_name','')}".strip()]+=1
        genes=read(gdir/"genes.tsv")
        no=[r for r in genes if (r.get("gene_id") or r.get("locus_tag")) not in has]
        total_no+=len(no)
        # operons
        opath=gdir/"operons.tsv"
        ops=read(opath) if opath.exists() else []
        if not ops:
            no_op+=len(no); continue
        by_gene={}
        for op in ops:
            ids=[x for x in (op.get("gene_ids") or "").split(",") if x]
            for gid in ids: by_gene[gid]=op
        for gene in no:
            gid=gene.get("gene_id") or gene.get("locus_tag") or ""
            op=by_gene.get(gid)
            if not op:
                singleton+=1; continue
            ids=[x for x in (op.get("gene_ids") or "").split(",") if x]
            ann=[x for x in ids if x in has]
            if not ann:
                continue
            c=Counter()
            for x in ann: c.update(terms.get(x,{}))
            top=[t for t,_ in c.most_common(5)]
            rows.append({"genome":g,"gene_id":gid,"operon_id":op.get("operon_id",""),
                         "context_type":"operon_with_annotated_neighbor",
                         "n_members":len(ids),"n_annotated":len(ann),
                         "neighbor_genes":",".join(ann[:20]),
                         "context_terms":";".join(top),
                         "source":"operon_context","evidence":"IGC","confidence":"low"})
            with_ctx+=1
    out=a.out/"gene_context.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"total_no_data":total_no,"with_context":with_ctx,"singleton_or_no_operon":singleton,
             "no_operon_file":no_op,"context_rows":len(rows)}
    (a.out/"context_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
