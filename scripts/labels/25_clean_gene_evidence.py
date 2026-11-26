#!/usr/bin/env python3
"""Deduplica no_data y elimina no_data de genes con evidencia."""
from __future__ import annotations
import csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def main():
    data=ROOT/"data"; p=data/"labels/gene_evidence.tsv"; tmp=p.with_suffix(".clean.tmp")
    evidence=set(); rows=0
    with p.open(encoding="utf-8",errors="replace") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            if r.get("source")!="no_data": evidence.add((r["genome"],r["gene_id"]))
    with p.open(encoding="utf-8",errors="replace") as fh, tmp.open("w",newline="") as of:
        rd=csv.DictReader(fh,delimiter="\t"); fields=rd.fieldnames
        w=csv.DictWriter(of,fieldnames=fields,delimiter="\t"); w.writeheader()
        for r in rd:
            if r.get("source")=="no_data": continue
            w.writerow(r); rows+=1
        for d in sorted((data/"regulation").iterdir()):
            if not d.is_dir() or not (d/"genes.tsv").exists(): continue
            g=d.name
            with (d/"genes.tsv").open(encoding="utf-8",errors="replace") as gf:
                for gene in csv.DictReader(gf,delimiter="\t"):
                    gid=gene.get("gene_id") or gene.get("locus_tag") or ""
                    if gid and (g,gid) not in evidence:
                        w.writerow({"genome":g,"gene_id":gid,"protein_id":gid,
                                    "source":"no_data","evidence":"no_data","aspect":"",
                                    "term_id":"","term_name":"","confidence":"unknown",
                                    "with_from":"","coverage":"none","notes":"sin evidencia GO"})
                        rows+=1
    tmp.replace(p)
    summary={"rows":rows,"genes_with_evidence":len(evidence),"no_data_rows":rows-len(evidence)}
    (data/"labels/gene_evidence_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
