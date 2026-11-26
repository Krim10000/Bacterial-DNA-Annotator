#!/usr/bin/env python3
"""Construye gene_evidence.tsv para clasificador/sigma/circuitos."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from webapp.gene_evidence import build

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/gene_evidence.tsv")
    ap.add_argument("--genomes",nargs="*",default=None)
    ap.add_argument("--per-genome-dir",type=Path,default=None)
    a=ap.parse_args()
    if a.per_genome_dir:
        import csv
        from webapp.gene_evidence import _protein_map
        a.per_genome_dir.mkdir(parents=True,exist_ok=True)
        genomes=a.genomes
        if not genomes:
            genomes=[d.name for d in sorted((a.data_root/"regulation").iterdir())
                     if d.is_dir() and (d/"genes.tsv").exists()]
        total={"genomes":0,"genes":0,"rows":0,"no_data_rows":0,"by_source":{}}
        for g in genomes:
            out=a.per_genome_dir/f"{g}.tsv"
            if out.exists() and out.stat().st_size>0:
                continue
            s=build(a.data_root,out,[g])
            total["genomes"]+=1; total["genes"]+=s["genes"]; total["rows"]+=s["rows"]; total["no_data_rows"]+=s["no_data_rows"]
            for k,v in s["by_source"].items(): total["by_source"][k]=total["by_source"].get(k,0)+v
            print(g,s["rows"],flush=True)
        (a.per_genome_dir/"summary.json").write_text(json.dumps(total,indent=1),encoding="utf-8")
        print(json.dumps(total,indent=1))
    else:
        s=build(a.data_root,a.out,a.genomes); print(json.dumps(s,indent=1))
if __name__=="__main__": raise SystemExit(main())
