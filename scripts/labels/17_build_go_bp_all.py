#!/usr/bin/env python3
"""BP oficial MF→BP para todos los genomas con gene_go.tsv."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from webapp.go_bp import build_for_genome

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/go_bp")
    a=ap.parse_args(); summary={}
    for d in sorted((a.data_root/"labels"/"go_mf_cc").iterdir()):
        if not d.is_dir() or not (d/"gene_go.tsv").exists(): continue
        r=build_for_genome(a.data_root,d.name,a.out/d.name); summary[d.name]=r["summary"]
    (a.out/"bp_index_summary.json").write_text(json.dumps({"genomes":len(summary),"per_genome":summary},indent=1),encoding="utf-8")
    print("genomes",len(summary),"rows",sum(v["n_bp_rows"] for v in summary.values()),
          "genes",sum(v["n_genes_bp"] for v in summary.values()))
    print({k:summary[k] for k in ("ecoli_k12","subtilis") if k in summary})
    return 0
if __name__=="__main__": raise SystemExit(main())
