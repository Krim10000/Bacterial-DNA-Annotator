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
    a=ap.parse_args(); s=build(a.data_root,a.out,a.genomes); print(json.dumps(s,indent=1))
if __name__=="__main__": raise SystemExit(main())
