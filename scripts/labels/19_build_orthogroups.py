#!/usr/bin/env python3
"""Construye FASTA único de proteomas y clusters DIAMOND para transferencia GO."""
from __future__ import annotations
import argparse, csv, gzip, json, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def fasta_path(data, g):
    cands=[data/"annotations_ab_initio"/g/"proteins.faa",
           data/"annotations"/f"{g}.faa.gz",
           data/"annotations"/f"{g}.faa",
           data/"regulation"/g/"proteins.faa"]
    if g=="subtilis": cands.insert(0, data.parent/"results/subtilis/subtilis.proteins.faa")
    return next((p for p in cands if p.exists()), None)


def iter_fasta(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",errors="replace") as fh:
        name=None; seq=[]
        for line in fh:
            line=line.strip()
            if line.startswith(">"):
                if name: yield name,"".join(seq)
                name=line[1:].split()[0]; seq=[]
            elif line and name is not None:
                seq.append(line)
        if name: yield name,"".join(seq)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/go_indirect")
    ap.add_argument("--threads",type=int,default=8)
    ap.add_argument("--cluster",action="store_true")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    fa=a.out/"proteins.faa"; mapping=a.out/"proteins_map.tsv"
    n=0; genomes=0
    with fa.open("w") as out, mapping.open("w") as mf:
        mf.write("protein_id\tgenome\tgene_id\n")
        for d in sorted((a.data_root/"regulation").iterdir()):
            if not d.is_dir() or not (d/"genes.tsv").exists(): continue
            g=d.name; p=fasta_path(a.data_root,g)
            if not p: continue
            genomes+=1
            for seqid,seq in iter_fasta(p):
                pid=f"{g}|{seqid}"
                out.write(f">{pid}\n{seq}\n"); mf.write(f"{pid}\t{g}\t{seqid}\n"); n+=1
    (a.out/"proteins_summary.json").write_text(json.dumps({"genomes":genomes,"proteins":n,"fasta":str(fa)},indent=1),encoding="utf-8")
    print("proteins",n,"genomes",genomes,flush=True)
    if a.cluster:
        cmd=["diamond","cluster","--db",str(fa),"--out",str(a.out/"clusters.tsv"),
             "--approx-id","30","--threads",str(a.threads),"--memory-limit","8G"]
        print(" ".join(cmd),flush=True)
        subprocess.run(cmd,check=True)
    return 0

if __name__=="__main__": raise SystemExit(main())
