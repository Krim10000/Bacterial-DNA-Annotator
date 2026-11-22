#!/usr/bin/env python3
"""Busca homólogos de los targets esenciales en la bolsa de proteínas."""
import csv, subprocess, sys, os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=Path(os.environ.get("TARGET_DIR", D/"data"/"metabolism"/"targets")); S=O/"search"; S.mkdir(exist_ok=True)
genomes=sorted({r["genome"] for r in csv.DictReader(open(O/"target_genomes.tsv"),delimiter="\t")})
db=O/"refs.dmnd"
if not db.exists():
    subprocess.run(["diamond","makedb","--in",str(O/"refs.faa"),"--db",str(db)],check=True)
def run(g):
    faa=D/"data"/"annotations_ab_initio"/g/"proteins.faa"
    out=S/f"{g}.diamond.tsv"; dom=S/f"{g}.domtbl"
    if not out.exists():
        subprocess.run(["diamond","blastp","--query",str(faa),"--db",str(db),"--outfmt","6","qseqid","sseqid","pident","qcovhsp","evalue","bitscore",
                        "--max-target-seqs","20","--fast","--threads","2","--out",str(out)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    if not dom.exists():
        subprocess.run(["hmmsearch","--domtblout",str(dom),"--noali","--cut_ga",str(O/"panel.hmm"),str(faa)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    return g
with ThreadPoolExecutor(max_workers=8) as ex:
    for i,g in enumerate(ex.map(run,genomes),1):
        if i%10==0: print(i,"/",len(genomes),flush=True)
print("done",len(genomes))
