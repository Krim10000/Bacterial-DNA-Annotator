#!/usr/bin/env python3
"""Co-ocurrencia de candidatos en los 400 genomas."""
import csv, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
D=Path(__file__).resolve().parents[2]; O=D/"data"/"metabolism"/"targets"; C=O/"comparative"; C.mkdir(exist_ok=True)
db=O/"candidates.dmnd"
if not db.exists():
    subprocess.run(["diamond","makedb","--in",str(O/"candidate_proteins.faa"),"--db",str(db)],check=True)
genomes=sorted(p.name for p in (D/"data"/"annotations_ab_initio").iterdir() if p.is_dir())
def run(g):
    faa=D/"data"/"annotations_ab_initio"/g/"proteins.faa"; out=C/f"{g}.tsv"
    if not faa.exists() or out.exists(): return g
    subprocess.run(["diamond","blastp","--query",str(faa),"--db",str(db),
                    "--outfmt","6","qseqid","sseqid","pident","qcovhsp","evalue",
                    "--max-target-seqs","20","--fast","--threads","2","--out",str(out)],
                   check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    return g
with ThreadPoolExecutor(max_workers=8) as ex:
    for i,g in enumerate(ex.map(run,genomes),1):
        if i%25==0: print(i,"/",len(genomes),flush=True)
print("done",len(genomes))
