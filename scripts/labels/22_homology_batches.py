#!/usr/bin/env python3
"""Transferencia GO por homología en lotes (ISS)."""
from __future__ import annotations
import argparse, csv, gzip, json, subprocess
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def fasta_records(path):
    with path.open(encoding="utf-8",errors="replace") as fh:
        name=None; seq=[]
        for line in fh:
            line=line.strip()
            if line.startswith(">"):
                if name: yield name,"".join(seq)
                name=line[1:].split()[0]; seq=[]
            elif line and name is not None: seq.append(line)
        if name: yield name,"".join(seq)


def split_batches(src, outdir, size):
    outdir.mkdir(parents=True,exist_ok=True)
    paths=[]; fh=None; n=0; idx=0
    for name,seq in fasta_records(src):
        if fh is None or n>=size:
            if fh: fh.close()
            idx+=1; p=outdir/f"query_{idx:04d}.faa"; fh=p.open("w"); paths.append(p); n=0
        fh.write(f">{name}\n{seq}\n"); n+=1
    if fh: fh.close()
    return paths


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir",type=Path,default=ROOT/"data/labels/go_indirect")
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--batch-size",type=int,default=20000)
    ap.add_argument("--max-batches",type=int,default=0)
    ap.add_argument("--threads",type=int,default=4)
    ap.add_argument("--run",action="store_true")
    a=ap.parse_args()
    queries=split_batches(a.indir/"unannotated.faa",a.indir/"batches",a.batch_size)
    if a.max_batches: queries=queries[:a.max_batches]
    print("batches",len(queries),flush=True)
    if not a.run:
        print("preparado; usa --run",flush=True); return 0
    ann=defaultdict(list)
    for p in (a.data_root/"labels/go_mf_cc").glob("*/gene_go.tsv"):
        g=p.parent.name
        for r in read_tsv(p):
            ann[f"{g}|{r['gene_id']}"].append((r["go_id"],r["aspect"],r["source"],r["evidence"]))
    for p in (a.data_root/"labels/go_bp").glob("*/gene_bp.tsv"):
        g=p.parent.name
        for r in read_tsv(p):
            ann[f"{g}|{r['gene_id']}"].append((r["bp_go_id"],"P",r["origin"],r["evidence"]))
    mapping={r["protein_id"]:r for r in read_tsv(a.indir/"proteins_map.tsv")}
    rows=[]; per=[]
    for p in queries:
        hits=p.with_suffix(".hits.tsv")
        cmd=["diamond","blastp","--query",str(p),"--db",str(a.indir/"donors.dmnd"),
             "--out",str(hits),"--outfmt","6","qseqid","sseqid","pident","length",
             "qcovhsp","scovhsp","evalue","bitscore","--threads",str(a.threads),
             "--max-target-seqs","1","--evalue","1e-5","--quiet"]
        subprocess.run(cmd,check=True)
        n=0
        for line in hits.read_text().splitlines():
            f=line.split("\t")
            if len(f)<8: continue
            q,s,pid,ln,qcov,scov,ev,bits=f[:8]
            if float(pid)<30 or float(qcov)<70: continue
            m=mapping.get(q)
            if not m: continue
            for go,asp,src,evd in ann.get(s,[]):
                rows.append({"genome":m["genome"],"gene_id":m["gene_id"],"go_id":go,"aspect":asp,
                             "source":"homology_transfer","evidence":"ISS","with_from":s,
                             "identity":pid,"coverage":qcov,"confidence":"low",
                             "donor_source":src,"donor_evidence":evd}); n+=1
        per.append({"batch":p.name,"transferred":n})
        print(p.name,n,flush=True)
    out=a.indir/"gene_go_homology_batches.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"batches":len(queries),"rows":len(rows),"genes":len({(r['genome'],r['gene_id']) for r in rows}),"per_batch":per}
    (a.indir/"homology_batches_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
    return 0

if __name__=="__main__": raise SystemExit(main())
