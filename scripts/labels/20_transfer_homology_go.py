#!/usr/bin/env python3
"""Transfiere GO por homología desde proteínas anotadas a no anotadas (ISS)."""
from __future__ import annotations
import argparse, csv, gzip, json, subprocess
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def fasta_iter(path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",errors="replace") as fh:
        name=None; seq=[]
        for line in fh:
            line=line.strip()
            if line.startswith(">"):
                if name: yield name,"".join(seq)
                name=line[1:].split()[0]; seq=[]
            elif line and name is not None: seq.append(line)
        if name: yield name,"".join(seq)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--indir",type=Path,default=ROOT/"data/labels/go_indirect")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/go_indirect")
    ap.add_argument("--threads",type=int,default=8)
    ap.add_argument("--min-pident",type=float,default=30)
    ap.add_argument("--min-qcov",type=float,default=70)
    ap.add_argument("--run",action="store_true")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    # annotated genes: direct MF/CC or BP
    ann=defaultdict(list)  # gene key -> (go_id,source,evidence)
    for p in (a.data_root/"labels/go_mf_cc").glob("*/gene_go.tsv"):
        g=p.parent.name
        for r in read_tsv(p):
            ann[f"{g}|{r['gene_id']}"].append((r["go_id"],r["aspect"],r["source"],r["evidence"]))
    for p in (a.data_root/"labels/go_bp").glob("*/gene_bp.tsv"):
        g=p.parent.name
        for r in read_tsv(p):
            ann[f"{g}|{r['gene_id']}"].append((r["bp_go_id"],"P",r["source"],r["evidence"]))
    mapping={r["protein_id"]:r for r in read_tsv(a.indir/"proteins_map.tsv")}
    donor=a.out/"donors.faa"; query=a.out/"unannotated.faa"
    nd=nq=0
    with donor.open("w") as df, query.open("w") as qf:
        for pid,seq in fasta_iter(a.indir/"proteins.faa"):
            m=mapping.get(pid); 
            if not m: continue
            key=f"{m['genome']}|{m['gene_id']}"
            if key in ann:
                df.write(f">{pid}\n{seq}\n"); nd+=1
            else:
                qf.write(f">{pid}\n{seq}\n"); nq+=1
    print("donors",nd,"queries",nq,flush=True)
    summary={"donors":nd,"queries":nq,"min_pident":a.min_pident,"min_qcov":a.min_qcov}
    if a.run:
        db=a.out/"donors.dmnd"
        subprocess.run(["diamond","makedb","--in",str(donor),"--db",str(db)],check=True)
        hits=a.out/"hits.tsv"
        subprocess.run(["diamond","blastp","--query",str(query),"--db",str(db),
                        "--out",str(hits),"--outfmt","6","qseqid","sseqid","pident",
                        "length","qcovhsp","scovhsp","evalue","bitscore",
                        "--threads",str(a.threads),"--max-target-seqs","1",
                        "--evalue","1e-5","--quiet"],check=True)
        # transfer
        rows=[]; byq=defaultdict(list)
        for line in hits.read_text().splitlines():
            f=line.split("\t")
            if len(f)<8: continue
            q,s,pident,length,qcov,scov,evalue,bits=f[:8]
            if float(pident)<a.min_pident or float(qcov)<a.min_qcov: continue
            byq[q].append((s,float(pident),float(qcov),float(bits)))
        for q,hs in byq.items():
            s,pident,qcov,bits=max(hs,key=lambda x:x[3])
            if s not in ann: continue
            m=mapping.get(q)
            if not m: continue
            for go,asp,src,ev in ann[s]:
                rows.append({"genome":m["genome"],"gene_id":m["gene_id"],"go_id":go,
                             "aspect":asp,"source":"homology_transfer","evidence":"ISS",
                             "with_from":s,"identity":pident,"coverage":qcov,
                             "confidence":"low","donor_source":src,"donor_evidence":ev})
        with (a.out/"gene_go_homology.tsv").open("w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["genome"],delimiter="\t")
            w.writeheader(); w.writerows(rows)
        summary["hits"]=sum(1 for _ in hits.read_text().splitlines())
        summary["transferred_rows"]=len(rows)
        summary["transferred_genes"]=len({(r["genome"],r["gene_id"]) for r in rows})
    (a.out/"homology_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
    return 0
if __name__=="__main__": raise SystemExit(main())
