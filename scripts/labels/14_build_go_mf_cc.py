#!/usr/bin/env python3
"""Indexa GO MF/CC por accesión UniProt en SQLite y anota todos los genomas."""
from __future__ import annotations
import argparse, csv, gzip, json, sqlite3, sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))


def read_obo(path: Path):
    names={}; ns={}; cur=None
    for line in path.read_text(encoding="utf-8",errors="replace").splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
        elif cur and line.startswith("namespace: "): ns[cur]=line.split("namespace: ",1)[1].strip()
    return names,ns


def acc_norm(a): return (a or "").strip().split("-")[0]
def class_ev(e):
    if e=="IEA": return "electronic"
    if e in {"EXP","IDA","IPI","IMP","IGI","IEP"}: return "experimental"
    return "curated"


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--out",type=Path,default=ROOT/"data/labels/go_mf_cc")
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    genomes={}
    for d in sorted((a.data_root/"regulation").iterdir()):
        if d.is_dir() and (d/"genes.tsv").exists():
            genomes[d.name]=list(csv.DictReader((d/"genes.tsv").open(encoding="utf-8",errors="replace"),delimiter="\t"))
    accs={acc_norm(r.get("acc","")) for rows in genomes.values() for r in rows if acc_norm(r.get("acc",""))}
    print("genomes",len(genomes),"genes",sum(len(v) for v in genomes.values()),"accs",len(accs),flush=True)
    names,ns=read_obo(a.data_root/"labels/external/go/go-basic.obo")
    db=a.out/"goa.sqlite"; db.unlink(missing_ok=True)
    con=sqlite3.connect(db); con.execute("PRAGMA journal_mode=DELETE"); con.execute("PRAGMA synchronous=OFF")
    con.execute("""CREATE TABLE go(acc TEXT, go_id TEXT, go_name TEXT, aspect TEXT,
                   namespace TEXT, evidence TEXT, evidence_class TEXT, source TEXT,
                   taxon TEXT, assigned_by TEXT, PRIMARY KEY(acc,go_id,evidence,source))""")
    files=0; total=0
    for p in sorted((a.data_root/"labels/external/goa/by_code").glob("*.gaf.gz")):
        files+=1; batch=[]
        with gzip.open(p,"rt",errors="replace") as fh:
            for line in fh:
                if line.startswith("!"): continue
                f=line.rstrip("\n").split("\t")
                if len(f)<15: continue
                acc=acc_norm(f[1])
                if acc not in accs: continue
                if f[8] not in ("F","C"): continue
                go=f[4]; ev=f[6]
                tax=next((x.split(":",1)[1] for x in f if x.startswith("taxon:")),"")
                batch.append((acc,go,names.get(go,""),f[8],ns.get(go,""),ev,class_ev(ev),p.name,tax,f[14]))
                if len(batch)>=5000:
                    con.executemany("INSERT OR IGNORE INTO go VALUES (?,?,?,?,?,?,?,?,?,?)",batch); total+=len(batch); batch=[]
        if batch:
            con.executemany("INSERT OR IGNORE INTO go VALUES (?,?,?,?,?,?,?,?,?,?)",batch); total+=len(batch)
        con.commit()
        if files%20==0: print("files",files,"rows",total,flush=True)
    con.execute("CREATE INDEX idx_acc ON go(acc)"); con.commit()
    print("indexed rows",total,"files",files,flush=True)
    by_acc=defaultdict(list)
    for r in con.execute("SELECT acc,go_id,go_name,aspect,namespace,evidence,evidence_class,source,taxon,assigned_by FROM go"):
        by_acc[r[0]].append(r)
    summary={}
    for g,rows in genomes.items():
        out=a.out/g; out.mkdir(parents=True,exist_ok=True); n=0; annotated=set(); terms=set()
        with (out/"gene_go.tsv").open("w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=["genome","gene_id","gene_name","acc","go_id","go_name","aspect","namespace","evidence","evidence_class","source","taxon","assigned_by"],delimiter="\t")
            w.writeheader()
            for r in rows:
                acc=acc_norm(r.get("acc",""))
                if not acc: continue
                for x in by_acc.get(acc,[]):
                    w.writerow({"genome":g,"gene_id":r.get("gene_id",""),"gene_name":r.get("gene_name",""),
                                "acc":acc,"go_id":x[1],"go_name":x[2],"aspect":x[3],"namespace":x[4],
                                "evidence":x[5],"evidence_class":x[6],"source":x[7],"taxon":x[8],"assigned_by":x[9]})
                    n+=1; annotated.add(r.get("gene_id")); terms.add(x[1])
        summary[g]={"n_genes":len(rows),"n_with_acc":sum(1 for r in rows if acc_norm(r.get("acc",""))),
                    "n_go_rows":n,"n_genes_annotated":len(annotated),"n_go_terms":len(terms),"n_gaf_files":files}
        (out/"go_summary.json").write_text(json.dumps(summary[g],indent=1),encoding="utf-8")
    (a.out/"go_index_summary.json").write_text(json.dumps({"genomes":len(genomes),"accs":len(accs),"indexed_rows":total,"gaf_files":files,"per_genome":summary},indent=1),encoding="utf-8")
    print(json.dumps({k:summary[k] for k in ("ecoli_k12","subtilis") if k in summary},indent=1))
    print("annotated genes",sum(v["n_genes_annotated"] for v in summary.values()),"of",sum(v["n_genes"] for v in summary.values()))
    return 0

if __name__=="__main__": raise SystemExit(main())
