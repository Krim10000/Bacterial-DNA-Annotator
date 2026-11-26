#!/usr/bin/env python3
"""Indexa GO BP directo (aspect P) por accesión UniProt desde los GAF."""
from __future__ import annotations
import csv, gzip, json, sqlite3, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def acc_norm(x): return (x or "").strip().split("-")[0]
def ev_class(e):
    if e=="IEA": return "electronic"
    if e in {"EXP","IDA","IPI","IMP","IGI","IEP"}: return "experimental"
    return "curated"


def main():
    data=ROOT/"data"; out=data/"labels/go_bp_direct"; out.mkdir(parents=True,exist_ok=True)
    accs=set()
    for d in (data/"regulation").iterdir():
        if not (d/"genes.tsv").exists(): continue
        for r in read_tsv(d/"genes.tsv"):
            a=acc_norm(r.get("acc",""))
            if a: accs.add(a)
    # B. subtilis functional
    fp=ROOT/"results/subtilis/subtilis.functional_annotation.json"
    if fp.exists():
        for h in json.loads(fp.read_text()).get("hits",{}).values():
            a=acc_norm(h.get("uniprot",""))
            if a: accs.add(a)
    names={}; ns={}; cur=None
    for line in (data/"labels/external/go/go.obo").read_text(encoding="utf-8",errors="replace").splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
        elif cur and line.startswith("namespace: "): ns[cur]=line.split("namespace: ",1)[1].strip()
    db=out/"goa_bp.sqlite"; db.unlink(missing_ok=True)
    con=sqlite3.connect(db); con.execute("PRAGMA journal_mode=DELETE"); con.execute("PRAGMA synchronous=OFF")
    con.execute("""CREATE TABLE bp(acc TEXT, go_id TEXT, go_name TEXT, aspect TEXT,
                   evidence TEXT, evidence_class TEXT, source TEXT, taxon TEXT,
                   assigned_by TEXT, qualifier TEXT, with_from TEXT,
                   PRIMARY KEY(acc,go_id,evidence,source))""")
    files=0; n=0
    for p in sorted((data/"labels/external/goa/by_code").glob("*.gaf.gz")):
        files+=1; batch=[]
        with gzip.open(p,"rt",errors="replace") as fh:
            for line in fh:
                if line.startswith("!"): continue
                f=line.rstrip("\n").split("\t")
                if len(f)<15 or f[8]!="P": continue
                acc=acc_norm(f[1])
                if acc not in accs: continue
                tax=next((x.split(":",1)[1] for x in f if x.startswith("taxon:")),"")
                batch.append((acc,f[4],names.get(f[4],""),f[8],f[6],ev_class(f[6]),p.name,tax,f[14],f[3],f[7]))
                if len(batch)>=5000:
                    con.executemany("INSERT OR IGNORE INTO bp VALUES (?,?,?,?,?,?,?,?,?,?,?)",batch); n+=len(batch); batch=[]
        if batch:
            con.executemany("INSERT OR IGNORE INTO bp VALUES (?,?,?,?,?,?,?,?,?,?,?)",batch); n+=len(batch)
        con.commit()
        if files%20==0: print("files",files,"rows",n,flush=True)
    con.execute("CREATE INDEX idx_acc ON bp(acc)"); con.commit()
    total=con.execute("select count(*) from bp").fetchone()[0]
    summary={"files":files,"accs":len(accs),"rows":total,"db":str(db)}
    (out/"bp_index_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
    return 0

if __name__=="__main__": raise SystemExit(main())
