#!/usr/bin/env python3
"""P2Rank sobre structures de sigma, reanudable y con baja prioridad."""
from __future__ import annotations
import argparse, csv, json, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]


def read_tsv(p):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def plddt_from_pdb(path):
    vals=[]
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ","HETATM")) and line[12:16].strip()=="CA":
            try: vals.append(float(line[60:66]))
            except ValueError: pass
    return round(sum(vals)/len(vals),3) if vals else None


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,default=ROOT/"data")
    ap.add_argument("--threads",type=int,default=2)
    ap.add_argument("--limit",type=int,default=0)
    ap.add_argument("--run",action="store_true")
    a=ap.parse_args()
    root=a.data_root
    afdb=root/"metabolism/targets/afdb25k"
    pdbdir=root/"sigma/structures/pdb"
    outroot=root/"sigma/structures/p2rank_sigma"; outroot.mkdir(parents=True,exist_ok=True)
    rows=[]; structures=[]
    for d in sorted((root/"sigma").iterdir()):
        f=d/"sigma_factors.tsv"; gp=root/"regulation"/d.name/"genes.tsv"
        if not f.exists() or not gp.exists(): continue
        g=d.name
        genes={r["gene_id"]:r.get("acc","") for r in read_tsv(gp)}
        for r in read_tsv(f):
            gid=r.get("gene_id",""); acc=genes.get(gid,"")
            if not acc: continue
            pdb=pdbdir/f"{acc}.pdb"
            apdb=afdb/f"{acc}.pdb"
            if pdb.exists(): path,origin=pdb,"PDB_experimental"
            elif apdb.exists(): path,origin=apdb,"AFDB"
            else: continue
            sid=f"{g}|{gid}"
            structures.append((sid,g,gid,acc,path,origin))
    if a.limit: structures=structures[:a.limit]
    print("structures",len(structures),flush=True)
    for sid,g,gid,acc,path,origin in structures:
        out=outroot/sid.replace("|","__")
        pred=out/f"{path.name}_predictions.csv"
        if pred.exists():
            pass
        elif a.run:
            out.mkdir(parents=True,exist_ok=True)
            subprocess.run(["nice","-n","19","p2rank","predict","-f",str(path),
                            "-o",str(out),"-threads",str(a.threads)],
                           check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        if not pred.exists(): continue
        with pred.open(errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.lower().startswith("name"): continue
                f=[x.strip() for x in line.split(",")]
                if len(f)<10: continue
                rows.append({"sigma_id":sid,"genome":g,"gene_id":gid,"accession":acc,
                             "structure_origin":origin,"structure_path":str(path),
                             "pocket_id":f[0],"pocket_rank":f[1],"score":f[2],
                             "probability":f[3],"sas_points":f[4],"surf_atoms":f[5],
                             "residues":f[9] if len(f)>9 else "",
                             "plddt":plddt_from_pdb(path)})
    out=a.out if hasattr(a,"out") else root/"sigma/structures/sigma_pockets.tsv"
    with out.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]) if rows else ["sigma_id"],delimiter="\t"); w.writeheader(); w.writerows(rows)
    summary={"structures":len(structures),"sigma_with_structure":len({x[0] for x in structures}),
             "pocket_rows":len(rows),"sigma_with_pocket":len({r['sigma_id'] for r in rows}),
             "out":str(out),"p2rank_threads":a.threads}
    (root/"sigma/structures/p2rank_sigma_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    print(json.dumps(summary,indent=1))
if __name__=="__main__": main()
