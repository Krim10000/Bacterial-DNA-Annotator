"""BP: directo desde GAF (aspect P) + propagación oficial MF→BP."""
from __future__ import annotations
import csv, json, sqlite3
from collections import defaultdict
from functools import lru_cache
from pathlib import Path


def read_tsv(path: Path):
    if not path.exists(): return []
    with path.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


@lru_cache(maxsize=4)
def load_direct_bp(data_root: Path):
    idx=defaultdict(list)
    db=data_root/"labels"/"go_bp_direct"/"goa_bp.sqlite"
    if not db.exists(): return idx
    con=sqlite3.connect(db)
    for r in con.execute("SELECT acc,go_id,go_name,evidence,evidence_class,source,taxon,assigned_by,qualifier,with_from FROM bp"):
        idx[r[0]].append(r)
    return idx


@lru_cache(maxsize=4)
def load_mf_to_bp(data_root: Path):
    names={}; ns={}; obsolete=set(); partof=defaultdict(set); cur=None
    path=data_root/"labels/external/go/go.obo"
    for line in path.read_text(encoding="utf-8",errors="replace").splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
        elif cur and line.startswith("namespace: "): ns[cur]=line.split("namespace: ",1)[1].strip()
        elif cur and line.startswith("is_obsolete: true"): obsolete.add(cur)
        elif cur and line.startswith("relationship: part_of GO:"):
            partof[cur].add(line.split("part_of ",1)[1].split()[0])
    mf_to_bp={}
    for mf,parents in partof.items():
        if ns.get(mf)!="molecular_function" or mf in obsolete: continue
        bps=[p for p in parents if ns.get(p)=="biological_process" and p not in obsolete]
        if bps: mf_to_bp[mf]=bps
    return names,ns,mf_to_bp


def build_for_genome(data_root: Path, genome: str, out_dir: Path):
    out_dir.mkdir(parents=True,exist_ok=True)
    names,ns,mf_to_bp=load_mf_to_bp(data_root)
    direct=load_direct_bp(data_root)
    genes=read_tsv(data_root/"regulation"/genome/"genes.tsv")
    bsu={}
    if genome=="subtilis":
        from webapp.go_mf_cc import _bsu_acc_by_gene
        bsu=_bsu_acc_by_gene(data_root)
    go_rows=read_tsv(data_root/"labels"/"go_mf_cc"/genome/"gene_go.tsv")
    rows=[]; seen=set(); n_direct=0; n_infer=0
    # directo: GAF aspect P
    for g in genes:
        acc=(g.get("acc") or bsu.get(g.get("gene_id"),"") or "").strip().split("-")[0]
        for r in direct.get(acc,[]):
            bp=r[1]; key=(g.get("gene_id"),bp)
            if key in seen: continue
            seen.add(key)
            rows.append({"genome":genome,"gene_id":g.get("gene_id",""),"gene_name":g.get("gene_name",""),
                         "acc":acc,"bp_go_id":bp,"bp_go_name":r[2],"mf_go_id":"","mf_go_name":"",
                         "origin":"direct_gaf","evidence":r[3],"source":r[5],
                         "with_from":r[9] or "","confidence":("high" if r[4]=="experimental" else "medium")})
            n_direct+=1
    # inferido: MF --part_of--> BP
    for r in go_rows:
        if r.get("aspect")!="F": continue
        for bp in mf_to_bp.get(r["go_id"],[]):
            key=(r["gene_id"],bp)
            if key in seen: continue
            seen.add(key)
            rows.append({"genome":genome,"gene_id":r["gene_id"],"gene_name":r.get("gene_name",""),
                         "acc":r.get("acc",""),"bp_go_id":bp,"bp_go_name":names.get(bp,""),
                         "mf_go_id":r["go_id"],"mf_go_name":r.get("go_name",""),
                         "origin":"mf_part_of","evidence":"IEA","source":"GO_REF:0000108",
                         "with_from":r["go_id"],"confidence":"automatic"})
            n_infer+=1
    with (out_dir/"gene_bp.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["genome","gene_id","gene_name","acc","bp_go_id","bp_go_name","mf_go_id","mf_go_name","origin","evidence","source","with_from","confidence"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"genome":genome,"n_genes":len(genes),"n_bp_rows":len(rows),
             "n_genes_bp":len({r["gene_id"] for r in rows}),
             "n_bp_terms":len({r["bp_go_id"] for r in rows}),
             "n_direct":n_direct,"n_mf_part_of":n_infer,
             "n_mf_links":len(mf_to_bp),
             "sources":{"direct_gaf":"GAF aspect P","mf_part_of":"GO_REF:0000108"}}
    (out_dir/"bp_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    return {"out_dir":str(out_dir),"summary":summary}
