"""BP por propagación oficial MF→BP (GO_REF:0000108) desde go.obo."""
from __future__ import annotations
import csv, json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path


def read_tsv(path: Path):
    if not path.exists(): return []
    with path.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


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
    go_rows=read_tsv(data_root/"labels/go_mf_cc"/genome/"gene_go.tsv")
    rows=[]; seen=set()
    for r in go_rows:
        if r.get("aspect")!="F": continue
        for bp in mf_to_bp.get(r["go_id"],[]):
            key=(r["gene_id"],bp)
            if key in seen: continue
            seen.add(key)
            rows.append({"genome":genome,"gene_id":r["gene_id"],"gene_name":r.get("gene_name",""),
                         "acc":r.get("acc",""),"bp_go_id":bp,"bp_go_name":names.get(bp,""),
                         "mf_go_id":r["go_id"],"mf_go_name":r.get("go_name",""),
                         "evidence":"IEA","source":"GO_REF:0000108","with_from":r["go_id"],
                         "confidence":"automatic"})
    with (out_dir/"gene_bp.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["genome","gene_id","gene_name","acc","bp_go_id","bp_go_name","mf_go_id","mf_go_name","evidence","source","with_from","confidence"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"genome":genome,"n_genes":len({r["gene_id"] for r in go_rows}),
             "n_bp_rows":len(rows),"n_genes_bp":len({r["gene_id"] for r in rows}),
             "n_bp_terms":len({r["bp_go_id"] for r in rows}),
             "n_mf_links":len(mf_to_bp),"source":"go.obo MF part_of BP; GO_REF:0000108"}
    (out_dir/"bp_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    return {"out_dir":str(out_dir),"summary":summary}
