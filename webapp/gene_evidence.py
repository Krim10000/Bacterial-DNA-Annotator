"""Contrato de evidencia gen→término para clasificador/sigma/circuitos."""
from __future__ import annotations
import csv, json
from collections import defaultdict
from functools import lru_cache
from pathlib import Path


def read_tsv(p: Path):
    if not p.exists(): return []
    with p.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


@lru_cache(maxsize=4)
def _protein_map(data_root: Path):
    out={}
    for r in read_tsv(data_root/"labels/go_indirect/proteins_map.tsv"):
        out[(r.get("genome",""),r.get("gene_id",""))]=r.get("protein_id","")
    return out


def build(data_root: Path, out_path: Path, genomes: list[str] | None = None) -> dict:
    pmap=_protein_map(data_root)
    if genomes is None:
        genomes=[]
        for d in sorted((data_root/"regulation").iterdir()):
            if d.is_dir() and (d/"genes.tsv").exists(): genomes.append(d.name)
    out_path.parent.mkdir(parents=True,exist_ok=True)
    n=0; n_genes=0; n_no=0; by_source=defaultdict(int)
    fields=["genome","gene_id","protein_id","source","evidence","aspect","term_id","term_name","confidence","with_from","coverage","notes"]
    with out_path.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t"); w.writeheader()
        for g in genomes:
            genes=read_tsv(data_root/"regulation"/g/"genes.tsv")
            seen=set(); rows=[]
            for r in read_tsv(data_root/"labels/go_mf_cc"/g/"gene_go.tsv"):
                gid=r.get("gene_id",""); term=r.get("go_id","")
                if not gid or not term: continue
                key=(gid,r.get("source",""),term,r.get("evidence",""))
                if key in seen: continue
                seen.add(key)
                rows.append({"genome":g,"gene_id":gid,"protein_id":pmap.get((g,gid),gid),
                             "source":r.get("source",""),"evidence":r.get("evidence",""),
                             "aspect":r.get("aspect",""),"term_id":term,"term_name":r.get("go_name",""),
                             "confidence":"high" if r.get("evidence") in {"EXP","IDA","IPI","IMP","IGI","IEP"} else "medium",
                             "with_from":"","coverage":"direct","notes":"GOA"})
            for r in read_tsv(data_root/"labels/go_bp"/g/"gene_bp.tsv"):
                gid=r.get("gene_id",""); term=r.get("bp_go_id","")
                if not gid or not term: continue
                key=(gid,r.get("source",""),term,r.get("evidence",""))
                if key in seen: continue
                seen.add(key)
                rows.append({"genome":g,"gene_id":gid,"protein_id":pmap.get((g,gid),gid),
                             "source":r.get("source",""),"evidence":r.get("evidence",""),
                             "aspect":"P","term_id":term,"term_name":r.get("bp_go_name",""),
                             "confidence":r.get("confidence","medium"),"with_from":r.get("with_from",""),
                             "coverage":"direct" if r.get("origin")=="direct_gaf" else "inferred",
                             "notes":r.get("origin","")})
            if not rows:
                for gene in genes:
                    gid=gene.get("gene_id") or gene.get("locus_tag") or ""
                    if not gid: continue
                    rows.append({"genome":g,"gene_id":gid,"protein_id":pmap.get((g,gid),gid),
                                 "source":"no_data","evidence":"no_data","aspect":"","term_id":"",
                                 "term_name":"","confidence":"unknown","with_from":"","coverage":"none",
                                 "notes":"sin GO directo"})
                n_no+=len(rows)
            for r in rows: w.writerow(r); n+=1; by_source[r["source"]]+=1
            n_genes+=len(genes)
    summary={"genomes":len(genomes),"genes":n_genes,"rows":n,"no_data_rows":n_no,
             "by_source":dict(sorted(by_source.items(),key=lambda kv:-kv[1])),
             "schema":"gene_evidence.tsv","note":"no_data se emite por gen sin ninguna evidencia"}
    (out_path.parent/"gene_evidence_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    return summary
