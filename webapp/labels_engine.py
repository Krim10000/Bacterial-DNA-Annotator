"""Motor de etiquetas anidadas gen→operón, sin reanotar.

Insumos: genes/operones existentes, léxico controlado, GOA y KEGG ya descargados.
Salidas: gene_labels.tsv, operon_labels.tsv, label_summary.json.
"""
from __future__ import annotations

import csv, gzip, json, re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Any, Iterable

from webapp.labels import labels_for_text, is_regulator, LABELS

GAF = {
    "ecoli_k12": "ECOLI-uniprot.gaf.gz",
    "subtilis": "BACSU-uniprot.gaf.gz",
}
GO_PATTERNS = {
    "heat_shock": r"response to heat|heat shock|protein folding|chaperon",
    "flagellum": r"flagell|chemotaxis",
    "nitrogen": r"nitrogen|nitrate|nitrite|glutamine",
    "stationary": r"stationary|starvation",
    "envelope": r"cell wall|outer membrane|envelope|peptidoglycan|lipopolysaccharide",
    "sporulation": r"sporulation|spore",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists(): return []
    with path.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_obo(path: Path):
    names, parents = {}, defaultdict(set); cur = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line == "[Term]": cur = None
        elif line.startswith("id: GO:"): cur = line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur] = line.split("name: ",1)[1].strip()
        elif cur and line.startswith("is_a: "): parents[cur].add(line.split("is_a: ",1)[1].split("!")[0].strip())
    return names, parents


def ancestors(go, parents, cache):
    if go in cache: return cache[go]
    out={go}; stack=[go]
    while stack:
        x=stack.pop()
        for p in parents.get(x,()):
            if p not in out: out.add(p); stack.append(p)
    cache[go]=out; return out


def read_gaf(path: Path, allow_iea: bool = False):
    out=defaultdict(set)
    if not path.exists(): return out
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<9: continue
            if not allow_iea and f[6] == "IEA": continue
            out[f[1]].add(f[4])
    return out


def _kegg_links(path: Path):
    out=defaultdict(set)
    if not path.exists(): return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "\t" not in line: continue
        a,b=line.split("\t",1)
        out[a.split(":",1)[-1]].add(b.split(":",1)[-1])
    return out


def _gene_rows(data_root: Path, genome: str):
    return read_tsv(data_root/"regulation"/genome/"genes.tsv")


def _operon_rows(data_root: Path, genome: str):
    p=data_root/"regulation"/genome/"operons.tsv"
    if p.exists(): return read_tsv(p), "gene_ids", ","
    p=data_root.parent/"results"/genome/f"{genome}.operons_auto.tsv"
    if p.exists(): return read_tsv(p), "genes", ";"
    return [], "genes", ";"


def _gaf_for(data_root: Path, genome: str):
    return data_root/"labels"/"external"/"goa"/GAF.get(genome,"")


def _kegg_code(data_root: Path, genome: str):
    if genome == "ecoli_k12":
        return "eco"
    if genome == "subtilis":
        return "bsu"
    for r in read_tsv(data_root/"labels"/"external"/"kegg"/"kegg_org_map.tsv"):
        if r.get("genome")==genome: return r.get("kegg_code","")
    return ""


def build(data_root: Path, genome: str, out_dir: Path, *,
          use_lexicon=True, use_go=True, use_kegg=True,
          allow_iea=False) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    genes=_gene_rows(data_root, genome); ops, op_col, sep=_operon_rows(data_root, genome)
    by_id={}
    for g in genes:
        gid=g.get("gene_id") or g.get("locus_tag") or ""
        if gid: by_id[gid]=g
    gene_labels: dict[str,set[str]]=defaultdict(set)
    provenance=[]
    for g in genes:
        gid=g.get("gene_id") or g.get("locus_tag") or ""
        if not gid: continue
        if use_lexicon and not is_regulator(g.get("gene_name",""), g.get("product","")):
            for lab in labels_for_text(f"{g.get('gene_name','')} {g.get('product','')}"):
                gene_labels[gid].add(lab); provenance.append((gid,lab,"lexicon"))
    # GO
    if use_go:
        names,parents=read_obo(data_root/"labels"/"external"/"go"/"go-basic.obo")
        gaf=read_gaf(_gaf_for(data_root,genome),allow_iea)
        cache={}
        for g in genes:
            gid=g.get("gene_id") or g.get("locus_tag") or ""
            acc=(g.get("acc") or g.get("uniprot") or "").strip()
            for go in gaf.get(acc,()):
                for gidgo in ancestors(go,parents,cache):
                    nm=names.get(gidgo,"").lower()
                    for lab,pat in GO_PATTERNS.items():
                        if nm and re.search(pat,nm):
                            gene_labels[gid].add(f"GO::{lab}")
                            provenance.append((gid,f"GO::{lab}",gidgo))
        # B. subtilis: los genes.tsv no traen acc; usar cds de functional hits.
        if genome == "subtilis":
            fp=data_root.parent/"results"/genome/f"{genome}.functional_annotation.json"
            if fp.exists():
                hits=json.loads(fp.read_text()).get("hits",{})
                for cid,h in hits.items():
                    acc=(h.get("uniprot") or "").strip()
                    for go in gaf.get(acc,()):
                        for gidgo in ancestors(go,parents,cache):
                            nm=names.get(gidgo,"").lower()
                            for lab,pat in GO_PATTERNS.items():
                                if nm and re.search(pat,nm):
                                    gene_labels[cid].add(f"GO::{lab}")
                                    provenance.append((cid,f"GO::{lab}",gidgo))
    # KEGG
    if use_kegg:
        code=_kegg_code(data_root,genome)
        if code:
            paths=_kegg_links(data_root/"labels"/"external"/"kegg"/code/"pathway.tsv")
            for g in genes:
                gid=g.get("gene_id") or g.get("locus_tag") or ""
                locus=(g.get("locus_tag") or "").replace("_","")
                for path in set(paths.get(locus,())) | set(paths.get(g.get("locus_tag",""), ())):
                    gene_labels[gid].add(f"KEGG::{path}")
                    provenance.append((gid,f"KEGG::{path}",code))
    # operons
    op_rows=[]
    for op in ops:
        oid=op.get("operon_id","")
        ids=[x for x in (op.get(op_col) or "").split(sep) if x]
        labs=sorted({l for gid in ids for l in gene_labels.get(gid,())})
        op_rows.append({"operon_id":oid,"label":op.get("label",""),
                        "n_genes":len(ids),"labels":";".join(labs)})
    with (out_dir/"gene_labels.tsv").open("w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t"); w.writerow(["genome","gene_id","label","source"])
        for gid,lab,src in provenance: w.writerow([genome,gid,lab,src])
    with (out_dir/"operon_labels.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["operon_id","label","n_genes","labels"],delimiter="\t")
        w.writeheader(); w.writerows(op_rows)
    summary={"genome":genome,"n_genes":len(genes),"n_operons":len(ops),
             "n_gene_labels":sum(1 for x in gene_labels.values() if x),
             "n_operons_with_label":sum(1 for x in op_rows if x["labels"]),
             "label_counts":dict(Counter(l for x in gene_labels.values() for l in x)),
             "sources":{"lexicon":use_lexicon,"go":use_go,"kegg":use_kegg,
                        "allow_iea":allow_iea,"kegg_code":_kegg_code(data_root,genome)}}
    (out_dir/"label_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    return {"out_dir":str(out_dir),"summary":summary,
            "gene_labels":str(out_dir/"gene_labels.tsv"),
            "operon_labels":str(out_dir/"operon_labels.tsv")}
