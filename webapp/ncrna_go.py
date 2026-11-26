"""GO para rRNA/tRNA/ncRNA desde Rfam2GO, independiente de proteínas."""
from __future__ import annotations
import csv, gzip, json, re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

# Rfam genéricos cuando el GFF no trae Dbxref.
TYPE_RFAM={"rRNA":"RF00001","tRNA":"RF00005","tmRNA":"RF00023"}
RRNA_RFAM={"16S":"RF00177","23S":"RF02541","5S":"RF00001"}


@lru_cache(maxsize=4)
def _go_meta(data_root: Path):
    names={}; ns={}; cur=None
    for line in (data_root/"labels/external/go/go.obo").read_text(encoding="utf-8",errors="replace").splitlines():
        if line=="[Term]": cur=None
        elif line.startswith("id: GO:"): cur=line.split("id: ",1)[1].strip()
        elif cur and line.startswith("name: "): names[cur]=line.split("name: ",1)[1].strip()
        elif cur and line.startswith("namespace: "): ns[cur]=line.split("namespace: ",1)[1].strip()
    return names,ns


@lru_cache(maxsize=4)
def load_rfam2go(data_root: Path):
    out=defaultdict(list)
    for line in (data_root/"labels/external/go/rfam2go").read_text(encoding="utf-8",errors="replace").splitlines():
        if line.startswith("!"): continue
        m=re.match(r"Rfam:(RF\d+)\s+.*?;\s*(GO:\d+)",line)
        if m: out[m.group(1)].append(m.group(2))
    return out


def _attrs(text):
    out={}
    for part in (text or "").split(";"):
        if "=" in part:
            k,v=part.split("=",1); out[k.strip()]=v.strip().strip('"')
    return out


def _gff_lines(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)==9: yield f


def _infer_rfam(typ, text):
    low=(text or "").lower()
    if typ=="rRNA":
        for k,rf in RRNA_RFAM.items():
            if k.lower() in low: return rf
        return TYPE_RFAM["rRNA"]
    if typ=="tRNA": return TYPE_RFAM["tRNA"]
    if typ=="tmRNA": return TYPE_RFAM["tmRNA"]
    return ""


def build_for_genome(data_root: Path, genome: str, out_dir: Path):
    out_dir.mkdir(parents=True,exist_ok=True)
    names,ns=_go_meta(data_root); rfam_map=load_rfam2go(data_root)
    rows=[]; source=""
    gffs=[data_root/"annotations"/f"{genome}.gff.gz", data_root/"annotations"/f"{genome}.gff"]
    for p in gffs:
        if not p.exists(): continue
        source=f"GFF:{p.name}"
        for f in _gff_lines(p):
            typ=f[2]
            if typ not in {"rRNA","tRNA","ncRNA","tmRNA"}: continue
            at=_attrs(f[8]); txt=" ".join([at.get("Name",""),at.get("product",""),at.get("comment","")])
            rf=""
            m=re.search(r"Rfam:(RF\d+)",at.get("Dbxref","")+at.get("Ontology_term",""))
            if m: rf=m.group(1)
            if not rf: rf=_infer_rfam(typ,txt)
            if not rf or rf not in rfam_map: continue
            fid=at.get("ID") or at.get("locus_tag") or f"{f[0]}:{f[3]}-{f[4]}"
            for go in rfam_map[rf]:
                rows.append({"genome":genome,"feature_id":fid,"contig":f[0],"start":f[3],"end":f[4],"strand":f[6],"type":typ,"family":at.get("Name","") or at.get("product",""),"rfam":rf,"go_id":go,"go_name":names.get(go,""),"aspect":ns.get(go,"")[:1].upper() if ns.get(go) else "","evidence":"IEA","source":"Rfam2GO:GO_REF:0000115","confidence":"medium"})
        break
    if not rows:
        gbff_candidates=[data_root/"annotations"/f"{genome}.gbff.gz"]
        if genome=="ecoli_k12":
            gbff_candidates.append(data_root/"test_ecoli.gbff")
        for p in [x for x in gbff_candidates if x.exists()]:
            from Bio import SeqIO
            source=f"GBFF:{p.name}"
            for rec in SeqIO.parse(str(p),"genbank"):
                for f in rec.features:
                    typ=f.type
                    if typ not in {"rRNA","tRNA","ncRNA","tmRNA"}: continue
                    q=f.qualifiers; txt=" ".join(q.get("product",[])+q.get("note",[]))
                    rf=""
                    m=re.search(r"Rfam:(RF\d+)"," ".join(q.get("db_xref",[])))
                    if m: rf=m.group(1)
                    if not rf: rf=_infer_rfam(typ,txt)
                    if not rf or rf not in rfam_map: continue
                    fid=(q.get("locus_tag") or [f"{rec.id}:{int(f.location.start)+1}"])[0]
                    for go in rfam_map[rf]:
                        rows.append({"genome":genome,"feature_id":fid,"contig":rec.id,"start":str(int(f.location.start)+1),"end":str(int(f.location.end)),"strand":"+" if f.location.strand==1 else "-","type":typ,"family":(q.get("product") or [""])[0],"rfam":rf,"go_id":go,"go_name":names.get(go,""),"aspect":ns.get(go,"")[:1].upper() if ns.get(go) else "","evidence":"IEA","source":"Rfam2GO:GO_REF:0000115","confidence":"medium"})
            break
    with (out_dir/"ncrna_go.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["genome","feature_id","contig","start","end","strand","type","family","rfam","go_id","go_name","aspect","evidence","source","confidence"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"genome":genome,"n_rows":len(rows),"n_features":len({r["feature_id"] for r in rows}),"by_type":{},"source":source}
    for r in rows: summary["by_type"][r["type"]]=summary["by_type"].get(r["type"],0)+1
    (out_dir/"ncrna_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    return {"out_dir":str(out_dir),"summary":summary}
