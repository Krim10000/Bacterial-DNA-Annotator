"""Anotación GO MF/CC por accesión UniProt (índice SQLite)."""
from __future__ import annotations
import csv, json, sqlite3
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path):
    if not path.exists(): return []
    with path.open(encoding="utf-8",errors="replace") as fh:
        return list(csv.DictReader(fh,delimiter="\t"))


def load_index(data_root: Path):
    """Devuelve {acc: [GO rows]} desde el SQLite local."""
    db=data_root/"labels"/"go_mf_cc"/"goa.sqlite"
    idx=defaultdict(list)
    if not db.exists(): return idx
    con=sqlite3.connect(db)
    for r in con.execute("SELECT acc,go_id,go_name,aspect,namespace,evidence,evidence_class,source,taxon,assigned_by FROM go"):
        idx[r[0]].append(r)
    return idx


def _bsu_acc_by_gene(data_root: Path):
    """B. subtilis: cds_* (functional) -> gen BSGatlas por coordenadas."""
    import json, re
    fp=data_root.parent/"results/subtilis/subtilis.functional_annotation.json"
    if not fp.exists(): return {}
    hits=json.loads(fp.read_text()).get("hits",{})
    gff=data_root.parent/"results/subtilis/subtilis.gff3"
    genes=read_tsv(data_root/"regulation/subtilis/genes.tsv")
    by_coord={(g.get("contig"),g.get("start"),g.get("end"),g.get("strand")):g.get("gene_id")
              for g in genes}
    out={}
    for line in gff.read_text(encoding="utf-8",errors="replace").splitlines():
        if line.startswith("#"): continue
        f=line.split("\t")
        if len(f)!=9 or f[2]!="CDS": continue
        m=re.search(r"ID=([^;]+)",f[8])
        if not m: continue
        gid=by_coord.get((f[0],f[3],f[4],f[6]))
        if gid:
            acc=(hits.get(m.group(1),{}) or {}).get("uniprot","")
            if acc: out[gid]=acc
    return out


def build_for_genome(data_root: Path, genome: str, out_dir: Path):
    out_dir.mkdir(parents=True,exist_ok=True)
    idx=load_index(data_root)
    genes=read_tsv(data_root/"regulation"/genome/"genes.tsv")
    bsu=_bsu_acc_by_gene(data_root) if genome=="subtilis" else {}
    rows=[]
    for g in genes:
        acc=(g.get("acc") or bsu.get(g.get("gene_id"),"") or "").strip().split("-")[0]
        for x in idx.get(acc,[]):
            rows.append({"genome":genome,"gene_id":g.get("gene_id",""),
                         "gene_name":g.get("gene_name",""),"acc":acc,
                         "go_id":x[1],"go_name":x[2],"aspect":x[3],"namespace":x[4],
                         "evidence":x[5],"evidence_class":x[6],"source":x[7],
                         "taxon":x[8],"assigned_by":x[9]})
    with (out_dir/"gene_go.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=["genome","gene_id","gene_name","acc","go_id","go_name",
                                        "aspect","namespace","evidence","evidence_class",
                                        "source","taxon","assigned_by"],delimiter="\t")
        w.writeheader(); w.writerows(rows)
    summary={"genome":genome,"n_genes":len(genes),
             "n_genes_annotated":len({r["gene_id"] for r in rows}),
             "n_rows":len(rows),"n_mf":sum(1 for r in rows if r["aspect"]=="F"),
             "n_cc":sum(1 for r in rows if r["aspect"]=="C"),
             "sources":dict(sorted(Counter(r["source"] for r in rows).items()))}
    (out_dir/"go_summary.json").write_text(json.dumps(summary,indent=1),encoding="utf-8")
    return {"out_dir":str(out_dir),"summary":summary,"gene_go":str(out_dir/"gene_go.tsv")}
