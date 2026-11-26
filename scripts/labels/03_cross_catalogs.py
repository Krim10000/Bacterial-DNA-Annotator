#!/usr/bin/env python3
"""Cruza EC de genes con Rhea, GO, KEGG reactions y KEGG pathways."""
from __future__ import annotations
import argparse, csv, json, re, sys
from collections import defaultdict, Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def read(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_maps(data: Path):
    ec2rhea = defaultdict(set)
    for r in read(data / "metabolism" / "tsv" / "rhea2ec.tsv"):
        ec = (r.get("ID") or "").strip()
        if ec:
            ec2rhea[ec].add((r.get("RHEA_ID") or "").strip())
    rhea2go = defaultdict(set)
    for r in read(data / "metabolism" / "tsv" / "rhea2go.tsv"):
        rhea2go[(r.get("RHEA_ID") or "").strip()].add((r.get("ID") or "").strip())
    rhea2kegg = defaultdict(set)
    for r in read(data / "metabolism" / "tsv" / "rhea2kegg_reaction.tsv"):
        rhea2kegg[(r.get("RHEA_ID") or "").strip()].add((r.get("ID") or "").strip())
    rhea2path = defaultdict(set)
    for r in read(data / "kegg" / "pathway_reactions.tsv"):
        rid = (r.get("rhea") or "").replace("RHEA:", "").strip()
        if rid:
            rhea2path[rid].add((r.get("pathway_name") or "").strip())
    return ec2rhea, rhea2go, rhea2kegg, rhea2path


def genes_k12(data: Path):
    out = []
    for r in read(data / "regulation" / "ecoli_k12" / "genes.tsv"):
        out.append({"id": r.get("gene_id", ""), "name": r.get("gene_name", ""),
                    "product": r.get("product", ""), "ec": r.get("ec", "")})
    return out


def _bsu_ec_map(data: Path) -> dict[str, str]:
    """Recupera EC del GBFF local y los mapea a IDs cds_* del GFF."""
    from Bio import SeqIO
    gbff = data / "dbtbs" / "subtilis.gbff"
    gff = data.parent / "results" / "subtilis" / "subtilis.gff3"
    ec_coord: dict[tuple[int, int, str], str] = {}
    if gbff.exists():
        for rec in SeqIO.parse(str(gbff), "genbank"):
            for f in rec.features:
                if f.type != "CDS":
                    continue
                ecs = [str(x) for x in (f.qualifiers.get("EC_number") or [])]
                if not ecs:
                    continue
                strand = "+" if f.location.strand == 1 else "-"
                ec_coord[(int(f.location.start) + 1, int(f.location.end), strand)] = ";".join(ecs)
    out: dict[str, str] = {}
    if gff.exists():
        for line in gff.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            m = re.search(r"ID=([^;]+)", parts[8])
            if not m:
                continue
            key = (int(parts[3]), int(parts[4]), parts[6])
            if key in ec_coord:
                out[m.group(1)] = ec_coord[key]
    return out


def genes_bsu(data: Path):
    p = data.parent / "results" / "subtilis" / "subtilis.functional_annotation.json"
    hits = json.loads(p.read_text()).get("hits", {}) if p.exists() else {}
    ecs = _bsu_ec_map(data)
    return [{"id": k, "name": v.get("gene", ""), "product": v.get("product", ""),
             "ec": v.get("ec", "") or ecs.get(k, "")} for k, v in hits.items()]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "reports" / "labels_pilot")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    maps = build_maps(a.data_root)
    summary = {}
    for genome, fn in (("ecoli_k12", genes_k12), ("subtilis", genes_bsu)):
        genes = fn(a.data_root)
        rows = []
        for g in genes:
            ecs = [x.strip() for x in (g["ec"] or "").replace(";", ",").split(",") if x.strip()]
            rheads = set(); gos = set(); kegg = set(); paths = set()
            for ec in ecs:
                for rid in maps[0].get(ec, set()):
                    rheads.add(rid); gos |= maps[1].get(rid, set())
                    kegg |= maps[2].get(rid, set()); paths |= maps[3].get(rid, set())
            if ecs or rheads or gos or paths:
                rows.append({**g, "ecs": ";".join(sorted(ecs)),
                             "rhea": ";".join(sorted(rheads)),
                             "go": ";".join(sorted(gos)),
                             "kegg_reaction": ";".join(sorted(kegg)),
                             "pathways": ";".join(sorted(paths)),
                             "n_go": len(gos), "n_pathways": len(paths)})
        summary[genome] = {"n_genes": len(genes), "n_with_ec": sum(1 for g in genes if (g["ec"] or "").strip()),
                           "n_crosswalk": len(rows),
                           "n_with_go": sum(1 for r in rows if r["go"]),
                           "n_with_pathway": sum(1 for r in rows if r["pathways"]),
                           "n_rhea": len({x for r in rows for x in r["rhea"].split(";") if x}),
                           "n_go": len({x for r in rows for x in r["go"].split(";") if x}),
                           "n_pathway": len({x for r in rows for x in r["pathways"].split(";") if x})}
        with (a.out / f"{genome}_catalog_crosswalk.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["id", "name", "product", "ec", "ecs", "rhea", "go", "kegg_reaction", "pathways", "n_go", "n_pathways"], delimiter="\t")
            w.writeheader(); w.writerows(rows)
    (a.out / "crosswalk_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
