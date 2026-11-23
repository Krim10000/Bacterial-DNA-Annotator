#!/usr/bin/env python3
"""Audita complejos por vecindad contra las brechas de los 400 genomas.

Uso:
    python tools/audit_complexes_400.py                 # todos
    python tools/audit_complexes_400.py --limit 20      # muestra
    python tools/audit_complexes_400.py --out reports/complex_audit_400.tsv

No modifica las redes ni cierra brechas. Genera una tabla de auditoría y, si
se pide, `data/complexes/<genoma>.tsv` para el navegador.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "webapp"))
from complexes import infer, write  # noqa: E402
from functional import parse_hit  # noqa: E402


def parse_gff(path: Path) -> list[dict]:
    out = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        attrs = {}
        for part in c[8].split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                attrs[k] = v
        out.append({
            "id": attrs.get("ID", ""), "type": "CDS",
            "contig": c[0], "start": int(c[3]), "end": int(c[4]),
            "strand": c[6],
        })
    return out


def parse_functional(path: Path) -> dict[str, dict]:
    hits = {}
    if not path.exists():
        return hits
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        c = line.rstrip("\n").split("\t")
        if len(c) < 7:
            continue
        stitle = c[-1] if len(c) >= 13 else ""
        try:
            h = parse_hit(c[0], c[1], stitle, float(c[2]), int(c[3]),
                          float(c[10]) if len(c) > 10 else 1.0,
                          float(c[11]) if len(c) > 11 else 0.0)
        except ValueError:
            continue
        hits[c[0]] = {
            "cds_id": h.cds_id, "uniprot": h.uniprot, "gene": h.gene,
            "product": h.product, "ec": h.ec, "organism": h.organism,
            "pident": h.pident, "aln_len": h.aln_len, "evalue": h.evalue,
            "bitscore": h.bitscore, "path": h.path,
        }
    return hits


def load_bridge(data_root: Path):
    with open(data_root / "metabolism" / "rhea_reactions.pkl", "rb") as fh:
        reactions = pickle.load(fh)
    ec_products: dict[str, set[str]] = defaultdict(set)
    with open(data_root / "reactions" / "tsv" / "rhea2ec.tsv",
              encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or not c[3]:
                continue
            rid = int(c[0]) if c[0].isdigit() else None
            r = reactions.get(rid or -1, {})
            ec_products[c[3]] |= set(r.get("products", set()))
    return ec_products


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "reports" / "complex_audit_400.tsv")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--write-complexes", action="store_true",
                    help="escribe data/complexes/<genoma>.tsv para el navegador")
    a = ap.parse_args()

    ec_products = load_bridge(a.data_root)
    annotations = a.data_root / "annotations_ab_initio"
    genomes = sorted(p.name for p in annotations.iterdir() if p.is_dir())
    if a.limit:
        genomes = genomes[:a.limit]

    gaps_by_genome: dict[str, set[str]] = defaultdict(set)
    gap_path = a.data_root / "metabolism" / "gap_interpretation.tsv"
    if gap_path.exists():
        with gap_path.open(encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("chebi"):
                    gaps_by_genome[row.get("genome", "")].add(row["chebi"])

    rows = []
    for i, g in enumerate(genomes, 1):
        gdir = annotations / g
        feats = parse_gff(gdir / "annotation.gff3") if (gdir / "annotation.gff3").exists() else []
        hits = parse_functional(gdir / "functional.tsv")
        if not feats:
            continue
        result = infer(feats, hits, g)
        if a.write_complexes and result["complexes"]:
            outdir = a.data_root / "complexes"
            write(result, outdir)
        gaps = gaps_by_genome.get(g, set())
        matched: set[str] = set()
        complete: set[str] = set()
        incomplete: set[str] = set()
        for c in result["complexes"]:
            products: set[str] = set()
            for ec in str(c.get("ecs", "")).split(";"):
                if ec:
                    products |= ec_products.get(ec, set())
            hit = gaps & products
            if hit:
                matched |= hit
                if c["status"] == "complex_complete":
                    complete |= hit
                elif c["status"] == "complex_incomplete":
                    incomplete |= hit
        rows.append({
            "genome": g, "gaps": len(gaps),
            "complexes": len(result["complexes"]),
            "complete": result["n_complete"],
            "incomplete": result["n_incomplete"],
            "gaps_with_complex": len(matched),
            "gaps_complete_complex": len(complete),
            "gaps_incomplete_complex": len(incomplete),
            "gaps_unassociated": len(gaps - matched),
        })
        if i % 25 == 0:
            print(f"{i}/{len(genomes)}", flush=True)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0]) if rows else ["genome"]
    with a.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    summary = {k: sum(r[k] for r in rows) for k in
               ("gaps", "complexes", "complete", "incomplete",
                "gaps_with_complex", "gaps_complete_complex",
                "gaps_incomplete_complex", "gaps_unassociated")}
    (a.out.with_suffix(".json")).write_text(
        json.dumps({"summary": summary, "genomes": len(rows)},
                   indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
