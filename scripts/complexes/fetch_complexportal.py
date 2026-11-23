#!/usr/bin/env python3
"""Descarga y mapea un catálogo ComplexPortal para otra especie.

Requiere:
  --taxid          taxid de la especie en ComplexPortal
  --functional     functional.tsv del genoma (qseqid, uniprot en sseqid)
  --genes          genes.tsv del genoma (locus_tag, gene_name)

Salida:
  data/complexes/raw/complexportal_<taxid>.tsv
  data/complexes/gold_<taxid>.tsv

No redistribuye el catálogo; se descarga bajo su licencia CC BY 4.0.
"""
from __future__ import annotations

import argparse
import csv
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import complexes_lib as L  # noqa: E402

ROOT = HERE.parents[1]
RAW = ROOT / "data" / "complexes" / "raw"
OUT = ROOT / "data" / "complexes"
CP_URL = ("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/"
          "complextab/{taxid}.tsv")


def uniprot_map(functional: Path) -> dict[str, str]:
    out = {}
    if not functional.exists():
        return out
    for line in functional.read_text(encoding="utf-8", errors="replace").splitlines():
        c = line.split("\t")
        if len(c) < 2 or "|" not in c[1]:
            continue
        parts = c[1].split("|")
        if len(parts) >= 2:
            out[parts[1]] = c[0]
    return out


def gene_names(genes: Path) -> dict[str, str]:
    out = {}
    if not genes.exists():
        return out
    with genes.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tag = r.get("locus_tag") or r.get("gene_id")
            if tag:
                out[tag] = r.get("gene_name", "")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxid", required=True)
    ap.add_argument("--functional", type=Path, required=True)
    ap.add_argument("--genes", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--allow-incomplete", action="store_true")
    a = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    cp = RAW / f"complexportal_{a.taxid}.tsv"
    if not cp.exists():
        if a.no_download:
            raise SystemExit(f"falta {cp} y --no-download está activo")
        print("descargando", CP_URL.format(taxid=a.taxid))
        urllib.request.urlretrieve(CP_URL.format(taxid=a.taxid), cp)
    acc2loc = uniprot_map(a.functional)
    if not acc2loc:
        # fallback: genes.tsv con columna acc (K-12)
        with a.genes.open(encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("acc"):
                    acc2loc[r["acc"].split("|")[0]] = (
                        r.get("locus_tag") or r.get("gene_id"))
    names = gene_names(a.genes)
    complexes, stats = L.read_complexportal(
        cp, acc2loc, drop_incomplete=not a.allow_incomplete)
    out = a.out or (OUT / f"gold_{a.taxid}.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["complex_id", "complex_name", "subunits", "n_genes",
                    "locus_tags", "gene_names", "evidence_code",
                    "evidence_class", "properties", "source"])
        for c in sorted(complexes, key=lambda x: (-len(x["locus_tags"]), x["complex_id"])):
            w.writerow([c["complex_id"], c["complex_name"], c["subunits"],
                        len(c["locus_tags"]), ",".join(c["locus_tags"]),
                        ",".join(names.get(t, t) for t in c["locus_tags"]),
                        c["evidence_code"], c["evidence_class"],
                        c["properties"], c["source"]])
    print(f"ComplexPortal taxid {a.taxid}: {stats['total']} total, "
          f"{stats['kept']} usables, {stats['unmapped']} sin mapear")
    print("escrito:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
