#!/usr/bin/env python3
"""fetch_bsgatlas.py — segunda especie curada: Bacillus subtilis 168 (BSGatlas).

Por qué BSGatlas y no DBTBS
---------------------------
DBTBS es la base curada clásica de *B. subtilis*, pero su descarga **no es
automatizable**: la página `/download.html` es un formulario que exige nombre,
email, institución y propósito, y lo envía a `/cgi-bin/download.cgi`; además
`robots.txt` prohíbe `/cgi-bin/`. Rellenar ese formulario con datos inventados
sería falsificar un registro, así que no se hace.

BSGatlas (Geissler et al.) es una integración abierta y publicada que **incluye
DBTBS** entre sus fuentes, y ofrece descarga directa en GFF3. Es una segunda
especie curada legítima y, de paso, trae terminadores con energía libre y TSS
con factor sigma, que es justo lo que faltaba para validar fuera de *E. coli*.

Fuente: https://rth.dk/resources/bsgatlas/  (v1.1)
Genoma: Bacillus subtilis 168, NCBI ASM904v1 (NC_000964.3)

Salidas en `data/dbtbs/`:
    BSGatlas_v1.1.gff      GFF3 original (caché)
    genome.fna             genoma
    gt_operons.tsv         operones con genes (mismo esquema que RegulonDB)
    gt_terminators.tsv     terminadores con energía
    gt_tss.tsv             sitios de inicio con sigma

Uso
---
    python scripts/validation/fetch_bsgatlas.py
    python scripts/validation/fetch_bsgatlas.py --offline   # solo parsea
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "dbtbs"
BASE = "https://rth.dk/resources/bsgatlas/annotation/"
FILES = {"gff": "BSGatlas_v1.1.gff", "fna": "genome.fna"}
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64)"}


def fetch(name: str, offline: bool) -> Path:
    path = OUT / name
    if path.exists() and path.stat().st_size > 1000:
        return path
    if offline:
        raise SystemExit(f"falta {path} y --offline está activo")
    OUT.mkdir(parents=True, exist_ok=True)
    url = BASE + name
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as fh, open(path, "wb") as out:
        out.write(fh.read())
    print(f"  descargado {name}: {path.stat().st_size / 1e6:.1f} MB")
    return path


def attrs(field: str) -> dict[str, str]:
    out = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out


def parse(gff: Path) -> tuple[list[dict], list[dict], list[dict]]:
    operons: dict[str, dict] = {}
    genes_of: dict[str, list[dict]] = defaultdict(list)
    terms, tss = [], []

    for line in open(gff):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        ftype, start, end, strand, a = p[2], int(p[3]), int(p[4]), p[6], attrs(p[8])
        if ftype == "operon":
            operons[a.get("ID", "")] = {
                "operon_id": a.get("ID", ""), "strand": strand,
                "start": start, "end": end,
                "locus": p[0], "name": a.get("Name", "")}
        elif ftype == "gene":
            for par in a.get("Parent", "").split(","):
                if par.startswith("BSGatlas-operon"):
                    genes_of[par].append({
                        # el atributo trae alias separados por coma:
                        # locus_tag=BSU_00010,BSU00010 -> se usa el primero
                        "locus_tag": a.get("locus_tag", "").split(",")[0],
                        "name": a.get("Name", ""),
                        "start": start, "strand": strand})
        elif ftype == "terminator":
            m = re.search(r"Energy:\s*(-?[\d.]+)", a.get("comment", ""))
            terms.append({"start": start, "end": end, "strand": strand,
                          "energy": m.group(1) if m else "",
                          "parent": a.get("Parent", ""),
                          "evidence": a.get("comment", "")[:160]})
        elif ftype == "TSS":
            m = re.search(r"Sigma:\s*([^,;]*)", a.get("comment", ""))
            tss.append({"pos": start, "strand": strand,
                        "sigma": (m.group(1).strip() if m else ""),
                        "evidence": a.get("comment", "")[:160]})

    rows = []
    for oid, op in operons.items():
        # BSGatlas trae dos features `gene` por gen (dos fuentes); se deduplica
        # por coordenada de inicio, que es única por gen
        seen, uniq = set(), []
        for g in sorted(genes_of.get(oid, []), key=lambda g: g["start"]):
            if g["start"] in seen:
                continue
            seen.add(g["start"]); uniq.append(g)
        genes = uniq
        rows.append({
            "operon_id": oid,
            "operon_name": op["name"] or oid,
            "strand": op["strand"], "start": op["start"], "end": op["end"],
            "n_genes": len(genes),
            "genes": ",".join(g["locus_tag"] for g in genes),
            "gene_names": ",".join(g["name"] for g in genes),
            "n_genes_tus": len(genes),
            "mapping_ok": "yes" if genes else "no",
            "confidence": "", "evidence": "BSGatlas v1.1 (integra DBTBS, BsubCyc, Rfam, RefSeq, SubtiWiki)",
        })
    rows.sort(key=lambda r: (r["strand"], r["start"]))
    return rows, terms, tss


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args()

    print("BSGatlas (Bacillus subtilis 168):")
    gff = fetch(FILES["gff"], a.offline)
    fetch(FILES["fna"], a.offline)

    ops, terms, tss = parse(gff)
    write(OUT / "gt_operons.tsv", ops)
    write(OUT / "gt_terminators.tsv", terms)
    write(OUT / "gt_tss.tsv", tss)

    multi = [o for o in ops if o["n_genes"] >= 2]
    print(f"  operones: {len(ops)} (con >=2 genes: {len(multi)}, "
          f"sin genes: {sum(1 for o in ops if not o['genes'])})")
    print(f"  terminadores: {len(terms)} (con energía: "
          f"{sum(1 for t in terms if t['energy'])})")
    print(f"  TSS: {len(tss)} (con sigma: {sum(1 for t in tss if t['sigma'])})")
    print(f"  -> {(OUT / 'gt_operons.tsv').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
