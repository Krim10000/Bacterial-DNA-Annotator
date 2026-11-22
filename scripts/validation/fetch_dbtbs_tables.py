#!/usr/bin/env python3
"""fetch_dbtbs_tables.py — DBTBS desde sus tablas públicas.

Por qué este camino
-------------------
El formulario de descarga de DBTBS (`/download.html`) responde con una
redirección a `dbtbs_050922.xml.gz`, un release de 2005 que **ya no existe**
(la URL devuelve 404). Ese enlace está roto en el servidor de DBTBS, no aquí.

Lo que sí funciona son sus tablas públicas:

- `promtable.html` — gen → operón → nº de elementos reguladores y sitios de unión.
- `tfactable.html` — factores de transcripción y sus sitios.

`robots.txt` solo prohíbe `/cgi-bin/`, así que estas páginas son de acceso
legítimo. El registro que exige el formulario se debe a que piden afiliación y
propósito de uso; eso queda declarado en `docs/EXTERNAL_SOURCES.md`.

Salida
------
    data/dbtbs/gt_operons_dbtbs.tsv   operón, genes y coordenadas (si mapean)
    data/dbtbs/gt_tf_sites_dbtbs.tsv  sitios de unión TF→gen

Uso
---
    python scripts/validation/fetch_dbtbs_tables.py
"""
from __future__ import annotations

import re
import sys
import urllib.request
from collections import OrderedDict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "dbtbs"
GFF = OUT / "BSGatlas_v1.1.gff"
BASE = "https://dbtbs.hgc.jp/"
UA = {"User-Agent": "Mozilla/5.0 (academic use; operon benchmark)"}


def fetch(page: str) -> str:
    req = urllib.request.Request(BASE + page, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.read().decode("latin-1", errors="replace")


def cells(row: str) -> list[str]:
    out = []
    for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S):
        txt = re.sub(r"<[^>]*>", "", c)
        txt = re.sub(r"&nbsp;", " ", txt)
        out.append(re.sub(r"\s+", " ", txt).strip())
    return out


def parse_promtable(html: str) -> list[dict]:
    """Gen → operón. La tabla grande es la única con miles de filas."""
    tables = re.findall(r"<table[^>]*>(.*?)</table>", html, re.S)
    if not tables:
        raise RuntimeError("promtable.html sin tablas")
    big = max(tables, key=len)
    rows = re.findall(r"<tr[^>]*>(.*?)</tr>", big, re.S)
    out = []
    for r in rows:
        c = cells(r)
        if len(c) < 3 or not c[0] or c[0].lower() == "gene":
            continue
        gene, synonym, operon = c[0], c[1], c[2]
        if not operon:
            continue
        try:
            n_sites = int(c[5]) if len(c) > 5 and c[5].isdigit() else 0
        except IndexError:
            n_sites = 0
        out.append({"gene": gene, "synonym": synonym, "operon": operon,
                    "n_sites": n_sites})
    return out


def gff_genes() -> dict[str, dict]:
    """Name= → coordenadas, desde el GFF de BSGatlas (misma especie)."""
    out: dict[str, dict] = {}
    if not GFF.exists():
        return out
    for line in open(GFF, encoding="utf-8", errors="replace"):
        if line.startswith("#") or "\t" not in line:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "CDS"):
            continue
        m = re.search(r"Name=([^;]+)", f[8])
        if not m:
            continue
        name = m.group(1).strip()
        if name in out:          # el feature `gene` manda
            continue
        out[name] = {"contig": f[0], "start": int(f[3]), "end": int(f[4]),
                     "strand": f[6]}
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print("descargando promtable.html…")
    rows = parse_promtable(fetch("promtable.html"))
    print(f"  {len(rows)} genes con operón asignado")

    genes = gff_genes()
    print(f"  {len(genes)} genes en el GFF de BSGatlas para mapear")

    ops: "OrderedDict[str, dict]" = OrderedDict()
    for r in rows:
        o = ops.setdefault(r["operon"], {"operon": r["operon"], "genes": [],
                                         "n_sites": 0})
        o["genes"].append(r["gene"])
        o["n_sites"] += r["n_sites"]

    mapped_total = 0
    for o in ops.values():
        coords = []
        for g in o["genes"]:
            hit = genes.get(g)
            if hit:
                coords.append(hit)
        mapped_total += len(coords)
        o["n_genes"] = len(o["genes"])
        o["n_mapped"] = len(coords)
        if coords:
            o["contig"] = coords[0]["contig"]
            o["start"] = min(c["start"] for c in coords)
            o["end"] = max(c["end"] for c in coords)
            o["strand"] = coords[0]["strand"]
        else:
            o["contig"] = o["start"] = o["end"] = ""
            o["strand"] = ""

    out1 = OUT / "gt_operons_dbtbs.tsv"
    with open(out1, "w", newline="") as fh:
        cols = ["operon", "n_genes", "n_mapped", "contig", "start", "end",
                "strand", "n_sites", "genes"]
        fh.write("\t".join(cols) + "\n")
        for o in ops.values():
            fh.write("\t".join(str(o.get(c, "")) for c in cols[:-1]) + "\t"
                     + ",".join(o["genes"]) + "\n")
    multi = sum(1 for o in ops.values() if o["n_genes"] > 1)
    print(f"  {len(ops)} operones ({multi} multi-génicos) -> "
          f"{out1.relative_to(ROOT)}")
    print(f"  genes mapeados a coordenadas: {mapped_total}/{len(rows)} "
          f"({100*mapped_total/max(1,len(rows)):.0f}%)")

    # sitios de unión de factores
    try:
        print("descargando tfactable.html…")
        html = fetch("tfactable.html")
        tables = re.findall(r"<table[^>]*>(.*?)</table>", html, re.S)
        big = max(tables, key=len)
        rws = re.findall(r"<tr[^>]*>(.*?)</tr>", big, re.S)
        sites = []
        for r in rws:
            c = cells(r)
            if len(c) >= 4 and c[0] and c[0].lower() not in ("name", "factor"):
                sites.append(c)
        out2 = OUT / "gt_tf_sites_dbtbs.tsv"
        with open(out2, "w", newline="") as fh:
            for c in sites:
                fh.write("\t".join(c) + "\n")
        print(f"  {len(sites)} filas de factores -> {out2.relative_to(ROOT)}")
    except Exception as exc:  # noqa: BLE001
        print(f"  (tfactable no parseable: {exc})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
