#!/usr/bin/env python3
"""Capa 1 (paso 1c) — mapeo matriz→TF resolviendo el sistema de coordenadas.

Problema: las posiciones de las matrices de Reg-Seq son coordenadas del **amplicón**
(posición 0 = extremo del fragmento mutado), mientras que el suplemento 3 del paper da
los sitios en coordenadas **relativas al TSS**. De ahí que el solapamiento directo fallara.

Solución: alinear la secuencia de diseño (`prior_designs/wtsequences.csv`) contra el
genoma de E. coli K-12 y calcular el desplazamiento por promotor. Empíricamente:
−113 pb (hebra +) y −116 pb (hebra −), consistente en 107 promotores.

Con eso: `coord_TSS = pos_matriz + offset`, y se asigna TF por solapamiento con los
sitios del suplemento 3.

Salida: data/allostery/regseq_tf_map.tsv + resumen JSON.

Uso: python scripts/allostery/15_regseq_map_offsets.py
"""
from __future__ import annotations

import collections
import csv
import json
import re
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery"
SUPP = OUT / "work" / "regseq" / "supp3_sites.xlsx"
WT = "https://raw.githubusercontent.com/RPGroup-PBoC/RegSeq/master/data/prior_designs/wtsequences.csv"
GBFF = ROOT / "data" / "test_ecoli.gbff"
COMP = str.maketrans("acgtn", "tgcan")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def genome() -> str:
    txt = GBFF.read_text(errors="ignore")
    i, j = txt.find("ORIGIN"), txt.find("//", txt.find("ORIGIN"))
    return re.sub(r"[^acgtn]", "", txt[i:j].lower())


def offsets(g: str) -> dict[str, int]:
    data = urlopen(WT, timeout=60).read().decode()
    out = {}
    for r in csv.DictReader(data.splitlines()):
        seq = (r.get("geneseq") or "").lower()
        if len(seq) < 20:
            continue
        ss = int(float(r["start_site"]))
        rev = r["rev"].strip().lower() == "rev"
        q = g.find((rc(seq) if rev else seq)[:40])
        if q < 0:
            continue
        out[r["name"].lower()] = (ss - q - len(seq) + 1) if rev else (q - ss)
    return out


def load_sites():
    df = pd.read_excel(SUPP, sheet_name="3_TF_binding_sites", header=None)
    hdr = next(i for i in range(20) if str(df.iloc[i, 0]).strip().lower() == "gene")
    df = df.iloc[hdr + 1:].copy()
    df.columns = ["gene", "start", "end", "type", "identity"]
    sites = collections.defaultdict(list)
    for r in df.itertuples():
        if str(r.identity) == "nan" or pd.isna(r.identity):
            continue
        sites[str(r.gene).lower()].append((int(float(r.start)), int(float(r.end)), r.identity))
    return sites


def main():
    g = genome()
    offs = offsets(g)
    sites = load_sites()
    rows = list(csv.DictReader(open(OUT / "regseq_matrices.tsv"), delimiter="\t"))
    out = []
    for m in rows:
        p = m["promotor"].lower()
        key = p if p in offs else next((k for k in offs if k in p or p in k), None)
        rec = {**m, "offset": None, "tf_suplemento": "", "n_tf": 0}
        if key and key in sites and key in offs:
            off = offs[key]
            lo, hi = int(m["from"]) + off, int(m["to"]) + off
            hs = [s for s in sites[key] if s[0] <= hi and s[1] >= lo]
            names = sorted({h[2] for h in hs})
            rec.update({"offset": off, "tf_suplemento": "|".join(names), "n_tf": len(names)})
        out.append(rec)

    with open(OUT / "regseq_tf_map.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    uniq = [x for x in out if x["n_tf"] == 1]
    summary = {
        "n_matrices": len(out),
        "con_sitio": sum(1 for x in out if x["n_tf"] > 0),
        "tf_unico": len(uniq),
        "ambiguo": sum(1 for x in out if x["n_tf"] > 1),
        "tfs_distintos": len({x["tf_suplemento"] for x in uniq}),
        "offsets": collections.Counter(x["offset"] for x in out if x["offset"] is not None).most_common(5),
        "top_tfs": collections.Counter(x["tf_suplemento"] for x in uniq).most_common(15),
    }
    (OUT / "regseq_tf_map.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
