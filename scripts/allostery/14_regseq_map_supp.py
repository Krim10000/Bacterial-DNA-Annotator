#!/usr/bin/env python3
"""Capa 1 (paso 1b) — tabla promotor→TF desde el suplemento de Reg-Seq.

El suplemento 3 del paper (Ireland et al. 2020) trae los sitios de TF identificados:
columnas `gene, start, end, type, identity`. Con eso se asigna cada matriz de energía a
su TF por **solapamiento de ventana**, que es la vía correcta (el nombre del archivo solo
trae promotor y coordenadas).

Salida: data/allostery/regseq_tf_map.tsv + resumen.

Uso: python scripts/allostery/14_regseq_map_supp.py
"""
from __future__ import annotations

import collections
import csv
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery"
SUPP = OUT / "work" / "regseq" / "supp3_sites.xlsx"


def load_sites() -> list[dict]:
    df = pd.read_excel(SUPP, sheet_name="3_TF_binding_sites", header=None)
    header_row = None
    for i in range(min(20, len(df))):
        if str(df.iloc[i, 0]).strip().lower() == "gene":
            header_row = i
            break
    if header_row is None:
        raise SystemExit("no se encontró la cabecera 'gene'")
    df = df.iloc[header_row + 1:].copy()
    df.columns = ["gene", "start", "end", "type", "identity"]
    sites = []
    for _, r in df.iterrows():
        try:
            sites.append({"gene": str(r["gene"]).strip(),
                          "start": int(float(r["start"])), "end": int(float(r["end"])),
                          "tipo": str(r["type"]).strip(),
                          "tf": (None if pd.isna(r["identity"]) else str(r["identity"]).strip())})
        except Exception:
            continue
    return sites


def match_promoter(name: str, genes: set[str]) -> str | None:
    n = name.lower()
    if n in genes:
        return n
    # nombres compuestos tipo 'WaaA-coaD': probar cada parte
    for part in n.replace("-", "_").split("_"):
        if part in genes:
            return part
    # último recurso: contención
    for g in genes:
        if g and (g in n or n in g):
            return g
    return None


def main():
    sites = load_sites()
    genes = {s["gene"].lower() for s in sites}
    print(f"sitios en el suplemento: {len(sites)} · promotores: {len(genes)}")

    rows = list(csv.DictReader(open(OUT / "regseq_matrices.tsv"), delimiter="\t"))
    out = []
    for r in rows:
        gene = match_promoter(r["promotor"], genes)
        hits = []
        if gene:
            lo, hi = int(r["from"]), int(r["to"])
            for s in sites:
                if s["gene"].lower() != gene or not s["tf"]:
                    continue
                if s["start"] <= hi and s["end"] >= lo:      # solapamiento
                    hits.append(s)
        tfs = sorted({h["tf"] for h in hits})
        out.append({**r, "gen_suplemento": gene,
                    "tf_suplemento": "|".join(tfs) if tfs else "",
                    "n_sitios_solapan": len(hits),
                    "n_tf_solapan": len(tfs)})

    with open(OUT / "regseq_tf_map.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    uniq = [x for x in out if x["n_tf_solapan"] == 1]
    amb = [x for x in out if x["n_tf_solapan"] > 1]
    none = [x for x in out if x["n_tf_solapan"] == 0]
    summary = {
        "n_matrices": len(out),
        "promotor_en_suplemento": sum(1 for x in out if x["gen_suplemento"]),
        "tf_unico": len(uniq),
        "tf_ambiguo": len(amb),
        "sin_tf_solapante": len(none),
        "tfs_distintos": len({x["tf_suplemento"] for x in uniq}),
        "top_tfs": collections.Counter(x["tf_suplemento"] for x in uniq).most_common(15),
    }
    (OUT / "regseq_tf_map.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
