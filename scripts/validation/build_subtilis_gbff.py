#!/usr/bin/env python3
"""build_subtilis_gbff.py — GBFF mínimo de *B. subtilis* desde BSGatlas.

El motor determinista (`rule_operons.py`) lee GenBank, no GFF+FASTA. Este script
convierte la anotación BSGatlas v1.1 en un GenBank mínimo que ese lector entiende
(formato: `LOCUS` + `FEATURES` + `ORIGIN`).

Detalle importante: en BSGatlas los `locus_tag` (BSU_XXXXX) están en las
features `gene`, no en las `CDS` (solo 64 de 4.332 los llevan en `synonyms`).
Se emparejan por coordenadas idénticas (start, end, hebra).

Uso
---
    python scripts/validation/build_subtilis_gbff.py
    python scripts/validation/build_subtilis_gbff.py --limit 100   # prueba
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GFF = ROOT / "data" / "dbtbs" / "BSGatlas_v1.1.gff"
FASTA = ROOT / "data" / "dbtbs" / "genome.fna"
OUT = ROOT / "data" / "dbtbs" / "subtilis.gbff"


def attrs(field: str) -> dict[str, str]:
    out = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out


def first_bsu(value: str) -> str:
    for tok in value.replace(";", ",").split(","):
        tok = tok.strip()
        if tok.startswith("BSU_"):
            return tok
    return ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gff", type=Path, default=GFF)
    ap.add_argument("--fasta", type=Path, default=FASTA)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--limit", type=int, default=0, help="solo N CDS (pruebas)")
    a = ap.parse_args()

    # --- coordenadas de la secuencia
    seq_lines = []
    for line in open(a.fasta):
        if line.startswith(">"):
            header = line[1:].split()[0]
            continue
        seq_lines.append(line.strip())
    seq = "".join(seq_lines).upper()
    print(f"secuencia: {header} · {len(seq):,} pb")

    # --- features
    genes_by_coords: dict[tuple[int, int, str], dict] = {}
    genes_by_id: dict[str, dict] = {}
    cds_rows = []
    for line in open(a.gff):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        at = attrs(p[8])
        kind = p[2]
        if kind == "gene":
            info = {"locus": first_bsu(at.get("locus_tag", "")),
                    "name": at.get("Name", ""), "id": at.get("ID", "")}
            genes_by_coords.setdefault((int(p[3]), int(p[4]), p[6]), info)
            if info["id"]:
                genes_by_id[info["id"]] = info
        elif kind == "CDS":
            cds_rows.append({"start": int(p[3]), "end": int(p[4]), "strand": p[6],
                             "name": at.get("Name", ""), "ec": at.get("ec", ""),
                             "syn": at.get("synonyms", ""),
                             "derives": at.get("derives_from", ""),
                             "id": at.get("ID", "")})

    if a.limit:
        cds_rows = cds_rows[:a.limit]

    rows = []
    for c in cds_rows:
        info = genes_by_coords.get((c["start"], c["end"], c["strand"]), {})
        locus = info.get("locus") or first_bsu(c["syn"]) or c["derives"] or c["id"]
        name = c["name"] or info.get("name", "")
        ec = c["ec"].split(",")[0].replace("EC-", "") if c["ec"] else ""
        rows.append({"start": c["start"], "end": c["end"], "strand": c["strand"],
                     "locus": locus, "name": name, "ec": ec})
    rows.sort(key=lambda r: r["start"])
    con_locus = sum(1 for r in rows if r["locus"].startswith("BSU_"))
    print(f"CDS: {len(rows)} · con locus_tag BSU_: {con_locus} "
          f"({100*con_locus/max(1,len(rows)):.0f}%) · con nombre: "
          f"{sum(1 for r in rows if r['name'])}")

    # --- escribir GBFF (formato que espera find_operons.read_gbff)
    out = [f"LOCUS       {header}     {len(seq)} bp  DNA  circular BCT 21-NOV-2026",
           "FEATURES             Location/Qualifiers"]
    for r in rows:
        loc = (f"complement({r['start']}..{r['end']})" if r["strand"] == "-"
               else f"{r['start']}..{r['end']}")
        out.append(f"     CDS             {loc}")
        out.append(f'                     /locus_tag="{r["locus"]}"')
        if r["name"]:
            out.append(f'                     /gene="{r["name"]}"')
        if r["ec"]:
            out.append(f'                     /EC_number="{r["ec"]}"')
    out.append("ORIGIN")
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        out.append(f"{i + 1:>9} {groups}")
    out.append("//")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text("\n".join(out) + "\n")
    print(f"-> {a.out.relative_to(ROOT)} ({a.out.stat().st_size/1e6:.1f} MB)")

    # --- verificación de coherencia
    bad = [r for r in rows if not (1 <= r["start"] < r["end"] <= len(seq))]
    print(f"coordenadas incoherentes: {len(bad)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
