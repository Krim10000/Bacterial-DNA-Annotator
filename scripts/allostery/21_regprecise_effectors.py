#!/usr/bin/env python3
"""Efectores desde RegPrecise para los organismos del estudio presentes allí.

RegPrecise no tiene exportación machine-readable: se scrapea `genome.jsp?genome_id=N`,
cuya tabla trae `TF | Effector | Pathway | ...`.

De los 400 organismos del estudio, **19** están en RegPrecise (501 genomas). La tabla
resultante es útil más allá de esos 19 porque se indexa por **nombre de TF**: si un
candidato de otro organismo tiene un homólogo llamado igual, se puede transferir con
`status = regprecise_tf_name`.

Salida: data/allostery/work/regprecise_effectors.tsv + resumen JSON.

Uso: python scripts/allostery/21_regprecise_effectors.py
"""
from __future__ import annotations

import json
import re
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "data" / "allostery" / "work"
OUT = WORK / "regprecise_effectors.tsv"
BASE = "https://regprecise.lbl.gov/genome.jsp?genome_id="


def fetch(genome_id: str) -> str:
    req = urllib.request.Request(BASE + genome_id, headers={"User-Agent": "alloreg"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return fh.read().decode("latin-1", errors="ignore")


def parse(html: str) -> list[dict]:
    rows = re.findall(r"<tr[^>]*>(.*?)</tr>", html, re.S)
    out = []
    for r in rows:
        cells = [re.sub(r"<[^>]+>", "", c).strip()
                 for c in re.findall(r"<td[^>]*>(.*?)</td>", r, re.S)]
        if len(cells) < 3 or cells[0] not in ("TF", "RNA"):
            continue
        tf, eff, pathway = cells[1], cells[2], cells[3] if len(cells) > 3 else ""
        if not tf or not eff:
            continue
        out.append({"tf": tf, "efector": eff, "ruta": pathway})
    return out


def main():
    match = json.loads((WORK / "regprecise_match.json").read_text())
    rows = []
    for m in match:
        g = m["rp"]
        try:
            html = fetch(g["genome_id"])
        except Exception as e:
            print(f"  error {g['genome_id']}: {str(e)[:60]}")
            continue
        for r in parse(html):
            rows.append({"organismo": m["org"]["organismo"],
                         "genome_id": g["genome_id"], **r,
                         "fuente": "regprecise"})
        time.sleep(0.3)

    with open(OUT, "w", newline="") as fh:
        cols = ["organismo", "genome_id", "tf", "efector", "ruta", "fuente"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ")
                               for c in cols) + "\n")

    summary = {"genomas_consultados": len(match), "pares": len(rows),
               "tfs_distintos": len({r["tf"].lower() for r in rows}),
               "efectores_distintos": len({r["efector"] for r in rows})}
    (WORK / "regprecise_effectors.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
