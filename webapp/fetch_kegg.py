#!/usr/bin/env python3
"""Descarga el catálogo de rutas KEGG y sus reacciones (una vez).

    python -m webapp.fetch_kegg --data-root data

Genera data/kegg/pathway_reactions.tsv con (pathway_id, name, rhea).
Uso académico; respeta los términos de KEGG.
"""
from __future__ import annotations

import argparse
import csv
import time
import urllib.request
from pathlib import Path

BASE = "https://rest.kegg.jp"


def get(path: str, retries: int = 3) -> str:
    for i in range(retries):
        try:
            with urllib.request.urlopen(BASE + path, timeout=30) as r:
                return r.read().decode()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", type=Path, default=None)
    a = ap.parse_args()
    out_dir = a.out_dir or (a.data_root / "kegg")
    out_dir.mkdir(parents=True, exist_ok=True)

    # reacción KEGG -> Rhea
    k2r: dict[str, str] = {}
    master: dict[str, str] = {}
    rhea2kegg = a.data_root / "reactions" / "tsv" / "rhea2kegg_reaction.tsv"
    for r in csv.DictReader(open(rhea2kegg), delimiter="\t"):
        k2r[r["ID"]] = f"RHEA:{r['RHEA_ID']}"
        master[f"RHEA:{r['RHEA_ID']}"] = f"RHEA:{r['MASTER_ID']}"

    pathways = []
    for line in get("/list/pathway").strip().split("\n"):
        p = line.split("\t")
        if len(p) >= 2 and p[0].startswith("map"):
            pathways.append((p[0], p[1]))
    # metabolismos (excluye mapas globales)
    keep = [x for x in pathways if x[0].startswith(("map00", "map01", "map012"))
            and x[0] not in ("map01100", "map01110", "map01120", "map01130",
                             "map01200", "map01210", "map01220", "map01230", "map01240")]
    print(f"rutas KEGG: {len(keep)} de {len(pathways)}")

    rows = []
    for i, (pid, name) in enumerate(keep, 1):
        try:
            text = get(f"/link/reaction/{pid}")
        except Exception as exc:
            print(f"  {pid}: {exc}")
            continue
        for line in text.strip().split("\n"):
            p = line.split("\t")
            if len(p) < 2:
                continue
            rxn = p[1].replace("rn:", "")
            rhea = k2r.get(rxn)
            if rhea:
                rows.append({"pathway_id": pid, "pathway_name": name,
                             "rhea": master.get(rhea, rhea)})
        if i % 25 == 0:
            print(f"  {i}/{len(keep)} · reacciones: {len(rows):,}")

    # deduplicar (ruta, reacción maestra)
    rows = list({(r["pathway_id"], r["rhea"]): r for r in rows}.values())
    out = out_dir / "pathway_reactions.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pathway_id", "pathway_name", "rhea"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"listo: {len(rows):,} reacciones mapeadas -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
