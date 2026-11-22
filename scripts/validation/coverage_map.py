#!/usr/bin/env python3
"""coverage_map.py — qué evidencia tiene cada genoma del proyecto.

Para cada genoma anotado: genes, operones, operones con un terminador clase A
cerca del final y ortólogos de K-12. Sirve para saber qué genomas son mejorables
y con qué, en vez de aplicar una receta única a todos.

Salida: data/operon_finder/coverage_map.tsv
"""
from __future__ import annotations
import csv, glob
from collections import Counter
from pathlib import Path
import bisect

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    hp: dict[str, list[int]] = {}
    for f in glob.glob(str(ROOT / "data/operon_finder/*_hairpins.tsv")):
        g = Path(f).name.replace("_hairpins.tsv", "")
        hp[g] = sorted(int(r["start"]) for r in csv.DictReader(open(f), delimiter="\t")
                       if r.get("class") == "A_fuerte")
    ort = Counter()
    hits = ROOT / ".cache/synteny/hits.tsv"
    if hits.exists():
        for line in open(hits):
            ort[line.split("|")[0]] += 1
    rows = []
    for d in sorted((ROOT / "data/regulation").iterdir()):
        if not d.is_dir() or not (d / "genes.tsv").exists():
            continue
        genes = list(csv.DictReader(open(d / "genes.tsv"), delimiter="\t"))
        ops_path = d / "operons.tsv"
        ops = list(csv.DictReader(open(ops_path), delimiter="\t")) if ops_path.exists() else []
        arr = hp.get(d.name, [])
        term = 0
        for o in ops:
            try:
                e = int(o["end"])
            except (KeyError, ValueError):
                continue
            i = bisect.bisect_left(arr, e - 200)
            if i < len(arr) and arr[i] <= e + 200:
                term += 1
        rows.append({"genome": d.name, "genes": len(genes), "operons": len(ops),
                     "con_terminador": term, "ortologos_k12": ort.get(d.name, 0)})
    out = ROOT / "data/operon_finder/coverage_map.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    tot = sum(r["operons"] for r in rows); ct = sum(r["con_terminador"] for r in rows)
    print(f"genomas: {len(rows)} | operones: {tot:,} | con terminador: {ct:,} ({100*ct/max(1,tot):.0f}%)")
    print(f"-> {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
