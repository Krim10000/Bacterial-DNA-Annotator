#!/usr/bin/env python3
"""compare_batter.py — BATTER-TPE frente a los terminadores curados de RegulonDB.

Qué mide y qué no
-----------------
- **Sensibilidad**: de los terminadores curados de K-12, cuántos tienen una
  predicción de BATTER dentro de ±tolerancia. Es una medida limpia: el curado es
  el denominador.
- **Precisión (cota inferior)**: de las predicciones de BATTER, cuántas caen
  cerca de un terminador curado. Ojo: RegulonDB **no es exhaustivo**, así que
  muchas predicciones "fuera" pueden ser terminadores reales no curados. El
  valor que sale es un suelo, no la precisión real.

Las posiciones curadas de RegulonDB pueden venir en coordenadas del genoma
completo; por eso se comparan directamente con el BED de BATTER, sin conversión.

Uso
---
    python scripts/regulation/compare_batter.py \\
        --bed data/terminators/batter_k12_300k.bed \\
        --gt data/regulondb/TerminatorSet.tsv --tol 25
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load_gt(path: Path, limit: int | None = None) -> list[dict]:
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                a, b = int(r["2)leftEndPos"]), int(r["3)rightEndPos"])
            except (KeyError, ValueError):
                continue
            if limit and a > limit:
                continue
            rows.append({"start": min(a, b), "end": max(a, b),
                         "strand": r.get("4)strand", ""),
                         "type": r.get("7)type", ""),
                         "id": r.get("1)id", "")})
    return rows


def load_bed(path: Path) -> list[dict]:
    """Acepta BED de 3 columnas o con score/strand; ignora comentarios."""
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith(("#", "track", "browser")):
            continue
        p = line.split("\t")
        if len(p) < 3:
            p = line.split()
        try:
            start, end = int(float(p[1])), int(float(p[2]))
        except (IndexError, ValueError):
            continue
        score = None
        if len(p) > 3:
            try:
                score = float(p[3])
            except ValueError:
                score = None
        rows.append({"start": min(start, end), "end": max(start, end), "score": score})
    return rows


def near(pred: dict, gt: dict, tol: int) -> bool:
    return not (pred["end"] < gt["start"] - tol or pred["start"] > gt["end"] + tol)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--gt", default=str(ROOT / "data" / "regulondb" / "TerminatorSet.tsv"))
    ap.add_argument("--tol", type=int, default=25, help="tolerancia en pb")
    ap.add_argument("--limit", type=int, help="solo curados por debajo de esta posición")
    a = ap.parse_args()

    gt = load_gt(Path(a.gt), a.limit)
    pred = load_bed(Path(a.bed))
    if not pred:
        print("sin predicciones en el BED")
        return 1

    hit = [g for g in gt if any(near(p, g, a.tol) for p in pred)]
    used = [p for p in pred if any(near(p, g, a.tol) for g in gt)]
    rho_ind = [g for g in gt if "rho-independent" in g["type"]]
    rho_dep = [g for g in gt if "rho-dependent" in g["type"]]
    hit_ind = [g for g in rho_ind if any(near(p, g, a.tol) for p in pred)]
    hit_dep = [g for g in rho_dep if any(near(p, g, a.tol) for p in pred)]

    print(f"BATTER-TPE vs RegulonDB  (tolerancia ±{a.tol} pb)")
    print(f"  curados:            {len(gt)}")
    print(f"  predicciones:       {len(pred)}")
    print(f"  sensibilidad:       {100*len(hit)/len(gt):.1f}%  ({len(hit)}/{len(gt)})")
    print(f"    rho-independiente: {100*len(hit_ind)/max(1,len(rho_ind)):.1f}% "
          f"({len(hit_ind)}/{len(rho_ind)})")
    print(f"    rho-dependiente:   {100*len(hit_dep)/max(1,len(rho_dep)):.1f}% "
          f"({len(hit_dep)}/{len(rho_dep)})")
    print(f"  precisión (suelo):  {100*len(used)/len(pred):.1f}%  "
          f"({len(used)}/{len(pred)})")
    print("  nota: el suelo subestima porque RegulonDB no es exhaustivo")
    return 0


if __name__ == "__main__":
    sys.exit(main())
