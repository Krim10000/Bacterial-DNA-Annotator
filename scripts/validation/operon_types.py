#!/usr/bin/env python3
"""operon_types.py — ¿cuántos TIPOS de operón hay según sus componentes?

Tratar los operones como un caso único es lo que hace que un umbral global
funcione mal. Este script agrupa los 2.605 operones curados de K-12 por su
**firma de componentes** en tipos interpretables, y perfila cada uno.

Tipos (reglas por prioridad, no clustering opaco)
-------------------------------------------------
  monogenico_sin_regulacion  1 gen, sin componentes anotados
  monogenico_con_promotor    1 gen, con promotor/sigma
  multigenico_sin_evidencia  >=2 genes, sin componentes
  multigenico_interno        >=2 genes con más de una unidad de transcripción
                             (promotores internos: el caso que rompe la
                             delimitación simple)
  multigenico_sigma_alt      >=2 genes gobernados por sigma != sigma70
  multigenico_simple         >=2 genes, una sola TU, sigma70

Salida: `data/regulondb/gt/gt_operon_types.tsv` y un perfil por tipo.

Uso
---
    python scripts/validation/operon_types.py
"""
from __future__ import annotations

import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GT = ROOT / "data" / "regulondb" / "gt"
COMP = GT / "gt_operon_components.tsv"
OPS = GT / "gt_operons.tsv"
TUS = GT / "gt_tus.tsv"
OUT = GT / "gt_operon_types.tsv"


def load(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path), delimiter="\t"))


def main() -> int:
    comp = {r["operon_id"]: r for r in load(COMP)}
    tus: dict[str, list[dict]] = defaultdict(list)
    for t in load(TUS):
        tus[t["operon_id"]].append(t)

    rows = []
    for op in load(OPS):
        c = comp.get(op["operon_id"])
        if not c:
            continue
        n_genes = int(op.get("n_genes") or 0)
        n_tus = len(tus.get(op["operon_id"], []))
        sigmas = {t.get("sigma_factor", "") for t in tus.get(op["operon_id"], [])}
        sigmas.discard("")
        sigma = "sigma70" if "sigma70" in sigmas else (sorted(sigmas)[0] if sigmas else "")
        has_comp = int(c["n_components"]) > 0

        if n_genes == 1:
            tipo = "monogenico_con_promotor" if has_comp else "monogenico_sin_regulacion"
        elif not has_comp:
            tipo = "multigenico_sin_evidencia"
        elif n_tus >= 2:
            tipo = "multigenico_interno"
        elif sigma and sigma != "sigma70":
            tipo = "multigenico_sigma_alt"
        else:
            tipo = "multigenico_simple"

        rows.append({
            "operon_id": op["operon_id"], "operon_name": op.get("operon_name", ""),
            "type": tipo, "n_genes": n_genes, "n_tus": n_tus, "sigma": sigma,
            "has_terminator": c["has_terminator"], "has_tf": c["has_tf"],
            "n_components": c["n_components"], "n_pmids": c["n_pmids"],
            "confidence": op.get("confidence", ""),
        })

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    by_type: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_type[r["type"]].append(r)

    print(f"operones clasificados: {len(rows)} -> {OUT.relative_to(ROOT)}")
    print(f"\n{'tipo':28s} {'n':>5s} {'%':>5s} {'genes med':>9s} "
          f"{'con term':>8s} {'con TF':>7s}")
    for t, sub in sorted(by_type.items(), key=lambda kv: -len(kv[1])):
        genes = statistics.median(r["n_genes"] for r in sub)
        term = sum(1 for r in sub if r["has_terminator"] == "1") / len(sub)
        tf = sum(1 for r in sub if r["has_tf"] == "1") / len(sub)
        print(f"{t:28s} {len(sub):5d} {100*len(sub)/len(rows):4.0f}% "
              f"{genes:9.1f} {100*term:7.0f}% {100*tf:6.0f}%")

    multi = [r for r in rows if r["n_genes"] >= 2]
    print(f"\nmultigénicos: {len(multi)}")
    print("  sigma:", Counter(r["sigma"] for r in multi).most_common(6))
    print("  con >1 TU (promotores internos):",
          sum(1 for r in multi if r["n_tus"] >= 2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
