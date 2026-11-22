#!/usr/bin/env python3
"""ensemble.py — ¿suman dos métodos independientes o solo se repiten?

L2 (`data/regulation/ecoli_k12/operons.tsv`) y `find_operons`
(`data/operon_finder/ecoli_k12_operons.tsv`) llegan a ~0.78–0.81 por caminos
distintos. Si sus errores son distintos, combinarlos supera a ambos; si son los
mismos, solo se confirman.

Se compara a nivel de **par de genes co-transcritos**, que es la métrica del
benchmark y la única comparable entre métodos:

    A  = pares que L2 declara juntos
    B  = pares que find_operons declara juntos
    AND = A ∩ B   (conservador: solo lo que ambos afirman)
    OR  = A ∪ B   (agresivo: cualquiera de los dos)

Además reporta el desglose de aciertos/fallos de cada uno por separado, que es
lo que dice si hay complementariedad real.

Uso
---
    python scripts/regulation/ensemble.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GT = ROOT / "data" / "regulondb" / "gt" / "gt_operons.tsv"
L2 = ROOT / "data" / "regulation" / "ecoli_k12" / "operons.tsv"
FINDER = ROOT / "data" / "operon_finder" / "ecoli_k12_operons.tsv"
OUT = ROOT / "data" / "operon_finder" / "ensemble_ecoli_k12.tsv"


def split(field: str) -> list[str]:
    for sep in (",", "|", ";"):
        if sep in field:
            return [x.strip() for x in field.split(sep) if x.strip()]
    return [field.strip()] if field.strip() else []


def pairs_of(path: Path, key: str = "genes") -> set[frozenset]:
    out: set[frozenset] = set()
    if not path.exists():
        return out
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes = split(row.get(key, ""))
            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    if genes[i] != genes[j]:
                        out.add(frozenset((genes[i], genes[j])))
    return out


def gt_pairs() -> set[frozenset]:
    out: set[frozenset] = set()
    with open(GT, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("mapping_ok") not in (None, "yes"):
                continue
            if int(row.get("n_genes") or 0) < 2:
                continue
            genes = split(row["genes"])
            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    if genes[i] != genes[j]:
                        out.add(frozenset((genes[i], genes[j])))
    return out


def prf(pred: set, gt: set) -> tuple[float, float, float]:
    tp = len(pred & gt)
    p = tp / len(pred) if pred else 0.0
    r = tp / len(gt) if gt else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f


def main() -> int:
    gt = gt_pairs()
    a = pairs_of(L2)
    b = pairs_of(FINDER)
    if not a or not b:
        print("faltan predicciones: ejecuta L2 y find_operons primero", file=sys.stderr)
        return 1

    rows = []
    for name, pred in (("L2", a), ("find_operons", b),
                       ("AND (intersección)", a & b), ("OR (unión)", a | b)):
        p, r, f = prf(pred, gt)
        rows.append((name, len(pred), p, r, f))
        print(f"{name:20s} pares={len(pred):5d}  P={p:.3f}  R={r:.3f}  F1={f:.3f}")

    print("\ncomplementariedad (sobre los pares verdaderos):")
    tp_a, tp_b = a & gt, b & gt
    print(f"  aciertan ambos: {len(tp_a & tp_b):4d}")
    print(f"  solo L2:        {len(tp_a - tp_b):4d}")
    print(f"  solo finder:    {len(tp_b - tp_a):4d}")
    print(f"  ninguno:        {len(gt - (tp_a | tp_b)):4d}")

    print("\nerrores de cada uno (sobre sus predicciones):")
    fp_a, fp_b = a - gt, b - gt
    print(f"  falsos positivos L2:     {len(fp_a):4d} (de {len(a)})")
    print(f"  falsos positivos finder: {len(fp_b):4d} (de {len(b)})")
    print(f"  falsos positivos comunes: {len(fp_a & fp_b):4d}")

    best = max(rows, key=lambda x: x[4])
    print(f"\nmejor: {best[0]} con F1 {best[4]:.3f}")

    # guarda la unión como propuesta (mayor cobertura), con su F1 real
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pair_a", "pair_b", "en_L2", "en_finder"])
        for pair in sorted(a | b, key=lambda s: sorted(s)):
            g = sorted(pair)
            if len(g) < 2:
                continue
            w.writerow([g[0], g[1], "yes" if pair in a else "no",
                        "yes" if pair in b else "no"])
    print(f"pares combinados -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
