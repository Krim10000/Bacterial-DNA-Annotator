#!/usr/bin/env python3
"""benchmark_operons.py — mide un predictor de operones contra RegulonDB.

Entrada: el `operons.tsv` que produzca L2 (o cualquier predictor) y la verdad
de referencia de `regulondb_gt.py`. Todo se compara en b-numbers.

Dos métricas, y la principal NO es la de conjuntos exactos:

* **pair-level** (estándar en la literatura de predicción de operones): de
  todos los pares de genes que un operón declara co-transcritos, ¿cuántos lo
  están de verdad? Es robusta a diferencias de borde (un operón partido en dos
  pierde un par, no todo) y es la que hay que reportar.
* **exact-set**: el conjunto de genes coincide exactamente. Es dura y sirve
  para detectar si el modelo está "casi" bien o simplemente bien.

Además reporta cobertura (% de genes de RegulonDB que el predictor menciona),
genes no resueltos, violaciones (gen en dos operones) y desglose por hebra y
tamaño.

Uso
---
    python scripts/validation/benchmark_operons.py \\
        --pred data/regulation/ecoli_k12/operons.tsv --json out.json

El TSV predicho necesita, como mínimo, una columna de genes y una de hebra.
Se aceptan alias de columna (genes/gene_ids/gene_list/members, strand/sign).
Separador de genes: coma, barra o punto y coma.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GT_DEFAULT = ROOT / "data" / "regulondb" / "gt" / "gt_operons.tsv"
GENES_KEYS = ("genes", "gene_ids", "gene_list", "members", "genelist")
STRAND_KEYS = ("strand", "sign", "dir")
ID_KEYS = ("operon_id", "operon", "id", "name")


def _detect(fieldnames: list[str], keys: tuple[str, ...]) -> str | None:
    lower = {f.lower().strip(): f for f in fieldnames}
    for k in keys:
        if k in lower:
            return lower[k]
    for k in keys:
        for low, orig in lower.items():
            if k in low:
                return orig
    return None


def _split(field: str) -> list[str]:
    for sep in (",", "|", ";"):
        if sep in field:
            return [g.strip() for g in field.split(sep) if g.strip()]
    return [field.strip()] if field.strip() else []


def load_gt(path: Path, min_genes: int) -> tuple[list[set[str]], set[str]]:
    """Operones GT (filtrados) y el universo de genes GT (sin filtrar).

    El universo completo evita reportar como "no resuelto" un b-number válido
    que solo aparece en operones de un gen.
    """
    ops, all_genes = [], set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            genes = set(_split(row.get("genes", "")))
            if not genes:
                continue
            all_genes |= genes
            if row.get("mapping_ok") not in (None, "yes"):
                continue
            if len(genes) >= min_genes:
                ops.append(genes)
    return ops, all_genes


def load_gt_levels(path: Path, min_genes: int) -> dict[str, list[set[str]]]:
    """Operones GT agrupados por nivel de confianza de RegulonDB.

    `S` = evidencia fuerte (experimental); `W` = débil; `''` = sin declarar.
    Sirve para evaluar contra los casos ya documentados sin que los débiles
    diluyan el número.
    """
    out: dict[str, list[set[str]]] = {"S": [], "W": [], "": []}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("mapping_ok") not in (None, "yes"):
                continue
            genes = set(_split(row.get("genes", "")))
            if len(genes) < min_genes:
                continue
            out.setdefault(row.get("confidence", ""), []).append(genes)
    return out


def load_aliases(path: Path) -> dict[str, str]:
    """nombre de gen -> b-number, desde gt_gene_aliases.tsv."""
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {r["name"]: r["bnumber"] for r in csv.DictReader(fh, delimiter="\t")
                if r.get("name") and r.get("bnumber")}


def load_pred(path: Path, alias: dict[str, str], gt_genes: set[str]):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"{path}: sin cabecera")
        gcol = _detect(reader.fieldnames, GENES_KEYS)
        scol = _detect(reader.fieldnames, STRAND_KEYS)
        icol = _detect(reader.fieldnames, ID_KEYS)
        if not gcol:
            raise SystemExit(f"{path}: no encuentro columna de genes "
                             f"entre {reader.fieldnames}")
        ops, unresolved, seen = [], defaultdict(int), {}
        for row in reader:
            genes = set()
            for token in _split(row[gcol]):
                if token in gt_genes:
                    genes.add(token)
                elif token in alias:
                    genes.add(alias[token])
                else:
                    unresolved[token] += 1
            if not genes:
                continue
            for g in genes:
                seen[g] = seen.get(g, 0) + 1
            ops.append({"genes": genes,
                        "strand": (row.get(scol, "") if scol else ""),
                        "id": (row.get(icol, "") if icol else "")})
    violations = sum(1 for g, n in seen.items() if n > 1)
    return ops, unresolved, violations


def pairs(ops: list[set[str]]) -> set[frozenset[str]]:
    out = set()
    for op in ops:
        for a, b in combinations(sorted(op), 2):
            out.add(frozenset((a, b)))
    return out


def prf(tp: int, n_pred: int, n_gt: int) -> dict:
    p = tp / n_pred if n_pred else 0.0
    r = tp / n_gt if n_gt else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"precision": round(p, 4), "recall": round(r, 4), "f1": round(f, 4),
            "n_pred": n_pred, "n_gt": n_gt}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred", type=Path, required=True)
    ap.add_argument("--gt", type=Path, default=GT_DEFAULT)
    ap.add_argument("--min-genes", type=int, default=2,
                    help="ignora operones de menos genes (default 2)")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--by-confidence", action="store_true",
                    help="recall contra los casos documentados (S) y débiles (W)")
    args = ap.parse_args()

    if not args.gt.exists():
        raise SystemExit(f"falta {args.gt}: corre regulondb_gt.py primero")

    gt_ops, gt_genes = load_gt(args.gt, args.min_genes)
    alias = load_aliases(args.gt.parent / "gt_gene_aliases.tsv")
    pred_ops, unresolved, violations = load_pred(args.pred, alias, gt_genes)

    gt_sets = [op for op in gt_ops]
    pred_sets = [op["genes"] for op in pred_ops if len(op["genes"]) >= args.min_genes]

    gt_pairs, pred_pairs = pairs(gt_sets), pairs(pred_sets)
    pair = prf(len(gt_pairs & pred_pairs), len(pred_pairs), len(gt_pairs))

    gt_exact = {frozenset(o) for o in gt_sets}
    pred_exact = {frozenset(o) for o in pred_sets}
    exact = prf(len(gt_exact & pred_exact), len(pred_exact), len(gt_exact))

    pred_genes = {g for op in pred_ops for g in op["genes"]}
    coverage = len(pred_genes & gt_genes) / len(gt_genes) if gt_genes else 0.0

    buckets = {"2-3": (2, 3), "4-6": (4, 6), "7+": (7, 10 ** 6)}
    by_size = {}
    for label, (lo, hi) in buckets.items():
        p = [o for o in pred_sets if lo <= len(o) <= hi]
        t = [o for o in gt_sets if lo <= len(o) <= hi]
        by_size[label] = prf(len({frozenset(o) for o in p} & {frozenset(o) for o in t}),
                             len({frozenset(o) for o in p}), len({frozenset(o) for o in t}))

    result = {
        "pred": str(args.pred), "gt": str(args.gt),
        "min_genes": args.min_genes,
        "pair_level": pair, "exact_set": exact,
        "gene_coverage": round(coverage, 4),
        "pred_operons": len(pred_sets), "gt_operons": len(gt_sets),
        "unresolved_gene_tokens": len(unresolved),
        "genes_in_multiple_operons": violations,
        "by_size": by_size,
        "top_unresolved": sorted(unresolved.items(), key=lambda kv: -kv[1])[:10],
    }

    print(f"\n=== L2 vs RegulonDB ({len(gt_sets)} operones >= {args.min_genes} genes) ===")
    print(f"  pair-level  P={pair['precision']:.3f}  R={pair['recall']:.3f}  "
          f"F1={pair['f1']:.3f}   (pred {pair['n_pred']} vs gt {pair['n_gt']} pares)")
    print(f"  exact-set   P={exact['precision']:.3f}  R={exact['recall']:.3f}  "
          f"F1={exact['f1']:.3f}")
    print(f"  cobertura de genes GT: {coverage:.1%}")
    print(f"  operones predichos: {len(pred_sets)} | GT: {len(gt_sets)}")
    if violations:
        print(f"  AVISO: {violations} genes aparecen en más de un operón")
    if unresolved:
        print(f"  AVISO: {len(unresolved)} tokens de gen no resueltos "
              f"(ej. {', '.join(list(unresolved)[:3])})")
    for label, m in by_size.items():
        print(f"  tamaño {label:4s} F1={m['f1']:.3f} (n_gt={m['n_gt']})")

    if args.by_confidence:
        levels = load_gt_levels(args.gt, args.min_genes)
        print("\n  recall contra casos documentados (sin penalizar precisión):")
        for level in ("S", "W", ""):
            ops = levels.get(level) or []
            if not ops:
                continue
            g_sets = {frozenset(o) for o in ops}
            g_pairs = pairs(ops)
            found = len(g_sets & pred_exact)
            rec_pairs = len(g_pairs & pred_pairs) / len(g_pairs) if g_pairs else 0
            label = level or "(sin nivel)"
            print(f"    {label:10s} n={len(ops):3d}  exacto {found}/{len(ops)} "
                  f"({found / len(ops):.0%})  pares {rec_pairs:.3f}")
            result.setdefault("by_confidence", {})[label] = {
                "n": len(ops), "exact_recovered": found,
                "exact_rate": round(found / len(ops), 4),
                "pair_recall": round(rec_pairs, 4)}

    if args.json:
        args.json.write_text(json.dumps(result, indent=1))
        print(f"  json -> {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
