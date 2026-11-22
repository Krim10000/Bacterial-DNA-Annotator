#!/usr/bin/env python3
"""benchmark_subtilis.py — ¿generaliza el buscador de operones fuera de *E. coli*?

Compara operones predichos contra el ground truth de BSGatlas para *B. subtilis*
(`data/dbtbs/gt_operons.tsv`) a nivel de **pares de genes contiguos**, que es la
métrica que no depende de cómo se parta un operón largo.

También reporta el acierto de **conjunto exacto** (un operón predicho igual a uno
del GT) y estratifica por tamaño y por promotores internos.

Uso
---
    python scripts/validation/benchmark_subtilis.py \
        --pred otsu:data/operon_finder/subtilis_otsu/subtilis_operons_rules.tsv \
        --pred fijo:data/operon_finder/subtilis_fixed/subtilis_operons_rules.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GT = ROOT / "data" / "dbtbs" / "gt_operons.tsv"
OUT = ROOT / "data" / "dbtbs" / "benchmark_subtilis.json"


def read_gt(path: Path) -> list[dict]:
    ops = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
        if not genes:
            continue
        ops.append({
            "id": r.get("operon_id", ""), "genes": genes,
            "n": len(genes),
            "n_tus": int(r.get("n_genes_tus") or 0),
            "strand": r.get("strand", ""),
            "start": int(r.get("start") or 0), "end": int(r.get("end") or 0),
        })
    return ops


def add_internal_promoters(ops: list[dict], tss_path: Path) -> None:
    """Marca qué operones tienen >1 TSS dentro: esos tienen promotor interno."""
    if not tss_path.exists():
        return
    tss = [r for r in csv.DictReader(open(tss_path), delimiter="\t")]
    for o in ops:
        n = 0
        for t in tss:
            try:
                pos = int(t["pos"])
            except (KeyError, ValueError):
                continue
            if t.get("strand") != o["strand"]:
                continue
            if o["start"] - 200 <= pos <= o["end"] + 200:
                n += 1
        o["n_tss"] = n
        o["internal_promoter"] = n > 1


def read_pred(path: Path) -> list[dict]:
    ops = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        genes = [g.strip() for g in (r.get("genes") or "").split(",") if g.strip()]
        if genes:
            ops.append({"id": r.get("operon_id", ""), "genes": genes, "n": len(genes)})
    return ops


def pairs_of(ops: list[dict]) -> set[tuple[str, str]]:
    """Pares consecutivos, normalizados a orden alfabético por gen."""
    out = set()
    for op in ops:
        g = op["genes"]
        for a, b in zip(g, g[1:]):
            out.add((a, b) if a < b else (b, a))
    return out


def score(pred: list[dict], gt: list[dict]) -> dict:
    pp, gp = pairs_of(pred), pairs_of(gt)
    tp = len(pp & gp)
    p = tp / len(pp) if pp else 0.0
    r = tp / len(gp) if gp else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    gsets = {frozenset(o["genes"]) for o in gt}
    exact = sum(1 for o in pred if frozenset(o["genes"]) in gsets)
    return {"pred_pairs": len(pp), "gt_pairs": len(gp), "tp": tp,
            "P": round(p, 3), "R": round(r, 3), "F1": round(f, 3),
            "exact": exact, "exact_rate": round(exact / len(pred), 3) if pred else 0,
            "pred_ops": len(pred)}


def stratified(pred: list[dict], gt: list[dict]) -> list[dict]:
    """F1 por tamaño de operón del GT y por promotor interno."""
    rows = []
    pp = pairs_of(pred)
    for label, cond in (
        ("GT 2 genes", lambda o: o["n"] == 2),
        ("GT 3 genes", lambda o: o["n"] == 3),
        ("GT 4-5 genes", lambda o: 4 <= o["n"] <= 5),
        ("GT 6+ genes", lambda o: o["n"] >= 6),
        ("un solo TSS (sin promotor interno)", lambda o: not o.get("internal_promoter")),
        ("2+ TSS (con promotor interno)", lambda o: o.get("internal_promoter")),
    ):
        sub = [o for o in gt if cond(o)]
        gp = pairs_of(sub)
        tp = len(pp & gp)
        # precisión del estrato: solo los pares predichos cuyos dos genes
        # pertenecen al universo de este grupo (si no, se subestima)
        universe = {g for o in sub for g in o["genes"]}
        pred_in = {pr for pr in pp if pr[0] in universe and pr[1] in universe}
        p = tp / len(pred_in) if pred_in else 0.0
        r = tp / len(gp) if gp else 0.0
        f = 2 * p * r / (p + r) if (p + r) else 0.0
        rows.append({"grupo": label, "operones_gt": len(sub), "pares_gt": len(gp),
                     "P": round(p, 3), "R": round(r, 3), "F1": round(f, 3)})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred", action="append", default=[],
                    help="etiqueta:ruta (repetible)")
    ap.add_argument("--gt", type=Path, default=GT)
    ap.add_argument("--gbff", type=Path, default=ROOT / "data" / "dbtbs" / "subtilis.gbff",
                    help="GBFF para restringir a genes anotados como CDS")
    a = ap.parse_args()

    gt = read_gt(a.gt)
    add_internal_promoters(gt, a.gt.parent / "gt_tss.tsv")
    gt_pairs = pairs_of(gt)
    multi = [o for o in gt if o["n"] > 1]
    print(f"GT: {len(gt)} operones · {len(multi)} multi-génicos · "
          f"{len(gt_pairs)} pares · genes: {len({g for o in gt for g in o['genes']})}")

    if not a.pred:
        print("sin predicciones: pasa --pred etiqueta:ruta")
        return 1

    # techo de recall: genes del GT que no son CDS en la anotación
    ann = set()
    if a.gbff.exists():
        import re as _re
        for line in open(a.gbff):
            m = _re.search(r'/locus_tag="([^"]+)"', line)
            if m:
                ann.add(m.group(1))
    gt_all = gt
    if ann:
        gt = [dict(o, genes=[g for g in o["genes"] if g in ann]) for o in gt]
        gt_pairs_eval = pairs_of([o for o in gt if o["n"] > 1])
        print(f"pares GT evaluables (ambos genes CDS): {len(gt_pairs_eval)} "
              f"de {len(gt_pairs)} ({100*len(gt_pairs_eval)/max(1,len(gt_pairs)):.1f}%)")

    gt_sets = {frozenset(x["genes"]) for x in gt_all}
    results = {}
    print(f"\n{'configuración':28s} {'ops':>5s} {'pares':>6s} {'P':>6s} {'R':>6s} "
          f"{'F1':>6s} {'exacto':>7s}")
    for item in a.pred:
        label, _, path = item.partition(":")
        pred = read_pred(Path(path))
        pred_full = pred
        if ann:  # solo pares con ambos genes anotados (comparación justa)
            pred = [dict(o, genes=[g for g in o["genes"] if g in ann]) for o in pred]
            pred = [o for o in pred if len(o["genes"]) > 1]
        sc = score(pred, gt)
        sc["exact_full"] = sum(1 for o in pred_full
                               if frozenset(o["genes"]) in gt_sets)
        sc["exact_full_rate"] = round(sc["exact_full"] / max(1, len(pred_full)), 3)
        results[label] = sc
        print(f"{label:28s} {sc['pred_ops']:5d} {sc['pred_pairs']:6d} "
              f"{sc['P']:6.3f} {sc['R']:6.3f} {sc['F1']:6.3f} {sc['exact_full']:7d}")

    # estratificado: sobre la configuración con mejor F1 (la más informativa)
    label0 = max(results, key=lambda k: results[k]["F1"])
    pred0 = read_pred(Path(dict(x.partition(":")[::2] for x in a.pred)[label0]))
    print(f"\nestratificado ({label0}):")
    print(f"{'grupo':34s} {'ops':>5s} {'pares':>6s} {'P':>6s} {'R':>6s} {'F1':>6s}")
    strat = stratified(pred0, gt)
    for r in strat:
        print(f"{r['grupo']:34s} {r['operones_gt']:5d} {r['pares_gt']:6d} "
              f"{r['P']:6.3f} {r['R']:6.3f} {r['F1']:6.3f}")

    # fp/fn más frecuentes por regla (si existe el detalle de fronteras)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"gt": {"operons": len(gt), "multi": len(multi),
                                      "pairs": len(gt_pairs)},
                               "configs": results, "stratified": strat}, indent=1))
    print(f"\n-> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
