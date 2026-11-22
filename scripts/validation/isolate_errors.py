#!/usr/bin/env python3
"""isolate_errors.py — aísla los errores del motor y los diagnostica.

Cada frontera entre genes contiguos es una decisión. Comparada con el ground
truth, cada decisión cae en una de cuatro casillas:

    TP  unió y debía unir      FN  cortó y no debía cortar
    TN  cortó y debía cortar   FP  unió y no debía unir

El motor escribe qué regla decidió (`*_boundary_rules.tsv`). Este script cruza
esa decisión con las señales disponibles **en ese hueco** —TSS, horquilla,
coherencia enzimática, tamaño del hueco— y escribe una tabla con una fila por
error. Asi se puede ver si los FN tienen TSS que el motor no miró, si los FP
carecen de señal, etc., en vez de adivinar.

Salidas
-------
    data/operon_finder/errors_<etiqueta>.tsv   una fila por frontera
    consola                                    diagnóstico cruzado

Uso
---
    python scripts/validation/isolate_errors.py \
        --boundary data/operon_finder/sub_tss/subtilis_boundary_rules.tsv \
        --gbff data/dbtbs/subtilis.gbff \
        --tss data/dbtbs/gt_tss.tsv \
        --hairpins data/operon_finder/subtilis_hairpins.tsv \
        --gt operon --tag subtilis_tu
"""
from __future__ import annotations

import argparse
import bisect
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))


def load_tss(path: Path | None) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {"+": [], "-": []}
    if not path:
        return out
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            pos = int(r.get("pos") or r.get("start") or 0)
        except ValueError:
            continue
        st = (r.get("strand") or "+").strip()
        out.setdefault(st if st in ("+", "-") else "+", []).append(pos)
    for v in out.values():
        v.sort()
    return out


def load_hairpins(path: Path | None) -> list[tuple[int, str]]:
    out = []
    if not path or not Path(path).exists():
        return out
    for r in csv.DictReader(open(path), delimiter="\t"):
        try:
            out.append((int(r["start"]), r.get("class", "")))
        except (ValueError, KeyError):
            continue
    out.sort()
    return out


def between(sorted_vals: list, a: int, b: int):
    i = bisect.bisect_left(sorted_vals, a)
    return sorted_vals[i] if i < len(sorted_vals) and sorted_vals[i] <= b else None


def load_gt(path: Path, genes_col: str = "genes") -> dict[str, str]:
    if not Path(path).exists():
        return {}
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        oid = r.get("operon_id") or r.get("operon") or r.get("tu_id") or ""
        for g in (r.get(genes_col) or "").split(","):
            if g.strip():
                out[g.strip()] = oid
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", required=True)
    ap.add_argument("--gbff", default=None, help="para gap; si falta, se omite")
    ap.add_argument("--tss", default=None)
    ap.add_argument("--hairpins", default=None)
    ap.add_argument("--gt", required=True, help="TSV con columna genes y operon_id/operon/tu")
    ap.add_argument("--genes-col", default="genes")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--ec", default=None, help="TSV de genes con locus_tag y ec")
    a = ap.parse_args()

    gt = load_gt(Path(a.gt), a.genes_col)
    tss = load_tss(Path(a.tss) if a.tss else None)
    hp = load_hairpins(Path(a.hairpins) if a.hairpins else None)
    ec = {}
    if a.ec and Path(a.ec).exists():
        for r in csv.DictReader(open(a.ec), delimiter="\t"):
            if r.get("ec"):
                ec[r["locus_tag"]] = r["ec"]

    coords = {}
    if a.gbff and Path(a.gbff).exists():
        import find_operons as F
        _, genes = F.read_gbff(Path(a.gbff))
        for g in genes:
            coords[g["locus_tag"]] = (int(g["start"]), int(g["end"]),
                                      "+" if str(g.get("strand")) in ("1", "+") else "-")

    rows = []
    for r in csv.DictReader(open(a.boundary), delimiter="\t"):
        l, rr = r["left"], r["right"]
        if l not in gt or rr not in gt:
            continue
        same = gt[l] == gt[rr]
        cut = str(r.get("cut", "0")) == "1"
        tipo = ("TP" if (same and not cut) else "TN" if (not same and cut)
                else "FN" if (same and cut) else "FP")
        st = coords.get(l, (0, 0, "+"))[2]
        gap = ""
        if l in coords and rr in coords:
            gap = coords[rr][0] - coords[l][1] - 1
        t = between(tss.get(st, []), (coords[l][1] + 1) if l in coords else 0,
                    (coords[rr][0] - 1) if rr in coords else 0) if tss else None
        hpos = between([x[0] for x in hp], (coords[l][1] + 1) if l in coords else 0,
                       (coords[rr][0] - 1) if rr in coords else 0) if hp else None
        hcls = next((c for p, c in hp if p == hpos), "") if hpos else ""
        e1, e2 = ec.get(l, ""), ec.get(rr, "")
        ec_same = bool(e1 and e2 and ".".join(e1.split(".")[:2]) == ".".join(e2.split(".")[:2]))
        rows.append({"type": tipo, "left": l, "right": rr, "gap": gap,
                     "tss": int(t is not None) if tss else "",
                     "hairpin": int(hpos is not None) if hp else "",
                     "hairpin_class": hcls, "ec_left": e1, "ec_right": e2,
                     "ec_same": int(ec_same) if ec else "",
                     "gt_operon": gt[l], "rule": r.get("rule", ""),
                     "why": r.get("why", "")})

    out = ROOT / "data" / "operon_finder" / f"errors_{a.tag}.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    cnt = Counter(x["type"] for x in rows)
    print(f"fronteras evaluadas: {len(rows)}  "
          f"TP={cnt['TP']} TN={cnt['TN']} FN={cnt['FN']} FP={cnt['FP']}")
    print(f"-> {out.relative_to(ROOT)}\n")

    for tipo in ("FN", "FP"):
        sub = [x for x in rows if x["type"] == tipo]
        if not sub:
            continue
        print(f"== {tipo} ({len(sub)}) ==")
        print("  por regla: " + ", ".join(f"{k}={v}" for k, v in
                                          Counter(x["rule"] for x in sub).most_common()))
        if tss:
            con = sum(1 for x in sub if x["tss"] == 1)
            print(f"  con TSS en el hueco: {con} ({100*con/len(sub):.0f}%)")
        if hp:
            con = sum(1 for x in sub if x["hairpin"] == 1)
            print(f"  con horquilla en el hueco: {con} ({100*con/len(sub):.0f}%)")
        if ec:
            con = sum(1 for x in sub if x["ec_same"] == 1)
            print(f"  con EC compartida: {con} ({100*con/len(sub):.0f}%)")
        gaps = [x["gap"] for x in sub if isinstance(x["gap"], int)]
        if gaps:
            gaps.sort()
            print(f"  hueco: mediana {gaps[len(gaps)//2]} pb, "
                  f"p10 {gaps[len(gaps)//10]} pb, p90 {gaps[9*len(gaps)//10]} pb")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
