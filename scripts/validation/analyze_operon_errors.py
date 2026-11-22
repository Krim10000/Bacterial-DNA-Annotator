#!/usr/bin/env python3
"""analyze_operon_errors.py — ¿qué reglas se equivocan, y en qué casos?

En vez de contar pares, se analizan las **decisiones de frontera**: cada par de
genes contiguos recibe una regla y una decisión (cortar / no cortar). Comparando
esa decisión con el ground truth de RegulonDB se obtiene el error atribuido a la
regla que lo causó:

- FN (falso negativo): la regla **cortó** donde el GT dice que siguen juntos.
- FP (falso positivo): la regla **no cortó** donde el GT dice que se separan.

Salida: tabla por regla, tabla por tamaño de operón y ejemplos.

Uso
---
    python scripts/validation/analyze_operon_errors.py \
        --boundary data/operon_finder/k12_p16b/ecoli_k12_boundary_rules.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))
import find_operons as F  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundary", required=True)
    ap.add_argument("--examples", type=int, default=4)
    ap.add_argument("--md", default=str(ROOT / "docs" / "OPERON_ERRORS.md"))
    a = ap.parse_args()

    op_of = F.load_operon_of()
    size_of = Counter(op_of.values())

    rows = list(csv.DictReader(open(a.boundary), delimiter="\t"))
    ok_merge = ok_cut = fp = fn = unknown = 0
    fp_by_rule: Counter = Counter()
    fn_by_rule: Counter = Counter()
    fn_by_size: Counter = Counter()
    fp_by_size: Counter = Counter()
    fn_total_by_size: Counter = Counter()
    examples_fp: defaultdict = defaultdict(list)
    examples_fn: defaultdict = defaultdict(list)

    for r in rows:
        l, rr = r["left"], r["right"]
        if l not in op_of or rr not in op_of:
            unknown += 1
            continue
        same = op_of[l] == op_of[rr]
        cut = str(r.get("cut", "0")) == "1"
        rule = r.get("rule", "?")
        why = r.get("why", "")
        gt_size = size_of.get(op_of[l], 0) if same else 0
        if same and not cut:
            ok_merge += 1
            fn_total_by_size[gt_size] += 1
        elif same and cut:
            fn += 1
            fn_by_rule[rule] += 1
            fn_by_size[gt_size] += 1
            fn_total_by_size[gt_size] += 0  # ya contado arriba
            if len(examples_fn[rule]) < a.examples:
                examples_fn[rule].append((l, rr, why))
        elif not same and cut:
            ok_cut += 1
        else:
            fp += 1
            fp_by_rule[rule] += 1
            fp_by_size[1] += 1
            if len(examples_fp[rule]) < a.examples:
                examples_fp[rule].append((l, rr, why))

    total = ok_merge + ok_cut + fp + fn
    print(f"decisiones evaluadas: {total} (sin GT: {unknown})")
    print(f"  aciertos al unir: {ok_merge}")
    print(f"  aciertos al cortar: {ok_cut}")
    print(f"  FN (cortó de más): {fn}")
    print(f"  FP (unió de más): {fp}")

    lines = ["# Dónde se equivoca el motor de operones", "",
             f"Analiza {total} decisiones de frontera contra RegulonDB.", "",
             "| decisión | n |", "|---|---:|",
             f"| acierto al unir | {ok_merge} |",
             f"| acierto al cortar | {ok_cut} |",
             f"| **FN: cortó de más** | **{fn}** |",
             f"| **FP: unió de más** | **{fp}** |", "",
             "## FN por regla (cortes que no debían hacerse)", "",
             "| regla | n | % de los FN |", "|---|---:|---:|"]
    for rule, n in fn_by_rule.most_common():
        lines.append(f"| `{rule}` | {n} | {100*n/max(1,fn):.1f}% |")
    lines += ["", "## FP por regla (fusiones que no debían hacerse)", "",
              "| regla | n | % de los FP |", "|---|---:|---:|"]
    for rule, n in fp_by_rule.most_common():
        lines.append(f"| `{rule}` | {n} | {100*n/max(1,fp):.1f}% |")

    lines += ["", "## FN según tamaño del operón real", "",
              "| genes en el operón | FN | aciertos | % de error |", "|---|---:|---:|---:|"]
    for size in sorted(fn_by_size):
        ac = fn_total_by_size.get(size, 0)
        lines.append(f"| {size} | {fn_by_size[size]} | {ac} | "
                     f"{100*fn_by_size[size]/max(1, ac):.1f}% |")

    for title, ex in (("Ejemplos de FN", examples_fn), ("Ejemplos de FP", examples_fp)):
        lines += ["", f"## {title}", ""]
        for rule, items in ex.items():
            for l, rr, why in items:
                lines.append(f"- `{rule}` · {l} → {rr}: {why}")

    Path(a.md).write_text("\n".join(lines) + "\n")
    print(f"-> {Path(a.md).relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
