#!/usr/bin/env python3
"""F1a — baseline de familia sobre los represores decisivos.

Pregunta: ¿la familia del TF (dominio Pfam de unión a ADN) por sí sola predice
si su efector es inductor o correpresor?

Dos evaluaciones, porque la puerta G1 del plan estaba mal definida:

1. LOFO (leave-one-family-out) con SOLO la familia es **indefinido**: al sacar una
   familia no queda ningún ejemplo de entrenamiento de esa familia. Se reporta
   cuántas familias son singleton para que quede documentado.
2. LOO dentro de familia (predice la clase mayoritaria de la familia sin el propio
   TF). Es la medida honesta de "¿cuánto sabe la familia?".

Salida: data/allostery/baseline_family.json + baseline_family.md
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import WORK, ROOT, load_decisive, load_families  # noqa: E402

OUT = ROOT / "data" / "allostery"


def metrics(y, pred):
    tp = sum(1 for a, b in zip(y, pred) if a == "correpresor" and b == "correpresor")
    tn = sum(1 for a, b in zip(y, pred) if a == "inductor" and b == "inductor")
    fp = sum(1 for a, b in zip(y, pred) if a == "inductor" and b == "correpresor")
    fn = sum(1 for a, b in zip(y, pred) if a == "correpresor" and b == "inductor")
    rec_c = tp / (tp + fn) if tp + fn else 0.0
    rec_i = tn / (tn + fp) if tn + fp else 0.0
    bal = (rec_c + rec_i) / 2
    den = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    mcc = (tp * tn - fp * fn) / den if den else 0.0
    return {"balanced_acc": bal, "mcc": mcc, "recall_correpresor": rec_c,
            "recall_inductor": rec_i,
            "confusion": {"tp_cor": tp, "tn_ind": tn, "fp_ind_as_cor": fp, "fn_cor_as_ind": fn}}


def main():
    data = [d for d in load_decisive() if d["seq"]]
    fams = load_families()
    for d in data:
        d["family"] = fams.get(d["tf"], "sin_dominio")

    by_fam = collections.defaultdict(list)
    for d in data:
        by_fam[d["family"]].append(d)

    rows = []
    for f, members in sorted(by_fam.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        c = collections.Counter(m["mode"] for m in members)
        rows.append({"family": f, "n": len(members),
                     "inductor": c["inductor"], "correpresor": c["correpresor"],
                     "majority": c.most_common(1)[0][0] if members else None})

    # LOO dentro de familia
    y, pred = [], []
    for d in data:
        others = [m for m in by_fam[d["family"]] if m["tf"] != d["tf"]]
        if others:
            c = collections.Counter(m["mode"] for m in others)
            guess = c.most_common(1)[0][0]
        else:
            guess = "inductor"  # mayoría global
        y.append(d["mode"])
        pred.append(guess)
    loo = metrics(y, pred)

    # mayoría global
    glob = collections.Counter(d["mode"] for d in data).most_common(1)[0][0]
    base = metrics([d["mode"] for d in data], [glob] * len(data))

    singletons = sum(1 for r in rows if r["n"] == 1)
    summary = {
        "n": len(data),
        "n_familias": len(rows),
        "familias_singleton": singletons,
        "fraccion_singleton": singletons / len(rows) if rows else None,
        "por_familia": rows,
        "loo_dentro_de_familia": loo,
        "mayoria_global": {"clase": glob, **base},
        "lofo_con_solo_familia": "indefinido (familia singleton no tiene train)",
    }
    (OUT / "baseline_family.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    L = ["# F1a — baseline de familia", "",
         f"n={len(data)} · familias={len(rows)} · singleton={singletons} "
         f"({100 * singletons / max(1, len(rows)):.0f}%)", "",
         "| familia | n | ind | cor | mayoritaria |", "|---|---:|---:|---:|---|"]
    for r in rows:
        L.append(f"| {r['family']} | {r['n']} | {r['inductor']} | {r['correpresor']} | {r['majority']} |")
    L += ["", "## Resultado", "",
          "| modelo | balanced acc | MCC | recall cor | recall ind |", "|---|---:|---:|---:|---:|",
          f"| LOO dentro de familia | {loo['balanced_acc']:.3f} | {loo['mcc']:.3f} | "
          f"{loo['recall_correpresor']:.3f} | {loo['recall_inductor']:.3f} |",
          f"| mayoría global ({glob}) | {base['balanced_acc']:.3f} | {base['mcc']:.3f} | "
          f"{base['recall_correpresor']:.3f} | {base['recall_inductor']:.3f} |",
          "", "**LOFO con solo la familia es indefinido**: "
          f"{singletons}/{len(rows)} familias tienen un único TF decisivo, así que al "
          "sacar la familia no queda entrenamiento.", ""]
    (OUT / "baseline_family.md").write_text("\n".join(L) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "por_familia"},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
