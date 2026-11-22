#!/usr/bin/env python3
"""Capa 1 (paso 2) — asignar cada matriz de Reg-Seq a su TF por motivo.

El mapeo por promotor falla mucho (54 % de los promotores de Reg-Seq no están en la red
de RegulonDB: eran justamente los poco caracterizados). Alternativa independiente:
comparar el consenso de la matriz con la biblioteca de PWMs de RegulonDB
(`data/regulondb/gt/tf_pwms.tsv`) y asignar el TF de mejor puntuación normalizada.

Salida: data/allostery/regseq_motif_map.tsv + resumen en JSON.

Uso: python scripts/allostery/12_regseq_motif_map.py
"""
from __future__ import annotations

import ast
import collections
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery"
PWM = ROOT / "data" / "regulondb" / "gt" / "tf_pwms.tsv"
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def load_pwms():
    out = []
    for r in csv.DictReader(open(PWM), delimiter="\t"):
        try:
            m = ast.literal_eval(r["pwm"])
        except Exception:
            continue
        if not m:
            continue
        cons = r["consensus"]
        self_score = sum(m[i][IDX[cons[i]]] for i in range(min(len(cons), len(m))))
        out.append({"tf": r["tf"], "pwm": m, "width": len(m),
                    "cons": cons, "self": self_score})
    return out


def best_match(seq: str, pwms: list[dict]):
    best, second = (-1e9, None), (-1e9, None)
    for p in pwms:
        w, m = p["width"], p["pwm"]
        if w > len(seq):
            continue
        for s in range(len(seq) - w + 1):
            score = 0.0
            for i, ch in enumerate(seq[s:s + w]):
                score += m[i][IDX.get(ch, 0)]
            norm = score / p["self"] if p["self"] > 0 else 0.0
            if norm > best[0]:
                second = best
                best = (norm, p["tf"])
            elif norm > second[0]:
                second = (norm, p["tf"])
    return best[0], best[1], second[0]


def main():
    pwms = load_pwms()
    rows = list(csv.DictReader(open(OUT / "regseq_matrices.tsv"), delimiter="\t"))
    out = []
    for r in rows:
        cons = r["consenso"]
        score, tf, second = best_match(cons, pwms)
        margin = score - second if second > -1e8 else None
        rd = r["tfs_candidatos"].split("|") if r["tfs_candidatos"] else []
        out.append({**r, "tf_motivo": tf, "score_norm": round(score, 3),
                    "margen": None if margin is None else round(margin, 3),
                    "coincide_regulondb": (tf in rd) if rd and tf else None})

    with open(OUT / "regseq_motif_map.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    conf = [o for o in out if o["score_norm"] >= 0.9 and (o["margen"] or 0) >= 0.03]
    rd_ok = [o for o in out if o["coincide_regulondb"] is True]
    rd_bad = [o for o in out if o["coincide_regulondb"] is False]
    summary = {
        "n_matrices": len(out),
        "con_match_convincente": len(conf),
        "fraccion": round(len(conf) / len(out), 3) if out else None,
        "tfs_distintos_asignados": len({o["tf_motivo"] for o in conf}),
        "acuerdo_con_regulondb": f"{len(rd_ok)}/{len(rd_ok) + len(rd_bad)}",
        "top_tfs": collections.Counter(o["tf_motivo"] for o in conf).most_common(12),
    }
    (OUT / "regseq_motif_map.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
