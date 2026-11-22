#!/usr/bin/env python3
"""Capa 1 (paso 2b) — asignar TF a cada matriz de Reg-Seq por núcleo del PWM.

Corrige dos fallos del intento anterior:
  1. Los PWM de RegulonDB miden 30-40 pb; las ventanas de matriz, 8-50. Se recorta cada
     PWM a su **núcleo de alta información** (ventana contigua de 12-16 pb que maximiza
     la suma de rangos por posición).
  2. El motivo puede estar en cualquiera de las dos hebras: se prueban ambas.

Puntuación normalizada: score(ventana) / score(consenso del núcleo), en [0, 1].

Validación: para las matrices cuyo promotor tiene **un solo** TF en RegulonDB, ¿coincide
la asignación por motivo?

Salida: data/allostery/regseq_motif_core.tsv + resumen JSON.
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
COMP = str.maketrans("ACGT", "TGCA")
CORE_W = 14


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def load_cores():
    out = []
    for r in csv.DictReader(open(PWM), delimiter="\t"):
        try:
            m = ast.literal_eval(r["pwm"])
        except Exception:
            continue
        if not m or len(m) < CORE_W:
            continue
        info = [max(row) - min(row) for row in m]
        best, best_i = -1e9, 0
        for i in range(len(m) - CORE_W + 1):
            s = sum(info[i:i + CORE_W])
            if s > best:
                best, best_i = s, i
        core = m[best_i:best_i + CORE_W]
        cons = "".join("ACGT"[max(range(4), key=lambda k: row[k])] for row in core)
        self_f = sum(core[i][IDX[cons[i]]] for i in range(CORE_W))
        rc = revcomp(cons)
        self_r = sum(core[i][IDX[rc[i]]] for i in range(CORE_W))
        out.append({"tf": r["tf"], "core": core, "cons": cons,
                    "cons_rc": rc, "self_f": self_f, "self_r": self_r,
                    "info": round(best, 2)})
    return out


def score(seq: str, p: dict):
    """Mejor score normalizado sobre ambas hebras (None si la ventana es corta)."""
    best = None
    for cons, self_s in ((p["cons"], p["self_f"]), (p["cons_rc"], p["self_r"])):
        if self_s <= 0:
            continue
        w = CORE_W
        if w > len(seq):
            continue
        for s in range(len(seq) - w + 1):
            sc = sum(p["core"][i][IDX.get(ch, 0)] for i, ch in enumerate(seq[s:s + w]))
            v = sc / self_s
            if best is None or v > best:
                best = v
    return best


def main():
    cores = load_cores()
    rows = list(csv.DictReader(open(OUT / "regseq_matrices.tsv"), delimiter="\t"))
    res = []
    for r in rows:
        seq = r["consenso"]
        scored = []
        for p in cores:
            v = score(seq, p)
            if v is not None:
                scored.append((v, p["tf"]))
        if not scored:
            res.append({**r, "tf_core": None, "score": None, "margen": None,
                        "coincide": None})
            continue
        scored.sort(reverse=True)
        best_v, best_tf = scored[0]
        second = scored[1][0] if len(scored) > 1 else -1
        rd = r["tfs_candidatos"].split("|") if r["tfs_candidatos"] else []
        res.append({**r, "tf_core": best_tf, "score": round(best_v, 3),
                    "margen": round(best_v - second, 3),
                    "coincide": (best_tf in rd) if rd else None})

    with open(OUT / "regseq_motif_core.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(res[0].keys()))
        w.writeheader()
        w.writerows(res)

    conf = [x for x in res if x["score"] is not None and x["score"] >= 0.85
            and (x["margen"] or 0) >= 0.05]
    val = [x for x in res if x["coincide"] is not None]
    agree = sum(1 for x in val if x["coincide"])
    summary = {
        "n_matrices": len(res),
        "con_score": sum(1 for x in res if x["score"] is not None),
        "confiables": len(conf),
        "fraccion_confiables": round(len(conf) / len(res), 3) if res else None,
        "tfs_asignados": len({x["tf_core"] for x in conf}),
        "validacion_regulondb": f"{agree}/{len(val)}",
        "top_tfs": collections.Counter(x["tf_core"] for x in conf).most_common(12),
    }
    (OUT / "regseq_motif_core.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
