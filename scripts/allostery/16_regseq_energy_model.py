#!/usr/bin/env python3
"""Capa 1 (paso 3) — modelo de energía y validación sobre los pares matriz→TF.

Para cada matriz de Reg-Seq ya mapeada a un TF (174 pares):
  1. Se parsea la matriz (ΔΔG por base y posición).
  2. Se toma la secuencia wild-type de la ventana desde `prior_designs/wtsequences.csv`.
  3. Se puntúa: ΔG(sitio) = Σ_i matriz[i][base_i]   (más bajo = mejor unión).
  4. Nulo: 500 permutaciones de la misma ventana (misma composición de bases).
  5. z = (ΔG_wt − media_azar) / sd_azar. Si la matriz es informativa, z << 0.

Se prueban las dos orientaciones de la secuencia (directa y reverso-complementaria) y se
reporta cuál funciona mejor: es un diagnóstico del convenio de hebra.

Salida: data/allostery/regseq_energy_model.tsv + resumen JSON.
"""
from __future__ import annotations

import collections
import csv
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery"
MAT = OUT / "work" / "regseq" / "matrices"
WT = "https://raw.githubusercontent.com/RPGroup-PBoC/RegSeq/master/data/prior_designs/wtsequences.csv"
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def parse_matrix(path: Path):
    rows = [l.split() for l in path.read_text(errors="ignore").splitlines()[1:] if l.strip()]
    pos, vals = [], []
    for f in rows:
        if len(f) < 5:
            continue
        try:
            pos.append(int(f[0]))
            vals.append([float(x) for x in f[1:5]])
        except ValueError:
            continue
    return pos, np.array(vals) if vals else None


def score(seq: str, vals: np.ndarray) -> float:
    if len(seq) != len(vals):
        return float("nan")
    return float(sum(vals[i][IDX[ch]] for i, ch in enumerate(seq) if ch in IDX))


def load_wt() -> dict[str, dict]:
    data = urlopen(WT, timeout=60).read().decode()
    out = {}
    for r in csv.DictReader(data.splitlines()):
        out[r["name"].lower()] = {"seq": (r.get("geneseq") or "").upper(),
                                  "rev": r["rev"].strip().lower() == "rev"}
    return out


def main():
    wt = load_wt()
    rows = list(csv.DictReader(open(OUT / "regseq_tf_map.tsv"), delimiter="\t"))
    rows = [r for r in rows if r["n_tf"] == "1"]
    rng = np.random.default_rng(0)

    out = []
    for r in rows:
        p = r["promotor"].lower()
        key = p if p in wt else next((k for k in wt if k in p or p in k), None)
        if not key:
            continue
        seq = wt[key]["seq"]
        pos, vals = parse_matrix(MAT / r["archivo"])
        if vals is None:
            continue
        lo, hi = int(r["from"]), int(r["to"])
        # las posiciones de la matriz son índices del amplicón
        idx = [i for i, pp in enumerate(pos) if lo <= pp <= hi]
        if len(idx) < 6 or hi >= len(seq):
            continue
        sub_vals = vals[idx]
        start = pos[idx[0]]
        window = seq[start:start + len(idx)]
        if len(window) != len(idx):
            continue
        for tag, w in (("fwd", window), ("rc", rc(window))):
            s0 = score(w, sub_vals)
            if np.isnan(s0):
                continue
            perm = []
            for _ in range(500):
                sh = list(w)
                rng.shuffle(sh)
                perm.append(score("".join(sh), sub_vals))
            perm = np.array(perm)
            z = float((s0 - perm.mean()) / perm.std()) if perm.std() > 0 else 0.0
            out.append({"archivo": r["archivo"], "tf": r["tf_suplemento"],
                        "promotor": r["promotor"], "hebra": tag, "n_pos": len(idx),
                        "dg_wt": round(s0, 4), "z": round(z, 2)})

    with open(OUT / "regseq_energy_model.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    summary = {}
    for tag in ("fwd", "rc"):
        sub = [x for x in out if x["hebra"] == tag]
        if not sub:
            continue
        zs = np.array([x["z"] for x in sub])
        summary[tag] = {"n": len(sub),
                        "z_mediana": round(float(np.median(zs)), 2),
                        "fraccion_sig_negativa": round(float((zs < -1.645).mean()), 3),
                        "fraccion_sig_positiva": round(float((zs > 1.645).mean()), 3)}
    best = max(summary, key=lambda k: summary[k]["fraccion_sig_negativa"]) if summary else None
    summary["mejor_hebra"] = best
    (OUT / "regseq_energy_model.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
