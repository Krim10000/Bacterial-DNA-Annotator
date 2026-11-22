#!/usr/bin/env python3
"""Capa 1 — tabla de afinidad relativa para los pares represor–sitio.

Para cada par (represor, sitio) con confianza S/C y evidencia experimental:
  - score PWM del núcleo de 14 pb (normalizado por el autoscore del consenso)
  - percentil frente a 10.000 ventanas aleatorias del genoma de E. coli
  - z frente a 200 permutaciones del propio sitio (misma composición)

Unidad: **relativa** (percentil / z). No es kcal/mol: no hay fuente pública que permita
esa conversión sin rehacer el pipeline de Reg-Seq (ver PLAN_ALOSTERIA §8.9).

Salida: data/allostery/affinity_table.tsv + resumen JSON.

Uso: python scripts/allostery/17_affinity_table.py
"""
from __future__ import annotations

import ast
import collections
import csv
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "allostery"
RD = ROOT / "data" / "regulondb"
GBFF = ROOT / "data" / "test_ecoli.gbff"
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
CORE_W = 14
N_BG = 10_000
N_SHUF = 200


def genome() -> str:
    """Secuencia del cromosoma vía Bio.SeqIO (el gbff puede traer varios registros)."""
    from Bio import SeqIO
    best = None
    for rec in SeqIO.parse(str(GBFF), "genbank"):
        if best is None or len(rec.seq) > len(best.seq):
            best = rec
    return str(best.seq).upper()


def load_cores():
    out = {}
    for r in csv.DictReader(open(RD / "gt" / "tf_pwms.tsv"), delimiter="\t"):
        try:
            m = ast.literal_eval(r["pwm"])
        except Exception:
            continue
        if not m or len(m) < CORE_W:
            continue
        info = [max(row) - min(row) for row in m]
        bi = max(range(len(m) - CORE_W + 1), key=lambda i: sum(info[i:i + CORE_W]))
        core = np.array(m[bi:bi + CORE_W])
        cons = "".join("ACGT"[int(np.argmax(row))] for row in core)
        self_s = float(sum(core[i][IDX[cons[i]]] for i in range(CORE_W)))
        if self_s <= 0:
            continue
        out[r["tf"]] = {"core": core, "self": self_s}
    return out


def encode(seq: str) -> np.ndarray:
    return np.array([IDX.get(c, -1) for c in seq])


def best_scores(core: np.ndarray, encoded: np.ndarray, self_s: float) -> np.ndarray:
    """Score máximo sobre desplazamientos, vectorizado sobre filas (ventanas)."""
    if encoded.ndim == 1:
        encoded = encoded[None, :]
    L = encoded.shape[1]
    best = np.full(encoded.shape[0], -1e9)
    for s in range(L - CORE_W + 1):
        acc = np.zeros(encoded.shape[0])
        ok = np.ones(encoded.shape[0], dtype=bool)
        for i in range(CORE_W):
            b = encoded[:, s + i]
            ok &= b >= 0
            acc += np.where(b >= 0, core[i][np.clip(b, 0, 3)], 0.0)
        acc = np.where(ok, acc, -1e9)
        best = np.maximum(best, acc)
    return best / self_s


def main():
    g = genome()
    cores = load_cores()
    rng = np.random.default_rng(0)

    rows = [r for r in csv.DictReader(open(RD / "TF-RISet.tsv"), delimiter="\t")
            if r["2)type"] == "TF-promoter" and r["11)riFunction"] == "repressor"
            and r["20)confidenceLevel"] in ("S", "C")
            and any(t.startswith("EXP") for t in r["21)tfrsEvidence"].split("|"))]
    print(f"pares represor S/C con evidencia experimental: {len(rows)}")

    # fondo: ventanas genómicas del largo típico de sitio
    Ls = [len(max(re.findall(r"[ACGT]{6,}", r["10)tfrsSeq"]), key=len, default="")) for r in rows]
    L = int(np.median([x for x in Ls if x])) or 20
    starts = rng.integers(0, len(g) - L, N_BG)
    bg = np.array([encode(g[s:s + L]) for s in starts])
    bg = bg[(bg >= 0).all(1)]
    print(f"ventanas de fondo: {bg.shape[0]} de {L} pb")

    out = []
    for r in rows:
        tf = r["4)regulatorName"]
        if tf not in cores:
            continue
        site = max(re.findall(r"[ACGT]{6,}", r["10)tfrsSeq"]), key=len, default="")
        if len(site) < CORE_W:
            continue
        c = cores[tf]
        s_site = float(best_scores(c["core"], encode(site), c["self"])[0])
        s_bg = best_scores(c["core"], bg, c["self"])
        pct = float((s_bg < s_site).mean())
        shuf_seq = []
        for _ in range(N_SHUF):
            sh = list(site)
            rng.shuffle(sh)
            shuf_seq.append(encode("".join(sh)))
        shuf = best_scores(c["core"], np.array(shuf_seq), c["self"])
        z = float((s_site - shuf.mean()) / shuf.std()) if shuf.std() > 0 else 0.0
        out.append({"tf": tf, "sitio": site, "len": len(site),
                    "score": round(s_site, 4),
                    "percentil_genoma": round(pct, 4),
                    "z_shuffle": round(z, 2),
                    "conf": r["20)confidenceLevel"],
                    "evidencia": r["21)tfrsEvidence"][:60],
                    "pmid": r["26)tfrsPMIDS"][:40]})

    with open(OUT / "affinity_table.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    pcts = np.array([o["percentil_genoma"] for o in out])
    zs = np.array([o["z_shuffle"] for o in out])
    per_tf = collections.defaultdict(list)
    for o in out:
        per_tf[o["tf"]].append(o["percentil_genoma"])
    tf_med = {t: float(np.median(v)) for t, v in per_tf.items() if len(v) >= 3}
    summary = {
        "n_pares": len(out),
        "n_tfs": len(per_tf),
        "percentil_mediano": round(float(np.median(pcts)), 3),
        "fraccion_top5pct": round(float((pcts >= 0.95).mean()), 3),
        "z_mediano": round(float(np.median(zs)), 2),
        "fraccion_z_mayor_2": round(float((zs > 2).mean()), 3),
        "tfs_con_mediana_top5pct": sum(1 for v in tf_med.values() if v >= 0.95),
        "mejores_tfs": sorted(tf_med.items(), key=lambda kv: -kv[1])[:8],
        "peores_tfs": sorted(tf_med.items(), key=lambda kv: kv[1])[:5],
    }
    (OUT / "affinity_table.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
