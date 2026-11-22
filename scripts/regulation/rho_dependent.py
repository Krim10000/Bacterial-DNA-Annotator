#!/usr/bin/env python3
"""rho_dependent.py — detector de terminadores dependientes de Rho.

Por qué una herramienta aparte
------------------------------
Un terminador dependiente de Rho **no tiene horquilla**: lo desenrolla la
proteína Rho, que se engancha en una región del ARN rica en C y pobre en G
(sitios *rut*) y sin estructura estable. Por eso el detector de horquillas no
puede encontrarlos por definición —41% de "aciertos" sobre los 37 curados es
puro solapamiento— y mezclar ambos en una sola métrica no dice nada.

Qué mide
--------
Sobre la ventana de 100 nt aguas arriba del punto de terminación:

- `c_frac`: fracción de C  (los sitios rut son C-ricos)
- `g_frac`: fracción de G  (Rho evita regiones G-ricas)
- `c_minus_g`: asimetría C−G
- `longest_c`: racha de C más larga
- `dg`: energía libre de plegado (menos negativo = menos estructura)

Las dos primeras son las que discriminan (AUC 0,70 y 0,60 medidos sobre la
curaduría); el resto afina.

Uso
---
    # evaluar contra la curaduría y calibrar
    python scripts/regulation/rho_dependent.py --mode evaluate

    # escanear regiones candidatas
    python scripts/regulation/rho_dependent.py --mode scan \
        --fasta data/terminators/regions_k12.fna --out data/terminators/rho_calls.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "regulation"))

WINDOW = 100


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def longest_run(s: str, base: str) -> int:
    best = cur = 0
    for ch in s:
        cur = cur + 1 if ch == base else 0
        best = max(best, cur)
    return best


def fold_energy(s: str) -> float:
    try:
        import RNA
        return float(RNA.fold(s)[1])
    except Exception:
        return 0.0


def features(win: str) -> dict:
    win = win.upper().replace("N", "")
    n = max(1, len(win))
    c = win.count("C") / n
    g = win.count("G") / n
    # motivos rut: Rho se une a regiones C-ricas y pobres en G. Los patrones
    # CAA/CAC y las rachas de C son los descritos en la literatura.
    cc = win.count("CC") / n
    caa = (win.count("CAA") + win.count("CAC")) / n
    a = win.count("A") / n
    return {"c_frac": round(c, 4), "g_frac": round(g, 4),
            "c_minus_g": round(c - g, 4),
            "longest_c": longest_run(win, "C"),
            "longest_g": longest_run(win, "G"),
            "cc_frac": round(cc, 4), "caa_cac_frac": round(caa, 4),
            "a_frac": round(a, 4),
            "dg": round(fold_energy(win), 2)}


def upstream(seq: str, pos: int, strand: str, w: int = WINDOW) -> str:
    """Ventana aguas arriba del TTS, en la hebra del transcrito."""
    if strand == "+":
        return seq[max(0, pos - w):pos]
    return revcomp(seq[pos:pos + w])


def load_curated():
    from Bio import SeqIO
    rec = next(SeqIO.parse(str(ROOT / "data" / "test_ecoli.gbff"), "genbank"))
    seq = str(rec.seq).upper()
    pos, neg = [], []
    for r in csv.DictReader(open(ROOT / "data" / "regulondb" / "TerminatorSet.tsv"),
                            delimiter="\t"):
        p = (r.get("2)leftEndPos") or "").strip()
        if not p.isdigit():
            continue
        st = (r.get("4)strand") or "").strip().lower()
        stn = "+" if st.startswith("f") else "-"
        tipo = (r.get("7)type") or "").strip()
        win = upstream(seq, int(p), stn)
        if len(win) < 60:
            continue
        f = features(win)
        if tipo == "rho-dependent":
            pos.append(f)
        elif tipo == "rho-independent":
            neg.append(f)
    return seq, pos, neg


def random_negatives(seq: str, n: int = 5000) -> list[dict]:
    """Negativos en regiones intergénicas al azar: la tasa real de FP."""
    random.seed(11)
    out = []
    while len(out) < n:
        p = random.randint(200, len(seq) - 200)
        st = random.choice("+-")
        win = upstream(seq, p, st)
        if len(win) >= 60:
            out.append(features(win))
    return out


def evaluate() -> dict:
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    seq, pos, neg_ind = load_curated()
    neg_rand = random_negatives(seq)
    print(f"positivos (rho-dependiente): {len(pos)}")
    print(f"negativos curados (rho-independiente): {len(neg_ind)}")
    print(f"negativos aleatorios (intergénico): {len(neg_rand)}")

    keys = ["c_frac", "g_frac", "c_minus_g", "longest_c", "longest_g",
            "cc_frac", "caa_cac_frac", "a_frac", "dg"]
    neg = neg_ind + neg_rand
    X = np.array([[f[k] for k in keys] for f in pos + neg])
    y = np.array([1] * len(pos) + [0] * len(neg))

    n_pos = len(pos)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
    aucs, sens = [], []
    for tr, te in skf.split(X, y):
        clf = LogisticRegression(max_iter=1000, class_weight="balanced").fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        # sensibilidad con 5% de FP: umbral en el cuantil 95 de los negativos de test
        neg_te = p[y[te] == 0]
        thr = np.quantile(neg_te, 0.95)
        aucs.append(roc_auc_score(y[te], p))
        sens.append(float((p[y[te] == 1] >= thr).mean()))
    clf = LogisticRegression(max_iter=1000, class_weight="balanced").fit(X, y)
    coef = {k: round(float(c), 3) for k, c in zip(keys, clf.coef_[0])}
    res = {"n_pos": len(pos), "n_neg_curados": len(neg_ind),
           "n_neg_aleatorios": len(neg_rand),
           "auc_cv": round(float(np.mean(aucs)), 3),
           "auc_cv_sd": round(float(np.std(aucs)), 3),
           "sens_at_5pct_fpr": round(float(np.mean(sens)), 3),
           "sens_sd": round(float(np.std(sens)), 3),
           "pesos": coef}
    print(json.dumps(res, indent=1, ensure_ascii=False))
    return res


def scan(fasta: Path, out: Path, threshold: float) -> int:
    """Puntúa cada región candidata con un score simple e interpretable."""
    n = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["region", "offset", "score", "c_frac", "g_frac", "longest_c", "call"])
        name = ""
        buf = []
        for line in open(fasta):
            if line.startswith(">"):
                if name and buf:
                    n += _emit(name, "".join(buf), w, threshold)
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if name and buf:
            n += _emit(name, "".join(buf), w, threshold)
    print(f"{n} llamadas rho-dependientes -> {out}")
    return 0


def _emit(name: str, seq: str, w, threshold: float) -> int:
    # score interpretable: asimetría C-G + racha de C normalizada
    f = features(seq)
    score = f["c_minus_g"] + 0.02 * f["longest_c"]
    if score >= threshold:
        w.writerow([name.split("|")[0], name.split("|")[1] if "|" in name else "",
                    round(score, 4), f["c_frac"], f["g_frac"], f["longest_c"],
                    "rho_dependent_candidate"])
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["evaluate", "scan"], default="evaluate")
    ap.add_argument("--fasta")
    ap.add_argument("--out")
    ap.add_argument("--threshold", type=float, default=0.10,
                    help="score mínimo en modo scan (asimetría C-G + racha)")
    a = ap.parse_args()
    if a.mode == "evaluate":
        evaluate()
        return 0
    if not (a.fasta and a.out):
        print("scan necesita --fasta y --out")
        return 1
    return scan(Path(a.fasta), Path(a.out), a.threshold)


if __name__ == "__main__":
    sys.exit(main())
