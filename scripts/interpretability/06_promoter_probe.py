#!/usr/bin/env python3
"""Sonda de promotor para embeddings pretrained.

Carga <run>/embeddings.npy + windows.tsv con columna `label` (1=TATAAT,
0=dinucleótido barajado), entrena regresión logística con CV estratificada y
reporta AUROC/accuracy. Incluye baseline GC.

Uso:
  python scripts/interpretability/06_promoter_probe.py \
      --run data/embeddings/pretrained/prokbert_promoter
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    run = Path(args.run)
    meta = pd.read_csv(run / "windows.tsv", sep="\t")
    X = np.load(run / "embeddings.npy")
    y = meta["label"].astype(int).values
    if len(meta) != len(X):
        raise SystemExit("Desalineado")

    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    clf = make_pipeline(StandardScaler(),
                        LogisticRegression(max_iter=2000, class_weight="balanced"))
    auc = cross_val_score(clf, X, y, cv=cv, scoring="roc_auc")
    acc = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")

    # Baseline GC (1 feature).
    seqs = []
    for r in meta.itertuples():
        seqs.append(0.0)
    gc = np.array([0.0])  # placeholder si no hay secuencia
    # Si el windows.tsv trae gc, usarlo.
    if "gc" in meta.columns:
        gc = meta["gc"].astype(float).values.reshape(-1, 1)
        auc_gc = cross_val_score(clf, gc, y, cv=cv, scoring="roc_auc")
    else:
        auc_gc = np.array([np.nan])

    res = {
        "run": str(run),
        "n": int(len(y)),
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auroc_mean": float(np.mean(auc)),
        "auroc_std": float(np.std(auc)),
        "accuracy_mean": float(np.mean(acc)),
        "auroc_gc_baseline": float(np.mean(auc_gc)) if np.isfinite(auc_gc).any() else None,
    }
    out = Path(args.out) if args.out else run
    (out / "promoter_probe.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
