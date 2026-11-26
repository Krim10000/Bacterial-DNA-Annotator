#!/usr/bin/env python3
"""Entrena un calibrador k-mer específico de especie (sin transferencia)."""
from __future__ import annotations
import argparse, importlib, json, sys
from pathlib import Path
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
E = importlib.import_module("scripts.sigma.50_eval_layers")
from webapp.sigma_kmer import KmerModel, save_model  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gold", type=Path, required=True)
    ap.add_argument("--source", default="auto")
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--up", type=int, default=40)
    ap.add_argument("--down", type=int, default=5)
    ap.add_argument("--ngram", type=int, default=4)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    rows = E.load_gold(a.gold, a.source)
    seqs = E.load_sequence(a.data_root, a.genome) or E._fallback_sequence(a.data_root, a.genome)
    if not seqs:
        raise SystemExit(f"sin secuencia para {a.genome}")
    labels = sorted({s for r in rows for s in r["sigmas"]})
    if not labels:
        raise SystemExit("gold sin etiquetas mapeadas")
    texts = []
    for r in rows:
        contig = r["contig"] or next(iter(seqs), "")
        texts.append(E._window(seqs.get(contig, ""), r["tss"], r["strand"], a.up, a.down))
    vec = CountVectorizer(analyzer="char", ngram_range=(a.ngram, a.ngram), min_df=2)
    X = vec.fit_transform(texts)
    vocab = {str(k): int(v) for k, v in vec.vocabulary_.items()}
    scaler = StandardScaler(with_mean=False).fit(X)
    Xs = scaler.transform(X)
    coef, intercept = [], []
    for s in labels:
        y = np.array([1 if s in r["sigmas"] else 0 for r in rows], dtype=int)
        if len(set(y)) < 2:
            coef.append([0.0] * Xs.shape[1]); intercept.append(-20.0)
            continue
        clf = LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0)
        clf.fit(Xs, y)
        coef.append([float(v) for v in clf.coef_[0]])
        intercept.append(float(clf.intercept_[0]))
    model = KmerModel(
        genome=a.genome, labels=tuple(labels), vocab=vocab,
        scale=tuple(float(v) for v in scaler.scale_),
        coef=tuple(tuple(row) for row in coef), intercept=tuple(intercept),
        ngram=a.ngram, window_up=a.up, window_down=a.down,
        trained_on=str(a.gold),
    )
    save_model(model, a.out)
    print(json.dumps({
        "genome": a.genome, "n": len(rows), "labels": labels,
        "vocab": len(vocab), "window": [-a.up, a.down], "ngram": a.ngram,
        "out": str(a.out),
    }, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
