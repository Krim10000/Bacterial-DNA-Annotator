#!/usr/bin/env python3
"""F1b — baseline ESM-2 con LOFO (leave-one-family-out).

Predice el modo del efector (inductor / correpresor) desde la secuencia de proteína.

Protocolo fijado antes de mirar resultados:
  - Representación: ESM-2 (8M local, CPU), media de residuos de la última capa.
  - Modelo: StandardScaler -> PCA(10) -> LogisticRegression(C=0.1, balanced).
  - Validación principal: LOFO por familia de dominio Pfam.
  - Contraste: 5-fold estratificado aleatorio (con fuga de familia) y LOO por TF.
  - Nulo: permutación de etiquetas con el mismo LOFO (500 réplicas).

Optimización: la proyección PCA de cada fold depende solo de X, no de y. Se
precomputa una vez y las permutaciones solo reajustan la regresión logística.

Puerta G1: si LOFO balanced acc >= 0.80, el bolsillo no se justifica para el modo.
Salida: data/allostery/baseline_esm.json + baseline_esm.md
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import ROOT, WORK, load_decisive, load_families  # noqa: E402

OUT = ROOT / "data" / "allostery"
COR, IND = "correpresor", "inductor"


def embed_cached(seqs, tag="esm8m"):
    cache = WORK / f"emb_{tag}.npz"
    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        if len(d["X"]) == len(seqs):
            return d["X"]
    import esm
    import torch

    torch.hub.set_dir(str(ROOT / ".cache" / "torch" / "hub"))
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    bc = alphabet.get_batch_converter()
    out = []
    with torch.no_grad():
        for i in range(0, len(seqs), 8):
            batch = [(f"s{j}", s) for j, s in enumerate(seqs[i:i + 8])]
            _, _, toks = bc(batch)
            rep = model(toks, repr_layers=[model.num_layers], return_contacts=False)
            h = rep["representations"][model.num_layers]
            for j, (_, s) in enumerate(batch):
                out.append(h[j, 1:len(s) + 1].mean(0).numpy())
    X = np.vstack(out)
    np.savez(cache, X=X)
    return X


def make_folds(X, groups, n_pca=10):
    """Proyecta cada fold una sola vez: (train_idx, test_idx, Ztr, Zte)."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    folds = []
    for g in sorted(set(groups)):
        te = np.array([i for i, gg in enumerate(groups) if gg == g])
        tr = np.array([i for i, gg in enumerate(groups) if gg != g])
        if len(tr) < 5 or len(te) == 0:
            continue
        sc = StandardScaler().fit(X[tr])
        k = max(2, min(n_pca, len(tr) - 1, X.shape[1]))
        pca = PCA(n_components=k, random_state=0).fit(sc.transform(X[tr]))
        folds.append((tr, te, pca.transform(sc.transform(X[tr])),
                      pca.transform(sc.transform(X[te]))))
    return folds


def make_folds_kfold(X, y, n_splits=5, seed=0):
    from sklearn.decomposition import PCA
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler

    folds = []
    for tr, te in StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X, y):
        sc = StandardScaler().fit(X[tr])
        pca = PCA(n_components=10, random_state=0).fit(sc.transform(X[tr]))
        folds.append((tr, te, pca.transform(sc.transform(X[tr])),
                      pca.transform(sc.transform(X[te]))))
    return folds


def make_folds_loo(X):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    folds = []
    for i in range(len(X)):
        tr = np.array([j for j in range(len(X)) if j != i])
        sc = StandardScaler().fit(X[tr])
        pca = PCA(n_components=10, random_state=0).fit(sc.transform(X[tr]))
        folds.append((tr, np.array([i]), pca.transform(sc.transform(X[tr])),
                      pca.transform(sc.transform(X[[i]]))))
    return folds


def predict(folds, y, C=0.1):
    from sklearn.linear_model import LogisticRegression

    pred = np.empty(len(y), dtype=object)
    for tr, te, Ztr, Zte in folds:
        lr = LogisticRegression(C=C, class_weight="balanced", solver="liblinear", max_iter=500)
        lr.fit(Ztr, y[tr])
        pred[te] = lr.predict(Zte)
    return pred


def metrics(y, pred):
    y = np.asarray(y); pred = np.asarray(pred)
    tp = int(np.sum((y == COR) & (pred == COR))); tn = int(np.sum((y == IND) & (pred == IND)))
    fp = int(np.sum((y == IND) & (pred == COR))); fn = int(np.sum((y == COR) & (pred == IND)))
    rec_c = tp / (tp + fn) if tp + fn else 0.0
    rec_i = tn / (tn + fp) if tn + fp else 0.0
    den = ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    return {"balanced_acc": (rec_c + rec_i) / 2,
            "mcc": (tp * tn - fp * fn) / den if den else 0.0,
            "recall_correpresor": rec_c, "recall_inductor": rec_i,
            "tp": tp, "tn": tn, "fp": fp, "fn": fn}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--C", type=float, default=0.1)
    args = ap.parse_args()

    data = [d for d in load_decisive() if d["seq"]]
    fams = load_families()
    for d in data:
        d["family"] = fams.get(d["tf"], "sin_dominio")
    y = np.array([d["mode"] for d in data])
    groups = [d["family"] for d in data]
    X = embed_cached([d["seq"] for d in data])
    print(f"embeddings {X.shape} · familias {len(set(groups))}", file=sys.stderr)

    folds = make_folds(X, groups)
    main_res = metrics(y, predict(folds, y, args.C))

    leak = metrics(y, predict(make_folds_kfold(X, y), y, args.C))
    loo = metrics(y, predict(make_folds_loo(X), y, args.C))

    rng = np.random.default_rng(0)
    null = []
    for i in range(args.perms):
        yp = rng.permutation(y)
        null.append(metrics(yp, predict(folds, yp, args.C))["balanced_acc"])
        if (i + 1) % 50 == 0:
            print(f"  perm {i + 1}/{args.perms}", file=sys.stderr)
    null = np.array(null)
    pval = float((np.sum(null >= main_res["balanced_acc"]) + 1) / (len(null) + 1))

    folds_report = []
    pred_main = predict(folds, y, args.C)
    for g in sorted(set(groups)):
        te = np.array([i for i, gg in enumerate(groups) if gg == g])
        tr = np.array([i for i, gg in enumerate(groups) if gg != g])
        if len(tr) < 5 or len(te) == 0:
            continue
        folds_report.append({"family": g, "n": int(len(te)),
                             "acc": float(np.mean(pred_main[te] == y[te]))})

    summary = {
        "n": len(data), "n_familias": len(set(groups)),
        "clases": dict(collections.Counter(y.tolist())),
        "modelo": {"esm": "esm2_t6_8M_UR50D", "pca": 10, "C": args.C,
                   "pooling": "mean_residuos"},
        "lofo": main_res, "lofo_folds": folds_report,
        "leaky_5fold": leak, "loo_tf": loo,
        "nulo_permutacion": {"media": float(null.mean()), "sd": float(null.std()),
                             "p95": float(np.percentile(null, 95)),
                             "p_valor": pval, "reps": len(null)},
        "puerta_G1": ("pasa — la secuencia sola predice el modo"
                      if main_res["balanced_acc"] >= 0.80
                      else "no pasa — F2 (bolsillo) se justifica"),
    }
    (OUT / "baseline_esm.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    L = ["# F1b — baseline ESM-2 (LOFO por familia)", "",
         f"n={len(data)} · familias={len(set(groups))} · clases={dict(collections.Counter(y.tolist()))}", "",
         "| validación | balanced acc | MCC | recall cor | recall ind |", "|---|---:|---:|---:|---:|"]
    for name, m in (("LOFO (familia)", main_res), ("5-fold aleatorio (fuga)", leak), ("LOO por TF", loo)):
        L.append(f"| {name} | {m['balanced_acc']:.3f} | {m['mcc']:.3f} | "
                 f"{m['recall_correpresor']:.3f} | {m['recall_inductor']:.3f} |")
    L += ["", f"Nulo por permutación ({len(null)} réplicas, mismo LOFO): "
          f"{null.mean():.3f} ± {null.std():.3f}, p95={np.percentile(null, 95):.3f}, p={pval:.3f}",
          "", "## Por familia (LOFO)", "", "| familia | n test | acc |", "|---|---:|---:|"]
    for f in folds_report:
        L.append(f"| {f['family']} | {f['n']} | {f['acc']:.2f} |")
    L += ["", f"**Puerta G1:** {summary['puerta_G1']}", ""]
    (OUT / "baseline_esm.md").write_text("\n".join(L) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "lofo_folds"},
                     indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
