#!/usr/bin/env python3
"""Firma multivariada por componente de ADN.

Idea: cada tipo de elemento (CDS, tRNA, rRNA, intergénico, upstream, …) no se
describe por una sola dimensión, sino por un **patrón de 384 valores**. Este
script:

  1. Asigna cada ventana a un componente.
  2. Calcula el centroide (media de 384 dims) y el patrón estandarizado
     (d de Cohen) de cada componente frente al resto.
  3. Identifica qué dimensiones del patrón son consistentes entre genomas y
     significativas (p de ventana, p meta por genoma, FDR, consistencia de
     signo).
  4. Evalúa si el patrón permite clasificar el componente con validación
     leave-one-genome-out (LOGO), comparado con baseline mayoritario y GC.
  5. Guarda el patrón disperso (dims seleccionadas) para usar como firma.

Salidas principales:
  class_centroids.tsv        media cruda de cada componente (384 dims)
  class_cohen_d.tsv          patrón estandarizado (384 dims por componente)
  class_signature_dims.tsv   tabla larga con estadística por dim
  class_signature_selected.tsv  dims que forman la firma de cada componente
  signature_scores.tsv       score de cada ventana contra cada firma
  classification.json        LOGO: accuracy, balanced accuracy, confusión
  heatmap_signatures.png     mapa de calor de las firmas
  report.md

Uso:
  python scripts/interpretability/03_component_signatures.py \
      --run data/embeddings/v2_epoch9_pilot6

  # ventana corta
  python scripts/interpretability/03_component_signatures.py \
      --run data/embeddings/v2_epoch9_pilot6_256
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
MPL_DIR = ROOT / ".cache" / "matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

BASE_CLASSES = ["CDS", "tRNA", "rRNA", "intergenic"]
EXTRA_CLASSES = {
    "upstream": "is_upstream",
    "transposon": "is_transposon",
    "prophage": "is_prophage",
    "crispr": "is_crispr",
}


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    q = np.empty(n, dtype=float)
    q[order] = np.minimum.accumulate((pv[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    out[ok] = np.clip(q, 0.0, 1.0)
    return out


def cohen_d(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ma, mb = a.mean(axis=0), b.mean(axis=0)
    va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    pooled = np.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) /
                     max(1, len(a) + len(b) - 2))
    out = np.zeros(a.shape[1], dtype=float)
    np.divide(ma - mb, pooled, out=out, where=pooled > 0)
    return out


def assign_components(df: pd.DataFrame) -> pd.Series:
    """Componente principal por contenido en el centro de la ventana."""
    comp = pd.Series("", index=df.index, dtype=object)
    cds = ((df["feature_type"] == "CDS") & (df["overlap_cds"] >= 0.5) &
           (df["overlap_trna"] < 0.10) & (df["overlap_rrna"] < 0.10))
    trna = (df["feature_type"] == "tRNA") & (df["overlap_trna"] >= 0.10)
    rrna = (df["feature_type"] == "rRNA") & (df["overlap_rrna"] >= 0.30)
    ig = df["feature_type"] == "intergenic"
    # Prioridad: rRNA > tRNA > CDS > intergénico.
    comp[ig] = "intergenic"
    comp[cds] = "CDS"
    comp[trna] = "tRNA"
    comp[rrna] = "rRNA"
    return comp


def clean_intergenic_mask(df: pd.DataFrame) -> pd.Series:
    return ((df["feature_type"] == "intergenic") &
            (df["overlap_cds"] < 0.20) &
            (df["overlap_trna"] < 0.05) &
            (df["overlap_rrna"] < 0.05) &
            (df["overlap_ncrna"] < 0.05))


def center_by_group(A: np.ndarray, groups: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=np.float64).copy()
    for g in np.unique(groups):
        idx = groups == g
        if idx.sum() == 0:
            continue
        if A.ndim == 1:
            A[idx] -= A[idx].mean()
        else:
            A[idx] -= A[idx].mean(axis=0, keepdims=True)
    return A


def adjust_matrix(X: np.ndarray, gc: np.ndarray, groups: np.ndarray,
                  mode: str) -> np.ndarray:
    """Elimina GC/GC^2 y, opcionalmente, la media por genoma.

    - none: sin ajuste.
    - gc: residualiza sobre [1, gc, gc^2] global.
    - gc_genome: centra por genoma y residualiza sobre [gc_c, gc_c^2].
    """
    X = X.astype(np.float64)
    if mode == "none":
        return X
    if mode == "gc":
        C = np.column_stack([np.ones(len(gc)), gc, gc * gc])
        beta, *_ = np.linalg.lstsq(C, X, rcond=None)
        return X - C @ beta
    if mode == "gc_genome":
        Xc = center_by_group(X, groups)
        gc_c = center_by_group(gc, groups)
        C = np.column_stack([gc_c, gc_c * gc_c])
        beta, *_ = np.linalg.lstsq(C, Xc, rcond=None)
        return Xc - C @ beta
    raise ValueError(f"adjust desconocido: {mode}")


# ---------------------------------------------------------------------------
# Firma de un componente
# ---------------------------------------------------------------------------
def signature_stats(X: np.ndarray, mask: np.ndarray, groups: np.ndarray,
                    name: str, d_min: float, consistency_min: float,
                    min_genomes: int = 4) -> dict | None:
    n = int(mask.sum())
    if n < 20 or (~mask).sum() < 20:
        return None
    Xc, Xr = X[mask].astype(np.float64), X[~mask].astype(np.float64)
    mean_c, mean_r = Xc.mean(axis=0), Xr.mean(axis=0)
    diff = mean_c - mean_r
    d = cohen_d(Xc, Xr)

    # p por ventana (Welch) y FDR dentro del componente.
    with np.errstate(invalid="ignore", divide="ignore"):
        _, p_win = stats.ttest_ind(Xc, Xr, axis=0, equal_var=False)
    p_win = np.nan_to_num(p_win, nan=1.0)
    q_win = bh(p_win)

    # Diferencia por genoma y meta-análisis (unidad = genoma).
    per_genome = []
    for g in np.unique(groups):
        gc = mask & (groups == g)
        gr = (~mask) & (groups == g)
        if gc.sum() >= 2 and gr.sum() >= 2:
            per_genome.append(X[gc].mean(axis=0) - X[gr].mean(axis=0))
    D = np.array(per_genome) if per_genome else np.zeros((0, X.shape[1]))
    n_genomes = len(D)
    meta_diff = D.mean(axis=0) if n_genomes else np.full(X.shape[1], np.nan)
    if n_genomes >= 3:
        meta_t, meta_p = stats.ttest_1samp(D, 0.0, axis=0)
        meta_p = np.nan_to_num(meta_p, nan=1.0)
        meta_q = bh(meta_p)
        sign = np.sign(meta_diff)
        consistency = np.where(sign >= 0, (D > 0).mean(axis=0),
                               (D < 0).mean(axis=0))
    else:
        meta_t = np.full(X.shape[1], np.nan)
        meta_p = np.ones(X.shape[1])
        meta_q = np.ones(X.shape[1])
        consistency = np.full(X.shape[1], np.nan)

    # Criterio principal: significancia por ventana + tamaño de efecto +
    # consistencia entre genomas. La meta-ANOVA por genoma (q_meta) se reporta
    # aparte porque con 6 genomas tiene poca potencia.
    selected_window = ((q_win < 0.01) & (np.abs(d) >= d_min) &
                       (consistency >= consistency_min) &
                       (n_genomes >= min_genomes))
    selected_meta = selected_window & (meta_q < 0.05)
    selected = selected_window

    return {
        "name": name, "n": n, "n_rest": int((~mask).sum()),
        "n_genomes": n_genomes,
        "mean_c": mean_c, "mean_r": mean_r, "diff": diff, "d": d,
        "p_win": p_win, "q_win": q_win,
        "meta_diff": meta_diff, "meta_t": meta_t, "meta_p": meta_p,
        "meta_q": meta_q, "consistency": consistency,
        "selected": selected, "selected_window": selected_window,
        "selected_meta": selected_meta,
    }


# ---------------------------------------------------------------------------
# Clasificación LOGO
# ---------------------------------------------------------------------------
def classify_logo(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                  classes: list[str], permutations: int = 0, seed: int = 42):
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                                 confusion_matrix, classification_report,
                                 roc_auc_score)
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    def make_pipeline():
        return Pipeline([
            ("scaler", StandardScaler()),
            ("select", SelectKBest(f_classif, k=min(50, X.shape[1]))),
            ("clf", LogisticRegression(max_iter=500, tol=1e-3,
                                       class_weight="balanced",
                                       multi_class="multinomial")),
        ])

    def run(y_use, collect=False):
        pred = np.empty(len(y_use), dtype=object)
        proba = np.zeros((len(y_use), len(classes)))
        for g in np.unique(groups):
            test = groups == g
            train = ~test
            if len(np.unique(y_use[train])) < 2:
                pred[test] = classes[0]
                continue
            pipe = make_pipeline()
            pipe.fit(X[train], y_use[train])
            pred[test] = pipe.predict(X[test])
            if collect:
                p = pipe.predict_proba(X[test])
                # Alinear columnas con `classes`.
                for j, c in enumerate(pipe.classes_):
                    proba[test, classes.index(c)] = p[:, j]
        pred = np.asarray(pred)
        acc = accuracy_score(y_use, pred)
        bal = balanced_accuracy_score(y_use, pred)
        if collect:
            return pred, proba, acc, bal
        return acc, bal

    pred, proba, acc, bal = run(y, collect=True)
    labels = classes
    cm = confusion_matrix(y, pred, labels=labels)
    rep = classification_report(y, pred, labels=labels, output_dict=True,
                                zero_division=0)
    auc = {}
    for j, c in enumerate(labels):
        y_bin = (y == c).astype(int)
        if len(np.unique(y_bin)) > 1:
            auc[c] = float(roc_auc_score(y_bin, proba[:, j]))

    # Baseline: clase mayoritaria del train en cada fold.
    base_pred = np.empty(len(y), dtype=object)
    for g in np.unique(groups):
        test = groups == g
        train = ~test
        vals, counts = np.unique(y[train], return_counts=True)
        base_pred[test] = vals[np.argmax(counts)]
    base_acc = accuracy_score(y, base_pred)
    base_bal = balanced_accuracy_score(y, base_pred)

    # Baseline GC-only.
    gc = None
    if "gc" in X.shape:  # placeholder; GC se pasa aparte si existe
        pass

    perm = None
    if permutations > 0:
        rng = np.random.default_rng(seed)
        null = []
        for _ in range(permutations):
            yp = y.copy()
            for g in np.unique(groups):
                idx = np.flatnonzero(groups == g)
                yp[idx] = rng.permutation(yp[idx])
            _, b = run(yp)
            null.append(b)
        null = np.array(null)
        perm = {
            "n": permutations,
            "balanced_accuracy_null_mean": float(null.mean()),
            "balanced_accuracy_null_p95": float(np.quantile(null, 0.95)),
            "p_value": float((1 + (null >= bal).sum()) / (1 + len(null))),
        }

    return {
        "classes": labels,
        "n": int(len(y)),
        "accuracy": float(acc),
        "balanced_accuracy": float(bal),
        "baseline_majority_accuracy": float(base_acc),
        "baseline_majority_balanced_accuracy": float(base_bal),
        "confusion_matrix": cm.tolist(),
        "classification_report": rep,
        "auc_ovr": auc,
        "permutation": perm,
        "predictions": pred,
        "probabilities": proba,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--d-min", type=float, default=0.25,
                    help="|d| mínimo para seleccionar una dim de la firma")
    ap.add_argument("--consistency-min", type=float, default=0.80,
                    help="Fracción mínima de genomas con el mismo signo")
    ap.add_argument("--min-genomes", type=int, default=4)
    ap.add_argument("--top-dims", type=int, default=60,
                    help="Dims a mostrar en el heatmap")
    ap.add_argument("--adjust", choices=["none", "gc", "gc_genome"],
                    default="gc_genome",
                    help="Ajuste de composición antes de la firma (default: GC + genoma)")
    ap.add_argument("--permutations", type=int, default=0,
                    help="Permutaciones del clasificador (0 = desactivar)")
    ap.add_argument("--include-extra", action="store_true", default=True,
                    help="Incluye upstream/transposón/profago/CRISPR como one-vs-rest")
    args = ap.parse_args()

    run = Path(args.run)
    out = Path(args.out) if args.out else run
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(run / "windows.tsv", sep="\t")
    X_raw = np.load(run / "embeddings.npy").astype(np.float32)
    if len(meta) != X_raw.shape[0]:
        raise SystemExit(f"Desalineado: {len(meta)} vs {X_raw.shape[0]}")
    dims = [f"dim_{i:03d}" for i in range(X_raw.shape[1])]
    groups = meta["genome"].values
    gc = meta["gc"].astype(float).values if "gc" in meta.columns else np.zeros(len(meta))
    X = adjust_matrix(X_raw, gc, groups, args.adjust)
    print(f"Ajuste: {args.adjust} (GC + genoma si aplica)")

    comp = assign_components(meta)
    meta["component"] = comp
    ig_clean = clean_intergenic_mask(meta)
    meta["component_clean"] = comp.where(~ig_clean, "intergenic_clean")

    print(f"Ventanas: {len(meta)} | genomas: {meta['genome'].nunique()} | dims: {X.shape[1]}")
    print("Componentes base:")
    for c, n in comp.value_counts().items():
        print(f"  {c:18s} {n:5d}  genomas={meta.loc[comp == c, 'genome'].nunique()}")
    print(f"  {'intergenic_clean':18s} {int(ig_clean.sum()):5d}  "
          f"genomas={meta.loc[ig_clean, 'genome'].nunique()}")

    # Lista de componentes a analizar.
    targets = []
    for c in BASE_CLASSES:
        mask = (comp == c).values
        if mask.sum() >= 20 and meta.loc[mask, "genome"].nunique() >= args.min_genomes:
            targets.append((c, mask))
    if ig_clean.sum() >= 20 and meta.loc[ig_clean, "genome"].nunique() >= 3:
        targets.append(("intergenic_clean", ig_clean.values))
    if args.include_extra:
        for name, col in EXTRA_CLASSES.items():
            if col in meta.columns:
                mask = (meta[col].astype(bool)).values
                if mask.sum() >= 20 and meta.loc[mask, "genome"].nunique() >= 3:
                    targets.append((name, mask))

    if not targets:
        raise SystemExit("No hay componentes con suficientes ventanas.")

    # Firmas.
    signatures = []
    for name, mask in targets:
        res = signature_stats(X, mask, groups, name, args.d_min,
                              args.consistency_min, args.min_genomes)
        if res is None:
            print(f"  [skip] {name}: muy pocas ventanas")
            continue
        signatures.append(res)
        print(f"  firma {name:18s} n={res['n']:5d} "
              f"genomas={res['n_genomes']} dims seleccionadas={int(res['selected'].sum())}")

    if not signatures:
        raise SystemExit("No se pudo calcular ninguna firma.")

    # Matrices.
    # Centroides crudos (sin ajustar) para inspección; el patrón usa X ajustada.
    raw_centroids = {}
    target_map = dict(targets)
    for sig in signatures:
        mask = target_map[sig["name"]]
        raw_centroids[sig["name"]] = X_raw[mask].astype(np.float64).mean(axis=0)
    pd.DataFrame(raw_centroids, index=dims).T.to_csv(
        out / "class_centroids_raw.tsv", sep="\t")
    centroids = pd.DataFrame({s["name"]: s["mean_c"] for s in signatures}, index=dims).T
    cohen = pd.DataFrame({s["name"]: s["d"] for s in signatures}, index=dims).T
    diffs = pd.DataFrame({s["name"]: s["diff"] for s in signatures}, index=dims).T
    centroids.to_csv(out / "class_centroids.tsv", sep="\t")
    cohen.to_csv(out / "class_cohen_d.tsv", sep="\t")
    diffs.to_csv(out / "class_diff.tsv", sep="\t")

    # Firma dispersa en formato ancho: d si la dim fue seleccionada, 0 si no.
    sparse = {}
    mask = {}
    for s in signatures:
        sparse[s["name"]] = np.where(s["selected"], s["d"], 0.0)
        mask[s["name"]] = s["selected"].astype(int)
    pd.DataFrame(sparse, index=dims).T.to_csv(
        out / "class_signature_sparse.tsv", sep="\t")
    pd.DataFrame(mask, index=dims).T.to_csv(
        out / "class_signature_mask.tsv", sep="\t")

    # Tabla larga.
    rows = []
    for s in signatures:
        for i, dim in enumerate(dims):
            rows.append({
                "component": s["name"], "dim": dim, "dim_index": i,
                "n": s["n"], "n_rest": s["n_rest"], "n_genomes": s["n_genomes"],
                "mean_component": s["mean_c"][i],
                "mean_rest": s["mean_r"][i],
                "diff": s["diff"][i],
                "cohen_d": s["d"][i],
                "p_window": s["p_win"][i],
                "q_window": s["q_win"][i],
                "meta_diff": s["meta_diff"][i],
                "meta_t": s["meta_t"][i],
                "meta_p": s["meta_p"][i],
                "q_meta": s["meta_q"][i],
                "consistency": s["consistency"][i],
                "selected": int(s["selected"][i]),
                "selected_meta": int(s["selected_meta"][i]),
            })
    long = pd.DataFrame(rows)
    long.to_csv(out / "class_signature_dims.tsv", sep="\t", index=False)

    # Firma dispersa: solo dims seleccionadas, ordenadas por |d|.
    selected = long[long["selected"] == 1].copy()
    selected["abs_d"] = selected["cohen_d"].abs()
    selected = selected.sort_values(["component", "abs_d"], ascending=[True, False])
    selected.drop(columns=["abs_d"]).to_csv(
        out / "class_signature_selected.tsv", sep="\t", index=False)

    # Scores por ventana contra cada firma (coseno con el patrón d).
    # Centrar por la media global para que el coseno no dependa del offset.
    Xc = X - X.mean(axis=0, keepdims=True)
    sig_matrix = np.stack([s["d"] for s in signatures], axis=1)  # 384 x C
    sig_norm = np.linalg.norm(sig_matrix, axis=0)
    x_norm = np.linalg.norm(Xc, axis=1, keepdims=True)
    denom = np.where((x_norm * sig_norm[None, :]) > 0,
                     x_norm * sig_norm[None, :], 1.0)
    scores = (Xc @ sig_matrix) / denom
    score_df = pd.DataFrame(scores, columns=[s["name"] for s in signatures])
    score_df.insert(0, "window_id", meta["window_id"].values)
    score_df.insert(1, "genome", meta["genome"].values)
    score_df.insert(2, "true_component", meta["component"].values)
    score_df["predicted_component"] = score_df[
        [s["name"] for s in signatures]].idxmax(axis=1)
    score_df.to_csv(out / "signature_scores.tsv", sep="\t", index=False)

    # Clasificación LOGO con las clases base presentes en >= min_genomes.
    base_targets = [(name, mask) for name, mask in targets
                    if name in BASE_CLASSES and meta.loc[mask, "genome"].nunique() >= 4]
    classification = None
    if len(base_targets) >= 2:
        mask = np.zeros(len(meta), dtype=bool)
        for _, m in base_targets:
            mask |= m
        y = meta.loc[mask, "component"].values
        cls = [name for name, _ in base_targets]
        print(f"Clasificación LOGO: {len(cls)} clases, {int(mask.sum())} ventanas")
        classification = classify_logo(X[mask], y, groups[mask], cls,
                                       permutations=args.permutations)
        pred = classification.pop("predictions")
        proba = classification.pop("probabilities")
        (out / "classification.json").write_text(
            json.dumps(classification, indent=2, ensure_ascii=False, default=float))
        pred_df = pd.DataFrame({
            "window_id": meta.loc[mask, "window_id"].values,
            "genome": groups[mask],
            "true_component": y,
            "predicted_component": pred,
        })
        for j, c in enumerate(cls):
            pred_df[f"proba_{c}"] = proba[:, j]
        pred_df.to_csv(out / "classification_predictions.tsv", sep="\t", index=False)

    # Heatmap de firmas (top dims por |d| máximo entre componentes).
    try:
        top_idx = np.argsort(np.abs(cohen.values).max(axis=0))[::-1][: args.top_dims]
        _plot_signatures(cohen.iloc[:, top_idx], out / "heatmap_signatures.png")
    except Exception as exc:
        print(f"  [warn] heatmap: {exc}")

    # Informe.
    _write_report(out, meta, signatures, selected, long, classification, args)

    print(f"\nOK -> {out}")
    print(f"  class_centroids.tsv          {centroids.shape}")
    print(f"  class_cohen_d.tsv            {cohen.shape}")
    print(f"  class_signature_sparse.tsv   {len(signatures)} x {len(dims)}")
    print(f"  class_signature_selected.tsv {len(selected)} dims")
    if classification:
        print(f"  classification.json          balanced_acc="
              f"{classification['balanced_accuracy']:.3f} "
              f"(baseline {classification['baseline_majority_balanced_accuracy']:.3f})")


def _plot_signatures(mat: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    values = mat.values.astype(float)
    fig, ax = plt.subplots(figsize=(max(8, 0.22 * values.shape[1] + 3), 0.7 * values.shape[0] + 3))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-1.0, vmax=1.0,
                   interpolation="nearest")
    ax.set_yticks(range(values.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=9)
    step = max(1, values.shape[1] // 30)
    ax.set_xticks(range(0, values.shape[1], step))
    ax.set_xticklabels([mat.columns[i] for i in range(0, values.shape[1], step)],
                       rotation=90, fontsize=6)
    ax.set_title("Firma de cada componente (d de Cohen vs resto)")
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="d")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_report(out: Path, meta: pd.DataFrame, signatures: list[dict],
                  selected: pd.DataFrame, long: pd.DataFrame,
                  classification: dict | None, args) -> None:
    lines = []
    lines.append("# Firmas multivariadas por componente de ADN\n")
    lines.append(f"Generado: {datetime.now(timezone.utc).isoformat()}\n")
    lines.append(f"- Run: `{args.run}`")
    lines.append(f"- Ventanas: **{len(meta)}**")
    lines.append(f"- Genomas: **{meta['genome'].nunique()}**")
    lines.append(f"- Dimensiones: **{len(long['dim'].unique())}**")
    lines.append(f"- Ajuste de composición: `{args.adjust}`")
    lines.append(f"- Criterio de firma: `q_ventana < 0.01`, `|d| >= {args.d_min}`, "
                 f"consistencia >= {args.consistency_min}. `q_meta` se reporta "
                 f"aparte (la meta por genoma tiene menos potencia que la prueba por ventana).\n")

    lines.append("## Conteo de componentes\n")
    lines.append("| componente | n | genomas | dims firma | dims meta q<0.05 |")
    lines.append("|---|---:|---:|---:|---:|")
    for s in signatures:
        lines.append(f"| {s['name']} | {s['n']} | {s['n_genomes']} | "
                     f"{int(s['selected'].sum())} | {int(s['selected_meta'].sum())} |")
    lines.append("")

    lines.append("## Patrón de cada componente (top dims)\n")
    for s in signatures:
        idx = np.argsort(np.abs(s["d"]))[::-1][:12]
        lines.append(f"### {s['name']}\n")
        lines.append("| dim | media componente | media resto | d | consistencia | q ventana | q meta |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for i in idx:
            lines.append(
                f"| dim_{i:03d} | "
                f"{s['mean_c'][i]:.3f} | {s['mean_r'][i]:.3f} | "
                f"{s['d'][i]:.3f} | {s['consistency'][i]:.2f} | "
                f"{s['q_win'][i]:.2e} | {s['meta_q'][i]:.2e} |"
            )
        lines.append("")

    if classification:
        lines.append("## Clasificación LOGO (leave-one-genome-out)\n")
        lines.append(f"- Clases: {', '.join(classification['classes'])}")
        lines.append(f"- n: {classification['n']}")
        lines.append(f"- Accuracy: **{classification['accuracy']:.3f}**")
        lines.append(f"- Balanced accuracy: **{classification['balanced_accuracy']:.3f}**")
        lines.append(f"- Baseline mayoritario (balanced): "
                     f"{classification['baseline_majority_balanced_accuracy']:.3f}")
        lines.append(f"- AUC one-vs-rest: "
                     f"{ {k: round(v, 3) for k, v in classification['auc_ovr'].items()} }")
        lines.append("\nMatriz de confusión (filas = verdad, columnas = predicción):\n")
        lines.append("| | " + " | ".join(classification["classes"]) + " |")
        lines.append("|---|" + "|".join(["---:"] * len(classification["classes"])) + "|")
        for i, c in enumerate(classification["classes"]):
            row = classification["confusion_matrix"][i]
            lines.append(f"| **{c}** | " + " | ".join(str(x) for x in row) + " |")
        lines.append("")

    lines.append("## Cómo leer la firma\n")
    lines.append("- `class_centroids.tsv` es el valor medio crudo de las 384 dims "
                 "para cada componente. Es el patrón en el espacio del modelo.")
    lines.append("- `class_cohen_d.tsv` es el mismo patrón estandarizado: positivo "
                 "= la dim es mayor en el componente que en el resto; negativo = menor.")
    lines.append("- `class_signature_selected.tsv` es la firma dispersa: solo las "
                 "dims que pasan FDR meta, tamaño de efecto y consistencia entre genomas.")
    lines.append("- `signature_scores.tsv` da el coseno de cada ventana contra la "
                 "firma de cada componente; `predicted_component` es el argmax.")
    lines.append("")
    lines.append("## Límites\n")
    n_gen = meta["genome"].nunique()
    lines.append(f"- Con {n_gen} genomas, la meta-análisis por genoma tiene menos "
                 f"potencia que la prueba por ventana; la consistencia se calcula "
                 f"sobre {n_gen} valores y debe leerse como exploratoria.")
    lines.append("- Las ventanas se solapan y algunas clases (rRNA, CRISPR) tienen "
                 "pocos ejemplos; los p-valores por ventana no corrigen dependencia "
                 "espacial.")
    lines.append("- El checkpoint es el epoch 9 (subentrenado). La firma describe lo "
                 "que el modelo actual codifica, no una verdad biológica.")
    (out / "report_signatures.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
