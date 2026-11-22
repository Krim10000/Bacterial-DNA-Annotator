#!/usr/bin/env python3
"""Correlaciona las 384 dimensiones del embedding con anotaciones de ADN.

Entrada: salida de 01_extract_embeddings.py (windows.tsv + embeddings.npy).

Para cada (dimensión, anotación) calcula:
  - r de Pearson y su p-valor
  - rho de Spearman
  - r parcial controlando GC (y GC^2), para separar "composición" de "anotación"
  - d de Cohen y AUC de una sola dimensión (descriptivos, no el objetivo)
  - q-valor FDR (Benjamini-Hochberg) dentro de cada anotación
  - p de permutación del estadístico máximo (opcional)

Dos alcances:
  - `global`: todas las ventanas juntas. Mezcla diferencias entre genomas.
  - `within_genome`: centra cada variable y cada dimensión por genoma (efectos
    fijos). Responde "dentro de un mismo genoma, ¿qué dimensión distingue esta
    anotación?", sin que la identidad del genoma explique la correlación.

Salidas:
  correlations.tsv                    tabla larga (scope x dim x anotación)
  correlations_raw_r_<scope>.tsv      matriz 384 x anotaciones (r de Pearson)
  correlations_partial_r_<scope>.tsv  matriz 384 x anotaciones (r parcial | GC)
  top_correlations.tsv                dims con q<0.05 y |r_parcial|>=0.2
  annotation_summary.tsv              resumen por anotación y scope
  heatmap_partial_r_<scope>.png       mapa de calor 384 x anotaciones
  report.md                           informe legible

Uso:
  python scripts/interpretability/02_correlate_annotations.py \
      --run data/embeddings/v2_epoch9_pilot --scope both
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
MPL_DIR = ROOT / ".cache" / "matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))

BINARY = [
    "is_cds", "is_trna", "is_rrna", "is_ncrna", "is_intergenic",
    "is_upstream", "is_oric", "is_crispr", "is_prophage", "is_transposon",
]
CONTINUOUS = [
    "gc", "gc_skew", "at_skew", "entropy", "n_frac",
    "overlap_cds", "overlap_trna", "overlap_rrna", "overlap_ncrna",
]
COMPOSITION = {"gc", "gc_skew", "at_skew", "entropy", "n_frac",
               "log_dist_to_feature"}


def _family(name: str) -> str:
    if name in COMPOSITION:
        return "composition"
    if name == "split_excluded":
        return "genome"
    if name.startswith("overlap_"):
        return "overlap"
    if name.startswith("is_") or name == "strand_plus":
        return "class"
    return "other"


# ---------------------------------------------------------------------------
# Estadística vectorizada (384 columnas a la vez)
# ---------------------------------------------------------------------------
def _pearson_cols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    Xc = X - X.mean(axis=0, keepdims=True)
    yc = y - y.mean()
    denom = np.sqrt((Xc * Xc).sum(axis=0) * (yc * yc).sum())
    num = (Xc * yc[:, None]).sum(axis=0)
    out = np.zeros(X.shape[1], dtype=np.float64)
    np.divide(num, denom, out=out, where=denom > 0)
    return out


def _p_from_r(r: np.ndarray, df: int) -> np.ndarray:
    if df <= 0:
        return np.full_like(r, np.nan, dtype=np.float64)
    r = np.clip(r, -0.999999, 0.999999)
    t = r * np.sqrt(df / (1.0 - r * r))
    return 2.0 * stats.t.sf(np.abs(t), df=df)


def _cohen_d(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    a, b = X[y == 1], X[y == 0]
    if len(a) < 2 or len(b) < 2:
        return np.full(X.shape[1], np.nan)
    ma, mb = a.mean(axis=0), b.mean(axis=0)
    va, vb = a.var(axis=0, ddof=1), b.var(axis=0, ddof=1)
    pooled = np.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2))
    out = np.zeros(X.shape[1], dtype=np.float64)
    np.divide(ma - mb, pooled, out=out, where=pooled > 0)
    return out


def _auc_cols(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.full(X.shape[1], np.nan)
    ranks = stats.rankdata(X, axis=0)
    return (ranks[y == 1].sum(axis=0) - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    q = np.empty(n, dtype=float)
    q[order] = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out[ok] = np.clip(q, 0.0, 1.0)
    return out


def _residualize(Y: np.ndarray, C: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(C, Y, rcond=None)
    return Y - C @ beta


def _center_by_group(A: np.ndarray, groups: np.ndarray | None) -> np.ndarray:
    A = np.asarray(A, dtype=np.float64).copy()
    if groups is None:
        return A
    for g in np.unique(groups):
        idx = groups == g
        if idx.sum() == 0:
            continue
        if A.ndim == 1:
            A[idx] -= A[idx].mean()
        else:
            A[idx] -= A[idx].mean(axis=0, keepdims=True)
    return A


def _center_by_group_masked(y: np.ndarray, groups: np.ndarray | None) -> np.ndarray:
    y = np.asarray(y, dtype=np.float64).copy()
    if groups is None:
        return y
    for g in np.unique(groups):
        idx = (groups == g) & np.isfinite(y)
        if idx.sum():
            y[idx] -= y[idx].mean()
    return y


def _perm_max_r(X: np.ndarray, y: np.ndarray, groups: np.ndarray | None,
                n_perm: int, seed: int) -> np.ndarray:
    """Nulo del máximo |r| entre las 384 dims.

    Si `groups` no es None, permuta dentro de cada genoma (nulo correcto para
    el alcance within-genome). Si es None, permuta globalmente.
    """
    if n_perm <= 0:
        return np.array([])
    rng = np.random.default_rng(seed)
    Xc = X - X.mean(axis=0, keepdims=True)
    out = np.empty(n_perm, dtype=np.float64)
    for p in range(n_perm):
        yp = y.copy()
        if groups is not None:
            for g in np.unique(groups):
                idx = np.flatnonzero(groups == g)
                if len(idx) > 1:
                    yp[idx] = rng.permutation(yp[idx])
        else:
            yp = rng.permutation(yp)
        yc = yp - yp.mean()
        denom = np.sqrt((Xc * Xc).sum(axis=0) * (yc * yc).sum())
        r = np.zeros(X.shape[1], dtype=np.float64)
        np.divide((Xc * yc[:, None]).sum(axis=0), denom, out=r, where=denom > 0)
        out[p] = np.nanmax(np.abs(r))
    return out


def analyze_variable(Xuse, y_raw, groups, name, kind, scope, n_perm, seed,
                     C, min_class, n_groups):
    """Devuelve vectores de 384 valores o None si la variable no es analizable."""
    y = _center_by_group_masked(y_raw, groups) if scope == "within_genome" else y_raw
    mask = np.isfinite(y_raw) & np.isfinite(y)
    if mask.sum() < 20:
        return None

    Xv = Xuse[mask]
    yv = y[mask]
    yraw = y_raw[mask]
    gv = groups[mask] if groups is not None else None
    n = len(yv)
    if np.nanstd(yv) < 1e-12:
        return None

    if kind == "binary":
        yb = (yraw > 0.5).astype(float)
        n_pos, n_neg = int(yb.sum()), int((1 - yb).sum())
        if min(n_pos, n_neg) < min_class:
            return None
    else:
        n_pos = n_neg = -1

    # Covariables de la correlación parcial: gc y gc^2 (2 en ambos alcances).
    k_cov = 2
    df = n - (n_groups if scope == "within_genome" else 0) - 2 - k_cov
    r = _pearson_cols(Xv, yv)
    p = _p_from_r(r, df)
    ranks_X = stats.rankdata(Xv, axis=0)
    ranks_y = stats.rankdata(yv)
    rho = _pearson_cols(ranks_X, ranks_y)
    p_rho = _p_from_r(rho, n - 2)

    Xr = _residualize(Xv, C[mask])
    yr = _residualize(yv, C[mask])
    r_part = _pearson_cols(Xr, yr)
    p_part = _p_from_r(r_part, df)

    if kind == "binary":
        d = _cohen_d(Xv, (yraw > 0.5).astype(float))
        auc = _auc_cols(Xv, (yraw > 0.5).astype(float))
    else:
        d = np.full(Xv.shape[1], np.nan)
        auc = np.full(Xv.shape[1], np.nan)

    perm_p = np.full(Xv.shape[1], np.nan)
    if n_perm > 0:
        null = _perm_max_r(Xv, yv, gv if scope == "within_genome" else None,
                           n_perm, seed)
        if len(null):
            perm_p = (1.0 + (null[:, None] >= np.abs(r)[None, :]).sum(axis=0)) / (1.0 + len(null))

    return {
        "name": name, "kind": kind, "scope": scope,
        "n": n, "n_pos": n_pos, "n_neg": n_neg, "df": df,
        "r": r, "p": p, "rho": rho, "p_rho": p_rho,
        "r_part": r_part, "p_part": p_part, "d": d, "auc": auc,
        "perm_p": perm_p,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="Directorio con windows.tsv y embeddings.npy")
    ap.add_argument("--windows", default="",
                    help="TSV de metadatos alternativo (mismo orden que embeddings.npy)")
    ap.add_argument("--out", default="", help="Directorio de salida (default: el mismo --run)")
    ap.add_argument("--scope", choices=["global", "within_genome", "both"], default="both")
    ap.add_argument("--permutations", type=int, default=100,
                    help="Permutaciones para el estadístico máximo (0 = desactivar)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--min-class", type=int, default=10)
    ap.add_argument("--top-r", type=float, default=0.2,
                    help="|r parcial| mínimo para top_correlations.tsv")
    ap.add_argument("--no-heatmap", action="store_true")
    args = ap.parse_args()

    run = Path(args.run)
    out = Path(args.out) if args.out else run
    out.mkdir(parents=True, exist_ok=True)

    meta_path = Path(args.windows) if args.windows else run / "windows.tsv"
    meta = pd.read_csv(meta_path, sep="\t")
    emb = np.load(run / "embeddings.npy")
    if len(meta) != emb.shape[0]:
        raise SystemExit(f"Desalineado: {len(meta)} filas vs {emb.shape[0]} embeddings")
    if emb.shape[1] != 384:
        raise SystemExit(f"Se esperaban 384 dims, hay {emb.shape[1]}")
    X = emb.astype(np.float32)
    dims = [f"dim_{i:03d}" for i in range(X.shape[1])]

    if "strand" in meta.columns:
        meta["strand_plus"] = (meta["strand"] == 1).astype(int)
    if "split" in meta.columns:
        meta["split_excluded"] = (meta["split"] == "excluded").astype(int)
    if "dist_to_feature" in meta.columns:
        d = meta["dist_to_feature"].astype(float)
        meta["log_dist_to_feature"] = np.where(d >= 0, np.log1p(d), np.nan)

    variables: list[tuple[str, str]] = []
    for c in BINARY:
        if c in meta.columns and meta[c].nunique(dropna=True) > 1:
            variables.append((c, "binary"))
    for c in CONTINUOUS:
        if c in meta.columns and meta[c].nunique(dropna=True) > 1:
            variables.append((c, "continuous"))
    for c in ("strand_plus", "split_excluded", "log_dist_to_feature"):
        if c in meta.columns and meta[c].nunique(dropna=True) > 1:
            kind = "binary" if c in {"strand_plus", "split_excluded"} else "continuous"
            variables.append((c, kind))

    groups = meta["genome"].values if "genome" in meta.columns else None
    n_groups = len(np.unique(groups)) if groups is not None else 0
    gc = meta["gc"].astype(float).values if "gc" in meta.columns else np.zeros(len(meta))

    X_within = _center_by_group(X, groups) if groups is not None else X
    gc_w = _center_by_group(gc, groups) if groups is not None else gc
    C_global = np.column_stack([np.ones(len(meta)), gc, gc * gc])
    C_within = np.column_stack([gc_w, gc_w * gc_w])

    scopes = ["global", "within_genome"] if args.scope == "both" else [args.scope]
    print(f"Ventanas: {len(meta)} | genomas: {meta['genome'].nunique()} | dims: {X.shape[1]}")
    print(f"Alcances: {scopes} | variables: {len(variables)}")

    all_results = []
    for scope in scopes:
        Xuse = X if scope == "global" else X_within
        C = C_global if scope == "global" else C_within
        print(f"\n[{scope}]")
        for name, kind in variables:
            y = meta[name].astype(float).values
            res = analyze_variable(Xuse, y, groups, name, kind, scope,
                                   args.permutations, args.seed, C,
                                   args.min_class, n_groups)
            if res is None:
                print(f"  [skip] {name}: clases insuficientes o constante")
                continue
            all_results.append(res)
            n_sig = int((_bh(res["p_part"]) < 0.05).sum())
            print(f"  {name:22s} n={res['n']:5d} "
                  f"max|r|={np.nanmax(np.abs(res['r'])):.3f} "
                  f"max|r_gc|={np.nanmax(np.abs(res['r_part'])):.3f} "
                  f"dims q<0.05={n_sig}")

    if not all_results:
        raise SystemExit("No hay variables analizables (clases insuficientes).")

    # Tabla larga.
    rows = []
    for res in all_results:
        q = _bh(res["p"])
        q_part = _bh(res["p_part"])
        for i, dim in enumerate(dims):
            rows.append({
                "scope": res["scope"], "dim": dim, "dim_index": i,
                "annotation": res["name"], "kind": res["kind"],
                "family": _family(res["name"]),
                "n": res["n"], "n_pos": res["n_pos"], "n_neg": res["n_neg"],
                "df": res["df"],
                "pearson_r": res["r"][i], "p_pearson": res["p"][i],
                "spearman_r": res["rho"][i], "p_spearman": res["p_rho"][i],
                "partial_r_gc": res["r_part"][i], "p_partial": res["p_part"][i],
                "cohen_d": res["d"][i], "auc_single_dim": res["auc"][i],
                "auc_strength": (max(res["auc"][i], 1.0 - res["auc"][i])
                                 if np.isfinite(res["auc"][i]) else np.nan),
                "q_pearson": q[i], "q_partial": q_part[i],
                "perm_p_max": res["perm_p"][i],
            })
    corr = pd.DataFrame(rows)

    # Medias por clase para anotaciones binarias (descriptivo) y tabla de
    # efecto por dimensión.
    means = {}
    for res in all_results:
        if res["kind"] != "binary":
            continue
        raw = meta[res["name"]].astype(float).values
        ok = np.isfinite(raw)
        yb = raw[ok] > 0.5
        if yb.sum() == 0 or (~yb).sum() == 0:
            continue
        means[res["name"]] = (X[ok][yb].mean(axis=0), X[ok][~yb].mean(axis=0))
    for name, (mp, mn) in means.items():
        for scope in scopes:
            sel = (corr["scope"] == scope) & (corr["annotation"] == name)
            if not sel.any():
                continue
            corr.loc[sel, "mean_pos"] = mp
            corr.loc[sel, "mean_neg"] = mn
            corr.loc[sel, "delta_mean"] = mp - mn
    corr.to_csv(out / "correlations.tsv", sep="\t", index=False)

    # Matrices por alcance.
    part_by_scope = {}
    raw_by_scope = {}
    for scope in scopes:
        raw = {r["name"]: r["r"] for r in all_results if r["scope"] == scope}
        part = {r["name"]: r["r_part"] for r in all_results if r["scope"] == scope}
        if not raw:
            continue
        pd.DataFrame(raw, index=dims).to_csv(
            out / f"correlations_raw_r_{scope}.tsv", sep="\t")
        pd.DataFrame(part, index=dims).to_csv(
            out / f"correlations_partial_r_{scope}.tsv", sep="\t")
        raw_by_scope[scope] = raw
        part_by_scope[scope] = part

    # Resumen.
    summary_rows = []
    for res in all_results:
        q = _bh(res["p"])
        q_part = _bh(res["p_part"])
        rr = np.abs(res["r"])
        rp = np.abs(res["r_part"])
        i_raw = int(np.nanargmax(rr))
        i_part = int(np.nanargmax(rp))
        summary_rows.append({
            "scope": res["scope"], "annotation": res["name"], "kind": res["kind"],
            "family": _family(res["name"]),
            "n": res["n"], "n_pos": res["n_pos"], "n_neg": res["n_neg"],
            "dims_q05_raw": int((q < 0.05).sum()),
            "dims_q05_partial": int((q_part < 0.05).sum()),
            "max_abs_r": float(rr[i_raw]), "top_dim_raw": dims[i_raw],
            "max_abs_r_partial": float(rp[i_part]),
            "top_dim_partial": dims[i_part],
            "top_r_partial_signed": float(res["r_part"][i_part]),
            "top_q_partial": float(q_part[i_part]),
        })
    summary = pd.DataFrame(summary_rows).sort_values(
        ["scope", "max_abs_r_partial"], ascending=[True, False])
    summary.to_csv(out / "annotation_summary.tsv", sep="\t", index=False)

    # Tabla por dimensión: qué anotación explica mejor cada una.
    if "within_genome" in scopes:
        w = corr[corr["scope"] == "within_genome"].copy()
        w["abs_partial"] = w["partial_r_gc"].abs()
        w = w.dropna(subset=["abs_partial"])
        if not w.empty:
            idx = w.groupby("dim")["abs_partial"].idxmax()
            best = w.loc[idx].copy()
            best = best.sort_values("partial_r_gc", key=lambda s: s.abs(),
                                    ascending=False)
            # Estadísticos de la dimensión a través de todas las ventanas.
            stats_dim = pd.DataFrame({
                "dim": dims,
                "emb_mean": X.mean(axis=0),
                "emb_std": X.std(axis=0),
                "emb_min": X.min(axis=0),
                "emb_max": X.max(axis=0),
            })
            best = best.merge(stats_dim, on="dim", how="left")
            cols = ["dim", "emb_mean", "emb_std", "emb_min", "emb_max",
                    "annotation", "family", "kind", "pearson_r",
                    "partial_r_gc", "q_partial", "cohen_d", "auc_single_dim",
                    "auc_strength", "mean_pos", "mean_neg", "delta_mean"]
            best[[c for c in cols if c in best.columns]].to_csv(
                out / "dimension_summary.tsv", sep="\t", index=False)
            # Mejor anotación biológica por dimensión (sin composición ni
            # genoma), elegida solo entre las biológicas.
            bio_all = w[w["family"] == "class"].copy()
            if not bio_all.empty:
                idx_bio = bio_all.groupby("dim")["abs_partial"].idxmax()
                bio_best = bio_all.loc[idx_bio].copy()
                bio_best = bio_best.merge(stats_dim, on="dim", how="left")
                bio_best = bio_best.sort_values(
                    "partial_r_gc", key=lambda s: s.abs(), ascending=False)
                bio_best[[c for c in cols if c in bio_best.columns]].to_csv(
                    out / "dimension_summary_bio.tsv", sep="\t", index=False)
            ov_all = w[w["family"] == "overlap"].copy()
            if not ov_all.empty:
                idx_ov = ov_all.groupby("dim")["abs_partial"].idxmax()
                ov_best = ov_all.loc[idx_ov].copy()
                ov_best = ov_best.merge(stats_dim, on="dim", how="left")
                ov_best = ov_best.sort_values(
                    "partial_r_gc", key=lambda s: s.abs(), ascending=False)
                ov_best[[c for c in cols if c in ov_best.columns]].to_csv(
                    out / "dimension_summary_overlap.tsv", sep="\t", index=False)

    # Top correlaciones.
    top = corr[(corr["q_partial"] < 0.05) &
               (corr["partial_r_gc"].abs() >= args.top_r)].copy()
    top = top.sort_values(["scope", "partial_r_gc"],
                          key=lambda s: s.abs() if s.name == "partial_r_gc" else s,
                          ascending=[True, False])
    top.to_csv(out / "top_correlations.tsv", sep="\t", index=False)
    top_bio = top[top["family"] == "class"].copy()
    top_bio.to_csv(out / "top_correlations_bio.tsv", sep="\t", index=False)

    # Heatmaps.
    if not args.no_heatmap:
        for scope, part in part_by_scope.items():
            mat = pd.DataFrame(part, index=dims)
            try:
                path = out / f"heatmap_partial_r_{scope}.png"
                _plot_heatmap(mat, path, vlim=0.6)
                print(f"  heatmap -> {path}")
            except Exception as exc:
                print(f"  [warn] heatmap {scope}: {exc}")
            bio_cols = [c for c in mat.columns
                        if c not in COMPOSITION and c != "split_excluded"]
            if bio_cols:
                try:
                    path = out / f"heatmap_partial_r_bio_{scope}.png"
                    _plot_heatmap(mat[bio_cols], path, vlim=0.35)
                    print(f"  heatmap -> {path}")
                except Exception as exc:
                    print(f"  [warn] heatmap bio {scope}: {exc}")

    _write_report(out, meta, summary, top, corr, run, args)

    print(f"\nOK -> {out}")
    print(f"  correlations.tsv            {len(corr)} filas")
    print(f"  top_correlations.tsv        {len(top)} filas")
    print(f"  top_correlations_bio.tsv    {len(top_bio)} filas")
    print(f"  annotation_summary.tsv      {len(summary)} filas")
    dim_summary_path = out / "dimension_summary.tsv"
    if dim_summary_path.exists():
        n_dims = len(pd.read_csv(dim_summary_path, sep="\t"))
        print(f"  dimension_summary.tsv       {n_dims} filas")
    bio_path = out / "dimension_summary_bio.tsv"
    if bio_path.exists():
        n_bio = len(pd.read_csv(bio_path, sep="\t"))
        print(f"  dimension_summary_bio.tsv   {n_bio} filas")


def _plot_heatmap(mat: pd.DataFrame, path: Path, vlim: float = 0.6) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, leaves_list

    values = mat.fillna(0.0).values.astype(float)
    if values.shape[0] > 2:
        Z = linkage(values, method="average", metric="euclidean")
        order = leaves_list(Z)
        values = values[order]
        ylabels = [mat.index[i] for i in order]
    else:
        ylabels = list(mat.index)

    fig, ax = plt.subplots(figsize=(max(6, 0.55 * mat.shape[1] + 3), 10))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vlim, vmax=vlim,
                   interpolation="nearest")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=8)
    step = max(1, len(ylabels) // 24)
    ax.set_yticks(range(0, len(ylabels), step))
    ax.set_yticklabels([ylabels[i] for i in range(0, len(ylabels), step)], fontsize=6)
    ax.set_ylabel("dimensión del embedding (384)")
    ax.set_title(f"r parcial (controlando GC) · {path.stem.replace('heatmap_partial_r_', '')}")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="r parcial")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_report(out: Path, meta: pd.DataFrame, summary: pd.DataFrame,
                  top: pd.DataFrame, corr: pd.DataFrame, run: Path,
                  args) -> None:
    lines = []
    lines.append("# Embeddings v2 vs anotaciones de ADN\n")
    lines.append(f"Generado: {datetime.now(timezone.utc).isoformat()}\n")
    lines.append(f"- Run: `{run}`")
    lines.append(f"- Ventanas: **{len(meta)}**")
    lines.append(f"- Genomas: **{meta['genome'].nunique()}**")
    if "window_bp" in meta.columns:
        lines.append(f"- Ventana: **{int(meta['window_bp'].iloc[0])} pb**")
    lines.append("- Dimensiones: **384**")
    lines.append(f"- Alcances: **{args.scope}**")
    lines.append(f"- Permutaciones: **{args.permutations}**\n")

    if "sample_class" in meta.columns:
        counts = meta["sample_class"].str.split("|").explode().value_counts()
        lines.append("## Composición de la muestra\n")
        lines.append("| clase de muestreo | n |")
        lines.append("|---|---:|")
        for k, v in counts.items():
            lines.append(f"| {k} | {v} |")
        lines.append("")

    for scope in summary["scope"].unique():
        sub = summary[summary["scope"] == scope]
        lines.append(f"## Resumen por anotación — `{scope}`\n")
        cols = ["annotation", "kind", "n", "n_pos", "n_neg",
                "dims_q05_raw", "dims_q05_partial",
                "max_abs_r", "top_dim_raw", "max_abs_r_partial", "top_dim_partial"]
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("|" + "|".join(["---"] * len(cols)) + "|")
        for _, r in sub.iterrows():
            vals = []
            for c in cols:
                v = r[c]
                vals.append(f"{v:.3f}" if isinstance(v, float) else str(v))
            lines.append("| " + " | ".join(vals) + " |")
        lines.append("")

    lines.append("## Top correlaciones (q parcial < 0.05)\n")
    if top.empty:
        lines.append("No hay dimensiones que superen el umbral tras corregir por GC. "
                     "Esto es un resultado, no un fallo: el embedding actual no "
                     "contiene señal lineal específica de estas anotaciones.\n")
    else:
        lines.append("| scope | dim | anotación | r | r parcial | q parcial | d | AUC 1D |")
        lines.append("|---|---|---|---:|---:|---:|---:|---:|")
        for _, r in top.head(30).iterrows():
            d = "" if pd.isna(r["cohen_d"]) else f"{r['cohen_d']:.2f}"
            auc = "" if pd.isna(r.get("auc_strength")) else f"{r['auc_strength']:.3f}"
            lines.append(
                f"| {r['scope']} | {r['dim']} | {r['annotation']} | "
                f"{r['pearson_r']:.3f} | {r['partial_r_gc']:.3f} | "
                f"{r['q_partial']:.2e} | {d} | {auc} |"
            )
        lines.append("")

    lines.append("## Top correlaciones biológicas (clases discretas)\n")
    bio = top[top["family"] == "class"].copy() if "family" in top.columns else top.iloc[0:0]
    if bio.empty:
        lines.append("Ninguna clase discreta (CDS/tRNA/rRNA/...) supera el umbral. "
                     "El embedding actual está dominado por composición.\n")
    else:
        lines.append("| scope | dim | anotación | r parcial | q parcial | d | AUC 1D |")
        lines.append("|---|---|---|---:|---:|---:|---:|")
        for _, r in bio.head(25).iterrows():
            d = "" if pd.isna(r["cohen_d"]) else f"{r['cohen_d']:.2f}"
            auc = "" if pd.isna(r.get("auc_strength")) else f"{r['auc_strength']:.3f}"
            lines.append(
                f"| {r['scope']} | {r['dim']} | {r['annotation']} | "
                f"{r['partial_r_gc']:.3f} | {r['q_partial']:.2e} | {d} | {auc} |"
            )
        lines.append("")

    lines.append("## Cómo leer esto\n")
    lines.append("- `global` mezcla diferencias entre genomas. Si un genoma tiene "
                 "un GC o un estilo de vida distinto, muchas dimensiones lo "
                 "capturan y cualquier anotación sesgada por genoma parecerá "
                 "correlacionada.")
    lines.append("- `within_genome` centra cada variable y cada dimensión por "
                 "genoma: elimina la identidad del genoma como explicación. Es "
                 "el alcance que hay que mirar para interpretar biología.")
    lines.append("- `r` es la correlación bruta; `r parcial` elimina GC y GC². "
                 "Si `r` es alta pero `r parcial` cae a ~0, la dimensión codifica "
                 "composición GC, no la anotación.")
    lines.append("- `q` es FDR de Benjamini-Hochberg dentro de cada anotación "
                 "(384 tests). `perm_p_max` es el p de permutación del máximo |r| "
                 "entre las 384 dimensiones (más conservador).")
    lines.append("- `AUC 1D` es descriptivo: cuánto separa esa única dimensión. "
                 "No es el AUROC del modelo ni una validación biológica.")
    lines.append("- Las ventanas pueden solaparse y varias provienen del mismo "
                 "genoma; los p-valores no corrigen dependencia espacial. Para "
                 "inferencia formal, usar split por genoma y ventanas no "
                 "solapadas.")
    lines.append("")

    lines.append("## Límites conocidos del checkpoint\n")
    lines.append("Este análisis describe lo que el checkpoint actual (epoch 9) "
                 "codifica; no demuestra que sea biológicamente correcto. "
                 "El checkpoint está subentrenado (perplejidad ~474 en el genoma "
                 "de entrenamiento, peor que un unigram) y fue entrenado con "
                 "32k nucleótidos (~7.6k tokens), no 32k tokens.\n")

    (out / "report_correlations.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
