#!/usr/bin/env python3
"""Estabilidad de las firmas por componente entre genomas.

Con 6 genomas una firma puede ser un artefacto de esos 6 linajes. Este script
parte un conjunto grande (20 genomas) en dos mitades al azar, calcula la firma
en cada mitad y mide:

  - Jaccard de las dims seleccionadas
  - correlación de los vectores d completos (384)
  - acuerdo de signo en las dims compartidas
  - dims que sobreviven en >=80% de las particiones con signo estable

También compara la firma de 6 genomas (piloto) con la de 20 genomas.

Salidas:
  signature_stability.tsv          una fila por partición × componente
  signature_stability_summary.tsv  agregado por componente
  signature_stable_dims.tsv        dims estables (frecuencia + signo)
  signature_6g_vs_20g.tsv          Jaccard piloto 6 vs 20
  report.md

Uso:
  python scripts/interpretability/04_signature_stability.py \
      --run data/embeddings/v2_epoch9_20g_256 \
      --pilot data/embeddings/v2_epoch9_pilot6_256 \
      --splits 20
"""
from __future__ import annotations

import argparse
import importlib.util
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

COMPONENTS = ["CDS", "tRNA", "rRNA", "intergenic"]


def load_sig03():
    path = Path(__file__).with_name("03_component_signatures.py")
    spec = importlib.util.spec_from_file_location("sig03", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--pilot", default="", help="Run de 6 genomas para comparar")
    ap.add_argument("--out", default="")
    ap.add_argument("--splits", type=int, default=20)
    ap.add_argument("--d-min", type=float, default=0.25)
    ap.add_argument("--consistency-min", type=float, default=0.80)
    ap.add_argument("--stable-freq", type=float, default=0.80)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--adjust", choices=["none", "gc", "gc_genome"],
                    default="gc_genome")
    args = ap.parse_args()

    sig03 = load_sig03()
    run = Path(args.run)
    out = Path(args.out) if args.out else run
    out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(run / "windows.tsv", sep="\t")
    X_raw = np.load(run / "embeddings.npy").astype(np.float32)
    if len(meta) != X_raw.shape[0]:
        raise SystemExit("Desalineado")
    groups = meta["genome"].values
    gc = meta["gc"].astype(float).values if "gc" in meta.columns else np.zeros(len(meta))
    X = sig03.adjust_matrix(X_raw, gc, groups, args.adjust)

    comp = sig03.assign_components(meta)
    meta["component"] = comp
    genomes = np.array(sorted(pd.unique(groups)))
    print(f"Run: {run}")
    print(f"Ventanas: {len(meta)} | genomas: {len(genomes)} | ajuste: {args.adjust}")

    # Máscaras por componente.
    masks = {}
    for c in COMPONENTS:
        m = (comp == c).values
        if m.sum() >= 20:
            masks[c] = m
    if not masks:
        raise SystemExit("Sin componentes analizables.")

    rng = np.random.default_rng(args.seed)
    rows = []
    # Guardamos por componente: lista de (selected, d) por partición.
    per_component = {c: [] for c in masks}

    for s in range(args.splits):
        perm = rng.permutation(genomes)
        half = len(perm) // 2
        A, B = set(perm[:half]), set(perm[half:])
        idxA = np.isin(groups, list(A))
        idxB = np.isin(groups, list(B))
        for c, mask in masks.items():
            mA, mB = mask[idxA], mask[idxB]
            if mA.sum() < 20 or (~mA).sum() < 20 or mB.sum() < 20 or (~mB).sum() < 20:
                continue
            rA = sig03.signature_stats(X[idxA], mA, groups[idxA], c,
                                       args.d_min, args.consistency_min,
                                       min_genomes=2)
            rB = sig03.signature_stats(X[idxB], mB, groups[idxB], c,
                                       args.d_min, args.consistency_min,
                                       min_genomes=2)
            if rA is None or rB is None:
                continue
            sA, sB = rA["selected"], rB["selected"]
            inter = int((sA & sB).sum())
            union = int((sA | sB).sum())
            jac = inter / union if union else float("nan")
            d_corr = (float(np.corrcoef(rA["d"], rB["d"])[0, 1])
                      if np.std(rA["d"]) > 0 and np.std(rB["d"]) > 0 else np.nan)
            if inter:
                sign_agree = float((np.sign(rA["d"][sA & sB]) ==
                                    np.sign(rB["d"][sA & sB])).mean())
            else:
                sign_agree = np.nan
            rows.append({
                "split": s, "component": c,
                "n_A": int(mA.sum()), "n_B": int(mB.sum()),
                "dims_A": int(sA.sum()), "dims_B": int(sB.sum()),
                "overlap": inter, "union": union, "jaccard": jac,
                "d_correlation": d_corr, "sign_agreement": sign_agree,
            })
            per_component[c].append({"selected": sA & sB,
                                     "selected_A": sA, "selected_B": sB,
                                     "d_A": rA["d"], "d_B": rB["d"],
                                     "q_A": rA["q_win"], "q_B": rB["q_win"]})
        if (s + 1) % 5 == 0:
            print(f"  partición {s + 1}/{args.splits}")

    det = pd.DataFrame(rows)
    if det.empty:
        raise SystemExit("No se pudo calcular ninguna partición.")
    det.to_csv(out / "signature_stability.tsv", sep="\t", index=False)

    summary = (det.groupby("component")
               .agg(splits=("split", "nunique"),
                    jaccard_mean=("jaccard", "mean"),
                    jaccard_median=("jaccard", "median"),
                    d_corr_mean=("d_correlation", "mean"),
                    d_corr_median=("d_correlation", "median"),
                    sign_agree_mean=("sign_agreement", "mean"),
                    overlap_mean=("overlap", "mean"))
               .reset_index()
               .sort_values("jaccard_mean", ascending=False))
    summary.to_csv(out / "signature_stability_summary.tsv", sep="\t", index=False)

    # Dims estables: frecuencia de selección en las mitades y signo.
    stable_rows = []
    for c, records in per_component.items():
        if not records:
            continue
        n = len(records)
        freq = np.zeros(384)
        sign_pos = np.zeros(384)
        d_all = []
        for r in records:
            sel = r["selected_A"] | r["selected_B"]
            freq += sel
            sign_pos += (r["d_A"] > 0) & sel
            d_all.append((r["d_A"] + r["d_B"]) / 2.0)
        freq /= n
        d_mean = np.mean(d_all, axis=0)
        d_sd = np.std(d_all, axis=0)
        # Consistencia de signo solo entre las veces que la dim fue seleccionada.
        count_sel = np.maximum(freq * n, 1.0)
        sign_cons = np.where(d_mean >= 0, sign_pos / count_sel,
                             1.0 - sign_pos / count_sel)
        stable = (freq >= args.stable_freq) & (sign_cons >= 0.90)
        for i in range(384):
            if stable[i]:
                stable_rows.append({
                    "component": c, "dim": f"dim_{i:03d}",
                    "frequency": freq[i], "sign_consistency": sign_cons[i],
                    "d_mean": d_mean[i], "d_sd": d_sd[i],
                })
    stable = pd.DataFrame(stable_rows)
    stable.to_csv(out / "signature_stable_dims.tsv", sep="\t", index=False)

    # Comparación 6 vs 20 genomas.
    comp_rows = []
    if args.pilot:
        pilot = Path(args.pilot)
        p_mask_path = pilot / "class_signature_mask.tsv"
        if p_mask_path.exists():
            p_mask = pd.read_csv(p_mask_path, sep="\t", index_col=0)
            # Recalcular la máscara del run grande con todos los genomas.
            full_masks = {}
            for c, mask in masks.items():
                r = sig03.signature_stats(X, mask, groups, c, args.d_min,
                                          args.consistency_min, min_genomes=4)
                if r is not None:
                    full_masks[c] = r["selected"]
            for c, sel in full_masks.items():
                if c not in p_mask.index:
                    continue
                sel6 = p_mask.loc[c].values.astype(bool)
                inter = int((sel & sel6).sum())
                union = int((sel | sel6).sum())
                comp_rows.append({
                    "component": c,
                    "dims_6g": int(sel6.sum()),
                    "dims_20g": int(sel.sum()),
                    "overlap": inter,
                    "jaccard_6g_vs_20g": inter / union if union else np.nan,
                })
            pd.DataFrame(comp_rows).to_csv(
                out / "signature_6g_vs_20g.tsv", sep="\t", index=False)

    _write_report(out, det, summary, stable, comp_rows, args, len(genomes))
    print(f"\nOK -> {out}")
    print(summary.to_string(index=False))
    print(f"  dims estables: {len(stable)}")
    if comp_rows:
        print(pd.DataFrame(comp_rows).to_string(index=False))


def _write_report(out: Path, det: pd.DataFrame, summary: pd.DataFrame,
                  stable: pd.DataFrame, comp_rows: list[dict],
                  args, n_genomes: int) -> None:
    lines = []
    lines.append("# Estabilidad de firmas por componente\n")
    lines.append(f"Generado: {datetime.now(timezone.utc).isoformat()}\n")
    lines.append(f"- Run: `{args.run}`")
    lines.append(f"- Genomas: **{n_genomes}**")
    lines.append(f"- Particiones: **{args.splits}**")
    lines.append(f"- Ajuste: `{args.adjust}`")
    lines.append(f"- Criterio de dim estable: frecuencia ≥ {args.stable_freq}, "
                 f"consistencia de signo ≥ 0.90\n")

    lines.append("## Resumen por componente\n")
    lines.append("| componente | particiones | Jaccard medio | Jaccard mediana | "
                 "corr(d) media | acuerdo de signo |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for _, r in summary.iterrows():
        lines.append(f"| {r['component']} | {int(r['splits'])} | "
                     f"{r['jaccard_mean']:.3f} | {r['jaccard_median']:.3f} | "
                     f"{r['d_corr_mean']:.3f} | {r['sign_agree_mean']:.3f} |")
    lines.append("")

    lines.append("## Interpretación\n")
    lines.append("- Jaccard mide cuánto se repiten las dims seleccionadas entre dos "
                 "mitades de genomas. 1 = idénticas; 0 = disjuntas.")
    lines.append("- `corr(d)` mide si el patrón completo de 384 dims es el mismo "
                 "en ambas mitades. Es la medida más informativa para un modelo "
                 "nuevo: si es alta, la firma es una propiedad estable del genoma; "
                 "si es baja, el modelo está capturando idiosincrasias de linaje.")
    lines.append("- Una firma con Jaccard bajo pero `corr(d)` alta puede tener "
                 "dims seleccionadas distintas pero el mismo patrón subyacente.")
    lines.append("")

    if comp_rows:
        lines.append("## Piloto 6 vs 20 genomas\n")
        lines.append("| componente | dims 6g | dims 20g | solapamiento | Jaccard |")
        lines.append("|---|---:|---:|---:|---:|")
        for r in comp_rows:
            lines.append(f"| {r['component']} | {r['dims_6g']} | {r['dims_20g']} | "
                         f"{r['overlap']} | {r['jaccard_6g_vs_20g']:.3f} |")
        lines.append("")

    if not stable.empty:
        lines.append("## Dims estables\n")
        lines.append("| componente | n dims estables | ejemplos |")
        lines.append("|---|---:|---|")
        for c, g in stable.groupby("component"):
            ex = ", ".join(g.head(5)["dim"].tolist())
            lines.append(f"| {c} | {len(g)} | {ex} |")
        lines.append("")

    lines.append("## Implicación para el modelo nuevo\n")
    lines.append("- Si tRNA/rRNA/CDS tienen `corr(d)` alta y Jaccard razonable, el "
                 "modelo nuevo debe conservar la capacidad de representar esos "
                 "componentes y no gastar capacidad en reconstruir composición GC.")
    lines.append("- Si intergénico/upstream no son estables, no conviene añadir una "
                 "cabeza de promotor sobre este embedding: primero hay que cambiar "
                 "tokenización/contexto/objetivo.")
    lines.append("- Una `corr(d)` baja entre mitades indica que el modelo actual "
                 "aprende firmas dependientes del linaje; el modelo nuevo necesita "
                 "más diversidad de genomas y regularización/objetivos auxiliares.")
    (out / "report_stability.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
