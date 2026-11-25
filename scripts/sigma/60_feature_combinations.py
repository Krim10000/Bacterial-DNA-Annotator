#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minar combinaciones de features que refinan la asignación de sigma.

No entrena un modelo: evalúa reglas interpretables sobre la tabla
`layer_sigma_features.tsv` (una fila por TSS × sigma) y reporta precisión,
recall, F1 y lift frente al top-1.

Features disponibles:
  z, pvalue, score, spacer_dev, dist_to_start, margin, inventory.
No usa `n_true`, `gold_confidence` ni `source`: son metadatos del gold, no
señales del pipeline.

Uso:
    python scripts/sigma/60_feature_combinations.py \
        --features data/sigma/layers/ecoli_k12/layer_sigma_features.tsv

    python scripts/sigma/60_feature_combinations.py --all
    python scripts/sigma/60_feature_combinations.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Iterable

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

VERSION = "1.0.0"
DEFAULT_MIN_SUPPORT = 5
DEFAULT_TOP = 10

#: Condiciones interpretables. Cada una devuelve True/False.
CONDITIONS: dict[str, tuple[str, Callable[[dict[str, str]], bool]]] = {
    "z>=2": ("z ≥ 2", lambda r: _f(r.get("z")) >= 2.0),
    "z>=3": ("z ≥ 3", lambda r: _f(r.get("z")) >= 3.0),
    "p<=0.05": ("p ≤ 0,05", lambda r: _f(r.get("pvalue"), 1.0) <= 0.05),
    "p<=0.01": ("p ≤ 0,01", lambda r: _f(r.get("pvalue"), 1.0) <= 0.01),
    "m>=1": ("margen ≥ 1", lambda r: _f(r.get("margin")) >= 1.0),
    "m>=2": ("margen ≥ 2", lambda r: _f(r.get("margin")) >= 2.0),
    "spacer_ok": ("spacer en rango", lambda r: _f(r.get("spacer_dev")) == 0.0),
    "dist60": ("distancia ≤ 60 pb", lambda r: 0 <= _f(r.get("dist_to_start"), 9999) <= 60),
    "inventory": ("inventario presente", lambda r: (r.get("inventory") or "").lower() == "yes"),
}


def _f(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_features(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _metrics(rows: list[dict[str, str]], rule: Callable[[dict[str, str]], bool]) -> dict[str, Any]:
    tp = fp = fn = 0
    n_pred = 0
    for r in rows:
        true = (r.get("true") or "").lower() == "yes"
        pred = bool(rule(r))
        if pred:
            n_pred += 1
            if true:
                tp += 1
            else:
                fp += 1
        elif true:
            fn += 1
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall else None)
    return {"n_pred": n_pred, "tp": tp, "fp": fp, "fn": fn,
            "precision": precision, "recall": recall, "f1": f1}


def evaluate_sigma(rows: list[dict[str, str]], sigma: str,
                   min_support: int = DEFAULT_MIN_SUPPORT,
                   top: int = DEFAULT_TOP) -> dict[str, Any]:
    sub = [r for r in rows if (r.get("sigma") or "") == sigma]
    if not sub:
        return {"sigma": sigma, "n": 0, "baseline": None, "rules": []}
    baseline = _metrics(sub, lambda r: (r.get("is_top1") or "").lower() == "yes")
    base_prec = baseline["precision"]
    results = []
    cond_names = list(CONDITIONS)
    for k in (1, 2, 3):
        for combo in itertools.combinations(cond_names, k):
            funcs = [CONDITIONS[c][1] for c in combo]
            rule = lambda r, fs=funcs: all(f(r) for f in fs)
            m = _metrics(sub, rule)
            if m["tp"] < min_support or m["n_pred"] < min_support:
                continue
            lift = (m["precision"] / base_prec) if base_prec and m["precision"] else None
            results.append({
                "sigma": sigma,
                "rule": " AND ".join(CONDITIONS[c][0] for c in combo),
                "keys": list(combo),
                "n_pred": m["n_pred"], "tp": m["tp"], "fp": m["fp"],
                "fn": m["fn"], "precision": m["precision"],
                "recall": m["recall"], "f1": m["f1"], "lift": lift,
            })
    results.sort(key=lambda r: (r["f1"] if r["f1"] is not None else -1,
                                r["precision"] if r["precision"] is not None else -1,
                                r["tp"]), reverse=True)
    # Quitar reglas idénticas en métricas, conservar la más simple.
    seen = set()
    unique = []
    for r in results:
        key = (r["n_pred"], r["tp"], r["fp"], r["fn"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(r)
    return {"sigma": sigma, "n": len(sub), "baseline": baseline,
            "rules": unique[:top]}


def evaluate_file(path: Path, min_support: int = DEFAULT_MIN_SUPPORT,
                  top: int = DEFAULT_TOP,
                  sigma_filter: Iterable[str] | None = None) -> dict[str, Any]:
    rows = read_features(path)
    genome = path.parent.name
    sigmas = sorted({r.get("sigma") or "" for r in rows if r.get("sigma")})
    if sigma_filter:
        wanted = set(sigma_filter)
        sigmas = [s for s in sigmas if s in wanted]
    return {
        "version": VERSION,
        "genome": genome,
        "input": str(path),
        "n_rows": len(rows),
        "n_tss": len({(r.get("contig"), r.get("tss"), r.get("strand")) for r in rows}),
        "per_sigma": [evaluate_sigma(rows, s, min_support, top) for s in sigmas],
    }


def _write_tsv(path: Path, results: list[dict[str, Any]]) -> None:
    cols = ["genome", "sigma", "rule", "n_pred", "tp", "fp", "fn",
            "precision", "recall", "f1", "lift",
            "baseline_precision", "baseline_recall", "baseline_f1"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for res in results:
            for sigma_res in res.get("per_sigma", []):
                base = sigma_res.get("baseline") or {}
                for rule in sigma_res.get("rules", []):
                    vals = {
                        "genome": res.get("genome", ""),
                        "sigma": sigma_res.get("sigma", ""),
                        "rule": rule.get("rule", ""),
                        "n_pred": rule.get("n_pred", ""),
                        "tp": rule.get("tp", ""), "fp": rule.get("fp", ""),
                        "fn": rule.get("fn", ""),
                        "precision": rule.get("precision", ""),
                        "recall": rule.get("recall", ""),
                        "f1": rule.get("f1", ""), "lift": rule.get("lift", ""),
                        "baseline_precision": base.get("precision", ""),
                        "baseline_recall": base.get("recall", ""),
                        "baseline_f1": base.get("f1", ""),
                    }
                    fh.write("\t".join(str(vals.get(c, "")) for c in cols) + "\n")


def _report(results: list[dict[str, Any]]) -> str:
    lines = ["# Combinaciones de features para sigma", "",
             "Reglas interpretables sobre `layer_sigma_features.tsv`. "
             "No es un modelo entrenado; es un análisis descriptivo.", ""]
    for res in results:
        lines += [f"## {res['genome']}", "",
                  f"TSS: {res['n_tss']} · filas TSS×sigma: {res['n_rows']}", ""]
        for s in res["per_sigma"]:
            base = s.get("baseline") or {}
            if not base:
                continue
            lines += [f"### {s['sigma']} (n={s['n']})", "",
                      f"Baseline top-1: precision {base.get('precision')}, "
                      f"recall {base.get('recall')}, F1 {base.get('f1')}", "",
                      "| regla | pred | TP | FP | precision | recall | F1 | lift |",
                      "|---|---:|---:|---:|---:|---:|---:|---:|"]
            for r in s.get("rules", [])[:5]:
                def f(x): return f"{x:.3f}" if isinstance(x, float) else "—"
                lines.append(f"| {r['rule']} | {r['n_pred']} | {r['tp']} | {r['fp']} | "
                             f"{f(r['precision'])} | {f(r['recall'])} | {f(r['f1'])} | "
                             f"{f(r['lift'])} |")
            lines.append("")
    lines += ["## Lectura", "",
              "Una regla con buen F1 y lift alto puede refinar un caso, pero si se",
              "elige sobre el mismo gold no es validación. Hay que probarla en",
              "held-out. Las reglas con `inventory` dependen de que exista el",
              "inventario de sigma del genoma.", ""]
    return "\n".join(lines)


def selftest() -> str:
    rows = [
        {"sigma": "sigma70", "true": "yes", "is_top1": "yes", "z": "4",
         "pvalue": "0.001", "margin": "3", "spacer_dev": "0",
         "dist_to_start": "30", "inventory": "yes"},
        {"sigma": "sigma70", "true": "no", "is_top1": "no", "z": "1",
         "pvalue": "0.2", "margin": "0", "spacer_dev": "2",
         "dist_to_start": "100", "inventory": "no"},
    ]
    res = evaluate_sigma(rows, "sigma70", min_support=1, top=3)
    assert res["baseline"]["tp"] == 1
    assert res["rules"]
    assert any("z ≥" in r["rule"] for r in res["rules"])
    return "métricas y combinaciones sigma OK"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features", type=Path,
                    help="Un layer_sigma_features.tsv.")
    ap.add_argument("--all", action="store_true",
                    help="Todos los data/sigma/layers/*/layer_sigma_features.tsv.")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "sigma" / "feature_combinations")
    ap.add_argument("--sigma", nargs="*", default=None,
                    help="Filtra por sigma (p. ej. sigma38 sigma54).")
    ap.add_argument("--min-support", type=int, default=DEFAULT_MIN_SUPPORT)
    ap.add_argument("--top", type=int, default=DEFAULT_TOP)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        print(selftest())
        return 0

    paths: list[Path] = []
    if args.all:
        paths = sorted((ROOT / "data" / "sigma" / "layers").glob(
            "*/layer_sigma_features.tsv"))
    elif args.features:
        paths = [args.features]
    if not paths:
        print("No hay features. Corre primero 50_eval_layers.py.", file=sys.stderr)
        return 2
    all_results = []
    for path in paths:
        res = evaluate_file(path, args.min_support, args.top, args.sigma)
        all_results.append(res)
        print(f"{res['genome']}: {res['n_rows']} filas, {res['n_tss']} TSS")
    _write_tsv(args.out / "feature_combinations.tsv", all_results)
    (args.out / "feature_combinations.json").write_text(
        json.dumps(all_results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (args.out / "feature_combinations.md").write_text(_report(all_results),
                                                       encoding="utf-8")
    print(f"salida: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
