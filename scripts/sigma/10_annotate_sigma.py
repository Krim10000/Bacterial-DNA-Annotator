#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Anota factores sigma y promotores sigma-dependientes.

Uso:
    python scripts/sigma/10_annotate_sigma.py --genome ecoli_k12
    python scripts/sigma/10_annotate_sigma.py --all
    python scripts/sigma/10_annotate_sigma.py --selftest

Salida:
    data/sigma/<genoma>/sigma_factors.tsv
    data/sigma/<genoma>/sigma_promoters.tsv
    data/sigma/<genoma>/sigma_summary.json
    data/sigma/<genoma>/sigma_report.md

No descarga datos, no usa GPU y no escribe fuera de `data/sigma/`.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.sigma.sigma_annotator import annotate_genome, selftest  # noqa: E402


def _genomes(data_root: Path) -> list[str]:
    reg = data_root / "regulation"
    if not reg.exists():
        return []
    return sorted(d.name for d in reg.iterdir()
                  if d.is_dir() and (d / "genes.tsv").exists()
                  and (d / "operons.tsv").exists())


def _write_batch_summary(data_root: Path, rows: list[dict]) -> Path:
    path = data_root / "sigma" / "summary_all.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["genome", "n_operons", "n_assigned", "n_ambiguous", "n_unresolved",
            "n_factors", "n_factors_countable", "n_factors_unresolved",
            "sospecha_fuerte", "sospecha_media", "sospecha_baja",
            "sin_sospecha", "sha256"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    return path


def _run_one(payload: tuple) -> dict:
    """Worker picklable para el lote paralelo."""
    data_root, genome, out_root, kwargs = payload
    res = annotate_genome(data_root, genome, Path(out_root) / genome, **kwargs)
    s = res["summary"]
    return {
        "genome": genome,
        "n_operons": s["n_operons"],
        "n_assigned": s["n_assigned"],
        "n_ambiguous": s["n_ambiguous"],
        "n_unresolved": s["n_unresolved"],
        "n_factors": s["n_factors"],
        "n_factors_countable": s.get("n_factors_countable", ""),
        "n_factors_unresolved": s.get("n_factors_unresolved", ""),
        "sospecha_fuerte": s["confidence_counts"].get("sospecha_fuerte", 0),
        "sospecha_media": s["confidence_counts"].get("sospecha_media", 0),
        "sospecha_baja": s["confidence_counts"].get("sospecha_baja", 0),
        "sin_sospecha": s["confidence_counts"].get("sin_sospecha", 0),
        "sha256": s["sha256"],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--out-root", type=Path, default=ROOT / "data" / "sigma")
    ap.add_argument("--window", type=int, default=180)
    ap.add_argument("--null-windows", type=int, default=250)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--max-p", type=float, default=0.05)
    ap.add_argument("--min-z", type=float, default=2.0)
    ap.add_argument("--min-margin", type=float, default=1.0)
    ap.add_argument("--no-expression", action="store_true")
    ap.add_argument("--pwm", type=Path, default=None,
                    help="PWM específica (p. ej. data/dbtbs/promoter_pwms_subtilis_all.tsv)")
    ap.add_argument("--policy", choices=("argmax", "two_stage"), default="argmax")
    ap.add_argument("--delta", type=float, default=1.0)
    ap.add_argument("--fallback", default="sigma70")
    ap.add_argument("--kmer-model", type=Path, default=None)
    ap.add_argument("--kmer-threshold", type=float, default=0.6)
    ap.add_argument("--labels", type=Path, default=None)
    ap.add_argument("--label-bonus", type=float, default=0.0)
    ap.add_argument("--jobs", type=int, default=1,
                    help="Procesos en paralelo para --all (por defecto 1).")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        print(selftest())
        return 0

    data_root = ROOT / "data"
    if args.all:
        genomes = _genomes(data_root)
        if args.offset:
            genomes = genomes[args.offset:]
        if args.limit:
            genomes = genomes[:args.limit]
        rows = []
        t0 = time.time()
        kwargs = dict(
            window=args.window, null_windows=args.null_windows,
            seed=args.seed, max_p=args.max_p, min_z=args.min_z,
            min_margin=args.min_margin,
            use_expression=not args.no_expression,
            pwm_path=args.pwm, policy=args.policy,
            delta=args.delta, fallback=args.fallback,
            kmer_model=args.kmer_model, kmer_threshold=args.kmer_threshold,
            labels_path=args.labels, label_bonus=args.label_bonus,
        )
        jobs = max(1, int(args.jobs))
        if jobs == 1:
            for i, genome in enumerate(genomes, 1):
                try:
                    rows.append(_run_one((data_root, genome, args.out_root, kwargs)))
                    s = rows[-1]
                    print(f"[{i}/{len(genomes)}] {genome}: "
                          f"{s['n_assigned']} asignados, {s['n_ambiguous']} ambiguos",
                          flush=True)
                except Exception as exc:  # noqa: BLE001
                    print(f"FALLO {genome}: {exc}", file=sys.stderr)
        else:
            payloads = [(data_root, g, args.out_root, kwargs) for g in genomes]
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                futures = {ex.submit(_run_one, p): p[1] for p in payloads}
                for i, fut in enumerate(as_completed(futures), 1):
                    genome = futures[fut]
                    try:
                        rows.append(fut.result())
                        s = rows[-1]
                        print(f"[{i}/{len(genomes)}] {genome}: "
                              f"{s['n_assigned']} asignados, {s['n_ambiguous']} ambiguos",
                              flush=True)
                    except Exception as exc:  # noqa: BLE001
                        print(f"FALLO {genome}: {exc}", file=sys.stderr)
        rows.sort(key=lambda r: r["genome"])
        path = _write_batch_summary(data_root, rows)
        print(f"lote: {len(rows)}/{len(genomes)} genomas en {time.time()-t0:.1f}s")
        print(f"resumen: {path}")
        return 0 if rows else 1

    if not args.genome:
        ap.error("indica --genome, --all o --selftest")
    res = annotate_genome(
        data_root, args.genome, args.out_root / args.genome,
        window=args.window, null_windows=args.null_windows, seed=args.seed,
        max_p=args.max_p, min_z=args.min_z, min_margin=args.min_margin,
        use_expression=not args.no_expression,
        pwm_path=args.pwm, policy=args.policy,
        delta=args.delta, fallback=args.fallback,
        kmer_model=args.kmer_model, kmer_threshold=args.kmer_threshold,
        labels_path=args.labels, label_bonus=args.label_bonus,
    )
    s = res["summary"]
    print(f"{args.genome}: {s['n_assigned']} asignados · "
          f"{s['n_ambiguous']} ambiguos · {s['n_unresolved']} sin resolver")
    print(f"salida: {res['out_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
