#!/usr/bin/env python3
"""Orquestador del pipeline completo para UN genoma.

Encadena las etapas ya existentes y deja un manifiesto con tiempos, códigos de
salida y artefactos:

    1. ghost.run              candidatos T1 (+ ESM-2 + conservación)
    2. ghost.gap_conservation brechas -> pérdida de gen vs artefacto
    3. ghost.gap_targets      enzimas objetivo y posibles cierres
    4. ghost.report           reporte Markdown

Uso:
    python -m ghost.analyze --genome <GCF> --data-root <data> --out ./out
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


def _step(name: str, args: list[str], log) -> float:
    t0 = time.time()
    log(f"\n=== {name} ===")
    subprocess.run([sys.executable, "-m", *args], check=True)
    dt = time.time() - t0
    log(f"--- {name}: {dt:.1f}s")
    return dt


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Pipeline completo de un genoma.")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--esm-top", type=int, default=100)
    ap.add_argument("--no-esm", action="store_true")
    ap.add_argument("--no-gap", action="store_true", help="omite las etapas metabólicas")
    ap.add_argument("--db-proteomes", type=Path, default=None)
    ap.add_argument("--db-swissprot", type=Path, default=None)
    ap.add_argument("--report-top", type=int, default=15)
    a = ap.parse_args(argv)

    out = a.out
    out.mkdir(parents=True, exist_ok=True)
    db_prot = a.db_proteomes or (a.data_root / "db" / "proteomes.dmnd")
    db_sp = a.db_swissprot or (a.data_root / "db" / "swissprot.dmnd")
    log_path = out / f"{a.genome}_analyze.log"

    def log(msg: str) -> None:
        print(msg, flush=True)
        with open(log_path, "a") as fh:
            fh.write(msg + "\n")

    log(f"# analyze {a.genome} · {time.strftime('%Y-%m-%d %H:%M:%S')}")
    steps: list[dict] = []

    # 1. candidatos + ESM + conservación
    run_args = [
        "ghost.run", "--genome", a.genome, "--data-root", str(a.data_root),
        "--out", str(out), "--diamond", "proteomes", "--db", str(db_prot),
        "--threads", str(a.threads),
    ]
    if not a.no_esm:
        run_args += ["--esm", "--esm-top", str(a.esm_top)]
    dt = _step("candidatos", run_args, log)
    steps.append({"step": "candidates", "seconds": round(dt, 1), "ok": True})

    if not a.no_gap:
        # 2. conservación de la reacción por filo
        dt = _step("gap_conservation", [
            "ghost.gap_conservation", "--genome", a.genome,
            "--data-root", str(a.data_root), "--out-dir", str(out),
        ], log)
        steps.append({"step": "gap_conservation", "seconds": round(dt, 1), "ok": True})

        # 3. cierre dirigido sobre las pérdidas reales
        cons = out / f"{a.genome}_gap_conservation.tsv"
        dt = _step("gap_targets", [
            "ghost.gap_targets", "--genome", a.genome,
            "--data-root", str(a.data_root), "--out-dir", str(out),
            "--db", str(db_sp), "--threads", str(a.threads),
            "--conservation-tsv", str(cons),
        ], log)
        steps.append({"step": "gap_targets", "seconds": round(dt, 1), "ok": True})

    # 4. reporte
    dt = _step("report", [
        "ghost.report", "--out-dir", str(out), "--data-root", str(a.data_root),
        "--genomes", a.genome, "--report", str(out / f"{a.genome}_report.md"),
        "--top", str(a.report_top),
    ], log)
    steps.append({"step": "report", "seconds": round(dt, 1), "ok": True})

    manifest = {
        "genome": a.genome,
        "data_root": str(a.data_root),
        "out": str(out),
        "esm": not a.no_esm,
        "gap": not a.no_gap,
        "steps": steps,
        "total_seconds": round(sum(s["seconds"] for s in steps), 1),
        "artifacts": sorted(p.name for p in out.glob(f"{a.genome}*")),
        "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    (out / f"{a.genome}_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=1)
    )
    log(f"\n# manifiesto -> {out / (a.genome + '_analysis_manifest.json')}")
    log(f"# total: {manifest['total_seconds']}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
