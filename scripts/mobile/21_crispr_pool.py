#!/usr/bin/env python3
"""21_crispr_pool.py — barrido de arrays CRISPR de los 400 genomas.

Reutiliza `detect_genome`/`write_detection` de `04_crispr_arrays.py`. Cada
genoma es independiente y reanudable (`crispr_stats.json` existente = saltar).
Los intervalos CDS para el filtro de repeticiones intragen se leen directamente
del GFF ab initio, sin depender de `01_extract_candidates.py`.

Uso:
    python scripts/mobile/21_crispr_pool.py --workers 4
    python scripts/mobile/21_crispr_pool.py --workers 4 --limit 5
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import importlib.util
import json
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data" / "mobile_typing"
WORK = BASE / "_work" / "full"
STATE = WORK / "state"
ANN = ROOT / "data" / "annotations_ab_initio"
RAW = ROOT / "data" / "raw"

_spec = importlib.util.spec_from_file_location(
    "crispr04", Path(__file__).resolve().parent / "04_crispr_arrays.py"
)
CRISPR = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CRISPR)


def cds_intervals(gff: Path):
    by_contig = defaultdict(list)
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9 and parts[2] == "CDS":
                by_contig[parts[0]].append((int(parts[3]), int(parts[4])))
    for contig in by_contig:
        by_contig[contig].sort()
    return by_contig


def one_genome(genome: str) -> dict:
    outdir = WORK / genome
    if (outdir / "crispr_stats.json").exists():
        return {"genome": genome, "status": "skip"}
    if not (RAW / f"{genome}_genomic.fna.gz").exists():
        return {"genome": genome, "status": "skip_no_fasta"}
    outdir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    try:
        fna = RAW / f"{genome}_genomic.fna.gz"
        gff = ANN / genome / "annotation.gff3"
        arrays, spacers, stats = CRISPR.detect_genome(genome, fna, cds_intervals(gff))
        CRISPR.write_detection(outdir, arrays, spacers, stats)
    except Exception as exc:  # noqa: BLE001
        return {"genome": genome, "status": "fail", "error": f"{type(exc).__name__}: {exc}"[:300]}
    return {"genome": genome, "status": "ok", "arrays": stats["arrays"],
            "spacers": stats["spacers"], "seconds": round(time.time() - started, 1)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    STATE.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((BASE / "full_genomes.tsv").open(), delimiter="\t"))
    todo = [r["genome"] for r in rows if not r["missing_inputs"]]
    if args.limit:
        todo = todo[: args.limit]
    print(f"CRISPR pool: {len(todo)} genomas, {args.workers} workers", flush=True)
    done = {"ok": 0, "skip": 0, "skip_no_fasta": 0, "fail": 0}
    failures = {}
    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(one_genome, todo):
            done[result["status"]] = done.get(result["status"], 0) + 1
            if result["status"] == "fail":
                failures[result["genome"]] = result.get("error", "")
            n = sum(done.values())
            if n % 20 == 0 or result["status"] == "fail":
                print(f"  [{n}/{len(todo)}] ok={done['ok']} skip={done['skip']} "
                      f"sin_fasta={done['skip_no_fasta']} fail={done['fail']} "
                      f"elapsed={(time.time()-start)/60:.1f} min", flush=True)
                (STATE / "crispr_progress.json").write_text(json.dumps(
                    {"done": n, **done, "failures": failures}, indent=2, sort_keys=True) + "\n")
    (STATE / "crispr_failures.json").write_text(json.dumps(failures, indent=2, sort_keys=True) + "\n")
    print(f"FIN: ok={done['ok']} skip={done['skip']} sin_fasta={done['skip_no_fasta']} "
          f"fail={done['fail']} "
          f"en {(time.time()-start)/60:.1f} min", flush=True)
    return 0 if done["fail"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
