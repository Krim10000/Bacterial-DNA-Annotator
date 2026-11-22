#!/usr/bin/env python3
"""20_pfam_pool.py — barrido Pfam de los 400 genomas, reanudable y paralelo.

Reutiliza el panel curado y el parser de `03_scan_pfam.py`. Cada genoma es
independiente: si `_work/full/<genoma>/pfam_hits.tsv` ya existe, se salta.
`--cpu 1` por proceso y N procesos en paralelo (paralelismo moderado).

Uso:
    python scripts/mobile/20_pfam_pool.py --workers 6
    python scripts/mobile/20_pfam_pool.py --workers 6 --limit 5   # prueba
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import importlib.util
import json
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data" / "mobile_typing"
WORK = BASE / "_work" / "full"
STATE = WORK / "state"
DEFAULT_PANEL = BASE / "_work" / "pfam_core" / "panel.hmm"
DEFAULT_MANIFEST = BASE / "_work" / "pfam_core" / "panel_manifest.tsv"
NUCLEUS = Path(__file__).resolve().parent / "03_scan_pfam.py"


def load_module():
    spec = importlib.util.spec_from_file_location("scan03", NUCLEUS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCAN = load_module()
PANEL = None
META = {}
LOCK = threading.Lock()
START = time.time()
DONE = {"ok": 0, "skip": 0, "fail": 0}
FAILURES = {}


def protein_count(faa: Path) -> int:
    with faa.open() as fh:
        return sum(1 for line in fh if line.startswith(">"))


def one_genome(genome: str, faa: Path, outdir: Path) -> dict:
    hits_tsv = outdir / "pfam_hits.tsv"
    stats_json = outdir / "pfam_scan_stats.json"
    if hits_tsv.exists() and stats_json.exists() and hits_tsv.stat().st_size > 0:
        return {"genome": genome, "status": "skip"}
    outdir.mkdir(parents=True, exist_ok=True)
    domtbl = outdir / "pfam.domtbl"
    cmd = [
        "hmmscan", "--cpu", "1", "--noali", "--domE", "1e-3",
        "--domtblout", str(domtbl), "-o", "/dev/null", str(PANEL), str(faa),
    ]
    started = time.time()
    try:
        subprocess.run(cmd, check=True, timeout=2400, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE)
    except Exception as exc:  # reintento único
        try:
            subprocess.run(cmd, check=True, timeout=2400, stdout=subprocess.DEVNULL,
                           stderr=subprocess.PIPE)
        except Exception as exc2:
            return {"genome": genome, "status": "fail", "error": f"{exc2}"[:300]}
    hits = SCAN.parse_domtbl(domtbl)
    rows = []
    n_ga = 0
    for hit in hits:
        info = META[hit["target_acc"]]
        passes = hit["dom_score"] >= info["ga_dom"] > 0
        n_ga += int(passes)
        rows.append(
            {
                "protein_id": hit["query"],
                "pfam_acc": hit["target_acc"],
                "pfam_name": info["name"],
                "pfam_desc": info["desc"],
                "categories": info["categories"],
                "passes_ga": "yes" if passes else "no",
                "ievalue": f"{hit['ievalue']:.3g}",
                "score": f"{hit['score']:.1f}",
                "hmm_from": hit["hmm_from"],
                "hmm_to": hit["hmm_to"],
                "ali_from": hit["ali_from"],
                "ali_to": hit["ali_to"],
                "env_from": hit["env_from"],
                "env_to": hit["env_to"],
            }
        )
    rows.sort(key=lambda r: (r["protein_id"], r["pfam_acc"]))
    fields = ["protein_id", "pfam_acc", "pfam_name", "pfam_desc", "categories",
              "passes_ga", "ievalue", "score", "hmm_from", "hmm_to", "ali_from",
              "ali_to", "env_from", "env_to"]
    tmp = hits_tsv.with_suffix(".tmp")
    with tmp.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(hits_tsv)
    stats = {
        "proteins": protein_count(faa),
        "hits": len(hits),
        "hits_passing_ga": n_ga,
        "models": len(META),
        "seconds": round(time.time() - started, 1),
    }
    stats_json.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    domtbl.unlink(missing_ok=True)
    return {"genome": genome, "status": "ok", "hits": len(hits),
            "seconds": stats["seconds"], "proteins": stats["proteins"]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--family", default="", help="filtra por prefijo de genoma")
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    global PANEL, META
    PANEL = args.panel.resolve()
    META = SCAN.load_manifest(args.manifest.resolve())

    STATE.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((BASE / "full_genomes.tsv").open(), delimiter="\t"))
    todo = []
    for row in rows:
        genome = row["genome"]
        if "proteins" in (row["missing_inputs"] or "").split(","):
            continue
        if args.family and not genome.startswith(args.family):
            continue
        todo.append((genome, ROOT / "data" / "annotations_ab_initio" / genome / "proteins.faa"))
    counts = {g: protein_count(faa) for g, faa in todo}
    todo.sort(key=lambda x: -counts[x[0]])  # los más grandes primero (mejor reparto)
    if args.limit:
        todo = todo[: args.limit]
    print(f"Pfam pool: {len(todo)} genomas, {args.workers} workers, "
          f"{sum(counts[g] for g, _ in todo)} proteínas", flush=True)

    def record(result):
        with LOCK:
            DONE[result["status"]] += 1
            if result["status"] == "fail":
                FAILURES[result["genome"]] = result.get("error", "")
            n = sum(DONE.values())
            if n % 10 == 0 or result["status"] == "fail":
                elapsed = time.time() - START
                print(f"  [{n}/{len(todo)}] ok={DONE['ok']} skip={DONE['skip']} "
                      f"fail={DONE['fail']} elapsed={elapsed/60:.1f} min", flush=True)
                (STATE / "pfam_progress.json").write_text(json.dumps(
                    {"done": n, "ok": DONE["ok"], "skip": DONE["skip"],
                     "fail": DONE["fail"], "failures": FAILURES}, indent=2, sort_keys=True) + "\n")

    if args.workers <= 1:
        for genome, faa in todo:
            record(one_genome(genome, faa, WORK / genome))
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(one_genome, genome, faa, WORK / genome): genome
                for genome, faa in todo
            }
            for future in concurrent.futures.as_completed(futures):
                record(future.result())
    (STATE / "pfam_failures.json").write_text(json.dumps(FAILURES, indent=2, sort_keys=True) + "\n")
    print(f"FIN: ok={DONE['ok']} skip={DONE['skip']} fail={DONE['fail']} "
          f"en {(time.time()-START)/60:.1f} min", flush=True)
    return 0 if DONE["fail"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
