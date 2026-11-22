#!/usr/bin/env python3
"""
38_batch_functional_annotation.py — Batch de anotacion funcional (37) sobre todo el corpus.

Corre scripts/metabolism/37_functional_annotation.py sobre cada directorio de
data/annotations_ab_initio/<prefix>/ que tenga proteins.faa, y produce
data/functional_annotation/batch_summary.json.

Caracteristicas (heredadas de 33):
  - Idempotente: si <genome>/functional_annotation.json existe, se omite (--force re-corre)
  - Resumable: el summary se guarda tras cada genoma
  - Ctrl+C: guarda progreso y sale limpio
  - Un genoma que falla no detiene el batch

Pre-requisito (una vez): DB DIAMOND creada (ver docstring de 37).

Uso:
    python scripts/metabolism/38_batch_functional_annotation.py
    python scripts/metabolism/38_batch_functional_annotation.py --only GCF_000006685
    python scripts/metabolism/38_batch_functional_annotation.py --limit 5
    python scripts/metabolism/38_batch_functional_annotation.py --force
"""
from __future__ import annotations

import argparse
import json
import logging
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("batch_func")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT_ROOT / "data" / "annotations_ab_initio"
DEFAULT_SUMMARY_DIR = PROJECT_ROOT / "data" / "functional_annotation"
ANNOTATOR = PROJECT_ROOT / "scripts" / "metabolism" / "37_functional_annotation.py"
SUMMARY_FILE = "batch_summary.json"
SUBPROCESS_TIMEOUT = 7200  # 2h por genoma (DIAMOND es rapido; margen holgado)


def format_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    if seconds < 3600:
        m, s = divmod(int(seconds), 60)
        return f"{m}m{s:02d}s"
    h, rem = divmod(int(seconds), 3600)
    m = rem // 60
    return f"{h}h{m:02d}m"


@dataclass
class Item:
    prefix: str
    genome_dir: str
    status: str  # success | skipped | failed
    elapsed_s: float = 0.0
    cds: int = 0
    hits: int = 0
    with_gene: int = 0
    error: str = ""


def load_func_json(genome_dir: Path) -> dict | None:
    p = genome_dir / "functional_annotation.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def run_one(genome_dir: Path, db: Path, threads: int, force: bool) -> Item:
    prefix = genome_dir.name

    if not force:
        existing = load_func_json(genome_dir)
        if existing is not None:
            rep = existing.get("report", {})
            return Item(
                prefix=prefix, genome_dir=str(genome_dir), status="skipped",
                cds=rep.get("n_cds", 0), hits=rep.get("n_hits", 0),
                with_gene=rep.get("with_gene", 0),
            )

    cmd = [sys.executable, str(ANNOTATOR), str(genome_dir),
           "--db", str(db), "--threads", str(threads)]
    t0 = time.time()
    try:
        res = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=SUBPROCESS_TIMEOUT)
    except subprocess.TimeoutExpired:
        return Item(prefix=prefix, genome_dir=str(genome_dir), status="failed",
                    elapsed_s=time.time() - t0,
                    error=f"timeout >{SUBPROCESS_TIMEOUT}s")

    elapsed = time.time() - t0
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "unknown")[-500:].strip()
        return Item(prefix=prefix, genome_dir=str(genome_dir), status="failed",
                    elapsed_s=elapsed, error=err)

    doc = load_func_json(genome_dir)
    if doc is None:
        return Item(prefix=prefix, genome_dir=str(genome_dir), status="failed",
                    elapsed_s=elapsed, error="functional_annotation.json no generado")

    rep = doc.get("report", {})
    return Item(
        prefix=prefix, genome_dir=str(genome_dir), status="success",
        elapsed_s=elapsed, cds=rep.get("n_cds", 0), hits=rep.get("n_hits", 0),
        with_gene=rep.get("with_gene", 0),
    )


def write_summary(path: Path, items: list[Item], started: str,
                  total_elapsed: float, interrupted: bool) -> None:
    ok = [i for i in items if i.status in ("success", "skipped")]
    doc = {
        "started": started,
        "updated": datetime.now().isoformat(timespec="seconds"),
        "interrupted": interrupted,
        "total_elapsed_s": round(total_elapsed, 1),
        "counts": {
            "total": len(items),
            "success": sum(1 for i in items if i.status == "success"),
            "skipped": sum(1 for i in items if i.status == "skipped"),
            "failed": sum(1 for i in items if i.status == "failed"),
        },
        "totals": {
            "cds": sum(i.cds for i in ok),
            "hits": sum(i.hits for i in ok),
            "with_gene": sum(i.with_gene for i in ok),
        },
        "items": [asdict(i) for i in items],
    }
    path.write_text(json.dumps(doc, indent=2))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Batch de anotacion funcional DIAMOND",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--summary-dir", type=Path, default=DEFAULT_SUMMARY_DIR)
    p.add_argument("--db", type=Path,
                   default=PROJECT_ROOT / "data" / "db" / "swissprot")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--only", default=None,
                   help="Solo genomas cuyo nombre contenga este substring")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.verbose:
        log.setLevel(logging.DEBUG)

    input_dir = args.input_dir.resolve()
    summary_dir = args.summary_dir.resolve()

    if not input_dir.is_dir():
        log.error("No existe: %s", input_dir)
        return 2
    if not ANNOTATOR.exists():
        log.error("Falta: %s", ANNOTATOR)
        return 2

    genome_dirs = sorted(
        p for p in input_dir.iterdir()
        if p.is_dir() and (p / "proteins.faa").exists()
    )
    if args.only:
        genome_dirs = [g for g in genome_dirs if args.only in g.name]
    if args.limit:
        genome_dirs = genome_dirs[: args.limit]

    if not genome_dirs:
        log.error("No hay genomas con proteins.faa en %s", input_dir)
        return 2

    summary_dir.mkdir(parents=True, exist_ok=True)

    pre_existing = sum(
        1 for g in genome_dirs if load_func_json(g) is not None
    )

    log.info("=== Batch anotacion funcional (DIAMOND vs Swiss-Prot) ===")
    log.info("Input  : %s", input_dir)
    log.info("DB     : %s", args.db)
    log.info("Threads: %d", args.threads)
    log.info("Genomas: %d (ya anotados: %d)", len(genome_dirs), pre_existing)
    log.info("Force  : %s", args.force)

    items: list[Item] = []
    summary_path = summary_dir / SUMMARY_FILE
    started = datetime.now().isoformat(timespec="seconds")
    t0 = time.time()
    interrupted = False

    def handle_sigint(sig, frame):
        nonlocal interrupted
        interrupted = True
        log.warning("Ctrl+C detectado: guardando progreso y saliendo...")

    signal.signal(signal.SIGINT, handle_sigint)

    try:
        for idx, gdir in enumerate(genome_dirs, 1):
            if interrupted:
                break
            log.info("[%d/%d] %s", idx, len(genome_dirs), gdir.name)

            item = run_one(gdir, args.db, args.threads, args.force)
            items.append(item)

            if item.status == "success":
                cov = 100 * item.hits / item.cds if item.cds else 0
                log.info("  OK %.1fs | CDS=%d hits=%d (%.1f%%) gene=%d",
                         item.elapsed_s, item.cds, item.hits, cov, item.with_gene)
            elif item.status == "skipped":
                log.info("  SKIP (ya existia)")
            else:
                log.warning("  FAIL: %s", item.error[:200])

            if idx > 1:
                avg = (time.time() - t0) / idx
                log.info("  ETA %s", format_duration(avg * (len(genome_dirs) - idx)))

            write_summary(summary_path, items, started,
                          time.time() - t0, interrupted)
    finally:
        write_summary(summary_path, items, started,
                      time.time() - t0, interrupted)

    n_ok = sum(1 for i in items if i.status == "success")
    n_skip = sum(1 for i in items if i.status == "skipped")
    n_fail = sum(1 for i in items if i.status == "failed")

    log.info("=== Fin ===")
    log.info("Success: %d | Skipped: %d | Failed: %d", n_ok, n_skip, n_fail)
    log.info("Tiempo total: %s", format_duration(time.time() - t0))
    log.info("Summary: %s", summary_path)

    if n_fail:
        log.warning("Genomas fallidos:")
        for i in items:
            if i.status == "failed":
                log.warning("  %s — %s", i.prefix, i.error[:150])
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
