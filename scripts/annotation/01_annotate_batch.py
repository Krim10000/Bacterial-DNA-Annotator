#!/usr/bin/env python3
"""
01_annotate_batch.py — Batch ab initio annotation over data/raw/.

Corre scripts/annotation/00_annotate_single.py sobre cada FASTA de un directorio,
escribe en data/annotations_ab_initio/<prefix>/, y produce batch_summary.json.

Caracteristicas:
  - Idempotente: si <out>/<prefix>/report.json existe, se omite (usar --force)
  - Resumable: el summary se guarda tras cada genoma
  - Ctrl+C: guarda progreso y sale limpiamente
  - Aisla fallos: un genoma que crashea no detiene el batch

Uso:
    python scripts/annotation/01_annotate_batch.py
    python scripts/annotation/01_annotate_batch.py --limit 5
    python scripts/annotation/01_annotate_batch.py --only GCF_000006685
    python scripts/annotation/01_annotate_batch.py --force
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
log = logging.getLogger("batch")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT_ROOT / "data" / "raw"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "annotations_ab_initio"
ANNOTATOR = PROJECT_ROOT / "scripts" / "annotation" / "00_annotate_single.py"
SUMMARY_FILE = "batch_summary.json"
SUBPROCESS_TIMEOUT = 3600


def prefix_from_filename(name: str) -> str:
    for suffix in ("_genomic.fna.gz", ".fna.gz", ".fasta.gz", ".fa.gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


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
    input: str
    status: str  # success | skipped | failed
    elapsed_s: float = 0.0
    cds: int = 0
    rrna: int = 0
    trna: int = 0
    error: str = ""


def load_report(out_dir: Path) -> dict | None:
    rp = out_dir / "report.json"
    if not rp.exists():
        return None
    try:
        return json.loads(rp.read_text())
    except Exception:
        return None


def run_one(input_path: Path, out_root: Path, force: bool) -> Item:
    prefix = prefix_from_filename(input_path.name)
    out_dir = out_root / prefix

    if not force:
        existing = load_report(out_dir)
        if existing is not None:
            return Item(
                prefix=prefix,
                input=str(input_path),
                status="skipped",
                elapsed_s=existing.get("elapsed_seconds", 0.0),
                cds=existing.get("cds_features", 0),
                rrna=existing.get("rrna_features", 0),
                trna=existing.get("trna_features", 0),
            )

    cmd = [sys.executable, str(ANNOTATOR), str(input_path), "--prefix", prefix]
    t0 = time.time()
    try:
        res = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=SUBPROCESS_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return Item(
            prefix=prefix, input=str(input_path), status="failed",
            elapsed_s=time.time() - t0,
            error=f"timeout >{SUBPROCESS_TIMEOUT}s",
        )

    elapsed = time.time() - t0
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "unknown")[-500:].strip()
        return Item(
            prefix=prefix, input=str(input_path), status="failed",
            elapsed_s=elapsed, error=err,
        )

    rep = load_report(out_dir)
    if rep is None:
        return Item(
            prefix=prefix, input=str(input_path), status="failed",
            elapsed_s=elapsed, error="report.json no generado",
        )

    return Item(
        prefix=prefix,
        input=str(input_path),
        status="success",
        elapsed_s=elapsed,
        cds=rep.get("cds_features", 0),
        rrna=rep.get("rrna_features", 0),
        trna=rep.get("trna_features", 0),
    )


def write_summary(
    path: Path, items: list[Item], started: str,
    total_elapsed: float, interrupted: bool,
) -> None:
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
            "rrna": sum(i.rrna for i in ok),
            "trna": sum(i.trna for i in ok),
            "bp": sum(0 for _ in ok),
        },
        "items": [asdict(i) for i in items],
    }
    path.write_text(json.dumps(doc, indent=2))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Batch ab initio annotation over data/raw/",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--pattern", default="*.fna.gz")
    p.add_argument("--limit", type=int, default=None,
                   help="Solo los primeros N genomas (para test)")
    p.add_argument("--only", default=None,
                   help="Solo genomas cuyo nombre contenga este substring")
    p.add_argument("--force", action="store_true",
                   help="Re-corre aunque exista report.json")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.verbose:
        log.setLevel(logging.DEBUG)

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not input_dir.is_dir():
        log.error("No existe: %s", input_dir)
        return 2
    if not ANNOTATOR.exists():
        log.error("Falta: %s", ANNOTATOR)
        return 2

    files = sorted(input_dir.glob(args.pattern))
    if args.only:
        files = [f for f in files if args.only in f.name]
    if args.limit:
        files = files[: args.limit]

    if not files:
        log.error("No hay archivos en %s con patron %s", input_dir, args.pattern)
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)

    pre_existing = sum(
        1 for f in files
        if load_report(output_dir / prefix_from_filename(f.name)) is not None
    )

    log.info("=== Batch ab_initio ===")
    log.info("Input  : %s", input_dir)
    log.info("Output : %s", output_dir)
    log.info("Files  : %d (ya anotados: %d)", len(files), pre_existing)
    log.info("Force  : %s", args.force)

    items: list[Item] = []
    summary_path = output_dir / SUMMARY_FILE
    started = datetime.now().isoformat(timespec="seconds")
    t0 = time.time()
    interrupted = False

    def handle_sigint(sig, frame):
        nonlocal interrupted
        interrupted = True
        log.warning("Ctrl+C detectado: guardando progreso y saliendo...")

    signal.signal(signal.SIGINT, handle_sigint)

    try:
        for idx, input_path in enumerate(files, 1):
            if interrupted:
                break

            prefix = prefix_from_filename(input_path.name)
            log.info("[%d/%d] %s", idx, len(files), prefix)

            item = run_one(input_path, output_dir, args.force)
            items.append(item)

            if item.status == "success":
                log.info("  OK %.1fs | CDS=%d rRNA=%d tRNA=%d",
                         item.elapsed_s, item.cds, item.rrna, item.trna)
            elif item.status == "skipped":
                log.info("  SKIP (ya existia, %.1fs)", item.elapsed_s)
            else:
                log.warning("  FAIL: %s", item.error[:200])

            if idx > 1:
                avg = (time.time() - t0) / idx
                eta = avg * (len(files) - idx)
                log.info("  ETA %s", format_duration(eta))

            write_summary(summary_path, items, started,
                          time.time() - t0, interrupted)
    finally:
        write_summary(summary_path, items, started,
                      time.time() - t0, interrupted)

    n_success = sum(1 for i in items if i.status == "success")
    n_skip = sum(1 for i in items if i.status == "skipped")
    n_fail = sum(1 for i in items if i.status == "failed")

    log.info("=== Fin ===")
    log.info("Success: %d | Skipped: %d | Failed: %d",
             n_success, n_skip, n_fail)
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
