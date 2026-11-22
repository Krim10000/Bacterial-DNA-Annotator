#!/usr/bin/env python3
"""Corre la capa de membrana en los 400 genomas, por lotes y con reanudacion.

Salidas:
  data/membrane/full/<genome>.tsv
  data/membrane/gff/full/<genome>.gff3
  data/membrane/full/genomes.tsv          (lista + metadatos)
  data/membrane/full/state.json           (estado reanudable por genoma)
  data/membrane/full/batches/batch_XXX.json
  data/membrane/full/run.log

No toca el piloto de 10 (data/membrane/*.tsv, data/membrane/gff/*.gff3).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCRIPT03 = HERE / "03_annotate_membrane.py"
HEADER = ["genome", "contig", "start", "end", "strand", "protein_id", "type",
          "family", "substrate", "evidence", "confidence"]
HEADER_LINE = "\t".join(HEADER)


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, path)


def load_meta() -> tuple[dict, dict, dict]:
    species = {}
    p = ROOT / "data/species_phylum.tsv"
    if p.exists():
        for r in csv.DictReader(open(p), delimiter="\t"):
            species[r["assembly_accession"]] = r
    asm = {}
    p = ROOT / "data/assembly_summary.txt"
    if p.exists():
        with open(p, errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) > 8:
                    asm[f[0]] = f[7]
    return species, asm, {}


def accession(genome: str) -> str:
    m = re.match(r"^(GCF_\d+\.\d+)", genome)
    return m.group(1) if m else genome


def fast_stats(annot_dir: Path) -> dict:
    contigs = {}
    cds = 0
    gff = annot_dir / "annotation.gff3"
    for line in open(gff, errors="replace"):
        if line.startswith("##sequence-region"):
            f = line.split()
            if len(f) >= 4:
                contigs[f[1]] = int(f[3])
        elif "\tCDS\t" in line:
            cds += 1
    nprot = sum(1 for line in open(annot_dir / "proteins.faa") if line.startswith(">"))
    return dict(contigs=len(contigs), bp=sum(contigs.values()), cds=cds, proteins=nprot)


def build_list(annot_root: Path, outdir: Path) -> list[dict]:
    species, asm, _ = load_meta()
    rows = []
    for d in sorted(annot_root.iterdir()):
        if not d.is_dir():
            continue
        g = d.name
        acc = accession(g)
        sp = species.get(acc, {})
        st = fast_stats(d)
        rows.append(dict(
            genome=g,
            phylum=sp.get("phylum", "unknown") or "unknown",
            organism=sp.get("organism_name") or asm.get(acc, g),
            contigs=st["contigs"], bp=st["bp"], cds=st["cds"], proteins=st["proteins"],
        ))
    with open(outdir / "genomes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["genome", "phylum", "organism", "contigs", "bp",
                                       "cds", "proteins"])
        w.writeheader()
        w.writerows(rows)
    return rows


def read_list(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path), delimiter="\t"))


def count_rows(path: Path) -> int:
    try:
        with open(path) as fh:
            next(fh)
            return sum(1 for _ in fh)
    except Exception:
        return 0


def is_complete(outdir: Path, gff_dir: Path, genome: str) -> bool:
    tsv = outdir / f"{genome}.tsv"
    gff = gff_dir / f"{genome}.gff3"
    if not (tsv.exists() and gff.exists()):
        return False
    if tsv.stat().st_size == 0 or gff.stat().st_size == 0:
        return False
    try:
        with open(tsv) as fh:
            return fh.readline().rstrip("\n") == HEADER_LINE
    except OSError:
        return False


def run_one(meta: dict, args, tmp_root: Path, outdir: Path, gff_dir: Path) -> dict:
    g = meta["genome"]
    t0 = time.time()
    tmpdir = Path(tempfile.mkdtemp(prefix=f"{g}_", dir=str(tmp_root)))
    try:
        cmd = [sys.executable, str(SCRIPT03), "--genome", g,
               "--annot-root", str(args.annot_root), "--panel-dir", str(args.panel_dir),
               "--outdir", str(tmpdir), "--cpu", str(args.cpu_per_worker)]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
        except subprocess.TimeoutExpired:
            return dict(status="failed", elapsed=round(time.time() - t0, 1),
                        rows=0, error=f"timeout {args.timeout}s")
        if res.returncode != 0:
            err = (res.stderr or res.stdout or "").strip()[-2000:]
            return dict(status="failed", elapsed=round(time.time() - t0, 1),
                        rows=0, error=err)
        tmp_tsv = tmpdir / f"{g}.tsv"
        tmp_gff = tmpdir / "gff" / f"{g}.gff3"
        if not tmp_tsv.exists() or not tmp_gff.exists():
            return dict(status="failed", elapsed=round(time.time() - t0, 1),
                        rows=0, error="salida incompleta (TSV/GFF)")
        outdir.mkdir(parents=True, exist_ok=True)
        gff_dir.mkdir(parents=True, exist_ok=True)
        os.replace(tmp_tsv, outdir / f"{g}.tsv")
        os.replace(tmp_gff, gff_dir / f"{g}.gff3")
        rows = count_rows(outdir / f"{g}.tsv")
        return dict(status="done", elapsed=round(time.time() - t0, 1), rows=rows, error="")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annot-root", default=str(ROOT / "data/annotations_ab_initio"))
    ap.add_argument("--outdir", default=str(ROOT / "data/membrane/full"))
    ap.add_argument("--gff-dir", default=str(ROOT / "data/membrane/gff/full"))
    ap.add_argument("--panel-dir", default=str(ROOT / "data/membrane/pfam_panel"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cpu-per-worker", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=50)
    ap.add_argument("--limit", type=int, default=0, help="solo los primeros N (prueba)")
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--retries", type=int, default=1)
    ap.add_argument("--reconcile-only", action="store_true",
                    help="solo reconciliar estado con salidas existentes; no anotar")
    ap.add_argument("--nice", type=int, default=10)
    ap.add_argument("--list", default=str(ROOT / "data/membrane/full/genomes.tsv"))
    ap.add_argument("--state", default=str(ROOT / "data/membrane/full/state.json"))
    ap.add_argument("--log", default=str(ROOT / "data/membrane/full/run.log"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    gff_dir = Path(args.gff_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    gff_dir.mkdir(parents=True, exist_ok=True)
    (outdir / "batches").mkdir(parents=True, exist_ok=True)
    tmp_root = outdir / ".tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)
    if args.nice:
        try:
            os.nice(args.nice)
        except OSError:
            pass

    log_path = Path(args.log)

    def log(msg: str) -> None:
        line = f"[{now()}] {msg}"
        print(line, flush=True)
        with open(log_path, "a") as fh:
            fh.write(line + "\n")

    if not Path(args.list).exists():
        log(f"generando lista {args.list}")
        metas = build_list(Path(args.annot_root), outdir)
    else:
        metas = read_list(Path(args.list))
        log(f"lista existente: {len(metas)} genomas")
    if args.limit:
        metas = metas[:args.limit]
    log(f"genomas a procesar: {len(metas)}; workers={args.workers}; cpu/worker={args.cpu_per_worker}")

    state_path = Path(args.state)
    state = {}
    if state_path.exists():
        try:
            state = json.load(open(state_path))
        except Exception:
            state = {}
    state.setdefault("schema_version", 1)
    state.setdefault("genomes", {})
    state.setdefault("batches", {})
    state["outdir"] = str(outdir)
    state["gff_dir"] = str(gff_dir)
    state["updated_at"] = now()

    # reconciliar estado con archivos ya presentes
    for m in metas:
        g = m["genome"]
        entry = state["genomes"].get(g, {})
        done = is_complete(outdir, gff_dir, g)
        if done:
            state["genomes"][g] = dict(status="done", elapsed_s=entry.get("elapsed_s", 0),
                                        rows=count_rows(outdir / f"{g}.tsv"),
                                        error="", batch=entry.get("batch", "pre"))
        elif entry.get("status") == "done":
            state["genomes"][g] = dict(status="pending", elapsed_s=0, rows=0,
                                        error="salida ausente", batch=entry.get("batch", "pre"))
    atomic_json(state_path, state)

    total = len(metas)
    done_count = sum(1 for m in metas if state["genomes"].get(m["genome"], {}).get("status") == "done")
    log(f"estado inicial: {done_count}/{total} ya completos")

    if args.reconcile_only:
        done = [m for m in metas if state["genomes"].get(m["genome"], {}).get("status") == "done"]
        failed = [m for m in metas if state["genomes"].get(m["genome"], {}).get("status") == "failed"]
        pending = [m for m in metas if m not in done and m not in failed]
        state["interrupted"] = True
        state["last_run_at"] = now()
        state["totals"] = dict(total=total, done=len(done), failed=len(failed),
                               pending=len(pending),
                               rows=sum(state["genomes"][m["genome"]]["rows"] for m in done),
                               updated_at=now())
        atomic_json(state_path, state)
        log(f"RECONCILIADO: {len(done)}/{total} OK; {len(failed)} fallos; {len(pending)} pendientes; filas={state['totals']['rows']}")
        return 0 if not failed else 2

    batch_size = max(1, args.batch_size)
    for b0 in range(0, total, batch_size):
        batch = metas[b0:b0 + batch_size]
        bnum = b0 // batch_size + 1
        pending = [m for m in batch if not is_complete(outdir, gff_dir, m["genome"])]
        if not pending:
            log(f"lote {bnum:03d}: ya completo ({len(batch)} genomas)")
            continue
        log(f"lote {bnum:03d}: {len(pending)} pendientes de {len(batch)}")
        bstart = now()
        t0 = time.time()
        results = []
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(run_one, m, args, tmp_root, outdir, gff_dir): m for m in pending}
            for fut in as_completed(futs):
                m = futs[fut]
                g = m["genome"]
                try:
                    r = fut.result()
                except Exception as exc:
                    r = dict(status="failed", elapsed=0, rows=0, error=f"{type(exc).__name__}: {exc}")
                r["batch"] = bnum
                state["genomes"][g] = dict(status=r["status"], elapsed_s=r["elapsed"],
                                            rows=r["rows"], error=r["error"], batch=bnum)
                results.append((g, r))
                log(f"  {g}: {r['status']} rows={r['rows']} {r['elapsed']}s"
                    + (f" error={r['error'][:200]}" if r["error"] else ""))
                atomic_json(state_path, state)
        bdone = sum(1 for _, r in results if r["status"] == "done")
        bfail = len(results) - bdone
        state["batches"][f"batch_{bnum:03d}"] = dict(
            started_at=bstart, finished_at=now(), total=len(batch), pending=len(pending),
            completed=bdone, failed=bfail, elapsed_s=round(time.time() - t0, 1),
            errors=[{"genome": g, "error": r["error"]} for g, r in results if r["status"] != "done"],
        )
        atomic_json(state_path, state)
        atomic_json(outdir / "batches" / f"batch_{bnum:03d}.json", state["batches"][f"batch_{bnum:03d}"])
        log(f"lote {bnum:03d}: {bdone} OK, {bfail} fallos, {state['batches'][f'batch_{bnum:03d}']['elapsed_s']}s")

    # un reintento para fallos
    for attempt in range(1, args.retries + 1):
        failed = [m for m in metas
                  if state["genomes"].get(m["genome"], {}).get("status") == "failed"]
        if not failed:
            break
        log(f"reintento {attempt}: {len(failed)} genomas")
        for m in failed:
            g = m["genome"]
            r = run_one(m, args, tmp_root, outdir, gff_dir)
            state["genomes"][g] = dict(status=r["status"], elapsed_s=r["elapsed"],
                                        rows=r["rows"], error=r["error"],
                                        batch=state["genomes"].get(g, {}).get("batch", "retry"))
            log(f"  retry {g}: {r['status']} rows={r['rows']} {r['elapsed']}s"
                + (f" error={r['error'][:200]}" if r["error"] else ""))
            atomic_json(state_path, state)

    done = [m for m in metas if state["genomes"].get(m["genome"], {}).get("status") == "done"]
    failed = [m for m in metas if state["genomes"].get(m["genome"], {}).get("status") != "done"]
    pending_n = total - len(done) - len(failed)
    state["totals"] = dict(total=total, done=len(done), failed=len(failed),
                           pending=pending_n,
                           rows=sum(state["genomes"][m["genome"]]["rows"] for m in done),
                           updated_at=now())
    atomic_json(state_path, state)
    log(f"FIN: {len(done)}/{total} OK; {len(failed)} fallos; filas={state['totals']['rows']}")
    if failed:
        log("fallos: " + ", ".join(m["genome"] for m in failed))
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
