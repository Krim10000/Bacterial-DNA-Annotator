#!/usr/bin/env python3
"""05_run_homology.py — DIAMOND all-hits de cada genoma contra los productores.

Para cada genoma del lote se buscan **todos** los hits (hasta 500 por proteína,
evalue <= 1e-3) contra la DB de enzimas productoras no eucariotas, en las dos
proteomas disponibles:

  * `ncbi`     : `data/annotations/<genome>.faa.gz` (proteínas anotadas)
  * `abinitio` : `data/annotations_ab_initio/<genome>/proteins.faa` (ORFs)

Los hits crudos quedan en `work/diamond/<genome>.tsv.gz` con la fuente como
primera columna; el mejor hit por accesión y brecha se resume en
`work/homology_hits.tsv`. Los umbrales solo etiquetan: la evidencia numérica
(pident/qcov/evalue/bitscore) se conserva completa.

Uso:
    python scripts/evidence/05_run_homology.py [--workers 8]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ABINITIO, ANNOT, MAX_EVALUE, MIN_PIDENT_STRONG,
                    MIN_PIDENT_WEAK, MIN_QCOV_STRONG, MIN_QCOV_WEAK, WORK,
                    read_tsv, write_json, write_tsv)

DB = WORK / "producer_db.dmnd"
RAW = WORK / "diamond"
RAW.mkdir(parents=True, exist_ok=True)
OUTFMT = ["6", "qseqid", "sseqid", "pident", "length", "qlen",
         "slen", "qcovhsp", "evalue", "bitscore"]
SOURCE_RANK = {"ncbi": 0, "abinitio": 1}


def plan() -> dict[str, list[dict]]:
    rows = read_tsv(WORK / "selected_gaps.tsv")
    acc_info = {r["acc"]: r for r in read_tsv(WORK / "acc_info.tsv")}
    by_genome: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        prok = [a for a in r["producer_accs"].split(";")
                if a and acc_info[a]["euk"] == "0"]
        by_genome[r["genome"]].append({"gap": r["gap"], "accs": prok})
    return by_genome


def _decompress(src: Path, dst: Path) -> None:
    with gzip.open(src, "rt") as fh, open(dst, "w") as out:
        out.write(fh.read())


def _scan_raw(genome: str) -> dict[str, bool]:
    """Ejecuta DIAMOND por fuente y deja el crudo etiquetado; True si corrió."""
    raw_path = RAW / f"{genome}.tsv.gz"
    if raw_path.exists() and raw_path.stat().st_size > 0:
        return {}
    sources = {"ncbi": ANNOT / f"{genome}.faa.gz",
               "abinitio": ABINITIO / genome / "proteins.faa"}
    ok: dict[str, bool] = {}
    with tempfile.TemporaryDirectory(dir=WORK) as td, \
            gzip.open(raw_path, "wt") as gz:
        for name, src in sources.items():
            if not src.exists():
                ok[name] = False
                continue
            q = src
            if str(q).endswith(".gz"):
                q = Path(td) / f"{name}.faa"
                _decompress(src, q)
            out = Path(td) / f"{name}.tsv"
            subprocess.run(
                ["diamond", "blastp", "--query", str(q), "--db", str(DB),
                 "--out", str(out), "--outfmt", *OUTFMT, "--evalue", "1e-3",
                 "--max-target-seqs", "500", "--threads", "1", "--quiet"],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with open(out) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 9:
                        gz.write(name + "\t" + "\t".join(p) + "\n")
            ok[name] = True
    return ok


def _scan_one(job: tuple[str, list[dict]]) -> tuple[str, list[dict]]:
    genome, specs = job
    _scan_raw(genome)
    # accesiones que interesan a este genoma
    wanted: set[str] = set()
    gap_accs: dict[str, list[str]] = {}
    for r in specs:
        gap_accs[r["gap"]] = r["accs"]
        wanted.update(r["accs"])

    # streaming del crudo: mejor hit y contadores por accesión
    per_acc: dict[str, dict] = defaultdict(
        lambda: {"hits": 0, "weak": 0, "strong": 0})
    with gzip.open(RAW / f"{genome}.tsv.gz", "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            source, qseqid, acc = p[0], p[1], p[2]
            if acc not in wanted:
                continue
            pident, qcov, ev, bits = (float(p[3]), float(p[7]) / 100.0,
                                      float(p[8]), float(p[9]))
            strong = (pident >= MIN_PIDENT_STRONG and qcov >= MIN_QCOV_STRONG
                      and ev <= MAX_EVALUE)
            weak = (pident >= MIN_PIDENT_WEAK and qcov >= MIN_QCOV_WEAK
                    and ev <= MAX_EVALUE)
            rec = per_acc[acc]
            rec["hits"] += 1
            rec["strong"] += int(strong)
            rec["weak"] += int(weak)
            score = (int(strong), int(weak), bits, pident, -ev,
                     -SOURCE_RANK.get(source, 9))
            if "best" not in rec or score > rec["best_score"]:
                rec["best_score"] = score
                rec["best"] = {"acc": acc, "source": source, "qseqid": qseqid,
                               "pident": pident, "qcov": qcov, "evalue": ev,
                               "bitscore": bits}

    rows = []
    for r in specs:
        accs = r["accs"]
        if not accs:
            rows.append(_empty(genome, r["gap"], "not_checked",
                               "sin productor no eucariota"))
            continue
        hit_accs = [a for a in accs if a in per_acc and "best" in per_acc[a]]
        if not hit_accs:
            rows.append(_empty(genome, r["gap"], "no_homolog", "sin hit"))
            continue
        best = max((per_acc[a]["best"] for a in hit_accs),
                   key=lambda b: (per_acc[b["acc"]]["best_score"]))
        rec = per_acc[best["acc"]]
        strong = any(per_acc[a]["strong"] for a in hit_accs)
        weak = any(per_acc[a]["weak"] for a in hit_accs)
        status = ("homolog_present" if strong else
                  "homolog_weak" if weak else "homolog_low")
        rows.append({
            "genome": genome, "gap": r["gap"], "best_acc": best["acc"],
            "source": best["source"], "qseqid": best["qseqid"],
            "pident": f"{best['pident']:.1f}", "qcov": f"{best['qcov']:.3f}",
            "evalue": f"{best['evalue']:.3g}",
            "bitscore": f"{best['bitscore']:.1f}",
            "n_hits": sum(per_acc[a]["hits"] for a in hit_accs),
            "n_accs_hit": len(hit_accs),
            "n_accs_strong": sum(1 for a in hit_accs if per_acc[a]["strong"]),
            "n_accs_weak": sum(1 for a in hit_accs if per_acc[a]["weak"]),
            "status": status,
        })
    return genome, rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    a = ap.parse_args()
    by_genome = plan()
    print(f"genomas: {len(by_genome)}; brechas: "
          f"{sum(len(v) for v in by_genome.values())}")
    rows: list[dict] = []
    done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_scan_one, (g, s)): g for g, s in by_genome.items()}
        for fut in as_completed(futs):
            _, part = fut.result()
            rows.extend(part)
            done += 1
            if done % 10 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} genomas", file=sys.stderr)

    rows.sort(key=lambda r: (r["genome"], r["gap"]))
    fields = ["genome", "gap", "best_acc", "source", "qseqid", "pident",
              "qcov", "evalue", "bitscore", "n_hits", "n_accs_hit",
              "n_accs_strong", "n_accs_weak", "status"]
    write_tsv(WORK / "homology_hits.tsv", rows, fields)
    stats = {
        "genomes": len(by_genome), "rows": len(rows),
        "status": dict(Counter(r["status"] for r in rows)),
        "with_hit": sum(1 for r in rows if r["best_acc"]),
        "thresholds": {
            "strong": f"pident>={MIN_PIDENT_STRONG} qcov>={MIN_QCOV_STRONG} evalue<={MAX_EVALUE}",
            "weak": f"pident>={MIN_PIDENT_WEAK} qcov>={MIN_QCOV_WEAK} evalue<={MAX_EVALUE}",
        },
        "evalue_search": 1e-3, "max_target_seqs": 500,
        "sources": ["ncbi", "abinitio"],
    }
    write_json(WORK / "homology.json", stats)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


def _empty(genome: str, gap: str, status: str, note: str) -> dict:
    return {"genome": genome, "gap": gap, "best_acc": "", "source": "",
            "qseqid": "", "pident": "", "qcov": "", "evalue": "",
            "bitscore": "", "n_hits": 0, "n_accs_hit": 0, "n_accs_strong": 0,
            "n_accs_weak": 0, "status": status, "note": note}


if __name__ == "__main__":
    raise SystemExit(main())
