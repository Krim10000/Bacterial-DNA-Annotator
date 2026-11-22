#!/usr/bin/env python3
"""Capa pesada: pliega el top de candidatos y busca homología estructural.

- Plegado con ESMFold (API de ESM Atlas). Devuelve PDB y pLDDT medio.
- Búsqueda con Foldseek contra la base de estructuras Swiss-Prot.

Uso:
    python3 -m ghost.structure --genome <GCF> --out-dir ./out --top 5 \
        --foldseek-db /ruta/data/db/foldseek_swissprot
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def fold_esmfold(seq: str, retries: int = 3, timeout: int = 300) -> str:
    """Devuelve el PDB. Lanza RuntimeError si la API no responde."""
    last = ""
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(
            ESMFOLD_URL, data=seq.encode(), method="POST",
            headers={"Content-Type": "text/plain"},
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read().decode()
            if body.startswith("ATOM") or body.startswith("HEADER"):
                return body
            last = body[:200]
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = str(exc)[:200]
        time.sleep(5 * attempt)
    raise RuntimeError(f"ESMFold no respondió: {last}")


def mean_plddt(pdb: str) -> float:
    vals = []
    for line in pdb.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return round(statistics.mean(vals), 2) if vals else 0.0


def foldseek_search(pdb: Path, db: Path, out: Path, threads: int = 4) -> list[dict]:
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [
                a.foldseek_bin, "easy-search", str(pdb), str(db), str(out), tmp,
                "--format-output", "query,target,pident,alnlen,evalue,bits",
                "--threads", str(threads), "--quiet",
            ],
            check=True,
        )
    rows = []
    with open(out) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6:
                rows.append(dict(zip(
                    ["query", "target", "pident", "alnlen", "evalue", "bits"], p
                )))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--foldseek-db", type=Path, required=True)
    ap.add_argument("--foldseek-bin", default="foldseek")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--skip-foldseek", action="store_true")
    a = ap.parse_args()

    cand_path = a.out_dir / f"{a.genome}_ghost_candidates.json"
    doc = json.load(open(cand_path))
    cands = sorted(doc["candidates"], key=lambda c: c["rank"], reverse=True)[: a.top]

    struct_dir = a.out_dir / "structures"
    struct_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for c in cands:
        locus = f"{c['start']}-{c['end']}{c['strand']}"
        pdb_path = struct_dir / f"{a.genome}_{locus}.pdb"
        try:
            if not pdb_path.exists():
                pdb = fold_esmfold(c["aa"])
                pdb_path.write_text(pdb)
            pdb = pdb_path.read_text()
        except RuntimeError as exc:
            print(f"  {locus}: plegado FALLÓ ({exc})")
            continue
        plddt = mean_plddt(pdb)
        hits: list[dict] = []
        if not a.skip_foldseek:
            try:
                hits = foldseek_search(
                    pdb_path, a.foldseek_db,
                    struct_dir / f"{a.genome}_{locus}_foldseek.tsv",
                    a.threads,
                )
            except subprocess.CalledProcessError as exc:
                print(f"  {locus}: foldseek FALLÓ ({exc})")
        results.append({
            "locus": locus, "aa_len": c["aa_len"], "rank": c["rank"],
            "esm_llr": c.get("esm_llr"), "plddt": plddt,
            "n_struct_hits": len(hits),
            "best_hit": hits[0]["target"] if hits else "",
            "best_evalue": hits[0]["evalue"] if hits else "",
        })
        print(f"  {locus}: pLDDT={plddt} hits={len(hits)} {hits[0]['target'] if hits else ''}")

    out_tsv = a.out_dir / f"{a.genome}_structure_summary.tsv"
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()) if results else ["locus"], delimiter="\t")
        w.writeheader()
        w.writerows(results)
    print(f"  -> {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
