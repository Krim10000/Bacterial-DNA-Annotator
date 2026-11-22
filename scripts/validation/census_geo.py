#!/usr/bin/env python3
"""census_geo.py — ¿qué especies del proyecto tienen datos de transcripción?

Consulta GEO (NCBI E-utilities) por cada binomio de especie representado en los
genomas anotados y cuenta estudios de dRNA-seq (TSS) y Term-seq (TTS).

Resultado del censo de 2026-11-21: **9 de 368 especies (2%)** tienen dRNA-seq.
Es el dato que justifica no intentar poblar la transcripción de los 400 por
transferencia: la información no existe.

Salida: data/operon_finder/geo_census.tsv
Requiere red. Con --limit se prueba un subconjunto.
"""
from __future__ import annotations
import argparse, csv, json, time, urllib.parse, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def count(term: str) -> int:
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmode": "json", "retmax": "1"})
    try:
        with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=25) as fh:
            return int(json.load(fh)["esearchresult"]["count"])
    except Exception:
        return -1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--delay", type=float, default=0.35)
    a = ap.parse_args()

    acc2org = {}
    summary = ROOT / "data/assembly_summary.txt"
    for line in open(summary, encoding="utf-8", errors="replace"):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 7:
            acc2org[f[0]] = f[7]
    bins = set()
    for d in (ROOT / "data/annotations_ab_initio").iterdir():
        if not d.is_dir():
            continue
        org = acc2org.get(d.name.split("_ASM")[0], "").split()
        if len(org) >= 2:
            bins.add(" ".join(org[:2]))
    bins = sorted(bins)
    if a.limit:
        bins = bins[:a.limit]
    print(f"binomios a consultar: {len(bins)}")
    rows = []
    for i, s in enumerate(bins, 1):
        rows.append({"species": s,
                     "drna": count(f'differential RNA-seq AND "{s}"'),
                     "termseq": count(f'Term-seq AND "{s}"')})
        if i % 50 == 0:
            print(f"  {i}/{len(bins)}", flush=True)
        time.sleep(a.delay)
    out = ROOT / "data/operon_finder/geo_census.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    for k in ("drna", "termseq"):
        n = sum(1 for r in rows if r[k] > 0)
        print(f"con {k}>0: {n}/{len(rows)} ({100*n/max(1,len(rows)):.0f}%)")
    print(f"-> {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
