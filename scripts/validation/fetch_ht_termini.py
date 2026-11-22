#!/usr/bin/env python3
"""fetch_ht_termini.py — extremos experimentales de transcripción de RegulonDB.

RegulonDB publica datasets de alto rendimiento con:
  - **TUS**: unidades de transcripción determinadas experimentalmente (genes,
    extremos 5' y 3', tipo de terminación).
  - **TTS**: sitios de terminación (extremo 3' real).
  - **TSS**: sitios de inicio (extremo 5' real).

Es el dato que no se puede derivar de la secuencia y que la estructura
secundaria no sustituye (`compare_structure.py`: ViennaRNA AUC 0.829 vs proxy
0.852, y los terminadores dependientes de Rho no tienen horquilla).

Salidas en `data/regulondb/ht/`:
    ht_tus.tsv   unidades de transcripción experimentales
    ht_tts.tsv   extremos 3'
    ht_tss.tsv   extremos 5'

Uso
---
    python scripts/validation/fetch_ht_termini.py
    python scripts/validation/fetch_ht_termini.py --offline   # usa caché
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "regulondb" / "ht"
URL = "https://regulondb.ccg.unam.mx/graphql"

KINDS = {
    "tus": ("TUS", "getAllTransUnitsOfDataset", "HtTranscriptionUnit",
            ["_id", "chromosome", "leftEndPosition", "rightEndPosition", "name",
             "strand", "length", "termType", "phantom", "pseudo", "temporalId"]),
    "tts": ("TTS", "getAllTTSOfDataset", "HTtranscriptionTerminationSite",
            ["_id", "chromosome", "leftEndPosition", "rightEndPosition", "name",
             "strand", "temporalId"]),
    "tss": ("TSS", "getAllTSSOfDataset", "HtTSS",
            ["_id", "chromosome", "leftEndPosition", "rightEndPosition", "pos_1",
             "strand"]),
}


def post(query: str) -> dict:
    req = urllib.request.Request(
        URL, method="POST", data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json"})
    try:
        return json.loads(urllib.request.urlopen(req, timeout=120).read())
    except urllib.error.HTTPError as e:
        return {"http": e.code, "body": e.read().decode()[:300]}


def datasets(dtype: str) -> list[str]:
    p = post('{getDatasetsWithMetadata(datasetType:"%s", source:"REGULONDB")'
             '{datasets{_id}}}' % dtype)
    ds = (p.get("data") or {}).get("getDatasetsWithMetadata", {}).get("datasets", [])
    return [d["_id"] for d in ds if d.get("_id")]


def fetch_all(kind: str) -> list[dict]:
    dtype, field, typename, cols = KINDS[kind]
    rows: list[dict] = []
    for ds in datasets(dtype):
        page = 1
        while True:
            q = ('{%s(datasetId:"%s", limit:500, page:%d){%s}}'
                 % (field, ds, page, " ".join(cols)))
            p = post(q)
            if "http" in p:
                print(f"  {kind}/{ds}: {p['http']} {p.get('body','')[:120]}", file=sys.stderr)
                break
            batch = (p.get("data") or {}).get(field) or []
            for b in batch:
                b["dataset_id"] = ds
            rows.extend(batch)
            if len(batch) < 500:
                break
            page += 1
            time.sleep(0.2)
        print(f"  {kind} {ds}: {len(rows)} acumuladas", file=sys.stderr)
    return rows


def write(rows: list[dict], name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    if not rows:
        path.write_text("")
        return path
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args()

    summary = {}
    for kind in ("tus", "tts", "tss"):
        path = OUT / f"ht_{kind}.tsv"
        if a.offline and path.exists():
            n = sum(1 for _ in open(path)) - 1
            summary[kind] = n
            print(f"{kind}: {n} (caché)")
            continue
        print(f"== {kind} ==", file=sys.stderr)
        rows = fetch_all(kind)
        path = write(rows, f"ht_{kind}.tsv")
        summary[kind] = len(rows)
        print(f"{kind}: {len(rows)} -> {path.relative_to(ROOT)}")

    print("\nresumen:", json.dumps(summary))
    return 0 if any(summary.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
