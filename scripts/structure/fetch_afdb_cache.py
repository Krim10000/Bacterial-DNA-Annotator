#!/usr/bin/env python3
"""Caché AFDB/PAE por accesión única (reanudable, sin duplicados).

- Reúne accesiones UniProt de los genes de los 400.
- Excluye las que ya tienen PDB/PAE local.
- Descarga PDB y PAE JSON de AlphaFold DB.
- Catálogo TSV con hash, tamaño y estado.
- ``--limit`` para pilotos; ``--pae`` para bajar también PAE.
- No descarga dos veces la misma accesión.

Ejemplo:
  python scripts/structure/fetch_afdb_cache.py --limit 25 --pae
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
CAT_COLS = ["accession", "status", "pdb", "pae", "sha256", "pdb_bytes",
            "pae_bytes", "pdb_url", "pae_url", "error", "updated"]


def accessions(data_root: Path) -> list[str]:
    seen = []
    for path in sorted((data_root / "regulation").glob("*/genes.tsv")):
        try:
            with path.open(encoding="utf-8", errors="replace") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    acc = (row.get("acc") or "").strip()
                    if acc:
                        seen.append(acc)
        except OSError:
            continue
    return sorted(set(seen))


def load_catalog(path: Path) -> dict[str, dict]:
    out = {}
    if path.exists():
        with path.open(encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("accession"):
                    out[row["accession"]] = row
    return out


def write_catalog(path: Path, rows: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CAT_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for acc in sorted(rows):
            w.writerow(rows[acc])
    tmp.replace(path)


def fetch(acc: str, targets: Path, catalog: dict, timeout: int, retries: int,
          want_pae: bool, delay: float = 0.0) -> dict:
    pdb = targets / f"{acc}.pdb"
    pae = targets / f"{acc}.pae.json"
    existing = catalog.get(acc, {})
    if pdb.exists() and (not want_pae or pae.exists()):
        return {**existing, "accession": acc, "status": "already_present",
                "pdb": str(pdb), "pae": str(pae) if pae.exists() else ""}
    last = ""
    for attempt in range(1, retries + 1):
        try:
            if delay:
                time.sleep(delay)
            meta = json.load(urllib.request.urlopen(
                API.format(acc=acc), timeout=timeout))
            if not meta:
                raise ValueError("sin entrada AFDB")
            item = meta[0]
            pdb_url = item.get("pdbUrl") or ""
            pae_url = item.get("paeDocUrl") or ""
            if pdb_url:
                data = urllib.request.urlopen(pdb_url, timeout=timeout * 3).read()
                tmp = pdb.with_suffix(".pdb.tmp")
                tmp.write_bytes(data)
                tmp.replace(pdb)
            if want_pae and pae_url:
                data = urllib.request.urlopen(pae_url, timeout=timeout * 3).read()
                tmp = pae.with_suffix(".json.tmp")
                tmp.write_bytes(data)
                tmp.replace(pae)
            sha = hashlib.sha256(pdb.read_bytes()).hexdigest() if pdb.exists() else ""
            return {"accession": acc, "status": "ok", "pdb": str(pdb),
                    "pae": str(pae) if pae.exists() else "",
                    "sha256": sha, "pdb_bytes": pdb.stat().st_size if pdb.exists() else 0,
                    "pae_bytes": pae.stat().st_size if pae.exists() else 0,
                    "pdb_url": pdb_url, "pae_url": pae_url, "error": "",
                    "updated": time.strftime("%Y-%m-%dT%H:%M:%S")}
        except urllib.error.HTTPError as exc:
            last = f"HTTPError: HTTP {exc.code}: {exc.reason}"
            if exc.code == 429:
                # AFDB limita por IP. Esperar más y reintentar; no tratarlo
                # como un 404 permanente.
                wait = int(exc.headers.get("Retry-After") or 0) or min(120, 20 * attempt)
                time.sleep(wait)
            elif exc.code == 404:
                return {"accession": acc, "status": "not_found", "pdb": "",
                        "pae": "", "sha256": "", "pdb_bytes": 0, "pae_bytes": 0,
                        "pdb_url": "", "pae_url": "", "error": last,
                        "updated": time.strftime("%Y-%m-%dT%H:%M:%S")}
            else:
                time.sleep(min(20, attempt * 3))
        except Exception as exc:  # noqa: BLE001
            last = f"{type(exc).__name__}: {exc}"
            time.sleep(min(20, attempt * 3))
    return {"accession": acc, "status": "error", "pdb": "", "pae": "",
            "sha256": "", "pdb_bytes": 0, "pae_bytes": 0, "pdb_url": "",
            "pae_url": "", "error": last, "updated": time.strftime("%Y-%m-%dT%H:%M:%S")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--targets", type=Path,
                    default=ROOT / "data" / "metabolism" / "targets" / "afdb25k")
    ap.add_argument("--catalog", type=Path,
                    default=ROOT / "data" / "metabolism" / "targets" / "afdb_catalog.tsv")
    ap.add_argument("--limit", type=int, default=0, help="0 = todas")
    ap.add_argument("--pae", action="store_true", help="descargar PAE JSON")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--delay", type=float, default=0.0,
                    help="segundos entre solicitudes por worker; usar >0 con 429")
    ap.add_argument("--timeout", type=int, default=25)
    ap.add_argument("--retries", type=int, default=2)
    a = ap.parse_args()
    a.targets.mkdir(parents=True, exist_ok=True)
    cat = load_catalog(a.catalog)
    accs = accessions(a.data_root)
    def _permanent(acc: str) -> bool:
        row = cat.get(acc) or {}
        return row.get("status") == "not_found" or "404" in str(row.get("error", ""))
    pending = [x for x in accs if not (
        (a.targets / f"{x}.pdb").exists()
        and (not a.pae or (a.targets / f"{x}.pae.json").exists()))
        and not _permanent(x)]
    if a.limit:
        pending = pending[:a.limit]
    print(f"accesiones totales={len(accs)} pendientes={len(pending)} "
          f"pae={a.pae} workers={a.workers}")
    done = 0
    with ThreadPoolExecutor(max_workers=max(1, a.workers)) as ex:
        futures = {ex.submit(fetch, acc, a.targets, cat, a.timeout,
                             a.retries, a.pae, a.delay): acc for acc in pending}
        for fut in as_completed(futures):
            row = fut.result()
            cat[row["accession"]] = row
            done += 1
            if done % 5 == 0 or done == len(pending):
                write_catalog(a.catalog, cat)
                print(f"  {done}/{len(pending)}")
    write_catalog(a.catalog, cat)
    ok = sum(1 for r in cat.values() if r.get("status") == "ok")
    err = sum(1 for r in cat.values() if r.get("status") == "error")
    print(f"catálogo={a.catalog} ok={ok} error={err}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
