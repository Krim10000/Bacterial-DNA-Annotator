#!/usr/bin/env python3
"""Descarga los .fna.gz desde ftp_path (Python puro, con reintentos)."""
import sys, time, urllib.request, urllib.error
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.io import read_config, setup_logging

def download_one(url, dest, retries=2, timeout=120):
    if dest.exists() and dest.stat().st_size > 0:
        return "exists"
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "dna_poc/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(dest, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return "ok"
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < retries:
                time.sleep(2 ** attempt)
                continue
            return f"fail: {e}"
    return "fail"

def main():
    cfg = read_config()
    log = setup_logging("04_download_fasta", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    data = root / "data"
    raw = data / "raw"
    raw.mkdir(parents=True, exist_ok=True)

    inp = data / "selected.tsv"
    if not inp.exists():
        log.error(f"No existe {inp}. Corre 03 primero.")
        sys.exit(1)

    ok = exists = fail = 0
    with open(inp) as f:
        header = f.readline()
        for i, line in enumerate(f, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            acc, _, org, ftp, phylum = parts[:5]
            basename = ftp.rstrip("/").split("/")[-1]
            fname = f"{basename}_genomic.fna.gz"
            url = f"{ftp}/{fname}"
            dest = raw / fname
            status = download_one(url, dest)
            if status == "ok":
                ok += 1
                print(f"[{i}] ✓ {org[:60]} ({phylum})")
            elif status == "exists":
                exists += 1
                print(f"[{i}] = {fname}")
            else:
                fail += 1
                print(f"[{i}] ✗ {status} → {url}")

    log.info(f"Descargas: {ok} nuevas, {exists} ya existían, {fail} fallidas")
    print(f"\nResumen: {ok} ok, {exists} ya existían, {fail} fallidas")

if __name__ == "__main__":
    main()
