#!/usr/bin/env python3
"""Descarga los archivos GFF3 y FAA de NCBI para nuestras 400 bacterias."""
import sys, urllib.request
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # scripts/ (lib)
sys.path.insert(0, str(Path(__file__).resolve().parent))        # esta carpeta
from lib.io import read_config, setup_logging

def download_one(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return "exists"
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "dna_poc/0.1"})
        with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
            f.write(r.read())
        return "ok"
    except Exception as e:
        return f"fail: {e}"

def main():
    cfg = read_config()
    root = Path.cwd()
    data_dir = root / "data"
    ann_dir = data_dir / "annotations"
    ann_dir.mkdir(exist_ok=True)

    sel_file = data_dir / "selected.tsv"
    if not sel_file.exists():
        print(f"No existe {sel_file}")
        return

    ok, exists, fail = 0, 0, 0
    with open(sel_file) as f:
        f.readline() # Header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4: continue
            ftp_path = parts[3]
            basename = ftp_path.rstrip("/").split("/")[-1]
            
            # NCBI usa .gff.gz y .faa.gz para las anotaciones
            gff_url = f"{ftp_path}/{basename}_genomic.gff.gz"
            faa_url = f"{ftp_path}/{basename}_protein.faa.gz"
            
            gff_dest = ann_dir / f"{basename}.gff.gz"
            faa_dest = ann_dir / f"{basename}.faa.gz"
            
            status_gff = download_one(gff_url, gff_dest)
            status_faa = download_one(faa_url, faa_dest)
            
            if status_gff == "ok" or status_faa == "ok":
                ok += 1
                print(f"✓ Descargado: {basename}")
            elif status_gff == "exists":
                exists += 1
            else:
                fail += 1
                print(f"✗ Falló: {basename}")

    print(f"\nResumen: {ok} nuevos, {exists} ya existían, {fail} fallidos.")
    print(f"Anotaciones guardadas en: {ann_dir}")

if __name__ == "__main__":
    main()
