#!/usr/bin/env python3
"""Descarga los recursos que la app necesita la primera vez.

    python -m webapp.setup_data --data-root data
    python -m webapp.setup_data --data-root data --pfam     # + Pfam (~600 MB)

Recursos:
  - Swiss-Prot FASTA + base DIAMOND  -> data/db/swissprot.dmnd
  - Rhea RDF y mapeos               -> data/metabolism/, data/reactions/tsv/
  - Pfam-A (opcional)               -> db/pfam/Pfam-A.hmm (presionado)
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

URLS = {
    "swissprot": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    "rhea_rdf": "https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz",
    "rhea2uniprot": "https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv",
    "rhea2ec": "https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv",
    "pfam": "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
}


def _progress(blocks: int, block_size: int, total: int) -> None:
    if total <= 0:
        return
    pct = min(100, blocks * block_size * 100 // total)
    print(f"\r  {pct:3d}%", end="", flush=True)


def download(url: str, dest: Path, force: bool = False) -> None:
    if dest.exists() and not force:
        print(f"  ya existe: {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  descargando {url}")
    urllib.request.urlretrieve(url, dest, reporthook=_progress)
    print()


def ensure_rhea(data_root: Path, force: bool = False) -> None:
    tsV = data_root / "metabolism" / "tsv"
    download(URLS["rhea_rdf"], tsV / "rhea.rdf.gz", force)
    download(URLS["rhea2uniprot"], data_root / "metabolism" / "rhea2uniprot_sprot.tsv", force)
    download(URLS["rhea2ec"], data_root / "reactions" / "tsv" / "rhea2ec.tsv", force)


def ensure_swissprot(data_root: Path, force: bool = False) -> None:
    if shutil.which("diamond") is None:
        print("ERROR: falta el binario 'diamond' (apt install diamond-aligner).")
        sys.exit(1)
    db = data_root / "db" / "swissprot.dmnd"
    if db.exists() and not force:
        print(f"  ya existe: {db}")
        return
    gz = data_root / "db" / "uniprot_sprot.fasta.gz"
    fasta = data_root / "db" / "uniprot_sprot.fasta"
    download(URLS["swissprot"], gz, force)
    print("  descomprimiendo…")
    with gzip.open(gz, "rb") as src, open(fasta, "wb") as dst:
        shutil.copyfileobj(src, dst)
    print("  creando base DIAMOND…")
    subprocess.run(["diamond", "makedb", "--in", str(fasta),
                    "-d", str(db.with_suffix("")), "--quiet"], check=True)


def ensure_pfam(pfam_dir: Path, force: bool = False) -> None:
    if shutil.which("hmmpress") is None:
        print("ERROR: falta 'hmmpress' (apt install hmmer).")
        sys.exit(1)
    gz = pfam_dir / "Pfam-A.hmm.gz"
    hmm = pfam_dir / "Pfam-A.hmm"
    download(URLS["pfam"], gz, force)
    if not hmm.exists() or force:
        print("  descomprimiendo Pfam…")
        with gzip.open(gz, "rb") as src, open(hmm, "wb") as dst:
            shutil.copyfileobj(src, dst)
    print("  presionando Pfam…")
    subprocess.run(["hmmpress", "-f", str(hmm)], check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--pfam", action="store_true", help="descargar también Pfam-A (~600 MB)")
    ap.add_argument("--only", choices=["all", "rhea", "swissprot"], default="all")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    data_root = a.data_root
    data_root.mkdir(parents=True, exist_ok=True)
    print("Recursos del anotador")
    if a.only in ("all", "swissprot"):
        print("\n[1/2] Swiss-Prot + DIAMOND")
        ensure_swissprot(data_root, a.force)
    if a.only in ("all", "rhea"):
        print("\n[2/2] Rhea")
        ensure_rhea(data_root, a.force)
    if a.pfam:
        print("\n[+] Pfam")
        ensure_pfam(Path("db") / "pfam", a.force)
    print("\nListo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
