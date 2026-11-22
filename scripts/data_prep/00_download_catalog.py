#!/usr/bin/env python3
"""Descarga assembly_summary.txt y taxdump de NCBI."""
import sys, urllib.request, tarfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.io import read_config, setup_logging

URLS = {
    "assembly_summary.txt": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
    "taxdump.tar.gz": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
}
def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0: return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  ↓ {url}")
    urllib.request.urlretrieve(url, dest)
    return True

def main():
    cfg = read_config()
    log = setup_logging("00_download_catalog", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    data, taxdump = root / "data", root / cfg["data"]["taxdump_dir"]
    for name, url in URLS.items():
        target = taxdump / name if "taxdump" in name else data / name
        try:
            if download(url, target): log.info(f"Descargado: {target}")
        except Exception as e:
            log.error(f"Falló {url}: {e}"); sys.exit(1)
    tgz = taxdump / "taxdump.tar.gz"
    needed = ["nodes.dmp", "names.dmp", "merged.dmp"]
    if not all((taxdump / n).exists() for n in needed):
        with tarfile.open(tgz, "r:gz") as tar:
            for member in tar.getmembers():
                if member.name in needed: tar.extract(member, taxdump)
    log.info("Catálogo y taxdump listos.")

if __name__ == "__main__": main()
