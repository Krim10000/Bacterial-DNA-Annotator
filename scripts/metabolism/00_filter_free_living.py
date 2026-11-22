#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_filter_free_living.py — Filtra genomas de vida libre.
Criterio: Genomas > 2.0 Mb (excluye parásitos obligados y simbiontes).
"""
import gzip
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("filter")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"
OUT_FILE = PROJECT_ROOT / "data" / "metabolism" / "free_living_genomes.json"
MIN_SIZE_MB = 2.0

def main():
    free_living = []
    excluded = []
    
    for fna in sorted(RAW_DIR.glob("*_genomic.fna.gz")):
        prefix = fna.name.replace("_genomic.fna.gz", "")
        
        # Calcular tamaño del genoma leyendo el FASTA
        size_bp = 0
        with gzip.open(fna, 'rt') as f:
            for line in f:
                if not line.startswith('>'):
                    size_bp += len(line.strip())
                    
        size_mb = size_bp / 1_000_000
        
        if size_mb >= MIN_SIZE_MB:
            free_living.append(prefix)
        else:
            excluded.append(prefix)
            
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(json.dumps(free_living, indent=2))
    
    log.info(f"Total evaluados: {len(free_living) + len(excluded)}")
    log.info(f"✅ Vida libre (> {MIN_SIZE_MB} Mb): {len(free_living)} genomas guardados en {OUT_FILE.name}")
    log.info(f"❌ Excluidos (< {MIN_SIZE_MB} Mb, parásitos/simbiontes): {len(excluded)}")

if __name__ == "__main__":
    main()
