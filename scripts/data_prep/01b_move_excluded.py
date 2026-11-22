#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
01b_move_excluded.py — Aísla los genomas excluidos (parásitos/simbiontes).
Mueve los archivos que NO están en free_living_genomes.json a data/raw_excluded/.
No borra nada y no sobreescribe.
"""
import json
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RAW_DIR = PROJECT_ROOT / "data" / "raw"
EXCLUDED_DIR = PROJECT_ROOT / "data" / "raw_excluded"
FREE_JSON = PROJECT_ROOT / "data" / "metabolism" / "free_living_genomes.json"

def main():
    if not FREE_JSON.exists():
        print(f"Error: No se encontró {FREE_JSON}. Corre 00_filter_free_living.py primero.")
        return
        
    free_living = set(json.loads(FREE_JSON.read_text()))
    EXCLUDED_DIR.mkdir(parents=True, exist_ok=True)
    
    moved = 0
    kept = 0
    
    for fna in RAW_DIR.glob("*_genomic.fna.gz"):
        prefix = fna.name.replace("_genomic.fna.gz", "")
        if prefix not in free_living:
            dest = EXCLUDED_DIR / fna.name
            shutil.move(str(fna), str(dest))
            moved += 1
        else:
            kept += 1
            
    print(f"✅ Genomas de vida libre mantenidos en data/raw/: {kept}")
    print(f"📦 Genomas parásitos/simbiontes movidos a data/raw_excluded/: {moved}")

if __name__ == "__main__":
    main()
