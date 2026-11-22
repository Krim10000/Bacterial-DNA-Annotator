#!/usr/bin/env python3
"""
Orquestador del Flujo 00_data_prep.
Ejecuta: catálogo -> filtro -> filo -> muestreo -> descarga.
Por defecto NO sobreescribe si los archivos ya existen (usar --force).
"""
import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("prep_orch")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
HERE = Path(__file__).resolve().parent

# (Script, Archivo de salida esperado, Descripción)
STEPS = [
    ("00_download_catalog.py", "data/assembly_summary.txt", "Descargando catálogo NCBI"),
    ("01_filter_candidates.py", "data/one_per_species.tsv", "Filtrando candidatos (excluyendo parásitos)"),
    ("02_annotate_phylum.py", "data/species_phylum.tsv", "Anotando filo taxonómico"),
    ("03_stratified_sample.py", "data/selected.tsv", "Muestreo estratificado de 400 especies"),
    ("04_download_fasta.py", "data/raw", "Descargando genomas FASTA")
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="Sobreescribir archivos existentes")
    args = ap.parse_args()

    for script_name, expected_out, desc in STEPS:
        out_path = PROJECT_ROOT / expected_out
        script_path = HERE / script_name

        # Lógica anti-sobreescritura
        if not args.force:
            if out_path.exists() and out_path.stat().st_size > 0:
                log.info(f"⏭️  Saltando {script_name} (El archivo {expected_out} ya existe).")
                continue
        
        log.info(f"▶ {desc} ({script_name})...")
        cmd = [sys.executable, str(script_path)]
        
        result = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
        if result.returncode != 0:
            log.error(f"❌ Falló el paso {script_name}.")
            sys.exit(1)
            
    log.info("✅ Flujo 00_data_prep completado.")

if __name__ == "__main__":
    main()
