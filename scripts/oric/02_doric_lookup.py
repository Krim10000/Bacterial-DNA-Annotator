#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_doric_lookup.py — Flujo ORIC (Paso 2).
Cruza los 400 genomas con la base de datos DoriC (formato v14).
Extrae las coordenadas reales de oriC (ground truth).
"""
import csv
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("doric")

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DORIC_DIR = PROJECT_ROOT / "data" / "doric"
OUT_FILE = PROJECT_ROOT / "data" / "doric_ground_truth.json"

def main():
    log.info("Buscando archivos DoriC en %s...", DORIC_DIR)
    files = list(DORIC_DIR.glob("*.csv")) + list(DORIC_DIR.glob("*.txt"))
    
    if not files:
        log.error("No se encontraron archivos. Descomprime DoriC en data/doric/")
        return
    
    doric_file = files[0]
    log.info("Usando archivo: %s", doric_file.name)
    
    gt = {}
    with open(doric_file, 'r', encoding='utf-8-sig') as f:
        first_line = f.readline()
        delimiter = ',' if ',' in first_line else '\t'
        f.seek(0)
        
        reader = csv.DictReader(f, delimiter=delimiter)
        
        count_genomes = 0
        for row in reader:
            refseq = row.get("refseq", "").strip()
            nc = row.get("nc", "").strip()
            orics_raw = row.get("orics", "")
            
            if not orics_raw or not (refseq or nc):
                continue
                
            try:
                # Parsear el string JSON a una lista de diccionarios
                orics_list = json.loads(orics_raw)
            except json.JSONDecodeError:
                continue
                
            # Guardar la lista completa de orígenes para este genoma
            info = {
                "organism": row.get("organism", ""),
                "orics": [
                    {
                        "start": int(o.get("oric_start", 0)),
                        "end": int(o.get("oric_end", 0)),
                        "length": int(o.get("oric_length", 0)),
                        "at_content": float(o.get("oric_at_content", 0.0))
                    } for o in orics_list
                ]
            }
            
            if refseq:
                gt[refseq] = info
            if nc:
                gt[nc] = info
                
            count_genomes += 1
            
    if not gt:
        log.warning("Se leyó el archivo pero no se extrajo ningún oriC.")
        return

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(json.dumps(gt, indent=2))
    log.info("✅ Guardados %d genomas (%d entradas con accesiones) en %s", 
             count_genomes, len(gt), OUT_FILE.name)

if __name__ == "__main__":
    main()
