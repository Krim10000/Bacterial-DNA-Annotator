#!/usr/bin/env python3
"""
Orquestador del Flujo 03_oric.
Ejecuta: validación de marcadores -> DoriC lookup -> Predicción ORCA.
"""
import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("oric_orch")

HERE = Path(__file__).resolve().parent

STEPS = [
    ("01_validate_markers.py", "Validando marcadores dnaA/dnaN"),
    ("02_doric_lookup.py", "Cruzando con base de datos DoriC"),
    ("03_orca_predict.py", "Ejecutando ORCA (ML)"),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="Pasar --force a los scripts")
    args = ap.parse_args()

    for script_name, desc in STEPS:
        log.info(f"▶ {desc} ({script_name})...")
        cmd = [sys.executable, str(HERE / script_name)]
        if args.force and script_name == "03_orca_predict.py":
            cmd.append("--force")
            
        result = subprocess.run(cmd)
        if result.returncode != 0:
            log.error(f"❌ Falló el paso {script_name}.")
            sys.exit(1)
            
    log.info("✅ Flujo 03_oric completado.")

if __name__ == "__main__":
    main()
