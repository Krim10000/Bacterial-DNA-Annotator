#!/usr/bin/env python3
"""Descarga la estructura 3D de una proteína usando la API de ESMFold."""
import requests, json
from pathlib import Path

def fetch_structure(sequence, out_file):
    print("Consultando ESMFold (esto tarda ~10 segundos)...")
    # API pública de ESMFold (Meta)
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    res = requests.post(url, data=sequence)
    
    if res.status_code == 200:
        with open(out_file, "w") as f:
            f.write(res.text)
        print(f"✅ Estructura 3D guardada en {out_file}")
    else:
        print(f"❌ Error: {res.status_code}")

def main():
    root = Path.cwd()
    out_dir = root / "data" / "structures"
    out_dir.mkdir(exist_ok=True)
    
    # Usamos una proteína de ejemplo (una toxina o enzima de tu bacteria)
    # En un caso real, esto vendría de tu archivo FAA.
    test_sequence = "MKKLLPTAALLLSACSTPAPAQETVKLQGSGKYTVTATGQTVTIKDLNGVQQVKATLDLQGQTVTLTGQTVTVKDLNG"
    fetch_structure(test_sequence, out_dir / "test_protein.pdb")

if __name__ == "__main__":
    main()
