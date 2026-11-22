#!/usr/bin/env python3
"""Busca transposasas y elementos móviles en los GFF de NCBI."""
import gzip, re, json
from pathlib import Path

MOBILE_KEYWORDS = [
    "transposase", "insertion sequence", "is element", "reverse transcriptase",
    "integrase", "resolvase", "relaxase", "conjugative"
]

def find_mobile_elements(gff_path):
    elements = []
    opener = gzip.open if str(gff_path).endswith('.gz') else open
    try:
        with opener(gff_path, 'rt') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS": continue
                
                match = re.search(r'product=([^;]+)', parts[8])
                if match:
                    prod = match.group(1).lower()
                    if any(kw in prod for kw in MOBILE_KEYWORDS):
                        elements.append({"start": int(parts[3]), "end": int(parts[4]), "name": match.group(1)})
    except Exception:
        pass
    return elements

def main():
    root = Path.cwd()
    out_file = root / "data" / "mobile_elements.json"
    gff_files = list((root / "data" / "annotations").glob("*.gff.gz"))
    
    print(f"Buscando elementos móviles en {len(gff_files)} genomas...")
    global_db = {}
    
    for i, gff in enumerate(gff_files):
        gid = gff.stem.replace(".gff", "")
        elems = find_mobile_elements(gff)
        if elems: global_db[gid] = elems
        if (i+1) % 100 == 0: print(f"  [{i+1}/{len(gff_files)}]...")

    with open(out_file, "w") as f:
        json.dump(global_db, f, indent=2)
    print(f"✅ Guardado en {out_file}. Se encontraron {sum(len(v) for v in global_db.values())} elementos.")

if __name__ == "__main__":
    main()
