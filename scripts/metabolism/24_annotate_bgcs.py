#!/usr/bin/env python3
"""Busca clústeres de biosíntesis (BGCs) en los GFF."""
import gzip, re, json
from pathlib import Path

BGC_KEYWORDS = [
    "polyketide synthase", "pks", "non-ribosomal peptide", "nrps", 
    "bacteriocin", "terpene", "lantibiotic", "siderophore"
]

def find_bgcs(gff_path):
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
                    if any(kw in prod for kw in BGC_KEYWORDS):
                        elements.append({"start": int(parts[3]), "end": int(parts[4]), "name": match.group(1)})
    except Exception: pass
    return elements

def main():
    root = Path.cwd()
    out_file = root / "data" / "bgcs.json"
    gff_files = list((root / "data" / "annotations").glob("*.gff.gz"))
    
    print(f"Buscando fábricas de antibióticos en {len(gff_files)} genomas...")
    global_db = {}
    for gff in gff_files:
        gid = gff.stem.replace(".gff", "")
        elems = find_bgcs(gff)
        if elems: global_db[gid] = elems

    with open(out_file, "w") as f:
        json.dump(global_db, f, indent=2)
    print(f"✅ Guardado en {out_file}. Se encontraron {sum(len(v) for v in global_db.values())} genes en BGCs.")

if __name__ == "__main__":
    main()
