#!/usr/bin/env python3
"""Muestreo estratificado por filo, reproducible con seed."""
import sys, random
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.io import read_config, setup_logging

def main():
    cfg = read_config()
    log = setup_logging("03_stratified_sample", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    data = root / "data"
    sel = cfg["selection"]
    cap = sel["cap_per_phylum"]
    random.seed(sel["seed"])

    inp = data / "species_phylum.tsv"
    if not inp.exists():
        log.error(f"No existe {inp}. Corre 02 primero.")
        sys.exit(1)

    by_phylum = {}
    with open(inp) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            phylum = parts[4] or "unclassified"
            by_phylum.setdefault(phylum, []).append(parts)

    selected = []
    print(f"Cap por filo: {cap}")
    for phylum, rows in by_phylum.items():
        random.shuffle(rows)
        take = min(len(rows), cap)
        selected.extend(rows[:take])
        print(f"  {phylum:30s} {len(rows):5d} disponibles -> {take:3d}")

    out = data / "selected.tsv"
    with open(out, "w") as f:
        f.write("assembly_accession\tspecies_taxid\torganism_name\tftp_path\tphylum\n")
        for r in selected:
            f.write("\t".join(r) + "\n")

    log.info(f"Seleccionados: {len(selected)} genomas de {len(by_phylum)} filos")
    print(f"\nTotal: {len(selected)} genomas -> {out}")

if __name__ == "__main__":
    main()
