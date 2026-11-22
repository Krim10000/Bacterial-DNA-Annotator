#!/usr/bin/env python3
"""Mapea species_taxid -> filo parseando taxdump.dmp directamente."""
import sys
from pathlib import Path
from collections import Counter
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.io import read_config, setup_logging

def load_taxdump(taxdump_dir):
    parent, rank = {}, {}
    nodes = taxdump_dir / "nodes.dmp"
    with open(nodes) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            tid, pid, r = int(parts[0]), int(parts[1]), parts[2]
            parent[tid] = pid
            rank[tid] = r
    names = {}
    names_path = taxdump_dir / "names.dmp"
    with open(names_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            tid, name, cls = int(parts[0]), parts[1], parts[3]
            if cls == "scientific name":
                names[tid] = name
    return parent, rank, names

def find_phylum(taxid, parent, rank, names):
    t = taxid
    seen = set()
    while t in parent and t != 1 and t not in seen:
        seen.add(t)
        if rank.get(t) == "phylum":
            return names.get(t, "unknown")
        t = parent[t]
    return "unknown"

def main():
    cfg = read_config()
    log = setup_logging("02_annotate_phylum", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    data = root / "data"
    taxdump = root / cfg["data"]["taxdump_dir"]

    inp = data / "one_per_species.tsv"
    if not inp.exists():
        log.error(f"No existe {inp}. Corre 01 primero.")
        sys.exit(1)

    print("Cargando taxdump (puede tardar 30-60s)...")
    parent, rank, names = load_taxdump(taxdump)
    log.info(f"Taxdump cargado: {len(parent):,} nodos")

    out = data / "species_phylum.tsv"
    counts = Counter()
    with open(inp) as f, open(out, "w") as o:
        f.readline()  # header
        o.write("assembly_accession\tspecies_taxid\torganism_name\tftp_path\tphylum\n")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            acc, sp_taxid, org, ftp = parts[:4]
            phylum = find_phylum(int(sp_taxid), parent, rank, names)
            o.write(f"{acc}\t{sp_taxid}\t{org}\t{ftp}\t{phylum}\n")
            counts[phylum] += 1

    log.info(f"Filos encontrados: {len(counts)}")
    print(f"\nDistribucion por filo (top 20):")
    for phylum, n in counts.most_common(20):
        print(f"  {n:5d}  {phylum}")

if __name__ == "__main__":
    main()
