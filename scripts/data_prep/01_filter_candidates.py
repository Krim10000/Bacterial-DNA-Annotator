#!/usr/bin/env python3
"""Filtra Complete Genome + representative. EXCLUYE parásitos intracelulares obligados."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lib.io import read_config, setup_logging

# Lista de géneros que arruinan el dataset de vida libre (genomas reducidos)
EXCLUDED_GENERA = {
    "Mycoplasma", "Ureaplasma", "Chlamydia", "Chlamydophila", "Rickettsia", 
    "Orientia", "Ehrlichia", "Anaplasma", "Buchnera", "Wigglesworthia", 
    "Blochmannia", "Sulcia", "Hodgkinia", "Tremblaya", "Carsonella", 
    "Nasuia", "Zinderia", "Kinetoplastibacterium", "Candidatus", "Symbiobacterium"
}

def main():
    cfg = read_config()
    log = setup_logging("01_filter_candidates", cfg["logging"]["dir"])
    root = Path(cfg["project_root"])
    data = root / "data"
    summary = data / "assembly_summary.txt"
    if not summary.exists():
        log.error(f"No existe {summary}. Corre 00 primero."); sys.exit(1)

    sel = cfg["selection"]
    wanted_level, wanted_cats = sel["assembly_level"], set(sel["refseq_categories"])

    candidates, seen_species = [], set()
    excluded_count = 0
    
    with open(summary) as f:
        f.readline()
        cols = f.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(cols)}
        for line in f:
            row = line.rstrip("\n").split("\t")
            if len(row) <= idx["ftp_path"]: continue
            if row[idx["assembly_level"]] != wanted_level: continue
            if row[idx["refseq_category"]] not in wanted_cats: continue
            
            # FILTRO DE VIDA LIBRE: Excluir parásitos por nombre
            org_name = row[idx["organism_name"]]
            genus = org_name.split()[0]
            if genus in EXCLUDED_GENERA or any(ex in org_name for ex in EXCLUDED_GENERA):
                excluded_count += 1
                continue

            sp = row[idx["species_taxid"]]
            if sp in seen_species: continue
            seen_species.add(sp)
            candidates.append({
                "assembly_accession": row[idx["assembly_accession"]],
                "species_taxid": sp,
                "organism_name": org_name,
                "ftp_path": row[idx["ftp_path"]],
            })

    out = data / "one_per_species.tsv"
    with open(out, "w") as f:
        f.write("assembly_accession\tspecies_taxid\torganism_name\tftp_path\n")
        for c in candidates:
            f.write(f"{c['assembly_accession']}\t{c['species_taxid']}\t{c['organism_name']}\t{c['ftp_path']}\n")

    log.info(f"Candidatos (vida libre): {len(candidates)}")
    log.info(f"Excluidos por ser parásitos/simbiontes: {excluded_count}")

if __name__ == "__main__": main()
