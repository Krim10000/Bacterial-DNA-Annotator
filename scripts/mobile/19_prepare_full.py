#!/usr/bin/env python3
"""19_prepare_full.py — inventario y lista de los 400 genomas.

Comprueba entradas mínimas (GFF ab initio, proteins.faa, GFF NCBI, FASTA),
calcula contigs/tamaño, cruza metadatos de especie/filo y escribe
`data/mobile_typing/full_genomes.tsv` + `_work/full/state/inputs.json`.

Uso: python scripts/mobile/19_prepare_full.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANN = ROOT / "data" / "annotations_ab_initio"
NCBI = ROOT / "data" / "annotations"
RAW = ROOT / "data" / "raw"
BASE = ROOT / "data" / "mobile_typing"
STATE = BASE / "_work" / "full" / "state"


def contig_lengths(gff: Path):
    out = []
    if not gff.exists():
        return out
    for line in gff.open():
        if line.startswith("##sequence-region"):
            parts = line.split()
            out.append((parts[1], int(parts[3])))
        elif not line.startswith("#"):
            break
    return out


def load_metadata():
    meta = {}
    for name in ("species_phylum.tsv", "selected.tsv", "one_per_species.tsv"):
        path = ROOT / "data" / name
        if not path.exists():
            continue
        with path.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                meta.setdefault(row["assembly_accession"], {})
                for key in ("organism_name", "phylum", "species_taxid"):
                    if row.get(key):
                        meta[row["assembly_accession"]].setdefault(key, row[key])
    return meta


def main() -> int:
    STATE.mkdir(parents=True, exist_ok=True)
    meta = load_metadata()
    genomes = sorted(d.name for d in ANN.iterdir() if d.is_dir())
    rows, problems = [], {}
    for genome in genomes:
        acc = "_".join(genome.split("_")[:2])
        checks = {
            "gff_ab_initio": (ANN / genome / "annotation.gff3").exists(),
            "proteins": (ANN / genome / "proteins.faa").exists(),
            "gff_ncbi": (NCBI / f"{genome}.gff.gz").exists(),
            "fasta": (RAW / f"{genome}_genomic.fna.gz").exists(),
        }
        missing = [k for k, ok in checks.items() if not ok]
        contigs = contig_lengths(ANN / genome / "annotation.gff3")
        m = meta.get(acc, {})
        rows.append(
            {
                "genome": genome,
                "phylum": m.get("phylum", "unknown"),
                "organism": m.get("organism_name", "unknown"),
                "contigs": len(contigs),
                "genome_bp": sum(length for _, length in contigs),
                "has_fasta": "yes" if checks["fasta"] else "no",
                "missing_inputs": ",".join(missing),
            }
        )
        if missing:
            problems[genome] = missing
    out = BASE / "full_genomes.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["genome", "phylum", "organism", "contigs", "genome_bp",
                        "has_fasta", "missing_inputs"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    (STATE / "inputs.json").write_text(
        json.dumps({"genomes": len(rows), "problems": problems}, indent=2, sort_keys=True) + "\n"
    )
    ready = sum(1 for r in rows if not r["missing_inputs"])
    print(f"OK: {len(rows)} genomas, {ready} con entradas completas, {len(problems)} con faltantes")
    for genome, missing in list(problems.items())[:10]:
        print(f"  {genome}: faltan {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
