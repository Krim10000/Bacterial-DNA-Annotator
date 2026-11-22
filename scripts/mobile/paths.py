"""Rutas del módulo móvil: piloto (10) o barrido completo (400).

El modo se activa con `--full` en la línea de comandos. El piloto conserva
exactamente sus rutas y nombres para no romper hashes ni el contrato existente.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data" / "mobile_typing"
FULL = "--full" in sys.argv

if FULL:
    TSV_DIR = BASE / "full"
    GFF_DIR = BASE / "gff" / "full"
    SPACER_DIR = BASE / "crispr_spacers" / "full"
    WORK = BASE / "_work" / "full"
    GENOMES = BASE / "full_genomes.tsv"
    SUMMARY = BASE / "summary_full.json"
    VALIDATION = BASE / "validation_full.json"
    README = BASE / "README_full.md"
else:
    TSV_DIR = BASE
    GFF_DIR = BASE / "gff"
    SPACER_DIR = BASE / "crispr_spacers"
    WORK = BASE / "_work"
    GENOMES = BASE / "pilot_genomes.tsv"
    SUMMARY = BASE / "summary.json"
    VALIDATION = BASE / "validation.json"
    README = BASE / "README.md"


def genome_list():
    """Lista de genomas (columna `genome`) del modo activo."""
    import csv

    with GENOMES.open() as fh:
        return [row["genome"] for row in csv.DictReader(fh, delimiter="\t")]
