#!/usr/bin/env python3
"""00_select_genomes.py — fija los 10 genomas piloto de la capa móvil.

Criterio (documentado y determinista, no aleatorio):
  1. Deben existir las tres entradas necesarias: GFF ab initio, GFF NCBI y FASTA
     genómico (para escanear arrays CRISPR de novo).
  2. Diez filos distintos, para cubrir arquitecturas de genoma muy diferentes
     (GC%, tamaño, número de contigs, carga de IS/profagos).
  3. Dentro de cada filo se elige un genoma con carga móvil/defensa previa
     suficiente para que el piloto sea informativo (≥10 elementos previos en
     `data/mobile_elements.json`) y, cuando es posible, multi-contig, porque el
     mapeo de coordenadas sin contig es el problema central del piloto.
  4. Se prefieren tamaños moderados (2–9 Mb) y especies de interés clínico o
     ambiental, con al menos un genoma por filo y ≥5 multi-contig en total.

Salida: `data/mobile_typing/pilot_genomes.tsv`.

Uso:
    python scripts/mobile/00_select_genomes.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANN = ROOT / "data" / "annotations_ab_initio"
NCBI_GFF = ROOT / "data" / "annotations"
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "mobile_typing"

# (GCF, filo, motivo de selección)
PILOT = [
    ("GCF_019048105.1_ASM1904810v1", "Pseudomonadota", "Providencia rettgeri: enterobacteria clínica, 3 contigs, carga AMR/móvil alta"),
    ("GCF_025558325.1_ASM2555832v1", "Bacillota", "Staphylococcus schleiferi: patógeno de piel, genoma pequeño con IS y cas"),
    ("GCF_001611795.1_ASM161179v1", "Actinomycetota", "Streptomyces qaidamensis: alto GC, 9 Mb, islas y elementos móviles"),
    ("GCF_018127985.1_ASM1812798v1", "Bacteroidota", "Prevotella multiformis: microbiota oral, 2 contigs, CRISPR previo"),
    ("GCF_000317105.1_ASM31710v1", "Cyanobacteriota", "Oxynema acuminatum: cianobacteria filamentosa, 3 contigs, CRISPR/Cas"),
    ("GCF_045784305.1_DRML_24019", "Deinococcota", "Deinococcus radiomollis: radiación, 9 contigs, muchos IS"),
    ("GCF_000147715.2_ASM14771v3", "Thermotogota", "Mesotoga prima: termófilo, 2 contigs, CRISPR previo"),
    ("GCF_017893965.1_ASM1789396v1", "Spirochaetota", "Treponema parvum: espiroqueta, defensa previa alta"),
    ("GCF_000092785.1_ASM9278v1", "Chlamydiota", "Waddlia chondrophila: intracellular, 2 contigs, móvil alto"),
    ("GCF_002305895.1_ASM230589v1", "Myxococcota", "Corallococcus macrosporus: 9 Mb, profagos y defensa"),
]


def read_sequence_regions(path: Path, gzipped: bool):
    """Contigs (accession, length) en el orden del fichero."""
    import gzip

    opener = gzip.open if gzipped else open
    out = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                parts = line.split()
                out.append((parts[1], int(parts[3])))
            elif not line.startswith("#"):
                break
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    problems = []
    for gcf, phylum, why in PILOT:
        gff_ab = ANN / gcf / "annotation.gff3"
        gff_nc = NCBI_GFF / f"{gcf}.gff.gz"
        fna = RAW / f"{gcf}_genomic.fna.gz"
        missing = [str(p) for p in (gff_ab, gff_nc, fna) if not p.exists()]
        if missing:
            problems.append(f"{gcf}: faltan {missing}")
            continue
        contigs = read_sequence_regions(gff_ab, gzipped=False)
        rows.append(
            {
                "genome": gcf,
                "phylum": phylum,
                "contigs": len(contigs),
                "genome_bp": sum(length for _, length in contigs),
                "contig_accessions": ",".join(acc for acc, _ in contigs),
                "reason": why,
            }
        )

    if problems:
        for p in problems:
            print(f"ERROR: {p}", file=sys.stderr)
        return 1
    if len({r["phylum"] for r in rows}) != len(rows):
        print("ERROR: filos repetidos en la selección", file=sys.stderr)
        return 1

    out_path = OUT / "pilot_genomes.tsv"
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["genome", "phylum", "contigs", "genome_bp", "contig_accessions", "reason"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"OK: {len(rows)} genomas piloto -> {out_path.relative_to(ROOT)}")
    for r in rows:
        print(f"  {r['genome']:36s} {r['phylum']:16s} contigs={r['contigs']:2d} {r['genome_bp']/1e6:6.2f} Mb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
