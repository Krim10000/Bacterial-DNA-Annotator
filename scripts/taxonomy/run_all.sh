#!/usr/bin/env bash
# Pipeline completo de taxonomia (ANI + linaje NCBI/GTDB) para los pilotos.
#
#   bash scripts/taxonomy/run_all.sh
#
# Variables utiles:
#   TAXONOMY_WORK        cache de trabajo (FASTA, referencias, GTDB)  [~/.cache/dna_poc/taxonomy]
#   FASTANI              ruta al binario fastANI                      [cache/bin/fastANI o PATH]
#   MAX_CONSPECIFIC      conspecificos por especie                    [4]
#   MAX_PER_GENUS        genomas de tipo por genero                   [200]
#   TAXONOMY_THREADS     hilos para fastANI                           [8]
set -euo pipefail

cd "$(dirname "$0")/../.."
HERE="scripts/taxonomy"

echo "== 00 herramientas =="
bash "$HERE/00_setup_fastani.sh"
python3 "$HERE/00_check_tools.py" > /dev/null

echo "== 01 panel de referencia =="
python3 "$HERE/01_build_refpanel.py"

echo "== 02 ANI (fastANI) =="
python3 "$HERE/02_run_fastani.py"

echo "== 03 asignacion de especie/linaje =="
python3 "$HERE/03_assign_taxonomy.py"

echo "== 04 validacion =="
python3 "$HERE/04_validate.py"

echo "OK -> data/taxonomy/"
