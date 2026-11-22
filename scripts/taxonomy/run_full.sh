#!/usr/bin/env bash
# Capa taxonómica completa (400 genomas). Reutiliza el cache del piloto.
#
#   bash scripts/taxonomy/run_full.sh
#
# Es resumible: 06 salta los grupos ya completados o con hits escritos, así que
# relanzar el script tras un corte continúa donde quedó.
#
# Variables:
#   TAXONOMY_WORK        cache (FASTA, panel, GTDB)   [~/.cache/dna_poc/taxonomy]
#   TAXONOMY_THREADS     hilos para fastANI           [8]
#   TAXONOMY_MAX_REFS    refs por filo en el barrido  [5]
set -euo pipefail

cd "$(dirname "$0")/../.."
HERE="scripts/taxonomy"

echo "== 05 metadata + panel + plan de alcance =="
python3 "$HERE/05_build_panel_full.py"

echo "== 06 ANI: barrido por filo (resumible) =="
python3 "$HERE/06_run_fastani_full.py"

echo "== 06b ANI: pasada dirigida (conspecíficas / familia / pilotos) =="
python3 "$HERE/06b_targeted_ani.py"

echo "== 07 asignación de especie/linaje (400 JSON) =="
python3 "$HERE/07_assign_taxonomy_full.py"

echo "== 08 validación =="
python3 "$HERE/08_validate_full.py"

echo "== 09 consistencia con el piloto =="
python3 "$HERE/09_compare_layers.py"

echo "== 10 estadísticas =="
python3 "$HERE/10_stats_full.py"

echo "OK -> data/taxonomy/full/"
