#!/usr/bin/env bash
# run_all.sh — cadena completa de la capa móvil para los 10 genomas piloto.
#
# Uso:
#   bash scripts/mobile/run_all.sh            # todo (incluye hmmscan del panel)
#   SKIP_SCAN=1 bash scripts/mobile/run_all.sh  # reutiliza pfam_hits.tsv existentes
#
# Requisitos: python3, hmmscan/hmmpress (HMMER 3.3) y `.cache/pfam/Pfam-A.hmm`.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python3}"

echo "== 00 selección de genomas piloto =="
"$PY" scripts/mobile/00_select_genomes.py
echo "== 01 extracción y mapeo de candidatos previos =="
"$PY" scripts/mobile/01_extract_candidates.py
echo "== 02 panel Pfam curado =="
"$PY" scripts/mobile/02_build_pfam_panel.py
if [ "${SKIP_SCAN:-0}" != "1" ]; then
  echo "== 03 hmmscan del panel =="
  "$PY" scripts/mobile/03_scan_pfam.py
fi
echo "== 04 arrays CRISPR de novo =="
"$PY" scripts/mobile/04_crispr_arrays.py
echo "== 05 tipado =="
"$PY" scripts/mobile/05_classify.py
echo "== 06 salidas TSV/GFF/summary =="
"$PY" scripts/mobile/06_write_outputs.py
echo "== 07 validación =="
"$PY" scripts/mobile/07_validate.py
