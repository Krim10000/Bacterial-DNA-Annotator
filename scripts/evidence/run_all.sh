#!/usr/bin/env bash
# run_all.sh — reconstruye la evidencia por brecha de punta a punta.
#
# Uso:
#   bash scripts/evidence/run_all.sh [workers]
#
# Reanudable: 03 usa caché de Pfam, 04 no re-descarga HMMs presentes y 05/06
# reutilizan los crudos de work/diamond y work/domtbl.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
WORKERS="${1:-7}"
S="scripts/evidence"

echo "== 01 selección del lote"
python3 "$S/01_select_worklist.py"
echo "== 02 DB de productores (taxonomía + DIAMOND)"
python3 "$S/02_build_producer_db.py"
echo "== 03 dominios Pfam esperados (UniProt + política core/majority/any)"
python3 "$S/03_resolve_domains.py"
echo "== 04 panel HMM Pfam (InterPro) + hmmpress"
python3 "$S/04_build_domain_panel.py"
echo "== 05 homología all-hits (DIAMOND, NCBI + ab initio)"
python3 "$S/05_run_homology.py" --workers "$WORKERS"
echo "== 06 dominios detectados (hmmscan --cut_ga)"
python3 "$S/06_run_domains.py" --workers "$WORKERS"
echo "== 07 tabla de evidencia + summary + README"
python3 "$S/07_build_evidence.py"
echo "== 08 validación"
python3 "$S/08_validate.py"
echo "OK: evidencia en data/metabolism/evidence/"
