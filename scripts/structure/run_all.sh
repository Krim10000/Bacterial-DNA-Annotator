#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python3 scripts/structure/00_select_structural.py
python3 scripts/structure/01_fetch_afdb.py
echo "Siguiente: autentica Modal (modal token new) y ejecuta:"
echo "modal run scripts/structure/modal_esmfold.py --fasta data/metabolism/targets/structural_candidates.faa --out data/metabolism/targets/structures/esmfold --chunk 5"
echo "Después: bash scripts/structure/02_foldseek.sh"
