#!/usr/bin/env bash
# run_batch_02.sh — lote 02 de evidencia por brecha.
#
# Universo: `data/metabolism/review/gap_review_current.tsv` (4.010 filas).
# Exclusiones: las 500 del lote 01 + accionables vigentes (54) + claves que ya
# no están en `real_gaps` (las 23 reannotated).
#
# Uso:
#   bash scripts/evidence/run_batch_02.sh [workers] [n]
#   (n por defecto 1000; el lote 01 usó 500)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
WORKERS="${1:-7}"
N="${2:-1000}"
S="scripts/evidence"

export EVIDENCE_BATCH=02
export EVIDENCE_REVIEW="$ROOT/data/metabolism/review/gap_review_current.tsv"
export EVIDENCE_N="$N"
export EVIDENCE_EXCLUDE="$ROOT/data/metabolism/evidence/work_batch_02/exclude_keys.tsv"

mkdir -p "$ROOT/data/metabolism/evidence/work_batch_02"
echo "== 00 exclusiones (lote 01)"
python3 "$S/00_make_exclusions.py" \
  --from "$ROOT/data/metabolism/evidence/gap_evidence.tsv" \
  --out "$EVIDENCE_EXCLUDE"
echo "== 01 selección del lote (n=$N)"
python3 "$S/01_select_worklist.py"
echo "== 02 DB de productores (taxonomía + DIAMOND)"
python3 "$S/02_build_producer_db.py"
echo "== 03 dominios Pfam esperados (caché compartida + UniProt)"
python3 "$S/03_resolve_domains.py"
echo "== 04 panel HMM Pfam (caché compartida + InterPro)"
python3 "$S/04_build_domain_panel.py"
echo "== 05 homología all-hits (DIAMOND, NCBI + ab initio)"
python3 "$S/05_run_homology.py" --workers "$WORKERS"
echo "== 06 dominios detectados (hmmscan --cut_ga, paneles mínimos)"
python3 "$S/06_run_domains.py" --workers "$WORKERS"
echo "== 07 tabla de evidencia + summary + README"
python3 "$S/07_build_evidence.py" --batch gap_evidence_batch_02
echo "== 08 validación"
python3 "$S/08_validate.py"
echo "OK: lote 02 en data/metabolism/evidence/"
