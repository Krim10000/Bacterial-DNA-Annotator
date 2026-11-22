#!/usr/bin/env bash
# run_batch_03.sh — lote 03 de evidencia por brecha.
#
# Universo: `real_gaps.tsv` (3.667 enzymatic_missing vigentes), no la revisión
# histórica de 4.021. Exclusiones (work_batch_03/exclude_keys.tsv):
#   * 1.500 ya cubiertas por los lotes 01+02 (gap_evidence_all.tsv);
#   * 77 accionables históricas;
#   * 343 promovidas a homology_present (ya dentro de las 1.500);
#   * 23 reannotated (ya dentro de las 77 accionables).
# El pool real es 2.456; el lote toma 1.000.
#
# Uso:
#   bash scripts/evidence/run_batch_03.sh [workers] [n]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
WORKERS="${1:-7}"
N="${2:-1000}"
S="scripts/evidence"

export EVIDENCE_BATCH=03
export EVIDENCE_REVIEW="$ROOT/data/metabolism/evidence/work_batch_03/universe.tsv"
export EVIDENCE_N="$N"
export EVIDENCE_EXCLUDE="$ROOT/data/metabolism/evidence/work_batch_03/exclude_keys.tsv"
# snapshot congelado: base y selección no dependen de cambios posteriores
export EVIDENCE_REAL="$ROOT/data/metabolism/evidence/work_batch_03/real_snapshot.tsv"
export EVIDENCE_QUEUE="$ROOT/data/metabolism/evidence/work_batch_03/queue_snapshot.tsv"
export EVIDENCE_REUSE_SELECTION=1

mkdir -p "$ROOT/data/metabolism/evidence/work_batch_03"
echo "== 09 universo y exclusiones (3.667 → pool 2.456)"
python3 "$S/09_prepare_batch03.py"
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
echo "== 07 tabla + no_bacterial_reference + summary + README"
python3 "$S/07_build_evidence.py" --batch gap_evidence_batch_03
echo "== 10 auditoría de deriva post-cierre"
python3 "$S/10_audit_batch03_drift.py"
echo "== 08 validación"
python3 "$S/08_validate.py"
echo "OK: lote 03 en data/metabolism/evidence/"
