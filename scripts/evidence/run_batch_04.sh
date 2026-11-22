#!/usr/bin/env bash
# run_batch_04.sh — lote 04 (cierre): todos los unknown restantes.
#
# Universo: `real_gaps.tsv` (3.433 enzymatic_missing). Exclusiones:
#   * 2.500 ya cubiertas por 01-03 (gap_evidence_all.tsv);
#   * 77 accionables históricas (54 siguen dentro);
#   * 23 reannotated;
#   * promovidas a homology_present (incluye 343 + 213 del contexto).
# Pool real: 1.456. Se procesan todas; las brechas sin productor bacteriano
# curado se separan en `no_bacterial_reference_batch_04.tsv` (no ausencia).
#
# Uso:
#   bash scripts/evidence/run_batch_04.sh [workers] [n]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
WORKERS="${1:-7}"
N="${2:-2000}"
S="scripts/evidence"
W="$ROOT/data/metabolism/evidence/work_batch_04"

export EVIDENCE_BATCH=04
export EVIDENCE_REVIEW="$W/universe.tsv"
export EVIDENCE_N="$N"
export EVIDENCE_EXCLUDE="$W/exclude_keys.tsv"
export EVIDENCE_REAL="$W/real_snapshot.tsv"
export EVIDENCE_QUEUE="$W/queue_snapshot.tsv"
export EVIDENCE_REUSE_SELECTION=1

mkdir -p "$W"
echo "== 11 universo y exclusiones (3.433 → pool 1.456)"
python3 "$S/11_prepare_batch04.py"
echo "== 01 selección del lote (n=$N; procesa todo el pool)"
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
python3 "$S/07_build_evidence.py" --batch gap_evidence_batch_04
echo "== 10 auditoría de deriva post-cierre"
python3 "$S/10_audit_batch03_drift.py"
echo "== 08 validación"
python3 "$S/08_validate.py"
echo "OK: lote 04 en data/metabolism/evidence/"
