#!/usr/bin/env bash
# 22_run_full.sh — cadena completa de los 400 genomas (reanudable por etapas).
#
#   bash scripts/mobile/22_run_full.sh prepare     # 19 inventario
#   bash scripts/mobile/22_run_full.sh extract     # 01 en 3 shards
#   bash scripts/mobile/22_run_full.sh pfam        # 20 (panel núcleo)
#   bash scripts/mobile/22_run_full.sh crispr      # 21
#   bash scripts/mobile/22_run_full.sh build       # 05 + 06 + 07
#   bash scripts/mobile/22_run_full.sh all
#
# Cada etapa salta lo ya hecho: 20/21 por ficheros de salida, 05/06/07
# reconstruyen desde los intermedios existentes.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python3}"
STAGE="${1:-all}"
PFAM_WORKERS="${PFAM_WORKERS:-5}"
CRISPR_WORKERS="${CRISPR_WORKERS:-5}"

case "$STAGE" in
  prepare) "$PY" scripts/mobile/19_prepare_full.py ;;
  extract)
    for i in 1 2 3; do
      "$PY" scripts/mobile/01_extract_candidates.py --full --shard "$i" --shards 3 &
    done
    wait ;;
  pfam)  "$PY" scripts/mobile/20_pfam_pool.py --workers "$PFAM_WORKERS" ;;
  crispr) "$PY" scripts/mobile/21_crispr_pool.py --workers "$CRISPR_WORKERS" ;;
  build)
    "$PY" scripts/mobile/05_classify.py --full
    "$PY" scripts/mobile/06_write_outputs.py --full
    "$PY" scripts/mobile/07_validate.py --full ;;
  all)
    "$0" prepare
    "$0" extract
    "$0" pfam & PFAM_PID=$!
    "$0" crispr & CRISPR_PID=$!
    wait "$PFAM_PID" "$CRISPR_PID"
    "$0" build ;;
  *)
    echo "Uso: $0 [prepare|extract|pfam|crispr|build|all]" >&2
    exit 1 ;;
esac
