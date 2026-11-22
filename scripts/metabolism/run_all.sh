#!/usr/bin/env bash
# Pipeline reproducible de la capa T1 bacteriana de transporte.
#
#   bash scripts/metabolism/run_all.sh [--threads N] [--force]
#
# Solo escribe en data/metabolism/transport_t1/. No hace commits.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

THREADS="${TRANSPORT_THREADS:-4}"
FORCE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE="--force"; shift ;;
    *) echo "argumento desconocido: $1" >&2; exit 2 ;;
  esac
done

PY="${PYTHON:-python3}"

echo "== transport_t1 $(date -Iseconds) =="
"$PY" "$HERE/01_build_reference.py"
"$PY" "$HERE/02_select_pilot.py"
"$PY" "$HERE/03_run_diamond.py" --threads "$THREADS" $FORCE
"$PY" "$HERE/04_classify.py"
"$PY" "$HERE/05_validate.py"

echo "== listo: data/metabolism/transport_t1/ =="
