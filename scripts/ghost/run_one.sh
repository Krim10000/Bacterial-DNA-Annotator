#!/usr/bin/env bash
# Corre el descubrimiento de genes fantasma para UN genoma.
# Uso: bash scripts/ghost/run_one.sh GCF_000006985.1_ASM698v1
set -euo pipefail

GENOME=${1:?uso: bash scripts/ghost/run_one.sh <GCF_...>}
PROJECT_ROOT=$(cd "$(dirname "$0")/../.." && pwd)
SCRIPT_DIR="$PROJECT_ROOT/scripts/ghost"
DATA="$PROJECT_ROOT/data"
PY="$PROJECT_ROOT/.venv/bin/python"
[ -x "$PY" ] || { echo "ERROR: no existe $PY" >&2; exit 1; }

export TORCH_HOME="$PROJECT_ROOT/.cache/torch"
export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$TORCH_HOME" "$PROJECT_ROOT/out"

cd "$PROJECT_ROOT"
exec "$PY" -m ghost.run \
  --genome "$GENOME" \
  --data-root "$DATA" \
  --out "$PROJECT_ROOT/out" \
  --esm --esm-top 100 \
  --diamond proteomes \
  --db "$DATA/db/proteomes.dmnd" \
  --threads "${THREADS:-4}"
