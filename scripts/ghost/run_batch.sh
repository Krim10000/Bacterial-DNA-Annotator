#!/usr/bin/env bash
# Corre el descubrimiento sobre todos los genomas del proyecto.
# Uso: bash scripts/ghost/run_batch.sh [N]
#   N = límite opcional de genomas (para probar)
set -uo pipefail

LIMIT=${1:-0}
PROJECT_ROOT=$(cd "$(dirname "$0")/../.." && pwd)
SCRIPT_DIR="$PROJECT_ROOT/scripts/ghost"
DATA="$PROJECT_ROOT/data"
PY="$PROJECT_ROOT/.venv/bin/python"
OUT="$PROJECT_ROOT/out"
[ -x "$PY" ] || { echo "ERROR: no existe $PY" >&2; exit 1; }

export TORCH_HOME="$PROJECT_ROOT/.cache/torch"
export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$TORCH_HOME" "$OUT"

ok=0; fail=0; skip=0; n=0
for d in "$DATA"/annotations_ab_initio/*/; do
  g=$(basename "$d")
  [ -f "$d/annotation.gff3" ] || continue
  if [ -f "$OUT/${g}_ghost_summary.json" ]; then skip=$((skip+1)); continue; fi
  if [ "$LIMIT" -gt 0 ] && [ "$n" -ge "$LIMIT" ]; then break; fi
  echo "=== $g ==="
  if "$PY" -m ghost.run --genome "$g" --data-root "$DATA" --out "$OUT" \
       --esm --esm-top 100 --diamond proteomes \
       --db "$DATA/db/proteomes.dmnd" --threads "${THREADS:-4}" \
       >> "$OUT/batch.log" 2>&1; then
    ok=$((ok+1))
  else
    fail=$((fail+1)); echo "FALLO $g" >> "$OUT/batch.log"
  fi
  n=$((n+1))
done
echo "Listo: ok=$ok fallos=$fail omitidos=$skip"
