#!/usr/bin/env bash
# Ejecución reproducible con --run para este gold.
# No descarga herramientas, imágenes ni bases. Bakta DB: exportar BAKTA_DB.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${PROJECT_ROOT:-}" ]]; then
    ROOT="$PROJECT_ROOT"
elif ROOT_GIT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null)"; then
    ROOT="$ROOT_GIT"
else
    ROOT="$(cd "$SCRIPT_DIR/../../../.." 2>/dev/null && pwd || true)"
fi
if [[ -z "${ROOT:-}" || ! -f "$ROOT/benchmarks/external/external_benchmark.py" ]]; then
    ROOT=/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc
fi
cd "$ROOT"
PYTHON="${PYTHON:-$ROOT/.venv/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
    PYTHON="${PYTHON_FALLBACK:-python3}"
fi
exec "$PYTHON" benchmarks/external/external_benchmark.py --genome benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna --reference benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 --tools all --run --threads 1 --timeout 3600.0 --save-gffs reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence --json reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.json --md reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.md "$@"
