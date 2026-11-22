#!/usr/bin/env bash
# Construye el grafo Rhea, busca cadenas y valida las salidas.
# Uso: bash scripts/regulation/run.sh [--include-unmapped-rhea] [--strict-cofactors]
set -euo pipefail
cd "$(dirname "$0")/../.."

python3 scripts/regulation/build_chains.py "$@"
python3 scripts/regulation/validate_chains.py
