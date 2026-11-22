#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python3 scripts/targets/00_prepare.py
python3 scripts/targets/01_extract_hmms.py
python3 scripts/targets/02_run_search.py
python3 scripts/targets/03_rank.py
