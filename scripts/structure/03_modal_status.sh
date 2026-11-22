#!/usr/bin/env bash
set -euo pipefail
echo "== Modal apps =="
modal app list || true
echo "== Modal containers =="
modal container list || true
echo "Recuerda: máximo 10 GPU activas por cuenta."
