#!/usr/bin/env bash
# Construye un panel ncRNA pequeño desde un Rfam.cm local.
# Uso: bash scripts/rrna/00_build_ncrna_panel.sh [Rfam.cm]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RFAM="${1:-$ROOT/Rfam.cm}"
OUT="$ROOT/data/rfam_models/bac_ncRNA.cm"
[ -f "$RFAM" ] || { echo "Falta $RFAM"; exit 1; }
command -v cmfetch >/dev/null || { echo "Falta cmfetch (Infernal)"; exit 1; }
mkdir -p "$(dirname "$OUT")"
: > "$OUT"
for name in RNaseP_bact_a RNaseP_bact_b Bacteria_small_SRP Bacteria_large_SRP tmRNA 6S CsrB CsrC; do
  cmfetch "$RFAM" "$name" >> "$OUT" 2>/dev/null || true
done
grep -c '^NAME' "$OUT"
echo "escrito: $OUT"
