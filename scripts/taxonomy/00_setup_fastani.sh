#!/usr/bin/env bash
# Descarga el binario estatico de fastANI (release pineada) al cache de trabajo.
# No instala nada en el sistema y no toca el repo.
set -euo pipefail

WORK="${TAXONOMY_WORK:-$HOME/.cache/dna_poc/taxonomy}"
VERSION="v1.34"          # release de GitHub; el binario reporta 1.33
URL="https://github.com/ParBLiSS/FastANI/releases/download/${VERSION}/fastANI-linux64-${VERSION}.zip"

mkdir -p "$WORK/bin"
if [ -x "$WORK/bin/fastANI" ]; then
  echo "fastANI ya presente: $WORK/bin/fastANI"
else
  echo "descargando $URL"
  tmp="$(mktemp -d)"
  curl -sSL --fail --retry 3 -o "$tmp/fastANI.zip" "$URL"
  unzip -q -o "$tmp/fastANI.zip" -d "$tmp"
  install -m 0755 "$tmp/fastANI" "$WORK/bin/fastANI"
  rm -rf "$tmp"
fi
"$WORK/bin/fastANI" --version
