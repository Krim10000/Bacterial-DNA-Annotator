#!/usr/bin/env bash
# annotate_all.sh — anota horquillas en todos los genomas con FASTA y GFF.
#
# Reanudable: salta los genomas que ya tienen `_hairpins.tsv`.
# Paraleliza por genoma (cada proceso abre ViennaRNA y su propio GFF).
#
#   bash scripts/regulation/annotate_all.sh          # 6 procesos
#   JOBS=8 bash scripts/regulation/annotate_all.sh
set -uo pipefail
cd "$(dirname "$0")/../.."

PY=.venv/bin/python
OUT=data/operon_finder
JOBS=${JOBS:-6}
mkdir -p "$OUT"

pending=()
for gff in data/annotations_ab_initio/*/annotation.gff3; do
  g=$(basename "$(dirname "$gff")")
  fa="data/raw/${g}_genomic.fna.gz"
  [ -f "$fa" ] || continue
  [ -s "$OUT/${g}_hairpins.tsv" ] && continue
  pending+=("$g")
done

if [ "${ORDER:-forward}" = "reverse" ]; then
  mapfile -t pending < <(printf '%s\n' "${pending[@]}" | tac)
  echo "orden inverso (para un segundo lote sin solaparse)"
fi
echo "pendientes: ${#pending[@]} genomas | ${JOBS} en paralelo"
[ ${#pending[@]} -eq 0 ] && exit 0

run_one() {
  local g=$1
  nice -n 10 "$PY" scripts/regulation/annotate_hairpins.py \
    --gff "data/annotations_ab_initio/$g/annotation.gff3" \
    --fasta "data/raw/${g}_genomic.fna.gz" \
    --genome "$g" > "logs/hairpins_$g.log" 2>&1
  if [ -s "$OUT/${g}_hairpins.tsv" ]; then
    echo "ok $g"
  else
    echo "FALLO $g (ver logs/hairpins_$g.log)"
  fi
}
export -f run_one
export PY OUT

printf '%s\n' "${pending[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}

echo "terminado: $(ls "$OUT"/*_hairpins.tsv 2>/dev/null | wc -l) genomas con horquillas"
