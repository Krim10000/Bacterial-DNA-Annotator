#!/usr/bin/env bash
# run_rules_all.sh — motor determinista de operones sobre los genomas del proyecto.
#
# Cada genoma usa SUS horquillas si existen; si no, corre sin evidencia de
# terminador (solo distancia). Reanudable: salta los que ya tienen salida.
#
#   JOBS=4 bash scripts/regulation/run_rules_all.sh
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
PY=".venv/bin/python"
[ -x "$PY" ] || PY="python"
JOBS="${JOBS:-4}"
OUT="data/operon_finder/rules_all"
mkdir -p "$OUT"

run_one() {
  g="$1"
  gbk="data/annotations_ab_initio/$g/annotation.gbk"
  [ -f "$gbk" ] || return 0
  dest="$OUT/${g}_operons_rules.tsv"
  [ -s "$dest" ] && return 0
  hp="data/operon_finder/${g}_hairpins.tsv"
  args=(--gbff "$gbk" --genome "$g" --out "$OUT")
  [ -f "$hp" ] && args+=(--hairpins "$hp")
  if "$PY" scripts/regulation/rule_operons.py "${args[@]}" > "logs/rules_${g}.log" 2>&1; then
    echo "ok $g"
  else
    echo "FALLO $g (logs/rules_${g}.log)"
  fi
}
export -f run_one
export PY OUT

# los genomas excluidos (parásitos/simbiontes) van a operon_finder/excluded/ y
# no entran en el conjunto libre: se saltan aquí
EXCLUIDOS="$(ls data/raw_excluded/*_genomic.fna.gz 2>/dev/null | xargs -r -n1 basename | sed 's/_genomic.fna.gz//')"
ls -d data/annotations_ab_initio/GCF_* 2>/dev/null | xargs -n1 basename | \
  grep -vxF -f <(printf '%s\n' $EXCLUIDOS) | \
  xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
TOTAL=$(ls -d data/annotations_ab_initio/GCF_* 2>/dev/null | xargs -n1 basename | \
        grep -vxF -f <(printf '%s\n' $EXCLUIDOS) | wc -l)
echo "listos: $(ls "$OUT"/*_operons_rules.tsv 2>/dev/null | wc -l) de ${TOTAL} libres"
