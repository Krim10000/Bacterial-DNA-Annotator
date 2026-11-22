#!/usr/bin/env bash
# run_closure_domain.sh — cadena completa del cierre de brechas por dominio.
#
#   worklist  ->  alcance, productores, dominios permitidos (core/majority/any)
#   scan      ->  hmmscan por genoma con panel mínimo (paralelo)
#   merge     ->  veredicto por brecha frente al catálogo previo
#   homology  ->  DIAMOND contra las enzimas productoras (método independiente)
#   final     ->  cruce dominio x homología -> validated / consensus / candidate
#
# Uso:
#   bash scripts/ghost/run_closure_domain.sh            # lote 1000, sin red
#   bash scripts/ghost/run_closure_domain.sh 50         # lote pequeño
#   RESOLVE=1 bash scripts/ghost/run_closure_domain.sh  # completa UniProt
#   FORCE=1   bash scripts/ghost/run_closure_domain.sh  # re-escanea todo
#
# Requisitos: data/db/pfam/panel_domain_closure.hmm (hmmfetch sobre Pfam-A) y
# hmmscan/hmmfetch/hmmpress en PATH. Sin la red el paso 1 usa la caché local.
set -euo pipefail
cd "$(dirname "$0")/../.."

PY=.venv/bin/python
TOP=${1:-1000}
WORKERS=${WORKERS:-8}

echo "== 1/3 worklist (top $TOP) =="
if [ "${RESOLVE:-0}" = "1" ]; then
  "$PY" scripts/ghost/build_domain_worklist.py --top "$TOP" --resolve
else
  "$PY" scripts/ghost/build_domain_worklist.py --top "$TOP"
fi

echo "== 2/3 scan ($WORKERS workers) =="
if [ "${FORCE:-0}" = "1" ]; then
  "$PY" scripts/ghost/scan_domain_closure.py --top "$TOP" \
      --workers "$WORKERS" --threads 1 --force
else
  "$PY" scripts/ghost/scan_domain_closure.py --top "$TOP" \
      --workers "$WORKERS" --threads 1
fi

echo "== 3/5 veredictos vs catálogo previo =="
"$PY" scripts/ghost/merge_closure_validation.py --top "$TOP"

echo "== 4/5 homología independiente (sin Pfam) =="
"$PY" scripts/ghost/homology_closure.py --top "$TOP" --build-db
"$PY" scripts/ghost/homology_closure.py --top "$TOP" --workers "$WORKERS"

echo "== 5/5 cruce dominio x homología =="
"$PY" scripts/ghost/homology_closure.py --top "$TOP" --merge

echo
echo "salidas:"
echo "  data/metabolism/domain_worklist_top$TOP.tsv"
echo "  data/metabolism/domain_closure_candidates_top$TOP.tsv"
echo "  data/metabolism/closure_validation_top$TOP.tsv"
echo "  data/metabolism/homology_closure_top$TOP.tsv"
echo "  data/metabolism/closure_final_top$TOP.tsv"
