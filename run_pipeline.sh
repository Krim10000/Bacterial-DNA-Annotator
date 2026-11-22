#!/usr/bin/env bash
# run_pipeline.sh — orquestador por etapas.
#
# Etapas de datos base y modelo, más el cierre de brechas por dominio:
#   bash run_pipeline.sh closure        # worklist -> scan -> veredictos
#   bash run_pipeline.sh closure 100    # lote de 100 brechas
#
# Nota: tras la reorganización de `scripts/` en subcarpetas, este fichero seguía
# apuntando a las rutas viejas y a `python` del sistema. Ambas cosas corregidas.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p logs reports data/tokenizer data/processed

PY="python"
[ -x ".venv/bin/python" ] && PY=".venv/bin/python"

# Sin argumentos NO se ejecuta todo: se pide la etapa explícita. Evita lanzar
# descargas completas por accidente (pasó dos veces durante el desarrollo).
STAGE="${1:-}"
ARG="${2:-}"
if [ -z "$STAGE" ]; then
  echo "Uso: $0 <etapa>  (usa 'all' para la cadena completa)"
  echo "  datos:  catalog|filter|phylum|sample|download|corpus|step1"
  echo "  modelo: tokenizer|precompute|verify|step2"
  echo "  opera:  tss|rules|rules_all|coverage|hairpins|suspects|closure [top]"
echo "  red:    census (consulta GEO, requiere internet)"
  echo "  todo:   all"
  exit 1
fi

case "$STAGE" in
  catalog)    "$PY" scripts/data_prep/00_download_catalog.py ;;
  filter)     "$PY" scripts/data_prep/01_filter_candidates.py ;;
  phylum)     "$PY" scripts/data_prep/02_annotate_phylum.py ;;
  sample)     "$PY" scripts/data_prep/03_stratified_sample.py ;;
  download)   "$PY" scripts/data_prep/04_download_fasta.py ;;
  corpus)     "$PY" scripts/data_prep/05_build_corpus.py ;;
  tokenizer)  "$PY" scripts/07_ml/06_train_tokenizer.py ;;
  precompute) "$PY" scripts/07_ml/07_precompute.py ;;
  verify)     "$PY" scripts/07_ml/08_verify_tokenizer.py ;;
  rna)        "$PY" scripts/annotation/18b_rebuild_rna_strand.py ;;
  parse)      "$PY" scripts/annotation/14_parse_gff.py ;;
  closure)    bash scripts/ghost/run_closure_domain.sh "${ARG:-1000}" ;;
  doctor)     "$PY" scripts/00_doctor.py ;;
  tss)
    # unir TSS curados y de alto rendimiento: la regla del promotor gana a todo
    # el resto cuando el dato existe (docs/TSS_FIRST.md)
    "$PY" scripts/validation/build_tss_union.py \
        --curated data/regulondb/tss_k12.tsv \
        --ht data/regulondb/ht/ht_tss.tsv \
        --out data/regulondb/tss_k12_union.tsv ;;
  rules)
    TSSARG=""
    [ -f data/regulondb/tss_k12_union.tsv ] && TSSARG="--tss data/regulondb/tss_k12_union.tsv --gap 99999"
    "$PY" scripts/regulation/rule_operons.py --prom 14 $TSSARG ;;
  suspects)   "$PY" scripts/regulation/operon_suspects.py ;;
  rules_all)  JOBS="${JOBS:-4}" bash scripts/regulation/run_rules_all.sh ;;
  coverage)   "$PY" scripts/validation/coverage_map.py ;;
  census)     "$PY" scripts/validation/census_geo.py ;;
  hairpins)   "$PY" scripts/regulation/annotate_hairpins.py --k12 ;;
  components)
    # verdad de referencia + inventario de piezas por operón (RegulonDB)
    "$PY" scripts/validation/regulondb_gt.py --offline
    "$PY" scripts/validation/operon_components.py ;;
  pwms)
    # PWMs de TF y de cajas -10/-35 por sigma desde sitios curados
    "$PY" scripts/validation/build_tf_pwms.py --test --promoters ;;
  step1)
    for s in catalog filter phylum sample download; do
      echo "=== $s ==="; "$PY" scripts/data_prep/0*_${s}*.py
    done ;;
  step2)
    echo "=== corpus ===";     "$PY" scripts/data_prep/05_build_corpus.py
    echo "=== tokenizer ===";  "$PY" scripts/07_ml/06_train_tokenizer.py
    echo "=== precompute ==="; "$PY" scripts/07_ml/07_precompute.py
    echo "=== verify ===";     "$PY" scripts/07_ml/08_verify_tokenizer.py ;;
  all)
    echo "=== step1 ==="; bash "$0" step1
    echo "=== step2 ==="; bash "$0" step2 ;;
  *)
    echo "Uso: $0 [catalog|filter|phylum|sample|download|corpus|tokenizer|"
    echo "          precompute|verify|rna|parse|closure [top]|components|pwms|"
    echo "          step1|step2|all]"
    exit 1 ;;
esac
