# scripts/mobile — capa móvil / defensa / CRISPR / AMR / virulencia

Cadena reproducible del piloto de 10 genomas. Todo escribe en
`data/mobile_typing/` (nunca en `data/annotations_ab_initio/`, `webapp/`,
`data/metabolism/`, `data/regulation/`).

## Orden

| script | función | salida principal |
|---|---|---|
| `00_select_genomes.py` | fija los 10 genomas piloto (10 filos, 6 multi-contig) | `data/mobile_typing/pilot_genomes.tsv` |
| `01_extract_candidates.py` | une JSON previos → CDS NCBI → CDS ab initio por solapamiento recíproco ≥0,8 | `_work/<genoma>/candidates.tsv`, `cds.tsv`, `contigs.tsv` |
| `02_build_pfam_panel.py` | extrae 339 modelos Pfam curados en una pasada del Pfam-A completo | `_work/pfam_panel/panel.hmm` |
| `03_scan_pfam.py` | `hmmscan --domE 1e-3` sobre `proteins.faa` | `_work/<genoma>/pfam_hits.tsv` |
| `04_crispr_arrays.py` | detector CRISPR tipo CRT sobre el FASTA | `_work/<genoma>/crispr_arrays.tsv` + espaciadores |
| `05_classify.py` | tipado por nombre + Pfam + previo; regiones de profago, RM y TA | `_work/<genoma>/elements.tsv`, `demoted.tsv` |
| `06_write_outputs.py` | TSV de contrato, GFF3, espaciadores, `summary.json` | `data/mobile_typing/` |
| `07_validate.py` | valida TSV/GFF/espaciadores y hashes | `validation.json` |

```bash
bash scripts/mobile/run_all.sh
SKIP_SCAN=1 bash scripts/mobile/run_all.sh   # reutiliza pfam_hits.tsv (evita ~25 min)
```

## Requisitos

- HMMER 3.3 (`hmmscan`, `hmmpress`) en `PATH`.
- `.cache/pfam/Pfam-A.hmm` (3,6 GB) y `.cache/Pfam-A.hmm.dat.gz`.
- `data/annotations_ab_initio/<genoma>/annotation.gff3`, `proteins.faa`,
  `functional_annotation.json`; `data/annotations/<genoma>.gff.gz`;
  `data/raw/<genoma>_genomic.fna.gz`; los tres JSON previos de `data/`.

## Panel Pfam

`pfam_panel_spec.json` lista las accesiones por categoría (transposasa,
integrasa, cas, restricción, toxina-antitoxina, fago, conjugación, defensa
otros, beta-lactamasa, aminoglucósido, cloranfenicol, efflux, tetraciclina,
virulencia). Se seleccionó por `NAME` de Pfam, no por descripción, para evitar
falsos amigos. Revisar a mano si se amplía.

## Notas

- Sin `--cut_ga` a propósito: el umbral GA se guarda por hit (`passes_ga`) y se
  usa para graduar la confianza; un hit débil solo apoya.
- El TSV final tiene exactamente las 10 columnas del contrato de `BACKLOG.md`.
- Las coordenadas son del GFF ab initio (1-based inclusivo); ver
  `data/mobile_typing/README.md` para el método, límites y hallazgos.
