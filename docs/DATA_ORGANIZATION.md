# Organización de los datos

Qué hay en cada carpeta y cuál usar según lo que se vaya a medir. La
exclusión de genomas es **deliberada**: se identificaron posibles parásitos o
simbiontes intracelulares, con biología distinta (genomas reducidos, dependencia
del huésped), y se apartaron para no distorsionar las métricas del conjunto
libre.

## Carpetas

| carpeta | contenido | n |
|---|---|---:|
| `data/raw/` | FASTA de los genomas **libres** | 310 |
| `data/raw_excluded/` | FASTA de los **excluidos** | 90 |
| `data/annotations_ab_initio/` | anotación de todos (genes, GFF, proteínas) | 400 |
| `data/regulation/` | capa L2 (genes, operones, promotores) de **todos** | 400 |
| `data/operon_finder/` | horquillas y salidas del motor de los libres + K-12 + *B. subtilis* | 312 horquillas, 310 operones |
| `data/operon_finder/rules_all/` | operones del motor, genomas libres | 310 |
| `data/operon_finder/excluded/` | operones y horquillas de los excluidos | 90 |
| `data/operon_finder/excluded_genomes.json` | **lista canónica** de los 90 | 90 |

## Qué usar según el objetivo

- **Métricas del conjunto libre**: `data/regulation/` **filtrando** con
  `excluded_genomes.json`, o directamente `data/operon_finder/rules_all/`.
- **Anotación completa (catálogo)**: `data/regulation/` sin filtrar, o
  `operon_finder/` + `operon_finder/excluded/`.
- **Genomas con evidencia de terminador**: los 311 (`*_hairpins.tsv` en
  `data/operon_finder/`) más los 90 excluidos en `excluded/`.

## Regla para no repetir el error

**Nunca inferir el conjunto desde la presencia de ficheros.** Un genoma puede
estar anotado (en `data/regulation`) y aun así excluido de las métricas. La
verdad está en `excluded_genomes.json`; el resto es consecuencia.

## Pendiente (fase 3.2)

`data/regulation/` no marca los excluidos. Propuesta para su responsable: añadir
`excluded=true` a `genes.tsv`/`operons.tsv`, o un `data/regulation/excluded.json`
que apunte a la lista canónica. Hasta entonces, filtrar con
`data/operon_finder/excluded_genomes.json`.
