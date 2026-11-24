# Cobertura móvil, defensa y AMR — estado real en disco

Auditoría de **solo lectura** generada por `tools/mobile_coverage.py`. No ejecuta HMMER, DefenseFinder, CRISPRCasFinder, AMRFinderPlus, ISEScan ni bases externas, y no crea ni completa anotaciones.

## Veredicto

**Estado global:** `covered` — El inventario está cubierto por TSV tipados. Las brechas restantes son de evaluabilidad (FASTA/CRISPR) o de métodos externos, no de filas.

| métrica | valor |
|---|---:|
| inventario | `mobile_typing/full_genomes.tsv` |
| genomas esperados (inventario) | 400 |
| genomas descubiertos por artefactos | 400 |
| genomas con TSV tipado | 400 / 400 (100,00 %) |
| genomas con GFF3 | 400 |
| genomas con capa previa | 400 |
| FASTA disponible | 400 / 400 |
| CRISPR evaluable (con FASTA) | 400; no evaluados: 0 |
| arrays CRISPR / espaciadores | 1.749 arrays en 245 genomas / 20.687 espaciadores |
| filas `cas` / genomas con `cas` | 3.288 / 245 |
| filas tipadas | 44.832 |

Cobertura por genoma: completa **400**, parcial **0**, sin datos tipados **0**.

## Inventario y fuentes

El total esperado **no está fijado en el código**: se lee de `mobile_typing/full_genomes.tsv`. Si ese archivo no existe, la cobertura se mide sobre los genomas descubiertos.

| fuente | genomas/archivos |
|---|---:|
| `mobile_typing_full_tsv` | 400 |
| `mobile_typing_pilot_tsv` | 10 |
| `mobile_typing_full_gff` | 400 |
| `mobile_typing_full_spacers` | 400 |
| `mobile_defense_features_genomes` | 400 |
| `annotations_ab_initio_genomes` | 400 |
| `raw_fasta_genomes` | 400 |

## Tipos agregados

Mapeo aplicado (no añade tipos, solo agrupa los existentes):

| categoría | `element_type` incluidos | filas tipadas | genomas tipados | filas previas | genomas previos |
|---|---|---:|---:|---:|---:|
| `mobile_element` | `mobile_element`, `transposon` | 25.142 | 400 | 22.866 | 400 |
| `prophage` | `prophage` | 858 | 265 | 6.981 | 392 |
| `defense` | `defense` | 2.939 | 385 | 33.754 | 400 |
| `crispr` | `crispr`, `cas` | 5.037 | 296 | 3.279 | 247 |
| `amr` | `amr` | 1.051 | 271 | 81.542 | 400 |
| `virulence` | `virulence` | 9.805 | 394 | 13.266 | 398 |

Tipos crudos en el TSV tipado: `transposon`: 25.142, `virulence`: 9.805, `cas`: 3.288, `defense`: 2.939, `crispr`: 1.749, `amr`: 1.051, `prophage`: 858

La capa previa (features.tsv) es keyword-based y no tiene contig ni hebra; sus recuentos se muestran aparte y **no** sustituyen al TSV.

## Pipeline que generó los datos

Atribución estática por scripts/mobile/22_run_full.sh, scripts/mobile/paths.py y los nombres de salida del repo. Los TSV no guardan metadatos del generador.

| cadena | orquestador | comando | activa |
|---|---|---|---|
| full | `scripts/mobile/22_run_full.sh` | `bash scripts/mobile/22_run_full.sh all` | sí |
| pilot | `scripts/mobile/run_all.sh` | `bash scripts/mobile/run_all.sh` | no |

| etapa | script | salida | encontrado | esperado | estado |
|---|---|---|---:|---:|---|
| prepare | `scripts/mobile/19_prepare_full.py` | `mobile_typing/full_genomes.tsv` (1) | 1 | 1 | `present` |
| extract | `scripts/mobile/01_extract_candidates.py` | `mobile_typing/_work/full/*/candidates.tsv` (400) | 400 | 400 | `present` |
| pfam | `scripts/mobile/20_pfam_pool.py` | `mobile_typing/_work/full/*/pfam_hits.tsv` (400) | 400 | 400 | `present` |
| crispr | `scripts/mobile/21_crispr_pool.py` | `mobile_typing/_work/full/*/crispr_arrays.tsv` (310) | 310 | 400 | `partial` |
| classify | `scripts/mobile/05_classify.py` | `mobile_typing/_work/full/*/elements.tsv` (400) | 400 | 400 | `present` |
| write | `scripts/mobile/06_write_outputs.py` | `mobile_typing/full/*.tsv` (400)<br>`mobile_typing/gff/full/*.gff3` (400)<br>`mobile_typing/crispr_spacers/full/*.tsv` (400) | 400/400/400 | 400 | `present` |
| validate | `scripts/mobile/07_validate.py` | `mobile_typing/validation_full.json` (1) | 1 | 1 | `present` |

**Contraste con artefactos declarados**

| artefacto | declarado | observado | coincide |
|---|---:|---:|---|
| `mobile_typing/summary_full.json` (genomas) | 400 | 400 | sí |
| `mobile_typing/summary_full.json` (filas) | 44.832 | 44.832 | sí |
| `mobile_typing/validation_full.json` | status=pass | errors=0 | — |

Fuentes previas (solo lectura; no son la capa tipada):

| sección | script | JSON | aplanado |
|---|---|---|---|
| mobile_element | `scripts/annotation/22_annotate_transposons.py` | `data/mobile_elements.json` | `data/mobile_defense/features.tsv` |
| defense / crispr / prophage | `scripts/annotation/20_annotate_defense.py` | `data/defense_systems.json` | `data/mobile_defense/features.tsv` |
| amr / virulence | `scripts/annotation/21_annotate_amr_virulence.py` | `data/defense_virulence.json` | `data/mobile_defense/features.tsv` |

No se encontró el script que aplana esos JSON en `data/mobile_defense/features.tsv`; se audita el archivo tal como está.

## Brechas

Alcance de las listas: **inventory**.

| brecha | n | detalle |
|---|---:|---|
| TSV tipado ausente | 0 | — |
| GFF3 ausente | 0 | — |
| Espaciadores ausentes con FASTA presente | 0 | — |
| FASTA ausente | 0 | — |
| CRISPR no evaluado | 0 | — |
| Sin capa previa | 0 | — |
| Sin ninguna salida (TSV/GFF/CRISPR/previo) | 0 | — |

Para cubrir el inventario faltan **0** TSV tipados.

## Cobertura por filo

| filo | genomas | con TSV | completa | parcial | sin datos | filas | FASTA |
|---|---:|---:|---:|---:|---:|---:|---:|
| unknown | 42 | 42 | 42 | 0 | 0 | 2.789 | 42 |
| Acidobacteriota | 15 | 15 | 15 | 0 | 0 | 2.055 | 15 |
| Actinomycetota | 15 | 15 | 15 | 0 | 0 | 1.239 | 15 |
| Bacillota | 15 | 15 | 15 | 0 | 0 | 1.596 | 15 |
| Bacteroidota | 15 | 15 | 15 | 0 | 0 | 1.649 | 15 |
| Campylobacterota | 15 | 15 | 15 | 0 | 0 | 933 | 15 |
| Cyanobacteriota | 15 | 15 | 15 | 0 | 0 | 4.931 | 15 |
| Deinococcota | 15 | 15 | 15 | 0 | 0 | 2.244 | 15 |
| Fusobacteriota | 15 | 15 | 15 | 0 | 0 | 1.056 | 15 |
| Myxococcota | 15 | 15 | 15 | 0 | 0 | 3.346 | 15 |
| Pseudomonadota | 15 | 15 | 15 | 0 | 0 | 2.775 | 15 |
| Spirochaetota | 15 | 15 | 15 | 0 | 0 | 651 | 15 |
| Thermodesulfobacteriota | 15 | 15 | 15 | 0 | 0 | 2.412 | 15 |
| Thermotogota | 15 | 15 | 15 | 0 | 0 | 1.089 | 15 |
| Planctomycetota | 14 | 14 | 14 | 0 | 0 | 2.904 | 14 |
| Aquificota | 13 | 13 | 13 | 0 | 0 | 676 | 13 |
| Chloroflexota | 13 | 13 | 13 | 0 | 0 | 1.179 | 13 |
| Verrucomicrobiota | 13 | 13 | 13 | 0 | 0 | 1.369 | 13 |
| Bdellovibrionota | 10 | 10 | 10 | 0 | 0 | 534 | 10 |
| Mycoplasmatota | 10 | 10 | 10 | 0 | 0 | 502 | 10 |
| Nitrospirota | 10 | 10 | 10 | 0 | 0 | 1.076 | 10 |
| Chlorobiota | 9 | 9 | 9 | 0 | 0 | 1.032 | 9 |
| Synergistota | 9 | 9 | 9 | 0 | 0 | 875 | 9 |
| Deferribacterota | 7 | 7 | 7 | 0 | 0 | 873 | 7 |
| Gemmatimonadota | 7 | 7 | 7 | 0 | 0 | 523 | 7 |
| Rhodothermota | 6 | 6 | 6 | 0 | 0 | 552 | 6 |
| Chlamydiota | 4 | 4 | 4 | 0 | 0 | 422 | 4 |
| Elusimicrobiota | 4 | 4 | 4 | 0 | 0 | 407 | 4 |
| Atribacterota | 3 | 3 | 3 | 0 | 0 | 294 | 3 |
| Armatimonadota | 2 | 2 | 2 | 0 | 0 | 316 | 2 |
| Chrysiogenota | 2 | 2 | 2 | 0 | 0 | 225 | 2 |
| Dictyoglomota | 2 | 2 | 2 | 0 | 0 | 89 | 2 |
| Ignavibacteriota | 2 | 2 | 2 | 0 | 0 | 144 | 2 |
| Kiritimatiellota | 2 | 2 | 2 | 0 | 0 | 300 | 2 |
| Minisyncoccota | 2 | 2 | 2 | 0 | 0 | 121 | 2 |
| Thermodesulfobiota | 2 | 2 | 2 | 0 | 0 | 98 | 2 |
| Thermomicrobiota | 2 | 2 | 2 | 0 | 0 | 101 | 2 |
| Abditibacteriota | 1 | 1 | 1 | 0 | 0 | 186 | 1 |
| Balneolota | 1 | 1 | 1 | 0 | 0 | 121 | 1 |
| Caldifarciminota | 1 | 1 | 1 | 0 | 0 | 78 | 1 |
| Caldisericota | 1 | 1 | 1 | 0 | 0 | 43 | 1 |
| Calditrichota | 1 | 1 | 1 | 0 | 0 | 204 | 1 |
| Candidatus Fervidibacterota | 1 | 1 | 1 | 0 | 0 | 96 | 1 |
| Candidatus Saccharimonadota | 1 | 1 | 1 | 0 | 0 | 30 | 1 |
| Coprothermobacterota | 1 | 1 | 1 | 0 | 0 | 23 | 1 |
| Fibrobacterota | 1 | 1 | 1 | 0 | 0 | 66 | 1 |
| Fidelibacterota | 1 | 1 | 1 | 0 | 0 | 51 | 1 |
| Lentisphaerota | 1 | 1 | 1 | 0 | 0 | 207 | 1 |
| Nitrospinota | 1 | 1 | 1 | 0 | 0 | 39 | 1 |
| Thermosulfidibacterota | 1 | 1 | 1 | 0 | 0 | 64 | 1 |
| Vulcanimicrobiota | 1 | 1 | 1 | 0 | 0 | 227 | 1 |
| Zhurongbacterota | 1 | 1 | 1 | 0 | 0 | 20 | 1 |

## Detalle por genoma

`M`=mobile_element, `P`=prophage, `D`=defense, `C`=crispr (incluye cas), `A`=amr, `V`=virulence. `no disponible` lista las capas `not_available`.

<details>
<summary>Ver los 400 genomas</summary>

| genoma | filo | estado | cobertura | filas | M | P | D | C | A | V | FASTA | no disponible |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| `GCF_000006685.1_ASM668v1` | unknown | `available` | `full` | 50 | 6 | 0 | 1 | 0 | 0 | 43 | sí | — |
| `GCF_000006985.1_ASM698v1` | Chlorobiota | `available` | `full` | 66 | 30 | 0 | 7 | 19 | 2 | 8 | sí | — |
| `GCF_000007205.1_ASM720v1` | unknown | `available` | `full` | 42 | 4 | 0 | 1 | 0 | 0 | 37 | sí | — |
| `GCF_000007605.1_ASM760v1` | unknown | `available` | `full` | 42 | 6 | 0 | 1 | 0 | 0 | 35 | sí | — |
| `GCF_000007985.2_ASM798v2` | Thermodesulfobacteriota | `available` | `full` | 126 | 59 | 2 | 6 | 19 | 2 | 38 | sí | — |
| `GCF_000008625.1_ASM862v1` | Aquificota | `available` | `full` | 43 | 9 | 0 | 1 | 20 | 1 | 12 | sí | — |
| `GCF_000008725.1_ASM872v1` | unknown | `available` | `full` | 26 | 4 | 0 | 0 | 0 | 0 | 22 | sí | — |
| `GCF_000010305.1_ASM1030v1` | Gemmatimonadota | `available` | `full` | 81 | 32 | 4 | 7 | 0 | 1 | 37 | sí | — |
| `GCF_000010985.1_ASM1098v1` | Deferribacterota | `available` | `full` | 75 | 38 | 1 | 1 | 18 | 0 | 17 | sí | — |
| `GCF_000011905.1_ASM1190v1` | Chloroflexota | `available` | `full` | 34 | 24 | 3 | 3 | 0 | 0 | 4 | sí | — |
| `GCF_000012485.1_ASM1248v1` | Chlorobiota | `available` | `full` | 69 | 24 | 1 | 13 | 9 | 1 | 21 | sí | — |
| `GCF_000014025.1_ASM1402v1` | Campylobacterota | `available` | `full` | 43 | 11 | 1 | 6 | 0 | 4 | 21 | sí | — |
| `GCF_000014125.1_ASM1412v1` | Bacillota | `available` | `full` | 79 | 47 | 4 | 2 | 11 | 2 | 13 | sí | — |
| `GCF_000015045.1_ASM1504v1` | Thermodesulfobacteriota | `available` | `full` | 182 | 124 | 4 | 13 | 10 | 6 | 25 | sí | — |
| `GCF_000015125.1_ASM1512v1` | Chlorobiota | `available` | `full` | 224 | 159 | 4 | 19 | 26 | 2 | 14 | sí | — |
| `GCF_000017685.1_ASM1768v1` | Spirochaetota | `available` | `full` | 49 | 15 | 0 | 7 | 1 | 2 | 24 | sí | — |
| `GCF_000017805.1_ASM1780v1` | Chloroflexota | `available` | `full` | 212 | 50 | 2 | 18 | 122 | 8 | 12 | sí | — |
| `GCF_000017865.1_ASM1786v1` | Thermotogota | `available` | `full` | 33 | 4 | 0 | 3 | 19 | 0 | 7 | sí | — |
| `GCF_000018605.1_ASM1860v1` | Thermotogota | `available` | `full` | 71 | 45 | 0 | 4 | 12 | 0 | 10 | sí | — |
| `GCF_000018865.1_ASM1886v1` | Chloroflexota | `available` | `full` | 160 | 53 | 2 | 16 | 63 | 6 | 20 | sí | — |
| `GCF_000018945.1_ASM1894v1` | Thermotogota | `available` | `full` | 29 | 4 | 0 | 5 | 11 | 0 | 9 | sí | — |
| `GCF_000019965.1_ASM1996v1` | Verrucomicrobiota | `available` | `full` | 87 | 46 | 1 | 8 | 0 | 4 | 28 | sí | — |
| `GCF_000020145.1_ASM2014v1` | Elusimicrobiota | `available` | `full` | 100 | 8 | 2 | 2 | 3 | 0 | 85 | sí | — |
| `GCF_000020465.1_ASM2046v1` | Chlorobiota | `available` | `full` | 135 | 79 | 3 | 19 | 22 | 3 | 9 | sí | — |
| `GCF_000020505.1_ASM2050v1` | Chlorobiota | `available` | `full` | 40 | 10 | 3 | 3 | 9 | 3 | 12 | sí | — |
| `GCF_000020525.1_ASM2052v1` | Chlorobiota | `available` | `full` | 78 | 25 | 4 | 9 | 28 | 0 | 12 | sí | — |
| `GCF_000020625.1_ASM2062v1` | Chlorobiota | `available` | `full` | 129 | 90 | 8 | 4 | 19 | 3 | 5 | sí | — |
| `GCF_000020645.1_ASM2064v1` | Chlorobiota | `available` | `full` | 175 | 105 | 3 | 31 | 22 | 2 | 12 | sí | — |
| `GCF_000020945.1_ASM2094v1` | Coprothermobacterota | `available` | `full` | 23 | 8 | 0 | 4 | 2 | 0 | 9 | sí | — |
| `GCF_000020965.1_ASM2096v1` | Dictyoglomota | `available` | `full` | 54 | 11 | 0 | 6 | 24 | 0 | 13 | sí | — |
| `GCF_000020985.1_ASM2098v1` | Nitrospirota | `available` | `full` | 91 | 17 | 0 | 12 | 46 | 1 | 15 | sí | — |
| `GCF_000021545.1_ASM2154v1` | Aquificota | `available` | `full` | 76 | 37 | 0 | 2 | 27 | 0 | 10 | sí | — |
| `GCF_000021565.1_ASM2156v1` | Aquificota | `available` | `full` | 51 | 14 | 1 | 3 | 12 | 0 | 21 | sí | — |
| `GCF_000021645.1_ASM2164v1` | Dictyoglomota | `available` | `full` | 35 | 6 | 0 | 0 | 16 | 0 | 13 | sí | — |
| `GCF_000021685.1_ASM2168v1` | Thermomicrobiota | `available` | `full` | 41 | 8 | 0 | 3 | 21 | 2 | 7 | sí | — |
| `GCF_000021905.1_ASM2190v1` | Thermodesulfobacteriota | `available` | `full` | 189 | 134 | 3 | 16 | 9 | 4 | 23 | sí | — |
| `GCF_000021945.1_ASM2194v1` | Chloroflexota | `available` | `full` | 154 | 39 | 6 | 17 | 75 | 4 | 13 | sí | — |
| `GCF_000022145.1_ASM2214v1` | Myxococcota | `available` | `full` | 122 | 60 | 2 | 8 | 12 | 6 | 34 | sí | — |
| `GCF_000022565.1_ASM2256v1` | Acidobacteriota | `available` | `full` | 124 | 73 | 2 | 13 | 12 | 1 | 23 | sí | — |
| `GCF_000023325.1_ASM2332v1` | Thermotogota | `available` | `full` | 126 | 85 | 1 | 7 | 25 | 2 | 6 | sí | — |
| `GCF_000024405.1_ASM2440v1` | Fusobacteriota | `available` | `full` | 91 | 41 | 11 | 2 | 5 | 4 | 28 | sí | — |
| `GCF_000024425.1_ASM2442v1` | Deinococcota | `available` | `full` | 106 | 39 | 1 | 14 | 39 | 2 | 11 | sí | — |
| `GCF_000024845.1_ASM2484v1` | Rhodothermota | `available` | `full` | 91 | 21 | 0 | 10 | 49 | 0 | 11 | sí | — |
| `GCF_000024905.1_ASM2490v1` | Synergistota | `available` | `full` | 24 | 7 | 2 | 5 | 0 | 0 | 10 | sí | — |
| `GCF_000024985.1_ASM2498v1` | Thermomicrobiota | `available` | `full` | 60 | 32 | 2 | 2 | 9 | 2 | 13 | sí | — |
| `GCF_000025005.1_ASM2500v1` | Chloroflexota | `available` | `full` | 51 | 10 | 0 | 0 | 27 | 5 | 9 | sí | — |
| `GCF_000025605.1_ASM2560v1` | Aquificota | `available` | `full` | 35 | 6 | 0 | 0 | 12 | 0 | 17 | sí | — |
| `GCF_000025725.1_ASM2572v1` | Deferribacterota | `available` | `full` | 191 | 159 | 1 | 6 | 4 | 2 | 19 | sí | — |
| `GCF_000025885.1_ASM2588v1` | Synergistota | `available` | `full` | 28 | 11 | 1 | 1 | 6 | 0 | 9 | sí | — |
| `GCF_000025905.1_ASM2590v1` | Verrucomicrobiota | `available` | `full` | 43 | 15 | 0 | 3 | 7 | 0 | 18 | sí | — |
| `GCF_000025945.1_ASM2594v1` | Thermodesulfobacteriota | `available` | `full` | 99 | 59 | 3 | 7 | 5 | 4 | 21 | sí | — |
| `GCF_000091545.1_ASM9154v1` | Deinococcota | `available` | `full` | 90 | 37 | 0 | 3 | 38 | 1 | 11 | sí | — |
| `GCF_000092425.1_ASM9242v1` | Deinococcota | `available` | `full` | 116 | 50 | 2 | 3 | 48 | 1 | 12 | sí | — |
| `GCF_000092785.1_ASM9278v1` | Chlamydiota | `available` | `full` | 138 | 100 | 2 | 4 | 0 | 3 | 29 | sí | — |
| `GCF_000143725.1_ASM14372v1` | Spirochaetota | `available` | `full` | 22 | 6 | 1 | 6 | 1 | 1 | 7 | sí | — |
| `GCF_000146065.2_ASM14606v1` | unknown | `available` | `full` | 76 | 19 | 0 | 8 | 43 | 0 | 6 | sí | — |
| `GCF_000146505.1_ASM14650v1` | Fibrobacterota | `available` | `full` | 66 | 16 | 0 | 11 | 21 | 1 | 17 | sí | — |
| `GCF_000147715.2_ASM14771v3` | Thermotogota | `available` | `full` | 92 | 42 | 0 | 3 | 35 | 1 | 11 | sí | — |
| `GCF_000164905.1_ASM16490v1` | Aquificota | `available` | `full` | 44 | 16 | 2 | 3 | 9 | 0 | 14 | sí | — |
| `GCF_000172155.1_ASM17215v1` | Verrucomicrobiota | `available` | `full` | 245 | 133 | 6 | 7 | 16 | 3 | 80 | sí | — |
| `GCF_000177635.2_ASM17763v2` | Chrysiogenota | `available` | `full` | 134 | 87 | 0 | 9 | 21 | 2 | 15 | sí | — |
| `GCF_000179915.2_ASM17991v2` | Acidobacteriota | `available` | `full` | 116 | 55 | 4 | 15 | 0 | 4 | 38 | sí | — |
| `GCF_000183405.1_ASM18340v1` | Deferribacterota | `available` | `full` | 72 | 27 | 0 | 7 | 22 | 1 | 15 | sí | — |
| `GCF_000185805.1_ASM18580v1` | Aquificota | `available` | `full` | 54 | 8 | 2 | 1 | 28 | 0 | 15 | sí | — |
| `GCF_000186245.1_ASM18624v1` | Campylobacterota | `available` | `full` | 89 | 50 | 6 | 7 | 14 | 0 | 12 | sí | — |
| `GCF_000191045.1_ASM19104v1` | Aquificota | `available` | `full` | 101 | 84 | 0 | 2 | 0 | 0 | 15 | sí | — |
| `GCF_000195335.1_ASM19533v1` | Deinococcota | `available` | `full` | 101 | 47 | 0 | 3 | 45 | 0 | 6 | sí | — |
| `GCF_000195975.1_ASM19597v1` | Cyanobacteriota | `available` | `full` | 38 | 22 | 2 | 4 | 0 | 0 | 10 | sí | — |
| `GCF_000196115.1_ASM19611v1` | Planctomycetota | `available` | `full` | 183 | 117 | 4 | 10 | 2 | 5 | 45 | sí | — |
| `GCF_000196175.1_ASM19617v1` | Bdellovibrionota | `available` | `full` | 38 | 10 | 1 | 2 | 0 | 4 | 21 | sí | — |
| `GCF_000199675.1_ASM19967v1` | Chloroflexota | `available` | `full` | 121 | 56 | 1 | 3 | 42 | 5 | 14 | sí | — |
| `GCF_000204255.1_ASM20425v1` | unknown | `available` | `full` | 49 | 6 | 0 | 1 | 0 | 0 | 42 | sí | — |
| `GCF_000210915.2_ASM21091v2` | Bdellovibrionota | `available` | `full` | 33 | 11 | 0 | 2 | 0 | 2 | 18 | sí | — |
| `GCF_000212395.1_ASM21239v1` | Thermodesulfobiota | `available` | `full` | 71 | 34 | 0 | 2 | 24 | 1 | 10 | sí | — |
| `GCF_000212675.2_ASM21267v2` | Thermodesulfobacteriota | `available` | `full` | 122 | 57 | 10 | 12 | 3 | 3 | 37 | sí | — |
| `GCF_000214155.1_ASM21415v1` | Actinomycetota | `available` | `full` | 65 | 32 | 0 | 10 | 3 | 6 | 14 | sí | — |
| `GCF_000215975.1_ASM21597v1` | Thermodesulfobacteriota | `available` | `full` | 52 | 16 | 0 | 3 | 17 | 2 | 14 | sí | — |
| `GCF_000218625.1_ASM21862v1` | Deferribacterota | `available` | `full` | 246 | 186 | 0 | 3 | 36 | 2 | 19 | sí | — |
| `GCF_000219725.1_ASM21972v1` | Spirochaetota | `available` | `full` | 87 | 57 | 1 | 8 | 8 | 1 | 12 | sí | — |
| `GCF_000225345.1_ASM22534v1` | Bacillota | `available` | `full` | 142 | 97 | 1 | 7 | 11 | 7 | 19 | sí | — |
| `GCF_000226295.1_ASM22629v1` | Acidobacteriota | `available` | `full` | 112 | 31 | 0 | 20 | 40 | 1 | 20 | sí | — |
| `GCF_000226315.1_ASM22631v1` | Bdellovibrionota | `available` | `full` | 48 | 20 | 3 | 4 | 0 | 3 | 18 | sí | — |
| `GCF_000235405.2_ASM23540v3` | Thermotogota | `available` | `full` | 93 | 55 | 2 | 3 | 21 | 0 | 12 | sí | — |
| `GCF_000237205.1_ASM23720v1` | Chlamydiota | `available` | `full` | 111 | 77 | 0 | 1 | 1 | 2 | 30 | sí | — |
| `GCF_000239475.1_ASM23947v1` | Spirochaetota | `available` | `full` | 9 | 2 | 0 | 1 | 0 | 0 | 6 | sí | — |
| `GCF_000252445.1_ASM25244v1` | Deinococcota | `available` | `full` | 169 | 109 | 6 | 8 | 19 | 6 | 21 | sí | — |
| `GCF_000253035.1_ASM25303v1` | Chlamydiota | `available` | `full` | 128 | 79 | 0 | 2 | 0 | 8 | 39 | sí | — |
| `GCF_000255135.1_ASM25513v1` | Thermotogota | `available` | `full` | 54 | 19 | 1 | 3 | 24 | 0 | 7 | sí | — |
| `GCF_000258405.1_ASM25840v1` | Ignavibacteriota | `available` | `full` | 84 | 47 | 2 | 3 | 20 | 1 | 11 | sí | — |
| `GCF_000259275.1_ASM25927v1` | Campylobacterota | `available` | `full` | 102 | 61 | 0 | 24 | 4 | 1 | 12 | sí | — |
| `GCF_000265295.1_ASM26529v1` | Campylobacterota | `available` | `full` | 50 | 26 | 0 | 7 | 8 | 1 | 8 | sí | — |
| `GCF_000265425.1_ASM26542v1` | Acidobacteriota | `available` | `full` | 76 | 29 | 4 | 8 | 4 | 3 | 28 | sí | — |
| `GCF_000279145.1_ASM27914v1` | Ignavibacteriota | `available` | `full` | 60 | 25 | 1 | 5 | 16 | 0 | 13 | sí | — |
| `GCF_000284115.1_ASM28411v1` | Planctomycetota | `available` | `full` | 100 | 37 | 0 | 2 | 26 | 0 | 35 | sí | — |
| `GCF_000284315.1_ASM28431v1` | Nitrospirota | `available` | `full` | 119 | 74 | 5 | 10 | 0 | 1 | 29 | sí | — |
| `GCF_000284335.1_ASM28433v1` | Caldisericota | `available` | `full` | 43 | 24 | 1 | 1 | 10 | 0 | 7 | sí | — |
| `GCF_000309885.1_ASM30988v1` | Deinococcota | `available` | `full` | 88 | 51 | 0 | 1 | 26 | 0 | 10 | sí | — |
| `GCF_000317105.1_ASM31710v1` | Cyanobacteriota | `available` | `full` | 266 | 135 | 2 | 26 | 73 | 5 | 25 | sí | — |
| `GCF_000317575.1_ASM31757v1` | Cyanobacteriota | `available` | `full` | 185 | 98 | 0 | 25 | 28 | 2 | 32 | sí | — |
| `GCF_000348725.1_ASM34872v1` | Bdellovibrionota | `available` | `full` | 41 | 15 | 0 | 3 | 4 | 3 | 16 | sí | — |
| `GCF_000392435.1_ASM39243v1` | Candidatus Saccharimonadota | `available` | `full` | 30 | 15 | 0 | 1 | 0 | 0 | 14 | sí | — |
| `GCF_000439435.1_ASM43943v1` | Mycoplasmatota | `available` | `full` | 13 | 9 | 0 | 0 | 4 | 0 | 0 | sí | — |
| `GCF_000455605.1_GPTSA100-9` | Bacteroidota | `available` | `full` | 52 | 31 | 0 | 11 | 2 | 1 | 7 | sí | — |
| `GCF_000471025.2_ASM47102v2` | unknown | `available` | `full` | 47 | 7 | 0 | 1 | 0 | 0 | 39 | sí | — |
| `GCF_000479335.1_ASM47933v1` | Bacillota | `available` | `full` | 70 | 27 | 1 | 1 | 4 | 8 | 29 | sí | — |
| `GCF_000484535.1_ASM48453v1` | Cyanobacteriota | `available` | `full` | 153 | 83 | 3 | 13 | 30 | 1 | 23 | sí | — |
| `GCF_000487995.2_ASM48799v2` | Deferribacterota | `available` | `full` | 126 | 76 | 2 | 16 | 6 | 1 | 25 | sí | — |
| `GCF_000504085.1_ASM50408v1` | Thermotogota | `available` | `full` | 37 | 5 | 0 | 2 | 23 | 0 | 7 | sí | — |
| `GCF_000508245.1_ASM50824v1` | unknown | `available` | `full` | 15 | 1 | 0 | 14 | 0 | 0 | 0 | sí | — |
| `GCF_000512735.1_ASM51273v1` | Aquificota | `available` | `full` | 47 | 13 | 1 | 1 | 20 | 0 | 12 | sí | — |
| `GCF_000513475.1_Bbla` | unknown | `available` | `full` | 14 | 7 | 0 | 0 | 0 | 2 | 5 | sí | — |
| `GCF_000522985.1_ASM52298v1` | Gemmatimonadota | `available` | `full` | 67 | 14 | 5 | 3 | 0 | 4 | 41 | sí | — |
| `GCF_000582535.1_ASM58253v1` | unknown | `available` | `full` | 8 | 3 | 0 | 2 | 3 | 0 | 0 | sí | — |
| `GCF_000695095.2_ASM69509v2` | Gemmatimonadota | `available` | `full` | 60 | 14 | 1 | 1 | 1 | 2 | 41 | sí | — |
| `GCF_000724625.1_ASM72462v1` | Armatimonadota | `available` | `full` | 140 | 26 | 1 | 14 | 0 | 3 | 96 | sí | — |
| `GCF_000767275.3_ASM76727v4` | Thermotogota | `available` | `full` | 130 | 98 | 2 | 5 | 14 | 0 | 11 | sí | — |
| `GCF_000800395.1_ASM80039v1` | Deinococcota | `available` | `full` | 90 | 64 | 2 | 5 | 0 | 9 | 10 | sí | — |
| `GCF_000803625.1_ASM80362v1` | Minisyncoccota | `available` | `full` | 21 | 7 | 1 | 3 | 0 | 0 | 10 | sí | — |
| `GCF_000828675.1_ASM82867v1` | Thermotogota | `available` | `full` | 55 | 4 | 0 | 3 | 39 | 0 | 9 | sí | — |
| `GCF_000829235.1_ASM82923v1` | Cyanobacteriota | `available` | `full` | 12 | 8 | 0 | 0 | 0 | 0 | 4 | sí | — |
| `GCF_000941075.1_ASM94107v1` | Mycoplasmatota | `available` | `full` | 84 | 46 | 0 | 3 | 4 | 0 | 31 | sí | — |
| `GCF_001017655.1_ASM101765v1` | Kiritimatiellota | `available` | `full` | 174 | 104 | 1 | 21 | 27 | 0 | 21 | sí | — |
| `GCF_001027025.1_ASM102702v1` | Thermotogota | `available` | `full` | 30 | 12 | 1 | 4 | 0 | 6 | 7 | sí | — |
| `GCF_001027285.1_ASM102728v1` | Myxococcota | `available` | `full` | 170 | 50 | 5 | 8 | 1 | 8 | 98 | sí | — |
| `GCF_001027545.1_ASM102754v1` | Elusimicrobiota | `available` | `full` | 40 | 9 | 0 | 3 | 13 | 0 | 15 | sí | — |
| `GCF_001273775.1_ASM127377v1` | Nitrospirota | `available` | `full` | 132 | 71 | 0 | 26 | 11 | 2 | 22 | sí | — |
| `GCF_001296145.1_ASM129614v1` | Fusobacteriota | `available` | `full` | 54 | 20 | 2 | 18 | 8 | 0 | 6 | sí | — |
| `GCF_001458695.1_NiCh1` | unknown | `available` | `full` | 95 | 24 | 2 | 11 | 30 | 2 | 26 | sí | — |
| `GCF_001547735.1_ASM154773v1` | Thermosulfidibacterota | `available` | `full` | 64 | 36 | 0 | 0 | 8 | 0 | 20 | sí | — |
| `GCF_001610955.1_ASM161095v1` | Bacillota | `available` | `full` | 161 | 113 | 3 | 15 | 6 | 3 | 21 | sí | — |
| `GCF_001611795.1_ASM161179v1` | Actinomycetota | `available` | `full` | 192 | 115 | 4 | 13 | 8 | 16 | 36 | sí | — |
| `GCF_001618865.1_ASM161886v1` | Acidobacteriota | `available` | `full` | 160 | 84 | 8 | 6 | 3 | 9 | 50 | sí | — |
| `GCF_001676785.2_ASM167678v2` | Spirochaetota | `available` | `full` | 41 | 23 | 1 | 7 | 2 | 1 | 7 | sí | — |
| `GCF_001701045.1_ASM170104v1` | Synergistota | `available` | `full` | 126 | 82 | 2 | 8 | 16 | 2 | 16 | sí | — |
| `GCF_001712875.1_ASM171287v1` | unknown | `available` | `full` | 18 | 11 | 0 | 5 | 0 | 1 | 1 | sí | — |
| `GCF_001747405.1_ASM174740v1` | Chlorobiota | `available` | `full` | 116 | 70 | 1 | 13 | 22 | 2 | 8 | sí | — |
| `GCF_001747425.1_ASM174742v1` | Actinomycetota | `available` | `full` | 82 | 15 | 1 | 7 | 28 | 13 | 18 | sí | — |
| `GCF_001886435.1_ASM188643v1` | Pseudomonadota | `available` | `full` | 86 | 29 | 0 | 6 | 11 | 2 | 38 | sí | — |
| `GCF_001886815.1_ASM188681v1` | Calditrichota | `available` | `full` | 204 | 137 | 1 | 10 | 45 | 0 | 11 | sí | — |
| `GCF_001889025.1_ASM188902v1` | Pseudomonadota | `available` | `full` | 48 | 18 | 1 | 5 | 2 | 5 | 17 | sí | — |
| `GCF_001953175.1_ASM195317v1` | Chloroflexota | `available` | `full` | 50 | 37 | 0 | 4 | 1 | 3 | 5 | sí | — |
| `GCF_002139875.1_ASM213987v1` | Campylobacterota | `available` | `full` | 50 | 20 | 1 | 4 | 13 | 3 | 9 | sí | — |
| `GCF_002157165.1_ASM215716v1` | Pseudomonadota | `available` | `full` | 153 | 108 | 3 | 12 | 2 | 5 | 23 | sí | — |
| `GCF_002269385.1_ASM226938v1` | Bacteroidota | `available` | `full` | 128 | 106 | 2 | 8 | 0 | 1 | 11 | sí | — |
| `GCF_002302565.1_ASM230256v1` | Bacteroidota | `available` | `full` | 55 | 36 | 0 | 5 | 6 | 0 | 8 | sí | — |
| `GCF_002305895.1_ASM230589v1` | Myxococcota | `available` | `full` | 254 | 134 | 7 | 14 | 21 | 7 | 71 | sí | — |
| `GCF_002310495.1_ASM231049v1` | Verrucomicrobiota | `available` | `full` | 50 | 21 | 0 | 4 | 1 | 2 | 22 | sí | — |
| `GCF_002355835.1_ASM235583v1` | unknown | `available` | `full` | 29 | 4 | 0 | 13 | 10 | 0 | 2 | sí | — |
| `GCF_002804025.1_ASM280402v1` | Mycoplasmatota | `available` | `full` | 38 | 36 | 1 | 1 | 0 | 0 | 0 | sí | — |
| `GCF_002813575.1_ASM281357v1` | Cyanobacteriota | `available` | `full` | 947 | 823 | 10 | 35 | 23 | 8 | 48 | sí | — |
| `GCF_002872415.1_ASM287241v1` | Bdellovibrionota | `available` | `full` | 54 | 21 | 0 | 5 | 0 | 3 | 25 | sí | — |
| `GCF_003019315.1_ASM301931v1` | Fusobacteriota | `available` | `full` | 79 | 48 | 4 | 3 | 12 | 1 | 11 | sí | — |
| `GCF_003019655.1_ASM301965v1` | Fusobacteriota | `available` | `full` | 52 | 25 | 3 | 2 | 13 | 0 | 9 | sí | — |
| `GCF_003019675.1_ASM301967v1` | Fusobacteriota | `available` | `full` | 62 | 28 | 4 | 5 | 13 | 2 | 10 | sí | — |
| `GCF_003057965.1_ASM305796v1` | Thermodesulfobiota | `available` | `full` | 27 | 13 | 0 | 1 | 2 | 1 | 10 | sí | — |
| `GCF_003173015.1_ASM317301v1` | Deinococcota | `available` | `full` | 78 | 43 | 6 | 3 | 0 | 11 | 15 | sí | — |
| `GCF_003258315.1_ASM325831v1` | unknown | `available` | `full` | 98 | 34 | 0 | 9 | 27 | 2 | 26 | sí | — |
| `GCF_003325455.1_ASM332545v1` | Actinomycetota | `available` | `full` | 96 | 48 | 0 | 14 | 14 | 4 | 16 | sí | — |
| `GCF_003330725.1_ASM333072v1` | Acidobacteriota | `available` | `full` | 290 | 224 | 7 | 12 | 2 | 5 | 40 | sí | — |
| `GCF_003351905.1_ASM335190v1` | Bacteroidota | `available` | `full` | 90 | 85 | 1 | 1 | 0 | 0 | 3 | sí | — |
| `GCF_003353065.1_ASM335306v1` | Balneolota | `available` | `full` | 121 | 69 | 0 | 20 | 18 | 0 | 14 | sí | — |
| `GCF_003355515.1_ASM335551v1` | Campylobacterota | `available` | `full` | 60 | 47 | 0 | 3 | 0 | 1 | 9 | sí | — |
| `GCF_003366055.1_ASM336605v1` | unknown | `available` | `full` | 13 | 8 | 0 | 0 | 0 | 3 | 2 | sí | — |
| `GCF_003491405.1_ASM349140v1` | Rhodothermota | `available` | `full` | 140 | 99 | 4 | 15 | 14 | 1 | 7 | sí | — |
| `GCF_003544815.1_ASM354481v1` | Campylobacterota | `available` | `full` | 53 | 27 | 1 | 6 | 1 | 3 | 15 | sí | — |
| `GCF_003574315.2_ASM357431v2` | Armatimonadota | `available` | `full` | 176 | 76 | 8 | 11 | 17 | 6 | 58 | sí | — |
| `GCF_003721155.3_ASM372115v4` | Actinomycetota | `available` | `full` | 36 | 10 | 1 | 3 | 2 | 3 | 17 | sí | — |
| `GCF_003815875.1_ASM381587v1` | Bacteroidota | `available` | `full` | 166 | 67 | 1 | 7 | 5 | 8 | 78 | sí | — |
| `GCF_003966735.1_ASM396673v1` | Thermodesulfobacteriota | `available` | `full` | 56 | 19 | 2 | 6 | 2 | 3 | 24 | sí | — |
| `GCF_003966975.1_ASM396697v1` | Chloroflexota | `available` | `full` | 36 | 15 | 0 | 1 | 10 | 1 | 9 | sí | — |
| `GCF_004087915.1_ASM408791v1` | Deferribacterota | `available` | `full` | 81 | 20 | 0 | 4 | 28 | 3 | 26 | sí | — |
| `GCF_004102945.1_ASM410294v1` | unknown | `available` | `full` | 37 | 13 | 2 | 2 | 0 | 0 | 20 | sí | — |
| `GCF_004103715.1_ASM410371v1` | Bacillota | `available` | `full` | 230 | 146 | 3 | 20 | 40 | 3 | 18 | sí | — |
| `GCF_004214795.1_ASM421479v1` | Campylobacterota | `available` | `full` | 124 | 101 | 1 | 4 | 5 | 3 | 10 | sí | — |
| `GCF_004376055.1_ASM437605v2` | unknown | `available` | `full` | 101 | 57 | 2 | 15 | 5 | 2 | 20 | sí | — |
| `GCF_004803815.1_ASM480381v1` | Campylobacterota | `available` | `full` | 43 | 15 | 1 | 3 | 9 | 1 | 14 | sí | — |
| `GCF_005144435.2_ASM514443v2` | Pseudomonadota | `available` | `full` | 199 | 99 | 6 | 14 | 0 | 13 | 67 | sí | — |
| `GCF_005697565.1_ASM569756v1` | unknown | `available` | `full` | 16 | 6 | 0 | 2 | 0 | 0 | 8 | sí | — |
| `GCF_006352065.1_ASM635206v1` | Actinomycetota | `available` | `full` | 69 | 34 | 1 | 8 | 0 | 5 | 21 | sí | — |
| `GCF_007475525.1_ASM747552v1` | Verrucomicrobiota | `available` | `full` | 29 | 10 | 0 | 1 | 12 | 1 | 5 | sí | — |
| `GCF_007747655.1_ASM774765v1` | Planctomycetota | `available` | `full` | 212 | 140 | 6 | 8 | 1 | 4 | 53 | sí | — |
| `GCF_007750655.1_ASM775065v1` | Planctomycetota | `available` | `full` | 36 | 16 | 1 | 2 | 0 | 1 | 16 | sí | — |
| `GCF_007904085.1_ASM790408v1` | Cyanobacteriota | `available` | `full` | 141 | 87 | 0 | 15 | 17 | 2 | 20 | sí | — |
| `GCF_007990365.1_ASM799036v1` | Fusobacteriota | `available` | `full` | 145 | 92 | 2 | 8 | 22 | 1 | 20 | sí | — |
| `GCF_007990505.1_ASM799050v1` | Fusobacteriota | `available` | `full` | 68 | 29 | 1 | 10 | 11 | 0 | 17 | sí | — |
| `GCF_007993755.1_ASM799375v1` | unknown | `available` | `full` | 74 | 19 | 0 | 17 | 8 | 2 | 28 | sí | — |
| `GCF_007993775.1_ASM799377v1` | unknown | `available` | `full` | 196 | 116 | 3 | 15 | 30 | 3 | 29 | sí | — |
| `GCF_008086545.1_ASM808654v1` | Mycoplasmatota | `available` | `full` | 15 | 14 | 0 | 1 | 0 | 0 | 0 | sí | — |
| `GCF_008087625.1_ASM808762v1` | Planctomycetota | `available` | `full` | 233 | 155 | 4 | 11 | 13 | 5 | 45 | sí | — |
| `GCF_008087665.1_ASM808766v1` | Planctomycetota | `available` | `full` | 149 | 100 | 0 | 9 | 11 | 1 | 28 | sí | — |
| `GCF_009002475.1_ASM900247v1` | unknown | `available` | `full` | 107 | 60 | 3 | 4 | 12 | 3 | 25 | sí | — |
| `GCF_009258225.1_ASM925822v1` | Campylobacterota | `available` | `full` | 54 | 13 | 2 | 3 | 4 | 2 | 30 | sí | — |
| `GCF_009688965.1_ASM968896v1` | Thermodesulfobacteriota | `available` | `full` | 270 | 191 | 12 | 18 | 0 | 4 | 45 | sí | — |
| `GCF_009689585.1_ASM968958v1` | Bacteroidota | `available` | `full` | 229 | 185 | 3 | 22 | 0 | 0 | 19 | sí | — |
| `GCF_009731575.1_ASM973157v1` | Verrucomicrobiota | `available` | `full` | 45 | 17 | 5 | 10 | 0 | 3 | 10 | sí | — |
| `GCF_009734425.1_ASM973442v1` | Actinomycetota | `available` | `full` | 95 | 65 | 1 | 8 | 0 | 5 | 16 | sí | — |
| `GCF_009905215.1_ASM990521v1` | Kiritimatiellota | `available` | `full` | 126 | 87 | 6 | 3 | 13 | 2 | 15 | sí | — |
| `GCF_010726765.1_ASM1072676v1` | Actinomycetota | `available` | `full` | 83 | 41 | 4 | 2 | 0 | 9 | 27 | sí | — |
| `GCF_011751885.2_ASM1175188v2` | Actinomycetota | `available` | `full` | 59 | 33 | 2 | 4 | 1 | 3 | 16 | sí | — |
| `GCF_013004105.1_ASM1300410v1` | Gemmatimonadota | `available` | `full` | 83 | 31 | 1 | 6 | 1 | 1 | 43 | sí | — |
| `GCF_013128195.2_ASM1312819v2` | Planctomycetota | `available` | `full` | 450 | 387 | 3 | 5 | 16 | 2 | 37 | sí | — |
| `GCF_013177655.1_ASM1317765v1` | Campylobacterota | `available` | `full` | 55 | 5 | 2 | 13 | 14 | 1 | 20 | sí | — |
| `GCF_013340845.1_ASM1334084v1` | Pseudomonadota | `available` | `full` | 54 | 36 | 1 | 3 | 0 | 4 | 10 | sí | — |
| `GCF_014274685.1_ASM1427468v1` | Acidobacteriota | `available` | `full` | 57 | 23 | 1 | 5 | 0 | 6 | 22 | sí | — |
| `GCF_014337215.1_ASM1433721v1` | Bacillota | `available` | `full` | 154 | 103 | 6 | 3 | 18 | 4 | 20 | sí | — |
| `GCF_014337255.1_ASM1433725v1` | Fusobacteriota | `available` | `full` | 82 | 47 | 1 | 3 | 12 | 10 | 9 | sí | — |
| `GCF_014529625.1_ASM1452962v2` | unknown | `available` | `full` | 29 | 14 | 2 | 0 | 0 | 0 | 13 | sí | — |
| `GCF_014672715.2_ASM1467271v2` | Chloroflexota | `available` | `full` | 38 | 23 | 3 | 2 | 0 | 3 | 7 | sí | — |
| `GCF_014701235.1_ASM1470123v1` | Nitrospirota | `available` | `full` | 113 | 45 | 1 | 12 | 43 | 0 | 12 | sí | — |
| `GCF_015277775.1_ASM1527777v1` | Acidobacteriota | `available` | `full` | 153 | 51 | 5 | 17 | 11 | 6 | 63 | sí | — |
| `GCF_015356815.2_ASM1535681v2` | unknown | `available` | `full` | 79 | 59 | 0 | 1 | 0 | 0 | 19 | sí | — |
| `GCF_015500095.1_ASM1550009v1` | Chloroflexota | `available` | `full` | 98 | 54 | 3 | 4 | 17 | 5 | 15 | sí | — |
| `GCF_015775515.1_ASM1577551v1` | Atribacterota | `available` | `full` | 143 | 107 | 0 | 11 | 18 | 0 | 7 | sí | — |
| `GCF_016027375.1_ASM1602737v1` | Bacillota | `available` | `full` | 54 | 24 | 1 | 6 | 0 | 1 | 22 | sí | — |
| `GCF_016446285.1_ASM1644628v1` | Thermodesulfobacteriota | `available` | `full` | 134 | 73 | 9 | 10 | 20 | 2 | 20 | sí | — |
| `GCF_016724785.1_ASM1672478v1` | Fusobacteriota | `available` | `full` | 74 | 54 | 0 | 5 | 9 | 0 | 6 | sí | — |
| `GCF_017498545.1_ASM1749854v1` | Acidobacteriota | `available` | `full` | 348 | 143 | 13 | 7 | 66 | 1 | 118 | sí | — |
| `GCF_017893965.1_ASM1789396v1` | Spirochaetota | `available` | `full` | 46 | 19 | 1 | 11 | 5 | 0 | 10 | sí | — |
| `GCF_018069705.1_ASM1806970v1` | Synergistota | `available` | `full` | 150 | 100 | 1 | 3 | 28 | 0 | 18 | sí | — |
| `GCF_018127985.1_ASM1812798v1` | Bacteroidota | `available` | `full` | 151 | 104 | 0 | 10 | 23 | 0 | 14 | sí | — |
| `GCF_018137965.1_ASM1813796v1` | Verrucomicrobiota | `available` | `full` | 79 | 43 | 1 | 5 | 6 | 4 | 20 | sí | — |
| `GCF_018243235.1_ASM1824323v1` | Pseudomonadota | `available` | `full` | 322 | 142 | 4 | 10 | 1 | 14 | 151 | sí | — |
| `GCF_018292205.1_ASM1829220v1` | Bacteroidota | `available` | `full` | 289 | 218 | 3 | 10 | 3 | 4 | 51 | sí | — |
| `GCF_018304765.1_ASM1830476v1` | Acidobacteriota | `available` | `full` | 95 | 28 | 0 | 19 | 26 | 1 | 21 | sí | — |
| `GCF_018343815.1_ASM1834381v1` | unknown | `available` | `full` | 49 | 6 | 0 | 1 | 0 | 0 | 42 | sí | — |
| `GCF_018398915.1_ASM1839891v1` | Synergistota | `available` | `full` | 166 | 107 | 2 | 8 | 32 | 0 | 17 | sí | — |
| `GCF_018771645.1_ASM1877164v1` | Deinococcota | `available` | `full` | 93 | 68 | 2 | 2 | 11 | 2 | 8 | sí | — |
| `GCF_018866205.1_ASM1886620v1` | unknown | `available` | `full` | 19 | 8 | 0 | 4 | 0 | 0 | 7 | sí | — |
| `GCF_019048065.1_ASM1904806v1` | Pseudomonadota | `available` | `full` | 227 | 42 | 4 | 9 | 3 | 35 | 134 | sí | — |
| `GCF_019048105.1_ASM1904810v1` | Pseudomonadota | `available` | `full` | 304 | 87 | 15 | 6 | 5 | 29 | 162 | sí | — |
| `GCF_019056495.1_ASM1905649v1` | unknown | `available` | `full` | 51 | 6 | 0 | 1 | 0 | 0 | 44 | sí | — |
| `GCF_019063885.1_ASM1906388v1` | Thermodesulfobacteriota | `available` | `full` | 97 | 29 | 4 | 6 | 10 | 7 | 41 | sí | — |
| `GCF_019552005.1_ASM1955200v1` | Fusobacteriota | `available` | `full` | 35 | 13 | 0 | 12 | 0 | 0 | 10 | sí | — |
| `GCF_020091585.1_ASM2009158v1` | Cyanobacteriota | `available` | `full` | 617 | 476 | 4 | 36 | 45 | 7 | 49 | sí | — |
| `GCF_020549965.1_ASM2054996v1` | unknown | `available` | `full` | 49 | 6 | 0 | 1 | 0 | 0 | 42 | sí | — |
| `GCF_020616335.1_ASM2061633v1` | Actinomycetota | `available` | `full` | 37 | 12 | 0 | 1 | 9 | 1 | 14 | sí | — |
| `GCF_021057185.1_ASM2105718v1` | unknown | `available` | `full` | 79 | 36 | 14 | 3 | 18 | 0 | 8 | sí | — |
| `GCF_021117425.1_ASM2111742v1` | Actinomycetota | `available` | `full` | 54 | 25 | 0 | 4 | 0 | 9 | 16 | sí | — |
| `GCF_021379005.1_ASM2137900v1` | Cyanobacteriota | `available` | `full` | 404 | 335 | 0 | 18 | 26 | 2 | 23 | sí | — |
| `GCF_022026155.1_ASM2202615v1` | unknown | `available` | `full` | 54 | 19 | 1 | 7 | 9 | 2 | 16 | sí | — |
| `GCF_022376835.1_ASM2237683v1` | Spirochaetota | `available` | `full` | 57 | 23 | 1 | 5 | 8 | 4 | 16 | sí | — |
| `GCF_022813675.1_ASM2281367v1` | Deferribacterota | `available` | `full` | 82 | 45 | 2 | 4 | 15 | 1 | 15 | sí | — |
| `GCF_022814725.1_ASM2281472v1` | Bdellovibrionota | `available` | `full` | 37 | 8 | 1 | 3 | 0 | 1 | 24 | sí | — |
| `GCF_022846135.1_ASM2284613v1` | Synergistota | `available` | `full` | 102 | 54 | 4 | 9 | 27 | 0 | 8 | sí | — |
| `GCF_022869165.1_ASM2286916v1` | Spirochaetota | `available` | `full` | 78 | 55 | 1 | 9 | 1 | 1 | 11 | sí | — |
| `GCF_022870745.1_ASM2287074v1` | Acidobacteriota | `available` | `full` | 98 | 41 | 4 | 10 | 2 | 4 | 37 | sí | — |
| `GCF_023035795.1_ASM2303579v1` | Spirochaetota | `available` | `full` | 45 | 7 | 4 | 7 | 0 | 0 | 27 | sí | — |
| `GCF_023035875.1_ASM2303587v1` | Spirochaetota | `available` | `full` | 53 | 12 | 7 | 17 | 0 | 0 | 17 | sí | — |
| `GCF_023169805.1_ASM2316980v1` | Abditibacteriota | `available` | `full` | 186 | 137 | 3 | 11 | 1 | 4 | 30 | sí | — |
| `GCF_023921165.1_ASM2392116v1` | Bacteroidota | `available` | `full` | 27 | 8 | 0 | 4 | 3 | 2 | 10 | sí | — |
| `GCF_023935825.1_ASM2393582v1` | Pseudomonadota | `available` | `full` | 114 | 58 | 3 | 1 | 1 | 5 | 46 | sí | — |
| `GCF_024205945.1_ASM2420594v1` | Bacillota | `available` | `full` | 85 | 37 | 3 | 9 | 14 | 9 | 13 | sí | — |
| `GCF_024338845.1_ASM2433884v1` | unknown | `available` | `full` | 22 | 9 | 0 | 5 | 0 | 0 | 8 | sí | — |
| `GCF_024584745.1_ASM2458474v1` | Cyanobacteriota | `available` | `full` | 268 | 122 | 1 | 36 | 85 | 3 | 21 | sí | — |
| `GCF_025148635.1_ASM2514863v1` | Bacteroidota | `available` | `full` | 203 | 140 | 0 | 12 | 6 | 2 | 43 | sí | — |
| `GCF_025273655.1_ASM2527365v1` | unknown | `available` | `full` | 30 | 11 | 0 | 3 | 0 | 0 | 16 | sí | — |
| `GCF_025558325.1_ASM2555832v1` | Bacillota | `available` | `full` | 111 | 75 | 3 | 5 | 7 | 7 | 14 | sí | — |
| `GCF_025725725.1_ASM2572572v1` | Actinomycetota | `available` | `full` | 98 | 54 | 3 | 3 | 10 | 9 | 19 | sí | — |
| `GCF_025997855.1_ASM2599785v1` | Deinococcota | `available` | `full` | 251 | 189 | 5 | 8 | 22 | 15 | 12 | sí | — |
| `GCF_025997995.1_ASM2599799v1` | Bacteroidota | `available` | `full` | 37 | 19 | 0 | 5 | 0 | 3 | 10 | sí | — |
| `GCF_026184435.1_ASM2618443v1` | Zhurongbacterota | `available` | `full` | 20 | 3 | 0 | 5 | 0 | 0 | 12 | sí | — |
| `GCF_026240635.1_ASM2624063v1` | Cyanobacteriota | `available` | `full` | 140 | 67 | 0 | 12 | 44 | 1 | 16 | sí | — |
| `GCF_026965555.1_ASM2696555v1` | Myxococcota | `available` | `full` | 297 | 94 | 7 | 12 | 102 | 6 | 76 | sí | — |
| `GCF_027594505.1_ASM2759450v1` | Chloroflexota | `available` | `full` | 35 | 6 | 1 | 5 | 10 | 2 | 11 | sí | — |
| `GCF_027625475.1_ASM2762547v1` | Pseudomonadota | `available` | `full` | 180 | 132 | 5 | 5 | 0 | 7 | 31 | sí | — |
| `GCF_027923555.1_ASM2792355v1` | Vulcanimicrobiota | `available` | `full` | 227 | 204 | 5 | 1 | 0 | 1 | 16 | sí | — |
| `GCF_028622095.1_ASM2862209v1` | Deinococcota | `available` | `full` | 95 | 59 | 3 | 1 | 12 | 5 | 15 | sí | — |
| `GCF_028728065.1_ASM2872806v1` | Lentisphaerota | `available` | `full` | 207 | 130 | 1 | 8 | 16 | 0 | 52 | sí | — |
| `GCF_028885455.1_ASM2888545v1` | Pseudomonadota | `available` | `full` | 334 | 235 | 5 | 10 | 25 | 4 | 55 | sí | — |
| `GCF_029023725.1_ASM2902372v1` | Verrucomicrobiota | `available` | `full` | 146 | 79 | 0 | 16 | 13 | 1 | 37 | sí | — |
| `GCF_029677705.1_ASM2967770v1` | unknown | `available` | `full` | 16 | 8 | 0 | 5 | 0 | 0 | 3 | sí | — |
| `GCF_029919255.1_ASM2991925v1` | Cyanobacteriota | `available` | `full` | 17 | 13 | 1 | 0 | 0 | 0 | 3 | sí | — |
| `GCF_030169785.1_ASM3016978v1` | Bdellovibrionota | `available` | `full` | 41 | 10 | 1 | 2 | 0 | 4 | 24 | sí | — |
| `GCF_030253535.1_ASM3025353v1` | unknown | `available` | `full` | 25 | 12 | 0 | 4 | 4 | 0 | 5 | sí | — |
| `GCF_030295385.1_ASM3029538v1` | Acidobacteriota | `available` | `full` | 47 | 10 | 1 | 9 | 1 | 8 | 18 | sí | — |
| `GCF_030389965.1_ASM3038996v1` | unknown | `available` | `full` | 310 | 258 | 8 | 13 | 14 | 2 | 15 | sí | — |
| `GCF_030438115.2_ASM3043811v2` | Spirochaetota | `available` | `full` | 40 | 10 | 10 | 11 | 0 | 0 | 9 | sí | — |
| `GCF_030439285.2_ASM3043928v2` | Spirochaetota | `available` | `full` | 29 | 10 | 3 | 2 | 0 | 0 | 14 | sí | — |
| `GCF_030444335.1_ASM3044433v1` | Bacteroidota | `available` | `full` | 25 | 13 | 0 | 3 | 0 | 1 | 8 | sí | — |
| `GCF_030503875.2_ASM3050387v2` | Bacteroidota | `available` | `full` | 39 | 21 | 0 | 5 | 6 | 0 | 7 | sí | — |
| `GCF_030520105.1_ASM3052010v1` | Candidatus Fervidibacterota | `available` | `full` | 96 | 19 | 1 | 13 | 39 | 0 | 24 | sí | — |
| `GCF_030585795.1_ASM3058579v1` | Bacillota | `available` | `full` | 77 | 52 | 4 | 5 | 0 | 2 | 14 | sí | — |
| `GCF_030623445.1_ASM3062344v1` | Gemmatimonadota | `available` | `full` | 102 | 71 | 0 | 13 | 0 | 0 | 18 | sí | — |
| `GCF_030866825.1_ASM3086682v1` | Actinomycetota | `available` | `full` | 102 | 64 | 1 | 11 | 7 | 4 | 15 | sí | — |
| `GCF_032050975.1_ASM3205097v1` | unknown | `available` | `full` | 42 | 13 | 1 | 6 | 1 | 3 | 18 | sí | — |
| `GCF_032051675.1_ASM3205167v1` | unknown | `available` | `full` | 108 | 57 | 8 | 15 | 1 | 3 | 24 | sí | — |
| `GCF_032102965.1_ASM3210296v1` | Actinomycetota | `available` | `full` | 33 | 17 | 1 | 3 | 0 | 3 | 9 | sí | — |
| `GCF_032595795.2_ASM3259579v2` | Spirochaetota | `available` | `full` | 21 | 5 | 3 | 3 | 0 | 0 | 10 | sí | — |
| `GCF_032698475.1_ASM3269847v1` | Pseudomonadota | `available` | `full` | 217 | 94 | 10 | 6 | 0 | 5 | 102 | sí | — |
| `GCF_033547035.1_ASM3354703v1` | Myxococcota | `available` | `full` | 243 | 79 | 8 | 13 | 13 | 13 | 117 | sí | — |
| `GCF_033878955.1_ASM3387895v1` | Verrucomicrobiota | `available` | `full` | 155 | 87 | 2 | 17 | 21 | 1 | 27 | sí | — |
| `GCF_033969665.2_ASM3396966v2` | Spirochaetota | `available` | `full` | 26 | 7 | 4 | 6 | 0 | 0 | 9 | sí | — |
| `GCF_033971325.1_ASM3397132v1` | Aquificota | `available` | `full` | 48 | 18 | 0 | 3 | 16 | 0 | 11 | sí | — |
| `GCF_034259205.1_ASM3425920v1` | Bdellovibrionota | `available` | `full` | 57 | 18 | 4 | 2 | 1 | 6 | 26 | sí | — |
| `GCF_034424215.1_ASM3442421v1` | Bacteroidota | `available` | `full` | 95 | 38 | 3 | 9 | 8 | 5 | 32 | sí | — |
| `GCF_037076345.1_ASM3707634v1` | Verrucomicrobiota | `available` | `full` | 83 | 45 | 2 | 6 | 1 | 3 | 26 | sí | — |
| `GCF_037076455.1_ASM3707645v1` | Thermodesulfobacteriota | `available` | `full` | 182 | 160 | 7 | 3 | 0 | 2 | 10 | sí | — |
| `GCF_037081835.1_ASM3708183v1` | Bdellovibrionota | `available` | `full` | 117 | 81 | 2 | 3 | 10 | 4 | 17 | sí | — |
| `GCF_037157775.1_ASM3715777v1` | Myxococcota | `available` | `full` | 202 | 68 | 8 | 10 | 3 | 12 | 101 | sí | — |
| `GCF_037477945.1_MHOCP` | Bacillota | `available` | `full` | 73 | 47 | 0 | 7 | 8 | 1 | 10 | sí | — |
| `GCF_038095605.1_ASM3809560v1` | Gemmatimonadota | `available` | `full` | 63 | 20 | 1 | 4 | 3 | 5 | 30 | sí | — |
| `GCF_038396615.1_ASM3839661v1` | Atribacterota | `available` | `full` | 81 | 34 | 3 | 4 | 35 | 0 | 5 | sí | — |
| `GCF_038400315.1_ASM3840031v1` | Planctomycetota | `available` | `full` | 276 | 226 | 2 | 8 | 1 | 3 | 36 | sí | — |
| `GCF_038443385.1_ASM3844338v1` | Planctomycetota | `available` | `full` | 269 | 170 | 3 | 23 | 27 | 2 | 44 | sí | — |
| `GCF_038725785.1_ASM3872578v1` | Verrucomicrobiota | `available` | `full` | 289 | 182 | 4 | 7 | 8 | 3 | 85 | sí | — |
| `GCF_039555295.1_ASM3955529v1` | Campylobacterota | `available` | `full` | 57 | 9 | 3 | 11 | 4 | 4 | 26 | sí | — |
| `GCF_039645985.1_ASM3964598v1` | Campylobacterota | `available` | `full` | 35 | 13 | 0 | 0 | 2 | 2 | 18 | sí | — |
| `GCF_039994875.1_ASM3999487v1` | Acidobacteriota | `available` | `full` | 87 | 50 | 0 | 1 | 1 | 3 | 32 | sí | — |
| `GCF_040256395.1_GCC_DsonKR87_v1` | Deinococcota | `available` | `full` | 70 | 43 | 0 | 4 | 1 | 12 | 10 | sí | — |
| `GCF_040450955.1_ASM4045095v1` | Nitrospirota | `available` | `full` | 62 | 20 | 0 | 2 | 28 | 0 | 12 | sí | — |
| `GCF_040451425.1_ASM4045142v1` | Nitrospirota | `available` | `full` | 69 | 15 | 0 | 6 | 34 | 0 | 14 | sí | — |
| `GCF_040564635.1_ASM4056463v1` | Cyanobacteriota | `available` | `full` | 270 | 120 | 1 | 26 | 85 | 4 | 34 | sí | — |
| `GCF_040833985.2_ASM4083398v1` | Spirochaetota | `available` | `full` | 48 | 18 | 0 | 7 | 0 | 1 | 22 | sí | — |
| `GCF_041026625.1_ASM4102662v1` | Mycoplasmatota | `available` | `full` | 14 | 5 | 0 | 3 | 3 | 0 | 3 | sí | — |
| `GCF_041154365.1_ASM4115436v1` | Fidelibacterota | `available` | `full` | 51 | 23 | 2 | 5 | 12 | 1 | 8 | sí | — |
| `GCF_041200005.1_ASM4120000v1` | Fusobacteriota | `available` | `full` | 61 | 24 | 2 | 11 | 17 | 0 | 7 | sí | — |
| `GCF_041505935.1_Fervidobacterium_thailandense_FC2004` | Thermotogota | `available` | `full` | 114 | 70 | 1 | 2 | 32 | 0 | 9 | sí | — |
| `GCF_042033325.1_ASM4203332v1` | Fusobacteriota | `available` | `full` | 25 | 14 | 0 | 3 | 0 | 0 | 8 | sí | — |
| `GCF_042847215.1_ASM4284721v1` | Pseudomonadota | `available` | `full` | 111 | 38 | 2 | 7 | 7 | 5 | 52 | sí | — |
| `GCF_044336245.1_ASM4433624v1` | Planctomycetota | `available` | `full` | 182 | 133 | 1 | 3 | 13 | 5 | 27 | sí | — |
| `GCF_044789115.1_ASM4478911v1` | Verrucomicrobiota | `available` | `full` | 45 | 26 | 2 | 5 | 0 | 1 | 11 | sí | — |
| `GCF_045784255.1_ASM4578425v1` | Deinococcota | `available` | `full` | 273 | 239 | 7 | 7 | 1 | 6 | 13 | sí | — |
| `GCF_045784305.1_DRML_24019` | Deinococcota | `available` | `full` | 510 | 472 | 6 | 8 | 3 | 6 | 15 | sí | — |
| `GCF_045864955.1_ASM4586495v1` | Atribacterota | `available` | `full` | 70 | 29 | 2 | 4 | 28 | 0 | 7 | sí | — |
| `GCF_046267995.1_ASM4626799v1` | Bacillota | `available` | `full` | 102 | 71 | 4 | 3 | 4 | 9 | 11 | sí | — |
| `GCF_046774955.1_ASM4677495v1` | Deinococcota | `available` | `full` | 114 | 74 | 3 | 6 | 0 | 6 | 25 | sí | — |
| `GCF_046845665.1_ASM4684566v1` | Rhodothermota | `available` | `full` | 48 | 15 | 2 | 4 | 0 | 1 | 26 | sí | — |
| `GCF_046845755.1_ASM4684575v1` | Rhodothermota | `available` | `full` | 46 | 12 | 2 | 6 | 0 | 2 | 24 | sí | — |
| `GCF_047159735.1_ASM4715973v1` | Minisyncoccota | `available` | `full` | 100 | 60 | 2 | 6 | 1 | 0 | 31 | sí | — |
| `GCF_048565275.1_ASM4856527v1` | Thermodesulfobacteriota | `available` | `full` | 307 | 174 | 11 | 52 | 18 | 4 | 48 | sí | — |
| `GCF_049669535.1_sGff_cult_auto_v1` | Mycoplasmatota | `available` | `full` | 233 | 218 | 7 | 1 | 0 | 0 | 7 | sí | — |
| `GCF_049676565.1_Dip510` | Elusimicrobiota | `available` | `full` | 142 | 7 | 4 | 0 | 0 | 1 | 130 | sí | — |
| `GCF_050039385.1_ASM5003938v1` | Myxococcota | `available` | `full` | 243 | 92 | 7 | 7 | 42 | 9 | 86 | sí | — |
| `GCF_050267035.1_ASM5026703v1` | unknown | `available` | `full` | 102 | 43 | 7 | 1 | 44 | 0 | 7 | sí | — |
| `GCF_050389785.1_ASM5038978v1` | Elusimicrobiota | `available` | `full` | 125 | 9 | 0 | 2 | 0 | 0 | 114 | sí | — |
| `GCF_050500235.1_ASM5050023v1` | Cyanobacteriota | `available` | `full` | 388 | 299 | 1 | 37 | 30 | 0 | 21 | sí | — |
| `GCF_050563455.1_ASM5056345v1` | Cyanobacteriota | `available` | `full` | 1.085 | 1.021 | 0 | 37 | 6 | 0 | 21 | sí | — |
| `GCF_050780155.1_QI0013` | Thermodesulfobacteriota | `available` | `full` | 119 | 37 | 14 | 3 | 6 | 3 | 56 | sí | — |
| `GCF_051175385.1_ASM5117538v1` | Planctomycetota | `available` | `full` | 198 | 162 | 4 | 6 | 3 | 2 | 21 | sí | — |
| `GCF_052161605.1_ASM5216160v1` | Planctomycetota | `available` | `full` | 188 | 122 | 6 | 17 | 9 | 1 | 33 | sí | — |
| `GCF_054054625.1_ASM5405462v1` | Bdellovibrionota | `available` | `full` | 68 | 33 | 1 | 2 | 0 | 5 | 27 | sí | — |
| `GCF_054917915.1_ASM5491791v1` | Pseudomonadota | `available` | `full` | 89 | 40 | 1 | 7 | 0 | 6 | 35 | sí | — |
| `GCF_055378285.1_ASM5537828v1` | Chrysiogenota | `available` | `full` | 91 | 57 | 1 | 8 | 9 | 2 | 14 | sí | — |
| `GCF_055378975.1_ASM5537897v1` | Bacillota | `available` | `full` | 132 | 96 | 4 | 7 | 2 | 4 | 19 | sí | — |
| `GCF_055379415.1_ASM5537941v1` | Gemmatimonadota | `available` | `full` | 67 | 13 | 2 | 1 | 1 | 3 | 47 | sí | — |
| `GCF_055379825.1_ASM5537982v1` | Rhodothermota | `available` | `full` | 135 | 85 | 0 | 19 | 25 | 1 | 5 | sí | — |
| `GCF_055379935.1_ASM5537993v1` | Rhodothermota | `available` | `full` | 92 | 48 | 0 | 21 | 17 | 1 | 5 | sí | — |
| `GCF_055380895.1_ASM5538089v1` | Aquificota | `available` | `full` | 35 | 10 | 1 | 1 | 5 | 0 | 18 | sí | — |
| `GCF_055383505.1_ASM5538350v1` | Fusobacteriota | `available` | `full` | 46 | 15 | 2 | 4 | 5 | 1 | 19 | sí | — |
| `GCF_055384375.1_ASM5538437v1` | Campylobacterota | `available` | `full` | 46 | 20 | 0 | 6 | 7 | 0 | 13 | sí | — |
| `GCF_055384745.1_ASM5538474v1` | Aquificota | `available` | `full` | 37 | 8 | 0 | 0 | 12 | 0 | 17 | sí | — |
| `GCF_055385515.1_ASM5538551v1` | Aquificota | `available` | `full` | 53 | 27 | 3 | 5 | 1 | 0 | 17 | sí | — |
| `GCF_055387075.1_ASM5538707v1` | Thermotogota | `available` | `full` | 68 | 47 | 0 | 5 | 9 | 0 | 7 | sí | — |
| `GCF_055387105.1_ASM5538710v1` | Caldifarciminota | `available` | `full` | 78 | 39 | 4 | 6 | 15 | 0 | 14 | sí | — |
| `GCF_055387525.1_ASM5538752v1` | Chloroflexota | `available` | `full` | 113 | 38 | 5 | 8 | 47 | 3 | 12 | sí | — |
| `GCF_055388355.1_ASM5538835v1` | Verrucomicrobiota | `available` | `full` | 73 | 39 | 0 | 6 | 2 | 5 | 21 | sí | — |
| `GCF_055390195.1_ASM5539019v1` | Aquificota | `available` | `full` | 52 | 17 | 0 | 2 | 18 | 1 | 14 | sí | — |
| `GCF_055390205.1_ASM5539020v1` | Thermotogota | `available` | `full` | 95 | 48 | 1 | 8 | 24 | 2 | 12 | sí | — |
| `GCF_055390565.1_ASM5539056v1` | Synergistota | `available` | `full` | 39 | 11 | 2 | 1 | 13 | 1 | 11 | sí | — |
| `GCF_055390845.1_ASM5539084v1` | Nitrospirota | `available` | `full` | 89 | 39 | 1 | 9 | 24 | 1 | 15 | sí | — |
| `GCF_055390895.1_ASM5539089v1` | Synergistota | `available` | `full` | 104 | 41 | 2 | 8 | 35 | 0 | 18 | sí | — |
| `GCF_055391325.1_ASM5539132v1` | Thermotogota | `available` | `full` | 62 | 17 | 1 | 2 | 30 | 0 | 12 | sí | — |
| `GCF_055396605.1_ASM5539660v1` | Bacillota | `available` | `full` | 82 | 40 | 0 | 2 | 20 | 3 | 17 | sí | — |
| `GCF_055694045.1_ASM5569404v1` | Fusobacteriota | `available` | `full` | 66 | 51 | 0 | 6 | 2 | 0 | 7 | sí | — |
| `GCF_055763695.1_ASM5576369v1` | Actinomycetota | `available` | `full` | 138 | 61 | 2 | 14 | 21 | 12 | 28 | sí | — |
| `GCF_056388415.1_ASM5638841v1` | Mycoplasmatota | `available` | `full` | 19 | 10 | 0 | 4 | 0 | 0 | 5 | sí | — |
| `GCF_058362665.1_ASM5836266v1` | Myxococcota | `available` | `full` | 120 | 29 | 3 | 9 | 10 | 6 | 63 | sí | — |
| `GCF_059632245.1_ASM5963224v1` | Planctomycetota | `available` | `full` | 102 | 55 | 5 | 7 | 5 | 0 | 30 | sí | — |
| `GCF_059696725.1_ASM5969672v1` | Chloroflexota | `available` | `full` | 77 | 14 | 2 | 5 | 44 | 2 | 10 | sí | — |
| `GCF_059705355.1_ASM5970535v1` | Campylobacterota | `available` | `full` | 72 | 53 | 1 | 7 | 0 | 1 | 10 | sí | — |
| `GCF_059726255.1_ASM5972625v1` | Acidobacteriota | `available` | `full` | 64 | 22 | 0 | 1 | 3 | 4 | 34 | sí | — |
| `GCF_059726345.1_ASM5972634v1` | Acidobacteriota | `available` | `full` | 228 | 127 | 3 | 13 | 16 | 7 | 62 | sí | — |
| `GCF_060026435.1_ASM6002643v1` | Myxococcota | `available` | `full` | 271 | 108 | 6 | 14 | 22 | 11 | 110 | sí | — |
| `GCF_060027055.1_ASM6002705v1` | Myxococcota | `available` | `full` | 266 | 77 | 5 | 16 | 41 | 3 | 124 | sí | — |
| `GCF_060030665.1_ASM6003066v1` | Myxococcota | `available` | `full` | 293 | 144 | 7 | 13 | 15 | 7 | 107 | sí | — |
| `GCF_060030705.1_ASM6003070v1` | Myxococcota | `available` | `full` | 136 | 87 | 1 | 7 | 0 | 4 | 37 | sí | — |
| `GCF_060030825.1_ASM6003082v1` | Myxococcota | `available` | `full` | 193 | 63 | 2 | 7 | 15 | 8 | 98 | sí | — |
| `GCF_060031025.1_ASM6003102v1` | Myxococcota | `available` | `full` | 333 | 60 | 2 | 13 | 63 | 6 | 189 | sí | — |
| `GCF_060138195.1_ASM6013819v1` | Myxococcota | `available` | `full` | 203 | 52 | 3 | 7 | 58 | 8 | 75 | sí | — |
| `GCF_900169565.1_NSJP_Ch1` | Nitrospirota | `available` | `full` | 66 | 25 | 1 | 1 | 1 | 7 | 31 | sí | — |
| `GCF_900184705.1_CAMBI-1` | unknown | `available` | `full` | 48 | 16 | 0 | 4 | 21 | 2 | 5 | sí | — |
| `GCF_900217235.1_DPRO_PRJEB22548_v1` | Thermodesulfobacteriota | `available` | `full` | 385 | 341 | 5 | 19 | 0 | 1 | 19 | sí | — |
| `GCF_900239975.1_Chlamydia_sp._S15-834K` | unknown | `available` | `full` | 50 | 6 | 0 | 1 | 0 | 0 | 43 | sí | — |
| `GCF_900475595.1_44343_G01` | Bacillota | `available` | `full` | 44 | 12 | 2 | 1 | 4 | 3 | 22 | sí | — |
| `GCF_900660485.1_50648_E01-3` | Mycoplasmatota | `available` | `full` | 32 | 16 | 0 | 10 | 0 | 0 | 6 | sí | — |
| `GCF_900660665.1_57612_D01-3` | Mycoplasmatota | `available` | `full` | 35 | 21 | 1 | 4 | 3 | 0 | 6 | sí | — |
| `GCF_901538265.1_ASM90153826v1` | Planctomycetota | `available` | `full` | 326 | 244 | 12 | 10 | 14 | 2 | 44 | sí | — |
| `GCF_902806995.1_15-2067_O50` | Chlamydiota | `available` | `full` | 45 | 5 | 0 | 1 | 0 | 0 | 39 | sí | — |
| `GCF_907165195.1_Cryomorphaceae_CECT30217` | Bacteroidota | `available` | `full` | 63 | 38 | 0 | 7 | 6 | 2 | 10 | sí | — |
| `GCF_918320885.1_Harlequin_endosymbiont1` | unknown | `available` | `full` | 258 | 203 | 4 | 5 | 0 | 2 | 44 | sí | — |
| `GCF_946900835.1_NSPWAT` | Nitrospinota | `available` | `full` | 39 | 8 | 1 | 4 | 10 | 1 | 15 | sí | — |
| `GCF_947241125.1_NsDNF` | Nitrospirota | `available` | `full` | 227 | 126 | 4 | 24 | 39 | 5 | 29 | sí | — |
| `GCF_949774925.1_MFUM-new` | unknown | `available` | `full` | 102 | 74 | 1 | 3 | 19 | 1 | 4 | sí | — |
| `GCF_963264975.1_NZ_220527_assembly_043` | Fusobacteriota | `available` | `full` | 116 | 79 | 3 | 5 | 20 | 0 | 9 | sí | — |
| `GCF_964030865.1_idDioLine1.Spiroplasma_sp_1.1` | Mycoplasmatota | `available` | `full` | 19 | 19 | 0 | 0 | 0 | 0 | 0 | sí | — |
| `GCF_964030885.1_icBarPell3.Wolbachia_sp_1.1` | Pseudomonadota | `available` | `full` | 337 | 321 | 5 | 3 | 0 | 3 | 5 | sí | — |
| `GCF_964642095.1_Ff3` | Synergistota | `available` | `full` | 136 | 76 | 2 | 14 | 28 | 0 | 16 | sí | — |
| `GCF_965152565.1_NS_bockiana` | unknown | `available` | `full` | 114 | 67 | 5 | 7 | 15 | 2 | 18 | sí | — |
| `GCF_976986435.1_NS_calida2` | Nitrospirota | `available` | `full` | 108 | 52 | 1 | 7 | 31 | 2 | 15 | sí | — |
| `GCF_978936825.1_DTS69` | Thermodesulfobacteriota | `available` | `full` | 92 | 49 | 0 | 8 | 19 | 2 | 14 | sí | — |

</details>

## Advertencias

- GCF_000006685.1_ASM668v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000007205.1_ASM720v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000007605.1_ASM760v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000008625.1_ASM862v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000008725.1_ASM872v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000011905.1_ASM1190v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000018945.1_ASM1894v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000020145.1_ASM2014v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000020945.1_ASM2094v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000020965.1_ASM2096v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000021545.1_ASM2154v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000021565.1_ASM2156v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000021645.1_ASM2164v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000024905.1_ASM2490v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000025605.1_ASM2560v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000025885.1_ASM2588v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000164905.1_ASM16490v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000185805.1_ASM18580v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000191045.1_ASM19104v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000204255.1_ASM20425v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000212395.1_ASM21239v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000215975.1_ASM21597v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000239475.1_ASM23947v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000259275.1_ASM25927v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000284335.1_ASM28433v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000392435.1_ASM39243v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000439435.1_ASM43943v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000471025.2_ASM47102v2: inventario dice has_fasta=no pero el archivo existe
- GCF_000479335.1_ASM47933v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000508245.1_ASM50824v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000512735.1_ASM51273v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000513475.1_Bbla: inventario dice has_fasta=no pero el archivo existe
- GCF_000582535.1_ASM58253v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000803625.1_ASM80362v1: inventario dice has_fasta=no pero el archivo existe
- GCF_000941075.1_ASM94107v1: inventario dice has_fasta=no pero el archivo existe
- GCF_001027545.1_ASM102754v1: inventario dice has_fasta=no pero el archivo existe
- GCF_001547735.1_ASM154773v1: inventario dice has_fasta=no pero el archivo existe
- GCF_001712875.1_ASM171287v1: inventario dice has_fasta=no pero el archivo existe
- GCF_002139875.1_ASM213987v1: inventario dice has_fasta=no pero el archivo existe
- GCF_002355835.1_ASM235583v1: inventario dice has_fasta=no pero el archivo existe
- GCF_002804025.1_ASM280402v1: inventario dice has_fasta=no pero el archivo existe
- GCF_003057965.1_ASM305796v1: inventario dice has_fasta=no pero el archivo existe
- GCF_003351905.1_ASM335190v1: inventario dice has_fasta=no pero el archivo existe
- GCF_003355515.1_ASM335551v1: inventario dice has_fasta=no pero el archivo existe
- GCF_003366055.1_ASM336605v1: inventario dice has_fasta=no pero el archivo existe
- GCF_004102945.1_ASM410294v1: inventario dice has_fasta=no pero el archivo existe
- GCF_005697565.1_ASM569756v1: inventario dice has_fasta=no pero el archivo existe
- GCF_008086545.1_ASM808654v1: inventario dice has_fasta=no pero el archivo existe
- GCF_013177655.1_ASM1317765v1: inventario dice has_fasta=no pero el archivo existe
- GCF_015356815.2_ASM1535681v2: inventario dice has_fasta=no pero el archivo existe
- GCF_018343815.1_ASM1834381v1: inventario dice has_fasta=no pero el archivo existe
- GCF_018866205.1_ASM1886620v1: inventario dice has_fasta=no pero el archivo existe
- GCF_019056495.1_ASM1905649v1: inventario dice has_fasta=no pero el archivo existe
- GCF_020549965.1_ASM2054996v1: inventario dice has_fasta=no pero el archivo existe
- GCF_021057185.1_ASM2105718v1: inventario dice has_fasta=no pero el archivo existe
- GCF_023035795.1_ASM2303579v1: inventario dice has_fasta=no pero el archivo existe
- GCF_023035875.1_ASM2303587v1: inventario dice has_fasta=no pero el archivo existe
- GCF_024338845.1_ASM2433884v1: inventario dice has_fasta=no pero el archivo existe
- GCF_025273655.1_ASM2527365v1: inventario dice has_fasta=no pero el archivo existe
- GCF_026184435.1_ASM2618443v1: inventario dice has_fasta=no pero el archivo existe
- GCF_029677705.1_ASM2967770v1: inventario dice has_fasta=no pero el archivo existe
- GCF_030253535.1_ASM3025353v1: inventario dice has_fasta=no pero el archivo existe
- GCF_030438115.2_ASM3043811v2: inventario dice has_fasta=no pero el archivo existe
- GCF_030439285.2_ASM3043928v2: inventario dice has_fasta=no pero el archivo existe
- GCF_030585795.1_ASM3058579v1: inventario dice has_fasta=no pero el archivo existe
- GCF_032595795.2_ASM3259579v2: inventario dice has_fasta=no pero el archivo existe
- GCF_033969665.2_ASM3396966v2: inventario dice has_fasta=no pero el archivo existe
- GCF_033971325.1_ASM3397132v1: inventario dice has_fasta=no pero el archivo existe
- GCF_040450955.1_ASM4045095v1: inventario dice has_fasta=no pero el archivo existe
- GCF_040451425.1_ASM4045142v1: inventario dice has_fasta=no pero el archivo existe
- GCF_041026625.1_ASM4102662v1: inventario dice has_fasta=no pero el archivo existe
- GCF_042033325.1_ASM4203332v1: inventario dice has_fasta=no pero el archivo existe
- GCF_046267995.1_ASM4626799v1: inventario dice has_fasta=no pero el archivo existe
- GCF_047159735.1_ASM4715973v1: inventario dice has_fasta=no pero el archivo existe
- GCF_049669535.1_sGff_cult_auto_v1: inventario dice has_fasta=no pero el archivo existe
- GCF_050389785.1_ASM5038978v1: inventario dice has_fasta=no pero el archivo existe
- GCF_055380895.1_ASM5538089v1: inventario dice has_fasta=no pero el archivo existe
- GCF_055384745.1_ASM5538474v1: inventario dice has_fasta=no pero el archivo existe
- GCF_055385515.1_ASM5538551v1: inventario dice has_fasta=no pero el archivo existe
- GCF_055390195.1_ASM5539019v1: inventario dice has_fasta=no pero el archivo existe
- GCF_055390845.1_ASM5539084v1: inventario dice has_fasta=no pero el archivo existe
- GCF_056388415.1_ASM5638841v1: inventario dice has_fasta=no pero el archivo existe
- GCF_900239975.1_Chlamydia_sp._S15-834K: inventario dice has_fasta=no pero el archivo existe
- GCF_900475595.1_44343_G01: inventario dice has_fasta=no pero el archivo existe
- GCF_900660485.1_50648_E01-3: inventario dice has_fasta=no pero el archivo existe
- GCF_900660665.1_57612_D01-3: inventario dice has_fasta=no pero el archivo existe
- GCF_902806995.1_15-2067_O50: inventario dice has_fasta=no pero el archivo existe
- GCF_918320885.1_Harlequin_endosymbiont1: inventario dice has_fasta=no pero el archivo existe
- GCF_964030865.1_idDioLine1.Spiroplasma_sp_1.1: inventario dice has_fasta=no pero el archivo existe
- GCF_964030885.1_icBarPell3.Wolbachia_sp_1.1: inventario dice has_fasta=no pero el archivo existe

## Método y límites

- Solo lectura: el informe se construye leyendo TSV/GFF/JSON/FASTA existentes.
- `not_available` significa que **no hay dato en disco**, no que la biología esté ausente. En particular, CRISPR sin FASTA es “no evaluado”.
- La atribución de cada artefacto a un script es estática (`22_run_full.sh`, `paths.py`, nombres de salida); no se ejecutó nada.
- Los recuentos de la capa previa y del TSV tipado no son equivalentes: la primera es keyword-based y sin coordenadas por contig.
- Para reproducir: `python tools/mobile_coverage.py`.
