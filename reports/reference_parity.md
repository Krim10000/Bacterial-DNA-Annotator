# Paridad de anotación: sistema vs GFF de referencia

**Generado:** `2026-11-25T16:39:00Z`  
**Manifiesto:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/reference_genomes.json`  
**Raíz de datos:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`

## Alcance

Compara los conteos de **CDS, tRNA y rRNA** del resumen del sistema (`results/<id>/<id>.summary.json`) contra los features del GFF de referencia (`annotation_path`) para cada genoma con `status` `annotated` o `available` y `system_summary` en el manifiesto.

Esta versión **no evalúa función, productos, coordenadas ni identidad de secuencia**. Tampoco modifica `results/`, los GFF ni el manifiesto.

## Cobertura

- Registros con status `annotated`/`available`: **20**
- Evaluados (con `system_summary`): **18** (ok **0**, partial **18**, not_available **0**)
- Omitidos: **24** (2 elegibles sin `system_summary`; 22 con otro status)

### Omitidos con status elegible

| genoma | status | motivo |
|---|---|---|
| `ecoli_k12` | available | sin system_summary en el manifiesto |
| `subtilis` | available | sin system_summary en el manifiesto |

## Resultados por genoma

| genoma | estado | CDS sis | CDS ref | Δ CDS | % CDS | tRNA sis | tRNA ref | Δ tRNA | % tRNA | rRNA sis | rRNA ref | Δ rRNA | % rRNA |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `mtb_h37rv` | partial | 3986 | 3906 | +80 | 102.05 | 45 | 45 | 0 | 100.00 | 2 | 3 | -1 | 66.67 |
| `cglutamicum` | partial | 3041 | 2999 | +42 | 101.40 | 60 | 59 | +1 | 101.69 | 12 | 18 | -6 | 66.67 |
| `pputida` | partial | 5507 | 5471 | +36 | 100.66 | 75 | 75 | 0 | 100.00 | 15 | 22 | -7 | 68.18 |
| `paeruginosa` | partial | 5715 | 5572 | +143 | 102.57 | 64 | 63 | +1 | 101.59 | 8 | 13 | -5 | 61.54 |
| `crescentus` | partial | 3738 | 3886 | -148 | 96.19 | 52 | 51 | +1 | 101.96 | 4 | 6 | -2 | 66.67 |
| `synechocystis` | partial | 3780 | 3698 | +82 | 102.22 | 41 | 41 | 0 | 100.00 | 2 | 6 | -4 | 33.33 |
| `saureus` | partial | 2675 | 2767 | -92 | 96.68 | 61 | 59 | +2 | 103.39 | 11 | 16 | -5 | 68.75 |
| `spneumoniae` | partial | 2168 | 2117 | +51 | 102.41 | 59 | 58 | +1 | 101.72 | 8 | 12 | -4 | 66.67 |
| `btheta` | partial | 4945 | 4801 | +144 | 103.00 | 71 | 70 | +1 | 101.43 | 10 | 15 | -5 | 66.67 |
| `scoelicolor` | partial | 7738 | 7755 | -17 | 99.78 | 72 | 66 | +6 | 109.09 | 12 | 18 | -6 | 66.67 |
| `mxanthus` | partial | 7325 | 7306 | +19 | 100.26 | 69 | 66 | +3 | 104.55 | 8 | 12 | -4 | 66.67 |
| `vcholerae` | partial | 3666 | 3665 | +1 | 100.03 | 102 | 102 | 0 | 100.00 | 21 | 31 | -10 | 67.74 |
| `cdifficile` | partial | 3700 | 3611 | +89 | 102.46 | 91 | 90 | +1 | 101.11 | 23 | 35 | -12 | 65.71 |
| `llactis` | partial | 2602 | 2567 | +35 | 101.36 | 71 | 70 | +1 | 101.43 | 13 | 19 | -6 | 68.42 |
| `lmonocytogenes` | partial | 2921 | 2867 | +54 | 101.88 | 67 | 67 | 0 | 100.00 | 12 | 18 | -6 | 66.67 |
| `efaecalis` | partial | 2730 | 2690 | +40 | 101.49 | 53 | 51 | +2 | 103.92 | 8 | 12 | -4 | 66.67 |
| `dradiodurans` | partial | 3150 | 3148 | +2 | 100.06 | 51 | 50 | +1 | 102.00 | 3 | 9 | -6 | 33.33 |
| `tthermophilus` | partial | 2227 | 2228 | -1 | 99.96 | 49 | 48 | +1 | 102.08 | 4 | 6 | -2 | 66.67 |

## Estado y motivos

| genoma | estado | motivo |
|---|---|---|
| `mtb_h37rv` | partial | delta CDS=+80, rRNA=-1 |
| `cglutamicum` | partial | delta CDS=+42, tRNA=+1, rRNA=-6 |
| `pputida` | partial | delta CDS=+36, rRNA=-7 |
| `paeruginosa` | partial | delta CDS=+143, tRNA=+1, rRNA=-5 |
| `crescentus` | partial | delta CDS=-148, tRNA=+1, rRNA=-2 |
| `synechocystis` | partial | delta CDS=+82, rRNA=-4 |
| `saureus` | partial | delta CDS=-92, tRNA=+2, rRNA=-5 |
| `spneumoniae` | partial | delta CDS=+51, tRNA=+1, rRNA=-4 |
| `btheta` | partial | delta CDS=+144, tRNA=+1, rRNA=-5 |
| `scoelicolor` | partial | delta CDS=-17, tRNA=+6, rRNA=-6 |
| `mxanthus` | partial | delta CDS=+19, tRNA=+3, rRNA=-4 |
| `vcholerae` | partial | delta CDS=+1, rRNA=-10 |
| `cdifficile` | partial | delta CDS=+89, tRNA=+1, rRNA=-12 |
| `llactis` | partial | delta CDS=+35, tRNA=+1, rRNA=-6 |
| `lmonocytogenes` | partial | delta CDS=+54, rRNA=-6 |
| `efaecalis` | partial | delta CDS=+40, tRNA=+2, rRNA=-4 |
| `dradiodurans` | partial | delta CDS=+2, tRNA=+1, rRNA=-6 |
| `tthermophilus` | partial | delta CDS=-1, tRNA=+1, rRNA=-2 |

## Entradas y metadatos

| genoma | especie | cepa | status manifiesto | GFF de referencia | resumen del sistema |
|---|---|---|---|---|---|
| `mtb_h37rv` | Mycobacterium tuberculosis | H37Rv | annotated | `data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.gff.gz` | `results/mtb_h37rv/mtb_h37rv.summary.json` |
| `cglutamicum` | Corynebacterium glutamicum | ATCC 13032 | annotated | `data/reference_genomes/cglutamicum/cglutamicum_genomic.gff.gz` | `results/cglutamicum/cglutamicum.summary.json` |
| `pputida` | Pseudomonas putida | KT2440 | annotated | `data/reference_genomes/pputida/pputida_genomic.gff.gz` | `results/pputida/pputida.summary.json` |
| `paeruginosa` | Pseudomonas aeruginosa | PAO1 | annotated | `data/reference_genomes/paeruginosa/paeruginosa_genomic.gff.gz` | `results/paeruginosa/paeruginosa.summary.json` |
| `crescentus` | Caulobacter crescentus | NA1000 | annotated | `data/reference_genomes/crescentus/crescentus_genomic.gff.gz` | `results/crescentus/crescentus.summary.json` |
| `synechocystis` | Synechocystis sp. | PCC 6803 | annotated | `data/reference_genomes/synechocystis/synechocystis_genomic.gff.gz` | `results/synechocystis/synechocystis.summary.json` |
| `saureus` | Staphylococcus aureus | NCTC 8325 | annotated | `data/reference_genomes/saureus/saureus_genomic.gff.gz` | `results/saureus/saureus.summary.json` |
| `spneumoniae` | Streptococcus pneumoniae | D39 | annotated | `data/reference_genomes/spneumoniae/spneumoniae_genomic.gff.gz` | `results/spneumoniae/spneumoniae.summary.json` |
| `btheta` | Bacteroides thetaiotaomicron | VPI-5482 | annotated | `data/reference_genomes/btheta/btheta_genomic.gff.gz` | `results/btheta/btheta.summary.json` |
| `scoelicolor` | Streptomyces coelicolor | A3(2) | annotated | `data/reference_genomes/scoelicolor/scoelicolor_genomic.gff.gz` | `results/scoelicolor/scoelicolor.summary.json` |
| `mxanthus` | Myxococcus xanthus | DK1622 | annotated | `data/reference_genomes/mxanthus/mxanthus_genomic.gff.gz` | `results/mxanthus/mxanthus.summary.json` |
| `vcholerae` | Vibrio cholerae | N16961 | annotated | `data/reference_genomes/vcholerae/vcholerae_genomic.gff.gz` | `results/vcholerae/vcholerae.summary.json` |
| `cdifficile` | Clostridioides difficile | 630 | annotated | `data/reference_genomes/cdifficile/cdifficile_genomic.gff.gz` | `results/cdifficile/cdifficile.summary.json` |
| `llactis` | Lactococcus lactis | IL1403 | annotated | `data/reference_genomes/llactis/llactis_genomic.gff.gz` | `results/llactis/llactis.summary.json` |
| `lmonocytogenes` | Listeria monocytogenes | EGD-e | annotated | `data/reference_genomes/lmonocytogenes/lmonocytogenes_genomic.gff.gz` | `results/lmonocytogenes/lmonocytogenes.summary.json` |
| `efaecalis` | Enterococcus faecalis | V583 | annotated | `data/reference_genomes/efaecalis/efaecalis_genomic.gff.gz` | `results/efaecalis/efaecalis.summary.json` |
| `dradiodurans` | Deinococcus radiodurans | R1 | annotated | `data/reference_genomes/dradiodurans/dradiodurans_genomic.gff.gz` | `results/dradiodurans/dradiodurans.summary.json` |
| `tthermophilus` | Thermus thermophilus | HB27 | annotated | `data/reference_genomes/tthermophilus/tthermophilus_genomic.gff.gz` | `results/tthermophilus/tthermophilus.summary.json` |

## Agregado (micro)

| clase | sistema | referencia | delta | % sistema/referencia |
|---|---:|---:|---:|---:|
| CDS | 71614 | 71054 | +560 | 100.79 |
| tRNA | 1153 | 1131 | +22 | 101.95 |
| rRNA | 176 | 271 | -95 | 64.94 |

El agregado suma solo genomas con conteo de sistema y referencia disponibles para esa clase.

## Criterio de estado

- **ok**: los tres conteos (CDS, tRNA y rRNA) coinciden exactamente con la referencia y la referencia tiene las tres clases.
- **partial**: la comparación fue posible pero al menos un conteo difiere, o la referencia carece de alguna de las tres clases.
- **not_available**: falta o es ilegible el resumen, el GFF, o el GFF no contiene ninguna de las tres clases.

En los recuentos de referencia, un CDS multiparte (varias filas con el mismo `ID`) cuenta como un solo feature. No se filtran pseudogenes por atributo: solo se usa la columna tipo del GFF.

## Límites

- Solo mide presencia y número de features; no mide exactitud de coordenadas, función, productos ni identidad de secuencia.
- El gold hereda los criterios y errores de la fuente RefSeq/NCBI; no es una verdad absoluta.
- Las diferencias de rRNA pueden reflejar que el pipeline no llama todas las copias o tipos (p. ej. 23S/5S); este informe no lo atribuye a una causa sin evidencia adicional.
- No se comparan ncRNA ni otros tipos, aunque estén presentes en el GFF de referencia.
- Los genomas sin `system_summary` en el manifiesto quedan fuera del cálculo; se listan arriba cuando su status era elegible.

## Reproducir

```bash
.venv/bin/python scripts/reference/reference_parity.py
```
