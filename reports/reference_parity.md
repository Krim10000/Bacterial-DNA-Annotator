# Paridad de anotación: sistema vs GFF de referencia

**Generado:** `2026-11-25T15:53:30Z`  
**Manifiesto:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/reference_genomes.json`  
**Raíz de datos:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`

## Alcance

Compara los conteos de **CDS, tRNA y rRNA** del resumen del sistema (`results/<id>/<id>.summary.json`) contra los features del GFF de referencia (`annotation_path`) para cada genoma con `status` `annotated` o `available` y `system_summary` en el manifiesto.

Esta versión **no evalúa función, productos, coordenadas ni identidad de secuencia**. Tampoco modifica `results/`, los GFF ni el manifiesto.

## Cobertura

- Registros con status `annotated`/`available`: **9**
- Evaluados (con `system_summary`): **7** (ok **0**, partial **7**, not_available **0**)
- Omitidos: **35** (2 elegibles sin `system_summary`; 33 con otro status)

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

## Agregado (micro)

| clase | sistema | referencia | delta | % sistema/referencia |
|---|---:|---:|---:|---:|
| CDS | 28442 | 28299 | +143 | 100.51 |
| tRNA | 398 | 393 | +5 | 101.27 |
| rRNA | 54 | 84 | -30 | 64.29 |

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
