# Paridad de anotación: sistema vs GFF de referencia

**Generado:** `2026-11-25T19:46:56Z`  
**Manifiesto:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/reference_genomes.json`  
**Raíz de datos:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`

## Alcance

Compara los conteos de **CDS, tRNA y rRNA** del resumen del sistema (`results/<id>/<id>.summary.json`) contra los features del GFF de referencia (`annotation_path`) para cada genoma con `status` `annotated` o `available` y `system_summary` en el manifiesto.

Esta versión **no evalúa función, productos, coordenadas ni identidad de secuencia**. Tampoco modifica `results/`, los GFF ni el manifiesto.

## Cobertura

- Registros con status `annotated`/`available`: **41**
- Evaluados (con `system_summary`): **39** (ok **0**, partial **39**, not_available **0**)
- Omitidos: **3** (2 elegibles sin `system_summary`; 1 con otro status)

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
| `hpylori` | partial | 1580 | 1584 | -4 | 99.75 | 37 | 36 | +1 | 102.78 | 2 | 6 | -4 | 33.33 |
| `cjejuni` | partial | 1673 | 1572 | +101 | 106.42 | 45 | 43 | +2 | 104.65 | 3 | 9 | -6 | 33.33 |
| `nmeningitidis` | partial | 2050 | 2142 | -92 | 95.70 | 61 | 59 | +2 | 103.39 | 8 | 12 | -4 | 66.67 |
| `hinfluenzae` | partial | 1819 | 1806 | +13 | 100.72 | 61 | 59 | +2 | 103.39 | 13 | 19 | -6 | 68.42 |
| `ctrachomatis` | partial | 899 | 887 | +12 | 101.35 | 37 | 37 | 0 | 100.00 | 4 | 6 | -2 | 66.67 |
| `tpallidum` | partial | 982 | 987 | -5 | 99.49 | 45 | 45 | 0 | 100.00 | 4 | 6 | -2 | 66.67 |
| `bburgdorferi` | partial | 1534 | 1528 | +6 | 100.39 | 32 | 32 | 0 | 100.00 | 3 | 5 | -2 | 60.00 |
| `mgenitalium` | partial | 494 | 524 | -30 | 94.27 | 36 | 36 | 0 | 100.00 | 2 | 3 | -1 | 66.67 |
| `abaumannii` | partial | 3801 | 3736 | +65 | 101.74 | 77 | 74 | +3 | 104.05 | 12 | 18 | -6 | 66.67 |
| `kpneumoniae` | partial | 5532 | 5779 | -247 | 95.73 | 91 | 62 | +29 | 146.77 | 17 | 25 | -8 | 68.00 |
| `styphimurium` | partial | 4619 | 4554 | +65 | 101.43 | 88 | 85 | +3 | 103.53 | 15 | 22 | -7 | 68.18 |
| `afabrum` | partial | 5330 | 5267 | +63 | 101.20 | 55 | 53 | +2 | 103.77 | 8 | 12 | -4 | 66.67 |
| `smeliloti` | partial | 6336 | 6287 | +49 | 100.78 | 56 | 55 | +1 | 101.82 | 6 | 9 | -3 | 66.67 |
| `mgryphiswaldense` | partial | 3945 | 3916 | +29 | 100.74 | 53 | 50 | +3 | 106.00 | 4 | 6 | -2 | 66.67 |
| `soneidensis` | partial | 4569 | 4491 | +78 | 101.74 | 106 | 104 | +2 | 101.92 | 19 | 28 | -9 | 67.86 |
| `cnecator` | partial | 7907 | 7887 | +20 | 100.25 | 73 | 69 | +4 | 105.80 | 10 | 15 | -5 | 66.67 |
| `zmobilis` | partial | 1859 | 1852 | +7 | 100.38 | 49 | 48 | +1 | 102.08 | 4 | 6 | -2 | 66.67 |
| `tmaritima` | partial | 1854 | 1866 | -12 | 99.36 | 47 | 46 | +1 | 102.17 | 2 | 3 | -1 | 66.67 |
| `neuropea` | partial | 2621 | 2631 | -10 | 99.62 | 42 | 41 | +1 | 102.44 | 2 | 3 | -1 | 66.67 |
| `dvulgaris` | partial | 3299 | 3249 | +50 | 101.54 | 72 | 67 | +5 | 107.46 | 11 | 16 | -5 | 68.75 |
| `rpalustris` | partial | 4947 | 4975 | -28 | 99.44 | 51 | 49 | +2 | 104.08 | 4 | 6 | -2 | 66.67 |

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
| `hpylori` | partial | delta CDS=-4, tRNA=+1, rRNA=-4 |
| `cjejuni` | partial | delta CDS=+101, tRNA=+2, rRNA=-6 |
| `nmeningitidis` | partial | delta CDS=-92, tRNA=+2, rRNA=-4 |
| `hinfluenzae` | partial | delta CDS=+13, tRNA=+2, rRNA=-6 |
| `ctrachomatis` | partial | delta CDS=+12, rRNA=-2 |
| `tpallidum` | partial | delta CDS=-5, rRNA=-2 |
| `bburgdorferi` | partial | delta CDS=+6, rRNA=-2 |
| `mgenitalium` | partial | delta CDS=-30, rRNA=-1 |
| `abaumannii` | partial | delta CDS=+65, tRNA=+3, rRNA=-6 |
| `kpneumoniae` | partial | delta CDS=-247, tRNA=+29, rRNA=-8 |
| `styphimurium` | partial | delta CDS=+65, tRNA=+3, rRNA=-7 |
| `afabrum` | partial | delta CDS=+63, tRNA=+2, rRNA=-4 |
| `smeliloti` | partial | delta CDS=+49, tRNA=+1, rRNA=-3 |
| `mgryphiswaldense` | partial | delta CDS=+29, tRNA=+3, rRNA=-2 |
| `soneidensis` | partial | delta CDS=+78, tRNA=+2, rRNA=-9 |
| `cnecator` | partial | delta CDS=+20, tRNA=+4, rRNA=-5 |
| `zmobilis` | partial | delta CDS=+7, tRNA=+1, rRNA=-2 |
| `tmaritima` | partial | delta CDS=-12, tRNA=+1, rRNA=-1 |
| `neuropea` | partial | delta CDS=-10, tRNA=+1, rRNA=-1 |
| `dvulgaris` | partial | delta CDS=+50, tRNA=+5, rRNA=-5 |
| `rpalustris` | partial | delta CDS=-28, tRNA=+2, rRNA=-2 |

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
| `hpylori` | Helicobacter pylori | 26695 | annotated | `data/reference_genomes/hpylori/hpylori_genomic.gff.gz` | `results/hpylori/hpylori.summary.json` |
| `cjejuni` | Campylobacter jejuni | NCTC 11168 | annotated | `data/reference_genomes/cjejuni/cjejuni_genomic.gff.gz` | `results/cjejuni/cjejuni.summary.json` |
| `nmeningitidis` | Neisseria meningitidis | MC58 | annotated | `data/reference_genomes/nmeningitidis/nmeningitidis_genomic.gff.gz` | `results/nmeningitidis/nmeningitidis.summary.json` |
| `hinfluenzae` | Haemophilus influenzae | Rd KW20 | annotated | `data/reference_genomes/hinfluenzae/hinfluenzae_genomic.gff.gz` | `results/hinfluenzae/hinfluenzae.summary.json` |
| `ctrachomatis` | Chlamydia trachomatis | D/UW-3/CX | annotated | `data/reference_genomes/ctrachomatis/ctrachomatis_genomic.gff.gz` | `results/ctrachomatis/ctrachomatis.summary.json` |
| `tpallidum` | Treponema pallidum | Nichols | annotated | `data/reference_genomes/tpallidum/tpallidum_genomic.gff.gz` | `results/tpallidum/tpallidum.summary.json` |
| `bburgdorferi` | Borrelia burgdorferi | B31 | annotated | `data/reference_genomes/bburgdorferi/bburgdorferi_genomic.gff.gz` | `results/bburgdorferi/bburgdorferi.summary.json` |
| `mgenitalium` | Mycoplasma genitalium | G37 | annotated | `data/reference_genomes/mgenitalium/mgenitalium_genomic.gff.gz` | `results/mgenitalium/mgenitalium.summary.json` |
| `abaumannii` | Acinetobacter baumannii | ATCC 17978 | annotated | `data/reference_genomes/abaumannii/abaumannii_genomic.gff.gz` | `results/abaumannii/abaumannii.summary.json` |
| `kpneumoniae` | Klebsiella pneumoniae | MGH 78578 | annotated | `data/reference_genomes/kpneumoniae/kpneumoniae_genomic.gff.gz` | `results/kpneumoniae/kpneumoniae.summary.json` |
| `styphimurium` | Salmonella enterica | Typhimurium LT2 | annotated | `data/reference_genomes/styphimurium/styphimurium_genomic.gff.gz` | `results/styphimurium/styphimurium.summary.json` |
| `afabrum` | Agrobacterium fabrum | C58 | annotated | `data/reference_genomes/afabrum/afabrum_genomic.gff.gz` | `results/afabrum/afabrum.summary.json` |
| `smeliloti` | Sinorhizobium meliloti | 1021 | annotated | `data/reference_genomes/smeliloti/smeliloti_genomic.gff.gz` | `results/smeliloti/smeliloti.summary.json` |
| `mgryphiswaldense` | Magnetospirillum gryphiswaldense | MSR-1 | annotated | `data/reference_genomes/mgryphiswaldense/mgryphiswaldense_genomic.gff.gz` | `results/mgryphiswaldense/mgryphiswaldense.summary.json` |
| `soneidensis` | Shewanella oneidensis | MR-1 | annotated | `data/reference_genomes/soneidensis/soneidensis_genomic.gff.gz` | `results/soneidensis/soneidensis.summary.json` |
| `cnecator` | Cupriavidus necator | H16 | annotated | `data/reference_genomes/cnecator/cnecator_genomic.gff.gz` | `results/cnecator/cnecator.summary.json` |
| `zmobilis` | Zymomonas mobilis | ZM4 | annotated | `data/reference_genomes/zmobilis/zmobilis_genomic.gff.gz` | `results/zmobilis/zmobilis.summary.json` |
| `tmaritima` | Thermotoga maritima | MSB8 | annotated | `data/reference_genomes/tmaritima/tmaritima_genomic.gff.gz` | `results/tmaritima/tmaritima.summary.json` |
| `neuropea` | Nitrosomonas europaea | ATCC 19718 | annotated | `data/reference_genomes/neuropea/neuropea_genomic.gff.gz` | `results/neuropea/neuropea.summary.json` |
| `dvulgaris` | Desulfovibrio vulgaris | Hildenborough | annotated | `data/reference_genomes/dvulgaris/dvulgaris_genomic.gff.gz` | `results/dvulgaris/dvulgaris.summary.json` |
| `rpalustris` | Rhodopseudomonas palustris | CGA009 | annotated | `data/reference_genomes/rpalustris/rpalustris_genomic.gff.gz` | `results/rpalustris/rpalustris.summary.json` |

## Agregado (micro)

| clase | sistema | referencia | delta | % sistema/referencia |
|---|---:|---:|---:|---:|
| CDS | 139264 | 138574 | +690 | 100.50 |
| tRNA | 2367 | 2281 | +86 | 103.77 |
| rRNA | 329 | 506 | -177 | 65.02 |

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
