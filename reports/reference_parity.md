# Paridad de anotación: sistema vs GFF de referencia

**Generado:** `2026-11-26T02:06:46Z`  
**Manifiesto:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/reference_genomes.json`  
**Raíz de datos:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`

## Alcance

Compara los conteos de **CDS, tRNA y rRNA** del resumen del sistema (`results/<id>/<id>.summary.json`) contra los features del GFF de referencia (`annotation_path`) para cada genoma con `status` `annotated` o `available` y `system_summary` en el manifiesto.

Esta versión **no evalúa función, productos, coordenadas ni identidad de secuencia**. Tampoco modifica `results/`, los GFF ni el manifiesto.

## Cobertura

- Registros con status `annotated`/`available`: **42**
- Evaluados (con `system_summary`): **40** (ok **0**, partial **40**, not_available **0**)
- Omitidos: **2** (2 elegibles sin `system_summary`; 0 con otro status)

### Omitidos con status elegible

| genoma | status | motivo |
|---|---|---|
| `ecoli_k12` | available | sin system_summary en el manifiesto |
| `subtilis` | available | sin system_summary en el manifiesto |

## Resultados por genoma

| genoma | estado | CDS sis | CDS ref | Δ CDS | % CDS | tRNA sis | tRNA ref | Δ tRNA | % tRNA | rRNA sis | rRNA ref | Δ rRNA | % rRNA |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `mtb_h37rv` | partial | 3986 | 3906 | +80 | 102.05 | 45 | 45 | 0 | 100.00 | 3 | 3 | 0 | 100.00 |
| `cglutamicum` | partial | 3041 | 2999 | +42 | 101.40 | 60 | 59 | +1 | 101.69 | 18 | 18 | 0 | 100.00 |
| `pputida` | partial | 5507 | 5471 | +36 | 100.66 | 75 | 75 | 0 | 100.00 | 22 | 22 | 0 | 100.00 |
| `paeruginosa` | partial | 5715 | 5572 | +143 | 102.57 | 64 | 63 | +1 | 101.59 | 12 | 13 | -1 | 92.31 |
| `crescentus` | partial | 3738 | 3886 | -148 | 96.19 | 52 | 51 | +1 | 101.96 | 6 | 6 | 0 | 100.00 |
| `synechocystis` | partial | 3780 | 3698 | +82 | 102.22 | 41 | 41 | 0 | 100.00 | 6 | 6 | 0 | 100.00 |
| `saureus` | partial | 2675 | 2767 | -92 | 96.68 | 61 | 59 | +2 | 103.39 | 16 | 16 | 0 | 100.00 |
| `spneumoniae` | partial | 2168 | 2117 | +51 | 102.41 | 59 | 58 | +1 | 101.72 | 12 | 12 | 0 | 100.00 |
| `btheta` | partial | 4945 | 4801 | +144 | 103.00 | 71 | 70 | +1 | 101.43 | 15 | 15 | 0 | 100.00 |
| `scoelicolor` | partial | 7738 | 7755 | -17 | 99.78 | 72 | 66 | +6 | 109.09 | 18 | 18 | 0 | 100.00 |
| `mxanthus` | partial | 7325 | 7306 | +19 | 100.26 | 69 | 66 | +3 | 104.55 | 12 | 12 | 0 | 100.00 |
| `vcholerae` | partial | 3666 | 3665 | +1 | 100.03 | 102 | 102 | 0 | 100.00 | 31 | 31 | 0 | 100.00 |
| `cdifficile` | partial | 3700 | 3611 | +89 | 102.46 | 91 | 90 | +1 | 101.11 | 35 | 35 | 0 | 100.00 |
| `llactis` | partial | 2602 | 2567 | +35 | 101.36 | 71 | 70 | +1 | 101.43 | 19 | 19 | 0 | 100.00 |
| `lmonocytogenes` | partial | 2921 | 2867 | +54 | 101.88 | 67 | 67 | 0 | 100.00 | 18 | 18 | 0 | 100.00 |
| `efaecalis` | partial | 2730 | 2690 | +40 | 101.49 | 53 | 51 | +2 | 103.92 | 12 | 12 | 0 | 100.00 |
| `dradiodurans` | partial | 3150 | 3148 | +2 | 100.06 | 51 | 50 | +1 | 102.00 | 9 | 9 | 0 | 100.00 |
| `tthermophilus` | partial | 2227 | 2228 | -1 | 99.96 | 49 | 48 | +1 | 102.08 | 6 | 6 | 0 | 100.00 |
| `hpylori` | partial | 1580 | 1584 | -4 | 99.75 | 37 | 36 | +1 | 102.78 | 6 | 6 | 0 | 100.00 |
| `cjejuni` | partial | 1673 | 1572 | +101 | 106.42 | 45 | 43 | +2 | 104.65 | 9 | 9 | 0 | 100.00 |
| `nmeningitidis` | partial | 2050 | 2142 | -92 | 95.70 | 61 | 59 | +2 | 103.39 | 12 | 12 | 0 | 100.00 |
| `hinfluenzae` | partial | 1819 | 1806 | +13 | 100.72 | 61 | 59 | +2 | 103.39 | 19 | 19 | 0 | 100.00 |
| `ctrachomatis` | partial | 899 | 887 | +12 | 101.35 | 37 | 37 | 0 | 100.00 | 6 | 6 | 0 | 100.00 |
| `tpallidum` | partial | 982 | 987 | -5 | 99.49 | 45 | 45 | 0 | 100.00 | 6 | 6 | 0 | 100.00 |
| `bburgdorferi` | partial | 1534 | 1528 | +6 | 100.39 | 32 | 32 | 0 | 100.00 | 5 | 5 | 0 | 100.00 |
| `mgenitalium` | partial | 494 | 524 | -30 | 94.27 | 36 | 36 | 0 | 100.00 | 3 | 3 | 0 | 100.00 |
| `abaumannii` | partial | 3801 | 3736 | +65 | 101.74 | 77 | 74 | +3 | 104.05 | 18 | 18 | 0 | 100.00 |
| `kpneumoniae` | partial | 5532 | 5779 | -247 | 95.73 | 91 | 62 | +29 | 146.77 | 25 | 25 | 0 | 100.00 |
| `styphimurium` | partial | 4619 | 4554 | +65 | 101.43 | 88 | 85 | +3 | 103.53 | 22 | 22 | 0 | 100.00 |
| `afabrum` | partial | 5330 | 5267 | +63 | 101.20 | 55 | 53 | +2 | 103.77 | 12 | 12 | 0 | 100.00 |
| `smeliloti` | partial | 6336 | 6287 | +49 | 100.78 | 56 | 55 | +1 | 101.82 | 9 | 9 | 0 | 100.00 |
| `rsphaeroides` | partial | 4417 | 4406 | +11 | 100.25 | 56 | 53 | +3 | 105.66 | 9 | 9 | 0 | 100.00 |
| `mgryphiswaldense` | partial | 3945 | 3916 | +29 | 100.74 | 53 | 50 | +3 | 106.00 | 6 | 6 | 0 | 100.00 |
| `soneidensis` | partial | 4569 | 4491 | +78 | 101.74 | 106 | 104 | +2 | 101.92 | 28 | 28 | 0 | 100.00 |
| `cnecator` | partial | 7907 | 7887 | +20 | 100.25 | 73 | 69 | +4 | 105.80 | 15 | 15 | 0 | 100.00 |
| `zmobilis` | partial | 1859 | 1852 | +7 | 100.38 | 49 | 48 | +1 | 102.08 | 6 | 6 | 0 | 100.00 |
| `tmaritima` | partial | 1854 | 1866 | -12 | 99.36 | 47 | 46 | +1 | 102.17 | 3 | 3 | 0 | 100.00 |
| `neuropea` | partial | 2621 | 2631 | -10 | 99.62 | 42 | 41 | +1 | 102.44 | 3 | 3 | 0 | 100.00 |
| `dvulgaris` | partial | 3299 | 3249 | +50 | 101.54 | 72 | 67 | +5 | 107.46 | 16 | 16 | 0 | 100.00 |
| `rpalustris` | partial | 4947 | 4975 | -28 | 99.44 | 51 | 49 | +2 | 104.08 | 6 | 6 | 0 | 100.00 |

## Estado y motivos

| genoma | estado | motivo |
|---|---|---|
| `mtb_h37rv` | partial | delta CDS=+80 |
| `cglutamicum` | partial | delta CDS=+42, tRNA=+1 |
| `pputida` | partial | delta CDS=+36 |
| `paeruginosa` | partial | delta CDS=+143, tRNA=+1, rRNA=-1 |
| `crescentus` | partial | delta CDS=-148, tRNA=+1 |
| `synechocystis` | partial | delta CDS=+82 |
| `saureus` | partial | delta CDS=-92, tRNA=+2 |
| `spneumoniae` | partial | delta CDS=+51, tRNA=+1 |
| `btheta` | partial | delta CDS=+144, tRNA=+1 |
| `scoelicolor` | partial | delta CDS=-17, tRNA=+6 |
| `mxanthus` | partial | delta CDS=+19, tRNA=+3 |
| `vcholerae` | partial | delta CDS=+1 |
| `cdifficile` | partial | delta CDS=+89, tRNA=+1 |
| `llactis` | partial | delta CDS=+35, tRNA=+1 |
| `lmonocytogenes` | partial | delta CDS=+54 |
| `efaecalis` | partial | delta CDS=+40, tRNA=+2 |
| `dradiodurans` | partial | delta CDS=+2, tRNA=+1 |
| `tthermophilus` | partial | delta CDS=-1, tRNA=+1 |
| `hpylori` | partial | delta CDS=-4, tRNA=+1 |
| `cjejuni` | partial | delta CDS=+101, tRNA=+2 |
| `nmeningitidis` | partial | delta CDS=-92, tRNA=+2 |
| `hinfluenzae` | partial | delta CDS=+13, tRNA=+2 |
| `ctrachomatis` | partial | delta CDS=+12 |
| `tpallidum` | partial | delta CDS=-5 |
| `bburgdorferi` | partial | delta CDS=+6 |
| `mgenitalium` | partial | delta CDS=-30 |
| `abaumannii` | partial | delta CDS=+65, tRNA=+3 |
| `kpneumoniae` | partial | delta CDS=-247, tRNA=+29 |
| `styphimurium` | partial | delta CDS=+65, tRNA=+3 |
| `afabrum` | partial | delta CDS=+63, tRNA=+2 |
| `smeliloti` | partial | delta CDS=+49, tRNA=+1 |
| `rsphaeroides` | partial | delta CDS=+11, tRNA=+3 |
| `mgryphiswaldense` | partial | delta CDS=+29, tRNA=+3 |
| `soneidensis` | partial | delta CDS=+78, tRNA=+2 |
| `cnecator` | partial | delta CDS=+20, tRNA=+4 |
| `zmobilis` | partial | delta CDS=+7, tRNA=+1 |
| `tmaritima` | partial | delta CDS=-12, tRNA=+1 |
| `neuropea` | partial | delta CDS=-10, tRNA=+1 |
| `dvulgaris` | partial | delta CDS=+50, tRNA=+5 |
| `rpalustris` | partial | delta CDS=-28, tRNA=+2 |

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
| `rsphaeroides` | Rhodobacter sphaeroides | 2.4.1 | annotated | `data/reference_genomes/rsphaeroides/rsphaeroides_genomic.gff.gz` | `results/rsphaeroides/rsphaeroides.summary.json` |
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
| CDS | 143681 | 142980 | +701 | 100.49 |
| tRNA | 2423 | 2334 | +89 | 103.81 |
| rRNA | 514 | 515 | -1 | 99.81 |

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
