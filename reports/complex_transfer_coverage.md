# Cobertura de la transferencia de complejos K-12

Generado: `2026-11-25T02:01:21Z`

Alcance: genomas con `data/complexes/transfer_*.tsv` o `data/complexes/model_*.complexes.tsv`. Solo lectura local: no se descargó ComplexPortal, no se consultó la red y no se modificó ningún complejo.

`not_available` significa que **no hay catálogo ComplexPortal independiente local** para ese genoma. No significa que ComplexPortal no tenga una entrada remota. La transferencia K-12 y el modelo K-12 son inferencias: no son validación.

## Resumen

| métrica | valor |
|---|---:|
| genomas auditados | 400 |
| con catálogo ComplexPortal independiente | 0 |
| sin catálogo independiente (`not_available`) | 400 |
| genomas con transferencia | 399 |
| genomas sin transferencia | 1 |
| genomas con modelo K-12 | 400 |
| genomas sin modelo K-12 | 0 |
| complejos transferidos | 19.085 |
| transferidos con `frac_mapped >= 0.8` (sospechosos) | 8.408 |
| transferidos con `frac_mapped < 0.8` (candidatos) | 10.677 |
| subunidades esperadas (suma `n_subunits`) | 95.278 |
| subunidades mapeadas (suma `n_mapped`) | 56.222 |
| cobertura global de subunidades | 0.5901 |
| candidatos del modelo K-12 | 281.090 |
| coincidencias con firmas curadas | 1.180 |
| firmas completas / incompletas | 247 / 933 |
| genomas con alguna firma | 382 |

## Estados

| estado | n | significado |
|---|---:|---|
| `transfer_suspected` | 8.408 | complejo transferido con frac_mapped >= 0.8; sospecha fuerte |
| `transfer_candidate` | 10.677 | complejo transferido con frac_mapped < 0.8; candidato |
| `transfer_stoich_complete` | 0 | observed_n >= expected_n en el TSV de transferencia |
| `transfer_stoich_incomplete` | 1.785 | expected_n presente y observed_n < expected_n |
| `transfer_stoich_unknown` | 17.300 | sin expected_n; estequiometría no evaluable |
| `transfer_full_subunits` | 7.398 | n_mapped == n_subunits del catálogo K-12 |
| `signature_complete` | 247 | todos los genes esenciales de la firma están en la transferencia |
| `signature_incomplete` | 933 | firma con >=2 genes pero sin todos los esenciales |
| `model_candidates` | 281.090 | filas del modelo K-12 con status complex_candidate |
| `model_complete` | 0 | filas del modelo K-12 con status complex_complete |
| `model_incomplete` | 0 | filas del modelo K-12 con status complex_incomplete |
| `model_other` | 0 | filas del modelo K-12 con otro status |

`transfer_full_subunits` cuenta filas con `n_mapped == n_subunits`; `transfer_stoich_complete` cuenta la estequiometría declarada (`observed_n >= expected_n`). No son lo mismo: el catálogo K-12 puede tener más subunidades que el `expected_n` de `properties`.

## Cobertura de subunidades

| métrica | valor |
|---|---:|
| subunidades esperadas | 95.278 |
| subunidades mapeadas | 56.222 |
| cobertura global | 0.5901 |
| media no ponderada de las medias de `frac_mapped` por genoma | 0.7041 |
| complejos con todas las subunidades del catálogo | 7.398 |

## Genomas con más complejos transferidos

| genome | transferidos | sospechosos | cobertura | firmas | modelo K-12 |
|---|---:|---:|---:|---:|---:|
| `GCF_019048065.1_ASM1904806v1` | 237 | 233 | 0.9912 | 13 | 892 |
| `GCF_019048105.1_ASM1904810v1` | 174 | 148 | 0.9197 | 12 | 846 |
| `GCF_005144435.2_ASM514443v2` | 108 | 61 | 0.7115 | 6 | 1.160 |
| `GCF_028885455.1_ASM2888545v1` | 105 | 79 | 0.8677 | 10 | 982 |
| `GCF_001886435.1_ASM188643v1` | 98 | 76 | 0.8733 | 9 | 874 |
| `GCF_018243235.1_ASM1824323v1` | 98 | 58 | 0.7623 | 10 | 927 |
| `GCF_032698475.1_ASM3269847v1` | 97 | 51 | 0.7400 | 8 | 972 |
| `GCF_002157165.1_ASM215716v1` | 96 | 59 | 0.7662 | 9 | 764 |
| `GCF_042847215.1_ASM4284721v1` | 85 | 59 | 0.8084 | 6 | 693 |
| `GCF_918320885.1_Harlequin_endosymbiont1` | 84 | 66 | 0.8761 | 5 | 363 |

## Genomas con menor cobertura de subunidades

| genome | transferidos | cobertura | subunidades mapeadas / esperadas |
|---|---:|---:|---:|
| `GCF_000508245.1_ASM50824v1` | 13 | 0.2545 | 28 / 110 |
| `GCF_900660665.1_57612_D01-3` | 13 | 0.3091 | 34 / 110 |
| `GCF_041026625.1_ASM4102662v1` | 14 | 0.3304 | 37 / 112 |
| `GCF_056388415.1_ASM5638841v1` | 14 | 0.3482 | 39 / 112 |
| `GCF_000941075.1_ASM94107v1` | 13 | 0.3545 | 39 / 110 |
| `GCF_900660485.1_50648_E01-3` | 16 | 0.3644 | 43 / 118 |
| `GCF_018866205.1_ASM1886620v1` | 20 | 0.3904 | 57 / 146 |
| `GCF_026184435.1_ASM2618443v1` | 34 | 0.4223 | 87 / 206 |
| `GCF_029677705.1_ASM2967770v1` | 21 | 0.4286 | 57 / 133 |
| `GCF_047159735.1_ASM4715973v1` | 21 | 0.4306 | 62 / 144 |

## Limitaciones

- **Sin catálogo independiente.** Los genomas del panel no tienen un catálogo ComplexPortal propio; el único catálogo local es el de K-12. Por eso la columna `independent_catalog` queda en `not_available` y no se usa como validación.
- **La transferencia no es validación.** Es RBH K-12, con identidad >= 30 % y evalue <= 1e-5, pero sin cobertura de alineamiento en los TSV actuales. `frac_mapped` bajo puede significar complejo ausente, anotación incompleta u ortología no detectada.
- **Estequiometría.** `expected_n` solo está presente en una parte de las filas; `stoich_status=unknown` no es un fallo del complejo, es un dato no evaluable. No hay ninguna fila `complete` en el panel actual.
- **Firmas curadas.** Se comparan solo nombres de gen K-12; no son un catálogo exhaustivo ni prueban interacción física. `signature_complete` significa que están los genes esenciales de la firma en la transferencia, no que el complejo exista.
- **Modelo K-12.** Cada fila es `complex_candidate`: una sospecha por vecindad/operón, sin estructura ni estequiometría. No cierra brechas.
- **Sin red.** No se descargó ComplexPortal ni se consultó el listado remoto. Un `not_available` es una ausencia local, no una afirmación sobre la base remota.

## Detalle por genoma

| genome | catálogo independiente | transferidos | sospechosos | cobertura | firmas (compl./inc.) | modelo K-12 | estoq. (C/I/U) |
|---|---|---:|---:|---:|---:|---:|---:|
| `GCF_000006685.1_ASM668v1` | `not_available` | 31 | 11 | 0.4432 | 0 / 1 | 198 | 0 / 4 / 27 |
| `GCF_000006985.1_ASM698v1` | `not_available` | 59 | 21 | 0.5654 | 1 / 2 | 519 | 0 / 4 / 55 |
| `GCF_000007205.1_ASM720v1` | `not_available` | 33 | 11 | 0.4866 | 0 / 1 | 225 | 0 / 4 / 29 |
| `GCF_000007605.1_ASM760v1` | `not_available` | 34 | 12 | 0.4921 | 0 / 1 | 223 | 0 / 4 / 30 |
| `GCF_000007985.2_ASM798v2` | `not_available` | 72 | 38 | 0.6780 | 2 / 5 | 627 | 0 / 6 / 66 |
| `GCF_000008625.1_ASM862v1` | `not_available` | 38 | 19 | 0.5744 | 1 / 2 | 403 | 0 / 4 / 34 |
| `GCF_000008725.1_ASM872v1` | `not_available` | 30 | 10 | 0.4368 | 0 / 1 | 193 | 0 / 4 / 26 |
| `GCF_000010305.1_ASM1030v1` | `not_available` | 59 | 26 | 0.6049 | 1 / 4 | 844 | 0 / 6 / 53 |
| `GCF_000010985.1_ASM1098v1` | `not_available` | 51 | 25 | 0.6221 | 2 / 5 | 416 | 0 / 4 / 47 |
| `GCF_000011905.1_ASM1190v1` | `not_available` | 32 | 13 | 0.5301 | 0 / 2 | 320 | 0 / 4 / 28 |
| `GCF_000012485.1_ASM1248v1` | `not_available` | 62 | 22 | 0.5828 | 2 / 2 | 457 | 0 / 4 / 58 |
| `GCF_000014025.1_ASM1402v1` | `not_available` | 49 | 19 | 0.5202 | 1 / 4 | 465 | 0 / 5 / 44 |
| `GCF_000014125.1_ASM1412v1` | `not_available` | 48 | 21 | 0.6198 | 0 / 1 | 347 | 0 / 4 / 44 |
| `GCF_000015045.1_ASM1504v1` | `not_available` | 74 | 34 | 0.6325 | 2 / 4 | 657 | 0 / 4 / 70 |
| `GCF_000015125.1_ASM1512v1` | `not_available` | 62 | 26 | 0.5754 | 2 / 2 | 629 | 0 / 4 / 58 |
| `GCF_000017685.1_ASM1768v1` | `not_available` | 48 | 22 | 0.5844 | 1 / 3 | 834 | 0 / 5 / 43 |
| `GCF_000017805.1_ASM1780v1` | `not_available` | 46 | 20 | 0.6404 | 2 / 3 | 1.061 | 0 / 5 / 41 |
| `GCF_000017865.1_ASM1786v1` | `not_available` | 35 | 13 | 0.5816 | 0 / 1 | 369 | 0 / 4 / 31 |
| `GCF_000018605.1_ASM1860v1` | `not_available` | 42 | 12 | 0.4652 | 0 / 2 | 314 | 0 / 4 / 38 |
| `GCF_000018865.1_ASM1886v1` | `not_available` | 52 | 19 | 0.6092 | 1 / 3 | 905 | 0 / 5 / 47 |
| `GCF_000018945.1_ASM1894v1` | `not_available` | 46 | 14 | 0.5331 | 0 / 2 | 335 | 0 / 4 / 42 |
| `GCF_000019965.1_ASM1996v1` | `not_available` | 59 | 29 | 0.6245 | 1 / 3 | 1.013 | 0 / 5 / 54 |
| `GCF_000020145.1_ASM2014v1` | `not_available` | 27 | 10 | 0.5032 | 0 / 1 | 311 | 0 / 4 / 23 |
| `GCF_000020465.1_ASM2046v1` | `not_available` | 63 | 25 | 0.5788 | 2 / 3 | 555 | 0 / 5 / 58 |
| `GCF_000020505.1_ASM2050v1` | `not_available` | 62 | 21 | 0.5734 | 3 / 2 | 438 | 0 / 4 / 58 |
| `GCF_000020525.1_ASM2052v1` | `not_available` | 49 | 22 | 0.6008 | 3 / 2 | 607 | 0 / 4 / 45 |
| `GCF_000020625.1_ASM2062v1` | `not_available` | 53 | 24 | 0.5771 | 2 / 2 | 474 | 0 / 4 / 49 |
| `GCF_000020645.1_ASM2064v1` | `not_available` | 62 | 22 | 0.5514 | 2 / 3 | 621 | 0 / 5 / 57 |
| `GCF_000020945.1_ASM2094v1` | `not_available` | 35 | 11 | 0.4873 | 0 / 1 | 246 | 0 / 4 / 31 |
| `GCF_000020965.1_ASM2096v1` | `not_available` | 43 | 17 | 0.5785 | 0 / 1 | 297 | 0 / 4 / 39 |
| `GCF_000020985.1_ASM2098v1` | `not_available` | 48 | 21 | 0.6154 | 0 / 3 | 349 | 0 / 4 / 44 |
| `GCF_000021545.1_ASM2154v1` | `not_available` | 38 | 19 | 0.5800 | 1 / 2 | 327 | 0 / 4 / 34 |
| `GCF_000021565.1_ASM2156v1` | `not_available` | 45 | 21 | 0.5819 | 1 / 2 | 403 | 0 / 4 / 41 |
| `GCF_000021645.1_ASM2164v1` | `not_available` | 42 | 19 | 0.5760 | 0 / 1 | 295 | 0 / 4 / 38 |
| `GCF_000021685.1_ASM2168v1` | `not_available` | 55 | 22 | 0.5751 | 0 / 5 | 564 | 0 / 5 / 50 |
| `GCF_000021905.1_ASM2190v1` | `not_available` | 64 | 29 | 0.6448 | 2 / 3 | 1.040 | 0 / 4 / 60 |
| `GCF_000021945.1_ASM2194v1` | `not_available` | 52 | 22 | 0.6187 | 1 / 3 | 836 | 0 / 5 / 47 |
| `GCF_000022145.1_ASM2214v1` | `not_available` | 64 | 34 | 0.6520 | 3 / 4 | 989 | 0 / 6 / 58 |
| `GCF_000022565.1_ASM2256v1` | `not_available` | 53 | 25 | 0.5992 | 1 / 4 | 754 | 0 / 5 / 48 |
| `GCF_000023325.1_ASM2332v1` | `not_available` | 34 | 12 | 0.4975 | 0 / 2 | 425 | 0 / 4 / 30 |
| `GCF_000024405.1_ASM2440v1` | `not_available` | 56 | 31 | 0.6506 | 0 / 1 | 707 | 0 / 4 / 52 |
| `GCF_000024425.1_ASM2442v1` | `not_available` | 53 | 25 | 0.6175 | 1 / 3 | 637 | 0 / 5 / 48 |
| `GCF_000024845.1_ASM2484v1` | `not_available` | 52 | 25 | 0.5833 | 1 / 3 | 626 | 0 / 5 / 47 |
| `GCF_000024905.1_ASM2490v1` | `not_available` | 36 | 13 | 0.5773 | 0 / 0 | 288 | 0 / 3 / 33 |
| `GCF_000024985.1_ASM2498v1` | `not_available` | 57 | 27 | 0.6066 | 0 / 3 | 740 | 0 / 5 / 52 |
| `GCF_000025005.1_ASM2500v1` | `not_available` | 50 | 23 | 0.6335 | 1 / 3 | 613 | 0 / 5 / 45 |
| `GCF_000025605.1_ASM2560v1` | `not_available` | 37 | 12 | 0.5308 | 2 / 2 | 329 | 0 / 4 / 33 |
| `GCF_000025725.1_ASM2572v1` | `not_available` | 60 | 25 | 0.6143 | 1 / 3 | 588 | 0 / 4 / 56 |
| `GCF_000025885.1_ASM2588v1` | `not_available` | 43 | 13 | 0.5596 | 0 / 0 | 306 | 0 / 3 / 40 |
| `GCF_000025905.1_ASM2590v1` | `not_available` | 46 | 21 | 0.5769 | 0 / 3 | 707 | 0 / 5 / 41 |
| `GCF_000025945.1_ASM2594v1` | `not_available` | 51 | 25 | 0.6109 | 0 / 2 | 659 | 0 / 4 / 47 |
| `GCF_000091545.1_ASM9154v1` | `not_available` | 43 | 22 | 0.6009 | 0 / 3 | 469 | 0 / 5 / 38 |
| `GCF_000092425.1_ASM9242v1` | `not_available` | 51 | 21 | 0.5673 | 1 / 3 | 627 | 0 / 5 / 46 |
| `GCF_000092785.1_ASM9278v1` | `not_available` | 37 | 15 | 0.5158 | 0 / 3 | 434 | 0 / 6 / 31 |
| `GCF_000143725.1_ASM14372v1` | `not_available` | 44 | 16 | 0.5152 | 0 / 1 | 529 | 0 / 4 / 40 |
| `GCF_000146065.2_ASM14606v1` | `not_available` | 33 | 13 | 0.5763 | 0 / 0 | 317 | 0 / 3 / 30 |
| `GCF_000146505.1_ASM14650v1` | `not_available` | 39 | 14 | 0.5385 | 1 / 1 | 722 | 0 / 3 / 36 |
| `GCF_000147715.2_ASM14771v3` | `not_available` | 45 | 13 | 0.4938 | 0 / 2 | 506 | 0 / 4 / 41 |
| `GCF_000164905.1_ASM16490v1` | `not_available` | 40 | 16 | 0.5117 | 1 / 2 | 402 | 0 / 4 / 36 |
| `GCF_000172155.1_ASM17215v1` | `not_available` | 55 | 24 | 0.6190 | 0 / 3 | 1.431 | 0 / 5 / 50 |
| `GCF_000177635.2_ASM17763v2` | `not_available` | 55 | 24 | 0.6391 | 2 / 2 | 508 | 0 / 4 / 51 |
| `GCF_000179915.2_ASM17991v2` | `not_available` | 48 | 21 | 0.5560 | 0 / 3 | 1.035 | 0 / 4 / 44 |
| `GCF_000183405.1_ASM18340v1` | `not_available` | 43 | 21 | 0.5991 | 1 / 3 | 319 | 0 / 4 / 39 |
| `GCF_000185805.1_ASM18580v1` | `not_available` | 41 | 21 | 0.5894 | 0 / 2 | 365 | 0 / 4 / 37 |
| `GCF_000186245.1_ASM18624v1` | `not_available` | 41 | 21 | 0.5625 | 1 / 3 | 395 | 0 / 4 / 37 |
| `GCF_000191045.1_ASM19104v1` | `not_available` | 42 | 21 | 0.5767 | 0 / 2 | 285 | 0 / 4 / 38 |
| `GCF_000195335.1_ASM19533v1` | `not_available` | 49 | 26 | 0.6087 | 1 / 3 | 464 | 0 / 5 / 44 |
| `GCF_000195975.1_ASM19597v1` | `not_available` | 41 | 16 | 0.5093 | 0 / 3 | 649 | 0 / 4 / 37 |
| `GCF_000196115.1_ASM19611v1` | `not_available` | 48 | 24 | 0.5957 | 1 / 3 | 1.283 | 0 / 6 / 42 |
| `GCF_000196175.1_ASM19617v1` | `not_available` | 49 | 24 | 0.6255 | 1 / 3 | 810 | 0 / 6 / 43 |
| `GCF_000199675.1_ASM19967v1` | `not_available` | 46 | 17 | 0.5537 | 0 / 2 | 668 | 0 / 4 / 42 |
| `GCF_000204255.1_ASM20425v1` | `not_available` | 35 | 13 | 0.4869 | 0 / 1 | 221 | 0 / 4 / 31 |
| `GCF_000210915.2_ASM21091v2` | `not_available` | 47 | 27 | 0.6202 | 1 / 3 | 756 | 0 / 5 / 42 |
| `GCF_000212395.1_ASM21239v1` | `not_available` | 40 | 17 | 0.5381 | 0 / 2 | 289 | 0 / 4 / 36 |
| `GCF_000212675.2_ASM21267v2` | `not_available` | 54 | 29 | 0.6292 | 1 / 2 | 785 | 0 / 4 / 50 |
| `GCF_000214155.1_ASM21415v1` | `not_available` | 38 | 14 | 0.5628 | 0 / 2 | 997 | 0 / 6 / 32 |
| `GCF_000215975.1_ASM21597v1` | `not_available` | 44 | 18 | 0.5672 | 0 / 2 | 326 | 0 / 4 / 40 |
| `GCF_000218625.1_ASM21862v1` | `not_available` | 46 | 24 | 0.6085 | 1 / 3 | 421 | 0 / 4 / 42 |
| `GCF_000219725.1_ASM21972v1` | `not_available` | 53 | 23 | 0.6145 | 1 / 0 | 557 | 0 / 3 / 50 |
| `GCF_000225345.1_ASM22534v1` | `not_available` | 54 | 19 | 0.5741 | 0 / 1 | 512 | 0 / 4 / 50 |
| `GCF_000226295.1_ASM22629v1` | `not_available` | 48 | 21 | 0.5680 | 0 / 2 | 658 | 0 / 4 / 44 |
| `GCF_000226315.1_ASM22631v1` | `not_available` | 53 | 27 | 0.6292 | 0 / 4 | 504 | 0 / 5 / 48 |
| `GCF_000235405.2_ASM23540v3` | `not_available` | 38 | 11 | 0.4932 | 0 / 2 | 391 | 0 / 4 / 34 |
| `GCF_000237205.1_ASM23720v1` | `not_available` | 39 | 13 | 0.5098 | 0 / 2 | 531 | 0 / 4 / 35 |
| `GCF_000239475.1_ASM23947v1` | `not_available` | 29 | 9 | 0.4971 | 0 / 0 | 183 | 0 / 3 / 26 |
| `GCF_000252445.1_ASM25244v1` | `not_available` | 53 | 24 | 0.6069 | 1 / 3 | 954 | 0 / 5 / 48 |
| `GCF_000253035.1_ASM25303v1` | `not_available` | 46 | 17 | 0.5146 | 1 / 4 | 576 | 0 / 5 / 41 |
| `GCF_000255135.1_ASM25513v1` | `not_available` | 38 | 11 | 0.5407 | 0 / 2 | 383 | 0 / 4 / 34 |
| `GCF_000258405.1_ASM25840v1` | `not_available` | 47 | 25 | 0.6092 | 1 / 3 | 598 | 0 / 6 / 41 |
| `GCF_000259275.1_ASM25927v1` | `not_available` | 34 | 13 | 0.4922 | 1 / 2 | 330 | 0 / 4 / 30 |
| `GCF_000265295.1_ASM26529v1` | `not_available` | 46 | 20 | 0.5607 | 0 / 2 | 488 | 0 / 4 / 42 |
| `GCF_000265425.1_ASM26542v1` | `not_available` | 50 | 21 | 0.5697 | 0 / 2 | 1.022 | 0 / 4 / 46 |
| `GCF_000279145.1_ASM27914v1` | `not_available` | 49 | 23 | 0.6286 | 2 / 3 | 500 | 0 / 5 / 44 |
| `GCF_000284115.1_ASM28411v1` | `not_available` | 44 | 19 | 0.5249 | 0 / 3 | 730 | 0 / 5 / 39 |
| `GCF_000284315.1_ASM28431v1` | `not_available` | 43 | 17 | 0.5808 | 3 / 2 | 472 | 0 / 4 / 39 |
| `GCF_000284335.1_ASM28433v1` | `not_available` | 31 | 9 | 0.4972 | 1 / 0 | 263 | 0 / 4 / 27 |
| `GCF_000309885.1_ASM30988v1` | `not_available` | 49 | 23 | 0.5934 | 1 / 3 | 544 | 0 / 5 / 44 |
| `GCF_000317105.1_ASM31710v1` | `not_available` | 51 | 24 | 0.6179 | 0 / 2 | 1.359 | 0 / 4 / 47 |
| `GCF_000317575.1_ASM31757v1` | `not_available` | 45 | 17 | 0.5593 | 0 / 2 | 1.161 | 0 / 4 / 41 |
| `GCF_000348725.1_ASM34872v1` | `not_available` | 45 | 22 | 0.6220 | 1 / 3 | 599 | 0 / 6 / 39 |
| `GCF_000392435.1_ASM39243v1` | `not_available` | 19 | 3 | 0.4437 | 0 / 1 | 250 | 0 / 3 / 16 |
| `GCF_000439435.1_ASM43943v1` | `not_available` | 20 | 6 | 0.4779 | 0 / 1 | 284 | 0 / 3 / 17 |
| `GCF_000455605.1_GPTSA100-9` | `not_available` | 37 | 15 | 0.5297 | 1 / 3 | 593 | 0 / 5 / 32 |
| `GCF_000471025.2_ASM47102v2` | `not_available` | 30 | 10 | 0.4682 | 0 / 1 | 209 | 0 / 4 / 26 |
| `GCF_000479335.1_ASM47933v1` | `not_available` | 43 | 14 | 0.5735 | 0 / 1 | 272 | 0 / 4 / 39 |
| `GCF_000484535.1_ASM48453v1` | `not_available` | 50 | 19 | 0.5637 | 0 / 3 | 1.094 | 0 / 4 / 46 |
| `GCF_000487995.2_ASM48799v2` | `not_available` | 47 | 21 | 0.5934 | 1 / 2 | 404 | 0 / 4 / 43 |
| `GCF_000504085.1_ASM50408v1` | `not_available` | 36 | 13 | 0.5829 | 0 / 1 | 370 | 0 / 4 / 32 |
| `GCF_000508245.1_ASM50824v1` | `not_available` | 13 | 2 | 0.2545 | 0 / 1 | 252 | 0 / 3 / 10 |
| `GCF_000512735.1_ASM51273v1` | `not_available` | 38 | 15 | 0.5556 | 2 / 2 | 355 | 0 / 4 / 34 |
| `GCF_000513475.1_Bbla` | `not_available` | 20 | 7 | 0.5115 | 0 / 1 | 218 | 0 / 3 / 17 |
| `GCF_000522985.1_ASM52298v1` | `not_available` | 51 | 24 | 0.6039 | 1 / 4 | 1.394 | 0 / 5 / 46 |
| `GCF_000582535.1_ASM58253v1` | `not_available` | 13 | 3 | 0.4364 | 0 / 1 | 138 | 0 / 3 / 10 |
| `GCF_000695095.2_ASM69509v2` | `not_available` | 55 | 27 | 0.6030 | 1 / 4 | 865 | 0 / 6 / 49 |
| `GCF_000724625.1_ASM72462v1` | `not_available` | 48 | 17 | 0.5476 | 0 / 2 | 1.128 | 0 / 4 / 44 |
| `GCF_000767275.3_ASM76727v4` | `not_available` | 39 | 11 | 0.4909 | 0 / 2 | 401 | 0 / 4 / 35 |
| `GCF_000800395.1_ASM80039v1` | `not_available` | 51 | 24 | 0.5896 | 1 / 3 | 751 | 0 / 5 / 46 |
| `GCF_000803625.1_ASM80362v1` | `not_available` | 20 | 3 | 0.4726 | 0 / 1 | 150 | 0 / 4 / 16 |
| `GCF_000828675.1_ASM82867v1` | `not_available` | 38 | 13 | 0.5598 | 0 / 1 | 365 | 0 / 4 / 34 |
| `GCF_000829235.1_ASM82923v1` | `not_available` | 37 | 13 | 0.5459 | 0 / 3 | 450 | 0 / 4 / 33 |
| `GCF_000941075.1_ASM94107v1` | `not_available` | 13 | 2 | 0.3545 | 0 / 1 | 245 | 0 / 3 / 10 |
| `GCF_001017655.1_ASM101765v1` | `not_available` | 42 | 15 | 0.5268 | 0 / 1 | 501 | 0 / 4 / 38 |
| `GCF_001027025.1_ASM102702v1` | `not_available` | 37 | 11 | 0.5071 | 0 / 2 | 387 | 0 / 4 / 33 |
| `GCF_001027285.1_ASM102728v1` | `not_available` | 66 | 33 | 0.6262 | 1 / 3 | 2.372 | 0 / 5 / 61 |
| `GCF_001027545.1_ASM102754v1` | `not_available` | 35 | 12 | 0.5048 | 0 / 2 | 169 | 0 / 4 / 31 |
| `GCF_001273775.1_ASM127377v1` | `not_available` | 52 | 28 | 0.6840 | 1 / 2 | 912 | 0 / 5 / 47 |
| `GCF_001296145.1_ASM129614v1` | `not_available` | 46 | 18 | 0.5776 | 0 / 1 | 383 | 0 / 4 / 42 |
| `GCF_001458695.1_NiCh1` | `not_available` | 51 | 28 | 0.6352 | 1 / 2 | 578 | 0 / 4 / 47 |
| `GCF_001547735.1_ASM154773v1` | `not_available` | 47 | 21 | 0.6229 | 1 / 3 | 312 | 0 / 4 / 43 |
| `GCF_001610955.1_ASM161095v1` | `not_available` | 64 | 27 | 0.6087 | 1 / 4 | 546 | 0 / 6 / 58 |
| `GCF_001611795.1_ASM161179v1` | `not_available` | 56 | 20 | 0.6180 | 0 / 4 | 1.828 | 0 / 5 / 51 |
| `GCF_001618865.1_ASM161886v1` | `not_available` | 61 | 32 | 0.6572 | 1 / 4 | 1.453 | 0 / 6 / 55 |
| `GCF_001676785.2_ASM167678v2` | `not_available` | 48 | 15 | 0.4980 | 0 / 1 | 585 | 0 / 4 / 44 |
| `GCF_001701045.1_ASM170104v1` | `not_available` | 36 | 13 | 0.5882 | 0 / 0 | 600 | 0 / 3 / 33 |
| `GCF_001712875.1_ASM171287v1` | `not_available` | 15 | 5 | 0.4909 | 0 / 0 | 66 | 0 / 3 / 12 |
| `GCF_001747405.1_ASM174740v1` | `not_available` | 61 | 22 | 0.5789 | 2 / 2 | 538 | 0 / 4 / 57 |
| `GCF_001747425.1_ASM174742v1` | `not_available` | 51 | 22 | 0.5929 | 1 / 4 | 1.218 | 0 / 6 / 45 |
| `GCF_001886435.1_ASM188643v1` | `not_available` | 98 | 76 | 0.8733 | 3 / 6 | 874 | 0 / 6 / 92 |
| `GCF_001886815.1_ASM188681v1` | `not_available` | 54 | 26 | 0.6412 | 2 / 4 | 695 | 0 / 6 / 48 |
| `GCF_001889025.1_ASM188902v1` | `not_available` | 51 | 27 | 0.6135 | 0 / 4 | 507 | 0 / 5 / 46 |
| `GCF_001953175.1_ASM195317v1` | `not_available` | 42 | 16 | 0.5322 | 0 / 2 | 450 | 0 / 4 / 38 |
| `GCF_002139875.1_ASM213987v1` | `not_available` | 46 | 20 | 0.5404 | 0 / 3 | 298 | 0 / 4 / 42 |
| `GCF_002157165.1_ASM215716v1` | `not_available` | 96 | 59 | 0.7662 | 1 / 8 | 764 | 0 / 6 / 90 |
| `GCF_002269385.1_ASM226938v1` | `not_available` | 36 | 16 | 0.5474 | 1 / 2 | 931 | 0 / 5 / 31 |
| `GCF_002302565.1_ASM230256v1` | `not_available` | 40 | 18 | 0.5477 | 1 / 2 | 448 | 0 / 5 / 35 |
| `GCF_002305895.1_ASM230589v1` | `not_available` | 70 | 33 | 0.6262 | 1 / 3 | 1.693 | 0 / 5 / 65 |
| `GCF_002310495.1_ASM231049v1` | `not_available` | 51 | 25 | 0.6327 | 0 / 3 | 839 | 0 / 5 / 46 |
| `GCF_002355835.1_ASM235583v1` | `not_available` | 30 | 11 | 0.5225 | 0 / 1 | 142 | 0 / 4 / 26 |
| `GCF_002804025.1_ASM280402v1` | `not_available` | 20 | 7 | 0.4779 | 0 / 1 | 195 | 0 / 3 / 17 |
| `GCF_002813575.1_ASM281357v1` | `not_available` | 57 | 26 | 0.6199 | 0 / 3 | 2.241 | 0 / 4 / 53 |
| `GCF_002872415.1_ASM287241v1` | `not_available` | 49 | 26 | 0.5962 | 1 / 3 | 868 | 0 / 5 / 44 |
| `GCF_003019315.1_ASM301931v1` | `not_available` | 48 | 19 | 0.6058 | 0 / 1 | 432 | 0 / 4 / 44 |
| `GCF_003019655.1_ASM301965v1` | `not_available` | 50 | 19 | 0.6008 | 0 / 1 | 539 | 0 / 4 / 46 |
| `GCF_003019675.1_ASM301967v1` | `not_available` | 50 | 19 | 0.6057 | 0 / 1 | 568 | 0 / 4 / 46 |
| `GCF_003057965.1_ASM305796v1` | `not_available` | 40 | 17 | 0.5430 | 0 / 2 | 279 | 0 / 4 / 36 |
| `GCF_003173015.1_ASM317301v1` | `not_available` | 54 | 23 | 0.6023 | 1 / 3 | 679 | 0 / 5 / 49 |
| `GCF_003258315.1_ASM325831v1` | `not_available` | 52 | 23 | 0.5779 | 1 / 3 | 831 | 0 / 5 / 47 |
| `GCF_003325455.1_ASM332545v1` | `not_available` | 44 | 15 | 0.5591 | 0 / 3 | 657 | 0 / 5 / 39 |
| `GCF_003330725.1_ASM333072v1` | `not_available` | 64 | 23 | 0.5947 | 0 / 4 | 1.538 | 0 / 5 / 59 |
| `GCF_003351905.1_ASM335190v1` | `not_available` | 25 | 8 | 0.5000 | 0 / 1 | 216 | 0 / 4 / 21 |
| `GCF_003353065.1_ASM335306v1` | `not_available` | 48 | 22 | 0.5794 | 1 / 2 | 688 | 0 / 5 / 43 |
| `GCF_003355515.1_ASM335551v1` | `not_available` | 43 | 16 | 0.4979 | 0 / 4 | 395 | 0 / 5 / 38 |
| `GCF_003366055.1_ASM336605v1` | `not_available` | 23 | 5 | 0.5000 | 0 / 1 | 235 | 0 / 4 / 19 |
| `GCF_003491405.1_ASM349140v1` | `not_available` | 46 | 19 | 0.5551 | 1 / 2 | 673 | 0 / 5 / 41 |
| `GCF_003544815.1_ASM354481v1` | `not_available` | 46 | 20 | 0.5330 | 1 / 4 | 514 | 0 / 5 / 41 |
| `GCF_003574315.2_ASM357431v2` | `not_available` | 54 | 24 | 0.6180 | 0 / 4 | 1.397 | 0 / 5 / 49 |
| `GCF_003721155.3_ASM372115v4` | `not_available` | 47 | 17 | 0.5776 | 0 / 3 | 713 | 0 / 6 / 41 |
| `GCF_003815875.1_ASM381587v1` | `not_available` | 46 | 20 | 0.5526 | 1 / 3 | 1.046 | 0 / 5 / 41 |
| `GCF_003966735.1_ASM396673v1` | `not_available` | 56 | 28 | 0.6227 | 1 / 2 | 640 | 0 / 4 / 52 |
| `GCF_003966975.1_ASM396697v1` | `not_available` | 52 | 17 | 0.5573 | 0 / 1 | 619 | 0 / 4 / 48 |
| `GCF_004087915.1_ASM408791v1` | `not_available` | 62 | 27 | 0.6246 | 3 / 3 | 515 | 0 / 4 / 58 |
| `GCF_004102945.1_ASM410294v1` | `not_available` | 37 | 17 | 0.5846 | 0 / 1 | 229 | 0 / 4 / 33 |
| `GCF_004103715.1_ASM410371v1` | `not_available` | 56 | 22 | 0.6097 | 1 / 0 | 526 | 0 / 3 / 53 |
| `GCF_004214795.1_ASM421479v1` | `not_available` | 41 | 21 | 0.5660 | 0 / 4 | 522 | 0 / 5 / 36 |
| `GCF_004376055.1_ASM437605v2` | `not_available` | 50 | 21 | 0.6285 | 1 / 2 | 670 | 0 / 4 / 46 |
| `GCF_004803815.1_ASM480381v1` | `not_available` | 37 | 13 | 0.5431 | 1 / 1 | 358 | 0 / 4 / 33 |
| `GCF_005144435.2_ASM514443v2` | `not_available` | 108 | 61 | 0.7115 | 2 / 4 | 1.160 | 0 / 5 / 103 |
| `GCF_005697565.1_ASM569756v1` | `not_available` | 17 | 3 | 0.5039 | 0 / 1 | 157 | 0 / 4 / 13 |
| `GCF_006352065.1_ASM635206v1` | `not_available` | 48 | 17 | 0.5672 | 0 / 3 | 794 | 0 / 6 / 42 |
| `GCF_007475525.1_ASM747552v1` | `not_available` | 46 | 19 | 0.5504 | 1 / 3 | 442 | 0 / 5 / 41 |
| `GCF_007747655.1_ASM774765v1` | `not_available` | 46 | 19 | 0.5477 | 0 / 3 | 1.632 | 0 / 5 / 41 |
| `GCF_007750655.1_ASM775065v1` | `not_available` | 43 | 16 | 0.4955 | 0 / 2 | 844 | 0 / 6 / 37 |
| `GCF_007904085.1_ASM790408v1` | `not_available` | 47 | 24 | 0.5895 | 0 / 2 | 787 | 0 / 4 / 43 |
| `GCF_007990365.1_ASM799036v1` | `not_available` | 48 | 23 | 0.6178 | 0 / 1 | 460 | 0 / 4 / 44 |
| `GCF_007990505.1_ASM799050v1` | `not_available` | 55 | 26 | 0.6078 | 0 / 1 | 342 | 0 / 4 / 51 |
| `GCF_007993755.1_ASM799375v1` | `not_available` | 64 | 28 | 0.6073 | 1 / 3 | 1.123 | 0 / 5 / 59 |
| `GCF_007993775.1_ASM799377v1` | `not_available` | 70 | 34 | 0.6204 | 1 / 4 | 1.380 | 0 / 6 / 64 |
| `GCF_008086545.1_ASM808654v1` | `not_available` | 21 | 7 | 0.4851 | 0 / 1 | 338 | 0 / 4 / 17 |
| `GCF_008087625.1_ASM808762v1` | `not_available` | 51 | 19 | 0.5620 | 0 / 4 | 1.869 | 0 / 6 / 45 |
| `GCF_008087665.1_ASM808766v1` | `not_available` | 40 | 16 | 0.5268 | 1 / 3 | 1.191 | 0 / 6 / 34 |
| `GCF_009002475.1_ASM900247v1` | `not_available` | 43 | 17 | 0.5409 | 0 / 2 | 1.553 | 0 / 5 / 38 |
| `GCF_009258225.1_ASM925822v1` | `not_available` | 40 | 16 | 0.5476 | 1 / 2 | 484 | 0 / 4 / 36 |
| `GCF_009688965.1_ASM968896v1` | `not_available` | 82 | 28 | 0.6083 | 2 / 5 | 1.311 | 0 / 6 / 76 |
| `GCF_009689585.1_ASM968958v1` | `not_available` | 47 | 19 | 0.5494 | 2 / 3 | 971 | 0 / 5 / 42 |
| `GCF_009731575.1_ASM973157v1` | `not_available` | 50 | 21 | 0.5643 | 1 / 3 | 538 | 0 / 5 / 45 |
| `GCF_009734425.1_ASM973442v1` | `not_available` | 40 | 17 | 0.6269 | 1 / 4 | 637 | 0 / 6 / 34 |
| `GCF_009905215.1_ASM990521v1` | `not_available` | 44 | 19 | 0.5385 | 0 / 1 | 649 | 0 / 4 / 40 |
| `GCF_010726765.1_ASM1072676v1` | `not_available` | 47 | 19 | 0.6123 | 0 / 2 | 1.464 | 0 / 4 / 43 |
| `GCF_011751885.2_ASM1175188v2` | `not_available` | 41 | 13 | 0.5721 | 0 / 3 | 809 | 0 / 5 / 36 |
| `GCF_013004105.1_ASM1300410v1` | `not_available` | 61 | 26 | 0.6048 | 1 / 4 | 935 | 0 / 6 / 55 |
| `GCF_013128195.2_ASM1312819v2` | `not_available` | 52 | 20 | 0.5444 | 0 / 4 | 1.878 | 0 / 6 / 46 |
| `GCF_013177655.1_ASM1317765v1` | `not_available` | 41 | 14 | 0.5067 | 0 / 2 | 347 | 0 / 4 / 37 |
| `GCF_013340845.1_ASM1334084v1` | `not_available` | 80 | 47 | 0.7701 | 2 / 6 | 467 | 0 / 5 / 75 |
| `GCF_014274685.1_ASM1427468v1` | `not_available` | 60 | 23 | 0.5876 | 0 / 4 | 848 | 0 / 4 / 56 |
| `GCF_014337215.1_ASM1433721v1` | `not_available` | 43 | 19 | 0.6029 | 0 / 0 | 564 | 0 / 3 / 40 |
| `GCF_014337255.1_ASM1433725v1` | `not_available` | 46 | 19 | 0.6000 | 0 / 1 | 298 | 0 / 4 / 42 |
| `GCF_014529625.1_ASM1452962v2` | `not_available` | 47 | 21 | 0.5720 | 0 / 3 | 874 | 0 / 6 / 41 |
| `GCF_014672715.2_ASM1467271v2` | `not_available` | 43 | 15 | 0.5256 | 0 / 2 | 423 | 0 / 4 / 39 |
| `GCF_014701235.1_ASM1470123v1` | `not_available` | 50 | 22 | 0.6398 | 2 / 4 | 432 | 0 / 5 / 45 |
| `GCF_015277775.1_ASM1527777v1` | `not_available` | 58 | 27 | 0.5949 | 1 / 3 | 1.742 | 0 / 5 / 53 |
| `GCF_015356815.2_ASM1535681v2` | `not_available` | 37 | 11 | 0.5000 | 1 / 1 | 324 | 0 / 4 / 33 |
| `GCF_015500095.1_ASM1550009v1` | `not_available` | 47 | 18 | 0.5783 | 1 / 4 | 1.061 | 0 / 5 / 42 |
| `GCF_015775515.1_ASM1577551v1` | `not_available` | 38 | 16 | 0.5950 | 0 / 2 | 437 | 0 / 4 / 34 |
| `GCF_016027375.1_ASM1602737v1` | `not_available` | 50 | 26 | 0.6500 | 0 / 1 | 451 | 0 / 4 / 46 |
| `GCF_016446285.1_ASM1644628v1` | `not_available` | 58 | 28 | 0.6439 | 1 / 4 | 544 | 0 / 4 / 54 |
| `GCF_016724785.1_ASM1672478v1` | `not_available` | 45 | 19 | 0.5911 | 0 / 1 | 388 | 0 / 4 / 41 |
| `GCF_017498545.1_ASM1749854v1` | `not_available` | 59 | 30 | 0.6450 | 1 / 4 | 1.592 | 0 / 6 / 53 |
| `GCF_017893965.1_ASM1789396v1` | `not_available` | 44 | 17 | 0.5551 | 0 / 0 | 489 | 0 / 3 / 41 |
| `GCF_018069705.1_ASM1806970v1` | `not_available` | 41 | 16 | 0.5973 | 0 / 1 | 453 | 0 / 3 / 38 |
| `GCF_018127985.1_ASM1812798v1` | `not_available` | 34 | 11 | 0.4826 | 0 / 2 | 551 | 0 / 4 / 30 |
| `GCF_018137965.1_ASM1813796v1` | `not_available` | 62 | 31 | 0.6489 | 0 / 3 | 937 | 0 / 6 / 56 |
| `GCF_018243235.1_ASM1824323v1` | `not_available` | 98 | 58 | 0.7623 | 3 / 7 | 927 | 0 / 6 / 92 |
| `GCF_018292205.1_ASM1829220v1` | `not_available` | 43 | 16 | 0.5223 | 0 / 2 | 943 | 0 / 4 / 39 |
| `GCF_018304765.1_ASM1830476v1` | `not_available` | 49 | 21 | 0.5731 | 0 / 4 | 655 | 0 / 4 / 45 |
| `GCF_018343815.1_ASM1834381v1` | `not_available` | 34 | 12 | 0.4868 | 0 / 1 | 219 | 0 / 4 / 30 |
| `GCF_018398915.1_ASM1839891v1` | `not_available` | 43 | 17 | 0.6018 | 0 / 1 | 448 | 0 / 3 / 40 |
| `GCF_018771645.1_ASM1877164v1` | `not_available` | 53 | 22 | 0.5758 | 1 / 3 | 544 | 0 / 5 / 48 |
| `GCF_018866205.1_ASM1886620v1` | `not_available` | 20 | 2 | 0.3904 | 0 / 1 | 223 | 0 / 4 / 16 |
| `GCF_019048065.1_ASM1904806v1` | `not_available` | 237 | 233 | 0.9912 | 7 / 6 | 892 | 0 / 6 / 231 |
| `GCF_019048105.1_ASM1904810v1` | `not_available` | 174 | 148 | 0.9197 | 7 / 5 | 846 | 0 / 6 / 168 |
| `GCF_019056495.1_ASM1905649v1` | `not_available` | 34 | 12 | 0.4815 | 0 / 1 | 219 | 0 / 4 / 30 |
| `GCF_019063885.1_ASM1906388v1` | `not_available` | 69 | 37 | 0.6877 | 3 / 5 | 773 | 0 / 6 / 63 |
| `GCF_019552005.1_ASM1955200v1` | `not_available` | 45 | 18 | 0.5833 | 0 / 1 | 326 | 0 / 4 / 41 |
| `GCF_020091585.1_ASM2009158v1` | `not_available` | 54 | 20 | 0.5970 | 0 / 3 | 1.977 | 0 / 4 / 50 |
| `GCF_020549965.1_ASM2054996v1` | `not_available` | 29 | 11 | 0.5000 | 0 / 1 | 211 | 0 / 4 / 25 |
| `GCF_020616335.1_ASM2061633v1` | `not_available` | 38 | 17 | 0.6033 | 0 / 2 | 485 | 0 / 6 / 32 |
| `GCF_021057185.1_ASM2105718v1` | `not_available` | 0 | 0 | 0.0000 | 0 / 0 | 293 | 0 / 0 / 0 |
| `GCF_021117425.1_ASM2111742v1` | `not_available` | 37 | 11 | 0.5528 | 0 / 2 | 873 | 0 / 4 / 33 |
| `GCF_021379005.1_ASM2137900v1` | `not_available` | 46 | 19 | 0.5917 | 0 / 2 | 1.204 | 0 / 4 / 42 |
| `GCF_022026155.1_ASM2202615v1` | `not_available` | 55 | 23 | 0.5884 | 2 / 2 | 353 | 0 / 4 / 51 |
| `GCF_022376835.1_ASM2237683v1` | `not_available` | 42 | 19 | 0.5869 | 1 / 0 | 680 | 0 / 3 / 39 |
| `GCF_022813675.1_ASM2281367v1` | `not_available` | 50 | 24 | 0.6061 | 1 / 3 | 373 | 0 / 4 / 46 |
| `GCF_022814725.1_ASM2281472v1` | `not_available` | 46 | 17 | 0.6000 | 1 / 3 | 747 | 0 / 6 / 40 |
| `GCF_022846135.1_ASM2284613v1` | `not_available` | 40 | 17 | 0.6000 | 0 / 0 | 431 | 0 / 3 / 37 |
| `GCF_022869165.1_ASM2286916v1` | `not_available` | 50 | 19 | 0.5691 | 0 / 1 | 646 | 0 / 4 / 46 |
| `GCF_022870745.1_ASM2287074v1` | `not_available` | 50 | 21 | 0.5945 | 1 / 2 | 1.000 | 0 / 4 / 46 |
| `GCF_023035795.1_ASM2303579v1` | `not_available` | 32 | 10 | 0.4865 | 0 / 0 | 308 | 0 / 3 / 29 |
| `GCF_023035875.1_ASM2303587v1` | `not_available` | 32 | 10 | 0.4811 | 0 / 0 | 346 | 0 / 3 / 29 |
| `GCF_023169805.1_ASM2316980v1` | `not_available` | 52 | 24 | 0.5977 | 0 / 3 | 1.473 | 0 / 6 / 46 |
| `GCF_023921165.1_ASM2392116v1` | `not_available` | 39 | 18 | 0.5507 | 0 / 3 | 489 | 0 / 5 / 34 |
| `GCF_023935825.1_ASM2393582v1` | `not_available` | 68 | 40 | 0.6711 | 1 / 4 | 848 | 0 / 5 / 63 |
| `GCF_024205945.1_ASM2420594v1` | `not_available` | 44 | 16 | 0.5936 | 0 / 3 | 385 | 0 / 6 / 38 |
| `GCF_024338845.1_ASM2433884v1` | `not_available` | 16 | 4 | 0.4508 | 0 / 1 | 204 | 0 / 3 / 13 |
| `GCF_024584745.1_ASM2458474v1` | `not_available` | 52 | 22 | 0.5914 | 0 / 3 | 1.193 | 0 / 4 / 48 |
| `GCF_025148635.1_ASM2514863v1` | `not_available` | 42 | 17 | 0.5291 | 0 / 2 | 849 | 0 / 4 / 38 |
| `GCF_025273655.1_ASM2527365v1` | `not_available` | 20 | 3 | 0.4658 | 0 / 2 | 251 | 0 / 3 / 17 |
| `GCF_025558325.1_ASM2555832v1` | `not_available` | 51 | 18 | 0.5935 | 0 / 2 | 381 | 0 / 6 / 45 |
| `GCF_025725725.1_ASM2572572v1` | `not_available` | 37 | 12 | 0.5820 | 0 / 4 | 643 | 0 / 5 / 32 |
| `GCF_025997855.1_ASM2599785v1` | `not_available` | 48 | 19 | 0.5837 | 1 / 3 | 1.056 | 0 / 5 / 43 |
| `GCF_025997995.1_ASM2599799v1` | `not_available` | 45 | 23 | 0.5570 | 1 / 3 | 584 | 0 / 5 / 40 |
| `GCF_026184435.1_ASM2618443v1` | `not_available` | 34 | 5 | 0.4223 | 0 / 2 | 276 | 0 / 4 / 30 |
| `GCF_026240635.1_ASM2624063v1` | `not_available` | 42 | 19 | 0.5734 | 0 / 2 | 879 | 0 / 4 / 38 |
| `GCF_026965555.1_ASM2696555v1` | `not_available` | 57 | 26 | 0.6180 | 1 / 3 | 2.537 | 0 / 5 / 52 |
| `GCF_027594505.1_ASM2759450v1` | `not_available` | 41 | 21 | 0.5759 | 1 / 2 | 648 | 0 / 5 / 36 |
| `GCF_027625475.1_ASM2762547v1` | `not_available` | 72 | 36 | 0.6759 | 2 / 5 | 940 | 0 / 5 / 67 |
| `GCF_027923555.1_ASM2792355v1` | `not_available` | 47 | 19 | 0.5934 | 1 / 3 | 778 | 0 / 6 / 41 |
| `GCF_028622095.1_ASM2862209v1` | `not_available` | 51 | 23 | 0.6129 | 1 / 3 | 815 | 0 / 5 / 46 |
| `GCF_028728065.1_ASM2872806v1` | `not_available` | 51 | 25 | 0.5983 | 0 / 2 | 903 | 0 / 5 / 46 |
| `GCF_028885455.1_ASM2888545v1` | `not_available` | 105 | 79 | 0.8677 | 4 / 6 | 982 | 0 / 6 / 99 |
| `GCF_029023725.1_ASM2902372v1` | `not_available` | 47 | 22 | 0.6379 | 0 / 3 | 1.003 | 0 / 5 / 42 |
| `GCF_029677705.1_ASM2967770v1` | `not_available` | 21 | 6 | 0.4286 | 0 / 1 | 234 | 0 / 3 / 18 |
| `GCF_029919255.1_ASM2991925v1` | `not_available` | 40 | 15 | 0.5479 | 0 / 3 | 503 | 0 / 4 / 36 |
| `GCF_030169785.1_ASM3016978v1` | `not_available` | 50 | 25 | 0.6312 | 1 / 3 | 793 | 0 / 6 / 44 |
| `GCF_030253535.1_ASM3025353v1` | `not_available` | 20 | 3 | 0.4452 | 0 / 1 | 178 | 0 / 4 / 16 |
| `GCF_030295385.1_ASM3029538v1` | `not_available` | 48 | 20 | 0.5926 | 0 / 3 | 648 | 0 / 4 / 44 |
| `GCF_030389965.1_ASM3038996v1` | `not_available` | 61 | 28 | 0.6162 | 2 / 2 | 1.100 | 0 / 5 / 56 |
| `GCF_030438115.2_ASM3043811v2` | `not_available` | 34 | 11 | 0.4896 | 0 / 0 | 280 | 0 / 3 / 31 |
| `GCF_030439285.2_ASM3043928v2` | `not_available` | 33 | 10 | 0.4842 | 0 / 0 | 244 | 0 / 3 / 30 |
| `GCF_030444335.1_ASM3044433v1` | `not_available` | 39 | 17 | 0.5528 | 1 / 2 | 525 | 0 / 5 / 34 |
| `GCF_030503875.2_ASM3050387v2` | `not_available` | 43 | 21 | 0.5640 | 1 / 3 | 751 | 0 / 5 / 38 |
| `GCF_030520105.1_ASM3052010v1` | `not_available` | 40 | 19 | 0.5668 | 0 / 2 | 766 | 0 / 4 / 36 |
| `GCF_030585795.1_ASM3058579v1` | `not_available` | 42 | 13 | 0.5519 | 0 / 1 | 286 | 0 / 3 / 39 |
| `GCF_030623445.1_ASM3062344v1` | `not_available` | 48 | 24 | 0.6160 | 1 / 4 | 594 | 0 / 5 / 43 |
| `GCF_030866825.1_ASM3086682v1` | `not_available` | 43 | 13 | 0.5660 | 0 / 5 | 972 | 0 / 6 / 37 |
| `GCF_032050975.1_ASM3205097v1` | `not_available` | 56 | 23 | 0.5889 | 1 / 2 | 779 | 0 / 4 / 52 |
| `GCF_032051675.1_ASM3205167v1` | `not_available` | 55 | 23 | 0.6008 | 1 / 2 | 832 | 0 / 4 / 51 |
| `GCF_032102965.1_ASM3210296v1` | `not_available` | 48 | 16 | 0.5546 | 0 / 4 | 614 | 0 / 5 / 43 |
| `GCF_032595795.2_ASM3259579v2` | `not_available` | 34 | 10 | 0.4948 | 0 / 0 | 222 | 0 / 3 / 31 |
| `GCF_032698475.1_ASM3269847v1` | `not_available` | 97 | 51 | 0.7400 | 3 / 5 | 972 | 0 / 6 / 91 |
| `GCF_033547035.1_ASM3354703v1` | `not_available` | 66 | 33 | 0.6230 | 1 / 3 | 2.429 | 0 / 5 / 61 |
| `GCF_033878955.1_ASM3387895v1` | `not_available` | 45 | 24 | 0.6063 | 0 / 3 | 977 | 0 / 5 / 40 |
| `GCF_033969665.2_ASM3396966v2` | `not_available` | 33 | 10 | 0.4947 | 0 / 0 | 247 | 0 / 3 / 30 |
| `GCF_033971325.1_ASM3397132v1` | `not_available` | 45 | 18 | 0.5240 | 1 / 2 | 417 | 0 / 4 / 41 |
| `GCF_034259205.1_ASM3425920v1` | `not_available` | 48 | 22 | 0.5827 | 1 / 3 | 956 | 0 / 5 / 43 |
| `GCF_034424215.1_ASM3442421v1` | `not_available` | 48 | 21 | 0.5560 | 1 / 3 | 1.039 | 0 / 5 / 43 |
| `GCF_037076345.1_ASM3707634v1` | `not_available` | 46 | 21 | 0.5498 | 0 / 3 | 1.048 | 0 / 5 / 41 |
| `GCF_037076455.1_ASM3707645v1` | `not_available` | 58 | 26 | 0.6642 | 1 / 2 | 755 | 0 / 4 / 54 |
| `GCF_037081835.1_ASM3708183v1` | `not_available` | 46 | 22 | 0.6196 | 2 / 2 | 579 | 0 / 5 / 41 |
| `GCF_037157775.1_ASM3715777v1` | `not_available` | 73 | 31 | 0.6055 | 1 / 2 | 2.532 | 0 / 6 / 67 |
| `GCF_037477945.1_MHOCP` | `not_available` | 54 | 25 | 0.6218 | 0 / 2 | 316 | 0 / 4 / 50 |
| `GCF_038095605.1_ASM3809560v1` | `not_available` | 48 | 22 | 0.6197 | 0 / 2 | 786 | 0 / 6 / 42 |
| `GCF_038396615.1_ASM3839661v1` | `not_available` | 42 | 19 | 0.6081 | 0 / 2 | 316 | 0 / 5 / 37 |
| `GCF_038400315.1_ASM3840031v1` | `not_available` | 41 | 17 | 0.5681 | 0 / 4 | 1.333 | 0 / 5 / 36 |
| `GCF_038443385.1_ASM3844338v1` | `not_available` | 41 | 19 | 0.5829 | 1 / 3 | 1.657 | 0 / 6 / 35 |
| `GCF_038725785.1_ASM3872578v1` | `not_available` | 51 | 26 | 0.6275 | 0 / 3 | 1.215 | 0 / 5 / 46 |
| `GCF_039555295.1_ASM3955529v1` | `not_available` | 32 | 14 | 0.5108 | 0 / 2 | 400 | 0 / 4 / 28 |
| `GCF_039645985.1_ASM3964598v1` | `not_available` | 45 | 20 | 0.5689 | 1 / 2 | 548 | 0 / 4 / 41 |
| `GCF_039994875.1_ASM3999487v1` | `not_available` | 54 | 27 | 0.6245 | 3 / 2 | 1.167 | 0 / 4 / 50 |
| `GCF_040256395.1_GCC_DsonKR87_v1` | `not_available` | 62 | 24 | 0.5918 | 1 / 3 | 999 | 0 / 5 / 57 |
| `GCF_040450955.1_ASM4045095v1` | `not_available` | 50 | 21 | 0.6230 | 0 / 3 | 338 | 0 / 4 / 46 |
| `GCF_040451425.1_ASM4045142v1` | `not_available` | 49 | 19 | 0.6151 | 0 / 3 | 315 | 0 / 4 / 45 |
| `GCF_040564635.1_ASM4056463v1` | `not_available` | 51 | 27 | 0.6179 | 0 / 2 | 1.324 | 0 / 4 / 47 |
| `GCF_040833985.2_ASM4083398v1` | `not_available` | 50 | 22 | 0.5785 | 1 / 3 | 845 | 0 / 5 / 45 |
| `GCF_041026625.1_ASM4102662v1` | `not_available` | 14 | 3 | 0.3304 | 0 / 1 | 197 | 0 / 3 / 11 |
| `GCF_041154365.1_ASM4115436v1` | `not_available` | 38 | 16 | 0.5792 | 1 / 1 | 421 | 0 / 3 / 35 |
| `GCF_041200005.1_ASM4120000v1` | `not_available` | 48 | 24 | 0.6211 | 0 / 1 | 374 | 0 / 4 / 44 |
| `GCF_041505935.1_Fervidobacterium_thailandense_FC2004` | `not_available` | 35 | 10 | 0.5222 | 0 / 2 | 378 | 0 / 4 / 31 |
| `GCF_042033325.1_ASM4203332v1` | `not_available` | 39 | 15 | 0.5441 | 0 / 1 | 167 | 0 / 4 / 35 |
| `GCF_042847215.1_ASM4284721v1` | `not_available` | 85 | 59 | 0.8084 | 1 / 5 | 693 | 0 / 6 / 79 |
| `GCF_044336245.1_ASM4433624v1` | `not_available` | 45 | 15 | 0.5353 | 0 / 4 | 1.231 | 0 / 6 / 39 |
| `GCF_044789115.1_ASM4478911v1` | `not_available` | 45 | 24 | 0.6137 | 0 / 3 | 646 | 0 / 5 / 40 |
| `GCF_045784255.1_ASM4578425v1` | `not_available` | 61 | 26 | 0.6042 | 1 / 3 | 1.015 | 0 / 5 / 56 |
| `GCF_045784305.1_DRML_24019` | `not_available` | 56 | 22 | 0.5818 | 0 / 3 | 1.245 | 0 / 5 / 51 |
| `GCF_045864955.1_ASM4586495v1` | `not_available` | 43 | 19 | 0.5983 | 0 / 2 | 327 | 0 / 5 / 38 |
| `GCF_046267995.1_ASM4626799v1` | `not_available` | 32 | 11 | 0.5691 | 0 / 1 | 272 | 0 / 4 / 28 |
| `GCF_046774955.1_ASM4677495v1` | `not_available` | 63 | 29 | 0.6117 | 1 / 3 | 982 | 0 / 5 / 58 |
| `GCF_046845665.1_ASM4684566v1` | `not_available` | 45 | 20 | 0.5708 | 1 / 2 | 767 | 0 / 5 / 40 |
| `GCF_046845755.1_ASM4684575v1` | `not_available` | 42 | 19 | 0.5759 | 1 / 2 | 744 | 0 / 5 / 37 |
| `GCF_047159735.1_ASM4715973v1` | `not_available` | 21 | 6 | 0.4306 | 0 / 1 | 202 | 0 / 4 / 17 |
| `GCF_048565275.1_ASM4856527v1` | `not_available` | 82 | 31 | 0.6074 | 2 / 3 | 1.667 | 0 / 6 / 76 |
| `GCF_049669535.1_sGff_cult_auto_v1` | `not_available` | 25 | 8 | 0.5033 | 0 / 1 | 417 | 0 / 5 / 20 |
| `GCF_049676565.1_Dip510` | `not_available` | 31 | 9 | 0.5057 | 0 / 1 | 427 | 0 / 4 / 27 |
| `GCF_050039385.1_ASM5003938v1` | `not_available` | 68 | 36 | 0.6307 | 1 / 3 | 2.037 | 0 / 5 / 63 |
| `GCF_050267035.1_ASM5026703v1` | `not_available` | 40 | 18 | 0.6019 | 0 / 1 | 356 | 0 / 4 / 36 |
| `GCF_050389785.1_ASM5038978v1` | `not_available` | 28 | 9 | 0.5092 | 0 / 1 | 333 | 0 / 4 / 24 |
| `GCF_050500235.1_ASM5050023v1` | `not_available` | 47 | 21 | 0.6000 | 0 / 2 | 1.136 | 0 / 4 / 43 |
| `GCF_050563455.1_ASM5056345v1` | `not_available` | 50 | 21 | 0.5732 | 0 / 3 | 1.608 | 0 / 4 / 46 |
| `GCF_050780155.1_QI0013` | `not_available` | 70 | 33 | 0.6540 | 2 / 1 | 751 | 0 / 4 / 66 |
| `GCF_051175385.1_ASM5117538v1` | `not_available` | 46 | 24 | 0.6099 | 0 / 4 | 1.414 | 0 / 5 / 41 |
| `GCF_052161605.1_ASM5216160v1` | `not_available` | 48 | 18 | 0.5708 | 0 / 4 | 1.346 | 0 / 6 / 42 |
| `GCF_054054625.1_ASM5405462v1` | `not_available` | 52 | 24 | 0.6387 | 1 / 3 | 719 | 0 / 5 / 47 |
| `GCF_054917915.1_ASM5491791v1` | `not_available` | 70 | 36 | 0.6677 | 1 / 4 | 856 | 0 / 5 / 65 |
| `GCF_055378285.1_ASM5537828v1` | `not_available` | 55 | 26 | 0.6453 | 1 / 2 | 514 | 0 / 4 / 51 |
| `GCF_055378975.1_ASM5537897v1` | `not_available` | 64 | 31 | 0.6633 | 1 / 2 | 619 | 0 / 4 / 60 |
| `GCF_055379415.1_ASM5537941v1` | `not_available` | 60 | 23 | 0.5671 | 1 / 4 | 1.056 | 0 / 5 / 55 |
| `GCF_055379825.1_ASM5537982v1` | `not_available` | 44 | 19 | 0.5625 | 1 / 2 | 648 | 0 / 5 / 39 |
| `GCF_055379935.1_ASM5537993v1` | `not_available` | 45 | 18 | 0.5595 | 1 / 2 | 691 | 0 / 6 / 39 |
| `GCF_055380895.1_ASM5538089v1` | `not_available` | 44 | 22 | 0.6009 | 0 / 2 | 293 | 0 / 4 / 40 |
| `GCF_055383505.1_ASM5538350v1` | `not_available` | 43 | 19 | 0.5860 | 0 / 1 | 372 | 0 / 4 / 39 |
| `GCF_055384375.1_ASM5538437v1` | `not_available` | 48 | 18 | 0.5490 | 1 / 3 | 399 | 0 / 4 / 44 |
| `GCF_055384745.1_ASM5538474v1` | `not_available` | 38 | 20 | 0.5928 | 0 / 2 | 327 | 0 / 4 / 34 |
| `GCF_055385515.1_ASM5538551v1` | `not_available` | 43 | 17 | 0.5333 | 2 / 2 | 446 | 0 / 4 / 39 |
| `GCF_055387075.1_ASM5538707v1` | `not_available` | 44 | 13 | 0.5316 | 0 / 2 | 395 | 0 / 4 / 40 |
| `GCF_055387105.1_ASM5538710v1` | `not_available` | 32 | 14 | 0.5632 | 0 / 1 | 393 | 0 / 4 / 28 |
| `GCF_055387525.1_ASM5538752v1` | `not_available` | 59 | 21 | 0.5882 | 1 / 2 | 860 | 0 / 4 / 55 |
| `GCF_055388355.1_ASM5538835v1` | `not_available` | 50 | 22 | 0.5896 | 0 / 3 | 1.184 | 0 / 5 / 45 |
| `GCF_055390195.1_ASM5539019v1` | `not_available` | 44 | 17 | 0.5263 | 1 / 2 | 407 | 0 / 4 / 40 |
| `GCF_055390205.1_ASM5539020v1` | `not_available` | 34 | 11 | 0.5000 | 0 / 2 | 523 | 0 / 4 / 30 |
| `GCF_055390565.1_ASM5539056v1` | `not_available` | 44 | 14 | 0.5636 | 0 / 0 | 367 | 0 / 3 / 41 |
| `GCF_055390845.1_ASM5539084v1` | `not_available` | 51 | 20 | 0.6240 | 0 / 3 | 378 | 0 / 4 / 47 |
| `GCF_055390895.1_ASM5539089v1` | `not_available` | 56 | 19 | 0.5484 | 0 / 1 | 574 | 0 / 3 / 53 |
| `GCF_055391325.1_ASM5539132v1` | `not_available` | 38 | 12 | 0.5024 | 0 / 1 | 338 | 0 / 4 / 34 |
| `GCF_055396605.1_ASM5539660v1` | `not_available` | 52 | 22 | 0.6349 | 0 / 1 | 436 | 0 / 4 / 48 |
| `GCF_055694045.1_ASM5569404v1` | `not_available` | 45 | 18 | 0.5808 | 0 / 1 | 374 | 0 / 4 / 41 |
| `GCF_055763695.1_ASM5576369v1` | `not_available` | 45 | 17 | 0.5721 | 0 / 2 | 1.241 | 0 / 6 / 39 |
| `GCF_056388415.1_ASM5638841v1` | `not_available` | 14 | 3 | 0.3482 | 0 / 1 | 263 | 0 / 3 / 11 |
| `GCF_058362665.1_ASM5836266v1` | `not_available` | 72 | 31 | 0.5994 | 1 / 3 | 1.780 | 0 / 5 / 67 |
| `GCF_059632245.1_ASM5963224v1` | `not_available` | 50 | 24 | 0.5878 | 1 / 3 | 1.137 | 0 / 5 / 45 |
| `GCF_059696725.1_ASM5969672v1` | `not_available` | 55 | 23 | 0.5993 | 1 / 3 | 911 | 0 / 6 / 49 |
| `GCF_059705355.1_ASM5970535v1` | `not_available` | 37 | 13 | 0.5220 | 0 / 2 | 483 | 0 / 4 / 33 |
| `GCF_059726255.1_ASM5972625v1` | `not_available` | 63 | 34 | 0.6357 | 0 / 4 | 931 | 0 / 6 / 57 |
| `GCF_059726345.1_ASM5972634v1` | `not_available` | 63 | 30 | 0.6053 | 0 / 3 | 1.636 | 0 / 4 / 59 |
| `GCF_060026435.1_ASM6002643v1` | `not_available` | 68 | 36 | 0.6472 | 1 / 3 | 2.557 | 0 / 5 / 63 |
| `GCF_060027055.1_ASM6002705v1` | `not_available` | 72 | 31 | 0.6159 | 1 / 2 | 2.030 | 0 / 5 / 67 |
| `GCF_060030665.1_ASM6003066v1` | `not_available` | 74 | 35 | 0.6277 | 1 / 3 | 2.750 | 0 / 5 / 69 |
| `GCF_060030705.1_ASM6003070v1` | `not_available` | 58 | 28 | 0.6093 | 1 / 3 | 1.464 | 0 / 5 / 53 |
| `GCF_060030825.1_ASM6003082v1` | `not_available` | 65 | 33 | 0.6367 | 1 / 3 | 2.132 | 0 / 5 / 60 |
| `GCF_060031025.1_ASM6003102v1` | `not_available` | 55 | 25 | 0.6066 | 1 / 3 | 2.100 | 0 / 5 / 50 |
| `GCF_060138195.1_ASM6013819v1` | `not_available` | 57 | 25 | 0.6007 | 1 / 3 | 2.096 | 0 / 5 / 52 |
| `GCF_900169565.1_NSJP_Ch1` | `not_available` | 46 | 21 | 0.6336 | 1 / 2 | 784 | 0 / 5 / 41 |
| `GCF_900184705.1_CAMBI-1` | `not_available` | 35 | 9 | 0.4975 | 1 / 1 | 436 | 0 / 4 / 31 |
| `GCF_900217235.1_DPRO_PRJEB22548_v1` | `not_available` | 67 | 29 | 0.6000 | 1 / 2 | 866 | 0 / 4 / 63 |
| `GCF_900239975.1_Chlamydia_sp._S15-834K` | `not_available` | 33 | 12 | 0.4839 | 0 / 1 | 220 | 0 / 4 / 29 |
| `GCF_900475595.1_44343_G01` | `not_available` | 48 | 14 | 0.5696 | 0 / 1 | 284 | 0 / 4 / 44 |
| `GCF_900660485.1_50648_E01-3` | `not_available` | 16 | 4 | 0.3644 | 0 / 1 | 283 | 0 / 3 / 13 |
| `GCF_900660665.1_57612_D01-3` | `not_available` | 13 | 2 | 0.3091 | 0 / 1 | 237 | 0 / 3 / 10 |
| `GCF_901538265.1_ASM90153826v1` | `not_available` | 48 | 18 | 0.5605 | 0 / 4 | 1.937 | 0 / 6 / 42 |
| `GCF_902806995.1_15-2067_O50` | `not_available` | 37 | 13 | 0.5149 | 0 / 1 | 221 | 0 / 4 / 33 |
| `GCF_907165195.1_Cryomorphaceae_CECT30217` | `not_available` | 41 | 18 | 0.5644 | 1 / 2 | 817 | 0 / 5 / 36 |
| `GCF_918320885.1_Harlequin_endosymbiont1` | `not_available` | 84 | 66 | 0.8761 | 3 / 2 | 363 | 0 / 6 / 78 |
| `GCF_946900835.1_NSPWAT` | `not_available` | 55 | 23 | 0.6023 | 1 / 3 | 587 | 0 / 5 / 50 |
| `GCF_947241125.1_NsDNF` | `not_available` | 63 | 27 | 0.5908 | 1 / 2 | 922 | 0 / 4 / 59 |
| `GCF_949774925.1_MFUM-new` | `not_available` | 48 | 21 | 0.5556 | 1 / 3 | 470 | 0 / 5 / 43 |
| `GCF_963264975.1_NZ_220527_assembly_043` | `not_available` | 48 | 20 | 0.5940 | 0 / 1 | 430 | 0 / 4 / 44 |
| `GCF_964030865.1_idDioLine1.Spiroplasma_sp_1.1` | `not_available` | 18 | 5 | 0.4609 | 0 / 1 | 275 | 0 / 3 / 15 |
| `GCF_964030885.1_icBarPell3.Wolbachia_sp_1.1` | `not_available` | 40 | 14 | 0.4836 | 0 / 4 | 385 | 0 / 5 / 35 |
| `GCF_964642095.1_Ff3` | `not_available` | 45 | 19 | 0.6304 | 0 / 0 | 436 | 0 / 3 / 42 |
| `GCF_965152565.1_NS_bockiana` | `not_available` | 51 | 21 | 0.6220 | 1 / 2 | 860 | 0 / 4 / 47 |
| `GCF_976986435.1_NS_calida2` | `not_available` | 58 | 23 | 0.5906 | 1 / 2 | 682 | 0 / 4 / 54 |
| `GCF_978936825.1_DTS69` | `not_available` | 48 | 21 | 0.5840 | 2 / 4 | 429 | 0 / 4 / 44 |

## Regenerar

```bash
python scripts/complexes/transfer_coverage_audit.py
```

El detalle completo de limitaciones y firmas está en `reports/complex_transfer_coverage.tsv`.
