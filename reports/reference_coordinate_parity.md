# Paridad por coordenadas: sistema vs GFF de referencia

**Generado:** `2026-11-25T19:46:46Z`  
**Manifiesto:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/reference_genomes.json`  
**Raíz de datos:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`  
**Tolerancia de bordes:** 5 pb  
**Solape recíproco mínimo:** 0.50

## Alcance

Compara los features **CDS, tRNA y rRNA** del GFF del sistema (`results/<id>/<id>.gff3`) contra el GFF de referencia (`annotation_path`) para cada registro con status `annotated` o `available` del manifiesto.

El emparejamiento es **1:1** y exige mismo contig, mismo tipo y misma hebra. Dos features son candidatos si se solapan recíprocamente o si sus dos bordes difieren como mucho **5 pb**. Si falta cualquiera de los dos GFF, el genoma queda como `not_available`; no se imputan coordenadas ni se completa con otros archivos.

Esta versión **no evalúa función, producto, identidad de secuencia ni calidad biológica** de las predicciones. Tampoco modifica `results/`, los GFF ni el manifiesto.

## Cobertura

- Registros con status elegible: **41**
- Comparables (ambos GFF legibles): **40** (ok **0**, partial **40**)
- No disponibles: **1**
- Omitidos por status: **1**

### Genomas no disponibles

| genoma | status | motivo |
|---|---|---|
| `ecoli_k12` | available | falta el GFF del sistema: results/ecoli_k12/ecoli_k12.gff3 |

## Resultados por genoma (todas las clases)

TP/FP/FN se calculan sobre CDS + tRNA + rRNA. Precisión, recall y F1 son fracciones (0–1); `% exactos` es un porcentaje.

| genoma | estado | sistema | referencia | TP | FP | FN | precisión | recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `ecoli_k12` | not_available | — | — | — | — | — | — | — | — |
| `subtilis` | partial | 4348 | 4448 | 4271 | 77 | 177 | 0.9823 | 0.9602 | 0.9711 |
| `mtb_h37rv` | partial | 4033 | 3954 | 3790 | 243 | 164 | 0.9397 | 0.9585 | 0.9490 |
| `cglutamicum` | partial | 3113 | 3076 | 2988 | 125 | 88 | 0.9598 | 0.9714 | 0.9656 |
| `pputida` | partial | 5597 | 5568 | 5480 | 117 | 88 | 0.9791 | 0.9842 | 0.9816 |
| `paeruginosa` | partial | 5787 | 5648 | 5616 | 171 | 32 | 0.9705 | 0.9943 | 0.9823 |
| `crescentus` | partial | 3794 | 3943 | 3755 | 39 | 188 | 0.9897 | 0.9523 | 0.9706 |
| `synechocystis` | partial | 3823 | 3745 | 3623 | 200 | 122 | 0.9477 | 0.9674 | 0.9574 |
| `saureus` | partial | 2747 | 2842 | 2644 | 103 | 198 | 0.9625 | 0.9303 | 0.9461 |
| `spneumoniae` | partial | 2235 | 2187 | 2056 | 179 | 131 | 0.9199 | 0.9401 | 0.9299 |
| `btheta` | partial | 5026 | 4886 | 4784 | 242 | 102 | 0.9519 | 0.9791 | 0.9653 |
| `scoelicolor` | partial | 7822 | 7839 | 7671 | 151 | 168 | 0.9807 | 0.9786 | 0.9796 |
| `mxanthus` | partial | 7402 | 7384 | 7232 | 170 | 152 | 0.9770 | 0.9794 | 0.9782 |
| `vcholerae` | partial | 3789 | 3798 | 3701 | 88 | 97 | 0.9768 | 0.9745 | 0.9756 |
| `cdifficile` | partial | 3814 | 3736 | 3693 | 121 | 43 | 0.9683 | 0.9885 | 0.9783 |
| `llactis` | partial | 2686 | 2656 | 2614 | 72 | 42 | 0.9732 | 0.9842 | 0.9787 |
| `lmonocytogenes` | partial | 3000 | 2952 | 2933 | 67 | 19 | 0.9777 | 0.9936 | 0.9856 |
| `efaecalis` | partial | 2791 | 2753 | 2722 | 69 | 31 | 0.9753 | 0.9887 | 0.9820 |
| `dradiodurans` | partial | 3204 | 3207 | 3153 | 51 | 54 | 0.9841 | 0.9832 | 0.9836 |
| `tthermophilus` | partial | 2280 | 2282 | 2225 | 55 | 57 | 0.9759 | 0.9750 | 0.9754 |
| `hpylori` | partial | 1619 | 1626 | 1562 | 57 | 64 | 0.9648 | 0.9606 | 0.9627 |
| `cjejuni` | partial | 1721 | 1624 | 1605 | 116 | 19 | 0.9326 | 0.9883 | 0.9596 |
| `nmeningitidis` | partial | 2119 | 2213 | 2010 | 109 | 203 | 0.9486 | 0.9083 | 0.9280 |
| `hinfluenzae` | partial | 1893 | 1884 | 1847 | 46 | 37 | 0.9757 | 0.9804 | 0.9780 |
| `ctrachomatis` | partial | 940 | 930 | 915 | 25 | 15 | 0.9734 | 0.9839 | 0.9786 |
| `tpallidum` | partial | 1031 | 1038 | 1007 | 24 | 31 | 0.9767 | 0.9701 | 0.9734 |
| `bburgdorferi` | partial | 1569 | 1565 | 1497 | 72 | 68 | 0.9541 | 0.9565 | 0.9553 |
| `mgenitalium` | partial | 532 | 563 | 524 | 8 | 39 | 0.9850 | 0.9307 | 0.9571 |
| `abaumannii` | partial | 3890 | 3828 | 3782 | 108 | 46 | 0.9722 | 0.9880 | 0.9800 |
| `kpneumoniae` | partial | 5640 | 5866 | 5438 | 202 | 428 | 0.9642 | 0.9270 | 0.9452 |
| `styphimurium` | partial | 4722 | 4661 | 4535 | 187 | 126 | 0.9604 | 0.9730 | 0.9667 |
| `afabrum` | partial | 5393 | 5332 | 5228 | 165 | 104 | 0.9694 | 0.9805 | 0.9749 |
| `smeliloti` | partial | 6398 | 6351 | 6110 | 288 | 241 | 0.9550 | 0.9621 | 0.9585 |
| `mgryphiswaldense` | partial | 4002 | 3972 | 3857 | 145 | 115 | 0.9638 | 0.9710 | 0.9674 |
| `soneidensis` | partial | 4694 | 4623 | 4504 | 190 | 119 | 0.9595 | 0.9743 | 0.9668 |
| `cnecator` | partial | 7990 | 7971 | 7685 | 305 | 286 | 0.9618 | 0.9641 | 0.9629 |
| `zmobilis` | partial | 1912 | 1906 | 1864 | 48 | 42 | 0.9749 | 0.9780 | 0.9764 |
| `tmaritima` | partial | 1903 | 1915 | 1882 | 21 | 33 | 0.9890 | 0.9828 | 0.9859 |
| `neuropea` | partial | 2665 | 2675 | 2542 | 123 | 133 | 0.9538 | 0.9503 | 0.9520 |
| `dvulgaris` | partial | 3382 | 3332 | 3258 | 124 | 74 | 0.9633 | 0.9778 | 0.9705 |
| `rpalustris` | partial | 5002 | 5030 | 4922 | 80 | 108 | 0.9840 | 0.9785 | 0.9812 |

### Bordes de los pares emparejados

Las diferencias son **sistema − referencia** en pb (GFF 1-based, inclusivo). `exactos` cuenta pares con inicio y fin idénticos.

| genoma | pares | exactos | % exactos | Δinicio medio | Δinicio mediano | Δfin medio | Δfin mediano | deriva máxima |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `ecoli_k12` | — | — | — | — | — | — | — | — |
| `subtilis` | 4271 | 3763 | 88.11% | -0.24 | +0.00 | -0.09 | +0.00 | 960 |
| `mtb_h37rv` | 3790 | 2855 | 75.33% | +4.85 | +0.00 | -4.46 | +0.00 | 1014 |
| `cglutamicum` | 2988 | 2512 | 84.07% | -1.20 | +0.00 | +0.04 | +0.00 | 1176 |
| `pputida` | 5480 | 5026 | 91.72% | -1.90 | +0.00 | +0.53 | +0.00 | 861 |
| `paeruginosa` | 5616 | 5091 | 90.65% | +0.52 | +0.00 | -0.41 | +0.00 | 639 |
| `crescentus` | 3755 | 2999 | 79.87% | +2.69 | +0.00 | -3.45 | +0.00 | 639 |
| `synechocystis` | 3623 | 3047 | 84.10% | +0.00 | +0.00 | -0.24 | +0.00 | 693 |
| `saureus` | 2644 | 2279 | 86.20% | -2.03 | +0.00 | +2.23 | +0.00 | 189 |
| `spneumoniae` | 2056 | 1833 | 89.15% | +4.40 | +0.00 | -4.83 | +0.00 | 894 |
| `btheta` | 4784 | 4237 | 88.57% | -0.93 | +0.00 | +0.55 | +0.00 | 1449 |
| `scoelicolor` | 7671 | 6524 | 85.05% | -0.12 | +0.00 | +0.87 | +0.00 | 1698 |
| `mxanthus` | 7232 | 5979 | 82.67% | -0.47 | +0.00 | +0.17 | +0.00 | 1134 |
| `vcholerae` | 3701 | 3394 | 91.70% | -0.43 | +0.00 | +0.06 | +0.00 | 586 |
| `cdifficile` | 3693 | 3333 | 90.25% | +0.02 | +0.00 | +0.30 | +0.00 | 691 |
| `llactis` | 2614 | 2366 | 90.51% | -0.43 | +0.00 | -0.24 | +0.00 | 692 |
| `lmonocytogenes` | 2933 | 2634 | 89.81% | -0.75 | +0.00 | +1.00 | +0.00 | 210 |
| `efaecalis` | 2722 | 2356 | 86.55% | -0.07 | +0.00 | +1.12 | +0.00 | 1677 |
| `dradiodurans` | 3153 | 2674 | 84.81% | +0.54 | +0.00 | -0.48 | +0.00 | 746 |
| `tthermophilus` | 2225 | 1946 | 87.46% | +2.29 | +0.00 | -1.14 | +0.00 | 585 |
| `hpylori` | 1562 | 1374 | 87.96% | +2.50 | +0.00 | -2.23 | +0.00 | 1062 |
| `cjejuni` | 1605 | 1488 | 92.71% | +1.17 | +0.00 | -1.48 | +0.00 | 1129 |
| `nmeningitidis` | 2010 | 1737 | 86.42% | +2.57 | +0.00 | -7.18 | +0.00 | 1290 |
| `hinfluenzae` | 1847 | 1705 | 92.31% | +0.69 | +0.00 | -2.10 | +0.00 | 866 |
| `ctrachomatis` | 915 | 811 | 88.63% | -0.35 | +0.00 | -1.04 | +0.00 | 204 |
| `tpallidum` | 1007 | 749 | 74.38% | -0.91 | +0.00 | -0.21 | +0.00 | 699 |
| `bburgdorferi` | 1497 | 1323 | 88.38% | +2.24 | +0.00 | -1.27 | +0.00 | 592 |
| `mgenitalium` | 524 | 467 | 89.12% | +4.02 | +0.00 | -2.25 | +0.00 | 381 |
| `abaumannii` | 3782 | 3426 | 90.59% | +0.09 | +0.00 | +0.42 | +0.00 | 1159 |
| `kpneumoniae` | 5438 | 4756 | 87.46% | -1.94 | +0.00 | +2.02 | +0.00 | 585 |
| `styphimurium` | 4535 | 3904 | 86.09% | +0.36 | +0.00 | -1.17 | +0.00 | 636 |
| `afabrum` | 5228 | 4622 | 88.41% | -1.05 | +0.00 | -0.29 | +0.00 | 1145 |
| `smeliloti` | 6110 | 5439 | 89.02% | +2.79 | +0.00 | -3.20 | +0.00 | 1355 |
| `mgryphiswaldense` | 3857 | 3336 | 86.49% | +0.80 | +0.00 | -1.01 | +0.00 | 672 |
| `soneidensis` | 4504 | 3995 | 88.70% | +3.41 | +0.00 | -2.06 | +0.00 | 765 |
| `cnecator` | 7685 | 6700 | 87.18% | +1.26 | +0.00 | -1.08 | +0.00 | 1008 |
| `zmobilis` | 1864 | 1683 | 90.29% | -0.91 | +0.00 | +0.78 | +0.00 | 858 |
| `tmaritima` | 1882 | 1551 | 82.41% | -0.48 | +0.00 | -0.31 | +0.00 | 1425 |
| `neuropea` | 2542 | 2183 | 85.88% | +0.76 | +0.00 | -2.88 | +0.00 | 1431 |
| `dvulgaris` | 3258 | 2839 | 87.14% | -0.28 | +0.00 | +1.32 | +0.00 | 624 |
| `rpalustris` | 4922 | 4381 | 89.01% | +1.35 | +0.00 | -0.67 | +0.00 | 548 |

## Resultados por clase

| genoma | clase | estado | sistema | referencia | TP | FP | FN | precisión | recall | F1 | exactos |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ecoli_k12` | CDS | not_available | — | — | — | — | — | — | — | — | — |
| `ecoli_k12` | tRNA | not_available | — | — | — | — | — | — | — | — | — |
| `ecoli_k12` | rRNA | not_available | — | — | — | — | — | — | — | — | — |
| `subtilis` | CDS | partial | 4242 | 4332 | 4165 | 77 | 167 | 0.9818 | 0.9614 | 0.9715 | 3743 |
| `subtilis` | tRNA | ok | 86 | 86 | 86 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 20 |
| `subtilis` | rRNA | partial | 20 | 30 | 20 | 0 | 10 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `mtb_h37rv` | CDS | partial | 3986 | 3906 | 3743 | 243 | 163 | 0.9390 | 0.9583 | 0.9486 | 2823 |
| `mtb_h37rv` | tRNA | ok | 45 | 45 | 45 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 32 |
| `mtb_h37rv` | rRNA | partial | 2 | 3 | 2 | 0 | 1 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `cglutamicum` | CDS | partial | 3041 | 2999 | 2917 | 124 | 82 | 0.9592 | 0.9727 | 0.9659 | 2453 |
| `cglutamicum` | tRNA | partial | 60 | 59 | 59 | 1 | 0 | 0.9833 | 1.0000 | 0.9916 | 59 |
| `cglutamicum` | rRNA | partial | 12 | 18 | 12 | 0 | 6 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `pputida` | CDS | partial | 5507 | 5471 | 5390 | 117 | 81 | 0.9788 | 0.9852 | 0.9820 | 4951 |
| `pputida` | tRNA | ok | 75 | 75 | 75 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 75 |
| `pputida` | rRNA | partial | 15 | 22 | 15 | 0 | 7 | 1.0000 | 0.6818 | 0.8108 | 0 |
| `paeruginosa` | CDS | partial | 5715 | 5572 | 5545 | 170 | 27 | 0.9703 | 0.9952 | 0.9826 | 5075 |
| `paeruginosa` | tRNA | partial | 64 | 63 | 63 | 1 | 0 | 0.9844 | 1.0000 | 0.9921 | 16 |
| `paeruginosa` | rRNA | partial | 8 | 13 | 8 | 0 | 5 | 1.0000 | 0.6154 | 0.7619 | 0 |
| `crescentus` | CDS | partial | 3738 | 3886 | 3700 | 38 | 186 | 0.9898 | 0.9521 | 0.9706 | 2948 |
| `crescentus` | tRNA | partial | 52 | 51 | 51 | 1 | 0 | 0.9808 | 1.0000 | 0.9903 | 51 |
| `crescentus` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `synechocystis` | CDS | partial | 3780 | 3698 | 3580 | 200 | 118 | 0.9471 | 0.9681 | 0.9575 | 3006 |
| `synechocystis` | tRNA | ok | 41 | 41 | 41 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 41 |
| `synechocystis` | rRNA | partial | 2 | 6 | 2 | 0 | 4 | 1.0000 | 0.3333 | 0.5000 | 0 |
| `saureus` | CDS | partial | 2675 | 2767 | 2574 | 101 | 193 | 0.9622 | 0.9302 | 0.9459 | 2223 |
| `saureus` | tRNA | partial | 61 | 59 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 56 |
| `saureus` | rRNA | partial | 11 | 16 | 11 | 0 | 5 | 1.0000 | 0.6875 | 0.8148 | 0 |
| `spneumoniae` | CDS | partial | 2168 | 2117 | 1990 | 178 | 127 | 0.9179 | 0.9400 | 0.9288 | 1775 |
| `spneumoniae` | tRNA | partial | 59 | 58 | 58 | 1 | 0 | 0.9831 | 1.0000 | 0.9915 | 58 |
| `spneumoniae` | rRNA | partial | 8 | 12 | 8 | 0 | 4 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `btheta` | CDS | partial | 4945 | 4801 | 4704 | 241 | 97 | 0.9513 | 0.9798 | 0.9653 | 4167 |
| `btheta` | tRNA | partial | 71 | 70 | 70 | 1 | 0 | 0.9859 | 1.0000 | 0.9929 | 70 |
| `btheta` | rRNA | partial | 10 | 15 | 10 | 0 | 5 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `scoelicolor` | CDS | partial | 7738 | 7755 | 7593 | 145 | 162 | 0.9813 | 0.9791 | 0.9802 | 6458 |
| `scoelicolor` | tRNA | partial | 72 | 66 | 66 | 6 | 0 | 0.9167 | 1.0000 | 0.9565 | 66 |
| `scoelicolor` | rRNA | partial | 12 | 18 | 12 | 0 | 6 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `mxanthus` | CDS | partial | 7325 | 7306 | 7158 | 167 | 148 | 0.9772 | 0.9797 | 0.9784 | 5913 |
| `mxanthus` | tRNA | partial | 69 | 66 | 66 | 3 | 0 | 0.9565 | 1.0000 | 0.9778 | 66 |
| `mxanthus` | rRNA | partial | 8 | 12 | 8 | 0 | 4 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `vcholerae` | CDS | partial | 3666 | 3665 | 3578 | 88 | 87 | 0.9760 | 0.9763 | 0.9761 | 3292 |
| `vcholerae` | tRNA | ok | 102 | 102 | 102 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 102 |
| `vcholerae` | rRNA | partial | 21 | 31 | 21 | 0 | 10 | 1.0000 | 0.6774 | 0.8077 | 0 |
| `cdifficile` | CDS | partial | 3700 | 3611 | 3580 | 120 | 31 | 0.9676 | 0.9914 | 0.9794 | 3243 |
| `cdifficile` | tRNA | partial | 91 | 90 | 90 | 1 | 0 | 0.9890 | 1.0000 | 0.9945 | 90 |
| `cdifficile` | rRNA | partial | 23 | 35 | 23 | 0 | 12 | 1.0000 | 0.6571 | 0.7931 | 0 |
| `llactis` | CDS | partial | 2602 | 2567 | 2531 | 71 | 36 | 0.9727 | 0.9860 | 0.9793 | 2296 |
| `llactis` | tRNA | partial | 71 | 70 | 70 | 1 | 0 | 0.9859 | 1.0000 | 0.9929 | 70 |
| `llactis` | rRNA | partial | 13 | 19 | 13 | 0 | 6 | 1.0000 | 0.6842 | 0.8125 | 0 |
| `lmonocytogenes` | CDS | partial | 2921 | 2867 | 2854 | 67 | 13 | 0.9771 | 0.9955 | 0.9862 | 2569 |
| `lmonocytogenes` | tRNA | ok | 67 | 67 | 67 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 65 |
| `lmonocytogenes` | rRNA | partial | 12 | 18 | 12 | 0 | 6 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `efaecalis` | CDS | partial | 2730 | 2690 | 2663 | 67 | 27 | 0.9755 | 0.9900 | 0.9827 | 2305 |
| `efaecalis` | tRNA | partial | 53 | 51 | 51 | 2 | 0 | 0.9623 | 1.0000 | 0.9808 | 51 |
| `efaecalis` | rRNA | partial | 8 | 12 | 8 | 0 | 4 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `dradiodurans` | CDS | partial | 3150 | 3148 | 3100 | 50 | 48 | 0.9841 | 0.9848 | 0.9844 | 2624 |
| `dradiodurans` | tRNA | partial | 51 | 50 | 50 | 1 | 0 | 0.9804 | 1.0000 | 0.9901 | 50 |
| `dradiodurans` | rRNA | partial | 3 | 9 | 3 | 0 | 6 | 1.0000 | 0.3333 | 0.5000 | 0 |
| `tthermophilus` | CDS | partial | 2227 | 2228 | 2173 | 54 | 55 | 0.9758 | 0.9753 | 0.9755 | 1898 |
| `tthermophilus` | tRNA | partial | 49 | 48 | 48 | 1 | 0 | 0.9796 | 1.0000 | 0.9897 | 48 |
| `tthermophilus` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `hpylori` | CDS | partial | 1580 | 1584 | 1524 | 56 | 60 | 0.9646 | 0.9621 | 0.9633 | 1338 |
| `hpylori` | tRNA | partial | 37 | 36 | 36 | 1 | 0 | 0.9730 | 1.0000 | 0.9863 | 36 |
| `hpylori` | rRNA | partial | 2 | 6 | 2 | 0 | 4 | 1.0000 | 0.3333 | 0.5000 | 0 |
| `cjejuni` | CDS | partial | 1673 | 1572 | 1559 | 114 | 13 | 0.9319 | 0.9917 | 0.9609 | 1445 |
| `cjejuni` | tRNA | partial | 45 | 43 | 43 | 2 | 0 | 0.9556 | 1.0000 | 0.9773 | 43 |
| `cjejuni` | rRNA | partial | 3 | 9 | 3 | 0 | 6 | 1.0000 | 0.3333 | 0.5000 | 0 |
| `nmeningitidis` | CDS | partial | 2050 | 2142 | 1943 | 107 | 199 | 0.9478 | 0.9071 | 0.9270 | 1678 |
| `nmeningitidis` | tRNA | partial | 61 | 59 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 59 |
| `nmeningitidis` | rRNA | partial | 8 | 12 | 8 | 0 | 4 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `hinfluenzae` | CDS | partial | 1819 | 1806 | 1775 | 44 | 31 | 0.9758 | 0.9828 | 0.9793 | 1646 |
| `hinfluenzae` | tRNA | partial | 61 | 59 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 59 |
| `hinfluenzae` | rRNA | partial | 13 | 19 | 13 | 0 | 6 | 1.0000 | 0.6842 | 0.8125 | 0 |
| `ctrachomatis` | CDS | partial | 899 | 887 | 874 | 25 | 13 | 0.9722 | 0.9853 | 0.9787 | 780 |
| `ctrachomatis` | tRNA | ok | 37 | 37 | 37 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 31 |
| `ctrachomatis` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `tpallidum` | CDS | partial | 982 | 987 | 958 | 24 | 29 | 0.9756 | 0.9706 | 0.9731 | 704 |
| `tpallidum` | tRNA | ok | 45 | 45 | 45 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 45 |
| `tpallidum` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `bburgdorferi` | CDS | partial | 1534 | 1528 | 1462 | 72 | 66 | 0.9531 | 0.9568 | 0.9549 | 1291 |
| `bburgdorferi` | tRNA | ok | 32 | 32 | 32 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 32 |
| `bburgdorferi` | rRNA | partial | 3 | 5 | 3 | 0 | 2 | 1.0000 | 0.6000 | 0.7500 | 0 |
| `mgenitalium` | CDS | partial | 494 | 524 | 486 | 8 | 38 | 0.9838 | 0.9275 | 0.9548 | 431 |
| `mgenitalium` | tRNA | ok | 36 | 36 | 36 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 36 |
| `mgenitalium` | rRNA | partial | 2 | 3 | 2 | 0 | 1 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `abaumannii` | CDS | partial | 3801 | 3736 | 3696 | 105 | 40 | 0.9724 | 0.9893 | 0.9808 | 3352 |
| `abaumannii` | tRNA | partial | 77 | 74 | 74 | 3 | 0 | 0.9610 | 1.0000 | 0.9801 | 74 |
| `abaumannii` | rRNA | partial | 12 | 18 | 12 | 0 | 6 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `kpneumoniae` | CDS | partial | 5532 | 5779 | 5359 | 173 | 420 | 0.9687 | 0.9273 | 0.9475 | 4695 |
| `kpneumoniae` | tRNA | partial | 91 | 62 | 62 | 29 | 0 | 0.6813 | 1.0000 | 0.8104 | 61 |
| `kpneumoniae` | rRNA | partial | 17 | 25 | 17 | 0 | 8 | 1.0000 | 0.6800 | 0.8095 | 0 |
| `styphimurium` | CDS | partial | 4619 | 4554 | 4435 | 184 | 119 | 0.9602 | 0.9739 | 0.9670 | 3904 |
| `styphimurium` | tRNA | partial | 88 | 85 | 85 | 3 | 0 | 0.9659 | 1.0000 | 0.9827 | 0 |
| `styphimurium` | rRNA | partial | 15 | 22 | 15 | 0 | 7 | 1.0000 | 0.6818 | 0.8108 | 0 |
| `afabrum` | CDS | partial | 5330 | 5267 | 5167 | 163 | 100 | 0.9694 | 0.9810 | 0.9752 | 4569 |
| `afabrum` | tRNA | partial | 55 | 53 | 53 | 2 | 0 | 0.9636 | 1.0000 | 0.9815 | 53 |
| `afabrum` | rRNA | partial | 8 | 12 | 8 | 0 | 4 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `smeliloti` | CDS | partial | 6336 | 6287 | 6049 | 287 | 238 | 0.9547 | 0.9621 | 0.9584 | 5384 |
| `smeliloti` | tRNA | partial | 56 | 55 | 55 | 1 | 0 | 0.9821 | 1.0000 | 0.9910 | 55 |
| `smeliloti` | rRNA | partial | 6 | 9 | 6 | 0 | 3 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `mgryphiswaldense` | CDS | partial | 3945 | 3916 | 3803 | 142 | 113 | 0.9640 | 0.9711 | 0.9675 | 3286 |
| `mgryphiswaldense` | tRNA | partial | 53 | 50 | 50 | 3 | 0 | 0.9434 | 1.0000 | 0.9709 | 50 |
| `mgryphiswaldense` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `soneidensis` | CDS | partial | 4569 | 4491 | 4381 | 188 | 110 | 0.9589 | 0.9755 | 0.9671 | 3891 |
| `soneidensis` | tRNA | partial | 106 | 104 | 104 | 2 | 0 | 0.9811 | 1.0000 | 0.9905 | 104 |
| `soneidensis` | rRNA | partial | 19 | 28 | 19 | 0 | 9 | 1.0000 | 0.6786 | 0.8085 | 0 |
| `cnecator` | CDS | partial | 7907 | 7887 | 7606 | 301 | 281 | 0.9619 | 0.9644 | 0.9631 | 6631 |
| `cnecator` | tRNA | partial | 73 | 69 | 69 | 4 | 0 | 0.9452 | 1.0000 | 0.9718 | 69 |
| `cnecator` | rRNA | partial | 10 | 15 | 10 | 0 | 5 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `zmobilis` | CDS | partial | 1859 | 1852 | 1812 | 47 | 40 | 0.9747 | 0.9784 | 0.9765 | 1635 |
| `zmobilis` | tRNA | partial | 49 | 48 | 48 | 1 | 0 | 0.9796 | 1.0000 | 0.9897 | 48 |
| `zmobilis` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `tmaritima` | CDS | partial | 1854 | 1866 | 1834 | 20 | 32 | 0.9892 | 0.9829 | 0.9860 | 1505 |
| `tmaritima` | tRNA | partial | 47 | 46 | 46 | 1 | 0 | 0.9787 | 1.0000 | 0.9892 | 46 |
| `tmaritima` | rRNA | partial | 2 | 3 | 2 | 0 | 1 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `neuropea` | CDS | partial | 2621 | 2631 | 2499 | 122 | 132 | 0.9535 | 0.9498 | 0.9516 | 2142 |
| `neuropea` | tRNA | partial | 42 | 41 | 41 | 1 | 0 | 0.9762 | 1.0000 | 0.9880 | 41 |
| `neuropea` | rRNA | partial | 2 | 3 | 2 | 0 | 1 | 1.0000 | 0.6667 | 0.8000 | 0 |
| `dvulgaris` | CDS | partial | 3299 | 3249 | 3180 | 119 | 69 | 0.9639 | 0.9788 | 0.9713 | 2772 |
| `dvulgaris` | tRNA | partial | 72 | 67 | 67 | 5 | 0 | 0.9306 | 1.0000 | 0.9641 | 67 |
| `dvulgaris` | rRNA | partial | 11 | 16 | 11 | 0 | 5 | 1.0000 | 0.6875 | 0.8148 | 0 |
| `rpalustris` | CDS | partial | 4947 | 4975 | 4869 | 78 | 106 | 0.9842 | 0.9787 | 0.9814 | 4332 |
| `rpalustris` | tRNA | partial | 51 | 49 | 49 | 2 | 0 | 0.9608 | 1.0000 | 0.9800 | 49 |
| `rpalustris` | rRNA | partial | 4 | 6 | 4 | 0 | 2 | 1.0000 | 0.6667 | 0.8000 | 0 |

## Agregado (micro)

Suma TP/FP/FN de todos los genomas comparables para cada clase.

| clase | TP | FP | FN | precisión | recall | F1 | exactos | % exactos | Δinicio medio | Δfin medio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CDS | 138809 | 4697 | 4097 | 0.9673 | 0.9713 | 0.9693 | 121173 | 87.29% | +0.47 | -0.65 |
| tRNA | 2367 | 86 | 0 | 0.9649 | 1.0000 | 0.9821 | 2144 | 90.58% | +0.05 | +0.14 |
| rRNA | 349 | 0 | 187 | 1.0000 | 0.6511 | 0.7887 | 0 | 0.00% | +2.70 | -2.86 |
| ALL | 141525 | 4783 | 4284 | 0.9673 | 0.9706 | 0.9689 | 123317 | 87.13% | +0.46 | -0.64 |

## Estado y motivos

| genoma | estado | motivo |
|---|---|---|
| `ecoli_k12` | not_available | falta el GFF del sistema: results/ecoli_k12/ecoli_k12.gff3 |
| `subtilis` | partial | TP=4271; FP=77; FN=177; deriva máxima=960 pb |
| `mtb_h37rv` | partial | TP=3790; FP=243; FN=164; deriva máxima=1014 pb |
| `cglutamicum` | partial | TP=2988; FP=125; FN=88; deriva máxima=1176 pb |
| `pputida` | partial | TP=5480; FP=117; FN=88; deriva máxima=861 pb |
| `paeruginosa` | partial | TP=5616; FP=171; FN=32; deriva máxima=639 pb |
| `crescentus` | partial | TP=3755; FP=39; FN=188; deriva máxima=639 pb |
| `synechocystis` | partial | TP=3623; FP=200; FN=122; deriva máxima=693 pb |
| `saureus` | partial | TP=2644; FP=103; FN=198; deriva máxima=189 pb |
| `spneumoniae` | partial | TP=2056; FP=179; FN=131; deriva máxima=894 pb |
| `btheta` | partial | TP=4784; FP=242; FN=102; deriva máxima=1449 pb |
| `scoelicolor` | partial | TP=7671; FP=151; FN=168; deriva máxima=1698 pb |
| `mxanthus` | partial | TP=7232; FP=170; FN=152; deriva máxima=1134 pb |
| `vcholerae` | partial | TP=3701; FP=88; FN=97; deriva máxima=586 pb |
| `cdifficile` | partial | TP=3693; FP=121; FN=43; deriva máxima=691 pb |
| `llactis` | partial | TP=2614; FP=72; FN=42; deriva máxima=692 pb |
| `lmonocytogenes` | partial | TP=2933; FP=67; FN=19; deriva máxima=210 pb |
| `efaecalis` | partial | TP=2722; FP=69; FN=31; deriva máxima=1677 pb |
| `dradiodurans` | partial | TP=3153; FP=51; FN=54; deriva máxima=746 pb |
| `tthermophilus` | partial | TP=2225; FP=55; FN=57; deriva máxima=585 pb |
| `hpylori` | partial | TP=1562; FP=57; FN=64; deriva máxima=1062 pb |
| `cjejuni` | partial | TP=1605; FP=116; FN=19; deriva máxima=1129 pb |
| `nmeningitidis` | partial | TP=2010; FP=109; FN=203; deriva máxima=1290 pb |
| `hinfluenzae` | partial | TP=1847; FP=46; FN=37; deriva máxima=866 pb |
| `ctrachomatis` | partial | TP=915; FP=25; FN=15; deriva máxima=204 pb |
| `tpallidum` | partial | TP=1007; FP=24; FN=31; deriva máxima=699 pb |
| `bburgdorferi` | partial | TP=1497; FP=72; FN=68; deriva máxima=592 pb |
| `mgenitalium` | partial | TP=524; FP=8; FN=39; deriva máxima=381 pb |
| `abaumannii` | partial | TP=3782; FP=108; FN=46; deriva máxima=1159 pb |
| `kpneumoniae` | partial | TP=5438; FP=202; FN=428; deriva máxima=585 pb |
| `styphimurium` | partial | TP=4535; FP=187; FN=126; deriva máxima=636 pb |
| `afabrum` | partial | TP=5228; FP=165; FN=104; deriva máxima=1145 pb |
| `smeliloti` | partial | TP=6110; FP=288; FN=241; deriva máxima=1355 pb |
| `mgryphiswaldense` | partial | TP=3857; FP=145; FN=115; deriva máxima=672 pb |
| `soneidensis` | partial | TP=4504; FP=190; FN=119; deriva máxima=765 pb |
| `cnecator` | partial | TP=7685; FP=305; FN=286; deriva máxima=1008 pb |
| `zmobilis` | partial | TP=1864; FP=48; FN=42; deriva máxima=858 pb |
| `tmaritima` | partial | TP=1882; FP=21; FN=33; deriva máxima=1425 pb |
| `neuropea` | partial | TP=2542; FP=123; FN=133; deriva máxima=1431 pb |
| `dvulgaris` | partial | TP=3258; FP=124; FN=74; deriva máxima=624 pb |
| `rpalustris` | partial | TP=4922; FP=80; FN=108; deriva máxima=548 pb |

## Entradas

| genoma | especie | cepa | status | GFF del sistema | GFF de referencia |
|---|---|---|---|---|---|
| `ecoli_k12` | Escherichia coli | K-12 MG1655 | available | `results/ecoli_k12/ecoli_k12.gff3` | `—` |
| `subtilis` | Bacillus subtilis | 168 | available | `results/subtilis/subtilis.gff3` | `data/annotations/subtilis.gff.gz` |
| `mtb_h37rv` | Mycobacterium tuberculosis | H37Rv | annotated | `results/mtb_h37rv/mtb_h37rv.gff3` | `data/reference_genomes/mtb_h37rv/mtb_h37rv_genomic.gff.gz` |
| `cglutamicum` | Corynebacterium glutamicum | ATCC 13032 | annotated | `results/cglutamicum/cglutamicum.gff3` | `data/reference_genomes/cglutamicum/cglutamicum_genomic.gff.gz` |
| `pputida` | Pseudomonas putida | KT2440 | annotated | `results/pputida/pputida.gff3` | `data/reference_genomes/pputida/pputida_genomic.gff.gz` |
| `paeruginosa` | Pseudomonas aeruginosa | PAO1 | annotated | `results/paeruginosa/paeruginosa.gff3` | `data/reference_genomes/paeruginosa/paeruginosa_genomic.gff.gz` |
| `crescentus` | Caulobacter crescentus | NA1000 | annotated | `results/crescentus/crescentus.gff3` | `data/reference_genomes/crescentus/crescentus_genomic.gff.gz` |
| `synechocystis` | Synechocystis sp. | PCC 6803 | annotated | `results/synechocystis/synechocystis.gff3` | `data/reference_genomes/synechocystis/synechocystis_genomic.gff.gz` |
| `saureus` | Staphylococcus aureus | NCTC 8325 | annotated | `results/saureus/saureus.gff3` | `data/reference_genomes/saureus/saureus_genomic.gff.gz` |
| `spneumoniae` | Streptococcus pneumoniae | D39 | annotated | `results/spneumoniae/spneumoniae.gff3` | `data/reference_genomes/spneumoniae/spneumoniae_genomic.gff.gz` |
| `btheta` | Bacteroides thetaiotaomicron | VPI-5482 | annotated | `results/btheta/btheta.gff3` | `data/reference_genomes/btheta/btheta_genomic.gff.gz` |
| `scoelicolor` | Streptomyces coelicolor | A3(2) | annotated | `results/scoelicolor/scoelicolor.gff3` | `data/reference_genomes/scoelicolor/scoelicolor_genomic.gff.gz` |
| `mxanthus` | Myxococcus xanthus | DK1622 | annotated | `results/mxanthus/mxanthus.gff3` | `data/reference_genomes/mxanthus/mxanthus_genomic.gff.gz` |
| `vcholerae` | Vibrio cholerae | N16961 | annotated | `results/vcholerae/vcholerae.gff3` | `data/reference_genomes/vcholerae/vcholerae_genomic.gff.gz` |
| `cdifficile` | Clostridioides difficile | 630 | annotated | `results/cdifficile/cdifficile.gff3` | `data/reference_genomes/cdifficile/cdifficile_genomic.gff.gz` |
| `llactis` | Lactococcus lactis | IL1403 | annotated | `results/llactis/llactis.gff3` | `data/reference_genomes/llactis/llactis_genomic.gff.gz` |
| `lmonocytogenes` | Listeria monocytogenes | EGD-e | annotated | `results/lmonocytogenes/lmonocytogenes.gff3` | `data/reference_genomes/lmonocytogenes/lmonocytogenes_genomic.gff.gz` |
| `efaecalis` | Enterococcus faecalis | V583 | annotated | `results/efaecalis/efaecalis.gff3` | `data/reference_genomes/efaecalis/efaecalis_genomic.gff.gz` |
| `dradiodurans` | Deinococcus radiodurans | R1 | annotated | `results/dradiodurans/dradiodurans.gff3` | `data/reference_genomes/dradiodurans/dradiodurans_genomic.gff.gz` |
| `tthermophilus` | Thermus thermophilus | HB27 | annotated | `results/tthermophilus/tthermophilus.gff3` | `data/reference_genomes/tthermophilus/tthermophilus_genomic.gff.gz` |
| `hpylori` | Helicobacter pylori | 26695 | annotated | `results/hpylori/hpylori.gff3` | `data/reference_genomes/hpylori/hpylori_genomic.gff.gz` |
| `cjejuni` | Campylobacter jejuni | NCTC 11168 | annotated | `results/cjejuni/cjejuni.gff3` | `data/reference_genomes/cjejuni/cjejuni_genomic.gff.gz` |
| `nmeningitidis` | Neisseria meningitidis | MC58 | annotated | `results/nmeningitidis/nmeningitidis.gff3` | `data/reference_genomes/nmeningitidis/nmeningitidis_genomic.gff.gz` |
| `hinfluenzae` | Haemophilus influenzae | Rd KW20 | annotated | `results/hinfluenzae/hinfluenzae.gff3` | `data/reference_genomes/hinfluenzae/hinfluenzae_genomic.gff.gz` |
| `ctrachomatis` | Chlamydia trachomatis | D/UW-3/CX | annotated | `results/ctrachomatis/ctrachomatis.gff3` | `data/reference_genomes/ctrachomatis/ctrachomatis_genomic.gff.gz` |
| `tpallidum` | Treponema pallidum | Nichols | annotated | `results/tpallidum/tpallidum.gff3` | `data/reference_genomes/tpallidum/tpallidum_genomic.gff.gz` |
| `bburgdorferi` | Borrelia burgdorferi | B31 | annotated | `results/bburgdorferi/bburgdorferi.gff3` | `data/reference_genomes/bburgdorferi/bburgdorferi_genomic.gff.gz` |
| `mgenitalium` | Mycoplasma genitalium | G37 | annotated | `results/mgenitalium/mgenitalium.gff3` | `data/reference_genomes/mgenitalium/mgenitalium_genomic.gff.gz` |
| `abaumannii` | Acinetobacter baumannii | ATCC 17978 | annotated | `results/abaumannii/abaumannii.gff3` | `data/reference_genomes/abaumannii/abaumannii_genomic.gff.gz` |
| `kpneumoniae` | Klebsiella pneumoniae | MGH 78578 | annotated | `results/kpneumoniae/kpneumoniae.gff3` | `data/reference_genomes/kpneumoniae/kpneumoniae_genomic.gff.gz` |
| `styphimurium` | Salmonella enterica | Typhimurium LT2 | annotated | `results/styphimurium/styphimurium.gff3` | `data/reference_genomes/styphimurium/styphimurium_genomic.gff.gz` |
| `afabrum` | Agrobacterium fabrum | C58 | annotated | `results/afabrum/afabrum.gff3` | `data/reference_genomes/afabrum/afabrum_genomic.gff.gz` |
| `smeliloti` | Sinorhizobium meliloti | 1021 | annotated | `results/smeliloti/smeliloti.gff3` | `data/reference_genomes/smeliloti/smeliloti_genomic.gff.gz` |
| `mgryphiswaldense` | Magnetospirillum gryphiswaldense | MSR-1 | annotated | `results/mgryphiswaldense/mgryphiswaldense.gff3` | `data/reference_genomes/mgryphiswaldense/mgryphiswaldense_genomic.gff.gz` |
| `soneidensis` | Shewanella oneidensis | MR-1 | annotated | `results/soneidensis/soneidensis.gff3` | `data/reference_genomes/soneidensis/soneidensis_genomic.gff.gz` |
| `cnecator` | Cupriavidus necator | H16 | annotated | `results/cnecator/cnecator.gff3` | `data/reference_genomes/cnecator/cnecator_genomic.gff.gz` |
| `zmobilis` | Zymomonas mobilis | ZM4 | annotated | `results/zmobilis/zmobilis.gff3` | `data/reference_genomes/zmobilis/zmobilis_genomic.gff.gz` |
| `tmaritima` | Thermotoga maritima | MSB8 | annotated | `results/tmaritima/tmaritima.gff3` | `data/reference_genomes/tmaritima/tmaritima_genomic.gff.gz` |
| `neuropea` | Nitrosomonas europaea | ATCC 19718 | annotated | `results/neuropea/neuropea.gff3` | `data/reference_genomes/neuropea/neuropea_genomic.gff.gz` |
| `dvulgaris` | Desulfovibrio vulgaris | Hildenborough | annotated | `results/dvulgaris/dvulgaris.gff3` | `data/reference_genomes/dvulgaris/dvulgaris_genomic.gff.gz` |
| `rpalustris` | Rhodopseudomonas palustris | CGA009 | annotated | `results/rpalustris/rpalustris.gff3` | `data/reference_genomes/rpalustris/rpalustris_genomic.gff.gz` |

## Criterio de estado

- **ok**: todas las clases con datos tienen `FP=0`, `FN=0`, al menos un par 1:1 y la deriva máxima de bordes no supera la tolerancia.
- **partial**: la comparación fue posible, pero hay FP/FN o algún par con deriva de bordes mayor que la tolerancia.
- **not_available**: falta o es ilegible el GFF del sistema o el de referencia, el manifiesto no declara `annotation_path`, o la referencia no contiene CDS/tRNA/rRNA. No se imputa nada.

## Límites

- **Depende del GFF del sistema**: se compara lo que el pipeline escribió en `results/<id>/<id>.gff3`, no el resumen de conteos.
- **Solo tres clases**: `ncRNA`, `tmRNA`, `pseudogene`, `gene`, `exon` y otros tipos quedan fuera aunque estén en el GFF.
- **Emparejamiento geométrico**: un par 1:1 no demuestra que sea el mismo gen ortólogo ni que la secuencia traducida coincida.
- **Contig y hebra exactos**: diferencias de nomenclatura de contig o hebras no normalizadas producen FP/FN masivos, no un empate silencioso.
- **Multiparte**: varias filas con el mismo `ID` en el mismo contig, tipo y hebra se colapsan a `[min(inicio), max(fin)]`; no se evalúan intrones ni fases.
- **El gold no es verdad absoluta**: hereda criterios y posibles errores de la anotación RefSeq/NCBI o de la fuente usada.
- **Tolerancia**: `ok` es paridad geométrica dentro de la tolerancia, no identidad exacta de bordes; el porcentaje de pares exactos se reporta aparte.
- **Sin causalidad**: un FP o FN no indica por sí solo un error del sistema ni de la referencia; no se usa para entrenamiento.

## Reproducir

```bash
.venv/bin/python scripts/reference/coordinate_parity.py
```

Opciones: `--root`, `--manifest`, `--tsv`, `--md`, `--tolerance`, `--min-reciprocal-overlap`.
