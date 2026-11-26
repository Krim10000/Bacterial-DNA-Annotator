# Paridad por coordenadas: sistema vs GFF de referencia

**Generado:** `2026-11-26T02:10:16Z`  
**Manifiesto:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc/data/reference_genomes.json`  
**Raíz de datos:** `/home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc`  
**Tolerancia de bordes:** 5 pb  
**Solape recíproco mínimo:** 0.50

## Alcance

Compara los features **CDS, tRNA y rRNA** del GFF del sistema (`results/<id>/<id>.gff3`) contra el GFF de referencia (`annotation_path`) para cada registro con status `annotated` o `available` del manifiesto.

El emparejamiento es **1:1** y exige mismo contig, mismo tipo y misma hebra. Dos features son candidatos si se solapan recíprocamente o si sus dos bordes difieren como mucho **5 pb**. Si falta cualquiera de los dos GFF, el genoma queda como `not_available`; no se imputan coordenadas ni se completa con otros archivos.

Esta versión **no evalúa función, producto, identidad de secuencia ni calidad biológica** de las predicciones. Tampoco modifica `results/`, los GFF ni el manifiesto.

## Cobertura

- Registros con status elegible: **42**
- Comparables (ambos GFF legibles): **41** (ok **0**, partial **41**)
- No disponibles: **1**
- Omitidos por status: **0**

### Genomas no disponibles

| genoma | status | motivo |
|---|---|---|
| `ecoli_k12` | available | falta el GFF del sistema: results/ecoli_k12/ecoli_k12.gff3 |

## Resultados por genoma (todas las clases)

TP/FP/FN se calculan sobre CDS + tRNA + rRNA. Precisión, recall y F1 son fracciones (0–1); `% exactos` es un porcentaje.

| genoma | estado | sistema | referencia | TP | FP | FN | precisión | recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `ecoli_k12` | not_available | — | — | — | — | — | — | — | — |
| `subtilis` | partial | 4358 | 4448 | 4281 | 77 | 167 | 0.9823 | 0.9625 | 0.9723 |
| `mtb_h37rv` | partial | 4034 | 3954 | 3791 | 243 | 163 | 0.9398 | 0.9588 | 0.9492 |
| `cglutamicum` | partial | 3119 | 3076 | 2994 | 125 | 82 | 0.9599 | 0.9733 | 0.9666 |
| `pputida` | partial | 5604 | 5568 | 5487 | 117 | 81 | 0.9791 | 0.9855 | 0.9823 |
| `paeruginosa` | partial | 5791 | 5648 | 5620 | 171 | 28 | 0.9705 | 0.9950 | 0.9826 |
| `crescentus` | partial | 3796 | 3943 | 3757 | 39 | 186 | 0.9897 | 0.9528 | 0.9709 |
| `synechocystis` | partial | 3827 | 3745 | 3627 | 200 | 118 | 0.9477 | 0.9685 | 0.9580 |
| `saureus` | partial | 2752 | 2842 | 2649 | 103 | 193 | 0.9626 | 0.9321 | 0.9471 |
| `spneumoniae` | partial | 2239 | 2187 | 2060 | 179 | 127 | 0.9201 | 0.9419 | 0.9309 |
| `btheta` | partial | 5031 | 4886 | 4789 | 242 | 97 | 0.9519 | 0.9801 | 0.9658 |
| `scoelicolor` | partial | 7828 | 7839 | 7677 | 151 | 162 | 0.9807 | 0.9793 | 0.9800 |
| `mxanthus` | partial | 7406 | 7384 | 7236 | 170 | 148 | 0.9770 | 0.9800 | 0.9785 |
| `vcholerae` | partial | 3799 | 3798 | 3711 | 88 | 87 | 0.9768 | 0.9771 | 0.9769 |
| `cdifficile` | partial | 3826 | 3736 | 3705 | 121 | 31 | 0.9684 | 0.9917 | 0.9799 |
| `llactis` | partial | 2692 | 2656 | 2620 | 72 | 36 | 0.9733 | 0.9864 | 0.9798 |
| `lmonocytogenes` | partial | 3006 | 2952 | 2939 | 67 | 13 | 0.9777 | 0.9956 | 0.9866 |
| `efaecalis` | partial | 2795 | 2753 | 2726 | 69 | 27 | 0.9753 | 0.9902 | 0.9827 |
| `dradiodurans` | partial | 3210 | 3207 | 3159 | 51 | 48 | 0.9841 | 0.9850 | 0.9845 |
| `tthermophilus` | partial | 2282 | 2282 | 2227 | 55 | 55 | 0.9759 | 0.9759 | 0.9759 |
| `hpylori` | partial | 1623 | 1626 | 1566 | 57 | 60 | 0.9649 | 0.9631 | 0.9640 |
| `cjejuni` | partial | 1727 | 1624 | 1611 | 116 | 13 | 0.9328 | 0.9920 | 0.9615 |
| `nmeningitidis` | partial | 2123 | 2213 | 2014 | 109 | 199 | 0.9487 | 0.9101 | 0.9290 |
| `hinfluenzae` | partial | 1899 | 1884 | 1853 | 46 | 31 | 0.9758 | 0.9835 | 0.9796 |
| `ctrachomatis` | partial | 942 | 930 | 917 | 25 | 13 | 0.9735 | 0.9860 | 0.9797 |
| `tpallidum` | partial | 1033 | 1038 | 1009 | 24 | 29 | 0.9768 | 0.9721 | 0.9744 |
| `bburgdorferi` | partial | 1571 | 1565 | 1499 | 72 | 66 | 0.9542 | 0.9578 | 0.9560 |
| `mgenitalium` | partial | 533 | 563 | 525 | 8 | 38 | 0.9850 | 0.9325 | 0.9580 |
| `abaumannii` | partial | 3896 | 3828 | 3788 | 108 | 40 | 0.9723 | 0.9896 | 0.9809 |
| `kpneumoniae` | partial | 5648 | 5866 | 5446 | 202 | 420 | 0.9642 | 0.9284 | 0.9460 |
| `styphimurium` | partial | 4729 | 4661 | 4542 | 187 | 119 | 0.9605 | 0.9745 | 0.9674 |
| `afabrum` | partial | 5397 | 5332 | 5232 | 165 | 100 | 0.9694 | 0.9812 | 0.9753 |
| `smeliloti` | partial | 6401 | 6351 | 6113 | 288 | 238 | 0.9550 | 0.9625 | 0.9587 |
| `rsphaeroides` | partial | 4482 | 4468 | 4367 | 115 | 101 | 0.9743 | 0.9774 | 0.9758 |
| `mgryphiswaldense` | partial | 4004 | 3972 | 3859 | 145 | 113 | 0.9638 | 0.9716 | 0.9677 |
| `soneidensis` | partial | 4703 | 4623 | 4513 | 190 | 110 | 0.9596 | 0.9762 | 0.9678 |
| `cnecator` | partial | 7995 | 7971 | 7690 | 305 | 281 | 0.9619 | 0.9647 | 0.9633 |
| `zmobilis` | partial | 1914 | 1906 | 1866 | 48 | 40 | 0.9749 | 0.9790 | 0.9769 |
| `tmaritima` | partial | 1904 | 1915 | 1883 | 21 | 32 | 0.9890 | 0.9833 | 0.9861 |
| `neuropea` | partial | 2666 | 2675 | 2543 | 123 | 132 | 0.9539 | 0.9507 | 0.9523 |
| `dvulgaris` | partial | 3387 | 3332 | 3263 | 124 | 69 | 0.9634 | 0.9793 | 0.9713 |
| `rpalustris` | partial | 5004 | 5030 | 4924 | 80 | 106 | 0.9840 | 0.9789 | 0.9814 |

### Bordes de los pares emparejados

Las diferencias son **sistema − referencia** en pb (GFF 1-based, inclusivo). `exactos` cuenta pares con inicio y fin idénticos.

| genoma | pares | exactos | % exactos | Δinicio medio | Δinicio mediano | Δfin medio | Δfin mediano | deriva máxima |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `ecoli_k12` | — | — | — | — | — | — | — | — |
| `subtilis` | 4281 | 3763 | 87.90% | -0.24 | +0.00 | -0.09 | +0.00 | 960 |
| `mtb_h37rv` | 3791 | 2855 | 75.31% | +4.85 | +0.00 | -4.45 | +0.00 | 1014 |
| `cglutamicum` | 2994 | 2530 | 84.50% | -1.21 | +0.00 | +0.05 | +0.00 | 1176 |
| `pputida` | 5487 | 5034 | 91.74% | -1.90 | +0.00 | +0.53 | +0.00 | 861 |
| `paeruginosa` | 5620 | 5091 | 90.59% | +0.52 | +0.00 | -0.41 | +0.00 | 639 |
| `crescentus` | 3757 | 2999 | 79.82% | +2.69 | +0.00 | -3.45 | +0.00 | 639 |
| `synechocystis` | 3627 | 3049 | 84.06% | +0.00 | +0.00 | -0.24 | +0.00 | 693 |
| `saureus` | 2649 | 2285 | 86.26% | -2.06 | +0.00 | +2.26 | +0.00 | 189 |
| `spneumoniae` | 2060 | 1837 | 89.17% | +4.38 | +0.00 | -4.82 | +0.00 | 894 |
| `btheta` | 4789 | 4242 | 88.58% | -0.93 | +0.00 | +0.55 | +0.00 | 1449 |
| `scoelicolor` | 7677 | 6530 | 85.06% | -0.12 | +0.00 | +0.87 | +0.00 | 1698 |
| `mxanthus` | 7236 | 5983 | 82.68% | -0.47 | +0.00 | +0.17 | +0.00 | 1134 |
| `vcholerae` | 3711 | 3405 | 91.75% | -0.44 | +0.00 | +0.06 | +0.00 | 586 |
| `cdifficile` | 3705 | 3344 | 90.26% | +0.01 | +0.00 | +0.31 | +0.00 | 691 |
| `llactis` | 2620 | 2373 | 90.57% | -0.44 | +0.00 | -0.22 | +0.00 | 692 |
| `lmonocytogenes` | 2939 | 2634 | 89.62% | -0.75 | +0.00 | +1.00 | +0.00 | 210 |
| `efaecalis` | 2726 | 2360 | 86.57% | -0.07 | +0.00 | +1.12 | +0.00 | 1677 |
| `dradiodurans` | 3159 | 2677 | 84.74% | +0.54 | +0.00 | -0.48 | +0.00 | 746 |
| `tthermophilus` | 2227 | 1948 | 87.47% | +2.29 | +0.00 | -1.14 | +0.00 | 585 |
| `hpylori` | 1566 | 1380 | 88.12% | +2.49 | +0.00 | -2.22 | +0.00 | 1062 |
| `cjejuni` | 1611 | 1488 | 92.36% | +1.13 | +0.00 | -1.51 | +0.00 | 1129 |
| `nmeningitidis` | 2014 | 1741 | 86.44% | +2.56 | +0.00 | -7.17 | +0.00 | 1290 |
| `hinfluenzae` | 1853 | 1712 | 92.39% | +0.68 | +0.00 | -2.09 | +0.00 | 866 |
| `ctrachomatis` | 917 | 811 | 88.44% | -0.51 | +0.00 | -1.02 | +0.00 | 204 |
| `tpallidum` | 1009 | 755 | 74.83% | -0.92 | +0.00 | -0.20 | +0.00 | 699 |
| `bburgdorferi` | 1499 | 1328 | 88.59% | +2.24 | +0.00 | -1.27 | +0.00 | 592 |
| `mgenitalium` | 525 | 470 | 89.52% | +4.01 | +0.00 | -2.18 | +0.00 | 381 |
| `abaumannii` | 3788 | 3432 | 90.60% | +0.09 | +0.00 | +0.42 | +0.00 | 1159 |
| `kpneumoniae` | 5446 | 4756 | 87.33% | -1.94 | +0.00 | +2.02 | +0.00 | 585 |
| `styphimurium` | 4542 | 3904 | 85.95% | +0.36 | +0.00 | -1.17 | +0.00 | 636 |
| `afabrum` | 5232 | 4626 | 88.42% | -1.05 | +0.00 | -0.30 | +0.00 | 1145 |
| `smeliloti` | 6113 | 5442 | 89.02% | +2.79 | +0.00 | -3.20 | +0.00 | 1355 |
| `rsphaeroides` | 4367 | 3985 | 91.25% | +0.83 | +0.00 | -0.32 | +0.00 | 699 |
| `mgryphiswaldense` | 3859 | 3338 | 86.50% | +0.80 | +0.00 | -1.01 | +0.00 | 672 |
| `soneidensis` | 4513 | 4005 | 88.74% | +3.40 | +0.00 | -2.05 | +0.00 | 765 |
| `cnecator` | 7690 | 6705 | 87.19% | +1.26 | +0.00 | -1.08 | +0.00 | 1008 |
| `zmobilis` | 1866 | 1685 | 90.30% | -0.90 | +0.00 | +0.76 | +0.00 | 858 |
| `tmaritima` | 1883 | 1552 | 82.42% | -0.48 | +0.00 | -0.31 | +0.00 | 1425 |
| `neuropea` | 2543 | 2184 | 85.88% | +0.76 | +0.00 | -2.88 | +0.00 | 1431 |
| `dvulgaris` | 3263 | 2845 | 87.19% | -0.28 | +0.00 | +1.32 | +0.00 | 624 |
| `rpalustris` | 4924 | 4383 | 89.01% | +1.35 | +0.00 | -0.67 | +0.00 | 548 |

## Resultados por clase

| genoma | clase | estado | sistema | referencia | TP | FP | FN | precisión | recall | F1 | exactos |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ecoli_k12` | CDS | not_available | — | — | — | — | — | — | — | — | — |
| `ecoli_k12` | tRNA | not_available | — | — | — | — | — | — | — | — | — |
| `ecoli_k12` | rRNA | not_available | — | — | — | — | — | — | — | — | — |
| `subtilis` | CDS | partial | 4242 | 4332 | 4165 | 77 | 167 | 0.9818 | 0.9614 | 0.9715 | 3743 |
| `subtilis` | tRNA | ok | 86 | 86 | 86 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 20 |
| `subtilis` | rRNA | ok | 30 | 30 | 30 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `mtb_h37rv` | CDS | partial | 3986 | 3906 | 3743 | 243 | 163 | 0.9390 | 0.9583 | 0.9486 | 2823 |
| `mtb_h37rv` | tRNA | ok | 45 | 45 | 45 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 32 |
| `mtb_h37rv` | rRNA | partial | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `cglutamicum` | CDS | partial | 3041 | 2999 | 2917 | 124 | 82 | 0.9592 | 0.9727 | 0.9659 | 2453 |
| `cglutamicum` | tRNA | partial | 60 | 59 | 59 | 1 | 0 | 0.9833 | 1.0000 | 0.9916 | 59 |
| `cglutamicum` | rRNA | ok | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 18 |
| `pputida` | CDS | partial | 5507 | 5471 | 5390 | 117 | 81 | 0.9788 | 0.9852 | 0.9820 | 4951 |
| `pputida` | tRNA | ok | 75 | 75 | 75 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 75 |
| `pputida` | rRNA | ok | 22 | 22 | 22 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 8 |
| `paeruginosa` | CDS | partial | 5715 | 5572 | 5545 | 170 | 27 | 0.9703 | 0.9952 | 0.9826 | 5075 |
| `paeruginosa` | tRNA | partial | 64 | 63 | 63 | 1 | 0 | 0.9844 | 1.0000 | 0.9921 | 16 |
| `paeruginosa` | rRNA | partial | 12 | 13 | 12 | 0 | 1 | 1.0000 | 0.9231 | 0.9600 | 0 |
| `crescentus` | CDS | partial | 3738 | 3886 | 3700 | 38 | 186 | 0.9898 | 0.9521 | 0.9706 | 2948 |
| `crescentus` | tRNA | partial | 52 | 51 | 51 | 1 | 0 | 0.9808 | 1.0000 | 0.9903 | 51 |
| `crescentus` | rRNA | partial | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `synechocystis` | CDS | partial | 3780 | 3698 | 3580 | 200 | 118 | 0.9471 | 0.9681 | 0.9575 | 3006 |
| `synechocystis` | tRNA | ok | 41 | 41 | 41 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 41 |
| `synechocystis` | rRNA | ok | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 2 |
| `saureus` | CDS | partial | 2675 | 2767 | 2574 | 101 | 193 | 0.9622 | 0.9302 | 0.9459 | 2223 |
| `saureus` | tRNA | partial | 61 | 59 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 56 |
| `saureus` | rRNA | ok | 16 | 16 | 16 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 6 |
| `spneumoniae` | CDS | partial | 2168 | 2117 | 1990 | 178 | 127 | 0.9179 | 0.9400 | 0.9288 | 1775 |
| `spneumoniae` | tRNA | partial | 59 | 58 | 58 | 1 | 0 | 0.9831 | 1.0000 | 0.9915 | 58 |
| `spneumoniae` | rRNA | ok | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 4 |
| `btheta` | CDS | partial | 4945 | 4801 | 4704 | 241 | 97 | 0.9513 | 0.9798 | 0.9653 | 4167 |
| `btheta` | tRNA | partial | 71 | 70 | 70 | 1 | 0 | 0.9859 | 1.0000 | 0.9929 | 70 |
| `btheta` | rRNA | partial | 15 | 15 | 15 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 5 |
| `scoelicolor` | CDS | partial | 7738 | 7755 | 7593 | 145 | 162 | 0.9813 | 0.9791 | 0.9802 | 6458 |
| `scoelicolor` | tRNA | partial | 72 | 66 | 66 | 6 | 0 | 0.9167 | 1.0000 | 0.9565 | 66 |
| `scoelicolor` | rRNA | ok | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 6 |
| `mxanthus` | CDS | partial | 7325 | 7306 | 7158 | 167 | 148 | 0.9772 | 0.9797 | 0.9784 | 5913 |
| `mxanthus` | tRNA | partial | 69 | 66 | 66 | 3 | 0 | 0.9565 | 1.0000 | 0.9778 | 66 |
| `mxanthus` | rRNA | ok | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 4 |
| `vcholerae` | CDS | partial | 3666 | 3665 | 3578 | 88 | 87 | 0.9760 | 0.9763 | 0.9761 | 3292 |
| `vcholerae` | tRNA | ok | 102 | 102 | 102 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 102 |
| `vcholerae` | rRNA | ok | 31 | 31 | 31 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 11 |
| `cdifficile` | CDS | partial | 3700 | 3611 | 3580 | 120 | 31 | 0.9676 | 0.9914 | 0.9794 | 3243 |
| `cdifficile` | tRNA | partial | 91 | 90 | 90 | 1 | 0 | 0.9890 | 1.0000 | 0.9945 | 90 |
| `cdifficile` | rRNA | ok | 35 | 35 | 35 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 11 |
| `llactis` | CDS | partial | 2602 | 2567 | 2531 | 71 | 36 | 0.9727 | 0.9860 | 0.9793 | 2296 |
| `llactis` | tRNA | partial | 71 | 70 | 70 | 1 | 0 | 0.9859 | 1.0000 | 0.9929 | 70 |
| `llactis` | rRNA | ok | 19 | 19 | 19 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 7 |
| `lmonocytogenes` | CDS | partial | 2921 | 2867 | 2854 | 67 | 13 | 0.9771 | 0.9955 | 0.9862 | 2569 |
| `lmonocytogenes` | tRNA | ok | 67 | 67 | 67 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 65 |
| `lmonocytogenes` | rRNA | partial | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `efaecalis` | CDS | partial | 2730 | 2690 | 2663 | 67 | 27 | 0.9755 | 0.9900 | 0.9827 | 2305 |
| `efaecalis` | tRNA | partial | 53 | 51 | 51 | 2 | 0 | 0.9623 | 1.0000 | 0.9808 | 51 |
| `efaecalis` | rRNA | ok | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 4 |
| `dradiodurans` | CDS | partial | 3150 | 3148 | 3100 | 50 | 48 | 0.9841 | 0.9848 | 0.9844 | 2624 |
| `dradiodurans` | tRNA | partial | 51 | 50 | 50 | 1 | 0 | 0.9804 | 1.0000 | 0.9901 | 50 |
| `dradiodurans` | rRNA | ok | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 3 |
| `tthermophilus` | CDS | partial | 2227 | 2228 | 2173 | 54 | 55 | 0.9758 | 0.9753 | 0.9755 | 1898 |
| `tthermophilus` | tRNA | partial | 49 | 48 | 48 | 1 | 0 | 0.9796 | 1.0000 | 0.9897 | 48 |
| `tthermophilus` | rRNA | ok | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 2 |
| `hpylori` | CDS | partial | 1580 | 1584 | 1524 | 56 | 60 | 0.9646 | 0.9621 | 0.9633 | 1338 |
| `hpylori` | tRNA | partial | 37 | 36 | 36 | 1 | 0 | 0.9730 | 1.0000 | 0.9863 | 36 |
| `hpylori` | rRNA | ok | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 6 |
| `cjejuni` | CDS | partial | 1673 | 1572 | 1559 | 114 | 13 | 0.9319 | 0.9917 | 0.9609 | 1445 |
| `cjejuni` | tRNA | partial | 45 | 43 | 43 | 2 | 0 | 0.9556 | 1.0000 | 0.9773 | 43 |
| `cjejuni` | rRNA | partial | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `nmeningitidis` | CDS | partial | 2050 | 2142 | 1943 | 107 | 199 | 0.9478 | 0.9071 | 0.9270 | 1678 |
| `nmeningitidis` | tRNA | partial | 61 | 59 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 59 |
| `nmeningitidis` | rRNA | ok | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 4 |
| `hinfluenzae` | CDS | partial | 1819 | 1806 | 1775 | 44 | 31 | 0.9758 | 0.9828 | 0.9793 | 1646 |
| `hinfluenzae` | tRNA | partial | 61 | 59 | 59 | 2 | 0 | 0.9672 | 1.0000 | 0.9833 | 59 |
| `hinfluenzae` | rRNA | ok | 19 | 19 | 19 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 7 |
| `ctrachomatis` | CDS | partial | 899 | 887 | 874 | 25 | 13 | 0.9722 | 0.9853 | 0.9787 | 780 |
| `ctrachomatis` | tRNA | ok | 37 | 37 | 37 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 31 |
| `ctrachomatis` | rRNA | partial | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `tpallidum` | CDS | partial | 982 | 987 | 958 | 24 | 29 | 0.9756 | 0.9706 | 0.9731 | 704 |
| `tpallidum` | tRNA | ok | 45 | 45 | 45 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 45 |
| `tpallidum` | rRNA | ok | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 6 |
| `bburgdorferi` | CDS | partial | 1534 | 1528 | 1462 | 72 | 66 | 0.9531 | 0.9568 | 0.9549 | 1291 |
| `bburgdorferi` | tRNA | ok | 32 | 32 | 32 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 32 |
| `bburgdorferi` | rRNA | ok | 5 | 5 | 5 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 5 |
| `mgenitalium` | CDS | partial | 494 | 524 | 486 | 8 | 38 | 0.9838 | 0.9275 | 0.9548 | 431 |
| `mgenitalium` | tRNA | ok | 36 | 36 | 36 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 36 |
| `mgenitalium` | rRNA | ok | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 3 |
| `abaumannii` | CDS | partial | 3801 | 3736 | 3696 | 105 | 40 | 0.9724 | 0.9893 | 0.9808 | 3352 |
| `abaumannii` | tRNA | partial | 77 | 74 | 74 | 3 | 0 | 0.9610 | 1.0000 | 0.9801 | 74 |
| `abaumannii` | rRNA | ok | 18 | 18 | 18 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 6 |
| `kpneumoniae` | CDS | partial | 5532 | 5779 | 5359 | 173 | 420 | 0.9687 | 0.9273 | 0.9475 | 4695 |
| `kpneumoniae` | tRNA | partial | 91 | 62 | 62 | 29 | 0 | 0.6813 | 1.0000 | 0.8104 | 61 |
| `kpneumoniae` | rRNA | ok | 25 | 25 | 25 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `styphimurium` | CDS | partial | 4619 | 4554 | 4435 | 184 | 119 | 0.9602 | 0.9739 | 0.9670 | 3904 |
| `styphimurium` | tRNA | partial | 88 | 85 | 85 | 3 | 0 | 0.9659 | 1.0000 | 0.9827 | 0 |
| `styphimurium` | rRNA | partial | 22 | 22 | 22 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 |
| `afabrum` | CDS | partial | 5330 | 5267 | 5167 | 163 | 100 | 0.9694 | 0.9810 | 0.9752 | 4569 |
| `afabrum` | tRNA | partial | 55 | 53 | 53 | 2 | 0 | 0.9636 | 1.0000 | 0.9815 | 53 |
| `afabrum` | rRNA | partial | 12 | 12 | 12 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 4 |
| `smeliloti` | CDS | partial | 6336 | 6287 | 6049 | 287 | 238 | 0.9547 | 0.9621 | 0.9584 | 5384 |
| `smeliloti` | tRNA | partial | 56 | 55 | 55 | 1 | 0 | 0.9821 | 1.0000 | 0.9910 | 55 |
| `smeliloti` | rRNA | partial | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 3 |
| `rsphaeroides` | CDS | partial | 4417 | 4406 | 4305 | 112 | 101 | 0.9746 | 0.9771 | 0.9758 | 3929 |
| `rsphaeroides` | tRNA | partial | 56 | 53 | 53 | 3 | 0 | 0.9464 | 1.0000 | 0.9725 | 53 |
| `rsphaeroides` | rRNA | partial | 9 | 9 | 9 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 3 |
| `mgryphiswaldense` | CDS | partial | 3945 | 3916 | 3803 | 142 | 113 | 0.9640 | 0.9711 | 0.9675 | 3286 |
| `mgryphiswaldense` | tRNA | partial | 53 | 50 | 50 | 3 | 0 | 0.9434 | 1.0000 | 0.9709 | 50 |
| `mgryphiswaldense` | rRNA | partial | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 2 |
| `soneidensis` | CDS | partial | 4569 | 4491 | 4381 | 188 | 110 | 0.9589 | 0.9755 | 0.9671 | 3891 |
| `soneidensis` | tRNA | partial | 106 | 104 | 104 | 2 | 0 | 0.9811 | 1.0000 | 0.9905 | 104 |
| `soneidensis` | rRNA | ok | 28 | 28 | 28 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 10 |
| `cnecator` | CDS | partial | 7907 | 7887 | 7606 | 301 | 281 | 0.9619 | 0.9644 | 0.9631 | 6631 |
| `cnecator` | tRNA | partial | 73 | 69 | 69 | 4 | 0 | 0.9452 | 1.0000 | 0.9718 | 69 |
| `cnecator` | rRNA | ok | 15 | 15 | 15 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 5 |
| `zmobilis` | CDS | partial | 1859 | 1852 | 1812 | 47 | 40 | 0.9747 | 0.9784 | 0.9765 | 1635 |
| `zmobilis` | tRNA | partial | 49 | 48 | 48 | 1 | 0 | 0.9796 | 1.0000 | 0.9897 | 48 |
| `zmobilis` | rRNA | partial | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 2 |
| `tmaritima` | CDS | partial | 1854 | 1866 | 1834 | 20 | 32 | 0.9892 | 0.9829 | 0.9860 | 1505 |
| `tmaritima` | tRNA | partial | 47 | 46 | 46 | 1 | 0 | 0.9787 | 1.0000 | 0.9892 | 46 |
| `tmaritima` | rRNA | ok | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 1 |
| `neuropea` | CDS | partial | 2621 | 2631 | 2499 | 122 | 132 | 0.9535 | 0.9498 | 0.9516 | 2142 |
| `neuropea` | tRNA | partial | 42 | 41 | 41 | 1 | 0 | 0.9762 | 1.0000 | 0.9880 | 41 |
| `neuropea` | rRNA | ok | 3 | 3 | 3 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 1 |
| `dvulgaris` | CDS | partial | 3299 | 3249 | 3180 | 119 | 69 | 0.9639 | 0.9788 | 0.9713 | 2772 |
| `dvulgaris` | tRNA | partial | 72 | 67 | 67 | 5 | 0 | 0.9306 | 1.0000 | 0.9641 | 67 |
| `dvulgaris` | rRNA | ok | 16 | 16 | 16 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 6 |
| `rpalustris` | CDS | partial | 4947 | 4975 | 4869 | 78 | 106 | 0.9842 | 0.9787 | 0.9814 | 4332 |
| `rpalustris` | tRNA | partial | 51 | 49 | 49 | 2 | 0 | 0.9608 | 1.0000 | 0.9800 | 49 |
| `rpalustris` | rRNA | partial | 6 | 6 | 6 | 0 | 0 | 1.0000 | 1.0000 | 1.0000 | 2 |

## Agregado (micro)

Suma TP/FP/FN de todos los genomas comparables para cada clase.

| clase | TP | FP | FN | precisión | recall | F1 | exactos | % exactos | Δinicio medio | Δfin medio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CDS | 143114 | 4809 | 4198 | 0.9675 | 0.9715 | 0.9695 | 125102 | 87.41% | +0.48 | -0.64 |
| tRNA | 2420 | 89 | 0 | 0.9645 | 1.0000 | 0.9819 | 2197 | 90.79% | +0.05 | +0.14 |
| rRNA | 544 | 0 | 1 | 1.0000 | 0.9982 | 0.9991 | 167 | 30.70% | +1.01 | -1.59 |
| ALL | 146078 | 4898 | 4199 | 0.9676 | 0.9721 | 0.9698 | 127466 | 87.26% | +0.47 | -0.63 |

## Estado y motivos

| genoma | estado | motivo |
|---|---|---|
| `ecoli_k12` | not_available | falta el GFF del sistema: results/ecoli_k12/ecoli_k12.gff3 |
| `subtilis` | partial | TP=4281; FP=77; FN=167; deriva máxima=960 pb |
| `mtb_h37rv` | partial | TP=3791; FP=243; FN=163; deriva máxima=1014 pb |
| `cglutamicum` | partial | TP=2994; FP=125; FN=82; deriva máxima=1176 pb |
| `pputida` | partial | TP=5487; FP=117; FN=81; deriva máxima=861 pb |
| `paeruginosa` | partial | TP=5620; FP=171; FN=28; deriva máxima=639 pb |
| `crescentus` | partial | TP=3757; FP=39; FN=186; deriva máxima=639 pb |
| `synechocystis` | partial | TP=3627; FP=200; FN=118; deriva máxima=693 pb |
| `saureus` | partial | TP=2649; FP=103; FN=193; deriva máxima=189 pb |
| `spneumoniae` | partial | TP=2060; FP=179; FN=127; deriva máxima=894 pb |
| `btheta` | partial | TP=4789; FP=242; FN=97; deriva máxima=1449 pb |
| `scoelicolor` | partial | TP=7677; FP=151; FN=162; deriva máxima=1698 pb |
| `mxanthus` | partial | TP=7236; FP=170; FN=148; deriva máxima=1134 pb |
| `vcholerae` | partial | TP=3711; FP=88; FN=87; deriva máxima=586 pb |
| `cdifficile` | partial | TP=3705; FP=121; FN=31; deriva máxima=691 pb |
| `llactis` | partial | TP=2620; FP=72; FN=36; deriva máxima=692 pb |
| `lmonocytogenes` | partial | TP=2939; FP=67; FN=13; deriva máxima=210 pb |
| `efaecalis` | partial | TP=2726; FP=69; FN=27; deriva máxima=1677 pb |
| `dradiodurans` | partial | TP=3159; FP=51; FN=48; deriva máxima=746 pb |
| `tthermophilus` | partial | TP=2227; FP=55; FN=55; deriva máxima=585 pb |
| `hpylori` | partial | TP=1566; FP=57; FN=60; deriva máxima=1062 pb |
| `cjejuni` | partial | TP=1611; FP=116; FN=13; deriva máxima=1129 pb |
| `nmeningitidis` | partial | TP=2014; FP=109; FN=199; deriva máxima=1290 pb |
| `hinfluenzae` | partial | TP=1853; FP=46; FN=31; deriva máxima=866 pb |
| `ctrachomatis` | partial | TP=917; FP=25; FN=13; deriva máxima=204 pb |
| `tpallidum` | partial | TP=1009; FP=24; FN=29; deriva máxima=699 pb |
| `bburgdorferi` | partial | TP=1499; FP=72; FN=66; deriva máxima=592 pb |
| `mgenitalium` | partial | TP=525; FP=8; FN=38; deriva máxima=381 pb |
| `abaumannii` | partial | TP=3788; FP=108; FN=40; deriva máxima=1159 pb |
| `kpneumoniae` | partial | TP=5446; FP=202; FN=420; deriva máxima=585 pb |
| `styphimurium` | partial | TP=4542; FP=187; FN=119; deriva máxima=636 pb |
| `afabrum` | partial | TP=5232; FP=165; FN=100; deriva máxima=1145 pb |
| `smeliloti` | partial | TP=6113; FP=288; FN=238; deriva máxima=1355 pb |
| `rsphaeroides` | partial | TP=4367; FP=115; FN=101; deriva máxima=699 pb |
| `mgryphiswaldense` | partial | TP=3859; FP=145; FN=113; deriva máxima=672 pb |
| `soneidensis` | partial | TP=4513; FP=190; FN=110; deriva máxima=765 pb |
| `cnecator` | partial | TP=7690; FP=305; FN=281; deriva máxima=1008 pb |
| `zmobilis` | partial | TP=1866; FP=48; FN=40; deriva máxima=858 pb |
| `tmaritima` | partial | TP=1883; FP=21; FN=32; deriva máxima=1425 pb |
| `neuropea` | partial | TP=2543; FP=123; FN=132; deriva máxima=1431 pb |
| `dvulgaris` | partial | TP=3263; FP=124; FN=69; deriva máxima=624 pb |
| `rpalustris` | partial | TP=4924; FP=80; FN=106; deriva máxima=548 pb |

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
| `rsphaeroides` | Rhodobacter sphaeroides | 2.4.1 | annotated | `results/rsphaeroides/rsphaeroides.gff3` | `data/reference_genomes/rsphaeroides/rsphaeroides_genomic.gff.gz` |
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
