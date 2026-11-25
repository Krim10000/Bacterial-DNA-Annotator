# Catálogo unificado de estructuras AFDB + ESMFold

Generado por `scripts/structure/unify_structure_catalog.py`.
Este reporte no es evidencia experimental: los modelos son predicciones computacionales.

## Resumen

| métrica | valor |
|---|---:|
| Filas AFDB de entrada | 620 |
| Filas ESMFold de entrada | 1565 |
| Accesiones unificadas | 2184 |
|   · fuente AFDB | 620 |
|   · fuente ESMFold | 1564 |
|   · estado ok | 2042 |
|   · estado missing | 0 |
|   · estado error | 142 |
|   · estado not_found | 0 |
| Duplicados por accesión | 1 |
| Duplicados por hash | 0 |
| PDB huérfanos | 1192 |
| Rutas no coincidentes | 0 |
| Archivos faltantes | 0 |

## Reglas aplicadas

- Una fila por accesión.
- Misma accesión en AFDB y ESMFold: se prefiere AFDB si su archivo es utilizable; si no, se conserva ESMFold.
- Mismo SHA-256 en accesiones distintas: canónica = archivo utilizable > estado > AFDB > mayor pLDDT > accesión menor.
- Los duplicados no se borran: quedan en el catálogo con `duplicate_of` y en este reporte.

## Duplicados por accesión

| accession | source | path | bytes | related | detail |
|---|---|---|---:|---|---|
| A0A089QKZ7 | ESMFold | data/metabolism/targets/esmfold_kaggle/A0A089QKZ7.pdb | 35520 | AFDB:data/metabolism/targets/afdb25k/A0A089QKZ7.pdb | AFDB preferido para la misma accesión |

## Duplicados por hash

Sin hallazgos.

## PDB huérfanos

| accession | source | path | bytes | related | detail |
|---|---|---|---:|---|---|
| A0A0K2JLU1 | AFDB | data/metabolism/targets/afdb25k/A0A0K2JLU1.pdb | 352673 |  | accesión A0A0K2JLU1 registrada sin ruta en AFDB (estado: error) |
| A0A0K2VM55 | AFDB | data/metabolism/targets/afdb25k/A0A0K2VM55.pdb | 283661 |  | accesión A0A0K2VM55 registrada sin ruta en AFDB (estado: error) |
| A0A0K3AUJ9 | AFDB | data/metabolism/targets/afdb25k/A0A0K3AUJ9.pdb | 133568 |  | accesión A0A0K3AUJ9 registrada sin ruta en AFDB (estado: error) |
| A0A0K6IWT9 | AFDB | data/metabolism/targets/afdb25k/A0A0K6IWT9.pdb | 288602 |  | accesión A0A0K6IWT9 registrada sin ruta en AFDB (estado: error) |
| A0QAB6 | AFDB | data/metabolism/targets/afdb25k/A0QAB6.pdb | 189701 |  | accesión A0QAB6 no registrada en el catálogo AFDB |
| A0R1B5 | AFDB | data/metabolism/targets/afdb25k/A0R1B5.pdb | 95660 |  | accesión A0R1B5 no registrada en el catálogo AFDB |
| A0R557 | AFDB | data/metabolism/targets/afdb25k/A0R557.pdb | 191483 |  | accesión A0R557 no registrada en el catálogo AFDB |
| A0RCV1 | AFDB | data/metabolism/targets/afdb25k/A0RCV1.pdb | 103598 |  | accesión A0RCV1 no registrada en el catálogo AFDB |
| A0RIY8 | AFDB | data/metabolism/targets/afdb25k/A0RIY8.pdb | 144584 |  | accesión A0RIY8 no registrada en el catálogo AFDB |
| A1JMC0 | AFDB | data/metabolism/targets/afdb25k/A1JMC0.pdb | 129599 |  | accesión A1JMC0 no registrada en el catálogo AFDB |
| A4IR67 | AFDB | data/metabolism/targets/afdb25k/A4IR67.pdb | 143774 |  | accesión A4IR67 no registrada en el catálogo AFDB |
| A6U5H5 | AFDB | data/metabolism/targets/afdb25k/A6U5H5.pdb | 91448 |  | accesión A6U5H5 no registrada en el catálogo AFDB |
| A7MKE8 | AFDB | data/metabolism/targets/afdb25k/A7MKE8.pdb | 126926 |  | accesión A7MKE8 no registrada en el catálogo AFDB |
| A8MGS9 | AFDB | data/metabolism/targets/afdb25k/A8MGS9.pdb | 148310 |  | accesión A8MGS9 no registrada en el catálogo AFDB |
| B0RVK2 | AFDB | data/metabolism/targets/afdb25k/B0RVK2.pdb | 460727 |  | accesión B0RVK2 no registrada en el catálogo AFDB |
| B7NAP9 | AFDB | data/metabolism/targets/afdb25k/B7NAP9.pdb | 242513 |  | accesión B7NAP9 no registrada en el catálogo AFDB |
| B7NSW2 | AFDB | data/metabolism/targets/afdb25k/B7NSW2.pdb | 127898 |  | accesión B7NSW2 no registrada en el catálogo AFDB |
| C0H3Q8 | AFDB | data/metabolism/targets/afdb25k/C0H3Q8.pdb | 65042 |  | accesión C0H3Q8 no registrada en el catálogo AFDB |
| C0H3V8 | AFDB | data/metabolism/targets/afdb25k/C0H3V8.pdb | 48842 |  | accesión C0H3V8 no registrada en el catálogo AFDB |
| C0H463 | AFDB | data/metabolism/targets/afdb25k/C0H463.pdb | 45764 |  | accesión C0H463 no registrada en el catálogo AFDB |
| C0QYX7 | AFDB | data/metabolism/targets/afdb25k/C0QYX7.pdb | 242270 |  | accesión C0QYX7 no registrada en el catálogo AFDB |
| C0SPB4 | AFDB | data/metabolism/targets/afdb25k/C0SPB4.pdb | 198935 |  | accesión C0SPB4 no registrada en el catálogo AFDB |
| C6Y4B8 | AFDB | data/metabolism/targets/afdb25k/C6Y4B8.pdb | 85697 |  | accesión C6Y4B8 no registrada en el catálogo AFDB |
| D4AK18 | AFDB | data/metabolism/targets/afdb25k/D4AK18.pdb | 130895 |  | accesión D4AK18 no registrada en el catálogo AFDB |
| O05083 | AFDB | data/metabolism/targets/afdb25k/O05083.pdb | 234656 |  | accesión O05083 no registrada en el catálogo AFDB |
| O05086 | AFDB | data/metabolism/targets/afdb25k/O05086.pdb | 150983 |  | accesión O05086 no registrada en el catálogo AFDB |
| O05235 | AFDB | data/metabolism/targets/afdb25k/O05235.pdb | 184760 |  | accesión O05235 no registrada en el catálogo AFDB |
| O05399 | AFDB | data/metabolism/targets/afdb25k/O05399.pdb | 168560 |  | accesión O05399 no registrada en el catálogo AFDB |
| O05405 | AFDB | data/metabolism/targets/afdb25k/O05405.pdb | 187028 |  | accesión O05405 no registrada en el catálogo AFDB |
| O05406 | AFDB | data/metabolism/targets/afdb25k/O05406.pdb | 138752 |  | accesión O05406 no registrada en el catálogo AFDB |
| O05408 | AFDB | data/metabolism/targets/afdb25k/O05408.pdb | 215135 |  | accesión O05408 no registrada en el catálogo AFDB |
| O05494 | AFDB | data/metabolism/targets/afdb25k/O05494.pdb | 156167 |  | accesión O05494 no registrada en el catálogo AFDB |
| O05502 | AFDB | data/metabolism/targets/afdb25k/O05502.pdb | 220805 |  | accesión O05502 no registrada en el catálogo AFDB |
| O05945 | AFDB | data/metabolism/targets/afdb25k/O05945.pdb | 273050 |  | accesión O05945 no registrada en el catálogo AFDB |
| O05972 | AFDB | data/metabolism/targets/afdb25k/O05972.pdb | 169289 |  | accesión O05972 no registrada en el catálogo AFDB |
| O05979 | AFDB | data/metabolism/targets/afdb25k/O05979.pdb | 368549 |  | accesión O05979 no registrada en el catálogo AFDB |
| O06012 | AFDB | data/metabolism/targets/afdb25k/O06012.pdb | 240569 |  | accesión O06012 no registrada en el catálogo AFDB |
| O06342 | AFDB | data/metabolism/targets/afdb25k/O06342.pdb | 671732 |  | accesión O06342 no registrada en el catálogo AFDB |
| O06391 | AFDB | data/metabolism/targets/afdb25k/O06391.pdb | 132515 |  | accesión O06391 no registrada en el catálogo AFDB |
| O06476 | AFDB | data/metabolism/targets/afdb25k/O06476.pdb | 413261 |  | accesión O06476 no registrada en el catálogo AFDB |
| O06481 | AFDB | data/metabolism/targets/afdb25k/O06481.pdb | 254744 |  | accesión O06481 no registrada en el catálogo AFDB |
| O06483 | AFDB | data/metabolism/targets/afdb25k/O06483.pdb | 266165 |  | accesión O06483 no registrada en el catálogo AFDB |
| O06747 | AFDB | data/metabolism/targets/afdb25k/O06747.pdb | 199016 |  | accesión O06747 no registrada en el catálogo AFDB |
| O06972 | AFDB | data/metabolism/targets/afdb25k/O06972.pdb | 108863 |  | accesión O06972 no registrada en el catálogo AFDB |
| O06978 | AFDB | data/metabolism/targets/afdb25k/O06978.pdb | 162566 |  | accesión O06978 no registrada en el catálogo AFDB |
| O07001 | AFDB | data/metabolism/targets/afdb25k/O07001.pdb | 132434 |  | accesión O07001 no registrada en el catálogo AFDB |
| O07016 | AFDB | data/metabolism/targets/afdb25k/O07016.pdb | 198044 |  | accesión O07016 no registrada en el catálogo AFDB |
| O07019 | AFDB | data/metabolism/targets/afdb25k/O07019.pdb | 131705 |  | accesión O07019 no registrada en el catálogo AFDB |
| O07079 | AFDB | data/metabolism/targets/afdb25k/O07079.pdb | 113156 |  | accesión O07079 no registrada en el catálogo AFDB |
| O07085 | AFDB | data/metabolism/targets/afdb25k/O07085.pdb | 229877 |  | accesión O07085 no registrada en el catálogo AFDB |
| O07522 | AFDB | data/metabolism/targets/afdb25k/O07522.pdb | 274589 |  | accesión O07522 no registrada en el catálogo AFDB |
| O07523 | AFDB | data/metabolism/targets/afdb25k/O07523.pdb | 265355 |  | accesión O07523 no registrada en el catálogo AFDB |
| O07541 | AFDB | data/metabolism/targets/afdb25k/O07541.pdb | 245672 |  | accesión O07541 no registrada en el catálogo AFDB |
| O07548 | AFDB | data/metabolism/targets/afdb25k/O07548.pdb | 134864 |  | accesión O07548 no registrada en el catálogo AFDB |
| O07556 | AFDB | data/metabolism/targets/afdb25k/O07556.pdb | 309014 |  | accesión O07556 no registrada en el catálogo AFDB |
| O07561 | AFDB | data/metabolism/targets/afdb25k/O07561.pdb | 316547 |  | accesión O07561 no registrada en el catálogo AFDB |
| O07577 | AFDB | data/metabolism/targets/afdb25k/O07577.pdb | 287549 |  | accesión O07577 no registrada en el catálogo AFDB |
| O07599 | AFDB | data/metabolism/targets/afdb25k/O07599.pdb | 285848 |  | accesión O07599 no registrada en el catálogo AFDB |
| O07601 | AFDB | data/metabolism/targets/afdb25k/O07601.pdb | 170099 |  | accesión O07601 no registrada en el catálogo AFDB |
| O07940 | AFDB | data/metabolism/targets/afdb25k/O07940.pdb | 291599 |  | accesión O07940 no registrada en el catálogo AFDB |
| O08455 | AFDB | data/metabolism/targets/afdb25k/O08455.pdb | 644435 |  | accesión O08455 no registrada en el catálogo AFDB |
| O25971 | AFDB | data/metabolism/targets/afdb25k/O25971.pdb | 209546 |  | accesión O25971 no registrada en el catálogo AFDB |
| O27068 | AFDB | data/metabolism/targets/afdb25k/O27068.pdb | 246806 |  | accesión O27068 no registrada en el catálogo AFDB |
| O27346 | AFDB | data/metabolism/targets/afdb25k/O27346.pdb | 88775 |  | accesión O27346 no registrada en el catálogo AFDB |
| O27383 | AFDB | data/metabolism/targets/afdb25k/O27383.pdb | 49814 |  | accesión O27383 no registrada en el catálogo AFDB |
| O28519 | AFDB | data/metabolism/targets/afdb25k/O28519.pdb | 139238 |  | accesión O28519 no registrada en el catálogo AFDB |
| O28550 | AFDB | data/metabolism/targets/afdb25k/O28550.pdb | 198449 |  | accesión O28550 no registrada en el catálogo AFDB |
| O28710 | AFDB | data/metabolism/targets/afdb25k/O28710.pdb | 94607 |  | accesión O28710 no registrada en el catálogo AFDB |
| O28724 | AFDB | data/metabolism/targets/afdb25k/O28724.pdb | 117773 |  | accesión O28724 no registrada en el catálogo AFDB |
| O28726 | AFDB | data/metabolism/targets/afdb25k/O28726.pdb | 187271 |  | accesión O28726 no registrada en el catálogo AFDB |
| O28865 | AFDB | data/metabolism/targets/afdb25k/O28865.pdb | 108539 |  | accesión O28865 no registrada en el catálogo AFDB |
| O28987 | AFDB | data/metabolism/targets/afdb25k/O28987.pdb | 113156 |  | accesión O28987 no registrada en el catálogo AFDB |
| O29056 | AFDB | data/metabolism/targets/afdb25k/O29056.pdb | 223559 |  | accesión O29056 no registrada en el catálogo AFDB |
| O29139 | AFDB | data/metabolism/targets/afdb25k/O29139.pdb | 95093 |  | accesión O29139 no registrada en el catálogo AFDB |
| O29470 | AFDB | data/metabolism/targets/afdb25k/O29470.pdb | 202823 |  | accesión O29470 no registrada en el catálogo AFDB |
| O29510 | AFDB | data/metabolism/targets/afdb25k/O29510.pdb | 188486 |  | accesión O29510 no registrada en el catálogo AFDB |
| O29828 | AFDB | data/metabolism/targets/afdb25k/O29828.pdb | 155033 |  | accesión O29828 no registrada en el catálogo AFDB |
| O29891 | AFDB | data/metabolism/targets/afdb25k/O29891.pdb | 162080 |  | accesión O29891 no registrada en el catálogo AFDB |
| O29973 | AFDB | data/metabolism/targets/afdb25k/O29973.pdb | 178442 |  | accesión O29973 no registrada en el catálogo AFDB |
| O30030 | AFDB | data/metabolism/targets/afdb25k/O30030.pdb | 316142 |  | accesión O30030 no registrada en el catálogo AFDB |
| O30272 | AFDB | data/metabolism/targets/afdb25k/O30272.pdb | 144341 |  | accesión O30272 no registrada en el catálogo AFDB |
| O31398 | AFDB | data/metabolism/targets/afdb25k/O31398.pdb | 164996 |  | accesión O31398 no registrada en el catálogo AFDB |
| O31448 | AFDB | data/metabolism/targets/afdb25k/O31448.pdb | 199340 |  | accesión O31448 no registrada en el catálogo AFDB |
| O31472 | AFDB | data/metabolism/targets/afdb25k/O31472.pdb | 102545 |  | accesión O31472 no registrada en el catálogo AFDB |
| O31480 | AFDB | data/metabolism/targets/afdb25k/O31480.pdb | 72008 |  | accesión O31480 no registrada en el catálogo AFDB |
| O31497 | AFDB | data/metabolism/targets/afdb25k/O31497.pdb | 91853 |  | accesión O31497 no registrada en el catálogo AFDB |
| O31500 | AFDB | data/metabolism/targets/afdb25k/O31500.pdb | 192536 |  | accesión O31500 no registrada en el catálogo AFDB |
| O31517 | AFDB | data/metabolism/targets/afdb25k/O31517.pdb | 237896 |  | accesión O31517 no registrada en el catálogo AFDB |
| O31530 | AFDB | data/metabolism/targets/afdb25k/O31530.pdb | 567647 |  | accesión O31530 no registrada en el catálogo AFDB |
| O31540 | AFDB | data/metabolism/targets/afdb25k/O31540.pdb | 210599 |  | accesión O31540 no registrada en el catálogo AFDB |
| O31560 | AFDB | data/metabolism/targets/afdb25k/O31560.pdb | 140615 |  | accesión O31560 no registrada en el catálogo AFDB |
| O31561 | AFDB | data/metabolism/targets/afdb25k/O31561.pdb | 259685 |  | accesión O31561 no registrada en el catálogo AFDB |
| O31563 | AFDB | data/metabolism/targets/afdb25k/O31563.pdb | 321083 |  | accesión O31563 no registrada en el catálogo AFDB |
| O31580 | AFDB | data/metabolism/targets/afdb25k/O31580.pdb | 72575 |  | accesión O31580 no registrada en el catálogo AFDB |
| O31603 | AFDB | data/metabolism/targets/afdb25k/O31603.pdb | 140777 |  | accesión O31603 no registrada en el catálogo AFDB |
| O31670 | AFDB | data/metabolism/targets/afdb25k/O31670.pdb | 74195 |  | accesión O31670 no registrada en el catálogo AFDB |
| O31674 | AFDB | data/metabolism/targets/afdb25k/O31674.pdb | 222749 |  | accesión O31674 no registrada en el catálogo AFDB |
| O31680 | AFDB | data/metabolism/targets/afdb25k/O31680.pdb | 155600 |  | accesión O31680 no registrada en el catálogo AFDB |
| O31711 | AFDB | data/metabolism/targets/afdb25k/O31711.pdb | 148715 |  | accesión O31711 no registrada en el catálogo AFDB |
| O31712 | AFDB | data/metabolism/targets/afdb25k/O31712.pdb | 245186 |  | accesión O31712 no registrada en el catálogo AFDB |
| O31716 | AFDB | data/metabolism/targets/afdb25k/O31716.pdb | 355751 |  | accesión O31716 no registrada en el catálogo AFDB |
| O31717 | AFDB | data/metabolism/targets/afdb25k/O31717.pdb | 197558 |  | accesión O31717 no registrada en el catálogo AFDB |
| O31723 | AFDB | data/metabolism/targets/afdb25k/O31723.pdb | 175364 |  | accesión O31723 no registrada en el catálogo AFDB |
| O31729 | AFDB | data/metabolism/targets/afdb25k/O31729.pdb | 63503 |  | accesión O31729 no registrada en el catálogo AFDB |
| O31761 | AFDB | data/metabolism/targets/afdb25k/O31761.pdb | 164672 |  | accesión O31761 no registrada en el catálogo AFDB |
| O31766 | AFDB | data/metabolism/targets/afdb25k/O31766.pdb | 286253 |  | accesión O31766 no registrada en el catálogo AFDB |
| O31767 | AFDB | data/metabolism/targets/afdb25k/O31767.pdb | 150497 |  | accesión O31767 no registrada en el catálogo AFDB |
| O31771 | AFDB | data/metabolism/targets/afdb25k/O31771.pdb | 190511 |  | accesión O31771 no registrada en el catálogo AFDB |
| O31834 | AFDB | data/metabolism/targets/afdb25k/O31834.pdb | 57023 |  | accesión O31834 no registrada en el catálogo AFDB |
| O31845 | AFDB | data/metabolism/targets/afdb25k/O31845.pdb | 122309 |  | accesión O31845 no registrada en el catálogo AFDB |
| O31859 | AFDB | data/metabolism/targets/afdb25k/O31859.pdb | 197963 |  | accesión O31859 no registrada en el catálogo AFDB |
| O32068 | AFDB | data/metabolism/targets/afdb25k/O32068.pdb | 158030 |  | accesión O32068 no registrada en el catálogo AFDB |
| O32078 | AFDB | data/metabolism/targets/afdb25k/O32078.pdb | 113885 |  | accesión O32078 no registrada en el catálogo AFDB |
| O32087 | AFDB | data/metabolism/targets/afdb25k/O32087.pdb | 52730 |  | accesión O32087 no registrada en el catálogo AFDB |
| O32103 | AFDB | data/metabolism/targets/afdb25k/O32103.pdb | 137618 |  | accesión O32103 no registrada en el catálogo AFDB |
| O32107 | AFDB | data/metabolism/targets/afdb25k/O32107.pdb | 103274 |  | accesión O32107 no registrada en el catálogo AFDB |
| O32152 | AFDB | data/metabolism/targets/afdb25k/O32152.pdb | 162485 |  | accesión O32152 no registrada en el catálogo AFDB |
| O32156 | AFDB | data/metabolism/targets/afdb25k/O32156.pdb | 280583 |  | accesión O32156 no registrada en el catálogo AFDB |
| O32161 | AFDB | data/metabolism/targets/afdb25k/O32161.pdb | 86912 |  | accesión O32161 no registrada en el catálogo AFDB |
| O32175 | AFDB | data/metabolism/targets/afdb25k/O32175.pdb | 83996 |  | accesión O32175 no registrada en el catálogo AFDB |
| O32207 | AFDB | data/metabolism/targets/afdb25k/O32207.pdb | 200231 |  | accesión O32207 no registrada en el catálogo AFDB |
| O32254 | AFDB | data/metabolism/targets/afdb25k/O32254.pdb | 217970 |  | accesión O32254 no registrada en el catálogo AFDB |
| O32256 | AFDB | data/metabolism/targets/afdb25k/O32256.pdb | 194156 |  | accesión O32256 no registrada en el catálogo AFDB |
| O32291 | AFDB | data/metabolism/targets/afdb25k/O32291.pdb | 230363 |  | accesión O32291 no registrada en el catálogo AFDB |
| O32965 | AFDB | data/metabolism/targets/afdb25k/O32965.pdb | 280745 |  | accesión O32965 no registrada en el catálogo AFDB |
| O33730 | AFDB | data/metabolism/targets/afdb25k/O33730.pdb | 182816 |  | accesión O33730 no registrada en el catálogo AFDB |
| O33752 | AFDB | data/metabolism/targets/afdb25k/O33752.pdb | 84320 |  | accesión O33752 no registrada en el catálogo AFDB |
| O34318 | AFDB | data/metabolism/targets/afdb25k/O34318.pdb | 79622 |  | accesión O34318 no registrada en el catálogo AFDB |
| O34319 | AFDB | data/metabolism/targets/afdb25k/O34319.pdb | 219428 |  | accesión O34319 no registrada en el catálogo AFDB |
| O34352 | AFDB | data/metabolism/targets/afdb25k/O34352.pdb | 256202 |  | accesión O34352 no registrada en el catálogo AFDB |
| O34382 | AFDB | data/metabolism/targets/afdb25k/O34382.pdb | 197477 |  | accesión O34382 no registrada en el catálogo AFDB |
| O34384 | AFDB | data/metabolism/targets/afdb25k/O34384.pdb | 124901 |  | accesión O34384 no registrada en el catálogo AFDB |
| O34416 | AFDB | data/metabolism/targets/afdb25k/O34416.pdb | 194318 |  | accesión O34416 no registrada en el catálogo AFDB |
| O34423 | AFDB | data/metabolism/targets/afdb25k/O34423.pdb | 184193 |  | accesión O34423 no registrada en el catálogo AFDB |
| O34424 | AFDB | data/metabolism/targets/afdb25k/O34424.pdb | 114938 |  | accesión O34424 no registrada en el catálogo AFDB |
| O34451 | AFDB | data/metabolism/targets/afdb25k/O34451.pdb | 220562 |  | accesión O34451 no registrada en el catálogo AFDB |
| O34458 | AFDB | data/metabolism/targets/afdb25k/O34458.pdb | 201041 |  | accesión O34458 no registrada en el catálogo AFDB |
| O34512 | AFDB | data/metabolism/targets/afdb25k/O34512.pdb | 344735 |  | accesión O34512 no registrada en el catálogo AFDB |
| O34524 | AFDB | data/metabolism/targets/afdb25k/O34524.pdb | 201689 |  | accesión O34524 no registrada en el catálogo AFDB |
| O34542 | AFDB | data/metabolism/targets/afdb25k/O34542.pdb | 146933 |  | accesión O34542 no registrada en el catálogo AFDB |
| O34574 | AFDB | data/metabolism/targets/afdb25k/O34574.pdb | 205577 |  | accesión O34574 no registrada en el catálogo AFDB |
| O34593 | AFDB | data/metabolism/targets/afdb25k/O34593.pdb | 83591 |  | accesión O34593 no registrada en el catálogo AFDB |
| O34611 | AFDB | data/metabolism/targets/afdb25k/O34611.pdb | 200717 |  | accesión O34611 no registrada en el catálogo AFDB |
| O34612 | AFDB | data/metabolism/targets/afdb25k/O34612.pdb | 110240 |  | accesión O34612 no registrada en el catálogo AFDB |
| O34643 | AFDB | data/metabolism/targets/afdb25k/O34643.pdb | 127574 |  | accesión O34643 no registrada en el catálogo AFDB |
| O34647 | AFDB | data/metabolism/targets/afdb25k/O34647.pdb | 78650 |  | accesión O34647 no registrada en el catálogo AFDB |
| O34648 | AFDB | data/metabolism/targets/afdb25k/O34648.pdb | 78326 |  | accesión O34648 no registrada en el catálogo AFDB |
| O34663 | AFDB | data/metabolism/targets/afdb25k/O34663.pdb | 242675 |  | accesión O34663 no registrada en el catálogo AFDB |
| O34719 | AFDB | data/metabolism/targets/afdb25k/O34719.pdb | 46979 |  | accesión O34719 no registrada en el catálogo AFDB |
| O34731 | AFDB | data/metabolism/targets/afdb25k/O34731.pdb | 166130 |  | accesión O34731 no registrada en el catálogo AFDB |
| O34733 | AFDB | data/metabolism/targets/afdb25k/O34733.pdb | 222587 |  | accesión O34733 no registrada en el catálogo AFDB |
| O34782 | AFDB | data/metabolism/targets/afdb25k/O34782.pdb | 166211 |  | accesión O34782 no registrada en el catálogo AFDB |
| O34805 | AFDB | data/metabolism/targets/afdb25k/O34805.pdb | 200393 |  | accesión O34805 no registrada en el catálogo AFDB |
| O34867 | AFDB | data/metabolism/targets/afdb25k/O34867.pdb | 65042 |  | accesión O34867 no registrada en el catálogo AFDB |
| O34886 | AFDB | data/metabolism/targets/afdb25k/O34886.pdb | 205739 |  | accesión O34886 no registrada en el catálogo AFDB |
| O34896 | AFDB | data/metabolism/targets/afdb25k/O34896.pdb | 173420 |  | accesión O34896 no registrada en el catálogo AFDB |
| O34897 | AFDB | data/metabolism/targets/afdb25k/O34897.pdb | 176741 |  | accesión O34897 no registrada en el catálogo AFDB |
| O34903 | AFDB | data/metabolism/targets/afdb25k/O34903.pdb | 153980 |  | accesión O34903 no registrada en el catálogo AFDB |
| O34908 | AFDB | data/metabolism/targets/afdb25k/O34908.pdb | 148391 |  | accesión O34908 no registrada en el catálogo AFDB |
| O34920 | AFDB | data/metabolism/targets/afdb25k/O34920.pdb | 214568 |  | accesión O34920 no registrada en el catálogo AFDB |
| O34948 | AFDB | data/metabolism/targets/afdb25k/O34948.pdb | 179900 |  | accesión O34948 no registrada en el catálogo AFDB |
| O34949 | AFDB | data/metabolism/targets/afdb25k/O34949.pdb | 107081 |  | accesión O34949 no registrada en el catálogo AFDB |
| O34969 | AFDB | data/metabolism/targets/afdb25k/O34969.pdb | 179090 |  | accesión O34969 no registrada en el catálogo AFDB |
| O34979 | AFDB | data/metabolism/targets/afdb25k/O34979.pdb | 150497 |  | accesión O34979 no registrada en el catálogo AFDB |
| O34984 | AFDB | data/metabolism/targets/afdb25k/O34984.pdb | 283823 |  | accesión O34984 no registrada en el catálogo AFDB |
| O35008 | AFDB | data/metabolism/targets/afdb25k/O35008.pdb | 215378 |  | accesión O35008 no registrada en el catálogo AFDB |
| O35014 | AFDB | data/metabolism/targets/afdb25k/O35014.pdb | 280988 |  | accesión O35014 no registrada en el catálogo AFDB |
| O35027 | AFDB | data/metabolism/targets/afdb25k/O35027.pdb | 89099 |  | accesión O35027 no registrada en el catálogo AFDB |
| O50310 | AFDB | data/metabolism/targets/afdb25k/O50310.pdb | 215378 |  | accesión O50310 no registrada en el catálogo AFDB |
| O51055 | AFDB | data/metabolism/targets/afdb25k/O51055.pdb | 257417 |  | accesión O51055 no registrada en el catálogo AFDB |
| O51468 | AFDB | data/metabolism/targets/afdb25k/O51468.pdb | 154385 |  | accesión O51468 no registrada en el catálogo AFDB |
| O53500 | AFDB | data/metabolism/targets/afdb25k/O53500.pdb | 178361 |  | accesión O53500 no registrada en el catálogo AFDB |
| O58549 | AFDB | data/metabolism/targets/afdb25k/O58549.pdb | 302048 |  | accesión O58549 no registrada en el catálogo AFDB |
| O58969 | AFDB | data/metabolism/targets/afdb25k/O58969.pdb | 286658 |  | accesión O58969 no registrada en el catálogo AFDB |
| O59672 | AFDB | data/metabolism/targets/afdb25k/O59672.pdb | 485351 |  | accesión O59672 no registrada en el catálogo AFDB |
| O59760 | AFDB | data/metabolism/targets/afdb25k/O59760.pdb | 156410 |  | accesión O59760 no registrada en el catálogo AFDB |
| O59832 | AFDB | data/metabolism/targets/afdb25k/O59832.pdb | 277181 |  | accesión O59832 no registrada en el catálogo AFDB |
| O66400 | AFDB | data/metabolism/targets/afdb25k/O66400.pdb | 120770 |  | accesión O66400 no registrada en el catálogo AFDB |
| O66404 | AFDB | data/metabolism/targets/afdb25k/O66404.pdb | 139238 |  | accesión O66404 no registrada en el catálogo AFDB |
| O66405 | AFDB | data/metabolism/targets/afdb25k/O66405.pdb | 225827 |  | accesión O66405 no registrada en el catálogo AFDB |
| O66407 | AFDB | data/metabolism/targets/afdb25k/O66407.pdb | 80594 |  | accesión O66407 no registrada en el catálogo AFDB |
| O66408 | AFDB | data/metabolism/targets/afdb25k/O66408.pdb | 220319 |  | accesión O66408 no registrada en el catálogo AFDB |
| O66409 | AFDB | data/metabolism/targets/afdb25k/O66409.pdb | 144584 |  | accesión O66409 no registrada en el catálogo AFDB |
| O66410 | AFDB | data/metabolism/targets/afdb25k/O66410.pdb | 80027 |  | accesión O66410 no registrada en el catálogo AFDB |
| O66411 | AFDB | data/metabolism/targets/afdb25k/O66411.pdb | 222425 |  | accesión O66411 no registrada en el catálogo AFDB |
| O66412 | AFDB | data/metabolism/targets/afdb25k/O66412.pdb | 153089 |  | accesión O66412 no registrada en el catálogo AFDB |
| O66413 | AFDB | data/metabolism/targets/afdb25k/O66413.pdb | 117935 |  | accesión O66413 no registrada en el catálogo AFDB |
| O66414 | AFDB | data/metabolism/targets/afdb25k/O66414.pdb | 185489 |  | accesión O66414 no registrada en el catálogo AFDB |
| O66416 | AFDB | data/metabolism/targets/afdb25k/O66416.pdb | 147581 |  | accesión O66416 no registrada en el catálogo AFDB |
| O66417 | AFDB | data/metabolism/targets/afdb25k/O66417.pdb | 134135 |  | accesión O66417 no registrada en el catálogo AFDB |
| O66418 | AFDB | data/metabolism/targets/afdb25k/O66418.pdb | 176255 |  | accesión O66418 no registrada en el catálogo AFDB |
| O66419 | AFDB | data/metabolism/targets/afdb25k/O66419.pdb | 105218 |  | accesión O66419 no registrada en el catálogo AFDB |
| O66423 | AFDB | data/metabolism/targets/afdb25k/O66423.pdb | 219509 |  | accesión O66423 no registrada en el catálogo AFDB |
| O66425 | AFDB | data/metabolism/targets/afdb25k/O66425.pdb | 101816 |  | accesión O66425 no registrada en el catálogo AFDB |
| O66427 | AFDB | data/metabolism/targets/afdb25k/O66427.pdb | 162890 |  | accesión O66427 no registrada en el catálogo AFDB |
| O66466 | AFDB | data/metabolism/targets/afdb25k/O66466.pdb | 143126 |  | accesión O66466 no registrada en el catálogo AFDB |
| O66516 | AFDB | data/metabolism/targets/afdb25k/O66516.pdb | 229148 |  | accesión O66516 no registrada en el catálogo AFDB |
| O66525 | AFDB | data/metabolism/targets/afdb25k/O66525.pdb | 89099 |  | accesión O66525 no registrada en el catálogo AFDB |
| O66536 | AFDB | data/metabolism/targets/afdb25k/O66536.pdb | 82457 |  | accesión O66536 no registrada en el catálogo AFDB |
| O66558 | AFDB | data/metabolism/targets/afdb25k/O66558.pdb | 170747 |  | accesión O66558 no registrada en el catálogo AFDB |

_992 filas adicionales en el TSV._

## Rutas no coincidentes

Sin hallazgos.

## Archivos faltantes

Sin hallazgos.

## Límites

- `plddt` se normaliza a 0-1. En AFDB se calcula como la media de B-factor de los CA del primer modelo; es una aproximación, no el valor oficial por residuo.
- `ptm` solo existe en el catálogo ESMFold; en AFDB queda vacío.
- `length` en AFDB se cuenta desde los residuos del PDB, no desde la secuencia UniProt.
- La deduplicación es por accesión y por SHA-256 del archivo PDB; no detecta la misma proteína con coordenadas distintas ni duplicados por secuencia.
- Solo se inspeccionan `*.pdb` en los directorios configurados; no hay mmCIF, ensamblados biológicos ni validación geométrica.
- No se borra ni se mueve ningún archivo, tampoco los huérfanos.
- El script no usa red.
