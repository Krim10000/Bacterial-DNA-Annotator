# Análisis set_004

## Integridad

| métrica | valor |
|---|---:|
| Proteínas | **127** |
| PDB | 127 |
| Errores | 0 |
| PDB con NaN | 0 |
| PDB sin ATOM | 0 |

## Tiempo y longitud

- Tiempo: **22,5 min**.
- Longitud: **187–209 aa**.
- Media: 199 aa.
- pLDDT medio: **0,802**.
- pTM medio: **0,819**.

## Calidad

| pLDDT | proteínas |
|---|---:|
| ≥0,80 | 86 |
| 0,70–0,80 | 28 |
| 0,60–0,70 | 8 |
| 0,50–0,60 | 3 |
| <0,50 | 2 |

## Acumulado

- **1.469** cortas completadas.
- **1.101** restantes.
- 57,2% de las 2.570.
- **4.338 clusters**.
- **9.070 proteínas**.
- **363 genomas**.

## Hallazgos

- 127/127 tienen hit no-self.
- **110/127** de alta confianza.

Categorías de los 110:

| categoría | cantidad |
|---|---:|
| otros | 43 |
| uncharacterized/DUF | 39 |
| oxidorreductasa | 8 |
| RNA metiltransferasa | 7 |
| peptidasa | 4 |
| regulador | 4 |
| transporte | 4 |
| aminotransferasa | 1 |

## Casos destacados

| query | hit | identidad | lectura |
|---|---|---:|---|
| P9WMC6 | regulador HTH Mb2936c | 72% | posible regulador |
| P9WMD6 | regulador HTH Mb0790c | 69% | posible regulador |
| P9WLE6 | Rv2466c | 50% | tiorredoxina reductasa |
| Q49YS0 | deglicasa 2 | 48% | posible deglicasa |
| P0ADI8 | YecD | 44% | isocorismatasa |
| P37872 | MJ0882 | 43% | SAM metiltransferasa |
| A6TAY4 | NudL | 43% | Nudix hidrolasa |
| Q09677 | MdaB | 42% | NADPH:quinona oxidorreductasa |
| P42398 | YciA | 40% | tioesterasa |
| P44137 | TruD | 40% | pseudouridina sintasa |

## Sobre subunidades

- Las proteínas cortas pueden ser subunidades.
- También pueden ser dominios pequeños, reguladores o proteínas independientes.
- El tamaño no prueba cuaternario.
- La evaluación de complejos sigue congelada hasta terminar los monómeros.
- Cuando toque: vecindad + pares curados + AlphaFold-Multimer/ipTM.

## Lectura

- set_004 es de mejor calidad que los sets 001–003.
- Tiene reguladores HTH y enzimas pequeñas interesantes.
- Las identidades son moderadas: familia, no función.
