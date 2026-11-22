# Dónde se equivoca el motor de operones

Analiza 2897 decisiones de frontera contra RegulonDB.

| decisión | n |
|---|---:|
| acierto al unir | 1409 |
| acierto al cortar | 1044 |
| **FN: cortó de más** | **370** |
| **FP: unió de más** | **74** |

## FN por regla (cortes que no debían hacerse)

| regla | n | % de los FN |
|---|---:|---:|
| `R4_distancia` | 295 | 79.7% |
| `R2_terminador` | 75 | 20.3% |

## FP por regla (fusiones que no debían hacerse)

| regla | n | % de los FP |
|---|---:|---:|
| `R5_defecto` | 70 | 94.6% |
| `R1_atenuador` | 4 | 5.4% |

## FN según tamaño del operón real

| genes en el operón | FN | aciertos | % de error |
|---|---:|---:|---:|
| 2 | 82 | 317 | 25.9% |
| 3 | 82 | 234 | 35.0% |
| 4 | 61 | 212 | 28.8% |
| 5 | 45 | 194 | 23.2% |
| 6 | 29 | 143 | 20.3% |
| 7 | 18 | 64 | 28.1% |
| 8 | 15 | 33 | 45.5% |
| 9 | 11 | 41 | 26.8% |
| 10 | 4 | 11 | 36.4% |
| 11 | 9 | 41 | 22.0% |
| 12 | 6 | 69 | 8.7% |
| 13 | 3 | 9 | 33.3% |
| 14 | 2 | 11 | 18.2% |
| 17 | 3 | 12 | 25.0% |

## Ejemplos de FN

- `R4_distancia` · b0014 → b0015: hueco 88 pb >= 42.5
- `R4_distancia` · b0027 → b0028: hueco 124 pb >= 42.5
- `R4_distancia` · b0036 → b0037: hueco 108 pb >= 42.5
- `R4_distancia` · b0038 → b0039: hueco 128 pb >= 42.5
- `R2_terminador` · b0019 → b0020: horquilla A a 15 pb, poli-U 4
- `R2_terminador` · b0037 → b0038: horquilla A a 45 pb, poli-U 5
- `R2_terminador` · b0042 → b0043: horquilla A a 7 pb, poli-U 4
- `R2_terminador` · b0053 → b0054: horquilla A a 18 pb, poli-U 5

## Ejemplos de FP

- `R5_defecto` · b0138 → b0139: hueco 16 pb, sin señal fuerte
- `R5_defecto` · b0146 → b0147: hueco 14 pb, sin señal fuerte
- `R5_defecto` · b0175 → b0176: hueco 11 pb, sin señal fuerte
- `R5_defecto` · b0176 → b0177: hueco 29 pb, sin señal fuerte
- `R1_atenuador` · b0553 → b0552: gen previo 19 aa (péptido líder); horquilla A_fuerte detrás = atenuador
- `R1_atenuador` · b4803 → b0694: gen previo 34 aa (péptido líder); horquilla B_media detrás = atenuador
- `R1_atenuador` · b4737 → b1015: gen previo 8 aa (péptido líder); horquilla B_media detrás = atenuador
- `R1_atenuador` · b4665 → b2913: gen previo 19 aa (péptido líder); horquilla B_media detrás = atenuador
