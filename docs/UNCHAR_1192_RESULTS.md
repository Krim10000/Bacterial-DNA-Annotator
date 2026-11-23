# Resultados de 1.192 accesiones con AFDB

## Cobertura

- 25.058 clusters.
- 3.931 accesiones únicas.
- 1.192 con AFDB descargable.
- 2.739 sin AFDB (404).
- 15.633 clusters cubiertos por las 1.192.
- 9.425 clusters sin modelo.

## Resultados

| filtro | pasa |
|---|---:|
| F1 Foldseek local | 1.192 |
| F2 EnzyMM/M-CSA | 195 |
| F3 P2Rank | 1.052 |
| F4 AlphaFill | 774 |
| F5 propio | 0 |
| F2+F3+F4 | **176** |

| prioridad | proteínas |
|---|---:|
| high | 0 |
| review | 1.110 |
| low | 82 |
| errores | 0 |

## Propagación

| prioridad | clusters |
|---|---:|
| review | 15.191 |
| low | 442 |
| con F2+F3+F4 | **3.374** |

## Tiempo

- Filtros: **6.314,5 s** (~1 h 45 min).
- Foldseek local: ~20 min adicionales.
- Total: ~2 h.

## Ejemplos

| proteína | plantilla | cofactor | orientación |
|---|---|---|---|
| P9WN14 | 1gro | Ca/Mg | glucósido hidrolasa |
| P0DUF3 | 3exr, 1aqj | MTA/SAH | metiltransferasa |
| P76585 | 1nww, 4pmh | Ca | posible enzima |
| P37443 | 4bzc, 2p18 | Zn/Mn | posible enzima |
| P71353 | 2om2 | ATP/VGL | posible ATPasa/ligasa |
| P9WIY6 | 1d4e | CMP/K | posible nucleótido-enzima |

## Límites

- F5 = 0: no hay mapas propios de sitio activo.
- F2 no confirma sustrato.
- F3 bolsillo no es sitio activo.
- F4 cofactor no define reacción.
- No hay cierres.
- 9.425 clusters siguen sin modelo AFDB.

## Archivos

- `UNCHAR_AFDB_PRESENT_filters.tsv`
- `UNCHAR_AFDB_PRESENT_filters.json`
- `UNCHAR_25K_RESULTS.tsv`
- `UNCHAR_UNIQUE_3931.tsv`
