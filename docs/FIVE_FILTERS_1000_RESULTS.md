# Resultados 1.000 proteínas poco documentadas

## Método

- 1.000 accesiones únicas, reparto por genomas.
- Modelos AFDB.
- F1 Foldseek local (`data/db/foldseek_swissprot`).
- F2 EnzyMM/M-CSA.
- F3 P2Rank.
- F4 AlphaFill.
- F5 sitio activo propio.

## Resultados

| filtro | pasa |
|---|---:|
| F1 | 1.000 |
| F2 | 180 |
| F3 | 890 |
| F4 | 680 |
| F5 | 0 |
| F2+F3+F4 | 157 |

| prioridad | cantidad |
|---|---:|
| high | 0 |
| review | 940 |
| low | 60 |
| errores | 0 |

## Tiempos

- Foldseek local, 1.000 vs 590.183: ~15 min.
- Filtros con caché: 431,6 s.
- Total local: ~22 min.

## Ejemplos

- P9WGQ6: oxidorreductasa probable; `3w8f`, `4k6f`, `1x7g`; MG/NA.
- Q57964: plantilla `1oie`; MG/UDP.
- P94396: hidrolasa probable; `3om8`, `1agy`, `1jkm`.
- P31447: sulfatasa probable; `4cyr`, `4mhx`.

## Límites

- F1 local detecta homología estructural; no confirma función.
- F2 da plantilla/EC; no confirma sustrato.
- F3 bolsillo no es sitio activo.
- F4 cofactor no define reacción.
- F5 = 0 porque no hay mapa propio para estas proteínas.
- Ninguna brecha cerrada.
- Ningún `high`.

## Archivos

- `FIVE_FILTERS_1000.tsv`
- `FIVE_FILTERS_1000.json`
- `FIVE_FILTERS_1000.foldseek.m8`
- `FIVE_FILTERS_1000_REVIEW.tsv`
