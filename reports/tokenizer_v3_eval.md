# Evaluación del tokenizador v3 vs v2

Generado: 2026-11-20T02:42:14.598757+00:00

- Genomas: 6 (4 no vistos + 2 vistos)

## Comparación global

| métrica | v2 | v3 |
|---|---:|---:|
| vocab | 4096 | 256 |
| nt/token | 4.004 | 2.969 |
| unk rate | 0.00000 | 0.00000 |
| roundtrip | 120/120 | 120/120 |
| tokens 1 kb (media) | 252.5 | 337.9 |
| tokens 1,386 Gbp | 346,134,714 | 466,778,782 |
| IUPAC completo | False | True |
| motivos 1 token | 7/7 | 7/7 |
| piezas ▁ | 2 | 0 |

## Motivos (token único y en contexto)

| motivo | v2 single | v2 contexto | v3 single | v3 contexto |
|---|---:|---:|---:|---:|
| TATAAT | True | 1.000 | True | 1.000 |
| TTGACA | True | 1.000 | True | 1.000 |
| AGGAGG | True | 0.696 | True | 0.703 |
| ATG | True | 0.989 | True | 0.993 |
| TAA | True | 1.000 | True | 1.000 |
| TAG | True | 1.000 | True | 1.000 |
| TGA | True | 0.978 | True | 0.658 |

## Motivos solapados

La protección atómica tiene un límite: si dos motivos comparten bases, el tokenizador elige uno y parte el otro. Casos medidos:

- `TGA` comparte la última base con `ATG`, `GTG` y `TTG`; en `ATGA` gana `ATG` y `TGA` queda partido. Por eso su tasa en contexto baja a ~0,66.
- `AGGAGG` contiene `GGAGG`; al proteger ambos, algunas apariciones se resuelven como `A` + `GGAGG`. La tasa ronda 0,70.
- `TATAAT` y `TTGACA` no se solapan con otros motivos protegidos: tasa 1,00.
Para el modelo nuevo: definir una **prioridad de solapamiento** o un pretokenizador que corte por límites de codón, en vez de confiar solo en `user_defined_symbols`.

## IUPAC

- v2: ['▁', 'ACGT', 'R', 'Y', 'K', 'M', 'S', 'W', 'BDHV', 'N', 'N'] (unk=1)
- v3: ['ACG', 'T', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', 'N'] (unk=0)

## Tokens especiales

| token | v2 | v3 |
|---|---:|---:|
| pad | -1 | 0 |
| unk | 0 | 1 |
| bos | 1 | 2 |
| eos | 2 | 3 |
| mask | 0 | 4 |

## Veredicto

- ✅ compresión 2-3 nt/token: 2.97
- ✅ IUPAC completo: True
- ✅ 7/7 motivos token único: 7/7
- ✅ sin ▁: 0
- ✅ special tokens distintos: True
- ✅ roundtrip: 120/120

Nota: es una evaluación **intrínseca** del tokenizador. El impacto en el modelo requiere entrenar un modelo con v3.