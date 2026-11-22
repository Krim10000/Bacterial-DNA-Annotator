# Corrección ORCA multi-configuración

Se reanalizaron los **95 genomas con confianza <0,5** con 4 configuraciones
(ambos marcadores, solo dnaA, solo dnaN, sin marcadores).

## Resultado

| fuente | genomas |
|---|---|
| DoriC-anchored | 48 |
| ORCA-consensus (dnaA) | 31 |
| ORCA (revisar) | 16 |

- **79/95 resueltos** (DoriC o consenso dnaA).
- De los 67 con DoriC: **63 a ≤50 kb**, mediana 0.
- 16 quedan marcados para revisión (sin DoriC y sin candidato cerca de dnaA).

## Hallazgo

Inyectar `dnaA`+`dnaN` juntos puede **empeorar** la predicción (caso 0,02).
Con **solo dnaA** subió a 0,42 y el error bajó de 1,3 Mb a 19 kb.

## Criterio de selección

1. Candidato a ≤50 kb de dnaA: elegir el de mayor confianza.
2. Si no, y hay DoriC: anclar a DoriC (etiquetado).
3. Si no: elegir el de mayor confianza y marcar "revisar".

## Próximo

- GC skew + densidad de cajas dnaA para los 16 restantes.
- Integrar la posición corregida al GFF, arquitectura y GUI.
