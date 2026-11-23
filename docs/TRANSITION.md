# Índice de traspaso

Documentación para retomar cada proceso sin depender del contexto de una
conversación. Cada entrada declara qué está implementado, qué falta y dónde
está la evidencia.

| proceso | documento principal | estado | siguiente |
|---|---|---|---|
| Visión y alcance | [`ROADMAP.md`](ROADMAP.md) / [`ROADMAP.en.md`](ROADMAP.en.md) | vivo | mantener estados honestos |
| Paridad 400 vs FASTA nuevo | [`PIPELINE_PARITY.md`](PIPELINE_PARITY.md) | auditoría 2026-11-23 | regulación/operones |
| oriC (ORCA/DoriC) | [`ORIC.md`](ORIC.md) | núcleo implementado | rotación circular, HMM RPA/Rep |
| Tres caminos A/B/C | [`THREE_PROTEIN_PATHS.md`](THREE_PROTEIN_PATHS.md) | A automático; B/C clasificados | estructura B/C automática |
| Reducción de desconocidas | [`UNKNOWN_REDUCTION_PLAN.md`](UNKNOWN_REDUCTION_PLAN.md) | filtros baratos + gating | categorías amplias |
| Filtros estructurales | [`FILTER_STRATEGY.md`](FILTER_STRATEGY.md) | pLDDT/gating implementado | colas separadas |
| Represores/alosteria | [`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md) | detectores; dirección bloqueada | presupuesto GPU/banco |
| Proteínas reguladoras | [`PROTEINAS_REGULADORAS.md`](PROTEINAS_REGULADORAS.md) | lecciones y límites | reabrir con datos |
| Integraciones | [`INTEGRATIONS.md`](INTEGRATIONS.md) | matriz y prioridades | Bakta/eggNOG/AlphaFold DB |
| GPU y Kaggle | [`GPU_STATUS.md`](GPU_STATUS.md) | política y límites | presupuesto |
| Benchmark | [`BENCHMARK.md`](BENCHMARK.md) | 12+ proyectos | benchmark R1 |
| Instalación | [`INSTALL.md`](INSTALL.md) / [`INSTALL_DOCKER.md`](INSTALL_DOCKER.md) | Ubuntu/Docker | macOS/WSL2 |
| Módulos | [`MODULES.md`](MODULES.md) | 10 módulos | `oric` como módulo |

## Regla de transición

Si un cambio altera una capacidad:

1. actualiza el documento del proceso;
2. actualiza [`PIPELINE_PARITY.md`](PIPELINE_PARITY.md) si cambia la paridad
   400 vs FASTA nuevo;
3. añade o ajusta el test correspondiente;
4. registra la fuente, versión y límite en el README cuando sea visible.

## Notas históricas que no se publican

- `data/orca_predictions/ORCA_FIX.md`
- `data/orca_predictions/ORCA_QUALITY.md`
- `data/orca_predictions/ORIGIN_PIPELINE.md`
- `data/orca_predictions_fix/ORCA_TRIAGE.md`

Son evidencia de la campaña de los 400. Sus conclusiones vigentes están
resumidas en [`ORIC.md`](ORIC.md); no se incluyen en el release porque `data/`
no se distribuye.
