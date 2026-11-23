# Visión, estado real y próximos pasos

🌐 **Español** · [English](ROADMAP.en.md)

## La idea

Desde una secuencia de ADN, el usuario debería poder responder, con evidencia y
límites explícitos:

```text
ADN
 → especie y taxonomía
 → genes y proteínas
 → función probable (UniProt/EC/Pfam)
 → rutas metabólicas y brechas
 → operones y promotores
 → reguladores y bucles de regulación
 → hipótesis estructurales
```

No es un “anotador que lo sabe todo”. Es un sistema que **reduce incertidumbre
por capas** y dice qué está medido, qué es inferencia y qué no se puede concluir.

## Qué está implementado

| capacidad | estado | evidencia |
|---|---|---|
| Especie/taxonomía | ✅ | ANI, GTDB/NCBI, confianza y conflictos |
| CDS, tRNA, rRNA | ✅ | pyrodigal, tRNAscan-SE, barrnap |
| Proteínas y ORFs no anotados | ✅ | `proteins.faa`, motor ghost |
| Función por homología | ✅ | DIAMOND + Swiss-Prot, EC, UniProt |
| Rutas y brechas metabólicas | ✅ | Rhea, transporte, cierres, evidencia |
| Operones y promotores | ✅ | reglas, TSS/TTS, PWMs, RegulonDB como validación |
| Navegador clicable | ✅ | ficha por gen/operón/brecha |
| Módulos y GUI ES/EN/ZH | ✅ | 10 módulos, autoprueba, resumen |
| Docker y DOI | ✅ | Compose, Zenodo `10.5281/zenodo.24469732` |

## Qué está parcial

### Represores y bucles de regulación

Implementado:

- pipeline genoma → candidatos a represor;
- familia por dominio, homólogo por DIAMOND y efector heredado;
- sitios PWM y puntuación relativa;
- recall 92 % en *E. coli* y 98,7 % de sitios frente al azar.

No implementado:

- dirección del efector (inductor/correpresor);
- magnitud del efecto (ΔΔG, afinidad absoluta);
- validación de que el bucle se cierre en un contexto real.

Motivo: requiere GPU (MM/GBSA, ABFE, co-folding) y, sobre todo, datos
experimentales (EMSA/ITC). El criterio de reapertura está en
[`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md): primero MM/GBSA de LacI,
PurR y TrpR; si no ordena, no se gasta más.

### Reducción de proteínas desconocidas en los 400 genomas

Implementado:

- deduplicación, clustering y manifiestos;
- filtros baratos antes de estructura;
- F1 Foldseek local;
- F2 EnzyMM/M-CSA;
- F3 P2Rank y F4 AlphaFill;
- propagación a clusters y priorización;
- **gating por pLDDT/pTM y early-exit** en `structural_filters.py`.

Pendiente:

- separar las colas baratas y caras en comandos independientes;
- completar las 2.570 estructuras ≤600 aa y las 169 >600 aa;
- asignar categorías amplias (`enzyme_hypothesis`, `transporter_hypothesis`,
  `regulator_hypothesis`, `unknown`) de forma automática;
- revisión manual y propagación final.

La regla es: **ningún filtro cierra una función**. Solo reduce candidatos.

### GPU: Kaggle como primer filtro

Kaggle se usa como **primer filtro de GPU**, no como dependencia de la app:

- permite probar ESMFold/co-folding y modelos grandes sin GPU local;
- las secuencias se suben manualmente y solo si el usuario acepta;
- los resultados se importan como tablas/estructuras;
- no se ejecuta en el arranque ni en la GUI;
- el pipeline local sigue funcionando en CPU.

Esto es deliberado: Kaggle no es reproducible como servicio permanente y no
debe ser un requisito para anotar un genoma.

## Próximamente

| prioridad | pieza | coste | depende de |
|---|---|---|---|
| 1 | `protein_query`: proteína → homólogos/función | 1–2 días | Swiss-Prot + DIAMOND |
| 2 | Reducción GPU-free completa de desconocidas | 3–5 días | tablas actuales |
| 3 | Categorías amplias automáticas | 2–3 días | filtros + contexto |
| 4 | Reapertura GPU de dirección de efectores | presupuesto | Kaggle/Modal + banco |
| 5 | Instalador macOS / guía WSL2 | 1–2 días | probar en limpio |
| 6 | PyPI / TestPyPI | 1 día | empaquetado |
| 7 | JOSS / bio.tools | 1–3 días | DOI ya disponible |
| 8 | Benchmark reproducible vs. Bakta/Prokka/DFAST | 2–4 días | 5–6 genomas |
| 9 | Exportar borrador SBML (CarveMe/COBRApy) | 1 semana | mapeo Rhea |

## Benchmark

La comparación con 12+ proyectos (Prokka, Bakta, DRAM, CarveMe, gapseq,
Operon-mapper, BV-BRC, etc.) está en [`BENCHMARK.md`](BENCHMARK.md).
Conclusiones: la anotación base es tabla stakes; nuestra ventaja es la
integración local con evidencia y la GUI ES/EN/ZH; la brecha principal es
base funcional curada, empaquetado y benchmark público.

## Lo que no prometemos

- Una brecha no prueba que falte una enzima.
- Un dominio HTH no confirma un regulón.
- ESMFold no confirma función.
- AlphaFill no confirma sustrato.
- P2Rank no confirma sitio activo.
- Ningún resultado es `confirmed` sin experimento, EC o contexto fuerte.
- El modelo de ADN custom es experimental y no está evaluado.

## Documentos relacionados

- [`UNKNOWN_REDUCTION_PLAN.md`](UNKNOWN_REDUCTION_PLAN.md)
- [`FILTER_STRATEGY.md`](FILTER_STRATEGY.md)
- [`ALOSTERIA_V1_CIERRE.md`](ALOSTERIA_V1_CIERRE.md)
- [`PROTEINAS_REGULADORAS.md`](PROTEINAS_REGULADORAS.md)
- [`PIPELINE_GENOMA.md`](PIPELINE_GENOMA.md)
- [`GPU_STATUS.md`](GPU_STATUS.md)
- [`BENCHMARK.md`](BENCHMARK.md)
- [`INTEGRATIONS.md`](INTEGRATIONS.md)
- [`ORIC.md`](ORIC.md)
- [`TRANSITION.md`](TRANSITION.md)
- [`PIPELINE_PARITY.md`](PIPELINE_PARITY.md)
- [`COMPLEX_GAPS.md`](COMPLEX_GAPS.md)
- [`COMPLEX_CPU.md`](COMPLEX_CPU.md)
- [`COMPLEX_TRANSFER.md`](COMPLEX_TRANSFER.md)
