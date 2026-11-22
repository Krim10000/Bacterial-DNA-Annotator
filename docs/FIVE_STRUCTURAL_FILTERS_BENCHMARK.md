# Benchmark de los 5 filtros

## Resultado

- 10 estructuras: **213,8 s** (~3,6 min).
- Promedio: **21,4 s/estructura**.
- En frío: ~26 s/estructura.
- Caché: 2 estructuras repitieron y bajaron a 0 s.
- CPU local; sin GPU; sin créditos.

## Tabla

Ver `data/metabolism/targets/FIVE_FILTERS_BENCHMARK.tsv`.

| # | PDB | F1 | F2 | F3 | F4 | F5 | pasan | prioridad |
|---|---|---|---:|---:|---:|---:|---:|---|
| 1 | test O06745 | sí | no | sí | sí | no | 3 | review |
| 2 | P77161 exacto | sí | no | sí | sí | sí | 4 | review |
| 3 | Q9KBE3 | sí | sí | sí | sí | no | 4 | review |
| 4 | Q12726 | sí | sí | sí | no | no | 3 | review |
| 5 | Q8VPD4 | sí | no | no | sí | no | 2 | review |
| 6 | 610 | sí | no | sí | sí | sí | 4 | review |
| 7 | 2055 | sí | no | sí | sí | sí | 4 | review |
| 8 | 2476 | sí | no | sí | sí | sí | 4 | review |
| 9 | 3962 | sí | no | sí | sí | sí | 4 | review |
| 10 | 1LYZ control | sí | sí | sí | sí | no | 4 | review |

- **high: 0**
- **review: 10**
- **low: 0**
- No hubo un caso con F2 + F5 que justificara `high`.

## Extrapolación

| supuesto | tiempo |
|---|---:|
| 1.000 proteínas, 1 núcleo | ~5,9 h |
| 647.588 proteínas, 1 núcleo | **~3.850 h** (~160 días) |
| 647.588 proteínas, 8 núcleos | ~481 h |
| 647.588 proteínas, 32 núcleos | ~120 h |
| 10 h, 1 núcleo | ~1.680 proteínas |
| 10 h, 8 núcleos | ~13.450 proteínas |

## Decisión

**No se inicia el barrido de 647.588 proteínas.**

Motivos:

1. Supera ampliamente 10 h incluso con 32 núcleos.
2. La mayoría de esas proteínas **no tienen PDB**; habría que plegarlas antes.
3. Plegar 647k con ESMFold/Neurosnap es inviable en tiempo y coste.
4. Los filtros deben aplicarse a candidatos priorizados, no a todo el proteoma.

## Estrategia

- Aplicar a los **9/10 ya estructurados**.
- Luego a los **5 estructurales restantes**.
- Luego a **17 P1/P2**.
- Luego a **51 P3** o una muestra.
- No lanzar 647k.

## Cómo reproducir

```bash
python scripts/structure/12_structural_filters.py \
  --pdb candidato.pdb --genome GCF_... --uniprot P77161 \
  --out resultados.tsv
```

GUI: **Ejecutar 5 filtros**.
Pipeline: `--structural-filters`.
