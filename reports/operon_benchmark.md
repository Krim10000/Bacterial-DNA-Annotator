# Benchmark de operones por reglas

**Generado:** 2026-11-24T02:28:59Z

**Alcance:** solo datos locales; sin red, sin reentrenar y sin modificar pipeline/GUI. No se inventa gold.

**Comando:** `python tools/benchmark_operons.py --gold data/regulondb/gt/gt_operons.tsv --pred data/operon_finder/ecoli_k12_operons_rules.tsv --genome ecoli_k12 --json reports/operon_benchmark.json --md reports/operon_benchmark.md`

**Genoma:** `ecoli_k12`

## Entradas

| artefacto | formato | ruta |
|---|---|---|
| gold | tsv | `data/regulondb/gt/gt_operons.tsv` |
| predicción | tsv | `data/operon_finder/ecoli_k12_operons_rules.tsv` |

**Estado:** evaluado (`evaluation_level=full`).

## Conteos

| campo | valor |
|---|---:|
| `gt_rows` | 2605 |
| `pred_rows` | 2490 |
| `gt_operons_filtered` | 827 |
| `pred_operons_filtered` | 847 |
| `pred_operons_with_mapped_genes` | 815 |
| `gt_genes` | 4509 |
| `pred_genes` | 4308 |
| `pred_genes_filtered` | 2665 |
| `shared_genes` | 4206 |
| `shared_genes_filtered` | 2599 |
| `pred_gene_mapping_coverage` | 0.9752 |
| `gt_without_coords` | 0 |
| `pred_without_coords` | 0 |
| `aliased_pred_genes` | 0 |

## Métricas por hebra

| hebra | nivel gen P | R | F1 | exact-set | nivel intervalo P | R | F1 | exact-intervalo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `all` | 0.828 | 0.770 | 0.798 | 0.561 | 0.761 | 0.780 | 0.771 | 0.559 |
| `+` | 0.868 | 0.740 | 0.799 | 0.574 | 0.754 | 0.795 | 0.774 | 0.576 |
| `-` | 0.796 | 0.799 | 0.797 | 0.549 | 0.769 | 0.765 | 0.767 | 0.542 |

### Hebra `all`

| métrica | P | R | F1 | TP | n_pred | n_gt | estado |
|---|---:|---:|---:|---:|---:|---:|---|
| Pares de genes | 0.828 | 0.770 | 0.798 | 3618 | 4368 | 4699 | evaluated |
| Conjuntos exactos | 0.569 | 0.561 | 0.565 | 464 | 815 | 827 | evaluated |
| Intervalos solapados | 0.761 | 0.780 | 0.771 | 645 | 847 | 827 | evaluated |
| Intervalos exactos | 0.545 | 0.559 | 0.552 | 462 | 847 | 827 | evaluated |

- exact-set: 464 de 827 operones gold (56.1%).
- exact-intervalo: 462 de 827 operones gold (55.9%).

### Hebra `+`

| métrica | P | R | F1 | TP | n_pred | n_gt | estado |
|---|---:|---:|---:|---:|---:|---:|---|
| Pares de genes | 0.868 | 0.740 | 0.799 | 1706 | 1966 | 2306 | evaluated |
| Conjuntos exactos | 0.565 | 0.574 | 0.569 | 230 | 407 | 401 | evaluated |
| Intervalos solapados | 0.754 | 0.795 | 0.774 | 319 | 423 | 401 | evaluated |
| Intervalos exactos | 0.546 | 0.576 | 0.561 | 231 | 423 | 401 | evaluated |

- exact-set: 230 de 401 operones gold (57.4%).
- exact-intervalo: 231 de 401 operones gold (57.6%).

### Hebra `-`

| métrica | P | R | F1 | TP | n_pred | n_gt | estado |
|---|---:|---:|---:|---:|---:|---:|---|
| Pares de genes | 0.796 | 0.799 | 0.797 | 1912 | 2402 | 2393 | evaluated |
| Conjuntos exactos | 0.574 | 0.549 | 0.561 | 234 | 408 | 426 | evaluated |
| Intervalos solapados | 0.769 | 0.765 | 0.767 | 326 | 424 | 426 | evaluated |
| Intervalos exactos | 0.545 | 0.542 | 0.543 | 231 | 424 | 426 | evaluated |

- exact-set: 234 de 426 operones gold (54.9%).
- exact-intervalo: 231 de 426 operones gold (54.2%).

## Avisos

- 21 filas con mapping_ok distinto de yes excluidas del gold evaluable.
- 3 filas sin lista de genes ignoradas.
- 66 genes de la predicción no están en el universo gold y se excluyen de las métricas por gen.
- El gold no declara contig en algunas filas: el contig se trata como comodín en el emparejamiento de intervalos.

## Promotores

- **Estado:** `not_evaluable`.
- **Razón:** Este benchmark mide estructura de operones. No se compara promoters_auto con gt_promoters porque una ventana upstream no es un TSS y no se define aquí un emparejamiento TSS→ventana; no se inventa gold.

## Límites

- No es un holdout ciego: las reglas de K-12 pueden haberse calibrado con el mismo gold; mide consistencia, no generalización.
- Nivel gen: los pares solo existen dentro de operones con ≥ 2 genes; IDs no mapeados se excluyen y eso puede favorecer la precisión.
- Nivel intervalo: solapamiento recíproco mínimo 0.5; coordenadas exactas son una métrica aparte y muy sensible a los bordes.
- Si el gold no declara contig, el contig se trata como comodín. Hebras `+` y `-` se evalúan por separado.
- Promotores/TSS no se evalúan aquí: una ventana upstream no es un TSS y no se define un emparejamiento TSS→ventana.
- Un `not_evaluable` no es cero: indica que falta gold, predicción o IDs/coordenadas comparables.
