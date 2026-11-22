# Anotación de TSS y TTS por regiones candidatas

Pipeline genérica para anotar **inicios y finales de transcripción** en cualquier
genoma bacteriano, incorporando desde el diseño la clasificación de los casos
anómalos.

## Por qué regiones y no el genoma entero

Un promotor o un terminador no puede estar dentro de un gen, así que analizar
los 4,6 Mb completos desperdicia ~80% del trabajo. Pero recortar el espacio
intergénico en seco **pierde terminadores reales**: de los 375 terminadores
curados de K-12, **48 (13%) están dentro de un CDS** y 74% a menos de 50 pb de
uno.

La solución: partir del intergénico y **entrar 50 pb en los genes vecinos**.

| configuración | genoma cubierto | terminadores curados cubiertos |
|---|---:|---:|
| intergénico estricto | 12,5% | ~85% |
| **intergénico + 50 pb** | **20,6%** | **92%** |
| intergénico + 100 pb | 28,0% | 94% |

## Qué anota

1. **TTS (Rho-independiente)**: escanea cada región buscando horquilla + tracto
   poli-U, y clasifica la estructura (A fuerte / B media / C débil).
2. **TSS**: no se predicen. Se cargan de un TSV (`pos, strand`) y se asignan al
   gen que empiezan. Si no hay datos, no se inventan.
3. **TTS externos**: `--batter` incorpora un BED de otro detector (BATTER-TPE),
   que tiene mayor recall y complementa a las horquillas.

## Casos anómalos: la parte importante

Cada detección se clasifica por dónde cae:

| clase | significado | qué hacer |
|---|---|---|
| `intergenic` | entre genes: lo esperado | usar |
| `gene_edge` | dentro de los 50 pb del margen | **legítimo**: los terminadores se solapan con el final del gen |
| `inside_gene` | más adentro que el margen | raro: antisentido mal etiquetado, atenuador interno o desajuste de coordenadas |
| `antisense` | dentro de un gen de la **hebra contraria** | **no es un operón**: termina un transcrito antisentido |
| `contig_edge` | a menos de 100 pb del borde | revisar; suele ser artefacto |

Esta clasificación no es decorativa. De los 375 terminadores curados de K-12,
**11 son antisentido** y **19 caen dentro de un gen** (probablemente por
desajuste de coordenadas entre RegulonDB y el GenBank). Mezclarlos con el resto
infla la cobertura y contamina cualquier benchmark.

## Uso

```bash
python scripts/regulation/pipeline.py \
    --gbff data/test_ecoli.gbff \
    --tss data/regulondb/tss_k12.tsv \
    --ref-tts data/regulondb/TerminatorSet.tsv \
    --out data/operon_finder/k12_tss_tts
```

Parámetros:

| opción | efecto |
|---|---|
| `--margin` | pb dentro de los genes vecinos (default 50) |
| `--min-region` | descarta regiones intergénicas menores (default 20) |
| `--max-region` | trunca regiones enormes al escanear (default 400) |
| `--batter` | BED de terminadores de BATTER-TPE |
| `--ref-tts` | conjunto de referencia para medir cobertura |
| `--tol` | tolerancia de comparación (default 25 pb) |

## Salidas

- `<out>.tsv` — una fila por detección: tipo, contig, posición, hebra, clase,
  estructura, ubicación y gen asignado.
- `<out>.json` — resumen por tipo y por clase de anomalía.

## Resultados

| especie | detecciones | TSS | TTS | cobertura de referencia |
|---|---:|---:|---:|---:|
| *E. coli* K-12 | 6.849 | 3.964 | 2.885 | 221/375 (58,9%) |
| *B. subtilis* 168 | 6.399 | 3.390 | 3.009 | 1.050/2.566 (40,9%) |

La diferencia entre especies no es un fallo del escaneo: BSGatlas incluye
terminadores dependientes de Rho y predichos por otros métodos, que una
horquilla no puede detectar por definición. Para igualar la comparación hay que
estratificar por mecanismo.

## Límites

- **La horquilla solo ve Rho-independientes.** Los dependientes de Rho no tienen
  estructura local; necesitan otro detector.
- **El TSS no se predice.** Sin datos experimentales, no hay TSS. Predecir con
  PWM da ~19% de sensibilidad con 3% de falsos positivos: insuficiente como
  cimiento.
- **La cobertura de referencia es un suelo**, no una precisión: los conjuntos
  curados no son exhaustivos.
