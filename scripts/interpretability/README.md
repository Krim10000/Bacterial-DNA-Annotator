# Interpretabilidad de embeddings v2

Extrae los **384 valores** del embedding por ventana de ADN y los cruza con
anotaciones genómicas. No reduce todo a un AUROC: produce la tabla de
dimensiones, el efecto de cada dimensión y su correlación con cada anotación.

## Flujo

```
FASTA + GFF
   │
   ├─ 01_extract_embeddings.py
   │     muestrea ventanas (CDS, tRNA, rRNA, intergénico, upstream, oriC,
   │     CRISPR/profago/transposón, aleatorias)
   │     → embeddings.npy [N, 384] + windows.tsv + embeddings.tsv
   │
   └─ 02_correlate_annotations.py
         para cada dim × anotación:
           r, rho, r parcial | GC, d de Cohen, AUC 1D, q FDR, p permutación
         → correlations.tsv, matrices, heatmap, dimension_summary.tsv
```

## Uso

```bash
cd dna_poc
source .venv/bin/activate

# 1) Piloto: 6 genomas, ~50 ventanas por clase, ventana de 1.024 pb
python scripts/interpretability/01_extract_embeddings.py \
    --n-genomes 6 --max-per-class 50 --window 1024 \
    --out data/embeddings/v2_epoch9_pilot6

# 2) Correlación global + within-genome
python scripts/interpretability/02_correlate_annotations.py \
    --run data/embeddings/v2_epoch9_pilot6 --scope both --permutations 100
```

Genomas concretos:

```bash
python scripts/interpretability/01_extract_embeddings.py \
    --genomes GCF_000006985.1_ASM698v1,GCF_000006685.1_ASM668v1 \
    --max-per-class 80 --window 1024 \
    --out data/embeddings/v2_epoch9_pilot2
```

## Salidas

| archivo | contenido |
|---|---|
| `windows.tsv` | una fila por ventana: genoma, contig, coordenadas, clase de muestreo, `feature_type`, flags `is_*`, GC, GC skew, entropía, solapamiento con CDS/tRNA/rRNA, distancia a feature |
| `embeddings.npy` | matriz `float32 [N, 384]`, mismo orden que `windows.tsv` |
| `embeddings.tsv` | `windows.tsv` + `dim_000` … `dim_383` (la tabla completa legible) |
| `manifest.json` | sha256 del checkpoint y tokenizer, genomas, parámetros, conteos |
| `correlations.tsv` | tabla larga `scope × dim × anotación` con r, r parcial, d, AUC 1D, q, p de permutación |
| `correlations_partial_r_<scope>.tsv` | matriz 384 × anotaciones (r parcial controlando GC) |
| `correlations_raw_r_<scope>.tsv` | matriz 384 × anotaciones (r bruta) |
| `dimension_summary.tsv` | para cada una de las 384 dims: la anotación que mejor la explica (incluye composición), r, q, efecto |
| `dimension_summary_bio.tsv` | para cada dim, la mejor **clase discreta** (CDS/tRNA/rRNA/…) |
| `dimension_summary_overlap.tsv` | para cada dim, el mejor solapamiento continuo (`overlap_cds`, etc.) |
| `annotation_summary.tsv` | por anotación: cuántas dims sobreviven FDR, máxima correlación, dim top |
| `top_correlations.tsv` | dims con `q_parcial < 0.05` y `|r_parcial| ≥ 0.2` (todas las familias) |
| `top_correlations_bio.tsv` | solo clases discretas (CDS/tRNA/rRNA/…) |
| `heatmap_partial_r_<scope>.png` | mapa de calor 384 × anotaciones |
| `report_correlations.md` | informe de correlaciones (paso 02) |
| `report_signatures.md` | informe de firmas (paso 03) |
| `report_stability.md` | informe de estabilidad (paso 04) |

## Cómo leer los dos alcances

- **`global`**: todas las ventanas juntas. Mezcla diferencias entre genomas.
  Un genoma con GC o estilo de vida distinto hace que muchas dimensiones lo
  capturen; cualquier anotación sesgada por genoma parecerá correlacionada.
  El piloto de 2 genomas mostró `split_excluded` con `r ≈ 0.96`: la identidad
  del genoma domina.
- **`within_genome`**: centra cada variable y cada dimensión por genoma
  (efectos fijos). Responde “dentro de un mismo genoma, ¿qué dimensión
  distingue esta anotación?”. Es el alcance que hay que mirar para biología.

## Columnas clave

| columna | significado |
|---|---|
| `pearson_r` | correlación bruta entre la dimensión y la anotación |
| `partial_r_gc` | correlación tras eliminar GC y GC²; si cae a ~0, la dim codifica composición |
| `cohen_d` | diferencia de medias estandarizada (solo binarias) |
| `auc_single_dim` | AUC de esa única dimensión en la dirección de la correlación (descriptivo) |
| `auc_strength` | `max(AUC, 1-AUC)`: separabilidad de esa dimensión, sin importar el signo |
| `q_partial` | FDR Benjamini-Hochberg dentro de cada anotación (384 tests) |
| `perm_p_max` | p de permutación del máximo `|r|` entre las 384 dims (conservador) |
| `delta_mean` | `mean_pos - mean_neg` de la dimensión (solo binarias) |

## Límites que hay que declarar en cualquier informe

1. **El checkpoint actual (epoch 9) está subentrenado.** Perplejidad ~474 en el
   genoma de entrenamiento, peor que un unigram; no generaliza promotores
   (AUROC ~0.57 en un genoma no visto). Este análisis describe lo que el
   embedding contiene, no que sea biológicamente correcto.
2. **Fuga de entrenamiento.** El split de entrenamiento es por contig, no por
   genoma. Cualquier genoma de `data/raw` estuvo en entrenamiento. Para una
   conclusión de generalización hace falta un hold-out por genoma.
3. **Ventanas solapadas.** Los p-valores no corrigen dependencia espacial.
   Para inferencia formal, usar ventanas no solapadas y split por genoma.
4. **Anotaciones en curso.** `windows.tsv` guarda `locus_tag`, `gene`,
   `product` y coordenadas para poder re-unir anotaciones nuevas sin
   re-extraer embeddings. El embedding es independiente de la anotación.
5. **Clases pequeñas.** `is_oric`, `is_prophage`, etc. se omiten si no llegan a
   `--min-class` (10 por defecto). En el piloto de 2 genomas, `is_oric` y
   `is_prophage` quedaron fuera.

## Escalar a los 400 genomas

La extracción es el cuello de botella: ~0.4 s/ventana en CPU con carga alta.
Para 400 genomas × ~300 ventanas ≈ 120.000 ventanas, usar GPU (Modal) o
reducir `--max-per-class`. El script ya soporta `--device cuda` y lotes.
Sugerencia: 8-12 genomas diversos para el análisis interpretativo; los 400
sirven para validar la estabilidad de las dimensiones top, no para descubrir
biología nueva.

## Reanudar cuando las anotaciones estén listas

```bash
# No re-extraer: los embeddings ya existen.
python scripts/interpretability/02_correlate_annotations.py \
    --run data/embeddings/v2_epoch9_pilot6 --scope both
```

Si las anotaciones nuevas añaden columnas a `windows.tsv`, el script las
detecta automáticamente si están en `BINARY` o `CONTINUOUS`; para columnas
nuevas, añadirlas a esas listas.

También se puede pasar una tabla de metadatos enriquecida sin re-extraer:

```bash
python scripts/interpretability/02_correlate_annotations.py     --run data/embeddings/v2_epoch9_pilot6     --windows data/embeddings/v2_epoch9_pilot6/windows_enriched.tsv     --scope both
```

La única condición es que el orden de filas coincida con `embeddings.npy`
(la extracción ya lo garantiza; para una tabla externa, conservar `window_id`).


## Resultado del piloto (6 genomas, 2.298 ventanas)

- Genomas: 5 de vida libre + 1 simbionte (`raw_excluded`), 6 filos.
- Ventana: 1.024 pb; 384 dims; 100 permutaciones.
- Composición: 403 CDS, 282 tRNA, 55 rRNA, 444 intergénico, 402 upstream,
  160 transposón, 77 profago, 44 CRISPR.

Correlaciones **within-genome** (la identidad del genoma no explica la señal):

| anotación | dim top | r parcial | q | d | AUC 1D |
|---|---|---:|---:|---:|---:|
| `is_trna` | dim_012 | -0.296 | 4.6e-45 | -0.79 | 0.72 |
| `is_rrna` | dim_038 | +0.278 | 3.3e-39 | +1.22 | 0.77 |
| `is_cds` | dim_018 | -0.250 | 1.6e-31 | +0.61 | 0.68 |
| `is_upstream` | dim_375 | +0.059 | 0.23 | — | — |
| `is_transposon` | dim_142 | -0.096 | 6.2e-4 | — | — |

Lectura:

- Hay señal débil pero reproducible para tRNA, rRNA y CDS.
- **No hay señal de promotor/upstream** (q=0.23), coherente con el AUROC ~0.57
  del checkpoint en genomas no vistos.
- La mayor parte de la varianza es composición: `entropy` |r|≈0.75,
  `gc_skew` |r|≈0.58. Tras controlar GC/GC², tRNA/rRNA/CDS sobreviven.
- En global, `split_excluded` llega a |r|≈0.51: la identidad del genoma domina.
  Por eso el alcance principal es `within_genome`.


## Firmas multivariadas por componente

`02_correlate_annotations.py` responde "qué dimensión correlaciona con qué
anotación". `03_component_signatures.py` responde otra pregunta: **¿qué patrón
de 384 valores caracteriza a cada componente?**

```bash
# Firma local (256 pb): recomendada para el patrón del componente
python scripts/interpretability/01_extract_embeddings.py     --n-genomes 6 --max-per-class 60 --window 256     --out data/embeddings/v2_epoch9_pilot6_256

python scripts/interpretability/03_component_signatures.py     --run data/embeddings/v2_epoch9_pilot6_256
```

Salidas:

| archivo | contenido |
|---|---|
| `class_centroids.tsv` | media cruda de las 384 dims por componente (el patrón en el espacio del modelo) |
| `class_cohen_d.tsv` | patrón estandarizado: `d` de cada dim, componente vs resto |
| `class_signature_selected.tsv` | firma dispersa: dims con `q_ventana<0.01`, `|d|>=0.25` y consistencia ≥0.8 |
| `class_signature_dims.tsv` | tabla larga con media, diferencia, `d`, `q_ventana`, `q_meta`, consistencia |
| `signature_scores.tsv` | coseno de cada ventana contra la firma de cada componente + predicción |
| `classification.json` | LOGO por genoma: accuracy, balanced accuracy, confusión, AUC one-vs-rest |
| `heatmap_signatures.png` | mapa de calor componentes × dims de la firma |

### Resultado (256 pb, 6 genomas, ajustado por GC y genoma)

Criterio: `q_ventana < 0.01`, `|d| >= 0.25`, consistencia >= 0.8. `q_meta` se
reporta aparte porque con 6 genomas tiene poca potencia.

| componente | n | dims firma | dims con q_meta<0.05 | patrón top (d) |
|---|---:|---:|---:|---|
| rRNA | 58 | 274 | 0 | dim_201 −1.72 · dim_298 +1.66 · dim_142 +1.64 |
| tRNA | 313 | 266 | 5 | dim_018 +1.21 · dim_380 +1.17 · dim_127 +1.17 |
| CDS | 1404 | 121 | 0 | dim_142 −0.61 · dim_324 +0.59 · dim_012 +0.59 |
| transposon | 180 | 8 | 0 | dim_061 +0.35 · dim_147 −0.33 |
| intergenic | 733 | 0 | 0 | sin firma tras ajustar |
| intergenic_clean | 307 | 0 | 0 | sin firma tras ajustar |
| upstream | 472 | 0 | 0 | sin firma |
| profago | 77 | 0 | 0 | sin firma |
| CRISPR | 44 | 0 | 0 | sin firma (solo 3 genomas) |

Lectura:

- **tRNA y rRNA tienen las firmas más fuertes y consistentes** (6/6 genomas,
  `|d|` > 1 en varias dims). CDS tiene una firma más débil (`|d|≈0.6`) pero
  clara (121 dims).
- **Intergénico pierde la firma al controlar GC y genoma.** Antes del ajuste
  aparecían 125 dims; eran composición/identidad de genoma, no una firma de
  "no gen". Es el resultado metodológico más importante.
- **Upstream sigue sin firma**, incluso con ventana local de 256 pb. No hay un
  patrón de promotor aprendido.
- **Transposón conserva una firma débil** (8 dims).
- **La clasificación LOGO mejora al ajustar**: balanced accuracy 0.414 vs 0.25
  de azar; AUC one-vs-rest tRNA 0.761, rRNA 0.707, CDS 0.623, intergénico 0.486.
  Es decir, el patrón existe y generaliza parcialmente entre genomas, pero el
  modelo no separa los cuatro componentes de forma limpia.

Con ventana de 1 kb y el mismo ajuste, el patrón es estable: CDS 131 dims,
tRNA 146, rRNA 252, intergénico 0; balanced accuracy 0.408. Eso sugiere que las
firmas de tRNA/rRNA/CDS no dependen del tamaño de ventana, mientras que la de
intergénico sí era un artefacto de contexto/composición.

### ¿Cómo se usa el patrón para un elemento nuevo?

1. Extraer el embedding de la ventana (384 dims).
2. Restarle la media global de la tabla.
3. Calcular el coseno contra cada columna de `class_cohen_d.tsv`.
4. El componente con mayor coseno es el candidato; `signature_scores.tsv`
   muestra ese cálculo para todas las ventanas del piloto.

La firma dispersa (`class_signature_selected.tsv`) es la versión interpretable:
solo las dims que pasaron significancia y consistencia, con su dirección.


## Escalado a 20 genomas y estabilidad

Run: `data/embeddings/v2_epoch9_20g_256` (20 filos, 9.849 ventanas de 256 pb,
16 de vida libre + 4 simbiontes).

```bash
python scripts/interpretability/03_component_signatures.py     --run data/embeddings/v2_epoch9_20g_256 --adjust gc_genome

python scripts/interpretability/04_signature_stability.py     --run data/embeddings/v2_epoch9_20g_256     --pilot data/embeddings/v2_epoch9_pilot6_256 --splits 20
```

| componente | dims firma | `|d|` máx | `q_meta` mín | Jaccard split-half | corr(d) | dims estables |
|---|---:|---:|---:|---:|---:|---:|
| tRNA | 249 | 1,00 | 2,5e−10 | 0,775 | 0,943 | 263 |
| rRNA | 63 | 1,30 | 1,8e−3 | 0,146 | 0,765 | 127 |
| CDS | 32 | 0,42 | 7,0e−4 | 0,149 | 0,366 | 52 |
| intergénico | 0 | — | — | 0,001 | −0,478 | 0 |
| upstream | 0 | — | — | — | — | — |

- **tRNA replica casi intacta del piloto de 6 genomas** (Jaccard 6g vs 20g =
  0,88): la firma es una propiedad estable, no un artefacto de 6 linajes.
- rRNA tiene patrón global estable (`corr(d)=0,76`) pero selección dispersa
  porque solo hay 143 ventanas; con más rRNA se estabilizaría.
- CDS es débil e inestable (`corr(d)=0,37`): el modelo distingue CDS, pero no
  de forma reproducible entre linajes.
- Intergénico/upstream no tienen firma; no conviene añadir una cabeza de
  promotor sobre este embedding.

Clasificación LOGO (20 genomas): balanced accuracy 0,407 vs 0,25 de azar;
AUC tRNA 0,748, rRNA 0,785, CDS 0,605, intergénico 0,474.

El dictamen de diseño del próximo modelo está en
[`reports/NEW_MODEL_DESIGN.md`](../docs/NEW_MODEL_DESIGN.md).


## Tokenizer v3 (para el próximo modelo)

Se entrenó y evaluó un barrido de tokenizers con las correcciones de diseño:
IUPAC completo, `<pad>/<unk>/<bos>/<eos>/<mask>` separados, sin `▁`, motivos
protegidos.

| tokenizer | vocab/max | nt/token | tokens 1,386 Gbp | IUPAC | motivos | `▁` |
|---|---:|---:|---:|---:|---:|---:|
| v2 | 4.096/6 | 4,00 | 346M | ❌ | 7/7 | 2 |
| v3a | 1.024/5 | 3,43 | 404M | ✅ | 7/7 | 0 |
| v3b | 512/5 | 3,28 | 423M | ✅ | 7/7 | 0 |
| **v3c** | **256/4** | **2,97** | **467M** | ✅ | 7/7 | 0 |

- **Recomendado: v3c 256/4.** Cumple 2–3 nt/token y todas las correcciones.
- **8k/3 es imposible:** con 4 letras solo hay 84 piezas de ≤3 nt.
- **Límite de solapamiento:** `TGA` en contexto 0,658 (choca con
  `ATG/GTG/TTG`); `AGGAGG` 0,703 (contiene `GGAGG`); `TATAAT`/`TTGACA` 1,00.
  El próximo modelo necesita una prioridad de solapamiento o un pretokenizador
  por límites de codón.

Ver [TOKENIZER_V3_REPORT.md](../docs/TOKENIZER_V3_REPORT.md) y
[NEW_MODEL_DESIGN.md](../docs/NEW_MODEL_DESIGN.md).

```bash
python scripts/06_train_tokenizer_v3.py --vocab-size 256 --max-piece 4     --out data/tokenizer_v3_256
python scripts/06_eval_tokenizer_v3.py     --v3 data/tokenizer_v3_256/dna_sp_v3.model     --out data/tokenizer_v3_256
```


## Benchmark de modelos pretrained

Se evaluaron ProkBERT-mini-long (27M, 4 kb) y NT-v2 50M (12 kb) con las mismas
ventanas/etiquetas que el modelo propio. Resultado (600 ventanas, 20 genomas,
LOGO, ajuste GC+genoma):

| modelo | balanced acc | CDS | tRNA | rRNA | intergénico | promotor |
|---|---:|---:|---:|---:|---:|---:|
| ProkBERT | **0,825** | 0,981 | 0,992 | 0,967 | 0,909 | **0,988** |
| NT-v2 50M | 0,510 | 0,818 | 0,840 | 0,663 | 0,626 | 0,908 |
| v2 epoch9 | 0,376 | 0,669 | 0,708 | 0,701 | 0,528 | 0,767 |
| 3-mer | 0,502 | 0,719 | 0,732 | 0,842 | 0,638 | 0,850* |

\* conteo de TATAAT. Detalle, capacidad del equipo y caveats en
[`docs/PRETRAINED_BENCHMARK.md`](../../docs/PRETRAINED_BENCHMARK.md).

```bash
export HF_HOME=$PWD/.cache/huggingface

python scripts/interpretability/05_pretrained_benchmark.py     --model prokbert     --windows data/embeddings/_pretrained_subset600/windows.tsv     --out data/embeddings/pretrained/prokbert_subset600

python scripts/interpretability/03_component_signatures.py     --run data/embeddings/pretrained/prokbert_subset600 --permutations 0

python scripts/interpretability/05_pretrained_benchmark.py     --model prokbert     --windows data/embeddings/_promoter_windows/windows.tsv     --out data/embeddings/pretrained/prokbert_promoter
python scripts/interpretability/06_promoter_probe.py     --run data/embeddings/pretrained/prokbert_promoter
```
