# Benchmark de validación del módulo sigma (TSS-level)

> **Estado:** benchmark reproducible, CPU, sin red. Fecha de la corrida: 2026-11-25.
>
> **Pregunta:** ¿el scanner de PWMs del módulo sigma recupera el factor sigma
> curado de un TSS? ¿En qué tipos falla? ¿La confianza del módulo ordena las
> predicciones? ¿Qué cambia al pasar de TSS a operón?
>
> **Gold curado y evaluación por capas:** [`SIGMA_LAYERED_EVAL.md`](SIGMA_LAYERED_EVAL.md).
> Con el gold nuevo, K-12 top-1 = 57,0 % (mayoría 62,6 %) y Subtilis 44,8 %
> (mayoría 87,6 %); este documento conserva el benchmark anterior como
> referencia del adaptador.
>
> **Respuesta corta:** a nivel TSS, el top-1 acierta **52.22%** en K-12 y
> **49.55%** en *B. subtilis*; el top-3 sube a **87.69%** y **87.60%**, pero
> está dominado por σ70. En K-12 el top-1 **no supera la mayoría trivial**
> (60.24%) y en *B. subtilis* tampoco (87.56%). La confianza del módulo sí
> ordena las predicciones (K-12: 65.83% en confianza media, 51.25% en baja,
> 36.15% sin sospecha), pero no es una probabilidad calibrada. La comparación
> por operón da **40.91%** de coincidencia exacta y expone una limitación
> estructural: una fila por operón no puede representar varios TSS/sigma.

Este documento describe el método, pega los números reales de la corrida y
separa explícitamente lo que el benchmark **sí** demuestra de lo que **no**.

---

## 1. Qué se construyó

| Archivo | Rol |
|---|---|
| [`scripts/sigma/30_validate_sigma.py`](../scripts/sigma/30_validate_sigma.py) | CLI del benchmark TSS-level, métricas, confusión y calibración |
| [`tests/test_sigma_validation.py`](../tests/test_sigma_validation.py) | tests de mapeo, ventana, métricas, calibración y smoke tests |
| [`data/sigma/benchmark/k12_tss_predictions.tsv`](../data/sigma/benchmark/k12_tss_predictions.tsv) | una fila por TSS de K-12, incluidos `no_pwm`, sin etiqueta y sin coordenada |
| [`data/sigma/benchmark/k12_confusion.json`](../data/sigma/benchmark/k12_confusion.json) | matriz de confusión TSS-level de K-12 |
| [`data/sigma/benchmark/k12_metrics.json`](../data/sigma/benchmark/k12_metrics.json) | cobertura, accuracy, P/R/F1, calibración y comparación por operón |
| [`data/sigma/benchmark/subtilis_tss_predictions.tsv`](../data/sigma/benchmark/subtilis_tss_predictions.tsv) | una fila por TSS de *B. subtilis* |
| [`data/sigma/benchmark/subtilis_confusion.json`](../data/sigma/benchmark/subtilis_confusion.json) | matriz de confusión de *B. subtilis* |
| [`data/sigma/benchmark/subtilis_metrics.json`](../data/sigma/benchmark/subtilis_metrics.json) | métricas y mapeo de etiquetas de *B. subtilis* |
| [`data/sigma/benchmark/benchmark_summary.md`](../data/sigma/benchmark/benchmark_summary.md) | resumen generado automáticamente por el CLI |

No se tocó la GUI, `app.py`, `browser.py`, `run_pipeline.sh` ni las salidas
existentes `data/sigma/<genoma>/`. El benchmark solo escribe bajo
`data/sigma/benchmark/`.

---

## 2. Datos

### 2.1 K-12

| Dato | Fuente | Filas | Uso |
|---|---|---:|---|
| TSS curados | [`data/regulondb/gt/gt_promoters.tsv`](../data/regulondb/gt/gt_promoters.tsv) | 4.057 | verdad de referencia |
| Secuencia | [`data/allostery/work/ecoli_k12.fna`](../data/allostery/work/ecoli_k12.fna) / [`data/test_ecoli.gbff`](../data/test_ecoli.gbff) | 4.641.652 pb | ventanas y fondo nulo |
| PWMs | [`data/regulondb/gt/promoter_pwms.tsv`](../data/regulondb/gt/promoter_pwms.tsv) | 6 sigma | scanner |
| Salida existente | [`data/sigma/ecoli_k12/sigma_promoters.tsv`](../data/sigma/ecoli_k12/sigma_promoters.tsv) | 2.544 operones | comparación operón/gen |

De las 4.057 filas, **3.964** tienen coordenada de TSS válida y **93** no
(68 de ellas sin etiqueta y 25 con etiqueta). Las 716 filas sin etiqueta se
reparten en 648 con coordenada válida y 68 sin coordenada.

### 2.2 *Bacillus subtilis*

| Dato | Fuente | Filas | Uso |
|---|---|---:|---|
| TSS curados | [`data/dbtbs/gt_tss.tsv`](../data/dbtbs/gt_tss.tsv) | 3.390 | verdad principal |
| TSS DBTBS | [`data/dbtbs/dbtbs_tss.tsv`](../data/dbtbs/dbtbs_tss.tsv) | 865 | solo enriquece la columna `operon` |
| Secuencia | [`data/dbtbs/genome.fna`](../data/dbtbs/genome.fna) | 4.215.606 pb | ventanas y fondo nulo |
| PWMs | [`data/regulondb/gt/promoter_pwms.tsv`](../data/regulondb/gt/promoter_pwms.tsv) | 6 sigma de *E. coli* | scanner transferido |

`dbtbs_tss.tsv` tiene 318 posiciones+hebra compartidas con `gt_tss.tsv` y 547
que solo aparecen en DBTBS. Esas 547 **no se agregaron** como verdad extra:
en los solapamientos hay discrepancias de etiqueta entre fuentes (p. ej.
`F`/`SigF`/`SigG`) y mezclar fuentes habría cambiado el denominador sin una
reconciliación explícita. Se usan solo para completar el nombre del operón.

---

## 3. Método

### 3.1 Ventana alrededor del TSS

Para un TSS en coordenada 1-based `p` y hebra `s`:

- hebra `+`: `seq[p-1-80 : p-1+20+1]` (101 nt);
- hebra `-`: `revcomp(seq[p-1-20 : p-1+80+1])`.

La ventana resultante tiene **101 nt** (`-80..+20`, TSS incluido) y el TSS queda
en el índice 0-based **80** en ambas hebras. La ventana se orienta en el sentido
de transcripción, de modo que las cajas -35/-10 (o -24/-12 para σ54) se buscan
en la misma orientación que el promotor.

Elección y límites:

- `-80/+20` cubre las arquitecturas de las seis PWMs (σ54 usa -24/-12 y
  espaciador 10–11; las demás -35/-10 con espaciadores 12–19).
- No se modelan elementos reguladores lejanos, UP-elements ni cambios de
  conformación del DNA.
- No se exige que la caja -10 esté a una distancia fija del TSS: `score_sigma`
  elige la mejor geometría dentro de la ventana y la limitación se discute en
  §6.

### 3.2 Scoring y fondo nulo

Para cada TSS evaluable:

1. se puntúan los seis sigma con `score_sigma(seq, rec, window_len)` sobre la
   misma ventana;
2. se calcula `z = _zscore(score, null[sigma])` y `p = _pvalue(score, null[sigma])`;
3. se ordenan los sigma por `(z, -p, score)`;
4. `top-1` es el primer sigma; `top-3` indica si la verdad está entre los tres
   primeros; `pred_rank` es la posición de la verdad (1–6);
5. `pred_margin = z(top-1) - z(segundo)`.

El fondo nulo se construye con `build_null` sobre **500 ventanas genómicas
aleatorias de 101 nt**, semilla **17**, en el mismo genoma y con las mismas
PWMs. El nulo es genómico, no un modelo de promotores: incluye promotores
reales y regiones no promotoras. Por eso el z-score es una referencia interna,
no un p-valor biológico.

### 3.3 Mapeo de etiquetas

**K-12** ya usa los nombres `sigma70`, `sigma54`, `sigma32`, `sigma38`,
`sigma28` y `sigma24`. `sigma19` no tiene PWM en
`promoter_pwms.tsv`; se reporta como `no_pwm` y **no se cuenta como fallo**.
Las filas sin etiqueta se reportan como `unlabeled`.

**B. subtilis** usa letras. Se mapean solo las etiquetas con una PWM
equivalente defendible:

| Etiqueta BSU | Sigma usado | Justificación | n |
|---|---|---|---:|
| `A` | `sigma70` | SigA, sigma housekeeping, ortólogo funcional de RpoD | 1.949 |
| `B` | `sigma38` | SigB, sigma de estrés general, análogo funcional de RpoS | 182 |
| `D` | `sigma28` | SigD, sigma flagelar, ortólogo funcional de FliA | 71 |
| `I` | `sigma32` | SigI, sigma de choque térmico, análogo funcional de RpoH | 13 |
| `L` | `sigma54` | SigL, familia σ54, relacionada con metabolismo de N | 11 |

No se mapean:

- esporulación `E`, `F`, `G`, `K` (no hay PWM equivalente en el repositorio);
- competencia/transición `H` (no hay equivalente fiable);
- ECF `W`, `X`, `M` (la PWM ECF de *E. coli* no es transferible);
- `YlaC` (etiqueta no reconocida) y `?`/vacío (`unlabeled`).

Todas las etiquetas sin PWM se reportan como `no_pwm` y quedan **fuera del
denominador** de accuracy; no son falsos negativos del scanner.

### 3.4 Confianza

Se usa `_confidence` del módulo con la evidencia disponible:

- `p`, `z` y `margin` del top-1;
- `inv=True` solo si el genoma tiene un factor del tipo predicho. En K-12 se
  usa [`data/sigma/ecoli_k12/sigma_factors.tsv`](../data/sigma/ecoli_k12/sigma_factors.tsv)
  (con el alias `sigmaE → sigma24`);
- `pathway=False` y `expr="not_available"`: el benchmark TSS-level no tiene
  contexto funcional ni expresión.

Con esa configuración, `sospecha_fuerte` es inalcanzable (requiere contexto o
expresión) y `sospecha_media` requiere inventario + umbrales estadísticos.
*B. subtilis* no tiene inventario sigma en el repositorio, así que solo puede
producir `sospecha_baja`/`sin_sospecha`. Esto es una limitación del benchmark,
no del genoma.

La calibración se reporta de tres formas:

- por etiqueta de confianza (`sospecha_media`, `sospecha_baja`, `sin_sospecha`);
- por bin de `confidence_score` (0.0–0.2, …, 0.8–1.0);
- por bin de `pred_z` (<1, 1–2, 2–3, 3–4, ≥4).

`confidence_score` **no es una probabilidad**. ECE y Brier se incluyen como
diagnósticos de ordenamiento, no como probabilidades calibradas.

### 3.5 Comparación K-12 por operón/gen

La salida existente tiene una fila por operón con `first_gene_name`. El GT tiene
TSS con `first_gene` (o un nombre de promotor como `yfeKp`, del que se deriva
`yfeK`). Para cada gen del GT:

- `exact`: el sigma predicho está entre los sigma curados para ese gen;
- `dominant`: el sigma predicho coincide con el sigma curado más frecuente;
- `multi-sigma`: el gen tiene ≥2 sigma curados distintos.

`exact` es más permisivo en genes con varios TSS; `dominant` es más estricto.
Ninguno de los dos puede representar varios sigma a la vez.

---

## 4. Resultados K-12

### 4.1 Cobertura y accuracy

| Métrica | Valor |
|---|---:|
| TSS en la tabla | 4.057 |
| TSS con coordenada válida | 3.964 |
| TSS sin coordenada | 93 |
| TSS evaluables (σ con PWM) | **3.315** |
| TSS `no_pwm` (σ19) | 1 |
| TSS sin etiqueta (con coordenada) | 648 |
| Cobertura sobre etiquetados con PWM | **99.97%** |
| Cobertura sobre TSS válidos | 83.63% |
| Cobertura sobre todas las filas | 81.71% |
| Tasa `no_pwm` sobre etiquetados | 0.03% |
| **Accuracy top-1** | **52.22%** (1.731/3.315) |
| **Accuracy top-3** | **87.69%** (2.907/3.315) |
| Baseline de mayoría (σ70) | 60.24% (1.997/3.315) |
| Ganancia sobre mayoría | **−8.02 pp** |
| Rank medio de la verdad | 1.92 |
| Macro-F1 (clases con apoyo) | 0.433 |
| Macro-F1 (6 clases, F1 ausente = 0) | 0.433 |
| Weighted-F1 | 0.558 |

Lectura: el top-3 es alto porque σ70 está casi siempre entre los tres primeros;
el top-1 no alcanza a la mayoría trivial. El módulo **no está recuperando bien
la clase mayoritaria** y, al mismo tiempo, no compensa con las clases
alternativas.

### 4.2 Por tipo de sigma

| sigma | support | predichos | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sigma70 | 1.997 | 1.119 | 959 | 160 | 1.038 | 0.857 | 0.480 | 0.616 |
| sigma54 | 99 | 204 | 35 | 169 | 64 | 0.172 | 0.354 | 0.231 |
| sigma32 | 312 | 367 | 192 | 175 | 120 | 0.523 | 0.615 | 0.566 |
| sigma38 | 241 | 696 | 126 | 570 | 115 | 0.181 | 0.523 | 0.269 |
| sigma28 | 144 | 261 | 66 | 195 | 78 | 0.253 | 0.458 | 0.326 |
| sigma24 | 522 | 668 | 353 | 315 | 169 | 0.528 | 0.676 | 0.593 |

σ70 tiene buena precision (0.857) pero recall bajo (0.480): pierde casi la
mitad de sus TSS, que se reparten sobre todo en σ38 (447), σ24 (228) y σ28
(131). σ38 es el falso positivo más claro: 696 predicciones para 241 TSS
reales. σ32 y σ24 son los mejores alternativos (F1 0.566 y 0.593); σ54 y σ38
son los peores (F1 0.231 y 0.269).

### 4.3 Matriz de confusión

Filas = etiqueta curada; columnas = top-1.

| verdad \ pred | σ70 | σ54 | σ32 | σ38 | σ28 | σ24 |
|---|---:|---:|---:|---:|---:|---:|
| **σ70** | **959** | 110 | 122 | 447 | 131 | 228 |
| **σ54** | 19 | **35** | 12 | 13 | 7 | 13 |
| **σ32** | 38 | 15 | **192** | 34 | 13 | 20 |
| **σ38** | 35 | 14 | 12 | **126** | 20 | 34 |
| **σ28** | 19 | 5 | 12 | 22 | **66** | 20 |
| **σ24** | 49 | 25 | 17 | 54 | 24 | **353** |

### 4.4 Calibración por confianza

| Etiqueta | n | correctos | accuracy | confianza media | gap |
|---|---:|---:|---:|---:|---:|
| `sospecha_media` | 796 | 524 | **65.83%** | 0.541 | +0.118 |
| `sospecha_baja` | 1.963 | 1.006 | 51.25% | 0.390 | +0.123 |
| `sin_sospecha` | 556 | 201 | 36.15% | 0.219 | +0.143 |

El orden es correcto: media > baja > sin sospecha. El subconjunto de confianza
media es el **24.01%** de los TSS evaluables y alcanza 65.83%, por encima del
baseline de mayoría (60.24%); el resto no. Es decir, la confianza sirve para
**priorizar**, pero no para afirmar que el top-1 es correcto.

Por bins de score:

| bin | n | accuracy | confianza media | gap |
|---|---:|---:|---:|---:|
| 0.0–0.2 | 109 | 33.95% | 0.182 | +0.158 |
| 0.2–0.4 | 1.482 | 46.76% | 0.279 | +0.188 |
| 0.4–0.6 | 1.523 | 58.63% | 0.497 | +0.089 |
| 0.6–0.8 | 201 | 53.73% | 0.627 | −0.090 |

ECE = **0.136**, Brier = **0.263**. El score subestima la accuracy en la mayoría
de los bins; el último bin se invierte (201 casos). No debe interpretarse como
probabilidad.

Por bins de z:

| bin de z | n | accuracy |
|---|---:|---:|
| <1 | 255 | 33.33% |
| 1–2 | 1.416 | 49.93% |
| 2–3 | 1.386 | 58.37% |
| 3–4 | 237 | 50.21% |
| ≥4 | 21 | 52.38% |

El z ordena hasta 2–3, pero se aplana/invierte por encima de 3: un z alto no
garantiza la etiqueta correcta cuando varias PWMs compiten con scores altos.

---

## 5. Resultados *B. subtilis*

### 5.1 Cobertura y accuracy

| Métrica | Valor |
|---|---:|
| TSS en la tabla | 3.390 |
| TSS con coordenada válida | 3.390 |
| TSS evaluables (σ mapeado con PWM) | **2.226** |
| TSS `no_pwm` | 848 |
| TSS sin etiqueta (`?`) | 316 |
| Cobertura sobre etiquetados con PWM | **72.41%** |
| Cobertura sobre TSS válidos | 65.66% |
| Tasa `no_pwm` sobre etiquetados | **27.59%** |
| Tasa `no_pwm` sobre todas las filas | 25.01% |
| **Accuracy top-1** | **49.55%** (1.103/2.226) |
| **Accuracy top-3** | **87.60%** (1.950/2.226) |
| Baseline de mayoría (A→σ70) | **87.56%** (1.949/2.226) |
| Ganancia sobre mayoría | **−38.01 pp** |
| Rank medio de la verdad | 1.95 |
| Macro-F1 (clases con apoyo) | 0.204 |
| Macro-F1 (6 clases, F1 ausente = 0) | 0.170 |
| Weighted-F1 | 0.608 |

El top-1 no supera la mayoría porque el 87.56% de los TSS mapeados son `A`
(σ70) y el scanner, entrenado con PWMs de *E. coli*, asigna σ alternativos a
muchos promotores `A`. El top-3 está inflado por la misma mayoría. Este
resultado **no mide la calidad biológica de la anotación de σ en *B. subtilis***:
mide la transferencia de PWMs de *E. coli* a etiquetas homologadas.

### 5.2 Por tipo de sigma

| sigma | support | predichos | TP | FP | FN | precision | recall | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sigma70 | 1.949 | 1.140 | 1.052 | 88 | 897 | 0.923 | 0.540 | 0.681 |
| sigma54 | 11 | 105 | 3 | 102 | 8 | 0.029 | 0.273 | 0.052 |
| sigma32 | 13 | 102 | 3 | 99 | 10 | 0.029 | 0.231 | 0.052 |
| sigma38 | 182 | 535 | 26 | 509 | 156 | 0.049 | 0.143 | 0.073 |
| sigma28 | 71 | 166 | 19 | 147 | 52 | 0.114 | 0.268 | 0.160 |
| sigma24 | 0 | 178 | 0 | 178 | 0 | 0.000 | not_available | not_available |

No hay TSS curados mapeados a σ24 en *B. subtilis*; las 178 predicciones σ24
son falsos positivos. Las clases alternativas tienen apoyo pequeño (11–182),
así que sus métricas tienen mucha incertidumbre; aun así, la precision es muy
baja y el patrón es consistente con un PWM no transferible.

### 5.3 Matriz de confusión

| verdad \ pred | σ70 | σ54 | σ32 | σ38 | σ28 | σ24 |
|---|---:|---:|---:|---:|---:|---:|
| **σ70 (A)** | **1.052** | 75 | 79 | 494 | 107 | 142 |
| **σ54 (L)** | 2 | **3** | 0 | 2 | 2 | 2 |
| **σ32 (I)** | 2 | 0 | **3** | 5 | 2 | 1 |
| **σ38 (B)** | 57 | 26 | 12 | **26** | 36 | 25 |
| **σ28 (D)** | 27 | 1 | 8 | 8 | **19** | 8 |
| **σ24** | 0 | 0 | 0 | 0 | 0 | 0 |

### 5.4 Calibración por confianza

| Etiqueta | n | correctos | accuracy | confianza media | gap |
|---|---:|---:|---:|---:|---:|
| `sospecha_baja` | 1.766 | 934 | **52.89%** | 0.264 | +0.265 |
| `sin_sospecha` | 460 | 169 | 36.74% | 0.068 | +0.299 |

No hay `sospecha_media` ni `sospecha_fuerte` porque el benchmark no tiene
inventario sigma de *B. subtilis* ni contexto funcional; no es que el genoma no
tenga esos factores. El orden baja > sin sospecha se mantiene, pero el score
subestima fuertemente la accuracy (gap +0.265/+0.299).

Por bins de z:

| bin de z | n | accuracy |
|---|---:|---:|
| <1 | 214 | 31.31% |
| 1–2 | 1.023 | 44.48% |
| 2–3 | 812 | 58.25% |
| 3–4 | 174 | 62.07% |
| ≥4 | 3 | 0.00% |

ECE = **0.272**, Brier = **0.314**. El z es ordinalmente informativo (salvo el
último bin, con n=3), pero el score no es una probabilidad.

### 5.5 `no_pwm` de *B. subtilis*

| Etiqueta | n | Motivo |
|---|---:|---|
| `E` | 278 | SigE, esporulación; sin PWM equivalente |
| `K` | 209 | SigK, esporulación; sin PWM equivalente |
| `F` | 189 | SigF, esporulación; sin PWM equivalente |
| `W` | 82 | SigW, ECF; PWM ECF no transferible |
| `H` | 77 | SigH, competencia/transición; sin PWM fiable |
| `G` | 6 | SigG, esporulación; sin PWM equivalente |
| `X` | 3 | SigX, ECF; PWM ECF no transferible |
| `M` | 3 | SigM, ECF; PWM ECF no transferible |
| `YlaC` | 1 | etiqueta no reconocida |
| **Total** | **848** | 27.59% de los TSS etiquetados |

Estas 848 filas **no son fallos del scanner**: son etiquetas que el benchmark
no puede evaluar con las PWMs disponibles. Reportarlas como `no_pwm` es una
decisión metodológica central para no inflar artificialmente el error.

---

## 6. Comparación K-12: TSS vs operón/gen

La salida existente `sigma_promoters.tsv` tiene **2.544 filas**, una por
operón, con un único `sigma`. El GT tiene **2.058 genes** con al menos un TSS
mapeable a las seis PWMs.

| Métrica | Valor |
|---|---:|
| Filas existentes | 2.544 |
| Filas existentes mapeadas a un gen del GT | 1.501 |
| Genes del GT con sigma mapeable | 2.058 |
| Genes del GT mapeados a la salida existente | 1.501 |
| Genes del GT sin fila existente | 557 |
| Filas existentes no asignables a un gen del GT con sigma mapeable | 1.043 |
| **Coincidencia exacta** | **40.91%** (614/1.501) |
| Coincidencia con el sigma dominante | 32.84% (493/1.501) |
| Genes con ≥2 TSS | 912 |
| Genes con ≥2 sigma distintos | 503 |
| % genes con ≥2 sigma | **24.44%** |

### 6.1 La limitación de un sigma por operón

| Subconjunto | n genes | Coincidencia | Tasa |
|---|---:|---:|---:|
| Un solo sigma curado | 1.125 | 391 | 34.76% |
| ≥2 sigma curados: coincide con alguno | 376 | 223 | 59.31% |
| ≥2 sigma curados: coincide con el dominante | 376 | 102 | 27.13% |
| ≥2 sigma curados: coincide con **todos** | 376 | **0** | 0.00% |

La coincidencia "con alguno" es más alta en genes multi-sigma (59.31%) que en
genes de un solo sigma (34.76%) simplemente porque hay más etiquetas contra las
que acertar. La coincidencia "con todos" es **0 por construcción**: una fila
con un solo `sigma` no puede expresar dos o más factores. No es un bug del
scanner; es una limitación del nivel operón.

Además, 912 genes del GT tienen ≥2 TSS (no necesariamente ≥2 sigma): el operón
puede tener promotores alternativos o TSS internos. La salida existente colapsa
esos casos en una sola asignación.

### 6.2 Dónde acierta la salida existente

| sigma predicho | n | exactos | tasa exacta |
|---|---:|---:|---:|
| sigma70 | 405 | 360 | 88.89% |
| sigma54 | 110 | 14 | 12.73% |
| sigma32 | 146 | 62 | 42.47% |
| sigma38 | 440 | 78 | 17.73% |
| sigma28 | 206 | 29 | 14.08% |
| sigma24 | 194 | 71 | 36.60% |

La salida existente es razonable cuando predice σ70 (88.89% exacto) y mucho
peor cuando predice alternativos, sobre todo σ54, σ38 y σ28. Esto es coherente
con el sesgo del scanner: σ70 tiene buena precision y bajo recall; σ38 se
sobrepredice.

**El mapeo es aproximado.** 1.043 filas existentes no se pudieron asignar a un
gen del GT con al menos una etiqueta mapeable: el GT no siempre trae
`first_gene` (1.771 filas con TSS válido; 1.864 sobre el total) y los nombres de promotor no siempre se pueden
parsear a un gen de `data/regulation/ecoli_k12/genes.tsv`. Por eso la
comparación es una cota, no un censo exacto.

---

## 7. Límites

1. **TSS vs operón.** El benchmark puntúa TSS aislados; el módulo anota
   operones con un sigma. Son unidades distintas. Un TSS puede ser interno o
   alternativo, y un operón puede tener varios TSS.
2. **Transferencia entre especies.** *B. subtilis* se evalúa con PWMs de
   *E. coli* sobre etiquetas homologadas. `A→sigma70` no es una PWM de SigA de
   *B. subtilis*; es un análogo funcional. El mal top-1 de *B. subtilis* es en
   buena parte un fallo de transferencia, no una medida de la biología del
   organismo.
3. **Etiquetas sin PWM.** `sigma19` en K-12 y `E/F/G/K/H/W/X/M/YlaC` en
   *B. subtilis* no tienen PWM equivalente. Se reportan como `no_pwm` y no se
   cuentan como fallo. El 27.59% de los TSS etiquetados de *B. subtilis* cae
   en esta categoría.
4. **Múltiples TSS.** 912 genes de K-12 tienen ≥2 TSS y 503 tienen ≥2 sigma
   distintos. Una predicción por operón no puede representarlos; el 0% de
   coincidencia "con todos" es estructural.
5. **Desbalance de clases.** σ70 es el 60.24% en K-12 y el 87.56% en
   *B. subtilis*. El top-1 y el top-3 no son comparables entre clases; por eso
   se reportan P/R/F1 por sigma, macro-F1 y baseline de mayoría.
6. **Confianza no calibrada.** `confidence_score` es una heurística del módulo;
   ECE/Brier son diagnósticos. En *B. subtilis* no hay inventario sigma ni
   contexto funcional, así que la confianza del módulo queda limitada a
   `sospecha_baja`/`sin_sospecha`.
7. **Fondo nulo genómico.** Las 500 ventanas aleatorias contienen promotores
   reales y regiones no promotoras; el z-score es una referencia interna, no un
   p-valor biológico.
8. **Posible circularidad en K-12.** Las PWMs del repositorio se derivaron de
   RegulonDB/GT y la verdad de K-12 también es RegulonDB/GT. El benchmark K-12
   mide **consistencia interna**, no generalización a un conjunto independiente.
   La evaluación de *B. subtilis* no tiene esa circularidad, pero sí el problema
   de transferencia entre especies.
9. **Ventana fija.** `-80/+20` no captura elementos reguladores lejanos,
   UP-elements ni la geometría completa de todos los promotores. Tampoco se
   modela la distancia -10/TSS como restricción dura.
10. **Fuente de *B. subtilis*.** `gt_tss.tsv` **no es DBTBS puro**: es el GFF
    de **BSGatlas v1.1**. De sus 2.226 filas, 1.777 (79,8 %) solo citan
    Nicolas et al. (tiling array), 379 (17,0 %) mencionan DBTBS y solo 430
    cumplen DBTBS/BsubCyc + `Resolution: 0`. El benchmark de B. subtilis, por
    tanto, está dominado por BSGatlas/Nicolas, no por DBTBS curado.
    Las etiquetas de las fuentes no siempre coinciden (`F` vs `SigF`/`SigG`).
11. **Comparación por operón aproximada.** El mapeo por `first_gene_name` deja
    1.043 filas existentes sin un gen del GT con etiqueta mapeable, y 557 genes
    del GT sin fila existente.

---

## 8. Conclusión honesta

**El módulo no es validado como anotador de sigma por este benchmark.** Lo que
los números muestran es:

- El scanner **ordena** razonablemente bien: la verdad está en top-3 el 87.69%
  (K-12) y 87.60% (*B. subtilis*), y la confianza del módulo separa
  predicciones (K-12: 65.83% media vs 36.15% sin sospecha).
- Pero el **top-1 no es fiable**: 52.22% en K-12 (por debajo del 60.24% de la
  mayoría) y 49.55% en *B. subtilis* (muy por debajo del 87.56%). El top-3 está
  inflado por la clase mayoritaria y no debe presentarse como accuracy.
- El sesgo principal es **σ38 sobrepredicho** y **σ70 con recall bajo** en
  K-12; en *B. subtilis* se suma la transferencia de PWMs de *E. coli*.
- A nivel operón, la salida existente coincide exactamente en **40.91%** de los
  genes mapeados; el **24.44%** de los genes tiene ≥2 sigma curados y una fila
  por operón no puede representarlos.
- La confianza es útil para **priorizar**, no para certificar. No es una
  probabilidad y en *B. subtilis* no puede alcanzar confianza media por falta
  de inventario/contexto.

Recomendaciones:

1. No presentar el top-1 como "sigma validado"; usar top-3/alternativas y la
   etiqueta de confianza como priorización.
2. Recalibrar o reemplazar las PWMs de σ38/σ54/σ28 (y las de *B. subtilis*)
   con datos de la especie; el benchmark muestra que las PWMs actuales no
   separan bien esas clases.
3. Implementar un nivel multi-TSS/multi-sigma en la salida por operón, o
   reportar explícitamente que el operón tiene promotores alternativos.
4. Mantener `no_pwm` fuera del denominador y reportarlo como cobertura, no como
   error.
5. Repetir el benchmark con un conjunto de TSS held-out no usado para construir
   las PWMs, para medir generalización y no solo consistencia interna.

---

## 9. Reproducción

```bash
cd /home/camilo/Documentos/GitHub/DNA_TOKEN_BAC/dna_poc

# Benchmark completo (K-12 + B. subtilis)
.venv/bin/python scripts/sigma/30_validate_sigma.py --all

# Por separado
.venv/bin/python scripts/sigma/30_validate_sigma.py --k12
.venv/bin/python scripts/sigma/30_validate_sigma.py --subtilis

# Autoprueba sin datos pesados
.venv/bin/python scripts/sigma/30_validate_sigma.py --selftest

# Tests
.venv/bin/python -m pytest tests/test_sigma_validation.py -q
```

Salidas: `data/sigma/benchmark/` (predicciones TSV, confusión JSON, métricas JSON
y `benchmark_summary.md`). Tiempo de la corrida completa en este equipo:
≈36 s.
