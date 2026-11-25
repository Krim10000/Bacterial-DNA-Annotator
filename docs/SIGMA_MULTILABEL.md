# Candidatos sigma multietiqueta por TSS

> **Estado:** v1.0.0 · experimental · CPU · sin red ni descargas.
>
> **Qué hace:** construye una tabla de candidatos `(genoma, TSS, sigma)` con
> ranking top-k y evalúa el conjunto contra el gold de TSS como problema
> multietiqueta.
>
> **Qué no hace:** no reemplaza la asignación de `sigma_promoters.tsv`, no
> cambia el núcleo (`sigma_annotator.py`), no toca la GUI y no inventa TSS.
>
> **Script:** [`scripts/sigma/74_multilabel.py`](../scripts/sigma/74_multilabel.py)
> **Pruebas:** [`tests/test_sigma_multilabel.py`](../tests/test_sigma_multilabel.py)

---

## 1. Por qué existe

`data/sigma/<genoma>/sigma_promoters.tsv` tiene **una fila por operón**. Esa
representación es útil para asignar un promotor, pero no representa los TSS
que tienen más de un sigma curado. En K-12 el gold tiene:

- 3.189 TSS,
- 3.314 etiquetas TSS×sigma,
- **124 TSS multi-sigma** (2 o 3 etiquetas),
- 249 etiquetas dentro de esos 124 TSS.

Un top-1 único puede acertar como máximo una etiqueta por TSS. Por eso el top-3
cualquiera llega a 90,2 %, pero el top-3 **completo** baja a 62,9 % en los TSS
multi-sigma (sección 5). La tabla multietiqueta hace visible esa diferencia.

---

## 2. Entradas y coordenada

| Genoma | Fuente usada | Filas | TSS | Coordenada | Evaluación |
|---|---|---:|---:|---|---|
| `ecoli_k12` | `layers/ecoli_k12/layer_sigma_features.tsv` | 19.134 | 3.189 | TSS 1-based | sí |
| `subtilis` | `layers/subtilis/layer_sigma_features.tsv` | 13.356 | 2.226 | TSS 1-based | sí |

Reglas de coordenada:

- **TSS real:** sale de `layer_sigma_features.tsv` (`contig`, `tss`, `strand`).
  Es 1-based, primera base transcrita. No se infiere, no se desplaza y no se
  redondea.
- **Respaldo de operón:** si no hay features de capa, se puede usar
  `sigma_promoters.tsv`. En ese caso `tss` queda **vacío**, la coordenada va a
  `operon_anchor` y `coord_type=operon_anchor`. Esa coordenada **no es un TSS**
  y no se evalúa contra el gold de TSS.
- En esta corrida ambos genomas tienen features de capa: las salidas son
  TSS-level.

Umbrales y top-k no cambian la coordenada; solo filtran u ordenan candidatos.

---

## 3. Construcción de candidatos

Para cada TSS/operón:

1. se leen todas las filas `TSS × sigma` disponibles;
2. se ordenan por `z` descendente; desempate por `pvalue`, orden canónico
   (`sigma70`, `sigma54`, `sigma32`, `sigma38`, `sigma28`, `sigma24`) y nombre;
3. se guardan los **top-k** (3 por defecto);
4. se calcula `confidence_label` solo con señal;
5. se calcula `refined_confidence` para σ54/σ38 si el refinamiento está activo.

Las columnas `true`, `is_top1`, `n_true`, `gold_confidence` y `source` de las
features **no se usan** para construir candidatos: son metadatos del gold. El
gold se lee después, solo para evaluar. Una prueba cambia esas columnas y
verifica que el ranking no cambie.

### 3.1 Confianza del candidato

Heurística, no calibrada:

| `confidence_label` | Condición |
|---|---|
| `alta` | `pvalue <= 0.05` y `margin >= 1.0` |
| `media` | `pvalue <= 0.05`, o (`margin >= 1.0` y `z >= 1.0`) |
| `baja` | `z >= 0.0`, o `pvalue <= 0.20` |
| `muy_baja` | resto |

No es una probabilidad. Sirve para priorizar revisión.

### 3.2 Refinamiento σ54/σ38

Se aplica **después** del ranking y **nunca cambia el sigma**:

| Sigma | Regla (`refined_confidence="alta"`) |
|---|---|
| σ54 | `pvalue <= 0.05` **y** `margin >= 1.0` **y** `0 <= dist_to_start <= 60` |
| σ38 | `margin >= 1.0` **y** `0 <= dist_to_start <= 60` |

Si una fila tenía `confidence_label=alta` pero no pasa la regla, se rebaja a
`media`. Si pasa, se escribe `alta`. El resto de los sigma conservan su
`confidence_label` en `refined_confidence`. Se puede desactivar con
`--no-refine`.

---

## 4. Evaluación

Contra `data/sigma/gold/<genoma>/tss_sigma.tsv`, agrupado por
`(contig, tss, strand)`:

| Métrica | Definición |
|---|---|
| `top1_any` | el rank 1 está en el conjunto verdadero |
| `top1_all` | **todos** los sigma verdaderos están en el rank 1 (solo posible en single) |
| `top3_any` | algún sigma del top-3 está en el conjunto verdadero |
| `top3_all` | todos los sigma verdaderos están en el top-3 |
| `jaccard_topk` | promedio de (intersección / unión) por TSS |
| `micro P/R/F1` | suma de TP/FP/FN sobre todos los TSS y sigmas |
| `macro P/R/F1` | promedio por sigma con `support > 0` |
| single / multi | TSS con 1 etiqueta / con más de 1 |

Supuestos y límites de la métrica:

- El gold es **positivo**: un candidato en un TSS que no está en el gold no se
  cuenta como falso positivo; se reporta aparte como “fuera del gold”.
- Se asume que el gold es exhaustivo **para los sigma de ese TSS**. Si una
  etiqueta falta en el gold, se contará como FP; por eso el gold es curado.
- La precisión top-k incluye k etiquetas por TSS. No es comparable con la
  precisión top-1 de una fila por operón.
- `top1_all` es 0 por construcción en multi-sigma; la métrica útil ahí es
  `top3_all`.
- `macro` excluye etiquetas sin soporte en el gold; en Subtilis, σ24 tiene
  soporte 0.

---

## 5. Resultados K-12

Corrida del 2026-11-25 con top-3, sin umbrales de señal y refinamiento activo.
Salida: `data/sigma/multilabel/ecoli_k12/`.

### 5.1 Global

| Métrica | Valor |
|---|---:|
| TSS gold | 3.189 |
| Etiquetas gold | 3.314 |
| Single / multi | 3.065 / 124 |
| Candidatos en la tabla | 9.567 (3 por TSS) |
| `top1_any` | 1.819/3.189 (**57,0 %**) |
| `top1_all` | 1.737/3.189 (54,5 %) |
| `top3_any` | 2.878/3.189 (**90,2 %**) |
| `top3_all` | 2.835/3.189 (**88,9 %**) |
| Jaccard top-1 | 0,558 |
| Jaccard top-3 | 0,308 |
| Micro P/R/F1 (top-3) | 0,309 / 0,892 / 0,459 |
| Macro P/R/F1 (top-3) | 0,246 / 0,820 / 0,345 |

Como referencia, top-1: micro P/R/F1 = 0,570 / 0,549 / 0,559; macro =
0,454 / 0,512 / 0,456. La ganancia de recall del top-3 se paga con una caída
fuerte de precisión: el conjunto predicho tiene tres etiquetas por TSS.

### 5.2 Single vs multi

| Subconjunto | n | `top1_any` | `top1_all` | `top3_any` | `top3_all` | Jaccard top-3 |
|---|---:|---:|---:|---:|---:|---:|
| single | 3.065 | 1.737 (56,7 %) | 1.737 (56,7 %) | 2.757 (90,0 %) | 2.757 (90,0 %) | 0,300 |
| multi | 124 | 82 (66,1 %) | 0 (0,0 %) | 121 (97,6 %) | **78 (62,9 %)** | 0,508 |

Lectura: en multi-sigma el top-3 contiene **alguna** etiqueta verdadera en 121
de 124 TSS (97,6 %), pero cubre **todas** en solo 78 (62,9 %). El número global
de 88,9 % está dominado por los 3.065 TSS single.

### 5.3 Por sigma (conjunto top-3)

| Sigma | Support | Pred | TP | FP | FN | P | R | F1 | Support multi | Cubiertos multi | Fuera multi |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| σ70 | 1.997 | 2.677 | 1.847 | 830 | 150 | 0,690 | 0,925 | 0,790 | 119 | 99 | 20 |
| σ54 | 99 | 834 | 63 | 771 | 36 | 0,076 | 0,636 | 0,135 | 4 | 1 | 3 |
| σ32 | 312 | 1.123 | 268 | 855 | 44 | 0,239 | 0,859 | 0,374 | 23 | 18 | 5 |
| σ38 | 241 | 2.169 | 204 | 1.965 | 37 | 0,094 | 0,846 | 0,169 | 93 | 75 | 18 |
| σ28 | 143 | 1.058 | 109 | 949 | 34 | 0,103 | 0,762 | 0,182 | 3 | 3 | 0 |
| σ24 | 522 | 1.706 | 466 | 1.240 | 56 | 0,273 | 0,893 | 0,418 | 7 | 4 | 3 |

La precisión top-3 es baja por diseño: cada TSS aporta tres etiquetas, dos de
ellas potencialmente falsas. σ70 conserva la mejor F1; σ54 y σ38 son los casos
donde el ranking top-3 incluye muchos rivales débiles.

### 5.4 Cobertura multi-sigma por k

“No cubiertos del todo” = TSS multi donde falta al menos una etiqueta.
“Fuera del top-k (TSS)” = TSS multi sin ninguna etiqueta en el top-k.
“Etiquetas fuera” = etiquetas verdaderas que no entran en el top-k.

| k | TSS multi | Alguno cubierto | Todos cubiertos | No cubiertos del todo | Fuera del top-k (TSS) | Etiquetas cubiertas | Etiquetas fuera |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 124 | 82 (66,1 %) | 0 (0,0 %) | 124 | 42 | 82/249 (32,9 %) | 167 |
| 2 | 124 | 110 (88,7 %) | 39 (31,5 %) | 85 | 14 | 149/249 (59,8 %) | 100 |
| **3\*** | 124 | 121 (97,6 %) | 78 (62,9 %) | 46 | 3 | 200/249 (80,3 %) | **49** |
| 4 | 124 | 124 (100,0 %) | 103 (83,1 %) | 21 | 0 | 228/249 (91,6 %) | 21 |
| 5 | 124 | 124 (100,0 %) | 119 (96,0 %) | 5 | 0 | 244/249 (98,0 %) | 5 |
| 6 | 124 | 124 (100,0 %) | 124 (100,0 %) | 0 | 0 | 249/249 (100,0 %) | 0 |

\* Fila que entrega la tabla top-3. Los k=4..6 son sensibilidad sobre el
ranking completo de features; **no** se guardan como candidatos en la salida
por defecto.

Con top-3 quedan **46 de 124 TSS multi incompletos** y **49 de 249 etiquetas
fuera**. Tres TSS no tienen **ninguna** etiqueta correcta en el top-3:

| Contig | TSS | Hebra | Verdad | Top-3 predicho |
|---|---|---:|---|---|
| NC_000913.3 | 1.330.789 | + | σ70, σ38 | σ54, σ32, σ28 |
| NC_000913.3 | 1.937.494 | + | σ70, σ38 | σ24, σ28, σ32 |
| NC_000913.3 | 3.039.851 | - | σ70, σ32 | σ24, σ54, σ28 |

Con top-6 se cubren los 124 TSS y las 249 etiquetas, pero se agregan tres
candidatos por TSS y la precisión baja aún más. Es un intercambio explícito,
no una mejora gratuita.

### 5.5 Refinamiento σ54/σ38

Se evalúa sobre **todas** las filas de señal del sigma (no solo las que entran
en top-3), contra el gold. Es un análisis del filtro, en el mismo gold; no es
validación *held-out*.

| Sigma | Regla | n pass | TP | FP | FN | P | R | F1 | Alta base P/R | Subidas | Bajadas |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| σ54 | p≤0,05 y margen≥1 y dist≤60 | 35 | 21 | 14 | 78 | **0,600** | 0,212 | 0,313 | 0,568 / 0,212 | 0 | 2 |
| σ38 | margen≥1 y dist≤60 | 191 | 60 | 131 | 181 | **0,314** | 0,249 | 0,278 | 0,309 / 0,245 | 7 | 7 |

σ54 sube precisión de 0,568 a 0,600 sin cambiar recall. σ38 apenas mueve la
precisión (0,309 → 0,314). Son filtros de revisión, no una nueva llamada de
sigma.

### 5.6 Sensibilidad al umbral de señal

Ejemplo con top-3, variando `--min-z` (no se escribe otra salida; es la misma
features evaluada en memoria):

| `min_z` | Candidatos | TSS con candidatos | `top3_any` | `top3_all` | Multi con todos | Etiquetas multi fuera | Micro P | Micro R |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ninguno | 9.567 | 3.189 | 2.878 | 2.835 | 78 | 49 | 0,309 | 0,892 |
| 0,0 | 9.014 | 3.182 | 2.828 | 2.784 | 73 | 58 | 0,322 | 0,876 |
| 0,5 | 7.689 | 3.135 | 2.691 | 2.642 | 60 | 79 | 0,358 | 0,830 |
| 1,0 | 5.740 | 2.937 | 2.410 | 2.355 | 41 | 112 | 0,427 | 0,740 |
| 2,0 | 1.876 | 1.536 | 1.055 | 1.015 | 7 | 195 | 0,566 | 0,320 |

Filtrar sube precisión, pero **pierde TSS y etiquetas verdaderas** porque
muchas etiquetas curadas tienen z bajo. Por eso la salida principal usa
`min_z=None` y deja la confianza como columna, no como filtro.

---

## 6. Resultados *B. subtilis*

Corrida del 2026-11-25 con top-3, sin umbrales y refinamiento activo.
Salida: `data/sigma/multilabel/subtilis/`.

### 6.1 Global

| Métrica | Valor |
|---|---:|
| TSS gold | 2.226 |
| Etiquetas gold | 2.226 |
| Single / multi | 2.226 / **0** |
| Candidatos en la tabla | 6.678 |
| `top1_any` = `top1_all` | 997/2.226 (**44,8 %**) |
| `top3_any` = `top3_all` | 1.904/2.226 (**85,5 %**) |
| Jaccard top-1 / top-3 | 0,448 / 0,285 |
| Micro P/R/F1 (top-3) | 0,285 / 0,855 / 0,428 |
| Macro P/R/F1 (top-3) | 0,212 / 0,674 / 0,237 |

El gold de Subtilis no tiene TSS multi-sigma, así que `top3_all` es igual a
`top3_any` y no hay cobertura multi que reportar.

### 6.2 Por sigma (conjunto top-3)

| Sigma | Support | Pred | TP | FP | FN | P | R | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| σ70 (A) | 1.949 | 1.926 | 1.740 | 186 | 209 | 0,903 | 0,893 | 0,898 |
| σ54 (L) | 11 | 565 | 9 | 556 | 2 | 0,016 | 0,818 | 0,031 |
| σ32 (I) | 13 | 731 | 6 | 725 | 7 | 0,008 | 0,462 | 0,016 |
| σ38 (B) | 182 | 1.661 | 105 | 1.556 | 77 | 0,063 | 0,577 | 0,114 |
| σ28 (D) | 71 | 631 | 44 | 587 | 27 | 0,070 | 0,620 | 0,125 |
| σ24 | 0 | 1.164 | 0 | 1.164 | 0 | 0,000 | — | — |

σ24 no tiene soporte en el gold de Subtilis; no entra en el macro. El patrón es
el esperado por transferencia: las PWMs son de *E. coli* y σ70 concentra casi
toda la señal. No es un problema del formato multietiqueta.

### 6.3 Refinamiento

| Sigma | n pass | TP | FP | FN | P | R | F1 | Alta base P/R | Subidas | Bajadas |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| σ54 | 15 | 1 | 14 | 10 | 0,067 | 0,091 | 0,077 | 0,056 / 0,091 | 0 | 3 |
| σ38 | 131 | 2 | 129 | 180 | 0,015 | 0,011 | 0,013 | 0,018 / 0,011 | 23 | 4 |

La regla **no mejora** en Subtilis. Las 23 “subidas” de σ38 a `alta` no se
traducen en precisión (P=0,015); son un artefacto de transferencia, no
evidencia. Para este genoma conviene `--no-refine` si se va a interpretar
`refined_confidence` como filtro.

---

## 7. Salidas

```text
data/sigma/multilabel/<genoma>/
  sigma_candidates.tsv   # una fila por (genoma, TSS/operón, sigma candidato)
  metrics.json           # configuración, conteos, métricas, cobertura, límites
  report.md              # reporte legible con tablas de esa corrida
```

### 7.1 `sigma_candidates.tsv`

| Columna | Contenido |
|---|---|
| `genome`, `contig`, `tss`, `strand` | TSS real (1-based) o vacío en respaldo de operón |
| `locus_id` | `contig:tss:strand` o `operon_id` |
| `coord_type` | `tss` o `operon_anchor` |
| `operon_anchor` | solo respaldo; **no es un TSS** |
| `rank` | 1..top-k, ordenado por `z` |
| `sigma` | etiqueta canónica |
| `z`, `pvalue`, `margin`, `score` | señal PWM |
| `spacer`, `spacer_dev`, `dist_to_start`, `inventory` | arquitectura/contexto disponible |
| `confidence_label` | heurística de la sección 3.1 |
| `refined_confidence` | resultado de la regla σ54/σ38; no cambia `sigma` |
| `refinement_rule` | regla aplicada, o vacío |
| `evidence` | texto de señal; no contiene etiquetas del gold |
| `source` | `layers` o `promoters` |

Tamaños de esta corrida: 9.568 líneas (con encabezado) en K-12 y 6.679 en
Subtilis.

### 7.2 `metrics.json`

Incluye `version`, `generated_at`, `genome`, `source`, `config`, `inputs`,
`candidate_stats`, `gold_stats`, `evaluation`, `refinement`, `limitations` y
`outputs`. Las métricas requeridas están en:

- `evaluation.overall`: `top1_any`, `top1_all`, `top3_any`, `top3_all`,
  Jaccard, micro, macro y `by_k`;
- `evaluation.single` / `evaluation.multi`;
- `evaluation.per_sigma`;
- `evaluation.multi_sigma_coverage`;
- `evaluation.topk_sensitivity`;
- `refinement.per_sigma`.

---

## 8. Uso

```bash
# Corrida principal (genera las salidas de K-12 y Subtilis)
.venv/bin/python scripts/sigma/74_multilabel.py \
  --genome ecoli_k12 --genome subtilis

# Otro top-k; la tabla guarda k candidatos y las métricas top-1/top-3 siguen
# disponibles si k >= 3
.venv/bin/python scripts/sigma/74_multilabel.py --genome ecoli_k12 --top-k 6

# Umbral de señal (pierde etiquetas verdaderas; usar con cuidado)
.venv/bin/python scripts/sigma/74_multilabel.py --genome ecoli_k12 \
  --min-z 0 --max-pvalue 0.10

# Respaldo de operón; usar otro out-root para no pisar la salida TSS
.venv/bin/python scripts/sigma/74_multilabel.py --genome ecoli_k12 \
  --source promoters --out-root data/sigma/multilabel_promoters

# Pruebas rápidas
.venv/bin/python -m pytest tests/test_sigma_multilabel.py -q
```

También acepta `--all` (todos los genomas con features de capa), `--no-refine`,
`--min-pwm-score`, `--features-root`, `--gold-root` y `--out-root`.

---

## 9. Límites

1. **Cobertura anclada al gold.** Las features de capa existen para los TSS del
   gold (3.189 y 2.226). La tabla no es un escaneo genómico independiente y no
   puede estimar precisión sobre TSS fuera del gold.
2. **Sin TSS nuevos.** No se generan coordenadas nuevas. El respaldo de
   `sigma_promoters.tsv` es a nivel de operón y no se evalúa como TSS.
3. **Mismo gold para diseñar y evaluar.** Los umbrales de confianza y las
   reglas σ54/σ38 se inspiraron en estas features/gold. Las métricas son
   descriptivas, no una validación *held-out*.
4. **Precisión top-k inflada de falsos positivos.** Cada TSS aporta k etiquetas.
   La precisión micro top-3 (0,309 en K-12) no es comparable con la precisión
   top-1 (0,570).
5. **El titular global esconde multi-sigma.** `top3_all` global = 88,9 %, pero
   en multi-sigma = 62,9 %. Reportar siempre ambos.
6. **Tres TSS multi sin ninguna etiqueta en top-3.** K-12: 1.330.789 (+),
   1.937.494 (+) y 3.039.851 (-). Recuperarlos exige otro ranking, features
   adicionales o top-6, no un umbral más laxo.
7. **Gold positivo y supuesto exhaustivo.** Un candidato fuera del gold no se
   cuenta como FP. Un sigma ausente en el gold se cuenta como FP si se predice;
   el gold curado reduce, pero no elimina, ese riesgo.
8. **Subtilis es transferencia.** Las PWMs son de *E. coli*; σ54/σ32/σ38
   tienen precisión muy baja. El refinamiento no ayuda allí.
9. **Confianza no calibrada.** `confidence_label` y `refined_confidence` son
   heurísticas, no probabilidades. `refined_confidence` no cambia `sigma` ni
   demuestra mecanismo.
10. **Sin red, sin instalación, sin GUI.** Este script no modifica el núcleo ni
    el pipeline; solo lee features/gold y escribe bajo
    `data/sigma/multilabel/`.

---

## 10. Recomendación

- Usar `sigma_candidates.tsv` **top-3** como conjunto de candidatos para
  revisión y para downstream multietiqueta. No reemplaza la asignación
  top-1 de `sigma_promoters.tsv`; la complementa.
- Si el objetivo es **cobertura exhaustiva** de multi-sigma, guardar top-6
  (124/124 TSS y 249/249 etiquetas) o, al menos, revisar los 46 TSS multi
  incompletos y los 3 TSS sin acierto del top-3. El costo es más FP.
- Usar `refined_confidence="alta"` para σ54/σ38 en **K-12** como filtro de
  revisión: σ54 sube a P=0,600; σ38 queda en P=0,314. No usarlo en Subtilis
  sin recalibrar.
- No filtrar por `min_z` en la salida principal: pierde TSS y etiquetas
  verdaderas (sección 5.6). Si se filtra, reportar cobertura y no solo precisión.
- Próximo paso real: evaluar en un conjunto **held-out** de TSS o en un escaneo
  genómico independiente. Mientras las features estén ancladas al gold, estas
  métricas no demuestran generalización.
