# Evaluación fina por capas del pipeline sigma

> **Principio:** un accuracy global mezcla datos, localización, motivo,
> competición, contexto e integración. Si algo “falla mucho”, puede ser que
> fallen tres capas a la vez o que haya procesos biológicos superpuestos.
>
> **Regla:** no se declara falso positivo sin contexto. Un TSS fuera del gold
> puede ser un TSS nuevo, un proceso superpuesto o una etiqueta incompleta.

---

## 1. Capas

| Capa | Pregunta | Qué mide | Qué NO mide |
|---|---|---|---|
| **0. Datos** | ¿hay TSS, secuencia, PWM y etiqueta? | cobertura y denominadores | acierto |
| **1. Localización** | ¿el TSS correcto tiene mejor señal que señuelos? | soporte del sitio | identidad del sigma |
| **2. Motivo** | ¿el PWM del sigma verdadero está soportado? | señal por sigma | competición entre sigma |
| **3. Asignación** | ¿el top-1/top-3 coincide? | acierto y multietiqueta | contexto biológico |
| **4. Contexto** | ¿hay inventario/ruta/expresión? | disponibilidad y calibración | verdad experimental |
| **5. Integración** | ¿el operón entra en el regulón? | enlace downstream | mecanismo |

El benchmark actual (`30_validate_sigma.py`) cubre sobre todo **3** con algo de
**0/2**. `50_eval_layers.py` separa las capas.

---

## 2. Gold esperado

Ruta:

```text
data/sigma/gold/<id>/tss_sigma.tsv
```

Columnas aceptadas (sin orden fijo):

```text
contig, tss/pos/position, strand, sigma/sigma_factor/factor,
operon/operon_id, gene/first_gene, evidence, confidence
```

- Una fila por TSS y sigma.
- Si un TSS tiene varios sigma, se agrupan y se evalúa como **multietiqueta**.
- `?`, `NA`, vacío → `unlabeled`; no cuentan como fallo.
- Sigma sin PWM → `no_pwm`; no cuentan como fallo.
- `<id>` debe ser el mismo ID de `data/regulation/<id>/` (p. ej. `ecoli_k12`,
  `subtilis`, `GCF_...`).

Estado del gold:

| Especie | Filas | Cobertura | Uso |
|---|---:|---|---|
| K-12 | 3.314 | 6 sigma | **Benchmark completo** |
| B. subtilis | 2.226 | 5 sigma | **Benchmark completo** |
| M. tuberculosis | 34 | 4/13 | Parcial; no accuracy global |
| S. coelicolor | 29 | 10 sigma | Parcial |
| S. aureus | 25 | 3 sigma | Parcial |
| P. aeruginosa PAO1 | 18 | 4 sigma | Exploratorio σE; no accuracy global |
| C. glutamicum | 32 | en curación | Parcial |

`--all` usa solo los genomas declarados completos en `schema.json`; los
parciales se omiten salvo `--include-partial`. No mezclar PAO1 con PA14 ni usar
las 18 filas de PAO1 para estimar σ54/σ38/σ70.

La evaluación oficial es:

```bash
python scripts/sigma/50_eval_layers.py --gold-root data/sigma/gold --all
```

Los adaptadores quedan solo como respaldo si falta el gold:

```bash
python scripts/sigma/50_eval_layers.py --genome ecoli_k12 --source k12
python scripts/sigma/50_eval_layers.py --genome subtilis --source subtilis
```

---

## 3. Resultados con el gold curado

Corrida del 2026-11-25 sobre `data/sigma/gold/`, ventana −80/+20, 200
señuelos, semilla 17. El gold ya está disponible; los adaptadores K-12/DBTBS
quedan solo como respaldo.

### 3.1 K-12

| Capa | Métrica | Valor |
|---|---|---:|
| 0 | TSS gold | 3.189 |
| 0 | Sin etiqueta | 0 |
| 0 | Sin PWM | 0 |
| 1 | TSS con soporte de motivo | 1.578/3.189 (**49,5 %**) |
| 3 | Evaluables | 3.189 |
| 3 | Single / multi | 3.065 / **124** |
| 3 | Top-1 | **57,0 %** |
| 3 | Top-3 cualquiera | **90,2 %** |
| 3 | Top-3 todos | 88,9 % |
| 3 | Mayoría trivial (σ70) | 62,6 % |
| 3 | `correct` / `partial_top3` / `missed` | 1.737 / 1.059 / 311 |
| 3 | `hit_one_of_multi` | 82 |

Por sigma (top-1):

| Sigma | Support | Precision | Recall | F1 |
|---|---:|---:|---:|---:|
| σ70 | 1.997 | 0,863 | 0,576 | 0,691 |
| σ24 | 522 | 0,565 | 0,625 | 0,593 |
| σ32 | 312 | 0,569 | 0,596 | 0,582 |
| σ28 | 143 | 0,268 | 0,497 | 0,348 |
| σ38 | 241 | 0,231 | 0,710 | 0,348 |
| σ54 | 99 | 0,169 | 0,333 | 0,224 |

### 3.2 *B. subtilis*

| Capa | Métrica | Valor |
|---|---|---:|
| 0 | TSS gold | 2.226 |
| 0 | Sin etiqueta | 0 |
| 0 | Sin PWM | 0 |
| 1 | TSS con soporte de motivo | 913/2.226 (**41,0 %**) |
| 3 | Evaluables | 2.226 |
| 3 | Single / multi | 2.226 / 0 |
| 3 | Top-1 | **44,8 %** |
| 3 | Top-3 cualquiera | 85,5 % |
| 3 | Mayoría trivial (A→σ70) | 87,6 % |
| 3 | `correct` / `partial_top3` / `missed` | 997 / 907 / 322 |

Por sigma (top-1):

| Sigma | Support | Precision | Recall | F1 |
|---|---:|---:|---:|---:|
| σ70 (A) | 1.949 | 0,929 | 0,487 | 0,639 |
| σ38 (B) | 182 | 0,050 | 0,154 | 0,075 |
| σ28 (D) | 71 | 0,122 | 0,211 | 0,155 |
| σ54 (L) | 11 | 0,027 | 0,273 | 0,048 |
| σ32 (I) | 13 | 0,017 | 0,154 | 0,030 |

**Lectura:** K-12 queda **5,6 pp por debajo** de la mayoría trivial (57,0 % vs
62,6 %). El top-3 es alto, pero está dominado por σ70. *B. subtilis* cae a
44,8 % porque las PWMs son de *E. coli*; no es un fallo del mismo tipo: es
transferencia entre especies.

> **Corrección (2026-11-25).** Confirmado: el 44,8 % de *B. subtilis* es
> transferencia de *E. coli* (el flag `--pwm` se ignoraba). Con la PWM
> propia el top-1 es **29,7 %** (top-3 80,7 %); con la capa de contexto
> asciende a **43,4 %** (top-3 86,9 %), todavía por debajo de la mayoría
> trivial. Detalle, ablación y límites en
> [`SIGMA_CONTEXT.md`](SIGMA_CONTEXT.md).

**Lo que el gold ya permite separar:**

- 124 TSS de K-12 tienen **más de un sigma**; un top-1 único no puede acertarlos
  todos. Por eso se reporta `hit_one_of_multi` (82) y `top3_all`.
- La capa 1 (49,5 %) es más baja que la capa 3 (57,0 %): el motivo no siempre
  supera a los señuelos, aunque la competición entre sigma acierte.
- La capa 2 muestra que σ54 y σ38 son los motivos con menor soporte
  (28,3 % y 49,4 %); no es lo mismo que “el asignador falle”.

## 4. Multietiqueta y procesos superpuestos

> Combinaciones de features para refinar: [`SIGMA_FEATURE_COMBINATIONS.md`](SIGMA_FEATURE_COMBINATIONS.md).
> Iteración 2 (localización, motivo, multietiqueta): [`SIGMA_ITERATION2.md`](SIGMA_ITERATION2.md).

El gold puede tener varios sigma por TSS. El script evalúa:

- `top1_any`: el top-1 está en el conjunto verdadero;
- `top3_any`: algún top-3 está en el conjunto;
- `top3_all`: todos los sigma verdaderos están en top-3;
- `jaccard_top3`: intersección/unión;
- por sigma: top-1 para single-label; top-3 para multi-label.

Un TSS con σ70 y σ32 no es un “error” si el modelo predice uno de los dos.
Tampoco es un acierto completo: hay que reportar `hit_any` y `hit_all` por
separado.

---

## 5. Falso positivo vs nuevo hallazgo

El gold suele contener **positivos**, no negativos. Por eso:

| Situación | Etiqueta |
|---|---|
| Predicción coincide con gold | `correct` |
| Coincide uno de varios sigma | `partial` / `hit_one_of_multi` |
| Gold tiene el TSS pero no ese sigma | `conflict` (solo si el gold es exhaustivo) |
| TSS fuera del gold, con señal fuerte y contexto | `novel_candidate` |
| TSS fuera del gold, sin contexto | `outside_gold` / `unresolved` |
| Gold cubre la región y hay evidencia negativa | `false_positive` |

**No se declara `false_positive` solo porque no esté en el gold.** Para eso
hace falta que el gold cubra esa región y que haya evidencia negativa
(secuencia, conservación, contexto). Si no, es `unresolved` o
`novel_candidate`.

El script actual reporta `correct`, `partial_top3`, `missed`, `unlabeled`,
`no_pwm`; la clasificación `novel_candidate`/`false_positive` queda para cuando
el gold traiga cobertura negativa.

---

## 6. Salidas

```text
data/sigma/layers/<genoma>/
  layer_metrics.json       # capas 0-4, por sigma, clases
  layer_predictions.tsv    # una fila por TSS gold
  layer_report.md          # informe legible
```

`layer_predictions.tsv` incluye:

```text
contig, tss, strand, true_sigmas, pred_top1, pred_top3, n_true,
class, jaccard_top3, no_pwm, has_pwm, operons, genes, evidence, confidence
```

---

## 7. Criterios para el gold nuevo

Cuando los agentes agreguen `tss_sigma.tsv`, verificar:

- [ ] `tss` es 1-based y la hebra es `+`/`-`.
- [ ] `sigma` usa nombres canónicos (`sigma70`, `sigma54`, …) o letras de
      *B. subtilis* documentadas.
- [ ] Varios sigma para un TSS van en filas separadas o separados por `;`.
- [ ] `unlabeled`/`?` se conserva, no se borra.
- [ ] `operon`/`gene` permiten cruzar con `data/regulation/`.
- [ ] Hay una fila `README` o `manifest` que declare fuente, versión y
      licencia.

---

## 8. Frase final

**La evaluación por capas separa “no hay datos” de “no hay PWM”, “el motivo no
se ve”, “el motivo se ve pero pierde la competición” y “el contexto no está
disponible”. Solo con esa separación se puede distinguir un fallo real de un
proceso superpuesto o de un hallazgo nuevo.**
