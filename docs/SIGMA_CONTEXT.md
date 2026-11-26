# Capa de contexto sigma: diseño, ablación y límites

> **Qué es:** una capa determinista y auditable que añade evidencias
> contextuales (operón, ruta/EC, dominio, conservación, expresión, geometría
> σ54) a la asignación de factor sigma. Cada evidencia declara fuente, peso y
> cobertura; la combinación es un **log-odds** transparente, no un modelo
> opaco.
>
> **Qué no es:** no sustituye al PWM, no usa etiquetas del gold para puntuar
> test, no rellena datos de otra especie y no convierte una ausencia de
> evidencia en una etiqueta.
>
> **Alcance:** la capa es **TSS-level** y vive en
> `scripts/sigma/sigma_context.py` + `scripts/sigma/51_eval_context.py`.
> No cambia la asignación por operón de `sigma_annotator.annotate_genome`
> (que sigue usando PWM + inventario + ruta/expresión). Se evalúa contra el
> gold y se deja como capa opcional: en K-12 el selector la apaga; en
> B. subtilis la enciende. Integrarla en producción exigiría fijar los
> pesos/escala y volver a correr el pipeline, que no se hace aquí.
>
> **Resultado corto**
>
> - **K-12:** la capa **no mejora** el top-1 fuera de muestra. El selector
>   out-of-fold la apaga sola (escala 0) y el resultado queda igual:
>   **56,98 % → 56,98 %** (top-3 90,25 % → 90,25 %). La única mejora
>   medible es de **calibración de confianza** (ECE 0,061 → 0,053) sin
>   cambiar top-1.
> - **B. subtilis:** con **PWMs propias** (sin transferir E. coli), la capa
>   sube el top-1 de **29,69 % → 43,40 %** (+13,70 pp; IC 95 % bootstrap
>   [11,86; 15,68] pp) y el top-3 de **80,68 % → 86,88 %**. Sigue **por
>   debajo de la mayoría trivial (87,56 %)** en top-1.
> - **Hallazgo de auditoría:** el **44,79 %** que se atribuía a las PWMs
>   propias de *B. subtilis* era en realidad una corrida con **PWMs de
>   E. coli**: `50_eval_layers.py` parseaba `--pwm` pero lo ignoraba. Con la
>   matriz propia el top-1 es **29,69 %**. Ver §7.

---

## 1. Capas y principio

La evaluación existente separa datos, localización, motivo, asignación,
contexto e integración ([`SIGMA_LAYERED_EVAL.md`](SIGMA_LAYERED_EVAL.md)).
Esta capa cubre la **capa 4** con una regla explícita:

```text
score_final(sigma) = z_PWM(sigma) + escala × Σ_e peso_e × LLR_e(sigma)
```

donde `LLR_e` es el log-likelihood ratio de la evidencia `e` para ese sigma,
estimado **solo con las particiones de entrenamiento**.

Tres guardarraíles evitan forzar la etiqueta:

1. **`escala = 0` reproduce el ranking del PWM** exactamente.
2. **La escala se elige en entrenamiento** (por fold); si el contexto no
   mejora el top-1 de entrenamiento en ≥ 1 pp, queda en 0 y la capa se apaga.
3. **Modo acotado (`gap = 1,0`):** el contexto solo puede promover candidatos
   que estén a ≤ 1 unidad z del mejor motivo. El top-3 se calcula con el
   ranking completo, no con el conjunto acotado.

---

## 2. Evidencias, fuentes y pesos

| Evidencia | Fuente local | Cobertura K-12 | Cobertura B. subtilis | Peso declarado |
|---|---|---:|---:|---:|
| `operon_position` | `data/regulation/<g>/genes.tsv` (TSS → gen aguas abajo) | 100 % | 100 % | 1,0 |
| `leader` | misma tabla, distancia al gen aguas abajo (bins) | 100 % | 100 % | 0,5 |
| `operon_size` | `data/regulation/<g>/operons.tsv`; B. subtilis: `data/dbtbs/gt_operons.tsv` | 100 % | 100 % | 0,5 |
| `pathway_ec` | `genes.tsv:ec` (primer dígito EC) | 27,5 % | 0 % | 1,0 |
| `domain` | `data/sigma/structures/structure_features.tsv` | 0 % | 0 % | 0,5 (no disponible) |
| `conservation` | `synteny_pairs.tsv` + `genes.tsv:acc`; B. subtilis: `bsu_vs_k12.tsv` | 23,4 % | 4,5 % (insuficiente) | 1,0 |
| `expression` | `data/expression/<g>.tsv` (opcional) | 0 % | 0 % | 1,0 (no disponible) |
| `sigma54_geometry` | PWM σ54 + ventana evaluada | 100 % (99 sitios) | 100 % (11 sitios) | 1,0 (σ54) |

Reglas de activación:

- Una evidencia sin fuente local o con **cobertura < 5 %** recibe peso 0:
  se conserva en la salida, pero no mueve el ranking.
- `domain` **no** usa `genes.tsv:fam`: en K-12 hay 4212 valores de `fam`
  para 4262 genes (casi todos singletons). No es un dominio reutilizable y
  sobreajustaría; se declara no disponible.
- `expression` no usa `expression_rules.tsv`: ese fichero contiene **reglas
  de regulación**, no niveles de expresión. No hay RNA-seq en el repo.
- `sigma54_geometry` solo se usa como prior si la PWM σ54 tiene **≥ 20
  sitios**. K-12: 96 → activa (pero el selector la deja en 0). B. subtilis:
  11 → **diagnóstico**, no mueve el ranking. La geometría de B. subtilis es
  una caja única; no hay dos cajas −24/−12 reales.

### 2.1 Por qué `none` es neutral

Si la ausencia de una evidencia se modela como una categoría más, el
clasificador aprende de **qué falta**, no de lo que se midió. En K-12 eso
producía una subida aparente de +1,1 pp que desaparecía al neutralizar
`none`. Por eso:

- `none` aporta `LLR = 0` en scoring;
- las categorías con soporte < `min_support` (5 sitios de entrenamiento)
  también aportan 0;
- la ablación usa pesos **excluyentes** (un dict explícito no vuelve a sumar
  las demás evidencias con su peso por defecto).

---

## 3. Matemática del log-odds

Para cada evidencia `e`, sigma `s` y valor `v` observado en entrenamiento:

```text
P(v | s) = (n_s(v) + α) / (N_s + α·(|V| + 1))
P(v)     = (n(v)   + α) / (N   + α·(|V| + 1))       α = 1 (Dirichlet)
LLR_e(s, v) = log( P(v | s) / P(v) )
```

- Un sitio multietiqueta cuenta para cada sigma verdadero en el numerador y
  una vez en el fondo.
- `N_s` es el número de sitios de entrenamiento con ese sigma; si es 0, la
  tabla de ese sigma queda vacía.
- La combinación asume independencia condicional entre evidencias (naive
  Bayes). No se estima una matriz de covarianza: sería un modelo opaco y
  sobreajustaría con 2.000–3.000 sitios.
- `α = 1` y `min_support = 5` son constantes declaradas, no ajustadas al
  test.

La Tabla de LLR se ajusta **por fold** (5 bloques contiguos por coordenada).
El modelo de test nunca ve el fold de test.

---

## 4. Protocolo de evaluación

| Elemento | K-12 | B. subtilis |
|---|---|---|
| Gold | `data/sigma/gold/ecoli_k12/tss_sigma.tsv` (3.314 filas, 3.189 TSS) | `data/sigma/gold/subtilis/tss_sigma.tsv` (2.226 filas/TSS) |
| PWM | `data/regulondb/gt/promoter_pwms.tsv` (RegulonDB) | `data/dbtbs/promoter_pwms_subtilis_all.tsv` (DBTBS propio) |
| Ventana | −80/+20 | −80/+20 (sensibilidad −40/+5) |
| Señuelos | 200 (sensibilidad 2.000) | 200 (sensibilidad 2.000) |
| Folds | 5 bloques contiguos | 5 bloques contiguos |
| Semilla | 17 | 17 |
| Escala | elegida por fold en entrenamiento | elegida por fold en entrenamiento |
| Calibración | temperatura por fold, ajustada en entrenamiento | ídem |

No se transfieren PWMs de E. coli a B. subtilis. El comando `--all` resuelve
la PWM propia de cada especie; si se pasa una PWM de `regulondb` a
`--genome subtilis`, el resultado se marca `ecoli_transfer=true` y no se usa
como número principal.

---

## 5. K-12: antes / después

Corrida oficial (`data/sigma/context/ecoli_k12/`):

| Capa | Top-1 | Top-3 cualquiera | Top-3 todos | Mayoría trivial |
|---|---:|---:|---:|---:|
| Motivo (línea base) | **56,98 %** | **90,25 %** | 88,90 % | 62,62 % (σ70) |
| Contexto (escala elegida = 0) | **56,98 %** | **90,25 %** | 88,90 % | 62,62 % (σ70) |
| Δ | 0,00 pp | 0,00 pp | 0,00 pp | — |

**Lectura:** el selector out-of-fold eligió **escala 0 en los cinco folds**.
No es que la capa no se haya probado: se probó y se apagó porque no mejoró
en entrenamiento. El top-1 y el top-3 quedan exactamente igual.

> La corrida de capas reporta 57,04 %; esta corrida reporta 56,98 %. La
> diferencia es de 2 TSS por desempate/redondeo, no por un cambio de datos.

### 5.1 Ablación (escala fija 1,0, fuera de muestra)

| Evidencia | Top-1 | Δ top-1 | Top-3 | Cambios | Aciertos | Errores |
|---|---:|---:|---:|---:|---:|---:|
| solo `operon_position` | 56,32 % | −0,66 pp | 90,34 % | 467 | 160 | 181 |
| solo `leader` | 56,88 % | −0,09 pp | 90,15 % | 156 | 51 | 54 |
| solo `operon_size` | 56,69 % | −0,28 pp | 90,44 % | 89 | 26 | 35 |
| solo `pathway_ec` | **57,04 %** | **+0,06 pp** | 90,37 % | 109 | 37 | 35 |
| solo `conservation` | 56,85 % | −0,13 pp | 90,25 % | 58 | 18 | 22 |
| solo `sigma54_geometry` | 56,98 % | 0,00 pp | 90,25 % | 1 | 0 | 0 |
| todas (1,0) | 56,19 % | −0,78 pp | 89,75 % | 594 | 199 | 224 |

Ninguna evidencia supera el ruido. `pathway_ec` es la única con signo
positivo (+0,06 pp), pero su cobertura es 27,5 % y el IC no se separa de 0.
Las evidencias con cobertura 100 % (`operon_position`, `leader`,
`operon_size`) **empeoran** solas y combinadas: la asociación con sigma es
débil y está confundida por el tipo de TSS.

### 5.2 Barrido de escala

| Escala | Top-1 | Δ top-1 | IC 95 % Δ | Top-3 |
|---:|---:|---:|---:|---:|
| 0,00 | 56,98 % | 0,00 pp | [0,00; 0,00] | 90,25 % |
| 0,25 | 56,76 % | −0,22 pp | [−0,91; +0,41] | 90,56 % |
| 0,50 | **57,07 %** | **+0,09 pp** | **[−0,82; +1,07]** | **90,56 %** |
| 0,75 | 56,66 % | −0,32 pp | [−1,47; +0,72] | 90,53 % |
| 1,00 | 56,19 % | −0,79 pp | [−2,04; +0,44] | 89,75 % |

El mejor punto (+0,09 pp) tiene un IC que incluye 0. No se elige 0,5 porque
el barrido se hace sobre las predicciones out-of-fold; elegirlo sería
ajustar sobre test. El selector de entrenamiento prefiere 0.

### 5.3 Calibración

| Capa | ECE | Brier | NLL | Confianza media |
|---|---:|---:|---:|---:|
| Motivo | 0,0529 | 0,1029 | 1,1498 | 0,549 |
| Contexto (escala 0) | 0,0529 | 0,1029 | 1,1498 | 0,549 |

No hay cambio en la vía principal. **Modo confianza-only** (exploratorio):
se conserva el top-1 del motivo y se usa el contexto solo para recalibrar su
confianza. La escala (0,5) y la temperatura se eligen en entrenamiento.

| Modo | ECE | Brier | NLL | Confianza media | Acierto |
|---|---:|---:|---:|---:|---:|
| Motivo | 0,0613 | 0,2359 | 0,6686 | 0,549 | 56,98 % |
| Confianza con contexto | **0,0531** | **0,2307** | **0,6580** | 0,555 | 56,98 % |

Es la única mejora K-12: **calibración de la confianza, no top-1**. El ECE
baja 0,8 pp sin tocar la predicción. Se reporta como exploratorio porque el
modo no se usa en la vía principal.

### 5.4 Multietiqueta

- Evaluables: 3.189 (3.065 single, **124 multi**).
- Top-1 cualquiera: 56,98 % (igual antes/después).
- Top-3 cualquiera: 90,25 %; top-3 todos: 88,90 %.
- Multi top-1 (alguno): 66,1 %; multi top-1 exacto: 0,0 % (un top-1 no puede
  cubrir varios sigma).
- Jaccard top-3: 0,308.
- Clases: 1.735 `correct`, 82 `hit_one_of_multi`, 1.061 `partial_top3`,
  311 `missed`.

### 5.5 Por sigma (top-1, línea base = contexto en K-12)

| Sigma | Support | Precisión | Recall | F1 |
|---|---:|---:|---:|---:|
| σ70 | 1.997 | 0,860 | 0,544 | 0,667 |
| σ24 | 522 | 0,610 | 0,621 | 0,615 |
| σ32 | 312 | 0,598 | 0,567 | 0,582 |
| σ38 | 241 | 0,188 | 0,531 | 0,278 |
| σ28 | 143 | 0,285 | 0,483 | 0,358 |
| σ54 | 99 | 0,182 | 0,323 | 0,233 |

La capa no cambia ninguna de estas filas (escala 0). Esta tabla usa
**top-1 estricto en todos los TSS**, incluidos los 124 multietiqueta. La
tabla de `SIGMA_LAYERED_EVAL.md` usa top-3 en los multietiqueta (para no
penalizar un top-1 único), por eso sus recalls de σ38/σ24 son mayores; no es
una discrepancia de datos. Los valores también difieren en <0,1 pp por
desempate/redondeo.

### 5.6 σ54

- σ54 verdaderos: 99; predichos top-1: 176.
- Precisión 18,2 %; recall top-1 32,3 %; recall top-3 63,6 %.
- Geometría: 99/99 verdaderos y 176/176 predichos caen en la geometría
  canónica (el escáner de dos cajas exige espaciador 10–11 nt, −24/−12).
- **Lectura:** la geometría no es el cuello de σ54 en K-12; el problema es
  la competición con los demás sigma (precisión 18 %). Añadir una evidencia
  de geometría no cambia el top-1 (1 cambio, 0 aciertos, 0 errores).

---

## 6. B. subtilis: antes / después

Corrida oficial (`data/sigma/context/subtilis/`), **PWM propia**, sin
transferencia de E. coli:

| Capa | Top-1 | Top-3 cualquiera | Top-3 todos | Mayoría trivial |
|---|---:|---:|---:|---:|
| Motivo (PWM propia) | **29,69 %** | **80,68 %** | 80,68 % | 87,56 % (σ70) |
| Contexto (escala 1,0) | **43,40 %** | **86,88 %** | 86,88 % | 87,56 % (σ70) |
| Δ | **+13,70 pp** | **+6,20 pp** | +6,20 pp | — |

IC bootstrap del Δ top-1: **[+11,86; +15,68] pp** (p bilateral ≈ 0,000,
1.000 remuestreos). Escalas elegidas: 1,0 en los cinco folds.

**Lectura:** el contexto recupera buena parte de lo que la PWM de una caja
no ve. Aun así, el top-1 (43,4 %) sigue **por debajo** de la mayoría trivial
(87,6 %): la capa mejora el sistema, no resuelve el problema.

### 6.1 Ablación (escala fija 1,0, fuera de muestra)

| Evidencia | Top-1 | Δ top-1 | Top-3 | Cambios | Aciertos | Errores |
|---|---:|---:|---:|---:|---:|---:|
| solo `operon_position` | 37,87 % | +8,18 pp | 84,50 % | 559 | 261 | 79 |
| solo `leader` | 33,60 % | +3,91 pp | 83,87 % | 327 | 134 | 47 |
| solo `operon_size` | 31,85 % | +2,16 pp | 84,50 % | 370 | 120 | 72 |
| todas (1,0) | **43,40 %** | **+13,70 pp** | **86,88 %** | 734 | 396 | 91 |
| sin `operon_position` | 36,57 % | +6,87 pp | 86,16 % | 497 | 225 | 72 |
| sin `leader` | 39,76 % | +10,06 pp | 86,21 % | 627 | 302 | 78 |
| sin `operon_size` | 41,73 % | +12,04 pp | 84,91 % | 671 | 357 | 89 |

Las tres evidencias aportan; `operon_position` es la principal. El cambio
neto es +305 TSS correctos (396 aciertos − 91 errores) sobre 2.226.

### 6.2 Barrido de escala

| Escala | Top-1 | Δ top-1 | IC 95 % Δ | Top-3 |
|---:|---:|---:|---:|---:|
| 0,00 | 29,69 % | 0,00 pp | [0,00; 0,00] | 80,68 % |
| 0,25 | 33,33 % | +3,64 pp | [+2,65; +4,72] | 84,01 % |
| 0,50 | 37,20 % | +7,50 pp | [+6,11; +9,03] | 85,71 % |
| 0,75 | 40,57 % | +10,87 pp | [+9,21; +12,67] | 86,66 % |
| 1,00 | **43,40 %** | **+13,70 pp** | **[+11,86; +15,68]** | **86,88 %** |

La ganancia es monótona y todos los IC se separan de 0. No hay un punto
óptimo ajustado: el selector de entrenamiento elige 1,0, el extremo
declarado de la rejilla.

### 6.3 Calibración

| Capa | ECE | Brier | NLL | Confianza media |
|---|---:|---:|---:|---:|
| Motivo | 0,0548 | 0,1271 | 1,4800 | 0,347 |
| Contexto | 0,0720 | **0,1175** | **1,3487** | 0,422 |

El Brier y el NLL mejoran; el ECE empeora (0,055 → 0,072). Es decir, las
probabilidades del contexto son mejores en conjunto pero **peor calibradas
en confianza top-1**: el modelo se vuelve más confiado y acierta más, pero
la confianza media (0,422) aún no alcanza el acierto (0,434); el ECE sube
por la forma de la curva, no porque el contexto sea inútil. Se reporta sin
maquillar.

Modo confianza-only (top-1 fijo del motivo):

| Modo | ECE | Brier | NLL | Confianza media | Acierto |
|---|---:|---:|---:|---:|---:|
| Motivo | 0,0533 | 0,2115 | 0,6148 | 0,274 | 29,69 % |
| Confianza con contexto | 0,0732 | **0,1962** | **0,5769** | 0,322 | 29,69 % |

Mismo patrón: Brier/NLL mejores, ECE peor.

### 6.4 Por sigma (top-1, línea base → contexto)

| Sigma | Support | Precisión | Recall | F1 |
|---|---:|---:|---:|---:|
| σ70 | 1.949 | 0,931 → 0,925 | 0,299 → **0,448** | 0,452 → **0,604** |
| σ38 | 182 | 0,173 → 0,193 | 0,308 → 0,368 | 0,221 → 0,253 |
| σ28 | 71 | 0,018 → 0,033 | 0,127 → 0,183 | 0,032 → 0,055 |
| σ32 | 13 | 0,014 → 0,022 | 0,308 → 0,385 | 0,027 → 0,041 |
| σ54 | 11 | 0,020 → 0,026 | 0,909 → 0,727 | 0,040 → 0,051 |

El contexto sube sobre todo el recall de σ70 (0,299 → 0,448), que es el
grueso del gold; a cambio baja el recall top-1 de σ54 (0,909 → 0,727)
porque deja de predecir σ54 en sitios que no lo son. La precisión σ54 sigue
siendo muy baja (2,6 %).

### 6.5 Multietiqueta y σ54

- B. subtilis no tiene TSS multietiqueta en el gold (2.226 single, 0 multi).
  El top-3 es la única vía para cubrir procesos superpuestos.
- Clases: 661 `correct`, 1.135 `partial_top3`, 430 `missed` (línea base) →
  966 `correct`, 968 `partial_top3`, 292 `missed` (contexto).
- σ54: 11 verdaderos, 303 predichos top-1. Precisión 2,6 %; recall top-1
  72,7 %; recall top-3 100 %. Geometría en verdaderos: 10 canónicos, 1
  desplazado; en predichos: 196 desplazados, 107 canónicos.
- La PWM σ54 de B. subtilis es de **una caja** y solo 11 sitios. La
  geometría se reporta como diagnóstico y **no** se usa como prior (regla
  ≥ 20 sitios). Aun así, la precisión σ54 es muy baja: el problema es la
  PWM, no la geometría.

### 6.6 Sensibilidad

| Configuración | Top-1 | Top-3 |
|---|---:|---:|
| PWM propia, 200 señuelos, −80/+20 | 29,69 % | 80,68 % |
| PWM propia, 2.000 señuelos, −80/+20 | 26,28 % | 80,05 % |
| PWM propia, 200 señuelos, −40/+5 | 33,11 % | 83,96 % |
| PWM propia, 2.000 señuelos, −40/+5 | 30,82 % | 83,87 % |
| + contexto, −80/+20 | **43,40 %** | **86,88 %** |
| + contexto, −40/+5 | 46,59 % | 87,83 % |
| + contexto con tabla curada DBTBS (`--curated-context`) | 46,63 % | 88,23 % |

La ganancia del contexto se mantiene en las tres variantes. La tabla curada
da más (+16,9 pp) porque `subtilis_tss_tts.tsv` comparte procedencia con el
gold; el número principal usa el mapeo genómico uniforme (`genes.tsv`) para
no inflar la mejora.

El cambio de 200 a 2.000 señuelos **baja** el top-1 de la PWM propia
(29,7 → 26,3 %): la señal de una caja es débil e inestable al tamaño del
nulo. Es una limitación de la PWM, no de la capa de contexto. En K-12 el
cambio es estable (56,98 → 57,07 %).

---

## 7. Auditoría: de dónde sale el 44,79 %

`docs/PWM_SUBTILIS.md` reportaba:

```text
top-1: 44,79 %
top-3: 85,53 %
mayoría trivial: 87,56 %
```

con el comando:

```bash
python scripts/sigma/50_eval_layers.py \
  --genome subtilis --source gold \
  --pwm data/dbtbs/promoter_pwms_subtilis_all.tsv
```

Pero `50_eval_layers.py` hacía `PWMS_PATH = None` y **nunca asignaba**
`args.pwm`; `load_pwm_library` seguía usando
`data/regulondb/gt/promoter_pwms.tsv` (E. coli). El 44,79 % era, por tanto,
**transferencia de E. coli**, no la PWM propia. `SIGMA_LAYERED_EVAL.md` ya
lo describía correctamente como transferencia entre especies.

Corrección aplicada: `global PWMS_PATH; PWMS_PATH = args.pwm` en `main`, con
test en `tests/test_sigma_context.py::test_eval_layers_honours_pwm_flag`.

Números tras la corrección:

| Configuración | Top-1 | Top-3 | Fuente |
|---|---:|---:|---|
| E. coli PWM sobre B. subtilis (bug `--pwm`) | 44,79 % | 85,53 % | `data/sigma/layers/subtilis` |
| PWM propia de B. subtilis (200 señuelos, −80/+20) | 29,69 % | 80,68 % | esta capa |
| PWM propia + contexto | **43,40 %** | **86,88 %** | esta capa |

**Conclusión de auditoría:** el 44,79 % no era un resultado de la PWM
específica de B. subtilis. La PWM propia es más débil (29,69 %), y el
contexto la lleva a 43,40 %, todavía por debajo de la mayoría trivial. No se
debe presentar el 44,79 % como “lo que da la PWM propia” ni como una
comparación válida contra el 87,56 %.

---

## 8. Límites

1. **No hay validación experimental.** El gold es RegulonDB (K-12) y
   DBTBS/BSGatlas (B. subtilis). No hay negativos: un TSS fuera del gold no
   es automáticamente un falso positivo.
2. **La validación es cruzada dentro de cada especie**, no una validación
   externa en un genoma nuevo. La mejora de B. subtilis es válida para la
   población y la anotación del gold, no para cualquier cepa.
3. **Procedencia compartida.** En B. subtilis, el gold y `genes.tsv` /
   `gt_operons.tsv` vienen de BSGatlas/DBTBS. La asociación
   “TSS interno → σ38” puede reflejar cómo se anotó el gold, no una regla
   universal. La tabla curada sube la ganancia a +16,9 pp, lo que sugiere
   parte de dependencia de la anotación.
4. **Cobertura de contexto desigual.** `pathway_ec` 27,5 % (solo K-12),
   `conservation` 23,4 % (K-12) / 4,5 % (B. subtilis, no usada), `domain`
   0 %, `expression` 0 %. No se inventan datos.
5. **`fam` no es dominio.** Se excluyó por alta cardinalidad (4212 valores
   para 4262 genes en K-12). `structure_features.tsv` cubre genomas GCF,
   no estos dos.
6. **σ54.** K-12 tiene geometría de dos cajas y 99 sitios; B. subtilis tiene
   una caja y 11 sitios. La precisión σ54 es baja (18,2 % / 2,6 %) y el
   recall top-3 es 63,6 % / 100 %. La capa de contexto no arregla la PWM σ54.
7. **Inestabilidad al nulo en B. subtilis.** 200 → 2.000 señuelos baja el
   top-1 de 29,7 % a 26,3 %. La señal de una caja es débil; hace falta más
   datos o una PWM de dos cajas para concluir.
8. **La calibración no siempre mejora.** En B. subtilis el ECE empeora
   (0,055 → 0,072) aunque Brier y NLL mejoran. En K-12 el ECE de la vía
   principal no cambia; el modo confianza-only sí mejora (0,061 → 0,053).
9. **Desempates.** La línea base K-12 de esta corrida es 56,98 %; la de
   `50_eval_layers.py` es 57,04 %. La diferencia (2 TSS) es de desempate,
   no de datos.
10. **No se toca el modelo ni los golds.** Todo es lectura + tablas
    derivadas en `data/sigma/context/`.

---

## 9. Reproducir

```bash
# autopruebas
python scripts/sigma/sigma_context.py --selftest
python scripts/sigma/51_eval_context.py --selftest
python -m pytest -q tests/test_sigma_context.py

# evaluación oficial (PWM propia por especie; K-12 y B. subtilis)
python scripts/sigma/51_eval_context.py --all

# sensibilidades de B. subtilis
python scripts/sigma/51_eval_context.py --genome subtilis --up 40 --down 5
python scripts/sigma/51_eval_context.py --genome subtilis --curated-context
python scripts/sigma/51_eval_context.py --genome subtilis --decoys 2000

# línea base de capas (con la corrección de --pwm)
python scripts/sigma/50_eval_layers.py \
  --genome subtilis --source gold \
  --pwm data/dbtbs/promoter_pwms_subtilis_all.tsv
```

Salidas:

```text
data/sigma/context/<genoma>/
  context_metrics.json            métricas, ablación, calibración, σ54
  context_sites.tsv               una fila por TSS
  context_predictions.tsv         una fila por TSS × sigma: z, LLR crudo,
                                  LLR ponderado, escala y score final
  context_evidence.tsv            tablas de LLR (ajuste completo, descriptivo)
  context_evidence_coverage.tsv   fuente, cobertura y peso por evidencia
  context_report.md               informe legible
```

---

## 10. Frase final

**La capa de contexto no rescata K-12: allí no hay señal contextual
out-of-fold y se apaga sola; su único aporte medible es calibrar la
confianza (ECE 0,061 → 0,053). En B. subtilis, con PWM propia y sin
transferir E. coli, sube el top-1 de 29,69 % a 43,40 % (+13,70 pp, IC 95 %
[11,86; 15,68]) y el top-3 de 80,68 % a 86,88 %, pero sigue por debajo de la
mayoría trivial del 87,56 %. El 44,79 % documentado era una corrida con
PWMs de E. coli: no se debe usar como resultado de la PWM propia.**
