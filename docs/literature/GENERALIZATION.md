# Generalización entre especies en modelos de regulación

> **Fecha:** 2026-11-26 · **Ámbito:** capa `sigma` de `dna_poc` (catálogo de ~400
> especies, 7 genomas con gold).
> **Pregunta:** ¿qué métodos permiten que un modelo de regulación entrenado en
> una especie funcione en otra, bajo qué supuestos, y cómo se valida sin
> engañarse?
> **Respuesta corta:** la generalización no se mide con exactitud top-1 sobre el
> gold disponible. Se mide con un diseño de partición que respete la filogenia
> (LOSO/LOCO), baselines calculados **por fold**, métricas por clase,
> calibración y abstención. En este repositorio el top-1 de K-12 es 57,0 %
> frente a 62,6 % de la mayoría trivial, y en *B. subtilis* la mejor política
> (híbrida) llega a 87,6 % —exactamente la mayoría trivial— mientras su PWM
> propia se queda en 29,7 %: el número que se cita como "precisión" mide sobre
> todo el prior de clase.

Este documento es una revisión metodológica con bibliografía verificada
(DOI/PMID resueltos el 2026-11-26, ver §14). No propone cambios de código; las
rutas y scripts que menciona son los existentes, y lo que no existe se marca
como **propuesto**.

---

## 0. Resumen ejecutivo

1. **Interpolar no es generalizar.** Un `GroupKFold` por bloques de 10–100 kb
   dentro del mismo genoma (lo que usa hoy `SIGMA_STRATEGIES.md`) mide
   interpolación espacial. Una partición por especie mide otra cosa. Una
   partición por filo mide una tercera. Las tres son legítimas, pero autorizan
   afirmaciones distintas (§1.2).
2. **El top-1 no distingue señal de prior.** Con 87,6 % de σ70 en *B. subtilis*,
   un modelo constante acierta 87,6 %. En K-12 la mayoría trivial (62,6 %) supera
   al asignador (57,0 %). Cualquier comparación entre especies con priors
   distintos es inválida sin macro-F1, recall por clase y curvas PR (§8).
3. **Transferir el motivo y transferir el umbral son operaciones distintas.** La
   forma del PWM viaja parcialmente entre ortólogos; el umbral, el fondo nulo y
   el prior **no** viajan. El propio repo documenta el caso: umbrales derivados
   de K-12 dan AUC < 0,5 para σ38 en *B. subtilis*
   ([`SIGMA_MOTIF_CALIBRATION.md`](../SIGMA_MOTIF_CALIBRATION.md)).
4. **El gold actual es solo de positivos.** No hay conjunto negativo, y "no está
   en el gold" no significa "no es promotor". La consecuencia práctica es que la
   **precisión no es estimable** en las 400; lo estimable es recall condicionado,
   cobertura y abstención (§6).
5. **La jerarquía es la herramienta natural, pero con 7 especies no se estima
   bien la varianza entre especies.** El *partial pooling* ayuda a estabilizar
   umbrales por especie; no sustituye datos ausentes (Gelman, 2006). Hay que
   reportar la incertidumbre, no solo el efecto (§3).
6. **Los 7 genomas con gold ya permiten un LOCO real de 3 pliegues** (tres filos:
   Pseudomonadota, Bacillota, Actinomycetota). Es la única medición de
   transferencia entre filos disponible hoy, y es descriptiva, no inferencial
   (§5.4, §13.3).
7. **Recomendación:** modelo jerárquico multi-tarea (motivo compartido por
   familia + intercepto por especie + cabeza multi-etiqueta), validado con tres
   diseños reportados por separado (interpolación, LOSO, LOCO), con compuertas de
   decisión numéricas y salida de tres estados (`assign` / `unresolved` /
   `no_model`). Si no pasa las compuertas, se publica inventario + evidencia de
   motivo, sin asignación de sigma (§13).

---

## 1. Generalización frente a interpolación: definición operativa

### 1.1 Los cuatro desplazamientos que importan

En genómica comparada, "otra especie" no es un solo tipo de cambio. Conviene
separarlos porque cada uno se corrige distinto (Pan y Yang, 2010; Ben-David et
al., 2010).

| Tipo de *shift* | Qué cambia entre especies | Ejemplo en `sigma` | Se corrige con |
|---|---|---|---|
| **Covariante** `P(x)` | Composición de secuencia: GC, uso de codones, sesgo de k-meros | Genomas AT-ricos vs GC-ricos; ventanas con distinto contenido de A/T | Fondo nulo específico de especie (`70_local_background.py`), normalización por GC, negativos apareados por composición |
| **Prior / etiqueta** `P(y)` | Frecuencia de uso de cada sigma | σ70 domina 62,6 % en K-12 y 87,6 % en *B. subtilis*; los regulones σ54 son pequeños | Recomputar el baseline por fold; corrección de *label shift* (Lipton et al., 2018); macro-F1 en vez de exactitud |
| **Concepto** `P(y\|x)` | El significado de la etiqueta cambia | σ38 de *E. coli* (RpoS) ≠ SigB de *B. subtilis*: el mapeo B→σ38 es "razonable" pero el motivo no coincide (AUC 0,276). La propia definición de "promotor" y "TSS" en bacterias es objeto de revisión conceptual (Mejía-Almonte et al., 2020) | Validar el mapeo de etiquetas con ortología + motivo; modelar por familia, no por nombre; declarar `no_model` |
| **Filogenético** | Correlación entre muestras por ancestría común | Dos cepas de *E. coli* no son dos muestras independientes | Particiones por clado, efectos aleatorios por especie/clado, intervalos que incorporen la filogenia |

La distinción entre *shift* de prior y de concepto es la que más se ignora. Un
modelo puede transferirse bien en forma de motivo y, aun así, producir métricas
peores porque el prior de la especie destino es distinto. Eso no es un fallo del
modelo: es un fallo de la métrica.

### 1.2 Criterio operativo: ¿qué afirmación autoriza cada partición?

Regla de trabajo: **una partición autoriza la afirmación más débil de la cadena
de dependencia que rompe**. Si el *holdout* comparte linaje cercano con el
entrenamiento, la afirmación es "tolera variación intra-género", no
"generaliza".

| Partición | Qué se deja fuera | Qué se mide | Afirmación permitida | Estado en el repo |
|---|---|---|---|---|
| Aleatoria por TSS | TSS sueltos del mismo genoma | Memorización | Ninguna (fuga por ventanas solapadas y operones vecinos) | No usar |
| `GroupKFold` por bloques de 10–100 kb | Regiones del mismo genoma | Interpolación espacial | "Funciona dentro de la especie" | **Sí**: `SIGMA_STRATEGIES.md` |
| Leave-one-operon-out | Operones completos del mismo genoma | Robustez a composición de operón | Interpolación | No implementado |
| LOSO, mismo género | Una cepa/especie hermana | Variación intra-género | "Tolera variación intra-género" | Posible con los 400 |
| LOSO, distinto género mismo filo | Un género no visto | Variación de composición y arquitectura | "Transfiere dentro del filo" | K-12 ↔ *B. subtilis* (sin reentrenar) |
| **LOCO por filo** | Un filo completo | Repertorio de sigma y gramática de promotor | "Generaliza a un filo no visto (con reservas)" | **Posible hoy**: 3 filos en el gold |
| Aleatorio sobre los 400 genomas | Genomas mezclados de todos los linajes | Prior y composición | **No es generalización** | Error común (§12) |
| Temporal (anotación antigua vs nueva) | Práctica de curación | Cambio de criterio de anotación | "Estable al cambio de versión de la base" | No implementado |

**Prueba de fuga cuantitativa (propuesta, no implementada).** Para cada genoma
*held-out*, calcular la similitud máxima contra los genomas de entrenamiento:
Jaccard de k-meros (5-meros, submuestreo) o ANI si hay genoma completo. Reportar
la distribución. Criterio propuesto: si la mediana del *held-out* supera ~0,5 de
Jaccard con algún genoma de entrenamiento, la partición es **interpolación** y no
debe etiquetarse como generalización. No es un umbral de la literatura; es una
convención del proyecto que hay que fijar **antes** de mirar resultados.

Para la distancia taxonómica conviene usar dos referencias complementarias: el
**ANI** con frontera de especie en 95–96 % (Richter y Rosselló-Móra, 2009) y la
taxonomía estandarizada por filogenia de genomas completos de **GTDB** (Parks et
al., 2018), que además da rangos comparables entre filos. El repo ya tiene
`data/taxonomy/ani_hits.tsv` y `data/taxonomy/ani_groups.json`: la infraestructura
para calcular esta distancia existe.

### 1.3 Por qué esto no es semántica: mecánica de la fuga

- **Fuga por estructura de datos.** Roberts et al. (2017) muestran que la
  validación cruzada aleatoria subestima el error cuando los datos tienen
  estructura temporal, espacial, jerárquica o filogenética, incluso si el modelo
  incluye términos autorregresivos o efectos mixtos. La solución es *block
  cross-validation*: partir por la unidad que induce dependencia (aquí, especie o
  clado).
- **Fuga por selección de modelo.** Varma y Simon (2006) cuantifican el sesgo
  optimista de elegir hiperparámetros y evaluar con las mismas particiones. Con
  ~30 scripts de evaluación en `scripts/sigma/`, este riesgo es real: cada
  variante probada sobre el mismo gold consume grados de libertad.
- **Fuga por práctica de investigación.** Kapoor y Narayanan (2023) documentan
  ocho tipos de fuga que explican fallos de reproducibilidad en ciencia basada en
  ML; Whalen et al. (2022) los aterrizan en genómica: estructura de datos que
  sesga la estimación de desempeño y baselines ausentes o mal calculados.

**Ejemplo del propio repo (y es el más instructivo).** Durante un tiempo
`50_eval_layers.py` parseaba `--pwm` pero lo ignoraba, de modo que la evaluación
"de *B. subtilis*" corría con la matriz de *E. coli* y reportaba 44,8 % de top-1.
Con la matriz propia el número real es **29,7 %**
([`SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md) §7). Dos lecciones:

1. Un flag que se ignora silenciosamente produce un resultado que *parece*
   transferencia y no lo es. La procedencia de la matriz debe ir en la salida
   (`model_source`), no en la línea de comandos.
2. La "transferencia" de *E. coli* **puntuó más alto** que la matriz nativa
   (44,8 % vs 29,7 %). La hipótesis más plausible es que la matriz ajena es menos
   específica y colapsa predicciones hacia σ70, que es la clase mayoritaria; es
   decir, ganó prior, no información. Es comprobable: descomponer la ganancia en
   cambios de prior vs cambios de verosimilitud (matriz de confusión por clase,
   calibración). **No se ha hecho todavía**; queda como primer experimento de
   auditoría (§11, fase B).

---

## 2. Transfer learning y domain adaptation en genómica

### 2.1 Taxonomía de métodos y qué exige cada uno

| Familia | Mecanismo | ¿Necesita etiquetas en destino? | Ejemplo genómico | Aplicabilidad a `sigma` |
|---|---|---|---|---|
| **Basada en features** | Representación compartida (ortología, embeddings, motivos) | No | SATL: integración entre especies **sin** ortología (Park et al., 2024); CNN de promotores procariotas (Umarov y Solovyev, 2017) | Alta: ya existe `data/sigma/sigma_orthologs.tsv` (2.432 grupos) |
| **Instancias / reweighting** | Reponderar muestras fuente por similitud con destino | No | Corrección de *label shift* con predictor caja negra (Lipton et al., 2018) | Media: requiere estimar el prior de sigma en destino, que es justo lo que falta |
| **Parámetros (preentrenar/afinar)** | Entrenar en fuente, afinar en destino | Sí (pocas) | Nucleotide Transformer, 850 genomas de 3.202 humanos (Dalla-Torre et al., 2025); Evo, 2,7 M genomas procariotas (Nguyen et al., 2024) | Media-baja: los genomas del catálogo son pocos para afinar por especie |
| **Alineación de momentos** | Igualar primer y segundo momento de la distribución entre especies | **No** | MORALE: fuente(s) con etiquetas → destino sin etiquetas; sin parámetros extra (Ebeid et al., 2026) | **Alta**: los 400 genomas sin gold aportan momentos no etiquetados |
| **Adversarial (GRL)** | Discriminador especie vs tarea | No | Baseline adversarial en MORALE; degradó frente a alineación de momentos en ese estudio | Baja: inestable con 2–5 dominios |
| **Aumento de features** | Concatenar features globales y por dominio ("frustratingly easy") | No | Nombre tomado del enfoque clásico; en genómica, MORALE es la instancia publicada | Media |
| **Multi-tarea** | Compartir tronco entre tareas/especies | Sí | PromoterAtlas (3.371 especies gammaproteobacterianas), DeepCROSS (*E. coli* + *P. aeruginosa*) | Alta (§4) |

### 2.2 Evidencia disponible (y qué mide realmente)

- **MORALE (Ebeid et al., 2026)** es el caso más cercano al problema de `sigma`:
  entrenar con una o varias especies fuente y predecir sitios de unión de TF en
  una especie destino **sin datos de unión**. Alinea momentos de la distribución
  y no añade parámetros. En el caso de cinco especies mamíferas mejoró la
  predicción en humano por encima de un modelo solo-humano. Es un resultado de
  un *preprint*/capítulo de congreso, con dos conjuntos de ChIP-seq; no es
  evidencia de que funcione en bacterias. Aun así, el principio —usar genomas no
  etiquetados para estimar la distribución destino— es el único que se puede
  aplicar hoy a las 400.
- **PromoterAtlas (Coppens y Ledesma-Amaro, 2026)**: transformer de 1,8 M de
  parámetros entrenado con 9 M de secuencias regulatorias de **3.371 especies**
  gammaproteobacterianas; reconoce RBS, promotores, sitios de TF y terminadores,
  y produce anotación genómica completa. Es la demostración a mayor escala de que
  un modelo multi-especie de regulación bacteriana es viable. Límite importante:
  cubre **una clase** (Gammaproteobacteria, dentro de Pseudomonadota); no dice
  nada sobre transferir a Bacillota o Actinomycetota.
- **DeepCROSS (Zhang et al., 2025)**: representación de 1,8 M de secuencias
  regulatorias de miles de genomas + optimización multi-tarea para diseñar
  secuencias en *E. coli* y *P. aeruginosa*. Aporta la idea de "meta-representación"
  como espacio común entre especies.
- **Promotores específicos de sigma con validación experimental (Van Brempt et
  al., 2020)**: bibliotecas de promotores σ70 de *E. coli* y σB/σF/σW de
  *B. subtilis*, con modelos que predicen frecuencia de inicio de transcripción y
  **ortogonalidad** entre sigmas. Es el precedente más directo de "modelo por
  sigma con validación cruzada de especie"; su límite es que son bibliotecas
  sintéticas, no genomas.
- **Elementos regulatorios selectivos de especie (Johns et al., 2018)**: minería
  metagenómica para obtener elementos con expresión selectiva de hospedero.
  Relevante como recordatorio de que "funciona en otra especie" puede ser un
  objetivo de diseño y no una propiedad del modelo.
- **Geneformer de ratón (Ito et al., 2025)**: modelo preentrenado en humano
  aplicado a ratón, con la observación de que los procesos conservados
  (housekeeping) transfieren mejor que los linajes-específicos. Es exactamente el
  patrón esperable en sigma: σ70 (housekeeping) transferirá mejor que σB/σF
  (respuesta específica).
- **Modelos de secuencia a función** (Enformer, Avsec et al., 2021; Basenji,
  Kelley et al., 2018; DeepSEA, Zhou y Troyanskaya, 2015; Sei, Chen et al., 2022;
  DNABERT, Ji et al., 2021) son el estado del arte eucariota. Se citan aquí por
  su metodología de evaluación (tareas múltiples, *held-out* por cromosoma), no
  porque sus arquitecturas transfieran a procariotas, donde la ventana útil es de
  decenas de pb, no de 100 kb.

### 2.3 Qué transfiere y qué no en `sigma`

| Componente | ¿Transfiere entre especies? | Evidencia / razón | Acción |
|---|---|---|---|
| Forma del motivo (PWM) | Parcialmente, dentro de familia | Feklistov et al. (2014): los dominios de reconocimiento de −10/−35 están conservados dentro de la familia σ70; σ54 usa −24/−12 y es otra máquina | PWM por familia ortóloga, no por nombre global |
| Ancho de caja y espaciador | No asumir | σ54: espaciador 10–11; σ70: 15–19; el repo ya separa arquitecturas (`SIGMA_AUDIT.md` §4.2) | Arquitectura como covariable discreta |
| Fondo nulo | **No** | Composición GC distinta por genoma | `70_local_background.py` por especie |
| Umbral de score | **No** | Umbrales de K-12 → *B. subtilis*: AUC σ38 = 0,276 (< 0,5) | Umbral por especie con *prior* jerárquico o abstención |
| Prior de clase | **No** | 62,6 % vs 87,6 % de σ70 | Recomputar por fold; nunca usar exactitud cruda |
| Mapeo de etiquetas | Parcialmente | σ38↔SigB: mapeo funcional razonable, motivo distinto; el inventario deja 62,8 % en `sigma_unresolved` | Validar mapeo con ortología + motivo; declarar ambigüedad |
| Inventario de genes | Sí, con ruido | Dominios Pfam conservados; pero 62,8 % sin familia y `rseD` mal etiquetado por nombre | Familia amplia (`sigma70_like`, `ecf_like`, `sigma54_like`) como unidad |

**Supuesto central que hay que hacer explícito:** la gramática del promotor
(qué posiciones importan y con qué peso) es aproximadamente compartida dentro de
una familia de sigma, y las diferencias entre especies son un desplazamiento de
fondo y de umbral, no un cambio de gramática. Este supuesto es **falso en los
bordes** (promotores con −10 extendida, UP elements, dependencia de bEBP para
σ54) y debe testarse por familia, no asumirse globalmente.

---

## 3. Modelos jerárquicos y filogenéticos

### 3.1 *Partial pooling*: la estructura mínima correcta

En vez de un umbral por especie (sin pooling: ruidoso con 11 sitios) o un umbral
global (pooling completo: ignora la especie), se estima un efecto por especie que
se encoge hacia la media de la familia:

```text
logit P(σ | x, s) = α_σ + u_{σ,s} + f_σ(x)
u_{σ,s} ~ Normal(0, τ²_σ)
```

- `α_σ`: intercepto global de la clase σ (el prior de familia).
- `u_{σ,s}`: desviación de la especie `s` respecto de la familia; se encoge según
  `τ_σ`.
- `f_σ(x)`: la señal de secuencia (log-odds del PWM o del modelo).

Si `τ_σ → 0`, el modelo decide que las especies son iguales (pooling completo).
Si `τ_σ → ∞`, cada especie va por su lado. **El dato decide**, y ese es el punto:
no hay que elegir a mano entre "un modelo global" y "un modelo por especie".

[`SIGMA_STRATEGIES.md`](../SIGMA_STRATEGIES.md) reporta un `jerárquico logístico`
como la mejor estrategia en K-12 (top-1 70,6 %, macro-F1 0,510). Dos
observaciones antes de reutilizarlo: (a) su definición no está en el documento ni
aparece implementada en `scripts/sigma/` (búsqueda de "jerárquico/hierarchical"
sin resultados), así que hay que recuperarla y documentarla; (b) por el contexto
del benchmark, lo más probable es que sea jerárquico **dentro** de una especie
(decidir primero σ70 vs alternativa), no entre especies. La generalización exige
añadir el nivel especie a la jerarquía, no reutilizar esta tal cual.

**Advertencia de Gelman (2006):** el modelado multinivel no sustituye información
ausente. Con 7 especies no se estima con precisión `τ_σ`; se obtiene una
estimación imprecisa que hay que reportar como tal (posterior o intervalo), no
como un número puntual. Harrison et al. (2018) ofrecen la guía práctica de
especificación y selección de modelos mixtos.

### 3.2 Encogimiento de PWM (empirical Bayes)

Problema concreto del repo: σ54 tiene 96–99 sitios en K-12 y **11** en
*B. subtilis* ([`SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md) §2). Estimar una PWM de
−24/−12 con 11 sitios produce sobreajuste.

Solución estándar: prior de Dirichlet hacia la PWM de familia (Sjölander et al.,
1996):

```text
p_{s,j} = (c_{s,j} + κ · m_fam,j) / (n_s + κ)
```

donde `c_{s,j}` son los conteos de la especie en la posición `j`, `m_fam,j` la
frecuencia de la familia y `κ` la fuerza del encogimiento. Con `n_s = 11`, la PWM
resultante es casi la de familia; con `n_s = 99`, domina la evidencia local. Esto
es preferible a "usar la PWM de *E. coli*" (lo que se hacía por accidente) y a
"usar la PWM propia con 11 sitios" (lo que sobreajusta).

### 3.3 Estructura filogenética

Tres niveles de uso, de menor a mayor compromiso:

1. **Particiones por clado** (obligatorio). Es el uso mínimo: agrupar especies por
   filo/clase antes de partir. No requiere modelo.
2. **Efectos aleatorios anidados** (recomendado). `u_{σ, filo} → u_{σ, especie}`
   con pooling parcial en cada nivel. Captura que especies del mismo filo se
   parecen más.
3. **Modelo filogenético explícito** (opcional). PGLS o regresión logística
   filogenética (Pagel, 1999; Ives y Garland, 2010) usando la matriz de
   covarianza inducida por el árbol. Requiere árbol y un modelo de evolución
   (Brownian motion por defecto), ambos discutibles en bacterias por transferencia
   horizontal.

Sobre el nivel 3 hay dos resultados que conviene tener presentes:

- **La señal filogenética se mide, no se asume.** Münkemüller et al. (2012)
  revisan las medidas (Blomberg's K, Pagel's λ) y sus nulos. Si λ ≈ 0 para una
  clase de sigma, imponer filogenia solo añade ruido.
- **Predecir en un taxón no muestreado no es evaluar la ecuación de regresión.**
  Gardner et al. (2025) muestran que los intervalos de predicción
  filogenéticamente informados mejoran 2–3× frente a usar la regresión OLS/PGLS
  como si el taxón nuevo fuera una observación más, y que el intervalo debe
  crecer con la longitud de rama terminal. Traducido: la incertidumbre de asignar
  sigma en un filo sin gold **debe** ser mayor que en una especie hermana de una
  especie con gold, y el reporte debe reflejarlo.

### 3.4 Qué compra y qué no compra la jerarquía

| Compra | No compra |
|---|---|
| Estabiliza umbrales con pocos sitios por especie | No crea señal donde el motivo no está |
| Permite préstamo de información entre especies del mismo clado | No corrige un mapeo de etiquetas equivocado |
| Hace explícita la incertidumbre por especie | No sustituye un conjunto negativo |
| Separa prior de clase de señal de secuencia | No justifica afirmar generalización sin LOSO/LOCO |

---

## 4. Multi-tarea y multi-especie

### 4.1 Qué aporta el multi-tarea (y cuándo falla)

Caruana (1997) establece el resultado básico: entrenar tareas relacionadas con
representación compartida mejora la generalización porque actúa como
regularizador y porque las tareas aportan información mutua. Falla cuando las
tareas compiten (capacidad limitada) o cuando la relación es espuria: es el
*negative transfer* que Pan y Yang (2010) documentan.

En `sigma` hay tres ejes de tarea que no compiten necesariamente:

| Eje | Tarea compartida | Por qué ayuda |
|---|---|---|
| Familia de sigma | Clasificación 6 clases (σ70, σ24, σ28, σ32, σ38, σ54) | Regulariza la representación de motivo |
| Especie | Intercepto/embedding por especie | Absorbe prior y composición sin tocar el motivo |
| Multi-etiqueta | ¿Qué sigmas pueden usar este TSS? | Refleja la biología (124 TSS con >1 sigma en K-12) |
| Auxiliares | RBS, terminador, operón, arquitectura | Añaden señal estructural y evitan que el modelo use composición como atajo |

Un **baseline fuerte y honesto** para cualquier modelo multi-tarea es ERM (riesgo
empírico minimizado) con validación en dominio de entrenamiento. Gulrajani y
Lopez-Paz (2021) mostraron que, con selección de modelo justa, ningún algoritmo
de *domain generalization* superaba consistentemente a ERM en su benchmark. La
consecuencia para este proyecto: **no adoptar un método complejo sin demostrar
que supera a PWM + umbral por especie en LOSO**. El método simple es el rival a
vencer, no el punto de partida.

### 4.2 Multi-especie en bacterias: qué existe

| Trabajo | Escala | Qué demuestra | Qué no demuestra |
|---|---|---|---|
| PromoterAtlas (2026) | 3.371 especies, 9 M secuencias, 1 clase | Un transformer multi-especie anota promotores y elementos | Transferencia fuera de Gammaproteobacteria |
| DeepCROSS (2025) | Miles de genomas, 2 especies objetivo | Representación común + optimización multi-tarea para diseño | Predicción de sigma nativo |
| Van Brempt et al. (2020) | *E. coli* + *B. subtilis*, 4 sigmas | Modelo por sigma con ortogonalidad validada experimentalmente | Escala genómica |
| MORALE (2026) | 2 y 5 especies mamíferas | Alineación de momentos sin etiquetas destino | Bacterias; sigma |
| SATL (2024) | Transcriptómica multi-especie | Transferencia sin ortología | Regulación a nivel de sitio |

No existe, hasta donde alcanza esta revisión, un modelo multi-especie de
**asignación de sigma a promotores** validado en varios filos. Eso no es un
inconveniente para el proyecto: es la contribución posible, siempre que la
validación sea la correcta.

### 4.3 El espacio de etiquetas: familia, ortólogo o nombre

Este es el punto donde el multi-especie se rompe en la práctica. El repo ya
documenta el problema:

- 5.993 filas de inventario, **62,8 % en `sigma_unresolved`**
  ([`SIGMA_AUDIT.md`](../SIGMA_AUDIT.md)).
- Etiquetas por nombre sin modelo asociado: `sigmaE`, `sigmaB`, `sigmaH`,
  `sigmaK`, `sigmaF`.
- Un falso positivo por regla de texto: `rseD` ("rpoE leader peptide") clasificado
  como `sigmaE`.
- En `data/sigma/sigma_orthologs.tsv`, el grupo **OG0001** (720 miembros, 381
  genomas) contiene a la vez `sigma70` (594), `sigma38` (18), `sigmaB` (23) y
  `sigma_unresolved` (85). Parte de esa mezcla es biología (la familia σ70 es
  amplia); parte es ruido de anotación. **No se puede distinguir sin auditarlo.**

Consecuencia de diseño: el objetivo de clasificación debe definirse en **dos
niveles** y no mezclarlos:

1. **Nivel familia** (`sigma70_like`, `ecf_like`, `sigma54_like`): es lo que
   transfiere entre filos, porque la maquinaria y la arquitectura se conservan.
2. **Nivel ortólogo/nombre** (RpoD, RpoS, SigB, SigF…): es lo que el usuario
   quiere, pero es lo que **no** transfiere de forma fiable entre filos.

La salida debe declarar en qué nivel está la evidencia. Un `sigma70_like` con
alta confianza y un `rpoS` con baja confianza son compatibles; un `rpoS` con alta
confianza en un genoma sin gold es una afirmación que hoy no está respaldada.

### 4.4 Multi-etiqueta: la biología no es top-1

Un TSS puede ser reconocido por más de un sigma in vitro y in vivo. En K-12, 124
TSS del gold tienen más de un sigma; el benchmark lo reporta como
`hit_one_of_multi` (82) y `top3_all`
([`SIGMA_LAYERED_EVAL.md`](../SIGMA_LAYERED_EVAL.md) §3.1). Un clasificador
top-1 **no puede** acertar todos por construcción. Métricas correctas:

- `hit_any` (¿al menos uno?), `hit_all` (¿todos?), Jaccard del conjunto predicho
  frente al verdadero;
- F1 micro y macro sobre conjuntos, no sobre etiquetas individuales;
- recall por clase calculado sobre TSS mono-etiqueta, para no mezclar regímenes.

---

## 5. Validación: LOSO, LOCO y sus parientes

### 5.1 Diseños y qué estima cada uno

| Diseño | Unidad de *holdout* | Estima | Varianza de la estimación |
|---|---|---|---|
| *k*-fold aleatorio | TSS | Sesgo optimista | Subestima |
| `GroupKFold` (bloques genómicos) | Región | Interpolación espacial | Apropiada dentro de especie |
| LOSO | Especie | Transferencia a especie no vista | 1 estimación por especie (7) |
| **LOCO** | Clado (filo/clase) | Transferencia a linaje no visto | 1 estimación por clado (3) |
| LOPO | Filo | Generalización profunda | Muy pocas unidades |
| Anidado | Selección + evaluación | Desempeño honesto | Requiere disciplina |

### 5.2 Protocolo mínimo (cinco reglas)

1. **Congelar las particiones antes de mirar resultados.** Un manifiesto
   `splits.tsv` con `genome → fold` versionado. Si cambia, se invalida el reporte.
2. **Selección dentro del entrenamiento.** Todo hiperparámetro (ventana, umbral,
   `κ`, `τ`) se elige con LOSO interno sobre las especies de entrenamiento.
   Varma y Simon (2006): elegir y evaluar con las mismas particiones infla el
   resultado.
3. **Baselines por fold, no globales.** La mayoría trivial debe recalcularse en
   cada fold (en *B. subtilis* es 87,6 %, en K-12 62,6 %). Un baseline global
   hace que un modelo mediocre parezca bueno en la especie con prior extremo.
4. **Reportar por fold, no solo agregado.** El agregado oculta que un modelo gana
   en 2 especies y pierde en 5. Con 7 folds, la mediana y el rango son más
   informativos que la media.
5. **Incertidumbre por *cluster*.** Los TSS de una especie no son independientes.
   El intervalo se calcula remuestreando **especies**, no TSS (Field y Welsh,
   2007), o con la corrección de Nadeau y Bengio (2003) para comparar modelos
   entre folds. En genómica microbiana esto es especialmente crítico: la
   estructura poblacional hace que las particiones aleatorias confluyan linaje con
   determinante, y las revisiones de ML aplicado a genomas bacterianos lo señalan
   como limitación práctica recurrente (Kim et al., 2022).

### 5.3 Comparar modelos entre folds

Con 7 folds, una diferencia de 2 pp puede ser ruido. Recomendaciones:

- test *t* corregido de Nadeau-Bengio (2003) para diferencias pareadas entre
  modelos;
- McNemar dentro de especie para comparaciones binarias (¿acierta A donde B
  falla?);
- reportar el número de folds en los que A > B, no solo la media.

### 5.4 Qué puede y qué no puede probar el LOSO de este proyecto

**Puede probar:** que un modelo entrenado sin una especie mantiene recall por
clase y calibración aceptables en esa especie, y cuánto cae respecto de la
evaluación intra-especie. Es la medición directa de "¿esto sirve para una especie
nueva?".

**No puede probar:**

- **Precisión en las 400.** El gold no tiene negativos; el LOSO con 7 especies no
  informa sobre la tasa de falsos positivos en un genoma sin curación.
- **Generalización a filos ausentes del gold.** Con 3 filos, LOCO da 3 números.
  Es descriptivo. No hay potencia para afirmar "funciona en cualquier filo".
- **Mecanismo.** Que el modelo acierte no implica que use el motivo correcto;
  para eso hace falta atribución (importancia por posición) y contraste con la
  PWM conocida. Sin eso, un modelo puede ganar por composición de GC.

**El número que sí es interpretable para las 400** no es exactitud: es
**cobertura × abstención × calibración**, más la comparación con el baseline
trivial dentro de cada clado.

---

## 6. Etiquetas: ruido, gold incompleto y negativos

### 6.1 Fuentes de ruido en `sigma` (todas verificables en el repo)

| Fuente | Mecanismo | Ejemplo local | Consecuencia |
|---|---|---|---|
| Anotación por nombre | Coincidencia de texto | `rseD` → `sigmaE` | Falso positivo sistemático |
| Anotación por dominio | Familias amplias | 62,8 % `sigma_unresolved` | Clase "basura" que el modelo aprende |
| Propagación entre especies | Transferencia de anotación | OG0001 con 4 etiquetas distintas | Ruido correlacionado con el clado |
| Curación con criterios distintos | RegulonDB vs DBTBS vs SubtiWiki | `confidence` S/W/C vs `curated` | El "gold" no es homogéneo |
| Coordenada de TSS | 1-based, hebra, TSS alternativos | gold de 3.314 filas vs 4.057 promotores curados | Desalineación de denominadores |
| Multi-etiqueta | Biología real | 124 TSS con >1 sigma | "Error" que no lo es |
| Ausencia de negativos | Solo positivos curados | `SIGMA_GOLD.md` §2: "No es un conjunto negativo" | Precisión no estimable |

Schnoes et al. (2009) documentan que el error de anotación en bases públicas es
sustancial y **se propaga** por transferencia; no es un problema hipotético.
Frénay y Verleysen (2014) revisan el efecto del ruido de etiqueta: no es
monótono —a tasas bajas puede regularizar, a tasas altas domina el sesgo— y
depende de la geometría de clases. La consecuencia práctica: **no se puede
asumir que el ruido es pequeño ni que es aleatorio**; hay que medirlo.

### 6.2 Aprendizaje con positivos y no etiquetados (PU)

El gold de `sigma` es exactamente un conjunto PU: hay TSS con sigma curado y
millones de posiciones sin etiquetar. El resultado clásico (Elkan y Noto, 2008)
es que se puede aprender un clasificador consistente desde positivos y no
etiquetados bajo el supuesto **SCAR** (*selected completely at random*): la
probabilidad de que un positivo esté etiquetado no depende de `x`.

Ese supuesto **falla** aquí, y hay que decirlo: la probabilidad de curación
depende de la especie (K-12 tiene 3.314 filas; *P. aeruginosa* 18), del tipo de
sigma (σ70 está sobrerrepresentado), del tipo de promotor (los canónicos se
curan más) y de la literatura disponible. Es **MNAR** (missing not at random).

Consecuencias:

- El recall medido sobre el gold es un recall **sobre lo curado**, no sobre la
  biología. Un TSS nuevo con señal fuerte no es un fallo.
- La precisión no se puede estimar sin negativos. Cualquier "precisión 0,86" del
  reporte es precisión **dentro del gold**, condicionada a que el TSS ya está en
  el gold.
- Los estimadores PU (Elkan-Noto) dan una cota, no un número exacto, y dependen
  de SCAR. Sirven para comparar modelos, no para publicar una precisión absoluta.

### 6.3 Diseño de negativos: cuatro opciones y sus trampas

| Tipo de negativo | Cómo se construye | Trampa | Uso recomendado |
|---|---|---|---|
| **Aleatorio genómico** | Ventanas al azar del mismo genoma | Puede caer dentro de un promotor no anotado | Solo como fondo nulo para z/p |
| **Apareado por composición** | *Shuffle* de dinucleótidos o coincidencia de GC | Destruye estructura de señal (CpG, oligo-A) y hace el problema más fácil de lo real | Para test de especificidad de motivo |
| **Regiones intergénicas lejanas** | > 500 pb de cualquier TSS/operón | Asume que "lejos" = no regulado | Negativo "duro" razonable |
| **Negativos experimentales** | TSS medidos sin sigma asignado o sin actividad | Requiere datos que el proyecto no tiene | Ideal, no disponible |

La literatura de evaluación de motivos ya advirtió que el desempeño depende del
conjunto de evaluación y de los controles usados, y que ningún método gana en
todos los conjuntos (Tompa et al., 2005). La misma advertencia aplica a `sigma`:
un AUC alto contra dinucleótidos barajados no predice AUC contra genoma real.

Para construir negativos de función (no de motivo) la referencia es NoGO (Youngs
et al., 2014): los negativos no se registran en las bases; hay que construirlos y
declarar el criterio. Es el análogo exacto del problema de `sigma`.

### 6.4 Auditoría de ruido (recomendada, no implementada)

1. **Consistencia dentro de ortogrupo.** Para cada OG con ≥ 5 miembros, medir la
   fracción de etiquetas de familia concordantes. OG0001 ya muestra discordancia.
2. **Concordancia entre fuentes.** Para K-12, comparar RegulonDB (Salgado et al.,
   2024) con EcoCyc; para *B. subtilis*, DBTBS (Sierro et al., 2008) con SubtiWiki
   (Elfmann et al., 2025). Reportar κ de Cohen por clase.
3. **Confident learning.** Northcutt et al. (2021) permiten estimar la matriz de
   ruido y rankear ejemplos sospechosos sin etiquetas limpias. Aplicable a nivel
   de clase de sigma.
4. **Etiquetado programático.** Snorkel (Ratner et al., 2017) formaliza combinar
   reglas débiles (nombre, dominio, ortología, sintenía) con un modelo de
   etiquetas que aprende pesos y permite que cada regla se abstenga. Es la vía
   natural para dar etiquetas débiles a las 400 y auditar su cobertura, en lugar
   de inventar un gold.
5. **Adjudicación humana de una muestra.** 100–200 TSS estratificados por especie
   y clase, con doble anotación y κ. Es lo único que calibra el resto.

---

## 7. Calibración y abstención

### 7.1 Definiciones y métricas

Un clasificador está **calibrado** si, entre las predicciones con confianza 0,8,
aproximadamente el 80 % son correctas. No es lo mismo que ser **discriminativo**
(ordenar bien). El repo ya mide ECE: la capa de contexto bajó ECE de 0,061 a
0,053 sin cambiar top-1 ([`SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md)). Ese es el
patrón típico: calibrar no mejora el ranking, mejora la interpretabilidad de la
confianza.

| Métrica | Qué mide | Cuándo usarla |
|---|---|---|
| ECE / MCE | Brecha media/máxima confianza−exactitud | Reporte por especie y por clase |
| Brier | Error cuadrático de la probabilidad | Comparar modelos probabilísticos |
| NLL | Verosimilitud | Entrenamiento y comparación |
| Curva de fiabilidad | Forma del desajuste | Diagnóstico (¿sobreconfianza en la cola?) |
| Riesgo-cobertura / AURC | Riesgo al abstenerse | Política de abstención |

Guo et al. (2017) mostraron que las redes modernas están **peor** calibradas que
las de hace una década y que *temperature scaling* corrige mucho con un solo
parámetro. Ovadia et al. (2019) añaden lo importante para este proyecto: **la
calibración se degrada bajo *dataset shift***, y las métricas de calibración
medidas en distribución no predicen el comportamiento fuera de distribución.

### 7.2 Calibración entre especies: lo que no se puede

- **La temperatura no transfiere.** Un modelo calibrado en K-12 puede estar
  sobreconfiado en *B. subtilis*. Hay que recalibrar por especie, y para eso hace
  falta un conjunto etiquetado por especie.
- **Sin etiquetas en destino, la calibración es una suposición.** Se puede
  intentar corrección de *label shift* (Lipton et al., 2018): estimar el prior de
  clases en destino a partir de las predicciones y la matriz de confusión de
  origen, y reponderar. Funciona si `P(x|y)` es estable y el prior cambia; es
  exactamente el escenario "misma biología, distinta frecuencia de uso de sigma".
  Falla si el cambio es de concepto (σ38 vs SigB).
- **Conformal prediction** (Angelopoulos y Bates, 2021) da cobertura marginal
  garantizada **bajo intercambiabilidad**. Las especies no son intercambiables:
  Barber et al. (2023) acotan la pérdida de cobertura en términos de la
  discrepancia entre la distribución de calibración y la de test. La versión
  práctica es *Mondrian* por clado (una calibración por filo) o conformal
  ponderado; ambas requieren datos de calibración por clado, que para los filos
  sin gold no existen.

### 7.3 Abstención: tres estados, no dos

El repo ya distingue estados que **no** deben colapsarse:

| Estado | Significado | ¿Es abstención del modelo? |
|---|---|---|
| `assign` | Hay evidencia suficiente y el margen supera el umbral | No |
| `unresolved` / `ambiguous` | El modelo tiene señal pero no alcanza el umbral, o dos sigmas compiten | **Sí** (abstención por incertidumbre) |
| `no_model` | No hay PWM/modelo para esa familia o esa especie | No: es falta de datos, no incertidumbre |
| `sin_sospecha` | No hay señal por encima del fondo | No: es una predicción negativa |

Mezclar `no_model` con `unresolved` infla artificialmente la "abstención" y
oculta que el problema es de cobertura de modelos. El reporte debe separarlos.

**Política de abstención (selective classification).** Geifman y El-Yaniv (2017)
formalizan el umbral de confianza que garantiza un riesgo objetivo con alta
probabilidad, **bajo iid**. Bajo *shift* entre especies la garantía no se
hereda; hay que re-estimar el umbral por especie o aceptar que la garantía es
solo nominal. La curva riesgo-cobertura es la herramienta de reporte: muestra
cuánto recall por clase se gana al precio de cuánta cobertura.

### 7.4 Lo que ya muestran los números del repo

| Configuración | Top-1 | Macro-F1 | Mayoría trivial | Lectura |
|---|---:|---:|---:|---|
| K-12, `argmax z` | 0,5667 | 0,4645 | 0,6127 | Pierde contra la mayoría |
| K-12, jerárquico logístico | 0,7057 | 0,5099 | 0,6127 | Gana, pero entrenado y evaluado en el mismo gold |
| *B. subtilis*, PWM propia | 0,2965 | 0,1536 | 0,8756 | La PWM sola no basta |
| *B. subtilis*, k-mer 4-mer | 0,8531 | 0,4191 | 0,8756 | Casi mayoría, mejor macro-F1 |
| *B. subtilis*, híbrido | 0,9003 | 0,5184 | 0,8756 | **Iguala** la mayoría en top-1; macro-F1 2,8× |

El caso del híbrido de *B. subtilis* es el más ilustrativo: **0,9003 de top-1
suena excelente y es 0,8756 de mayoría**. La mejora real está en macro-F1
(0,1867 → 0,5184): el modelo sí distingue alternativas, pero el top-1 no lo
refleja porque σ70 es el 87,6 % del gold. Si se reportara solo top-1, se
concluiría "casi perfecto" o "no sirve" según el baseline elegido.

> **Nota sobre denominadores.** Los valores de esta tabla son de
> [`SIGMA_STRATEGIES.md`](../SIGMA_STRATEGIES.md) (GroupKFold de 5 folds). Los de
> §0 y §8.1 (57,0 % / 62,6 %) son de
> [`SIGMA_LAYERED_EVAL.md`](../SIGMA_LAYERED_EVAL.md), que evalúa sobre el gold
> completo (3.189 TSS evaluables). Las diferencias de 0,3–1,3 pp son de
> denominador y de partición, no una contradicción; pero es exactamente el tipo
> de detalle que hay que fijar en un reporte: **qué conjunto, qué folds y qué
> denominador**. Un número sin esas tres cosas no es comparable.

---

## 8. Por qué la exactitud top-1 no basta

### 8.1 El argumento en una tabla (datos del repo)

| Especie | σ70 en gold | Mayoría top-1 | Modelo top-1 | Macro-F1 modelo | Veredicto con top-1 | Veredicto con macro-F1 |
|---|---:|---:|---:|---:|---|---|
| K-12 | 1.997/3.189 (62,6 %) | 0,613–0,626 | 0,567–0,570 | 0,465 | "Falla" | Falla, pero informa |
| *B. subtilis* | 1.949/2.226 (87,6 %) | 0,876 | 0,900 | 0,518 | "Casi perfecto" | Distingue alternativas |

Los rangos de K-12 corresponden a las dos corridas disponibles (GroupKFold de 5
folds y gold completo); ver la nota de §7.4. En ambas el modelo queda por debajo
de la mayoría trivial, así que la conclusión no depende de la partición.

Con dos especies y priors distintos, el mismo 0,90 significa cosas opuestas. La
exactitud top-1 es una función del prior y del umbral, no solo de la calidad del
modelo.

### 8.2 Los cinco mecanismos

1. **Dominancia del prior.** Con 6 clases y una que ocupa 60–90 %, la exactitud
   es casi la frecuencia de la clase mayoritaria. El *accuracy paradox* es
   estructural, no un accidente de este dataset (Saito y Rehmsmeier, 2015).
   Además, el desbalance interactúa con el diseño del entrenamiento (tamaño de
   batch, esquema de muestreo) y no se corrige solo reponderando la pérdida
   (Buda et al., 2018).
2. **Oculta el fallo por clase.** σ54 tiene precisión 0,169 y recall 0,333 en
   K-12; σ38 precisión 0,231 y recall 0,710. Un global de 57 % no dice nada sobre
   si σ54 es usable.
3. **No admite multi-etiqueta.** 124 TSS tienen más de un sigma; el top-1 no
   puede acertarlos todos y el agregado lo penaliza como error.
4. **Ignora el costo asimétrico.** Asignar σ70 a un promotor σ54 no es un error
   equivalente a lo contrario: la ruta de rescate busca activadores en
   [−100,−35] y no en −24/−12 (`SIGMA_AUDIT.md` §3). El error tiene consecuencia
   downstream, no solo estadística.
5. **No dice nada fuera del gold.** En las 400 no hay gold: la exactitud no es
   calculable, y sustituirla por "cobertura" sin decirlo es un cambio de
   definición silencioso.

### 8.3 Suite de métricas mínima

| Métrica | Qué responde | Cómo se reporta |
|---|---|---|
| **Macro-F1** | ¿Distingue clases raras? | Por especie y por fold |
| **Recall por clase** | ¿Qué sigma se detecta? | Con n y IC |
| **PR-AUC por clase** | Ranking independiente del umbral | Curvas, no un número |
| **Balanced accuracy** | Sensibilidad media por clase | Complemento de macro-F1 |
| **Top-k (k=3) `hit_any`/`hit_all`** | Ranking y multi-etiqueta | Separados, no sumados |
| **Brier / ECE** | ¿La confianza significa algo? | Por especie |
| **Riesgo-cobertura (AURC)** | ¿Cuánto cuesta abstenerse? | Curva + punto operativo |
| **Cobertura de inventario** | ¿Hay datos para puntuar? | % de genomas con modelo |
| **Tasa de abstención por causa** | ¿Falta modelo o falta señal? | `unresolved` vs `no_model` |

Reglas de reporte: (a) nunca una exactitud sin su mayoría trivial al lado; (b)
nunca un agregado sin desagregado por clase; (c) nunca una comparación entre
especies sin declarar el prior; (d) nunca un intervalo calculado sobre TSS
independientes (Johnson y Khoshgoftaar, 2019, para el tratamiento de desbalance;
Davis y Goadrich, 2006, para la relación PR/ROC).

---

## 9. Métodos aplicables a `sigma`, por capa

| Capa | Método recomendado | Qué problema resuelve | Supuesto | Costo |
|---|---|---|---|---|
| **0. Inventario** | Familia amplia como unidad + PU/multi-etiqueta; OG como unidad de transferencia | 62,8 % `unresolved`; ruido dentro de OG | El dominio define la familia mejor que el nombre | Bajo |
| **0. Inventario** | Etiquetado programático (Snorkel) con reglas nombre/dominio/ortología/sintenía | No hay gold en 400 | Las reglas son condicionalmente independientes (falso; usar modelo de etiquetas) | Medio |
| **1. Motivo** | PWM por familia ortóloga + fondo por especie (ya existe) | Composición distinta | El motivo se conserva dentro de familia | Ya hecho |
| **1. Motivo** | Encogimiento jerárquico de PWM (Dirichlet, `κ`) | σ54 con 11 sitios en *B. subtilis* | La familia es un prior razonable | Bajo |
| **1. Motivo** | Arquitectura como covariable (σ54 −24/−12; ECF; −10 extendida) | Un solo scanner no cubre todas | La arquitectura es discreta y conocida | Ya identificado |
| **2. Competición** | Softmax multinomial con intercepto por especie (*partial pooling*) | Umbrales que no transfieren | `P(x\|y)` estable entre especies | Medio |
| **2. Competición** | Margen top-2 + umbral calibrado por especie | Ambigüedad real entre sigmas | El margen separa ambiguos de claros | Bajo |
| **3. Asignación** | Multietiqueta + `hit_any`/`hit_all` | 124 TSS multi-sigma | Un TSS puede tener varios sigma | Medio |
| **4. Contexto** | Log-odds con pesos por especie estimados por pooling parcial | El contexto ayuda en una especie y no en otra | Las evidencias son informativas donde existen | Medio |
| **4. Contexto** | Manejo explícito de modalidad faltante (peso 0, no imputar) | `pathway_ec` 27,5 % vs 0 %; conservación 23,4 % vs 4,5 % | Falta ≠ ausencia | Ya implementado |
| **5. Calibración** | Temperatura por especie con prior jerárquico | Confianza no comparable entre especies | El modelo es el mismo salvo escala | Bajo |
| **5. Abstención** | Riesgo-cobertura + tres estados | No se puede asignar en todo | El costo del error es asimétrico | Bajo |
| **6. Validación** | LOSO (7) + LOCO (3 filos) + GroupKFold intra-especie | Confundir interpolación con generalización | Las particiones son las correctas | Medio |

**Nota sobre σ54.** No es "una clase más". Usa −24/−12, depende de bEBP y de IHF,
y su geometría es distinta. Tratarlo en el mismo softmax que σ70 obliga al modelo
a comparar scores de arquitecturas distintas; la alternativa es un modelo
jerárquico en dos etapas (¿arquitectura −24/−12? → clase dentro de la
arquitectura). El repo ya detectó esto (`SIGMA_AUDIT.md` §4.2: "No reutilizar el
scanner σ70 para σ54").

---

## 10. Supuestos, con su prueba y su plan B

| # | Supuesto | Por qué se necesita | Cómo se testea | Si se viola |
|---|---|---|---|---|
| 1 | La gramática del promotor se conserva dentro de la familia | Permite transferir el PWM | AUC del PWM de familia en especie held-out vs PWM propia | Modelo por clado; declarar `no_model` |
| 2 | El fondo nulo depende de la especie pero no del TSS | z y p comparables | Comparar distribución de z en regiones intergénicas por especie | Fondo empírico por región |
| 3 | El mapeo de etiquetas entre especies es correcto | Multi-tarea y evaluación conjunta | Concordancia de etiquetas dentro de OG; κ entre fuentes | Familia amplia; `unresolved` |
| 4 | La ausencia de etiqueta no es información | PU/negativos | Tasa de curación por especie/clase | Modelo MNAR explícito |
| 5 | `P(x\|y)` es estable entre especies | Corrección de *label shift* | Rendimiento de la corrección en LOSO | Reentrenar por clado |
| 6 | Las especies son intercambiables (supuesto de CP estándar) | Intervalos y cobertura | Distancia filogenética / Jaccard entre especies | Calibración por clado; intervalos más anchos |
| 7 | El árbol filogenético es correcto | PGLS/mixed models | Comparar árbol 16S vs GTDB vs ANI | Efectos aleatorios por género/filo, sin árbol |
| 8 | Los TSS de una especie son independientes entre sí | Intervalos y *p*-valores | Autocorrelación espacial de residuos | Bootstrap por bloque/especie |
| 9 | El gold es completo en la región evaluada | Declarar falsos positivos | Cobertura negativa del gold | Solo `outside_gold`/`novel_candidate` |
| 10 | La práctica de curación es estable en el tiempo | Comparar versiones | Diff entre versiones de RegulonDB/DBTBS | Congelar *snapshot* y declarar fecha |
| 11 | El prior de clase se puede estimar en destino | Umbrales y métricas | BBSE en especies con gold | No reportar exactitud; solo cobertura |
| 12 | La ventana evaluada captura el promotor | Comparabilidad | Barrido de ventana (sensibilidad) | Reportar sensibilidad a la ventana |

Los supuestos 1, 3 y 9 son los que más afectan las conclusiones de las 400. Los
tres son testeables con datos que el proyecto ya tiene.

---

## 11. Cómo validar: protocolo en ocho fases

**Fase A — Congelar el diseño (1 día).**
Manifiesto de particiones (**propuesto**; hoy no existe `data/sigma/splits/`):

```text
data/sigma/splits/
  loso.tsv        # 7 folds, un genoma con gold por fold
  loco.tsv        # 3 folds por filo (Pseudomonadota, Bacillota, Actinomycetota)
  intra.tsv       # GroupKFold por bloques de 10 kb (control de interpolación)
  interp_check.tsv# Jaccard/ANI de cada held-out contra train
```

**Fase B — Baselines por fold (1 día).**
Cuatro baselines, todos recalculados dentro de cada fold: (1) mayoría σ70; (2)
prior de familia; (3) PWM global argmax; (4) PWM propia argmax. Sin estos cuatro
números, ningún resultado posterior es interpretable. **Aquí se resuelve también
la anomalía del 44,8 %**: descomponer la ganancia de la matriz ajena en cambio de
prior vs cambio de verosimilitud.

**Fase C — Modelo jerárquico (3–5 días).**
Softmax multinomial con `α_σ + u_{σ,s} + f_σ(x)`, PWM con encogimiento, y
selección de `κ`, `τ` y ventana por LOSO interno. Comparar contra los cuatro
baselines.

**Fase D — Multietiqueta y abstención (2 días).**
Umbral de margen por especie, tres estados de salida, curvas riesgo-cobertura por
clase.

**Fase E — Calibración (1 día).**
Temperatura por especie con prior; ECE y Brier por especie y clase; comparar con
ECE del modelo sin calibrar.

**Fase F — Auditoría de etiquetas y negativos (2–3 días).**
Consistencia por OG, κ entre fuentes, confianza de etiquetas, y una muestra
adjudicada. Construir los dos conjuntos negativos (aleatorio y apareado) y
reportar AUC contra ambos por separado.

**Fase G — Estrés (1 día).**
Genomas AT-ricos, genomas pequeños, contigs fragmentados, especies sin TSS
(ventana proxy de inicio de operón), especies sin ningún candidato sigma (40
genomas del inventario). Reportar caída de cobertura y de calibración.

**Fase H — Reporte (1 día).**
Una tabla por fase, con n, IC por *cluster* y el número de folds donde el modelo
supera al baseline. Declarar explícitamente qué afirmación autoriza cada tabla.

**Criterios de aceptación (propuestos, no de la literatura):**

- [ ] Las 7 particiones LOSO y las 3 LOCO están congeladas y versionadas.
- [ ] Cada métrica tiene su baseline por fold al lado.
- [ ] La anomalía 44,8 % vs 29,7 % está explicada con matriz de confusión.
- [ ] Se reporta macro-F1 y recall por clase; no solo top-1.
- [ ] Se reporta ECE por especie; no solo global.
- [ ] La abstención se desglosa por causa (`unresolved` vs `no_model`).
- [ ] Ninguna afirmación de generalización se basa en un `GroupKFold` intra-especie.
- [ ] El 62,8 % `sigma_unresolved` se reporta como límite del espacio de etiquetas.

---

## 12. Errores comunes

| # | Error | Por qué engaña | Dónde aparece | Corrección |
|---|---|---|---|---|
| 1 | Partir aleatoriamente por TSS | Ventanas solapadas y operones vecinos filtran información | Cualquier `train_test_split` ingenuo | `GroupKFold` por bloque; LOSO para generalizar |
| 2 | Llamar "generalización" a un `GroupKFold` intra-especie | Mide interpolación espacial | `SIGMA_STRATEGIES.md` (correctamente rotulado como evaluación interna) | Nombrar cada partición por lo que rompe |
| 3 | Ajustar el umbral en la especie de test | El umbral deja de ser transferible | Riesgo en `72_motif_calibration.py` | Selección anidada |
| 4 | Comparar exactitud entre especies con priors distintos | 0,90 en *B. subtilis* es la mayoría; 0,57 en K-12 no | Reportes de top-1 | Macro-F1 + mayoría por fold |
| 5 | Usar "no está en el gold" como negativo | El gold es PU, no exhaustivo | Clasificación `false_positive` | Exigir cobertura negativa (`SIGMA_LAYERED_EVAL.md` §5) |
| 6 | Transferir el umbral junto con el PWM | El umbral depende del fondo y del prior | Umbrales K-12 → *B. subtilis* | Umbral por especie con pooling |
| 7 | Ignorar la multi-etiqueta | 124 TSS con >1 sigma | Agregados de top-1 | `hit_any`/`hit_all`, Jaccard |
| 8 | Confiar en un flag sin verificar procedencia | El caso 44,8 %: `--pwm` ignorado | `50_eval_layers.py` (corregido) | `model_source` en la salida |
| 9 | Calibrar con el mismo gold que se evalúa | ECE optimista | Cualquier `temperature scaling` ingenuo | Calibración dentro del fold |
| 10 | Reportar un agregado sin desagregar | Un modelo puede ganar en 2 especies y perder en 5 | Reportes de "media de folds" | Por fold y por clase |
| 11 | Tratar 400 genomas como muestras iid | Están correlacionados por filo | `--all` sobre genomas | Bootstrap por clado |
| 12 | Asumir que el vecino k-mer identifica especie | `NEW_SPECIES.md`: vecino ≠ ANI | Asignación de especie nueva | ANI o GTDB; declarar aproximación |
| 13 | Contar `no_model` como abstención | Infla la abstención y oculta falta de cobertura | Reportes de cobertura | Tres estados separados |
| 14 | Usar el nombre del gen como etiqueta de familia | `rseD` → `sigmaE` | Inventario | Dominio + ortología; `unresolved` |
| 15 | Declarar "modelo multi-especie" con un calibrador por especie | `SIGMA_STRATEGIES.md`: "No se declara modelo multi-especie" | Tentación de marketing | Un modelo con efectos por especie, o nombrarlo "ensemble por especie" |

---

## 13. Recomendación concreta para las 400 especies

### 13.1 Tres niveles de ambición, explícitos

| Nivel | Qué se entrega | Qué se afirma | Requisito |
|---|---|---|---|
| **A. Inventario + evidencia** | Gen → familia; promotor candidato con z, p, percentil, arquitectura | "Hay señal compatible con la familia X" | Ya disponible para ~400 |
| **B. Asignación con abstención** | σ asignado o `unresolved`, con margen y calibración | "El modelo asigna σ con esta cobertura y esta calibración" | Modelo jerárquico + LOSO/LOCO |
| **C. Asignación calibrada por especie** | Probabilidad calibrada y conjunto conformal | "Cobertura 1−α bajo el supuesto de intercambiabilidad dentro del clado" | Datos de calibración por clado (hoy: 3 filos) |

La propuesta es entregar **A para las 400, B para las 400 con las compuertas
superadas, C solo para los filos con gold**. No prometer C donde no hay datos.

### 13.2 Modelo recomendado (especificación)

```text
Entrada: ventana −80/+20 (o la que fije la ablación), codificación one-hot
  ↓
Tronco compartido: CNN pequeña (3–4 capas) o PWM aprendida por familia
  ↓
Cabeza 1 (familia):      softmax sobre {σ70_like, ecf_like, σ54_like}
Cabeza 2 (arquitectura): −35/−10 | −24/−12 | −10 extendida
Cabeza 3 (multi-etiqueta): sigmas compatibles (sigmoid por clase)
  ↓
Intercepto por especie: u_{σ,s} ~ N(0, τ²_σ)   (partial pooling)
  ↓
Calibración: temperatura por especie con prior jerárquico
  ↓
Abstención: margen top-2 + umbral por especie → assign | unresolved | no_model
```

Por qué así y no de otra forma:

- **Tres cabezas y no una**: la arquitectura es discreta y conocida; forzar un
  softmax único sobre 6 clases compara scores de geometrías distintas (σ54).
- **Tronco pequeño**: con 5.540 pares TSS–sigma y 6 clases, un modelo grande
  memoriza la especie. PromoterAtlas usa 1,8 M de parámetros con 9 M de
  secuencias; aquí hay tres órdenes de magnitud menos datos.
- **Intercepto por especie**: es la forma mínima de permitir que el prior cambie
  sin reentrenar el motivo.
- **Sin adversarial**: con 7 dominios etiquetados y 400 sin etiquetas, GRL es
  inestable y MORALE (moment matching) es más simple y no añade parámetros.

**Alternativa de bajo riesgo (si el tiempo aprieta):** quedarse con PWM por
familia + fondo por especie + encogimiento jerárquico + abstención. Es menos
ambicioso, pero cada componente es auditable y no requiere GPU. Dado que el
modelo de k-meros ya iguala la mayoría trivial en *B. subtilis* y el jerárquico
gana en K-12, el listón para un modelo más complejo es **superar a ambos en
LOSO**, no en su propia especie.

### 13.3 Plan de validación con los datos de hoy

Los 7 genomas con gold cubren tres filos, lo que permite tres mediciones
distintas con distinto poder:

| Diseño | Pliegues | Qué mide | Cómo se reporta |
|---|---|---|---|
| Intra-especie (`GroupKFold` 10 kb) | 5 × 7 | Techo de interpolación | "Techo" — no es generalización |
| LOSO | 7 | Transferencia a especie no vista | Mediana, rango y tabla por fold |
| **LOCO por filo** | 3 | Transferencia a linaje no visto | 3 números; descriptivo, con IC amplio |
| **Interpolación controlada** | LOSO mismo género (si hay pares) | ¿Cuánto de la caída es filogenia? | Diferencia LOSO − intra |

Composición de los pliegues LOCO (verificable en `data/sigma/gold/`):

- **Pseudomonadota**: `ecoli_k12`, `paeruginosa`.
- **Bacillota**: `subtilis`, `saureus`.
- **Actinomycetota**: `mtb_h37rv`, `scoelicolor`, `cglutamicum`.

Con 2–3 especies por filo, el LOCO estima **cuánto cae** el desempeño al salir del
filo, no la varianza poblacional. Es suficiente para decidir si el modelo se
puede aplicar a un filo nuevo con `confidence=low`, y insuficiente para afirmar
generalización universal. Escribirlo así en el reporte es parte del resultado.

### 13.4 Compuertas de decisión (política del proyecto)

Propuesta de reglas numéricas; son política interna, no umbrales de la
literatura. Se fijan **antes** de ver LOSO/LOCO.

| Compuerta | Criterio | Si no se cumple |
|---|---|---|
| G1. Señal | Macro-F1 LOSO > macro-F1 de la mayoría por fold en ≥ 5/7 | No publicar asignación; solo inventario + evidencia |
| G2. Clado | Macro-F1 LOCO > mayoría en ≥ 2/3 filos | Limitar la asignación a los filos donde sí pasa |
| G3. Clases raras | Recall σ54 ≥ 0,30 y σ32/σ38 ≥ 0,40 con abstención ≤ 40 % | Marcar esas clases como `unresolved` por defecto |
| G4. Calibración | ECE ≤ 0,05 dentro de especie tras calibrar | No publicar probabilidades; solo ordinales |
| G5. Cobertura | ≥ 80 % de los genomas con al menos una familia modelada | Declarar `no_model` explícito por genoma |

Si G1 falla, el resultado honesto sigue siendo útil: inventario completo,
evidencia de promotor y abstención. Es preferible a un número de exactitud que no
significa nada.

Los umbrales están calibrados al desempeño actual (σ54 recall 0,333 en K-12;
σ70 recall 0,576): son un **piso para no retroceder**, no una meta. Si al
implementar la jerarquía el recall de σ54 baja de 0,30, el modelo nuevo es peor
que el `argmax z` vigente en la clase que más cuesta, y eso invalida el cambio
aunque el macro-F1 global suba.

### 13.5 Qué reportar a usuarios y revisores

Permitido:

- "Cobertura del inventario: X genomas, Y candidatos, Z sin candidato."
- "Cobertura de modelos: N familias con PWM en M genomas."
- "En las 7 especies con gold, recall por clase y macro-F1 con IC por clúster."
- "En las 400 sin gold: distribución de scores, tasa de abstención y su causa."
- "El modelo no se calibró en especies sin gold; la confianza es ordinal."

Prohibido:

- "Exactitud del X % en 400 especies."
- "Generaliza a bacterias" a partir de 3 filos.
- "Precisión" sin conjunto negativo.
- Cualquier comparación entre especies sin su mayoría trivial al lado.

### 13.6 Orden de ejecución (2–3 semanas, sin GPU)

1. **Día 1** — Manifiesto de particiones y control de fuga (Fase A).
2. **Días 2–3** — Cuatro baselines por fold; explicar la anomalía 44,8 %/29,7 %
   (Fase B). Este paso solo puede cambiar conclusiones previas.
3. **Días 4–8** — Modelo jerárquico + encogimiento de PWM; LOSO; comparación con
   baselines (Fase C).
4. **Días 9–10** — Multietiqueta y abstención (Fase D).
5. **Día 11** — Calibración por especie (Fase E).
6. **Días 12–14** — Auditoría de etiquetas y negativos (Fase F).
7. **Día 15** — Estrés y reporte (Fases G–H).

Si solo hay una semana: fases A, B y E. Los baselines y la calibración son los
que cambian decisiones; el modelo nuevo es lo último que hay que añadir.

### 13.7 Lo que no hay que hacer

- No entrenar con las 400 y evaluar con un subconjunto aleatorio de las 400: es
  la misma fuga con otro nombre.
- No usar el gold de K-12 para elegir el modelo y reportar *B. subtilis* como
  validación externa limpia: *B. subtilis* ya se usó para decidir políticas
  (`SIGMA_STRATEGIES.md`). Está contaminado como test.
- No afirmar que un modelo "aprende biología" sin atribución por posición.
- No extrapolar el resultado de Gammaproteobacteria (PromoterAtlas) a los filos
  del catálogo.

---

## 14. Referencias verificadas

Verificación: DOI resueltos vía Crossref (`api.crossref.org/works/{doi}`) el
2026-11-26; PMID vía Europe PMC/NCBI eutils; arXiv vía `export.arxiv.org`.
Todas las entradas de esta tabla devolvieron título, venue y año coherentes con
la cita. Las que no se pudieron resolver se omitieron en lugar de citarse de
memoria.

### Transferencia y adaptación de dominio

| Referencia | Año | Venue | DOI / URL | PMID | Aporta |
|---|---|---|---|---|---|
| Pan y Yang | 2010 | IEEE TKDE | [10.1109/TKDE.2009.191](https://doi.org/10.1109/TKDE.2009.191) | — | Taxonomía de transferencia; *negative transfer* |
| Ben-David et al. | 2010 | Machine Learning | [10.1007/s10994-009-5152-4](https://doi.org/10.1007/s10994-009-5152-4) | — | Cota de error en dominio destino (divergencia + λ) |
| Ebeid et al. (MORALE) | 2026 | RECOMB-CG (LNCS) | [10.1007/978-3-032-26891-4_19](https://doi.org/10.1007/978-3-032-26891-4_19) | — | Alineación de momentos sin etiquetas destino; preprint [10.1101/2025.05.21.655414](https://doi.org/10.1101/2025.05.21.655414) |
| Park et al. (SATL) | 2024 | Brief Bioinform | [10.1093/bib/bbae004](https://doi.org/10.1093/bib/bbae004) | 38305455 | Transferencia entre especies sin ortología |
| Coppens y Ledesma-Amaro (PromoterAtlas) | 2026 | Nat Commun | [10.1038/s41467-026-72837-3](https://doi.org/10.1038/s41467-026-72837-3) | 42140912 | 3.371 especies, 9 M secuencias, un filo |
| Zhang et al. (DeepCROSS) | 2025 | Nat Commun | [10.1038/s41467-025-57031-1](https://doi.org/10.1038/s41467-025-57031-1) | 39971994 | Meta-representación + multi-tarea cross-species |
| Van Brempt et al. | 2020 | Nat Commun | [10.1038/s41467-020-19446-w](https://doi.org/10.1038/s41467-020-19446-w) | 33199691 | Modelos por sigma; ortogonalidad validada |
| Johns et al. | 2018 | Nat Methods | [10.1038/nmeth.4633](https://doi.org/10.1038/nmeth.4633) | — | Elementos con expresión selectiva de especie |

### Modelos de secuencia y multi-tarea

| Referencia | Año | Venue | DOI / URL | PMID | Aporta |
|---|---|---|---|---|---|
| Avsec et al. (Enformer) | 2021 | Nat Methods | [10.1038/s41592-021-01252-x](https://doi.org/10.1038/s41592-021-01252-x) | — | Secuencia→expresión con atención de largo alcance |
| Kelley et al. (Basenji) | 2018 | Genome Res | [10.1101/gr.227819.117](https://doi.org/10.1101/gr.227819.117) | — | Predicción multi-tarea por cromosoma |
| Zhou y Troyanskaya (DeepSEA) | 2015 | Nat Methods | [10.1038/nmeth.3547](https://doi.org/10.1038/nmeth.3547) | — | Efectos de variantes no codificantes |
| Chen et al. (Sei) | 2022 | Nat Genet | [10.1038/s41588-022-01102-2](https://doi.org/10.1038/s41588-022-01102-2) | 35817977 | 21.907 perfiles; clases de actividad regulatoria |
| Ji et al. (DNABERT) | 2021 | Bioinformatics | [10.1093/bioinformatics/btab083](https://doi.org/10.1093/bioinformatics/btab083) | — | Transformer preentrenado en DNA |
| Dalla-Torre et al. (Nucleotide Transformer) | 2025 | Nat Methods | [10.1038/s41592-024-02523-z](https://doi.org/10.1038/s41592-024-02523-z) | 39609566 | 850 genomas; *fine-tuning* con pocos datos |
| Nguyen et al. (Evo) | 2024 | Science | [10.1126/science.ado9336](https://doi.org/10.1126/science.ado9336) | 39541441 | 2,7 M genomas procariotas |
| Ito et al. (Mouse-Geneformer) | 2025 | PLoS Genet | [10.1371/journal.pgen.1011420](https://doi.org/10.1371/journal.pgen.1011420) | — | Transferencia humano→ratón; lo conservado transfiere mejor |
| Caruana | 1997 | Machine Learning | [10.1023/A:1007379606734](https://doi.org/10.1023/A:1007379606734) | — | Fundamento de multi-tarea |
| Gulrajani y Lopez-Paz (DomainBed) | 2021 | ICLR | [arXiv:2007.01434](https://arxiv.org/abs/2007.01434) | — | ERM como rival fuerte; selección de modelo honesta |

### Jerárquicos y filogenéticos

| Referencia | Año | Venue | DOI / URL | PMID | Aporta |
|---|---|---|---|---|---|
| Gelman | 2006 | Technometrics | [10.1198/004017005000000661](https://doi.org/10.1198/004017005000000661) | — | Qué puede y qué no el modelado multinivel |
| Harrison et al. | 2018 | PeerJ | [10.7717/peerj.4794](https://doi.org/10.7717/peerj.4794) | — | Guía práctica de modelos mixtos |
| Pagel | 1999 | Nature | [10.1038/44766](https://doi.org/10.1038/44766) | — | λ filogenético |
| Ives y Garland | 2010 | Syst Biol | [10.1093/sysbio/syp074](https://doi.org/10.1093/sysbio/syp074) | — | Regresión logística filogenética |
| Münkemüller et al. | 2012 | Methods Ecol Evol | [10.1111/j.2041-210X.2012.00196.x](https://doi.org/10.1111/j.2041-210X.2012.00196.x) | — | Medir señal filogenética |
| Gardner et al. | 2025 | Nat Commun | [10.1038/s41467-025-61036-1](https://doi.org/10.1038/s41467-025-61036-1) | 40610418 | Intervalos de predicción filogenéticos |
| Sjölander et al. | 1996 | CABIOS | [10.1093/bioinformatics/12.4.327](https://doi.org/10.1093/bioinformatics/12.4.327) | — | Priors de Dirichlet para PWM |
| Field y Welsh | 2007 | JRSS-B | [10.1111/j.1467-9868.2007.00593.x](https://doi.org/10.1111/j.1467-9868.2007.00593.x) | — | Bootstrap por clúster |

### Diseño de validación y fuga

| Referencia | Año | Venue | DOI / URL | PMID | Aporta |
|---|---|---|---|---|---|
| Roberts et al. | 2017 | Ecography | [10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881) | — | *Block CV* con estructura filogenética |
| Varma y Simon | 2006 | BMC Bioinformatics | [10.1186/1471-2105-7-91](https://doi.org/10.1186/1471-2105-7-91) | — | Sesgo de CV al seleccionar y evaluar |
| Nadeau y Bengio | 2003 | Machine Learning | [10.1023/A:1024068626366](https://doi.org/10.1023/A:1024068626366) | — | Inferencia sobre el error de CV |
| Kapoor y Narayanan | 2023 | Patterns | [10.1016/j.patter.2023.100804](https://doi.org/10.1016/j.patter.2023.100804) | 37720327 | Ocho tipos de fuga |
| Whalen et al. | 2022 | Nat Rev Genet | [10.1038/s41576-021-00434-9](https://doi.org/10.1038/s41576-021-00434-9) | 34837041 | Trampas de ML en genómica |
| Kim et al. | 2022 | Clin Microbiol Rev | [10.1128/cmr.00179-21](https://doi.org/10.1128/cmr.00179-21) | 35612324 | ML en genómica microbiana: límites de benchmark |

### Etiquetas, negativos y calibración

| Referencia | Año | Venue | DOI / URL | PMID | Aporta |
|---|---|---|---|---|---|
| Elkan y Noto | 2008 | KDD | [10.1145/1401890.1401920](https://doi.org/10.1145/1401890.1401920) | — | Aprendizaje PU |
| Youngs et al. (NoGO) | 2014 | PLoS Comput Biol | [10.1371/journal.pcbi.1003644](https://doi.org/10.1371/journal.pcbi.1003644) | 24922051 | Los negativos no se registran; hay que construirlos |
| Northcutt et al. | 2021 | JAIR | [10.1613/jair.1.12125](https://doi.org/10.1613/jair.1.12125) | — | *Confident learning* |
| Frénay y Verleysen | 2014 | IEEE TNNLS | [10.1109/tnnls.2013.2292894](https://doi.org/10.1109/tnnls.2013.2292894) | 24808033 | Efecto no monótono del ruido |
| Schnoes et al. | 2009 | PLoS Comput Biol | [10.1371/journal.pcbi.1000605](https://doi.org/10.1371/journal.pcbi.1000605) | — | Error de anotación que se propaga |
| Ratner et al. (Snorkel) | 2017 | VLDB | [10.14778/3157794.3157797](https://doi.org/10.14778/3157794.3157797) | — | Etiquetado programático + abstención |
| Tompa et al. | 2005 | Nat Biotechnol | [10.1038/nbt1053](https://doi.org/10.1038/nbt1053) | — | El desempeño depende de negativos y fondo |
| Guo et al. | 2017 | ICML | [arXiv:1706.04599](https://arxiv.org/abs/1706.04599) | — | Redes mal calibradas; temperature scaling |
| Ovadia et al. | 2019 | NeurIPS | [arXiv:1906.02530](https://arxiv.org/abs/1906.02530) | — | Calibración se degrada bajo *shift* |
| Geifman y El-Yaniv | 2017 | NeurIPS | [arXiv:1705.08500](https://arxiv.org/abs/1705.08500) | — | Clasificación selectiva con riesgo objetivo |
| Angelopoulos y Bates | 2021 | arXiv | [arXiv:2107.07511](https://arxiv.org/abs/2107.07511) | — | Introducción a conformal prediction |
| Barber et al. | 2023 | Ann Stat | [10.1214/23-AOS2276](https://doi.org/10.1214/23-AOS2276) | — | Conformal sin intercambiabilidad |
| Lipton et al. | 2018 | ICML | [arXiv:1802.03916](https://arxiv.org/abs/1802.03916) | — | Corrección de *label shift* con predictor caja negra |

### Métricas y biología de sigma

| Referencia | Año | Venue | DOI / URL | PMID | Aporta |
|---|---|---|---|---|---|
| Saito y Rehmsmeier | 2015 | PLoS ONE | [10.1371/journal.pone.0118432](https://doi.org/10.1371/journal.pone.0118432) | — | PR > ROC con desbalance |
| Davis y Goadrich | 2006 | ICML | [10.1145/1143844.1143874](https://doi.org/10.1145/1143844.1143874) | — | Relación PR/ROC |
| Johnson y Khoshgoftaar | 2019 | J Big Data | [10.1186/s40537-019-0192-5](https://doi.org/10.1186/s40537-019-0192-5) | — | Desbalance en deep learning |
| Buda et al. | 2018 | Neural Networks | [10.1016/j.neunet.2018.07.011](https://doi.org/10.1016/j.neunet.2018.07.011) | — | Estudio sistemático de desbalance en CNN |
| Feklistov et al. | 2014 | Annu Rev Microbiol | [10.1146/annurev-micro-092412-155737](https://doi.org/10.1146/annurev-micro-092412-155737) | — | Biología de factores sigma |
| Mejía-Almonte et al. | 2020 | Nat Rev Genet | [10.1038/s41576-020-0254-8](https://doi.org/10.1038/s41576-020-0254-8) | — | Redefinición de promotor/TSS en bacterias |
| Salgado et al. (RegulonDB v12) | 2024 | Nucleic Acids Res | [10.1093/nar/gkad1072](https://doi.org/10.1093/nar/gkad1072) | 37971353 | Gold de K-12 |
| Sierro et al. (DBTBS) | 2008 | Nucleic Acids Res | [10.1093/nar/gkm910](https://doi.org/10.1093/nar/gkm910) | 17962296 | Gold de *B. subtilis* |
| Elfmann et al. (SubtiWiki) | 2025 | Nucleic Acids Res | [10.1093/nar/gkae957](https://doi.org/10.1093/nar/gkae957) | 39441067 | Fuente alternativa para concordancia |
| Richter y Rosselló-Móra | 2009 | PNAS | [10.1073/pnas.0906412106](https://doi.org/10.1073/pnas.0906412106) | — | ANI 95–96 % como frontera de especie |
| Parks et al. (GTDB) | 2018 | Nat Biotechnol | [10.1038/nbt.4229](https://doi.org/10.1038/nbt.4229) | — | Taxonomía estandarizada por filogenia |
| Umarov y Solovyev | 2017 | PLoS ONE | [10.1371/journal.pone.0171410](https://doi.org/10.1371/journal.pone.0171410) | — | CNN para promotores procariotas/eucariotas |

---

## 15. Apéndice A — Glosario mínimo

| Término | Definición operativa en este documento |
|---|---|
| **Interpolación** | Evaluar en una región/genoma que comparte linaje y práctica de anotación con el entrenamiento. Mide memoria y suavizado. |
| **Generalización** | Evaluar en una unidad (especie, clado) cuya distribución no estuvo representada en el entrenamiento. |
| **LOSO** | *Leave-one-species-out*: una especie fuera por fold. |
| **LOCO** | *Leave-one-clade-out*: un clado (filo/clase) fuera por fold. |
| **PU** | Positivos y no etiquetados; el régimen real de este gold. |
| **SCAR** | *Selected completely at random*: supuesto de que la probabilidad de etiquetado no depende de `x`. Falso aquí. |
| **Partial pooling** | Estimar efectos por grupo encogiéndolos hacia la media global. |
| **ECE** | *Expected calibration error*: brecha media entre confianza y exactitud. |
| **AURC** | Área bajo la curva riesgo-cobertura. |
| **BBSE** | *Black-box shift estimation*: estimar el prior de clases en destino desde predicciones. |

## 16. Apéndice B — Métricas, con fórmula

```text
macro-F1      = (1/K) Σ_k 2 P_k R_k / (P_k + R_k)
balanced acc  = (1/K) Σ_k R_k
ECE           = Σ_b (n_b / N) · |acc(b) − conf(b)|
Brier         = (1/N) Σ_i (p_i − y_i)²
AURC          = ∫₀¹ riesgo(cobertura) d(cobertura)
Jaccard@k     = |pred_k ∩ true| / |pred_k ∪ true|
```

Interpretación: macro-F1 y balanced accuracy responden "¿distingue clases raras?";
ECE y Brier, "¿la confianza significa algo?"; AURC, "¿cuánto cuesta abstenerse?".
Ninguna responde "¿cuántas especies generalizan?" — para eso están LOSO y LOCO.

## 17. Apéndice C — Checklist de construcción de particiones

- [ ] La unidad de *holdout* está declarada (TSS, región, especie, clado).
- [ ] Se verificó que no hay genomas casi idénticos entre train y test (Jaccard/ANI).
- [ ] El baseline trivial se calcula **dentro** del fold.
- [ ] Los hiperparámetros se eligen con folds internos.
- [ ] La calibración se ajusta en train, no en test.
- [ ] El reporte incluye por fold y por clase, con n.
- [ ] Los intervalos se calculan por *cluster* (especie/clado).
- [ ] Se declara qué afirmación autoriza la partición (tabla §1.2).
- [ ] Las particiones están versionadas y son reproducibles por hash.

---

**Frase final.** La generalización entre especies no es una propiedad que se
declare: es una medición con un diseño. En `sigma`, el diseño disponible permite
tres mediciones honestas (interpolación, LOSO de 7 especies, LOCO de 3 filos) y
una conclusión que no depende de ninguna de ellas: mientras el top-1 de K-12 esté
por debajo de la mayoría trivial y el gold no tenga negativos, la afirmación
defendible para las 400 especies es **cobertura con abstención**, no exactitud.
