# Predicción computacional de promotores bacterianos y asignación de sigma

> **Tipo:** revisión de literatura y métodos, con verificación de referencias.
> **Fecha de revisión:** 2026-11-26 (America/Santiago).
> **Alcance:** promotores bacterianos, predicción específica de sigma y
> asignación de sigma; evaluación, negativos, transferencia entre especies y
> fugas de datos.
> **No incluye:** resultados nuevos de benchmark, ni código. Este documento no
> modifica el pipeline.
>
> **Convención de verificación.** «Verificado» significa: el DOI resuelve en
> Crossref con el título/revista/año indicados; el PMID se recuperó por
> E-utilities de PubMed; la URL se comprobó con una petición HTTP el
> 2026-11-26 (código de estado en el Apéndice A). Las afirmaciones sobre
> método, datos y métricas provienen del abstract y, cuando se indica, de texto
> completo accesible o de la documentación del servidor. Si un dato no está
> declarado en la fuente, se escribe **no declarado**; no se infiere.

---

## 0. Resumen ejecutivo

1. **Predecir un promotor no es asignar un sigma.** La mayor parte de las
   herramientas publicadas responde una pregunta binaria («¿esta ventana es un
   promotor σ70?»), no «¿qué sigma reconoce este TSS?». Solo un grupo pequeño
   hace clasificación multi-sigma explícita: bTSSfinder, BacPP,
   iPromoter-2L/2L2.0, PromoterLCNN y, a mayor escala, PromoterAtlas.
2. **σ54 es un caso aparte.** Su arquitectura (−24/−12, dependiente de
   activadores) hace que los modelos entrenados en σ70 no transfieran. La
   herramienta dedicada actual es ProPr54 (2025), entrenada con 446 sitios de
   33 especies.
3. **Las cifras publicadas no son comparables entre sí.** Cambian la ventana,
   el conjunto de negativos, la especie, la definición de positivo y el split.
   Un 97 % de accuracy en un artículo y un 90 % en otro pueden representar
   tareas distintas.
4. **Los negativos mandan.** Los negativos aleatorios con composición similar
   (usados por el benchmark de mSystems) premian a los modelos que explotan
   contenido AT; los negativos de secuencias codificantes sesgan aún más. Un
   buen AUC en ventanas no implica buena precisión a escala genómica.
5. **La transferencia entre especies casi no se evalúa.** E. coli domina el
   entrenamiento; Pseudomonas, Bacillus y cianobacterias aparecen sobre todo
   como tests. La evidencia disponible (SAPPHIRE, bTSSfinder, nuestro propio
   benchmark K-12 → B. subtilis) muestra que la transferencia directa es débil.
6. **La fuga de datos es estructural, no anecdótica.** Entrenar y testear con
   la misma base curada (RegulonDB/DBTBS), con ventanas solapadas, con PWMs
   estimadas en el mismo gold o con versiones consecutivas de la misma base
   produce estimaciones optimistas aunque el split se llame «independiente».
7. **Deep learning no elimina esos problemas.** Mejora la clasificación en
   ventana, pero sigue dependiendo de los negativos, del tamaño de datos y de
   la especie. El salto reciente es de escala (PromoterAtlas: transformer de
   1,8 M de parámetros, 9 M de secuencias, 3371 especies de
   Gammaproteobacteria), no una solución general.
8. **La disponibilidad es frágil.** El 2026-11-26: bTSSfinder responde 503,
   BacPP no resuelve, iPromoter-2L no responde, biosapphire.com (SAPPHIRE) no
   responde y pepper.molgenrug.nl ya no sirve PePPER sino ProPr. Pocos tienen
   código publicado (Sigma70Pred, PromoterAtlas, SigmoID; bTSSfinder era
   standalone); la mayoría son solo web.
9. **Para nuestro benchmark**, la comparación honesta es por niveles: T1
   (promotor binario), T2 (promotor σ-específico), T3 (asignación multi-clase
   en TSS) y T4 (anotación genómica). Nuestro módulo actual hace T2 (scoring
   PWM) + T3 (argmax) y se evalúa en T3 a nivel TSS; no debe presentarse como
   T4 ni como validación biológica.
10. **Recomendación central:** benchmark TSS-level, multi-clase, con negativos
    múltiples y declarados, split por locus y por especie, métricas
    top-k/macro-F1/MCC, calibración solo si hay probabilidades, y auditoría de
    fuga entre PWMs y gold. Detalle en §8.

---

## 1. Definiciones: cinco niveles que se confunden

La literatura usa «promoter prediction» para cosas distintas. Para comparar
hay que fijar la tarea, la unidad de evaluación y la etiqueta.

| Nivel | Pregunta | Salida típica | Unidad | Ejemplos |
|---|---|---|---|---|
| **T0. Mapeo de TSS** | ¿dónde empieza la transcripción? | coordenada + hebra + score | base genómica | dRNA-seq, Cappable-seq, bTSSfinder |
| **T1. Detección de promotor** | ¿esta secuencia es promotora? | binaria (promotor / no promotor) | ventana | BPROM, CNNProm, Sigma70Pred, SAPPHIRE |
| **T2. Promotor σ-específico** | ¿esta secuencia es un promotor del sigma X? | binaria por sigma (one-vs-rest) | ventana | ProPr54 (σ54), iPro70-FMWin (σ70) |
| **T3. Asignación de sigma** | dado un TSS/promotor, ¿qué sigma lo reconoce? | clase (σ70/σ38/σ32/σ28/σ24/σ54) o distribución | TSS o ventana | bTSSfinder, BacPP, iPromoter-2L, PromoterLCNN |
| **T4. Anotación genómica** | ¿qué promotores, sigma, RBS y terminadores hay en el genoma? | features coordinates + tipo | genoma | PromoterAtlas, ProPr, SigmoID |

Consecuencias:

- **T1 no implica T3.** Un detector binario de σ70 puede tener gran accuracy y
  no saber nada de σ38, σ32 o σ54. Muchas «herramientas de sigma» son en
  realidad detectores σ70 (Sigma70Pred, 70ProPred, iPro70-FMWin, CNNProm,
  SAPPHIRE, BPROM).
- **T2 no es T3.** Un score alto para σ32 en una ventana no dice que ese TSS
  sea de σ32 in vivo; dice que la secuencia se parece al motivo. La decisión
  T3 requiere comparar contra los otros sigma y contra un modelo de «ninguno».
- **T3 es biológicamente multi-etiqueta.** Un mismo TSS puede ser reconocido
  por más de un sigma in vitro, y la selección in vivo depende de
  disponibilidad, competencia por el core, activadores y factores de
  elongación (Gruber & Gross, 2003; Feklistov et al., 2014; Helmann, 2019).
  Nuestro propio gold de K-12 contiene TSS multi-sigma (124 en la evaluación
  de `SIGMA_LAYERED_EVAL.md`), y una fila por operón no puede representarlos.
- **T4 no se evalúa con accuracy de ventanas.** Se evalúa con precisión/recall
  a escala genómica, distancia al TSS real, y consistencia con la anotación de
  operones.

> **Regla para este documento:** cada vez que se cita una cifra se indica la
> tarea y el conjunto de evaluación. Las cifras de T1 no se usan para rankear
> herramientas de T3.

---

## 2. Representación y arquitectura

### 2.1 PWM/PSSM y escaneo de motivos

La representación clásica es una **matriz de pesos por posición (PWM/PSSM)**:
para cada posición de la caja se estima un log-odds frente a un modelo de
fondo, y se escanea la ventana sumando contribuciones independientes
(Stormo, 2000). En bacterias, el escaneo no es «una caja»: es una
**arquitectura** con restricciones geométricas.

| Elemento | Sigma típico | Geometría | Nota |
|---|---|---|---|
| Caja −35 y −10 | σ70 (grupo 1), σ38, σ32, σ28 | espaciador ~15–19 nt; −10 a ~7 nt del TSS | la distancia y la fase importan |
| Extended −10 (`TGn`) | σ70 | sustituye funcionalmente a −35 en algunos promotores | no modelada por PWMs de dos cajas |
| UP element | σ70 | ~−40 a −60, A/T rica | recluta la subunidad α |
| −24/−12 | σ54 (RpoN) | espaciador 10–11 nt; requiere activador ATPasa | arquitectura distinta, no «σ70 con otro motivo» |
| ECF (grupo 4) | σ24/σE, σW, σX… | −35/−10 atípicas, con motivos propios | los modelos σ70 no transfieren |
| Discriminador | varios | entre −10 y el TSS | afecta la fuerza de apertura |

Referencias de arquitectura: Shultzaberger et al. (2007) para la anatomía de
los promotores σ70 de E. coli; Barrios et al. (1999) para la compilación de
promotores σ54; Paget (2015) y Feklistov et al. (2014) para la clasificación
de familias; Helmann (2019) para selectividad de la iniciación.

Limitaciones del PWM puro:

- **Independencia entre posiciones.** Las cajas −35 y −10 no son
  independientes; la información mutua y la geometría acoplan sus
  contribuciones. Un PWM por caja y suma de scores ignora esa correlación.
- **Rigidez de espaciador.** Un PWM fijo no modela bien inserciones/deleciones
  ni la flexibilidad de fase; por eso muchos escáneres prueban varios
  espaciadores y se quedan con el mejor score.
- **Fondo y umbral.** El score no es un p-valor biológico. El umbral se elige
  sobre un conjunto (a menudo el mismo gold), y a escala genómica el número de
  ventanas multiplica los falsos positivos.
- **PWMs con pocos sitios.** Las matrices por sigma de RegulonDB/DBTBS se
  estiman con decenas o cientos de sitios; las menos representadas (σ54,
  σ28, ECF) son las más ruidosas.

### 2.2 De PWM a modelos aprendidos

| Familia | Representación | Fortalezas | Modos de fallo |
|---|---|---|---|
| PWM + geometría | matrices por caja, espaciador restringido | interpretable, pocos datos, inspeccionable | rígida; depende de la matriz; no capta contexto |
| HMM / perfiles | estados para cajas y espaciador | modela longitud variable | poca señal para motivos cortos; sobreajuste |
| SVM + k-mer/PseKNC | composición de k-tuplas, autocorrelación, propiedades fisicoquímicas | fuerte en datasets pequeños; rápido | explota composición (AT/GC); no localiza cajas; difícil de interpretar |
| CNN | filtros convolucionales sobre one-hot | aprende motivos y su disposición; buena en ventanas | necesita datos; memoriza; atribución difícil |
| CNN + RNN | convolución + recurrencia | orden y dependencias largas | más parámetros; menos interpretable |
| Transformer | atención sobre la secuencia | escala y transferencia entre especies | datos masivos; coste; riesgo de exposición en preentrenamiento |

La progresión histórica es visible en las herramientas: BacPP (reglas
extraídas de redes), 70ProPred/iPro70-FMWin/iPromoter-2L (SVM con features),
CNNProm/PromoterLCNN/ProPr54 (CNN), PromoterAtlas (transformer). Ninguna
familia elimina el problema de los negativos ni el de la etiqueta.

### 2.3 Qué optimiza cada modelo

- Un **clasificador binario** optimiza separar positivos de negativos; su
  score no está calibrado como probabilidad de que el TSS sea usado.
- Un **clasificador multi-clase** optimiza la clase correcta entre las
  presentes; si se entrena sin clase «ninguno/desconocido», fuerza una
  etiqueta para cualquier ventana.
- Un **anotador genómico** optimiza precisión a escala de genoma y suele
  incluir post-procesado (supresión de solapamientos, NMS, contexto de
  operón). Sus métricas no son comparables con las de un clasificador de
  ventanas.

---

## 3. Herramientas y métodos

> **Cómo leer las tablas.** «Validación» resume lo que el artículo declara
> (CV, conjunto independiente, genoma, benchmark externo). «Negativos» se
> detalla en §5.2 y en la Tabla 4. «Código/disponibilidad» se comprobó el
> 2026-11-26; los estados HTTP están en el Apéndice A.

### 3.1 Predicción genérica de promotor (T1), casi siempre σ70

| Herramienta | Método | Entrada | Salida | Especie | Validación reportada | Código / disponibilidad |
|---|---|---|---|---|---|---|
| **BPROM** | scoring de cajas −35/−10; método descrito en Solovyev & Salamov (2011), sin artículo dedicado | FASTA | promotores σ70 putativos, cajas, score | E. coli (modelo); uso extendido en otras bacterias | incluida en el benchmark de mSystems (2020), donde fue la peor del panel | solo web Softberry; sin código; URL no alcanzable desde este entorno |
| **CNNProm** | CNN sobre ventanas (Keras/Theano) | FASTA | probabilidad/score de promotor | E. coli σ70 y B. subtilis | Sn 0,90 / Sp 0,96 (E. coli σ70); Sn 0,91 / Sp 0,95 (B. subtilis) en sus datos | web Softberry (`CNNPromoter_b`); sin código |
| **Sigma70Pred** | SVM sobre ~8.000 features → 200; RFE | FASTA (una o muchas) | promotor σ70 / no promotor + score | E. coli K-12 | 741 σ70 / 1.400 no promotores (RegulonDB 9.0); 90,41 % acc, AUROC 0,953, MCC 0,794 en «independiente» RegulonDB 10.8 (1.134/638); 81,46 % en promotores constitutivos | web + standalone (GitHub `raghavagps/sigma70pred`, Zenodo); activo |
| **70ProPred** | SVM con PSTNPss + PseEIIP | FASTA | promotor σ70 / no promotor | E. coli (σ70 procariota) | 95,56 % acc, MCC 0,90 (jackknife) en sus datos | web (`server.malab.cn/70ProPred`); sin código verificado |
| **iPro70-FMWin** | múltiples ventanas + features mínimas | FASTA | promotor σ70 / no promotor | E. coli | mejor herramienta del benchmark mSystems en la mayoría de métricas | web; código no verificado |
| **SAPPHIRE** | red neuronal artificial sobre cajas −35/−10 | FASTA | probabilidad de promotor σ70 | *Pseudomonas aeruginosa* y *P. putida* | supera a herramientas establecidas en clasificación de promotores σ70 de Pseudomonas (según el artículo) | script Python 3 + web `biosapphire.com`; el dominio no responde (2026-11-26) |
| **SAPPHIRE.CNN** | CNN por especie, dirigida por dRNA-seq | FASTA | promotor σ70 por especie | especies con dRNA-seq | validación específica de especie (según el artículo) | disponibilidad no verificada; el dominio SAPPHIRE no responde |

**Notas.**

- **BPROM no tiene artículo metodológico propio.** Se describe en el capítulo
  de Solovyev & Salamov (2011) y en la documentación de Softberry. Es el
  ejemplo más claro de herramienta ubicua sin validación primaria publicada;
  el benchmark de mSystems la midió y quedó última.
- **CNNProm** se entrenó con negativos de secuencias codificantes según el
  análisis de mSystems; eso favorece el sesgo AT y hay que tenerlo en cuenta
  al comparar.
- **Sigma70Pred** es el caso más limpio de «independiente» temporal: train
  RegulonDB 9.0, test RegulonDB 10.8. Sigue siendo la misma base y la misma
  lógica de curación; la independencia es de versión, no de proceso (§5.5).
- **SAPPHIRE** y **SAPPHIRE.CNN** son aportaciones distintas (2020 y 2022) y
  no deben citarse como una sola herramienta.

### 3.2 Predicción σ-específica y asignación multi-sigma (T2/T3)

| Herramienta | Método | Entrada | Salida | Especie | Validación reportada | Código / disponibilidad |
|---|---|---|---|---|---|---|
| **bTSSfinder** | modelos de promotor por clase de σ + predicción de TSS; 5 clases E. coli (σ70, σ38, σ32, σ28, σ24) y 5 cianobacterias (σA, σC, σH, σG, σF) | genoma/secuencia FASTA | TSS putativos + clase de σ | E. coli + 3 cianobacterias | MCC 0,86 / F1 0,93 frente al siguiente mejor (MCC 0,59 / F1 0,79) en sus datos; incluida en mSystems | standalone (descarga KAUST) + web; web 503, descarga no verificada |
| **BacPP** | reglas ponderadas extraídas de redes neuronales; una secuencia → promotor/no + σ | una secuencia | promotor/no promotor + σ asignado (σ24, σ28, σ32, σ38, σ54, σ70) | E. coli / enterobacterias | accuracy por σ: 86,9 / 92,8 / 91,5 / 89,3 / 97,0 / 83,6 %; 76 % en otras enterobacterias | web (2016) + versión desktop original; servidor no resuelve (2026-11-26) |
| **iPromoter-2L** | dos capas: promotor/no promotor; luego tipo σ (σ24, σ28, σ32, σ38, σ54, σ70) con PseKNC multi-ventana + SVM | FASTA | promotor/no + tipo σ | E. coli (seis clases) | validación interna en sus datos; incluida en mSystems | web `bioinformatics.hitsz.edu.cn/iPromoter-2L`; no responde (2026-11-26) |
| **iPromoter-2L2.0** | smoothing cutting window + k-mer/PseKNC + ensemble de SVMs | FASTA | promotor/no + tipo σ | E. coli (seis clases) | mejora sobre 2L según el artículo | web; disponibilidad no verificada |
| **PromoterLCNN** | CNN ligera, dos etapas, multiclase | FASTA | promotor/no + clase σ (σ70, σ24, σ32, σ38, σ28, σ54) | E. coli | comparación con clasificadores CNN; ahorro 30–90 % de tiempo sin pérdida de calidad (según el artículo) | código no verificado |
| **MULTiPly** | multicapa: k-mer, autocovarianza, bi-profile Bayes, k-NN; F-score | FASTA | promotor/no + tipo específico | no declarada en el abstract; evaluada en E. coli por mSystems | 5-fold CV, jackknife y test independiente en sus datos | web `flagshipnt.erc.monash.edu/MULTiPly` |
| **ProPr54** | CNN híbrida (convolucional + recurrente) | genoma/secuencia | promotores σ54 + regulón | 33 especies bacterianas (446 sitios validados) | entrenamiento y test en regiones intergénicas; generalización a genomas sin sitios validados (según el artículo) | web `propr54.molgenrug.nl`; activo |
| **PromoterAtlas** | transformer de 1,8 M de parámetros, 9 M de secuencias, 3371 especies | genoma GenBank | GenBank anotado: tipo de σ (σ70, σ38, σ54, σ32, σ28, σ24), cajas −10/−35 (o −12/−24), RBS, terminadores | Gammaproteobacteria | «varios niveles de validación»; embeddings que agrupan promotores por σ más que por especie (según el artículo) | código + pesos en GitHub (`LucasCoppens/PromoterAtlas`) y Zenodo; activo |

**Notas.**

- **bTSSfinder** es la referencia canónica de T3 multi-sigma en E. coli. Su
  limitación es explícita: solo E. coli y cianobacterias; no es un modelo
  general.
- **BacPP** es el clásico que asigna σ con una salida explícita por factor y
  reporta accuracy por sigma; su debilidad es que el método son reglas y que
  el servidor ya no está.
- **iPromoter-2L** y **PromoterLCNN** son los comparadores naturales de
  nuestro T3 en E. coli: seis clases y la misma tarea; PromoterLCNN declara
  datos derivados de RegulonDB. Pero sus splits son CV dentro de la misma
  distribución (§5.5).
- **ProPr54** es el único caso verificado, entre las herramientas revisadas,
  de entrenamiento multi-especie explícito para un sigma alternativo (33
  especies). Es el modelo a imitar para σ54, no para σ70.
- **PromoterAtlas** es el único modelo verificado a escala de filo
  (Gammaproteobacteria). Su propia escala implica que E. coli y muchas especies
  de evaluación probablemente estuvieron en el entrenamiento; la validación
  «cross-species» debe leerse con cuidado.

### 3.3 Escaneo de motivos y frameworks de anotación

| Herramienta | Método | Entrada | Salida | Especie | Validación reportada | Código / disponibilidad |
|---|---|---|---|---|---|---|
| **Virtual Footprint** | matching de PWMs contra la base PRODORIC | secuencia/genoma + matriz | sitios TFBS/promotores con score | según matriz (PRODORIC) | incluida en mSystems; sin accuracy única | web `prodoric.de`; sin código |
| **PePPER** | predicción de promotores + TFBS + comparación de regulones | genoma GenBank / secuencia | promotores, TFBS, regulones | cualquier bacteria; conocimiento concentrado en *L. lactis* | validación biológica en *L. lactis* (WalRK, CodY); sin benchmark cuantitativo | web; el dominio ahora sirve **ProPr** (2026-11-26) |
| **ProPr** | modelos CNN por especie/GC (p. ej. `FGC3545_*`, `Bacillus`, `Acinetobacter`), predicción en regiones intergénicas + terminadores | genoma/GFF | promotores en tabla/GFF, terminadores | especies cubiertas por sus modelos | sin artículo peer-reviewed localizado (Crossref, 2026-11-26); modelos actualizados en 2025 según el servidor | web `pepper.molgenrug.nl` (sirve ProPr); sin código |
| **SigmoID** | aplicación de escritorio que integra MEME/FIMO, navegador genómico y bases (RegPrecise, RegulonDB) | genoma + matrices | promotores, TFBS, terminadores en contexto | bacterias; validado en *Pectobacterium atrosepticum* | >900 sitios inferidos y corrección de >300 anotaciones (según el artículo) | open source, multiplataforma |

**Notas.**

- **Virtual Footprint** y **SigmoID** no son modelos entrenados: son
  infraestructura de escaneo. Su resultado depende de la matriz y del umbral;
  no tienen una «accuracy» propia comparable.
- **PePPER** fue un servidor de referencia para regulones de bacterias
  Gram-positivas, pero su Interfaz original ya no está; el dominio sirve ProPr.
  Cualquier comparación con PePPER debe declarar que se usó ProPr.
- **ProPr** es el ejemplo de herramienta viva sin publicación indexada
  verificable. Es útil como anotador externo, no como referencia metodológica.

### 3.4 Herramientas relacionadas que conviene conocer

- **iPromoter-BnCNN** (Amin et al., 2020): CNN ramificada para identificar y
  clasificar promotores sigma.
- **iPro54-PseKNC** (Lin et al., 2014): predictor binario previo de promotores
  σ54 con PseKNC.
- **iProm-Sigma54** (Shujaat et al., 2023; DOI
  [10.3390/cells12060829](https://doi.org/10.3390/cells12060829)): CNN para
  σ54; otra alternativa binaria al problema que ProPr54 aborda con más datos.
- **G4PromFinder** (Di Salvo et al., 2018): promotores en genomas GC-ricos
  usando elementos AT-ricos y G-cuádruplex; ejemplo de modelo dependiente de
  composición, no de sigma.
- **FIMO/MEME** (Grant et al., 2011): escaneo PWM estándar; línea base
  obligatoria para cualquier comparación con PWMs propias.

---

## 4. Machine learning y deep learning para promotores bacterianos

### 4.1 Qué se ha hecho

1. **Features de secuencia + SVM/ensembles (2011–2019).** BacPP (reglas de
   redes), 70ProPred, iPro70-FMWin, iPromoter-2L, MULTiPly. Funcionan bien en
   datasets pequeños y son rápidos; su punto débil es que las features de
   composición (k-mer, PseKNC, autocorrelación) pueden separar positivos de
   negativos sin aprender la geometría del promotor.
2. **CNN (2017–2023).** CNNProm, PromoterLCNN, iPromoter-BnCNN, ProPr54,
   SAPPHIRE.CNN. Aprenden motivos y su disposición; mejoran la clasificación
   en ventana. Requieren más datos y son sensibles al split.
3. **Transformer y modelos de escala (2026).** PromoterAtlas: 1,8 M de
   parámetros, 9 M de secuencias de 3371 especies de Gammaproteobacteria,
   anotación de genoma completo (promotores por σ, RBS, TFBS, terminadores).
   Es el primer modelo verificado que integra reconocimiento de elementos y
   transferencia entre especies en un solo marco.

### 4.2 Qué significa «mejor» en cada caso

- **En ventana (T1/T2):** mejor AUC/AUPRC con los mismos negativos. No
  garantiza nada a escala genómica.
- **En T3:** mejor top-1/top-3 y macro-F1 con las mismas clases y el mismo
  gold. La accuracy global está dominada por σ70.
- **En T4:** mejor precisión a recall fijo en genomas completos, con distancia
  al TSS real y sin duplicar promotores.

Los artículos suelen reportar la primera; nuestro benchmark necesita la
segunda y, si se puede, la tercera.

### 4.3 Límites que el DL no resuelve

- **Datos pequeños y sesgados.** 741 promotores σ70 (Sigma70Pred), 446 sitios
  σ54 (ProPr54), 3.390 TSS de B. subtilis (DBTBS). Un transformer puede
  superarlos en su test, pero no crea diversidad de especies.
- **Negativos.** Un modelo más profundo explota mejor cualquier sesgo de los
  negativos; el problema es el diseño del negativo, no la capacidad del
  modelo.
- **Etiqueta.** La «verdad» de sigma en RegulonDB/DBTBS es una mezcla de
  experimentos, inferencia de regulón y asignación por motivo. El modelo
  aprende esa mezcla, no la biología.
- **Interpretabilidad.** Los filtros/motivos aprendidos no siempre
  corresponden a cajas conocidas; la atribución es útil para generar
  hipótesis, no para asignar sigma.
- **Coste y reproducibilidad.** Pesos, versiones y entornos cambian; sin
  código y versión fijados, la comparación no es reproducible.

---

## 5. Evaluación: TSS, promotor, sigma, negativos, cross-species

### 5.1 Qué se mide exactamente

| Objeto | Definición operativa | Trampa frecuente |
|---|---|---|
| **TSS** | posición 1-based de inicio de transcripción, con hebra | el TSS experimental depende del método (dRNA-seq, Cappable-seq, 5′-RACE); tolerancia ±0/±5/±10 cambia el resultado |
| **Promotor** | región −80/+20 (o −50/+10, etc.) alrededor del TSS | la «región» es una convención; dos artículos con ventanas distintas no son comparables |
| **Sigma** | factor que reconoce el promotor | se asigna por motivo, por regulón o por experimento; no son equivalentes |
| **Negativo** | ventana no promotora | su definición (codificante, aleatoria, intergénica, otro sigma) determina el resultado |
| **Operón** | unidad transcripcional con uno o varios TSS | agregar a una fila por operón oculta multi-sigma y TSS alternativos |

En nuestro benchmark (`docs/SIGMA_BENCHMARK.md`): ventana de 101 nt (−80/+20),
TSS en el índice 80, seis PWMs, nulo genómico de 500 ventanas aleatorias,
top-1/top-3, comparación por operón con `exact` y `dominant`. Es una elección
defendible para T3, pero debe declararse que el nulo es genómico (incluye
promotores reales) y que el z-score no es un p-valor biológico.

### 5.2 Negativos: seis diseños y qué premian

| Tipo de negativo | Qué tests | Sesgo que introduce |
|---|---|---|
| **Secuencias codificantes** | discriminación ORF vs promotor | composición y estructura de ORF; el modelo aprende a detectar código, no promotores |
| **Aleatorias con composición similar** | discriminación de motivo con composición controlada | premia a los modelos que explotan AT-riqueza; es el negativo del benchmark mSystems |
| **Genómicas aleatorias sin control** | detección a escala | incluye promotores reales no anotados; label noise |
| **Intergénicas no-TSS** | localización realista | puede contener TSS no anotados y promotores silentes |
| **TSS de otro sigma** | T3/T2 «duro» | es el negativo correcto para asignación; casi nunca se reporta |
| **Motivos barajados / dinucleótidos preservados** | control de motivo vs composición | demasiado fácil; no representa el genoma |

Evidencia:

- El benchmark de mSystems (Cassiano & Silva-Rocha, 2020) usó promotores
  validados de E. coli y un control de secuencias aleatorias con distribución
  de nucleótidos similar; encontró que varias herramientas están sesgadas
  hacia regiones AT-ricas y que BPROM fue la peor del panel, mientras
  CNNProm, iPro70-FMWin, 70ProPred e iPromoter-2L rindieron alto.
- El mismo trabajo señala que CNNProm e iPromoter-2L se entrenaron con
  secuencias codificantes como negativos, lo que explica parte de su
  comportamiento en secuencias aleatorias.
- Zhang et al. (2022) curaron 58 datasets de 7 especies y compararon 19
  predictores con servidor o herramienta local; concluyen que DL y ML
  superan a los métodos de scoring, pero que la comparación solo es válida
  dentro del mismo dataset y la misma tarea.

### 5.3 Métricas: qué reportar y qué no

**No usar accuracy sola cuando σ70 es mayoría.** En nuestro gold, la mayoría
trivial es 62,6 % (K-12) y 87,6 % (B. subtilis). Un modelo con 57 % o 44,8 %
de top-1 está por debajo de no hacer nada.

Reportar como mínimo:

- **Top-1 y top-3** (y rango medio de la verdad) para T3.
- **Macro-F1** y **balanced accuracy** para no premiar la clase mayoritaria.
- **MCC** (útil con desbalance; el benchmark mSystems lo usa).
- **Precisión/recall por sigma**, con support; reportar σ54 y ECF aparte.
- **AUPRC** si hay scores; **AUROC solo como complemento** (puede ser
  optimista con negativos fáciles).
- **Calibración** (ECE, Brier, reliability bins) solo si el modelo emite
  probabilidades; un z-score no es una probabilidad.
- **Precisión a recall fijo a escala genómica** (p. ej. recall 50 %, 80 %),
  y número de falsos positivos por Mb.
- **Cobertura**: proporción de TSS etiquetados con PWM disponible; en
  B. subtilis, 27,6 % de los TSS etiquetados no tienen PWM equivalente y
  quedan fuera del denominador. Reportarlo evita inflar la accuracy.
- **Intervalos de confianza** por bootstrap agrupado por locus/especie, no
  por ventana.
- **Multi-etiqueta**: subset accuracy/Hamming si se permiten varios sigma.

### 5.4 Cross-species: la evidencia disponible

| Evidencia | Qué muestra |
|---|---|
| SAPPHIRE (2020) | las herramientas históricas, centradas en E. coli, rinden por debajo en *Pseudomonas*; un modelo entrenado en TSS experimentales de *P. aeruginosa*/*P. putida* mejora |
| bTSSfinder (2017) | solo E. coli y cianobacterias; no es un modelo general |
| ProPr54 (2025) | 446 sitios de 33 especies; entrenamiento y test en regiones intergénicas; generaliza a genomas sin sitios validados (según el artículo) |
| PromoterAtlas (2026) | 3371 especies de un solo filo (Gammaproteobacteria); los embeddings agrupan por σ más que por especie |
| Nuestro benchmark (interno) | PWMs de E. coli transferidas a B. subtilis: top-1 44,8 % frente a mayoría 87,6 %; K-12 57,0 % frente a 62,6 % |

Lectura: **la transferencia entre géneros no está demostrada** para σ70 en
general; para σ54 hay un caso multi-especie creíble (ProPr54); para
Gammaproteobacteria hay un modelo a gran escala (PromoterAtlas) pero sin
held-out de filo. Cualquier afirmación de generalidad debe especificar la
distancia taxonómica y el número de especies held-out.

### 5.5 Fugas de datos: mecanismos y cómo detectarlas

| Mecanismo | Cómo se manifiesta | Cómo detectarlo / mitigarlo |
|---|---|---|
| **Ventanas solapadas** | TSS vecinos comparten 101 nt; el mismo promotor cae en train y test | agrupar por locus/región (GroupKFold 10–100 kb); deduplicar por identidad ≥90 % |
| **Misma base, distinta versión** | train RegulonDB 9.0, test 10.8; la mayoría de entradas persisten | split temporal real (fecha de curación) o holdout por locus, no por versión |
| **PWMs estimadas en el gold** | la matriz se construye con los mismos promotores que se evalúan | estimar PWMs solo en el fold de entrenamiento; dejar el fold de test fuera |
| **Etiqueta circular** | el sigma del gold se infirió por motivo/regulón, no por experimento | auditar provenance; separar «TSS experimental» de «asignación inferida» |
| **Negativos fáciles** | el modelo aprende composición, no motivo | negativos duros (otro sigma, TSS no usado), negativos intergénicos, controles de composición |
| **Reutilización de benchmark** | todos entrenan y testean en RegulonDB; el test deja de ser independiente | held-out de especie/filo; congelar una partición no usada por nadie |
| **Optimización en test** | elegir política/umbral mirando el gold completo | CV anidada; congelar hiperparámetros antes de mirar el test |
| **Exposición en preentrenamiento** | un modelo de 3371 especies probablemente vio el genoma de test | held-out por especie/tiempo; no evaluar en genomas del entrenamiento |

Evidencia interna de este riesgo: `docs/SIGMA_STRATEGIES.md` documenta que el
mejor resultado de K-12 (70,6 % top-1, modelo jerárquico) usa entrenamiento
sobre el mismo gold y PWMs de RegulonDB, y por tanto **no es generalización
independiente**. Nuestro `SIGMA_BENCHMARK.md` también usa PWMs de RegulonDB y
gold de RegulonDB; la evaluación es útil para diagnóstico, no para afirmar
generalidad.

---

## 6. Qué métodos son comparables

La comparabilidad no es una propiedad de la herramienta, sino de la terna
**(tarea, datos, métrica)**. Esta es la clasificación que propongo para
nuestro trabajo.

| Nivel | Pregunta | Herramientas comparables | Requisito para que la comparación sea válida |
|---|---|---|---|
| **T1** | promotor vs no promotor | BPROM, CNNProm, Sigma70Pred, 70ProPred, iPro70-FMWin, SAPPHIRE, SAPPHIRE.CNN | mismos negativos, misma ventana, misma especie, mismo umbral o misma curva |
| **T2** | promotor de σX vs no promotor | ProPr54 (σ54), iPro70-FMWin (σ70), Sigma70Pred (σ70), BacPP por sigma | mismo σ, mismos negativos, misma definición de positivo |
| **T3** | asignación de σ en TSS | bTSSfinder, BacPP, iPromoter-2L/2L2.0, PromoterLCNN, MULTiPly | mismo gold multi-clase, mismas clases, misma cobertura, mismo tratamiento de «sin PWM» |
| **T4** | anotación genómica | PromoterAtlas, ProPr, SigmoID, PePPER (histórico) | mismos genomas, misma tolerancia de TSS, misma métrica de precisión/recall |

**Comparaciones que sí tienen sentido hoy:**

1. **Nuestro T3 en K-12** contra bTSSfinder (5 clases), iPromoter-2L (6),
   PromoterLCNN (6) y BacPP (6), re-ejecutados sobre las mismas ventanas
   −80/+20 y el mismo gold, con y sin clase «desconocido».
2. **Nuestro T2 σ54** contra ProPr54, sobre genomas con TSS σ54 validados y
   sobre genomas sin ellos (transferencia).
3. **Nuestro T1** contra CNNProm/Sigma70Pred usando exactamente los mismos
   negativos (los cinco sets biológicos de §5.2), no los negativos de cada
   artículo.
4. **Nuestro T4** contra PromoterAtlas en genomas held-out de
   Gammaproteobacteria, midiendo precisión a recall fijo y distancia al TSS.

**Comparaciones que no tienen sentido:**

- Rankear BPROM (T1, σ70) contra bTSSfinder (T3, 5 clases) por accuracy.
- Comparar un CNN entrenado con negativos codificantes contra un PWM escaneado
  contra genoma completo.
- Usar cifras de artículos distintos como si fueran un leaderboard.
- Comparar modelos eucariotas con bacterianos.
- Comparar un modelo de fuerza/expresión con uno de detección.

---

## 7. Qué no se puede concluir

1. **No se puede concluir un ranking de herramientas a partir de las cifras
   publicadas.** Cada artículo usa su dataset, su ventana y sus negativos.
2. **No se puede concluir que una herramienta «asigna sigma» si es binaria.**
   Sigma70Pred, 70ProPred, iPro70-FMWin, CNNProm, SAPPHIRE y BPROM son
   detectores σ70; no asignan los otros sigma.
3. **No se puede concluir generalidad entre especies** sin un held-out de
   especie o filo. La evidencia actual es favorable solo para casos concretos
   (ProPr54 en σ54; PromoterAtlas dentro de Gammaproteobacteria).
4. **No se puede concluir que un hit de PWM sea un promotor activo in vivo.**
   El motivo es necesario pero no suficiente; la accesibilidad, los
   activadores y la competencia por el core deciden.
5. **No se puede concluir que «no hay hit» sea «no hay promotor».** Los
   promotores sin caja canónica, los dependientes de activador y los ECF
   escapan a las PWMs.
6. **No se puede concluir dirección, magnitud ni fuerza de expresión** a
   partir de la detección o de la asignación de sigma. Eso requiere un modelo
   distinto (fuerza de promotor, TRS, datos de expresión).
7. **No se puede concluir que la accuracy en ventanas prediga la precisión
   genómica.** El número de ventanas y la prevalencia cambian drásticamente
   la precisión.
8. **No se puede concluir que los negativos de un artículo representen el
   genoma.** Los negativos aleatorios/codificantes son los más cuestionables.
9. **No se puede concluir que un modelo sea temporalmente robusto** porque
   use una versión posterior de la misma base. La curación cambia, pero el
   proceso y los promotores se solapan.
10. **No se puede concluir superioridad general del DL.** Supera a scoring en
    varios benchmarks, pero no en todos los regímenes de datos ni en todas
    las especies, y no resuelve el sesgo de negativos.
11. **No se puede concluir que un modelo σ70 sirva para σ54 o ECF.** La
    arquitectura y los motivos son distintos; ProPr54 existe precisamente
    porque no hay transferencia.
12. **No se puede llamar «confirmada» a una asignación computacional.** Como
    máximo, «hipótesis con evidencia de motivo/contexto». La confirmación
    requiere experimento (p. ej. dRNA-seq, mutante de sigma, ChIP-seq).

---

## 8. Recomendaciones para nuestro benchmark

> Contexto: el módulo `sigma` hace scoring PWM por sigma (T2) y decide por
> `argmax`/`two_stage`/`hybrid` (T3). El benchmark actual
> (`docs/SIGMA_BENCHMARK.md`, `docs/SIGMA_LAYERED_EVAL.md`) evalúa T3 a nivel
> TSS con 101 nt y gold de RegulonDB/DBTBS. Las recomendaciones buscan que la
> comparación con la literatura sea defendible y que las conclusiones no
> excedan lo medido.

### P0 — imprescindible antes de comparar con herramientas externas

1. **Contrato de tarea explícito.** Cada tabla de resultados declara:
   tarea (T1/T2/T3/T4), unidad (ventana/TSS/genoma/operón), definición de
   positivo, definición de negativo, especie(s), versión del gold, ventana,
   umbral y política de decisión. Sin eso, no se mezclan filas.
2. **Split sin fuga.**
   - GroupKFold por locus/región de 10–100 kb (ya usado en
     `SIGMA_STRATEGIES.md`; conservarlo).
   - Held-out de especie completa (B. subtilis fuera de K-12 y viceversa).
   - Split temporal por versión/fecha de curación del gold, además del split
     por locus.
   - Deduplicar promotores con identidad ≥90 % en la ventana antes de
     partir; si no, el mismo promotor aparece en train y test.
   - Estimar las PWMs **solo** en el fold de entrenamiento; el fold de test
     no puede contribuir a la matriz.
3. **Cinco conjuntos de negativos biológicos, reportados por separado.**
   Codificantes, aleatorios con composición similar, genómicos aleatorios,
   intergénicos no-TSS y TSS de otro sigma; el control de motivos barajados se
   reporta aparte, como sanity check. El set «otro sigma» es el único que mide
   asignación de verdad; los otros miden detección.
4. **Líneas base obligatorias.** Mayoría σ70; PWM argmax (actual); PWM con
   matrices barajadas; k-mer/logística; composición GC/AT sola; detector
   binario promotor/no (para separar T1 de T3). Si el modelo no supera a
   «composición sola», no hay señal de motivo.
5. **Métricas mínimas.** Top-1, top-3, macro-F1, balanced accuracy, MCC,
   precisión/recall por sigma, cobertura (`no_pwm`), y precisión genómica a
   recall fijo. Intervalos por bootstrap agrupado por locus/especie.
   Calibración solo si se emiten probabilidades; un z-score no lo es.

### P1 — para que el benchmark sea útil fuera de nuestro repo

6. **Protocolo cross-species.** Evaluar K-12 → B. subtilis y B. subtilis →
   K-12; añadir un holdout de otro género (Pseudomonas o cianobacteria) sin
   reentrenar. Comparar PWM transferida vs PWM específica de especie.
   Para σ54, usar ProPr54 como referencia externa.
7. **Multi-sigma y multi-etiqueta.** Reportar los TSS con ≥2 sigma (124 en
   K-12) con subset accuracy/Hamming y permitir `ambiguous`/`unresolved` en
   la salida. No forzar una única etiqueta.
8. **Nivel operón con agregación explícita.** Reportar `exact`, `dominant` y
   `multi-sigma` por separado (ya existe en el benchmark); no usar una fila
   por operón para estimar T3.
9. **Comparación con herramientas externas reproducible.** Congelar
   versión/fecha/URL; usar las mismas ventanas FASTA; probar defaults y un
   umbral ajustado solo en train; registrar fallos y límites de envío
   (herramientas web con lotes pequeños → muestrear y declarar).

### P2 — calidad y mantenimiento

10. **Auditoría de provenance del gold.** Separar TSS experimentales de
    asignaciones inferidas por regulón/motivo; guardar el campo de
    provenance (el módulo ya lo tiene) y reportar métricas por subconjunto.
11. **Controles de sanidad.** Barajar etiquetas (debe caer a la mayoría),
    barajar dinucleótidos, intercambiar negativos, holdout de genoma
    completo. Si el modelo sigue rindiendo, hay fuga.
12. **Release versionado del benchmark.** Checksums de gold/PWMs/ventanas,
    entorno y comandos; una partición «sellada» que no se use para
    desarrollo.

### Receta mínima viable (CPU, sin red)

- Gold: K-12 4.057 TSS (3.315 evaluables con PWM) y B. subtilis 3.390
  (2.226 evaluables), con el mapeo de etiquetas ya documentado.
- Cinco folds GroupKFold (10–100 kb) + un holdout de especie.
- Cinco sets de negativos; el set «otro sigma» es obligatorio.
- Baselines: mayoría, PWM argmax, PWM barajada, k-mer, composición.
- Métricas: top-1/top-3, macro-F1, balanced accuracy, MCC, P/R por sigma,
  cobertura, y precisión a recall fijo.
- Salida: una tabla por tarea y por set de negativos; nunca una sola cifra
  agregada.

---

## 9. Referencias verificadas

> Todas las entradas tienen DOI resuelto en Crossref el 2026-11-26. El PMID
> se incluye solo cuando se recuperó por E-utilities. Los estados de los
> servidores están en el Apéndice A.

### 9.1 Biología, arquitectura y motivos

- **Barrios, H.; Valderrama, B.; Morett, E. (1999).** Compilation and analysis
  of σ54-dependent promoter sequences. *Nucleic Acids Research* 27(22):
  4305–4313. DOI: [10.1093/nar/27.22.4305](https://doi.org/10.1093/nar/27.22.4305).
- **Feklistov, A.; Sharon, B. D.; Darst, S. A.; Gross, C. A. (2014).**
  Bacterial sigma factors: a historical, structural, and genomic perspective.
  *Annual Review of Microbiology* 68: 357–376. DOI:
  [10.1146/annurev-micro-092412-155737](https://doi.org/10.1146/annurev-micro-092412-155737).
- **Gruber, T. M.; Gross, C. A. (2003).** Multiple sigma subunits and the
  partitioning of bacterial transcription space. *Annual Review of
  Microbiology* 57: 441–466. DOI:
  [10.1146/annurev.micro.57.030502.090913](https://doi.org/10.1146/annurev.micro.57.030502.090913).
- **Helmann, J. D. (2019).** Where to begin? Sigma factors and the selectivity
  of transcription initiation in bacteria. *Molecular Microbiology* 112(2):
  335–347. DOI: [10.1111/mmi.14309](https://doi.org/10.1111/mmi.14309).
- **Paget, M. S. (2015).** Bacterial sigma factors and anti-sigma factors:
  structure, function and distribution. *Biomolecules* 5(3): 1245–1265. DOI:
  [10.3390/biom5031245](https://doi.org/10.3390/biom5031245).
- **Shultzaberger, R. K.; Chen, Z.; Lewis, K. A.; Schneider, T. D. (2007).**
  Anatomy of *Escherichia coli* σ70 promoters. *Nucleic Acids Research* 35(3):
  771–788. DOI: [10.1093/nar/gkl956](https://doi.org/10.1093/nar/gkl956).
- **Stormo, G. D. (2000).** DNA binding sites: representation and discovery.
  *Bioinformatics* 16(1): 16–23. DOI:
  [10.1093/bioinformatics/16.1.16](https://doi.org/10.1093/bioinformatics/16.1.16).
- **Grant, C. E.; Bailey, T. L.; Noble, W. S. (2011).** FIMO: scanning for
  occurrences of a given motif. *Bioinformatics* 27(7): 1017–1018. DOI:
  [10.1093/bioinformatics/btr064](https://doi.org/10.1093/bioinformatics/btr064).

### 9.2 Herramientas

- **Shahmuradov, I. A.; Mohamad Razali, R.; Bougouffa, S.; Radovanovic, A.;
  Bajic, V. B. (2017).** bTSSfinder: a novel tool for the prediction of
  promoters in cyanobacteria and *Escherichia coli*. *Bioinformatics* 33(3):
  334–340. DOI: [10.1093/bioinformatics/btw629](https://doi.org/10.1093/bioinformatics/btw629).
  PMID: 27694198.
- **de Avila e Silva, S.; Echeverrigaray, S.; Gerhardt, G. J. L. (2011).**
  BacPP: bacterial promoter prediction — a tool for accurate sigma-factor
  specific assignment in enterobacteria. *Journal of Theoretical Biology*
  287: 92–99. DOI:
  [10.1016/j.jtbi.2011.07.017](https://doi.org/10.1016/j.jtbi.2011.07.017).
  PMID: 21827769.
- **de Avila e Silva, S.; Notari, D. L.; Neis, F. A.; Ribeiro, H. G.;
  Echeverrigaray, S. (2016).** BacPP: a web-based tool for Gram-negative
  bacterial promoter prediction. *Genetics and Molecular Research* 15(2):
  gmr.15027973. DOI:
  [10.4238/gmr.15027973](https://doi.org/10.4238/gmr.15027973).
  PMID: 27173187.
- **de Jong, A.; Pietersma, H.; Cordes, M.; Kuipers, O. P.; Kok, J. (2012).**
  PePPER: a webserver for prediction of prokaryote promoter elements and
  regulons. *BMC Genomics* 13: 299. DOI:
  [10.1186/1471-2164-13-299](https://doi.org/10.1186/1471-2164-13-299).
- **Coppens, L.; Lavigne, R. (2020).** SAPPHIRE: a neural network based
  classifier for σ70 promoter prediction in *Pseudomonas*. *BMC
  Bioinformatics* 21: 415. DOI:
  [10.1186/s12859-020-03730-z](https://doi.org/10.1186/s12859-020-03730-z).
  PMID: 32962628.
- **Coppens, L.; Wicke, L.; Lavigne, R. (2022).** SAPPHIRE.CNN: implementation
  of dRNA-seq-driven, species-specific promoter prediction using convolutional
  neural networks. *Computational and Structural Biotechnology Journal* 20:
  4969–4974. DOI:
  [10.1016/j.csbj.2022.09.006](https://doi.org/10.1016/j.csbj.2022.09.006).
- **Patiyal, S.; Singh, N.; Ali, M. Z.; Pundir, D. S.; Raghava, G. P. S.
  (2022).** Sigma70Pred: a highly accurate method for predicting sigma70
  promoter in *Escherichia coli* K-12 strains. *Frontiers in Microbiology*
  13: 1042127. DOI:
  [10.3389/fmicb.2022.1042127](https://doi.org/10.3389/fmicb.2022.1042127).
  PMID: 36452927.
- **Umarov, R. Kh.; Solovyev, V. V. (2017).** Recognition of prokaryotic and
  eukaryotic promoters using convolutional deep learning neural networks.
  *PLoS ONE* 12(2): e0171410. DOI:
  [10.1371/journal.pone.0171410](https://doi.org/10.1371/journal.pone.0171410).
  PMID: 28158264.
- **He, W.; Jia, C.; Duan, Y.; Zou, Q. (2018).** 70ProPred: a predictor for
  discovering sigma70 promoters based on combining multiple features. *BMC
  Systems Biology* 12(Suppl 4): 44. DOI:
  [10.1186/s12918-018-0570-1](https://doi.org/10.1186/s12918-018-0570-1).
  PMID: 29745856.
- **Rahman, M. S.; Aktar, U.; Jani, M. R.; Shatabda, S. (2019).**
  iPro70-FMWin: identifying Sigma70 promoters using multiple windowing and
  minimal features. *Molecular Genetics and Genomics* 294: 69–84. DOI:
  [10.1007/s00438-018-1487-5](https://doi.org/10.1007/s00438-018-1487-5).
  PMID: 30187132.
- **Liu, B.; Yang, F.; Huang, D.-S.; Chou, K.-C. (2018).** iPromoter-2L: a
  two-layer predictor for identifying promoters and their types by
  multi-window-based PseKNC. *Bioinformatics* 34(1): 33–40. DOI:
  [10.1093/bioinformatics/btx579](https://doi.org/10.1093/bioinformatics/btx579).
  PMID: 28968797.
- **Liu, B.; Li, K. (2019).** iPromoter-2L2.0: identifying promoters and their
  types by combining smoothing cutting window algorithm and sequence-based
  features. *Molecular Therapy — Nucleic Acids* 18: 80–87. DOI:
  [10.1016/j.omtn.2019.08.008](https://doi.org/10.1016/j.omtn.2019.08.008).
  PMID: 31536883.
- **Hernández, D.; Jara, N.; Araya, M.; Durán, R. E.; Buil-Aranda, C. (2022).**
  PromoterLCNN: a light CNN-based promoter prediction and classification
  model. *Genes* 13(7): 1126. DOI:
  [10.3390/genes13071126](https://doi.org/10.3390/genes13071126).
  PMID: 35885909.
- **Zhang, M.; Li, F.; Marquez-Lago, T. T.; Leier, A.; Fan, C.; Kwoh, C. K.;
  Chou, K.-C.; Song, J.; Jia, C. (2019).** MULTiPly: a novel multi-layer
  predictor for discovering general and specific types of promoters.
  *Bioinformatics* 35(17): 2957–2965. DOI:
  [10.1093/bioinformatics/btz016](https://doi.org/10.1093/bioinformatics/btz016).
  PMID: 30649179.
- **Achterberg, T.; de Jong, A. (2025).** ProPr54 web server: predicting σ54
  promoters and regulon with a hybrid convolutional and recurrent deep neural
  network. *NAR Genomics and Bioinformatics* 7(1): lqae188. DOI:
  [10.1093/nargab/lqae188](https://doi.org/10.1093/nargab/lqae188).
  PMID: 39781509.
- **Coppens, L.; Ledesma-Amaro, R. (2026).** PromoterAtlas: decoding regulatory
  sequences across Gammaproteobacteria using a transformer model. *Nature
  Communications* 17: 6451. DOI:
  [10.1038/s41467-026-72837-3](https://doi.org/10.1038/s41467-026-72837-3).
- **Münch, R.; Hiller, K.; Grote, A.; Scheer, M.; Klein, J.; Schobert, M.;
  Jahn, D. (2005).** Virtual Footprint and PRODORIC: an integrative framework
  for regulon prediction in prokaryotes. *Bioinformatics* 21(22): 4187–4189.
  DOI: [10.1093/bioinformatics/bti635](https://doi.org/10.1093/bioinformatics/bti635).
- **Nikolaichik, Y.; Damienikan, A. U. (2016).** SigmoID: a user-friendly tool
  for improving bacterial genome annotation through analysis of transcription
  control signals. *PeerJ* 4: e2056. DOI:
  [10.7717/peerj.2056](https://doi.org/10.7717/peerj.2056).
- **Amin, R.; Rahman, C. R.; Ahmed, S.; Sifat, M. H. R.; Liton, M. N. K.;
  Rahman, M. M.; Khan, M. Z. H.; Shatabda, S. (2020).** iPromoter-BnCNN: a
  novel branched CNN-based predictor for identifying and classifying sigma
  promoters. *Bioinformatics* 36(19): 4869–4875. DOI:
  [10.1093/bioinformatics/btaa609](https://doi.org/10.1093/bioinformatics/btaa609).
- **Lin, H.; Deng, E.-Z.; Ding, H.; Chen, W.; Chou, K.-C. (2014).**
  iPro54-PseKNC: a sequence-based predictor for identifying sigma-54 promoters
  in prokaryote with pseudo k-tuple nucleotide composition. *Nucleic Acids
  Research* 42(21): 12961–12972. DOI:
  [10.1093/nar/gku1019](https://doi.org/10.1093/nar/gku1019).
- **Shujaat, M.; Kim, H.; Tayara, H.; Chong, K. T. (2023).** iProm-Sigma54: a
  CNN base prediction tool for σ54 promoters. *Cells* 12(6): 829. DOI:
  [10.3390/cells12060829](https://doi.org/10.3390/cells12060829).
- **Di Salvo, M.; Pinatel, E.; Talà, A.; Fondi, M.; Peano, C.; Alifano, P.
  (2018).** G4PromFinder: an algorithm for predicting transcription promoters
  in GC-rich bacterial genomes based on AT-rich elements and G-quadruplex
  motifs. *BMC Bioinformatics* 19: 36. DOI:
  [10.1186/s12859-018-2049-x](https://doi.org/10.1186/s12859-018-2049-x).
- **Solovyev, V.; Salamov, A. (2011).** Automatic annotation of microbial
  genomes and metagenomic sequences. En: Li, R. W. (ed.), *Metagenomics and
  its Applications in Agriculture, Biomedicine and Environmental Studies*.
  Nova Science Publishers, Hauppauge, NY, pp. 61–78. Sin DOI. Descripción de
  BPROM citada por la documentación de Softberry y por bTSSfinder.

### 9.3 Evaluación, benchmarks y datos

- **Cassiano, M. H. A.; Silva-Rocha, R. (2020).** Benchmarking bacterial
  promoter prediction tools: potentialities and limitations. *mSystems* 5(4):
  e00439-20. DOI:
  [10.1128/mSystems.00439-20](https://doi.org/10.1128/mSystems.00439-20).
  PMID: 32843538.
- **Zhang, M.; Jia, C.; Li, F.; Li, C.; Zhu, Y.; Akutsu, T.; Webb, G. I.;
  Zou, Q.; Coin, L. J. M.; Song, J. (2022).** Critical assessment of
  computational tools for prokaryotic and eukaryotic promoter prediction.
  *Briefings in Bioinformatics* 23(1): bbab551. DOI:
  [10.1093/bib/bbab551](https://doi.org/10.1093/bib/bbab551).
- **Sharma, C. M.; Hoffmann, S.; Darfeuille, F.; Reignier, J.; Findeiß, S.;
  Sittka, A.; Chabas, S.; Reiche, K.; Hackermüller, J.; Reinhardt, R.;
  Stadler, P. F.; Vogel, J. (2010).** The primary transcriptome of the major
  human pathogen *Helicobacter pylori*. *Nature* 464: 250–255. DOI:
  [10.1038/nature08756](https://doi.org/10.1038/nature08756).
- **Ettwiller, L.; Buswell, J.; Yigit, E.; Schildkraut, I. (2016).** A novel
  enrichment strategy reveals unprecedented number of novel transcription
  start sites at single base resolution in a model prokaryote and the gut
  microbiome. *BMC Genomics* 17: 199. DOI:
  [10.1186/s12864-016-2539-z](https://doi.org/10.1186/s12864-016-2539-z).
- **Tierrafría, V. H.; Rioualen, C.; Salgado, H.; et al. (2022).** RegulonDB
  11.0: comprehensive high-throughput datasets on transcriptional regulation
  in *Escherichia coli* K-12. *Microbial Genomics* 8: mgen000833. DOI:
  [10.1099/mgen.0.000833](https://doi.org/10.1099/mgen.0.000833).
- **Sierro, N.; Makita, Y.; de Hoon, M.; Nakai, K. (2008).** DBTBS: a database
  of transcriptional regulation in *Bacillus subtilis* containing upstream
  intergenic conservation information. *Nucleic Acids Research* 36(Database
  issue): D93–D96. DOI:
  [10.1093/nar/gkm910](https://doi.org/10.1093/nar/gkm910).
- **Kılıç, S.; White, E. R.; Sagitova, D. M.; Cornish, J. P.; Erill, I.
  (2014).** CollecTF: a database of experimentally validated transcription
  factor-binding sites in Bacteria. *Nucleic Acids Research* 42(D1):
  D156–D160. DOI:
  [10.1093/nar/gkt1123](https://doi.org/10.1093/nar/gkt1123).

### 9.4 Documentos internos citados (no son literatura)

- [`docs/SIGMA_BENCHMARK.md`](../SIGMA_BENCHMARK.md) — benchmark TSS-level,
  métricas y limitaciones.
- [`docs/SIGMA_LAYERED_EVAL.md`](../SIGMA_LAYERED_EVAL.md) — evaluación por
  capas y gold curado.
- [`docs/SIGMA_STRATEGIES.md`](../SIGMA_STRATEGIES.md) — estrategias y
  resultados GroupKFold (incluye la advertencia de no generalización).
- [`docs/SIGMA_400_PLAN.md`](../SIGMA_400_PLAN.md) — regla de no transferir
  calibradores entre especies.

---

## Apéndice A. Registro de verificación (2026-11-26)

### A.1 DOI resueltos en Crossref

| DOI | Título / revista (abreviado) | Año |
|---|---|---|
| 10.1093/bioinformatics/btw629 | bTSSfinder / Bioinformatics | 2017 |
| 10.1016/j.jtbi.2011.07.017 | BacPP (original) / J Theor Biol | 2011 |
| 10.4238/gmr.15027973 | BacPP (web) / Genet Mol Res | 2016 |
| 10.1186/1471-2164-13-299 | PePPER / BMC Genomics | 2012 |
| 10.1186/s12859-020-03730-z | SAPPHIRE / BMC Bioinformatics | 2020 |
| 10.1016/j.csbj.2022.09.006 | SAPPHIRE.CNN / Comput Struct Biotechnol J | 2022 |
| 10.3389/fmicb.2022.1042127 | Sigma70Pred / Front Microbiol | 2022 |
| 10.1371/journal.pone.0171410 | CNNProm / PLoS ONE | 2017 |
| 10.1093/bioinformatics/btx579 | iPromoter-2L / Bioinformatics | 2018 |
| 10.1016/j.omtn.2019.08.008 | iPromoter-2L2.0 / Mol Ther Nucleic Acids | 2019 |
| 10.3390/genes13071126 | PromoterLCNN / Genes | 2022 |
| 10.1093/nargab/lqae188 | ProPr54 / NAR Genom Bioinform | 2025 |
| 10.1038/s41467-026-72837-3 | PromoterAtlas / Nat Commun | 2026 |
| 10.1186/s12918-018-0570-1 | 70ProPred / BMC Syst Biol | 2018 |
| 10.1007/s00438-018-1487-5 | iPro70-FMWin / Mol Genet Genomics | 2019 |
| 10.1093/bioinformatics/btz016 | MULTiPly / Bioinformatics | 2019 |
| 10.1093/bioinformatics/bti635 | Virtual Footprint / Bioinformatics | 2005 |
| 10.7717/peerj.2056 | SigmoID / PeerJ | 2016 |
| 10.1093/bioinformatics/btaa609 | iPromoter-BnCNN / Bioinformatics | 2020 |
| 10.1093/nar/gku1019 | iPro54-PseKNC / Nucleic Acids Res | 2014 |
| 10.3390/cells12060829 | iProm-Sigma54 / Cells | 2023 |
| 10.1186/s12859-018-2049-x | G4PromFinder / BMC Bioinformatics | 2018 |
| 10.1128/mSystems.00439-20 | Benchmark mSystems | 2020 |
| 10.1093/bib/bbab551 | Critical assessment / Brief Bioinform | 2022 |
| 10.1038/nature08756 | dRNA-seq *H. pylori* / Nature | 2010 |
| 10.1186/s12864-016-2539-z | Cappable-seq / BMC Genomics | 2016 |
| 10.1099/mgen.0.000833 | RegulonDB 11.0 / Microb Genom | 2022 |
| 10.1093/nar/gkm910 | DBTBS / Nucleic Acids Res | 2008 |
| 10.1093/nar/gkt1123 | CollecTF / Nucleic Acids Res | 2014 |
| 10.1093/nar/27.22.4305 | Barrios σ54 / Nucleic Acids Res | 1999 |
| 10.1146/annurev-micro-092412-155737 | Feklistov / Annu Rev Microbiol | 2014 |
| 10.1146/annurev.micro.57.030502.090913 | Gruber & Gross / Annu Rev Microbiol | 2003 |
| 10.1111/mmi.14309 | Helmann / Mol Microbiol | 2019 |
| 10.3390/biom5031245 | Paget / Biomolecules | 2015 |
| 10.1093/nar/gkl956 | Shultzaberger / Nucleic Acids Res | 2007 |
| 10.1093/bioinformatics/16.1.16 | Stormo / Bioinformatics | 2000 |
| 10.1093/bioinformatics/btr064 | FIMO / Bioinformatics | 2011 |

### A.2 PMID recuperados por E-utilities

| PMID | Artículo |
|---|---|
| 27694198 | bTSSfinder |
| 21827769 | BacPP (2011) |
| 27173187 | BacPP (web, 2016) |
| 32962628 | SAPPHIRE |
| 36452927 | Sigma70Pred |
| 28158264 | CNNProm |
| 28968797 | iPromoter-2L |
| 31536883 | iPromoter-2L2.0 |
| 35885909 | PromoterLCNN |
| 39781509 | ProPr54 |
| 29745856 | 70ProPred |
| 30187132 | iPro70-FMWin |
| 30649179 | MULTiPly |
| 32843538 | Benchmark mSystems |

### A.3 Estado de servidores y código (HTTP, 2026-11-26)

| Recurso | URL | Estado |
|---|---|---|
| bTSSfinder (web) | `https://www.cbrc.kaust.edu.sa/btssfinder/` | **503** |
| bTSSfinder (descarga) | `http://www.cbrc.kaust.edu.sa/download/` | no verificada (mismo host) |
| BacPP | `http://www.bacpp.bioinfoucs.com/home` | **sin respuesta** |
| PePPER / ProPr | `http://pepper.molgenrug.nl/` | **200**, sirve «ProPr – Promoter Prediction Prokaryotes» |
| ProPr54 | `http://propr54.molgenrug.nl/` | **200** |
| Sigma70Pred (web) | `https://webs.iiitd.edu.in/raghava/sigma70pred/` | **200** |
| Sigma70Pred (código) | `https://github.com/raghavagps/sigma70pred` | **200** |
| iPromoter-2L | `http://bioinformatics.hitsz.edu.cn/iPromoter-2L/` | **sin respuesta** |
| SAPPHIRE | `https://www.biosapphire.com/` | **sin respuesta** |
| PromoterAtlas (código) | `https://github.com/LucasCoppens/PromoterAtlas` | **200** |
| Softberry (BPROM/CNNProm) | `http://www.softberry.com/` | no alcanzable desde este entorno (TLS/red) |

> **Advertencia de reproducibilidad.** Los servidores web cambian sin aviso
> (PePPER → ProPr es el caso documentado aquí). Cualquier comparación con una
> herramienta web debe registrar fecha, URL, versión y salida cruda.

---

## Apéndice B. Ambigüedades y erratas detectadas

1. **DOI de bTSSfinder.** Es `10.1093/bioinformatics/btw629`, no `btw620`.
   Verificado en PubMed/Crossref.
2. **«SAPPHIRE» no es un solo método.** SAPPHIRE (2020, ANN, *Pseudomonas*)
   y SAPPHIRE.CNN (2022, CNN por especie, dRNA-seq) son publicaciones
   distintas; no deben citarse como la misma herramienta.
3. **«PePPER» y «ProPr» no son lo mismo.** El dominio de PePPER sirve ahora
   ProPr, con modelos CNN por especie/GC y actualizaciones de 2025. ProPr no
   tiene artículo peer-reviewed localizado en Crossref (2026-11-26).
4. **«iPromoter» es una familia.** Existen iPromoter-2L, iPromoter-2L2.0,
   iPromoter-BnCNN, iPro70-FMWin, iPro54-PseKNC, iProm-Sigma54. Decir
   «iPromoter» sin sufijo no identifica una herramienta.
5. **Sigma70Pred ≠ 70ProPred ≠ iPro70-FMWin.** Los tres predicen σ70 en
   E. coli, con métodos y datasets distintos; sus cifras no son comparables.
6. **BPROM no tiene artículo dedicado.** Se cita por el capítulo de Solovyev
   & Salamov (2011) y por la documentación de Softberry; el benchmark de
   mSystems es la validación independiente más útil.
7. **Cappable-seq.** El DOI correcto es `10.1186/s12864-016-2539-z`; el DOI
   `10.1093/nar/gkv1320` corresponde a otro artículo (G-quadruplex), no a
   Cappable-seq.
8. **«Independiente» no siempre significa independiente.** Sigma70Pred llama
   «independiente» a RegulonDB 10.8 tras entrenar en 9.0; es independencia de
   versión, no de proceso de curación ni de solapamiento de promotores.
9. **El «mejor» resultado de nuestro propio benchmark no es generalización.**
   `SIGMA_STRATEGIES.md` lo dice explícitamente: el modelo jerárquico usa el
   mismo gold y PWMs de RegulonDB. No debe presentarse como evidencia externa.
10. **Los TSS multi-sigma existen.** En K-12 hay 124 TSS con más de un sigma
    curado; una evaluación top-1 los penaliza por diseño. Reportar también
    métricas multi-etiqueta.
