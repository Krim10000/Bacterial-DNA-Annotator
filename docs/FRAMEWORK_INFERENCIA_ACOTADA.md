# Framework de Inferencia Acotada Represor–Operador–Efector

### Pipeline computacional a partir de genoma anotado + estructuras 3D

**Versión:** 1.0  
**Contexto:** bacteria completa anotada (genes, proteínas, estructuras 3D, contexto génico)  
**Objetivo:** identificar, para cada represor anotado, (i) el operador más probable en el genoma y (ii) el efector alostérico más probable, sin validación experimental, mediante inferencia física bajo restricciones biológicas.

> **Alcance del proyecto (2026-11-24):** el objetivo es un **censo de todos los
> represores**, no la revisión de un caso. La estructura 3D es evidencia
> auxiliar; el motor principal es **especificidad de familia + conservación
> evolutiva + contexto génico**. Ver
> [`CENSUS_REPRESSORES.md`](CENSUS_REPRESSORES.md).
>
> **Estado en este repositorio:** documento de referencia conceptual. No describe
> código ya integrado. La evaluación de qué existe, qué falta y qué se puede
> aplicar hoy está en
> [`FRAMEWORK_GAP_ANALYSIS.md`](FRAMEWORK_GAP_ANALYSIS.md).
>
> **Regla de lectura:** este marco produce **hipótesis priorizadas con
> incertidumbre**, no funciones confirmadas. No sustituye un ensayo de unión
> (EMSA/ITC/reportero); sirve para decidir qué medir primero.

---

## 0. Principio rector

El problema no es "estructura → diana única" en abstracto (mal condicionado,
degenerado, subdeterminado). El problema real es:

> **estructura 3D + geometría del DBD + contexto genómico + operón regulado +
> química del bolsillo → ranking probabilístico de operadores y efectores dentro
> de un espacio biológicamente acotado.**

La estructura no "deduce" la diana. La estructura **acota la geometría**; el
contexto genómico **acota el espacio de búsqueda**; el cálculo de energías libres
relativas **ordena los candidatos**. El output es una tabla priorizada con
incertidumbre, no una certeza.

### Modelo físico de dos estados (ON/OFF)

El represor se trata como un sistema efectivo de dos cuencas de energía libre
separadas por una barrera:

- **Estado ON (unido a ADN):** DBD en registro con el surco mayor, dímero con
  simetría compatible con el operador.
- **Estado OFF (libre o unido a efector):** DBD desplazado, bolsillo regulador
  ocupado o interfaz alostérica estabilizada.

La ocupación del operador se modela como:

\[
\theta_i = \frac{[R]\,K_i}{1+\sum_j [R]\,K_j},
\qquad
K_i = e^{-\Delta G_i/RT}
\]

y la probabilidad de estado:

\[
P_{\text{ON}} = \frac{1}{1+e^{\Delta G_{\text{ON}\to\text{OFF}}/RT}}
\]

Una sola estructura da una cuenca. El segundo estado se modela por homología,
docking, dinámica molecular (DM) o red de acoplamiento alostérico. **No se
necesita el paisaje completo; se necesitan los dos estados funcionales y su
diferencia de energía libre.**

> **Convención:** \(K_i\) es una constante relativa adimensional, definida
> respecto de un estado de referencia. Sin esa referencia, \(\Delta G_i\) solo
> ordena candidatos; no es una afinidad absoluta en kcal/mol.

### Nota sobre difusión facilitada (por qué la proximidad importa)

El represor no viaja como un monorriel procesivo. Su búsqueda es difusión 3D +
deslizamiento 1D intermitente + *hopping* + transferencia intersegmental. La
consecuencia operativa es que la distribución espacial del represor está sesgada
por la posición de su gen codificante, y la represión es más fuerte sobre dianas
próximas al gen del represor. **Esto convierte la proximidad genómica en un
criterio de priorización físicamente fundamentado, no en una heurística.**

> En este documento "difusión facilitada" se usa en el sentido de **búsqueda
> facilitada de la diana**, no en el sentido de transporte de membrana. La
> proximidad es un *prior*; su peso debe calibrarse contra un nulo genómico y
> contra la orientación del operón (divergón, mismo sentido, etc.).

---

## 1. Inputs requeridos

| Input | Formato | Fuente |
|---|---|---|
| Genoma completo | FASTA | Anotación |
| Anotación de genes | GFF3/BED | Anotación |
| Regiones intergénicas / promotores | BED | Derivado de GFF3 |
| Operones predichos | Tabla | Predicción (p. ej. DOOR, Rockhopper) |
| Estructuras 3D de proteínas | PDB/mmCIF | AlphaFold DB / modelado |
| Familias de reguladores | Tabla | Pfam/InterPro/NCBI CDD |
| Red metabólica / funciones génicas | Tabla | KEGG/UniProt/GO |

**Filtro inicial:** seleccionar solo proteínas anotadas como reguladores
transcripcionales (DBD detectado por HMM: HTH, zinc finger, etc.). El resto no
entra al pipeline.

**Inputs implícitos que el marco necesita y suele omitir:**

- Estado oligomérico funcional (monómero, dímero, tetrámero) y simetría del
  operador.
- Correspondencia `locus_tag` → UniProt → PDB/AFDB → estructura.
- Complejo proteína–ADN, no solo la proteína aislada.
- Biblioteca de ligandos contextualizada al genoma (metabolitos, iones,
  xenobióticos) con identificadores químicos estables (ChEBI/InChIKey).
- Verdad de referencia para calibrar: sitios de unión, efectores, dirección y,
  cuando exista, afinidad medida.

---

## 2. Arquitectura del pipeline

```text
[Genoma anotado + estructuras]
              │
              ▼
┌───────────────────────────────┐
│ MÓDULO 1: DBD → Motivo        │
│ (especificidad de secuencia)  │
└───────────────────────────────┘
              │ PWM / matriz de especificidad
              ▼
┌───────────────────────────────┐
│ MÓDULO 2: Motivo → Operador   │
│ (escaneo genómico acotado)    │
└───────────────────────────────┘
              │ sitios candidatos rankeados
              ▼
┌───────────────────────────────┐
│ MÓDULO 3: Bolsillo → Efector  │
│ (docking + acoplamiento)      │
└───────────────────────────────┘
              │ ligandos candidatos rankeados
              ▼
┌───────────────────────────────┐
│ MÓDULO 4: Integración         │
│ (score compuesto + output)    │
└───────────────────────────────┘
```

Los Módulos 1 y 3 corren en paralelo. El Módulo 2 depende del 1. El Módulo 4
integra 1 + 2 + 3.

---

## 3. MÓDULO 1 — De la estructura 3D al motivo de unión

**Objetivo:** producir una matriz de especificidad (PWM o equivalente) para cada
represor.

### 3.1 Anotación del DBD

- Buscar dominios de unión a ADN con HMM (Pfam/InterPro).
- Identificar motivos estructurales: HTH, *winged*-HTH, zinc finger, etc.
- Localizar la hélice de reconocimiento y los residuos que contactan bases.

### 3.2 Generación de la geometría del complejo

| Situación | Herramienta |
|---|---|
| Homólogo con complejo proteína–ADN en PDB | Modelado por homología (MODELLER, SWISS-MODEL) |
| Sin homólogo de complejo | **RoseTTAFold2NA** o **AlphaFold3** para predecir el complejo |
| Estructura experimental del complejo | Usar directamente |

**Advertencia:** AF3 puede recapitular interacción genérica polinucleótido–proteína
y fallar en especificidad fina. Úsalo para geometría, no como árbitro de secuencia.

### 3.3 Predicción de especificidad (motor principal)

| Método | Tipo | Uso en el pipeline |
|---|---|---|
| **DeepPBS** | DL geométrico sobre interfaz | **Motor principal**. Funciona sobre estructuras predichas. Salida interpretable por residuo. |
| **NA-MPNN** | Inversión de plegamiento de ácidos nucleicos | Alternativa/complemento. Supera métodos de docking fijo en especificidad. |
| **IDEA** | Optimización estructura–secuencia | Complemento. ~80 % de acierto en distinguir unidores fuertes/débiles en benchmark HT-SELEX. |
| **3M-BSIP** | Multimodal, sin estructura de ADN | Filtro de residuos de interfaz cuando la geometría del complejo es incierta. |
| ~~TF2DNA~~ | Homología + potencial empírico | **Descartado.** Tasa de éxito 45–81 %, inaceptable para pipeline determinista. |
| ~~footprintDB~~ | Base de datos de motivos | Solo como *lookup* de familias ya caracterizadas, no como predictor primario. |

### 3.4 Refinamiento energético (opcional, para candidatos top)

- Modelar el complejo con las secuencias top del PWM.
- Calcular \(\Delta\Delta G\) relativo con **Rosetta** (protocolo DNA-interface)
  entre candidatos.
- Rosetta correlaciona con fuerza de represión observada en sistemas como LacI.
  No da precisión absoluta; ordena.

### 3.5 Salida del Módulo 1

- PWM por represor (formato TRANSFAC/JASPAR).
- Tabla de residuos de interfaz y su contribución a la especificidad.
- Score de confianza del motivo (basado en identidad con plantilla y
  consistencia entre métodos).

---

## 4. MÓDULO 2 — Del motivo al operador en el genoma

**Objetivo:** escanear solo las regiones biológicamente plausibles y rankear
sitios candidatos.

### 4.1 Acotación del espacio de búsqueda

No se escanea el genoma entero. Se restringe a:

1. **Regiones intergénicas** (50–300 pb corriente arriba del TSS; región
   promotora del operón).
2. **Adyacencia al gen del represor** (transcripción divergente o proximidad
   física). **Peso alto** por difusión facilitada.
3. **Promotores de operones funcionalmente coherentes** con la anotación del
   regulador (p. ej. represor de metabolismo de azúcares → operones de
   transporte/catabolismo de azúcares).

### 4.2 Escaneo de motivos

| Herramienta | Función |
|---|---|
| **FIMO** (MEME Suite) | Escaneo estándar de PWM, restringido a coordenadas BED de regiones intergénicas |
| **PREDetector 4.0** | Escaneo sobre genoma bacteriano con PWMs personalizadas, umbrales dinámicos, integración de genes co-transcritos y enlaces KEGG/UniProt |
| **COMMBAT** | Score compuesto: *interaction score* (motivo) + *target score* (contexto genómico y coherencia funcional). Reduce falsos positivos que el PWM solo no filtra. |

### 4.3 Criterios de priorización de sitios

- Puntuación del PWM (FIMO/PREDetector).
- Contexto genómico (COMMBAT *target score*).
- Proximidad al gen del represor.
- Conservación de sintenia (si hay genomas cercanos).
- Coherencia funcional del gen/operón corriente abajo.

### 4.4 Salida del Módulo 2

- Tabla de sitios candidatos por represor: coordenadas, secuencia, score PWM,
  score de contexto, score compuesto.
- Top 3–5 operadores candidatos por represor.

---

## 5. MÓDULO 3 — Del bolsillo al efector

**Objetivo:** identificar el ligando que estabiliza el estado OFF (o ON) y
rankearlo.

### 5.1 Detección del bolsillo regulador

- **fpocket** o **SiteMap** sobre el dominio regulador (típicamente C-terminal en
  TetR; interfaz de dímero en otras familias).
- Comparar con homólogos co-cristalizados con ligando en PDB para localizar el
  sitio.

### 5.2 Generación de candidatos contextuales

No se busca en todo el espacio químico. Candidatos:

- Metabolitos del operón regulado (sustrato, intermediarios).
- Compuestos de la misma vía metabólica (KEGG).
- Antibióticos o iones si la familia lo sugiere (TetR → tetraciclina;
  MerR → metales; MarR → salicilato).

### 5.3 Docking

- **AutoDock Vina** o **Glide** para docking ciego/contextual.
- Rankear ligandos por energía de docking y complementariedad geométrica.
- Identificar residuos de la cavidad que contactan el ligando (validación interna
  con mutagénesis *in silico*).

### 5.4 Acoplamiento alostérico

El docking dice si el ligando encaja, **no** si estabiliza el estado OFF. Para
eso:

- **Dinámica molecular comparativa** apo vs. holo (100–500 ns).
- **Análisis de modos normales** para detectar acoplamiento bolsillo ↔ DBD.
- **Red de correlación** de residuos: ¿la unión en el bolsillo altera la
  geometría de la interfaz DBD–ADN?

AlphaFold3 puede predecir el complejo proteína–ADN, pero la transición inducida
por efector requiere simulación explícita.

### 5.5 Salida del Módulo 3

- Top 3–5 efectores candidatos por represor: energía de docking, coherencia
  funcional, evidencia de acoplamiento alostérico.

---

## 6. MÓDULO 4 — Integración y output

### 6.1 Score compuesto por par (regulador, operador)

\[
S_{\text{sitio}} = w_1 \cdot \text{PWM} + w_2 \cdot \text{contexto} +
w_3 \cdot \text{proximidad} + w_4 \cdot \Delta\Delta G_{\text{Rosetta}}
\]

Pesos ajustables según disponibilidad de datos. Reportar sensibilidad a los
pesos.

### 6.2 Score compuesto por par (regulador, efector)

\[
S_{\text{efector}} = v_1 \cdot \text{docking} +
v_2 \cdot \text{coherencia funcional} +
v_3 \cdot \text{acoplamiento alostérico}
\]

### 6.3 Formato de output

Para cada represor anotado:

| Campo | Contenido |
|---|---|
| Regulador | ID, familia, DBD, gen |
| Operador top-1 | Coordenadas, secuencia, score, genes corriente abajo |
| Operadores top-2…5 | Ídem, con score |
| Efector top-1 | Nombre, energía de docking, vía metabólica |
| Efectores top-2…5 | Ídem |
| Confianza | Alta / Media / Baja según consistencia entre métodos |
| Incertidumbre | Rango de score, sensibilidad a pesos |

**El output NO es** "este represor se une a esta secuencia y responde a este
metabolito".  
**El output ES** "para este represor, el conjunto de operadores/efectores con
mayor probabilidad bajo el modelo es X, con confianza Y".

---

## 7. Estrategia de escala (cientos de reguladores)

Ejecución jerárquica para evitar cuello de botella:

1. **Fase rápida (todo el proteoma):**
   - HMM de DBD → identificar reguladores.
   - DeepPBS/NA-MPNN sobre estructuras AF3/RF2NA → PWM por regulador.
   - FIMO + PREDetector sobre regiones intergénicas → sitios candidatos.
   - Filtro de contexto (COMMBAT) → top candidatos.

2. **Fase profunda (solo pares que pasan el filtro):**
   - Modelado de complejo + Rosetta \(\Delta\Delta G\).
   - Docking de efectores + DM de acoplamiento.
   - Solo decenas de casos, no miles.

---

## 8. Qué NO afirma este framework

- No afirma deducir la diana desde una estructura aislada sin contexto.
- No afirma precisión absoluta en \(\Delta G\); usa energías relativas para
  ranking.
- No sustituye validación experimental como criterio epistemológico; la
  reemplaza por **inferencia física acotada y cuantificada** cuando el
  experimento no está disponible.
- No asume que el operador top-1 es el funcional *in vivo*; asume que es el más
  probable bajo el modelo.
- No modela el paisaje conformacional completo; usa el modelo efectivo de dos
  estados.

> **Matiz necesario:** "reemplazar la validación experimental" es demasiado
> fuerte. La inferencia acotada **prioriza** hipótesis y **cuantifica su
> incertidumbre**; no convierte una hipótesis en un hecho. Un ranking top-1 sin
> validación sigue siendo una hipótesis top-1.

---

## 9. Glosario operativo

| Término | Definición en este framework |
|---|---|
| **Inferencia acotada** | Predicción probabilística dentro de un espacio restringido por geometría, genómica y función |
| **Modelo de dos estados** | Aproximación ON/OFF del paisaje de energía libre del represor |
| **Difusión facilitada** | Búsqueda diana por 3D + deslizamiento 1D + *hopping* + transferencia intersegmental |
| **Degeneración** | Muchas secuencias/ligandos compatibles con una misma estructura; se acota por contexto |
| **PWM** | Matriz de pesos posición-específica que representa la especificidad de unión |
| **Acoplamiento alostérico** | Comunicación física bolsillo regulador ↔ DBD que altera la afinidad por ADN |

---

## 10. Herramientas de referencia

| Categoría | Herramientas |
|---|---|
| Predicción de complejos | RoseTTAFold2NA, AlphaFold3, MODELLER, SWISS-MODEL |
| Especificidad estructura→secuencia | **DeepPBS**, **NA-MPNN**, **IDEA**, 3M-BSIP |
| Escaneo genómico | FIMO, PREDetector 4.0, COMMBAT |
| Energía de unión | Rosetta (DNA-interface), FoldX, MM/PBSA |
| Bolsillos y docking | fpocket, SiteMap, AutoDock Vina, Glide |
| Dinámica y alostería | GROMACS, AMBER, modos normales, redes de correlación |
| Anotación funcional | Pfam, InterPro, KEGG, UniProt, DOOR/Rockhopper |

---

## 11. Resumen en una frase

**Con estructura 3D + geometría del DBD + contexto genómico + operón regulado +
química del bolsillo, el espacio de soluciones se reduce a un conjunto pequeño y
manejable; dentro de ese conjunto, el cálculo de energías libres relativas
produce un ranking probabilístico de operadores y efectores con confianza
cuantificable, sin necesidad de validación experimental.**

---

## 12. Trazabilidad con este repositorio

| Pieza del marco | Dónde vive hoy (si existe) | Estado |
|---|---|---|
| Detección de DBD y familia | `scripts/allostery/18_genome_pipeline.py`; panel de 299 HMM | Prototipo |
| PWM por TF | `data/regulondb/gt/tf_pwms.tsv` (82 PWM de *E. coli*) | Parcial, solo K-12 |
| Estructura → motivo | No integrado | Ausente |
| Escaneo de sitios | `18_genome_pipeline.py` (ventana −250 pb del gen del TF) | Prototipo, no operón-diana |
| Operones | `data/regulation/<genoma>/operons.tsv`; gate L2 18/20 | Implementado |
| Regulón TF→operón | `data/regulation/ecoli_k12/regulon.json` (65 aristas, 7 TF con aristas; 8 modelados) | Parcial, solo K-12 |
| Bolsillo y docking | `scripts/allostery/03_pocket_specificity.py` | Benchmark, no pipeline |
| Acoplamiento alostérico | `scripts/allostery/06_anm_modal.py`, `07_laci_oligomer.py` | Experimental |
| Integración y GUI | `scripts/allostery/app_alloreg.py` | Prototipo |

La evaluación detallada, con números y prioridades, está en
[`FRAMEWORK_GAP_ANALYSIS.md`](FRAMEWORK_GAP_ANALYSIS.md).
