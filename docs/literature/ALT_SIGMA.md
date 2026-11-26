# Promotores de sigma alternativos: arquitecturas, evidencia y por qué una PWM única no basta

> **Qué es:** una revisión de los motivos de promotor de las sigma alternativas
> (σ54, ECF, σ28, σ32, σ38/σB) y de los elementos accesorios (-10 extendido,
> elemento UP, espaciador, discriminador), con la evidencia experimental que
> sostiene cada consenso y los límites de esa evidencia.
>
> **Qué no es:** no propone consensos nuevos ni matrices nuevas; no afirma que
> las PWM sean inútiles; no toca código ni datos. Es insumo para decidir qué
> parte del problema de asignación de sigma es modelable con PWM y qué parte no.
>
> **Alcance:** promotores bacterianos, con foco en *E. coli* K-12 y
> *B. subtilis*. La parte de familias, dominios y ortología está en
> [`SIGMA_BIOLOGY.md`](SIGMA_BIOLOGY.md); la de métodos y benchmark, en
> [`SIGMA_PREDICTION.md`](SIGMA_PREDICTION.md). Este documento se concentra en
> **los motivos y en su poder explicativo**.
>
> **Verificación:** todas las citas se comprobaron el 2026-11-26 contra PubMed
> (E-utilities) y Europe PMC. DOI, PMID y estado de cada servidor están en el
> [Apéndice A](#apéndice-a-registro-de-verificación-2026-11-26).

---

## 0. Resumen ejecutivo

1. Las sigma alternativas **no comparten una arquitectura única de promotor**.
   σ54 usa -24/-12 y no -35/-10; las ECF usan -35/-10 pero con motivos
   específicos de grupo; σ28, σ32 y σ38 usan el par -35/-10 con énfasis y
   tolerancias distintas.
2. Un **consenso es el resumen de un alineamiento**, no una regla de unión.
   Depende del organismo, del N y del criterio de inclusión de promotores.
3. Hay **secuencias genómicas muy parecidas al consenso que no funcionan** como
   promotores: el caso de σ28 es el más explícito (Ide et al. 1999).
4. Parte de la información de reconocimiento **no está en la identidad de las
   bases**: en σE el AA central del -35 importa por la geometría que induce, y
   en las ECF la apertura de -10 implica el volteo de una sola base.
5. Los elementos **se compensan entre sí**: un -35 ausente puede sustituirse por
   un -10 extendido más fuerte; un espaciador desviado puede compensarse con un
   espaciador rico en A/T; un UP fuerte cambia la fuerza sin cambiar el núcleo.
6. La actividad depende de **factores trans** (activadores AAA+ de σ54,
   anti-sigma, competencia por el core), que una PWM no ve.
7. Existen **promotores solapados y multi-sigma**; asignar una sola etiqueta es
   una decisión de modelado, no un hecho biológico.
8. Con N pequeño por sigma, una PWM puede **ajustarse al fondo del conjunto**
   (contenido GC/AT) en lugar del motivo. Las matrices de *B. subtilis* para
   σ38/σ28/σ32/σ54 del proyecto tienen N bajo y deben tratarse como hipótesis
   (ver [`../PWM_SUBTILIS.md`](../PWM_SUBTILIS.md)).

---

## 1. Ocho razones por las que una PWM única no alcanza

La clase σ70 se delimitó por comparación de secuencias (Lonetto et al. 1992);
dentro de ella, las ECF corresponden al grupo IV (Lane & Darst 2006), y σ54 es
una familia aparte. Una PWM (position weight matrix) asume, para cada posición, una contribución
independiente de la base, y una ventana fija. Esas dos suposiciones se rompen de
formas distintas en cada familia. Las razones, con su evidencia:

| # | Razón | Evidencia | Sección |
|---|---|---|---|
| 1 | Arquitecturas incompatibles entre familias | σ54 usa -24/-12; σ70-like usa -35/-10 (Bush & Dixon 2012) | §2, §3 |
| 2 | El consenso es un promedio, no una regla | Ide et al. 1999: homólogos genómicos no funcionales | §4, §9.2 |
| 3 | La información incluye geometría y mecánica de apertura | Lane & Darst 2006 (AA central por forma); Campagne et al. 2014 (base volteada) | §3, §9.3 |
| 4 | Los elementos son condicionales y se compensan | Mitchell et al. 2003 (ext -10 vs -35); Typas & Hengge 2006 (espaciador) | §7, §9.4 |
| 5 | La fuerza depende del UP y de [RNAP] | Rhodius et al. 2012: dos clases de promotor σE | §7.2, §10 |
| 6 | La actividad depende de factores trans | bEBP + ATP en σ54; anti-sigma; competencia por core | §2, §11 |
| 7 | Hay promotores multi-sigma solapados | Qiu & Helmann 2001; Huang et al. 1998; Typas & Hengge 2006 | §11 |
| 8 | Con N pequeño, la PWM aprende el fondo | N bajo en matrices de *B. subtilis*; HMM en Price et al. 2001 | §12.3 |

La conclusión práctica no es «descartar la PWM», sino **separar dos preguntas**:
¿hay un sitio de unión compatible? (la PWM sirve, con umbral y calibración) y
¿este promotor se transcribe, con qué fuerza y con qué sigma? (la PWM sola no
alcanza).

---

## 2. σ54: -24/-12, activadores y ATP

### 2.1 Arquitectura

σ54 es una familia aparte, sin similitud de secuencia relevante con la clase
σ70 (Bush & Dixon 2012). Reconoce dos elementos:

- **-24**: GG.
- **-12**: TGC/GC.

El consenso extendido que recoge Bush & Dixon (2012) es:

```text
YTGGCACGrNNNTTGCW
```

donde las mayúsculas son residuos muy conservados, las minúsculas débilmente
conservados, `N` no conservado, `Y` = C/T, `R` = A/G, `W` = A/T. El espaciado
entre ambos elementos es corto y **crítico**: las deleciones en la región entre
-12 y -24 eliminan la función, y los elementos deben quedar sobre la misma cara
de la hélice (Bush & Dixon 2012). El registro exacto del consenso varía entre
compilaciones; conviene no tratar la cadena como posiciones absolutas.

La compilación de referencia es la de Barrios et al. (1999): **186 secuencias**
-24/-12, con un consenso extendido que **aumenta la probabilidad de identificar
promotores genuinos** frente al consenso corto. Es decir: el consenso ampliado
mejora la detección, pero sigue siendo una descripción estadística de un
conjunto heterogéneo.

### 2.2 Mecanismo: por qué la secuencia no alcanza

A diferencia de Eσ70, Eσ54 forma un **complejo cerrado estable que no se
isomeriza espontáneamente**. La apertura requiere una proteína activadora
(bEBP) que hidroliza ATP (Bush & Dixon 2012):

- los bEBP son **AAA+**; el motivo conservado **GAFTGA** contacta σ54;
- el dominio C-terminal reconoce secuencias **UAS/enhancer** situadas
  típicamente **80-150 pb** aguas arriba del TSS;
- el DNA se curva (a menudo con IHF o HU) para poner el activador en contacto
  con la holoenzima;
- la región I de σ54 se une a -12 y bloquea la carga del DNA en el sitio activo;
  el activador remodea ese bloqueo.

Ejemplos de bEBP: NtrC, NifA, PspF, XylR, DmpR, DctD, NorR, ZraR, FhlA
(Bush & Dixon 2012, Tabla 1). Una revisión temprana del factor y de su
dependencia de activador es Buck et al. (2000); el detalle del mecanismo de
remodelación está en Wigneshweraraj et al. (2008). En *E. coli* K-12, el σ54 es RpoN; en
*B. subtilis*, el equivalente es **SigL (σL)**, identificado por su homología con
σ54 de gramnegativos y por el fenotipo de mutantes nulos (Débarbouillé et al.
1991).

### 2.3 Qué implica para un modelo de secuencia

- **Detección:** una PWM sobre la ventana -24/-12 funciona razonablemente bien,
  porque el motivo es más conservado que el de σ70 (Bush & Dixon 2012). Es el
  caso más favorable de toda la familia.
- **Actividad:** no es predecible desde la secuencia del promotor. Depende de la
  presencia, posición y estado del activador, de la geometría del DNA y del ATP.
- **Falsos positivos:** un sitio -24/-12 perfecto sin UAS ni activador es un
  promotor silencioso. Un escáner que solo mire la ventana del promotor no puede
  distinguirlos.

---

## 3. ECF (grupo IV): σE, σW, σX y compañía

### 3.1 Qué son

Las ECF (extracytoplasmic function) son el subgrupo **más grande y diverso** de
sigma alternativas; responden a estrés de envoltura y señalización de superficie
(Lane & Darst 2006). La subfamilia se definió originalmente al analizar `sigE`
de *Streptomyces coelicolor* (Lonetto et al. 1994). Se clasifican como **grupo
IV de la clase σ70** y se distinguen por su mecanismo de reconocimiento, no solo
por su secuencia.

### 3.2 Arquitectura y consenso: el caso σE

Para **σE de *E. coli***, el elemento -35 de consenso es:

```text
GGAACTT
```

La estructura cristalográfica del dominio σE4 unido a ese elemento (2,3 Å)
mostró algo importante para el modelado: **el AA central es esencial aunque no
hay contactos directos proteína-base con esas bases**. El AA induce una geometría
de DNA tipo tracto AA/TT (eje recto, surco menor estrecho), y esa geometría es la
que se lee. Los autores señalan que estos principios se aplican a un rango amplio
de ECF (Lane & Darst 2006).

Consecuencia directa: una PWM puede aprender que AA es frecuente en esa posición,
pero **no puede representar el mecanismo** (la forma), y por lo tanto no
generaliza bien a variantes que preservan la geometría sin preservar la
identidad exacta.

### 3.3 La apertura en -10 también es no canónica

En las ECF, la apertura en -10 ocurre por el **volteo de una única base** de la
cadena no molde, seleccionada por un loop proteico dinámico que puede sustituirse
para cambiar la especificidad (Campagne et al. 2014). Esto refuerza que la
posición -10 no es solo una columna de PWM: es un evento mecanístico con
requisitos de secuencia contextuales.

### 3.4 Especificidad fina: σX vs σW en *B. subtilis*

*B. subtilis* tiene al menos **17 sigma, siete de ellas ECF** (Huang et al.
1998). σX se identificó por búsqueda dirigida por consenso (Huang & Helmann
1998). Dos de las mejor caracterizadas:

| Sigma | -10 que reconoce | Referencia |
|---|---|---|
| σX | CGAC (y potencialmente CGTC) | Qiu & Helmann 2001 |
| σW | CGTA (y potencialmente CGTC) | Qiu & Helmann 2001 |

Lo importante es el experimento: **cambiar solo el elemento -10 basta para
cambiar un promotor de regulón** (σX→σW o σW→σX) in vitro e in vivo, aunque los
efectos de contexto modulan la fuerza (Qiu & Helmann 2001). Es decir, la
diferencia entre dos regulones puede residir en **una sola posición** del
promotor, y esa posición es precisamente donde una PWM entrenada con pocos
sitios tiene más ruido.

Otros hechos verificados de este sistema:

- σW está **autoregulada** (un promotor dependiente de σW precede a `sigW`) y su
  regulón **solapa parcialmente** con el de σX (Huang et al. 1998).
- El regulón de σW se definió comparando **búsqueda por consenso, ROMA y
  perfilado transcripcional**. El conjunto mínimo consensuado es de **≥30 sitios
  (~60 genes)**, y **ningún método solo recupera más del ~80 %** del regulón;
  algunos promotores σW no se detectan por perfilado por enmascaramiento de otros
  elementos o por reconocimiento solapado con otras ECF (Cao et al. 2002). Es un
  buen recordatorio de que «regulón» depende del método.
- En *B. subtilis*, σW se asocia a agentes que dañan la membrana, σX a
  resistencia a péptidos antimicrobianos catiónicos y σV a lisozima; los
  regulones se solapan (Helmann 2016).

### 3.5 Un modelo de promotor ECF completo

Rhodius et al. (2012) presentaron el primer modelo de promotor bacteriano **de
longitud completa**, aplicado a σE con 60 promotores: combinaron un modelo de
**UP-element** (descrito por longitud y frecuencia de tractos A/T) con un modelo
PWM del núcleo.
El resultado:

- predijo fuerza in vivo e in vitro, y distinguió promotores activos de
  débiles/inactivos;
- la contribución del UP **varía con la disponibilidad de RNAP**;
- el regulón σE contiene **dos clases de promotor**, una activa solo a
  concentraciones altas de RNAP.

Esto es el argumento más limpio contra «una PWM única»: aun restringido a un
solo sigma, el modelo necesario tiene **dos módulos** (UP + núcleo) y una
dependencia del contexto celular.

---

## 4. σ28 / flagelar

### 4.1 Arquitectura

σ28 (FliA en *E. coli*; σD en *B. subtilis*) transcribe los genes flagelares de
**clase 3**. El consenso clásico era:

```text
TAAA N15 GCCGATAA
```

Ide et al. (1999) revisaron la estructura y propusieron:

```text
TAAAGTTT N11 GCCGATAA
```

Obsérvese que el elemento -10 (`GCCGATAA`) tiene **9 nt**, más largo que el
hexámero canónico de σ70. El cambio de N15 a N11 no es cosmético: el -35 pasa a
definirse por 8 nt (`TAAAGTTT`) y el espaciador se acorta, de modo que cambia qué
secuencias cuentan como «parecidas al consenso», aunque el total entre elementos
se conserve.

### 4.2 La evidencia que importa

El mismo trabajo (Ide et al. 1999) encontró que el genoma de *E. coli* contiene
**muchas secuencias homólogas al consenso que no confieren actividad dependiente
de σ28**. Es una de las demostraciones más directas de que *hit de consenso* ≠
*promotor funcional*, y de que el umbral de una PWM no puede fijarse solo por
parecido al consenso.

### 4.3 Regulación

σ28 está inhibida por el anti-sigma **FlgM**; FlgM se secreta cuando se completa
el cuerpo basal del flagelo (hook-basal body), lo que acopla la transcripción de
clase 3 al estado de ensamblado (Chevance & Hughes 2008). De nuevo: la secuencia
del promotor no informa si σ28 está disponible.

### 4.4 Solapamiento con σ32

Briat et al. (1985) mostraron que **σ28 de *B. subtilis* y σ32 de *E. coli*
tienen especificidad de promotor solapada**. Es un ejemplo temprano de que
sigma de especies distintas, con papeles celulares distintos, pueden reconocer
secuencias parecidas, lo que complica cualquier asignación por matriz única.

---

## 5. σ32 / heat shock

### 5.1 Arquitectura y consenso

El consenso propuesto por Cowing et al. (1985), a partir de los promotores de
`dnaK`, `groE` y el gen de C62.5, es:

```text
-35: T-N-t-C-N-C-c-C-T-T-G-A-A
-10: C-C-C-C-A-t-T-a
```

(Se reproduce con la notación del artículo: mayúsculas, posiciones conservadas;
minúsculas, tendencia; `N`, cualquiera.)

Dos observaciones sobre la evidencia:

1. **La muestra era pequeña**: tres loci (los operones `dnaK` y `groE` y el gen
   de C62.5). Es el tipo de consenso que hoy llamaríamos «de baja cobertura».
2. Los propios autores propusieron un modelo en el que **σ32 reconocería
   principalmente la región -10**, con un -35 menos determinante. El modelo fue
   superado, pero ilustra cuán pronto la evidencia admitía más de una lectura.

### 5.2 Lo que agregó la genómica

Nonaka et al. (2006) combinaron análisis de expresión y validación de promotores
para identificar genes regulados directamente por σ32, y **triplicaron el número
de genes validados**. Es decir, el regulón validado es mucho mayor que el
conjunto que dio origen al consenso: el consenso clásico no puede tomarse como
una descripción completa del regulón.

### 5.3 Consecuencia

El caso σ32 muestra el patrón: **consenso derivado de pocos promotores fuertes →
regulón validado mucho más amplio y heterogéneo**. Una PWM calibrada con los
promotores canónicos tenderá a puntuar bajo los miembros débiles del regulón, que
pueden ser funcionales igual.

---

## 6. σ38 / stationary (y σB)

### 6.1 El caso incómodo: RpoS reconoce lo mismo que σ70

En *E. coli*, σS (RpoS) reconoce **in vitro casi los mismos consensos -35 y -10
que σ70** (Typas & Hengge 2006). No existe, por lo tanto, un «consenso σ38»
nítido que la distinga: la especificidad in vivo se logra por la **capacidad
diferencial de tolerar desviaciones** del consenso.

Los detalles verificados de Typas & Hengge (2006):

- muchos promotores σS naturales **sí poseen un -35**, a menudo mal espaciado
  (contra la idea de que el -35 no está conservado en σS);
- una desviación de **±2 pb** respecto del espaciador óptimo de **17 pb**, o la
  ausencia de -35, **reduce la actividad total pero favorece a EσS** en su
  competencia con Eσ70;
- un espaciador **rico en A/T** puede compensar la pérdida de actividad sin
  perder selectividad por σS;
- las regiones 2.5 y 4 de σS estarían implicadas en leer -35 subóptimos.

### 6.2 Por qué esto rompe una PWM

Si σS y σ70 comparten el núcleo, una PWM entrenada con promotores σS **no puede
separar** σS de σ70: la información discriminante no está en las posiciones
canónicas, sino en **combinaciones débiles y distribuidas** (espaciador,
tolerancia a desviaciones, composición del espaciador, UP, discriminador). Es el
caso paradigmático de «una PWM única no basta», y el más relevante para K-12.

Una revisión del mecanismo y de las bases estructurales de la selectividad está
en Typas et al. (2007).

El control de RpoS es además multinivel (transcripción, traducción, estabilidad y
actividad), revisado por Hengge-Aronis (2002) y Battesti et al. (2011); el
regulón se ha analizado a escala genómica (Weber et al. 2005). Nada de
eso es visible en la secuencia del promotor.

### 6.3 σB en *B. subtilis*: análogo funcional, no equivalente

En *B. subtilis* el papel de respuesta general al estrés lo cumple **σB**:

- ~**150 genes** de estrés general dependientes de SigB (Hecker et al. 2007);
- la actividad de σB está controlada por una cascada con **al menos tres vías**
  (estrés ambiental, depleción de energía, baja temperatura) (Hecker et al.
  2007);
- Price et al. (2001) usaron macroarrays de >3.900 ORFs, identificaron 127
  candidatos (30 ya conocidos) y construyeron **HMM** para posibles secuencias de
  reconocimiento σB;
- Petersohn et al. (2001) confirmaron 62 genes previos y detectaron 63 nuevos
  (125 en total), y al menos **24 parecían sujetos a un segundo mecanismo de
  inducción independiente de σB**.

Ese último dato es relevante para el solapamiento: parte del «regulón σB» se
induce también sin σB. Un modelo de promotor no puede decidir eso.

**Advertencia de transferencia:** σB y RpoS cumplen papeles análogos en grupos
distintos de bacterias, pero su regulación y su repertorio de promotores difieren
(Hecker et al. 2007; Hengge-Aronis 2002). Transferir una PWM de σ38 de *E. coli*
a *B. subtilis* no está justificado por la literatura; el proyecto ya construye
matrices propias para *B. subtilis*, y debe seguir haciéndolo.

---

## 7. Elementos accesorios: -10 extendido, UP, espaciador, discriminador

Estos elementos son la razón por la que el «núcleo -35/-10» no agota el
promotor. Todos quedan **fuera de una PWM de dos hexámeros**.

### 7.1 -10 extendido (extended -10)

- Motivo **5'-TG-3'** en la región -10 extendida (Mitchell et al. 2003).
- En una compilación de **554 promotores** de *E. coli*, los promotores con
  -10 extendido tienden a tener **espaciadores más largos**, **menos
  coincidencias con el -35** y **tramos cortos de T** en el espaciador
  (Mitchell et al. 2003).
- El TG contribuyó a la actividad en **siete promotores** analizados, y su
  importancia **dependía de la secuencia de los otros elementos** (Mitchell et
  al. 2003).

Consecuencia: el -10 extendido **sustituye funcionalmente** al -35 en algunos
promotores. Un modelo que exija -35 los pierde; uno que ignore la condicionalidad
los sobre-puntúa.

### 7.2 Elemento UP

- Secuencias **ricas en A/T** aguas arriba del -35 que se unen a la **αCTD** de
  la RNAP.
- Estrem et al. (1998) combinaron selección in vitro y cribado in vivo:
  obtuvieron **31 secuencias upstream** que aumentaban la actividad del promotor
  `rrnB P1` más que cualquier UP conocido; un híbrido representativo llegó a
  ~**100 veces** la actividad del promotor core.
- La subestructura del UP y su interacción con αCTD se analizaron en Estrem et
  al. (1999).
- Para σE, el UP se describe mejor por **longitud y frecuencia de tractos A/T**
  que por una matriz posición a posición (Rhodius et al. 2012).

Consecuencia: el UP es un **módulo estadístico distinto**, no una columna más de
la PWM. Y su contribución depende de la concentración de RNAP (Rhodius et al.
2012).

### 7.3 Espaciador

- El óptimo es **17 pb** entre -35 y -10 (Typas & Hengge 2006); la compilación
  clásica de Hawley & McClure (1983) es la referencia habitual para la
  arquitectura -35/-10 de σ70.
- Typas & Hengge (2006) usan 17 pb como óptimo y muestran que **±2 pb** cambia
  actividad y selectividad.
- Los promotores con -10 extendido tienden a espaciadores más largos (Mitchell
  et al. 2003).

Consecuencia: el espaciador es una **feature de longitud**, no una posición de
matriz. Un PWM de ventana fija lo colapsa o lo ignora.

### 7.4 Discriminador

- Región **rica en GC entre -10 y +1**; se asocia al control stringent y a la
  regulación por tasa de crecimiento (Travers 1980; Lamond & Travers 1985).
- Los efectos de ppGpp sobre la iniciación dependen de la secuencia del promotor
  (Barker et al. 2001).

Consecuencia: el discriminador está **downstream del -10**, típicamente fuera de
la ventana que se escanea. Su aporte a la fuerza y a la selección por condiciones
no es capturable con una PWM de -35/-10.

---

## 8. Tabla maestra: sigma, arquitectura, consenso y evidencia

> **Ningún consenso de esta tabla es universal.** Cada uno proviene de un
> alineamiento con un organismo, un N y un criterio de inclusión concretos. La
> columna «base» indica de dónde sale; la columna «límite» indica qué no puede
> concluirse de él.

| Sigma (ejemplos) | Organismo | Arquitectura | Consenso citado (base) | Evidencia principal | Límite del motivo |
|---|---|---|---|---|---|
| **σ70/RpoD** | K-12 | -35 + -10; espaciador 17 ± 1; -10 ext; UP; discriminador | TTGACA / TATAAT (Bush & Dixon 2012) | compilaciones; información del sitio 6,5 ± 0,1 bits (Shultzaberger et al. 2007) | muchos promotores sin -35; UP y discriminador fuera de la ventana |
| **σ54/RpoN** (K-12), **SigL/σL** (*B. subtilis*) | ambos | -24 (GG) y -12 (TGC); espaciador corto; UAS 80-150 pb | YTGGCACGrNNNTTGCW (Bush & Dixon 2012; 186 promotores en Barrios et al. 1999) | activador AAA+ obligatorio; complejo cerrado estable (Bush & Dixon 2012) | la actividad; el activador, el ATP y la geometría |
| **σE/RpoE** | K-12 | -35 GGAACTT + -10; UP | GGAACTT (Lane & Darst 2006) | estructura 2,3 Å; AA central por geometría; modelo UP (Rhodius et al. 2012) | lectura por forma; UP; dos clases de promotor |
| **σW** | *B. subtilis* | -10 CGTA/CGTC; -35 tipo ECF | -10 CGTA (Qiu & Helmann 2001) | autoregulación; ROMA; solape con σX (Huang et al. 1998; Cao et al. 2002) | competencia con σX; -35 específico de grupo |
| **σX** | *B. subtilis* | -10 CGAC/CGTC | -10 CGAC (Qiu & Helmann 2001) | cambio de un solo nucleótido cambia de regulón | idem |
| **σ28/FliA** (K-12), **σD** (*B. subtilis*) | ambos | -35 TAAAGTTT + -10 GCCGATAA (9 nt), N11 | TAAAGTTT-N11-GCCGATAA (Ide et al. 1999) | homólogos genómicos no funcionales; FlgM (Chevance & Hughes 2008) | hits no funcionales; contexto; solape con σ32 (Briat et al. 1985) |
| **σ32/RpoH** | K-12 | -35 + -10 | -35 TNtCNCcCTTGAA; -10 CCCCATtTa (Cowing et al. 1985) | 3 loci originales; genómica posterior triplica el regulón (Nonaka et al. 2006) | muestra pequeña; miembros débiles del regulón |
| **σ38/RpoS** | K-12 | mismos elementos que σ70; espaciador y tolerancias | sin consenso propio: casi idéntico a σ70 in vitro (Typas & Hengge 2006) | compite con σ70; ±2 pb; espaciador A/T compensa | indistinguible del núcleo σ70; UP y discriminador |
| **σB** | *B. subtilis* | -35/-10 tipo σ70, degenerados | sin consenso único en las fuentes revisadas (Price et al. 2001; Petersohn et al. 2001) | ~150 genes; HMM; 125 genes | no transferible desde RpoS; regulación por cascada |

---

## 9. Evidencia experimental y límites de los motivos

### 9.1 Mutagénesis de un solo nucleótido

Los experimentos de intercambio de -10 entre σX y σW (Qiu & Helmann 2001)
muestran que **una posición puede decidir el regulón**. Esto tiene dos lecturas:
la buena (la biología es específica y modelable) y la mala (la señal es frágil y
una PWM con pocos datos puede no capturarla).

### 9.2 Homólogos genómicos que no funcionan

Ide et al. (1999) encontraron muchas secuencias del genoma de *E. coli*
homólogas al consenso σ28 que **no** tienen actividad dependiente de σ28. Es
evidencia directa de que la puntuación de matriz no es suficiente para declarar
promotor. Cualquier umbral de PWM debería calibrarse contra **actividad**, no
contra parecido al consenso.

### 9.3 Estructura: forma y base volteada

- En σE, el AA central del -35 es esencial **sin contactos directos con las
  bases**: actúa por la geometría del DNA (Lane & Darst 2006).
- En ECF, la apertura de -10 implica el **volteo de una sola base**, seleccionada
  por un loop que puede cambiarse para alterar la especificidad (Campagne et al.
  2014).

Ambos resultados apuntan a lo mismo: la PWM, que solo codifica identidad de base
por posición, **no representa el mecanismo**. Funciona como aproximación
estadística, no como modelo físico. La lectura de la forma del DNA y las
desviaciones de la aditividad están documentadas de forma general (Siggers &
Gordân 2014; Slattery et al. 2014).

### 9.4 Compensación entre elementos

- Un -35 ausente puede sustituirse por un -10 extendido fuerte; la importancia
  del TG depende de los otros elementos (Mitchell et al. 2003).
- Un espaciador desviado reduce actividad pero puede compensarse con un
  espaciador rico en A/T (Typas & Hengge 2006).
- El UP puede cambiar la fuerza sin cambiar el núcleo, y su efecto depende de
  [RNAP] (Rhodius et al. 2012).

Esto viola la **aditividad** que asume una PWM de posiciones independientes: el
efecto de una posición depende de las demás.

### 9.5 Información y tamaño de muestra

Shultzaberger et al. (2007) estimaron la conservación de los sitios σ70 en
**6,5 ± 0,1 bits**, y mostraron que algunos promotores **sin -35** tienen un -10
extendido de **6,7 ± 0,2 bits**: casi la misma información que un promotor
bipartito. Es decir, hay **arquitecturas alternativas con información
equivalente**; una PWM bipartita no las ve.

Además, el contenido de información de un motivo pone un techo a lo que se puede
aprender con pocos ejemplos: con N pequeño, la matriz estimada tiene alta
varianza y puede reflejar el fondo del conjunto (sesgo AT/GC) más que el motivo.
El caso de *B. subtilis* en el proyecto es elocuente: las matrices σ38, σ28, σ32 y
σ54 se construyeron con 182, 71, 13 y 11 TSS respectivamente, y el propio
documento las declara hipótesis de baja cobertura
([`../PWM_SUBTILIS.md`](../PWM_SUBTILIS.md)).

### 9.6 Competencia y disponibilidad

- Las sigma alternativas **compiten por el core** de la RNAP; su disponibilidad
  depende de anti-sigma y de otras moléculas (Österberg et al. 2011).
- La partición del «espacio de transcripción» entre sigma es una propiedad del
  sistema, no del promotor aislado (Gruber & Gross 2003).
- En σ54, la actividad requiere activador y ATP (Bush & Dixon 2012).

Una perspectiva histórica y estructural de la familia, con la partición de
funciones entre sigma, está en Feklístov et al. (2014). Ninguno de estos
factores está en la secuencia del promotor.

---

## 10. Qué se puede y qué no se puede modelar con PWM

### 10.1 Supuestos de una PWM

1. **Independencia posicional**: cada base contribuye por separado.
2. **Aditividad**: el efecto de una posición no depende de las otras.
3. **Ventana fija**: el motivo tiene una longitud y un registro definidos.
4. **Una secuencia, una puntuación**: no modela competencia ni multi-etiqueta.
5. **Fondo conocido**: la puntuación se interpreta contra un modelo de fondo.

Los puntos 1-3 son razonables para el **núcleo** de σ70 en promotores canónicos;
fallan en los casos de §9.3-§9.4. Los puntos 4-5 son decisiones de calibración.

### 10.2 Tabla de modelabilidad

| Señal | ¿PWM? | Cómo tratarla |
|---|---|---|
| Hexámero -10 | Sí, parcial | PWM; recordar que la apertura implica volteo de base (Campagne et al. 2014) |
| Hexámero -35 | Sí, con umbral | PWM con umbral calibrado; débil en σ32, σ38 y ECF |
| Espaciador -35/-10 | No | feature explícita de longitud (17 óptimo; ±2 tolerado) |
| -10 extendido (TG) | Sí, como motivo aparte | PWM separada y **condicional** a la ausencia de -35 |
| Elemento UP | No | estadísticos de tractos A/T (longitud, frecuencia) |
| Discriminador | No | composición GC entre -10 y +1 |
| -24/-12 (σ54) | Sí para detección | PWM + feature de espaciador; la actividad requiere activador |
| Geometría/forma del DNA | No | modelo de shape o reglas (tractos AA/TT) |
| Activador / anti-sigma | No | capa de contexto/trans |
| Competencia por RNAP | No | priors por condición o fase de crecimiento |
| Promotores solapados | No | salida multi-etiqueta calibrada |

### 10.3 Detectabilidad vs. fuerza

Rhodius & Mutalik (2010) mostraron, para σE, que imponer **puntuaciones mínimas
por motivo** permitía discriminar promotores activos de débiles/inactivos, y
presentaron esto como una vía para mejorar la precisión de las predicciones
basadas en PWM. La lección es que el umbral **por motivo** importa tanto como la
matriz, y que un modelo de fuerza requiere más que el núcleo.

Rhodius et al. (2012) añaden la pieza que falta: **dos clases de promotor** en el
mismo regulón, una de ellas activa solo a alta concentración de RNAP. La fuerza
no es una función solo de la secuencia.

---

## 11. Promotores con múltiples sigma y solapamiento

| Ejemplo | Qué muestra | Referencia |
|---|---|---|
| σX / σW (*B. subtilis*) | solape de regulones; -10 decide especificidad; un cambio de base cambia el regulón | Qiu & Helmann 2001; Huang et al. 1998 |
| σS / σ70 (*E. coli*) | mismo núcleo in vitro; la selectividad viene de tolerar desviaciones (espaciador, -35, A/T) | Typas & Hengge 2006 |
| σ28 (*B. subtilis*) / σ32 (*E. coli*) | especificidad de promotor solapada entre sigma de especies distintas | Briat et al. 1985 |
| σ55 / σ37 (*B. subtilis*) | promotores solapados transcritos por dos holoenzimas durante crecimiento y fase estacionaria | Wang & Doi 1984 |
| Dos sigma de *B. subtilis* (gen regulado durante el desarrollo) | promotores solapados y discriminación entre sigma | Johnson et al. 1983 |
| σB (*B. subtilis*) | ≥24 de 125 genes inducidos por estrés parecen tener un segundo mecanismo independiente de σB | Petersohn et al. 2001 |
| Mapa de crosstalk ECF | subconjuntos de sigma/promotores/anti-sigma ortogonales; base para switches sintéticos | Rhodius et al. 2013 |

**Consecuencia para la anotación:** si el mismo fragmento puede ser promotor de
dos sigma, la salida correcta no es una etiqueta, sino un **ranking con
alternativas y una medida de ambigüedad**. El proyecto ya reporta top-3 y
alternativas en la capa sigma ([`../SIGMA_LAYERED_EVAL.md`](../SIGMA_LAYERED_EVAL.md));
la literatura respalda esa decisión de diseño.

---

## 12. Implicaciones para *B. subtilis* y K-12

### 12.1 Los repertorios no son comparables

| Función | K-12 | *B. subtilis* | Nota |
|---|---|---|---|
| Housekeeping | RpoD (σ70) | SigA (σA) | ortólogos |
| Estrés general / estacionaria | RpoS (σ38) | SigB (σB) | análogos funcionales; regulación distinta (Hecker et al. 2007) |
| Choque térmico | RpoH (σ32) | sin ortólogo directo claro; σI implicada en *Bacillus* (Kim & Wilson 2016, en *B. anthracis*) | **no transferir** la PWM de σ32 |
| Flagelar | FliA (σ28) | SigD (σD) | ortólogos; σD aislada por Helmann et al. (1988) |
| ECF | RpoE (σE) | siete ECF, entre ellas σM, σV, σW, σX, σY y σZ | uno-a-muchos; crosstalk (Helmann 2016) |
| σ54 | RpoN | SigL (σL) | mismos elementos; **distintos activadores y regulones** (Débarbouillé et al. 1991; Calogero et al. 1994) |
| Esporulación | no aplica | σH (transición), σE, σF, σG, σK | específicas de *B. subtilis*; sin contraparte en K-12 |

### 12.2 Colisiones de nombre que afectan el modelado

- **σE**: en *E. coli* es RpoE (ECF); en *B. subtilis* suele ser SigE
  (esporulación). Son factores distintos.
- **σF**: en *E. coli*/Salmonella es σ28 (FliA); en *B. subtilis* es SigF
  (esporulación).
- **σH**: en *E. coli* es RpoH (σ32); en *B. subtilis* es SigH (transición).
- **σ38/σB**: la etiqueta «σ38» se usa para RpoS y, en algunos contextos, para
  σB. No son ortólogos.

Un panel de seis sigma con etiquetas genéricas (`sigma70`, `sigma54`, `sigma32`,
`sigma38`, `sigma28`, `sigma24`) es una **simplificación de familia funcional**,
no un mapa de ortología. El propio proyecto lo advierte en
[`../SIGMA_GROUPS.md`](../SIGMA_GROUPS.md) y [`../SIGMA_MODULE.md`](../SIGMA_MODULE.md);
este documento aporta la razón biológica: los motivos y los regulones no se
corresponden uno a uno entre especies.

### 12.3 Consecuencias prácticas

1. **Matrices por especie, sin transferencia.** Ya es la práctica del proyecto
   ([`../PWM_SUBTILIS.md`](../PWM_SUBTILIS.md)). La literatura la respalda: los
   consensos son específicos de organismo y de conjunto.
2. **No mapear por nombre.** σ32/σ38/σE/σF/σH significan cosas distintas en cada
   especie (§12.2).
3. **Tratar las matrices de N bajo como hipótesis.** Con 11-71 TSS, la matriz
   estimada es inestable; el consenso de literatura puede usarse como **prior
   declarado**, no como matriz aprendida.
4. **Separar detección de actividad.** Para σ54, la detección por PWM es el caso
   más favorable; la actividad es el caso más desfavorable. No prometer lo
   segundo con lo primero.
5. **Modelar la ambigüedad.** En *B. subtilis*, σW/σX solapan; en K-12, σS/σ70.
   La salida multi-etiqueta no es un lujo.

### 12.4 Sobre la capa de contexto del proyecto

La capa de contexto descrita en [`../SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md)
intenta incorporar evidencia trans (operón, ruta, conservación, geometría σ54).
Este documento no la evalúa; solo aporta el fundamento de por qué una capa así es
necesaria: la PWM no puede ver activadores, anti-sigma ni competencia. Cualquier
ganancia de esa capa debe medirse fuera de muestra, como ya se hace allí.

---

## 13. Recomendaciones (sin tocar código)

1. **Un modelo por sigma y por especie.** Nunca una matriz global reutilizada
   entre K-12 y *B. subtilis*.
2. **Features explícitas y auditables** además de la PWM: longitud del espaciador,
   presencia de -10 extendido, tractos A/T del UP, GC del discriminador.
3. **Condicionalidad**: permitir que el -10 extendido sustituya al -35, y que un
   espaciador desviado se compense con composición A/T. Modelar como reglas
   separadas, no como columnas de la misma matriz.
4. **Umbral por motivo y por sigma**, calibrado contra actividad o contra
   negativos reales, no contra parecido al consenso (Ide et al. 1999).
5. **Fondo dependiente del GC** y p-valor empírico por genoma; el proyecto ya lo
   hace ([`../SIGMA_MODULE.md`](../SIGMA_MODULE.md)).
6. **Reportar N, cobertura e intervalo** de cada matriz; con N pequeño, declarar
   el prior de literatura usado.
7. **Salida multi-etiqueta/top-k** con alternativas y un indicador de ambigüedad
   cuando dos sigma compiten por la misma ventana.
8. **No modelar fuerza** con la PWM sola. Si se reporta fuerza, que sea con un
   modelo de dos módulos (núcleo + UP) y con la advertencia de que depende de
   [RNAP] y de condiciones (Rhodius et al. 2012).
9. **Para σ54**, añadir una feature de geometría/espaciador y, si hubiera datos,
   señales de UAS/activador. La secuencia del promotor no decide la actividad.
10. **Documentar la etiqueta**: ¿es un ortólogo o una familia funcional? La
    diferencia cambia cómo se interpreta el resultado.

---

## 14. Referencias verificadas

> Cada entrada se comprobó el 2026-11-26. «Usado en» indica las secciones donde
> se apoya en esa referencia. Los DOI enlazan a `doi.org`; los PMID a PubMed.

### 14.1 σ54

- **Barrios, H.; Valderrama, B.; Morett, E. (1999).** Compilation and analysis of
  σ54-dependent promoter sequences. *Nucleic Acids Research* 27(22): 4305-4313.
  DOI: [10.1093/nar/27.22.4305](https://doi.org/10.1093/nar/27.22.4305).
  PMID: [10536136](https://pubmed.ncbi.nlm.nih.gov/10536136/). Usado en §2, §8.
- **Buck, M.; Gallegos, M. T.; Studholme, D. J.; Guo, Y.; Gralla, J. D. (2000).**
  The bacterial enhancer-dependent σ54 (σN) transcription factor. *Journal of
  Bacteriology* 182(15): 4129-4136. DOI:
  [10.1128/JB.182.15.4129-4136.2000](https://doi.org/10.1128/JB.182.15.4129-4136.2000).
  PMID: [10894718](https://pubmed.ncbi.nlm.nih.gov/10894718/). Usado en §2.2.
- **Bush, M.; Dixon, R. (2012).** The role of bacterial enhancer binding proteins
  as specialized activators of σ54-dependent transcription. *Microbiology and
  Molecular Biology Reviews* 76(3): 497-529. DOI:
  [10.1128/MMBR.00006-12](https://doi.org/10.1128/MMBR.00006-12).
  PMID: [22933558](https://pubmed.ncbi.nlm.nih.gov/22933558/). Usado en §1, §2,
  §8, §9.6, §10.
- **Débarbouillé, M.; Martin-Verstraete, I.; Kunst, F.; Rapoport, G. (1991).**
  The *Bacillus subtilis* `sigL` gene encodes an equivalent of σ54 from
  gram-negative bacteria. *PNAS* 88(20): 9092-9096. DOI:
  [10.1073/pnas.88.20.9092](https://doi.org/10.1073/pnas.88.20.9092).
  PMID: [1924373](https://pubmed.ncbi.nlm.nih.gov/1924373/). Usado en §2.2, §12.1.
- **Calogero, S.; Gardan, R.; Glaser, P.; Schweizer, J.; Rapoport, G.;
  Débarbouillé, M. (1994).** RocR, a novel regulatory protein controlling
  arginine utilization in *Bacillus subtilis*, belongs to the NtrC/NifA family of
  transcriptional activators. *Journal of Bacteriology* 176(5): 1234-1241. DOI:
  [10.1128/jb.176.5.1234-1241.1994](https://doi.org/10.1128/jb.176.5.1234-1241.1994).
  PMID: [8113162](https://pubmed.ncbi.nlm.nih.gov/8113162/). Usado en §12.1.
- **Wigneshweraraj, S.; Bose, D.; Burrows, P. C.; Joly, N.; Schumacher, J.;
  Rappas, M.; Pape, T.; Zhang, X.; et al. (2008).** Modus operandi of the
  bacterial RNA polymerase containing the σ54 promoter-specificity factor.
  *Molecular Microbiology* 68(3): 538-546. DOI:
  [10.1111/j.1365-2958.2008.06181.x](https://doi.org/10.1111/j.1365-2958.2008.06181.x).
  PMID: [18331472](https://pubmed.ncbi.nlm.nih.gov/18331472/). Usado en §2.2.

### 14.2 ECF

- **Campagne, S.; Marsh, M. E.; Capitani, G.; Vorholt, J. A.; Allain, F. H.
  (2014).** Structural basis for -10 promoter element melting by environmentally
  induced sigma factors. *Nature Structural & Molecular Biology* 21(3): 269-276.
  DOI: [10.1038/nsmb.2777](https://doi.org/10.1038/nsmb.2777).
  PMID: [24531660](https://pubmed.ncbi.nlm.nih.gov/24531660/). Usado en §3.3,
  §9.3, §10.2.
- **Cao, M.; Kobel, P. A.; Morshedi, M. M.; Wu, M. F.; Paddon, C.; Helmann, J. D.
  (2002).** Defining the *Bacillus subtilis* σW regulon: a comparative analysis
  of promoter consensus search, run-off transcription/macroarray analysis
  (ROMA), and transcriptional profiling approaches. *Journal of Molecular
  Biology* 316(3): 443-457. DOI:
  [10.1006/jmbi.2001.5372](https://doi.org/10.1006/jmbi.2001.5372).
  PMID: [11866510](https://pubmed.ncbi.nlm.nih.gov/11866510/). Usado en §3.4, §8.
- **Helmann, J. D. (2016).** *Bacillus subtilis* extracytoplasmic function (ECF)
  sigma factors and defense of the cell envelope. *Current Opinion in
  Microbiology* 30: 122-132. DOI:
  [10.1016/j.mib.2016.02.002](https://doi.org/10.1016/j.mib.2016.02.002).
  PMID: [26901131](https://pubmed.ncbi.nlm.nih.gov/26901131/). Usado en §3.4,
  §12.1.
- **Huang, X.; Fredrick, K. L.; Helmann, J. D. (1998).** Promoter recognition by
  *Bacillus subtilis* σW: autoregulation and partial overlap with the σX regulon.
  *Journal of Bacteriology* 180(15): 3765-3770. DOI:
  [10.1128/JB.180.15.3765-3770.1998](https://doi.org/10.1128/JB.180.15.3765-3770.1998).
  PMID: [9683469](https://pubmed.ncbi.nlm.nih.gov/9683469/). Usado en §3.4, §11.
- **Huang, X.; Helmann, J. D. (1998).** Identification of target promoters for
  the *Bacillus subtilis* σX factor using a consensus-directed search. *Journal
  of Molecular Biology* 279(1): 165-173. DOI:
  [10.1006/jmbi.1998.1765](https://doi.org/10.1006/jmbi.1998.1765).
  PMID: [9636707](https://pubmed.ncbi.nlm.nih.gov/9636707/). Usado en §3.4.
- **Lane, W. J.; Darst, S. A. (2006).** The structural basis for promoter -35
  element recognition by the group IV sigma factors. *PLoS Biology* 4(9): e269.
  DOI: [10.1371/journal.pbio.0040269](https://doi.org/10.1371/journal.pbio.0040269).
  PMID: [16903784](https://pubmed.ncbi.nlm.nih.gov/16903784/). Usado en §3, §8,
  §9.3.
- **Qiu, J.; Helmann, J. D. (2001).** The -10 region is a key promoter
  specificity determinant for the *Bacillus subtilis* extracytoplasmic-function
  sigma factors σX and σW. *Journal of Bacteriology* 183(6): 1921-1927. DOI:
  [10.1128/JB.183.6.1921-1927.2001](https://doi.org/10.1128/JB.183.6.1921-1927.2001).
  PMID: [11222589](https://pubmed.ncbi.nlm.nih.gov/11222589/). Usado en §3.4,
  §8, §9.1, §11.
- **Rhodius, V. A.; Mutalik, V. K.; Gross, C. A. (2012).** Predicting the
  strength of UP-elements and full-length *E. coli* σE promoters. *Nucleic Acids
  Research* 40(7): 2907-2924. DOI:
  [10.1093/nar/gkr1190](https://doi.org/10.1093/nar/gkr1190).
  PMID: [22156164](https://pubmed.ncbi.nlm.nih.gov/22156164/). Usado en §3.5,
  §7.2, §10, §13.
- **Rhodius, V. A.; Mutalik, V. K. (2010).** Predicting strength and function for
  promoters of the *Escherichia coli* alternative sigma factor, σE. *PNAS*
  107(7): 2854-2859. DOI:
  [10.1073/pnas.0915066107](https://doi.org/10.1073/pnas.0915066107).
  PMID: [20133665](https://pubmed.ncbi.nlm.nih.gov/20133665/). Usado en §10.3.
- **Rhodius, V. A.; Segall-Shapiro, T. H.; Sharon, B. D.; Ghodasara, A.;
  Orlova, E.; Tabakh, H.; Burkhardt, D. H.; Clancy, K.; et al. (2013).** Design
  of orthogonal genetic switches based on a crosstalk map of σs, anti-σs, and
  promoters. *Molecular Systems Biology* 9: 702. DOI:
  [10.1038/msb.2013.58](https://doi.org/10.1038/msb.2013.58).
  PMID: [24169405](https://pubmed.ncbi.nlm.nih.gov/24169405/). Usado en §11.
- **Lonetto, M. A.; Brown, K. L.; Rudd, K. E.; Buttner, M. J. (1994).** Analysis
  of the *Streptomyces coelicolor* `sigE` gene reveals the existence of a
  subfamily of eubacterial RNA polymerase sigma factors involved in the
  regulation of extracytoplasmic functions. *PNAS* 91(16): 7573-7577. DOI:
  [10.1073/pnas.91.16.7573](https://doi.org/10.1073/pnas.91.16.7573).
  PMID: [8052622](https://pubmed.ncbi.nlm.nih.gov/8052622/). Usado en §3.1.

### 14.3 σ28, σ32, σ38/σB

- **Briat, J. F.; Gilman, M. Z.; Chamberlin, M. J. (1985).** *Bacillus subtilis*
  σ28 and *Escherichia coli* σ32 (htpR) are minor sigma factors that display an
  overlapping promoter specificity. *Journal of Biological Chemistry* 260(4):
  2038-2041. PMID:
  [3918995](https://pubmed.ncbi.nlm.nih.gov/3918995/). Usado en §4.4, §8, §11.
- **Cowing, D. W.; Bardwell, J. C.; Craig, E. A.; Woolford, C.; Hendrix, R. W.;
  Gross, C. A. (1985).** Consensus sequence for *Escherichia coli* heat shock
  gene promoters. *PNAS* 82(9): 2679-2683. DOI:
  [10.1073/pnas.82.9.2679](https://doi.org/10.1073/pnas.82.9.2679).
  PMID: [3887408](https://pubmed.ncbi.nlm.nih.gov/3887408/). Usado en §5, §8.
- **Chevance, F. F. V.; Hughes, K. T. (2008).** Coordinating assembly of a
  bacterial macromolecular machine. *Nature Reviews Microbiology* 6(6): 455-465.
  DOI: [10.1038/nrmicro1887](https://doi.org/10.1038/nrmicro1887).
  PMID: [18483484](https://pubmed.ncbi.nlm.nih.gov/18483484/). Usado en §4.3, §8.
- **Helmann, J. D.; Márquez, L. M.; Chamberlin, M. J. (1988).** Isolation and
  characterization of the *Bacillus subtilis* σ28 factor. *Journal of
  Bacteriology* 170(4): 1560-1567. DOI:
  [10.1128/jb.170.4.1560-1567.1988](https://doi.org/10.1128/jb.170.4.1560-1567.1988).
  PMID: [3127378](https://pubmed.ncbi.nlm.nih.gov/3127378/). Usado en §12.1.
- **Ide, N.; Ikebe, T.; Kutsukake, K. (1999).** Reevaluation of the promoter
  structure of the class 3 flagellar operons of *Escherichia coli* and
  *Salmonella*. *Genes & Genetic Systems* 74(3): 113-116. DOI:
  [10.1266/ggs.74.113](https://doi.org/10.1266/ggs.74.113).
  PMID: [10586520](https://pubmed.ncbi.nlm.nih.gov/10586520/). Usado en §0, §1,
  §4, §8, §9.2, §13.
- **Nonaka, G.; Blankschien, M.; Herman, C.; Gross, C. A.; Rhodius, V. A.
  (2006).** Regulon and promoter analysis of the *E. coli* heat-shock factor,
  σ32, reveals a multifaceted cellular response to heat stress. *Genes &
  Development* 20(13): 1776-1789. DOI:
  [10.1101/gad.1428206](https://doi.org/10.1101/gad.1428206).
  PMID: [16818608](https://pubmed.ncbi.nlm.nih.gov/16818608/). Usado en §5.2, §8.
- **Typas, A.; Hengge, R. (2006).** Role of the spacer between the -35 and -10
  regions in σS promoter selectivity in *Escherichia coli*. *Molecular
  Microbiology* 59(3): 1037-1051. DOI:
  [10.1111/j.1365-2958.2005.04998.x](https://doi.org/10.1111/j.1365-2958.2005.04998.x).
  PMID: [16420370](https://pubmed.ncbi.nlm.nih.gov/16420370/). Usado en §1, §6.1,
  §7.3, §8, §9.4, §11.
- **Typas, A.; Becker, G.; Hengge, R. (2007).** The molecular basis of selective
  promoter activation by the σS subunit of RNA polymerase. *Molecular
  Microbiology* 63(5): 1296-1306. PMID:
  [17302812](https://pubmed.ncbi.nlm.nih.gov/17302812/). Usado en §6.1.
- **Hengge-Aronis, R. (2002).** Signal transduction and regulatory mechanisms
  involved in control of the σS (RpoS) subunit of RNA polymerase. *Microbiology
  and Molecular Biology Reviews* 66(3): 373-395. DOI:
  [10.1128/MMBR.66.3.373-395.2002](https://doi.org/10.1128/MMBR.66.3.373-395.2002).
  PMID: [12208995](https://pubmed.ncbi.nlm.nih.gov/12208995/). Usado en §6.2,
  §6.3.
- **Battesti, A.; Majdalani, N.; Gottesman, S. (2011).** The RpoS-mediated
  general stress response in *Escherichia coli*. *Annual Review of Microbiology*
  65: 189-213. DOI:
  [10.1146/annurev-micro-090110-102946](https://doi.org/10.1146/annurev-micro-090110-102946).
  PMID: [21639793](https://pubmed.ncbi.nlm.nih.gov/21639793/). Usado en §6.2.
- **Weber, H.; Polen, T.; Heuveling, J.; Wendisch, V. F.; Hengge, R. (2005).**
  Genome-wide analysis of the general stress response network in *Escherichia
  coli*: σS-dependent genes, promoters, and sigma factor selectivity. *Journal of
  Bacteriology* 187(5): 1591-1603. DOI:
  [10.1128/JB.187.5.1591-1603.2005](https://doi.org/10.1128/JB.187.5.1591-1603.2005).
  PMID: [15716429](https://pubmed.ncbi.nlm.nih.gov/15716429/). Usado en §6.2.
- **Price, C. W.; Fawcett, P.; Cérémonie, H.; Su, N.; Murphy, C. K.; Youngman, P.
  (2001).** Genome-wide analysis of the general stress response in *Bacillus
  subtilis*. *Molecular Microbiology* 41(4): 757-774. DOI:
  [10.1046/j.1365-2958.2001.02534.x](https://doi.org/10.1046/j.1365-2958.2001.02534.x).
  PMID: [11532142](https://pubmed.ncbi.nlm.nih.gov/11532142/). Usado en §1, §6.3,
  §8.
- **Petersohn, A.; Brigulla, M.; Haas, S.; Hoheisel, J. D.; Völker, U.; Hecker,
  M. (2001).** Global analysis of the general stress response of *Bacillus
  subtilis*. *Journal of Bacteriology* 183(19): 5617-5631. DOI:
  [10.1128/JB.183.19.5617-5631.2001](https://doi.org/10.1128/JB.183.19.5617-5631.2001).
  PMID: [11544224](https://pubmed.ncbi.nlm.nih.gov/11544224/). Usado en §6.3,
  §8, §11.
- **Hecker, M.; Pané-Farré, J.; Völker, U. (2007).** SigB-dependent general
  stress response in *Bacillus subtilis* and related gram-positive bacteria.
  *Annual Review of Microbiology* 61: 215-236. DOI:
  [10.1146/annurev.micro.61.080706.093445](https://doi.org/10.1146/annurev.micro.61.080706.093445).
  PMID: [18035607](https://pubmed.ncbi.nlm.nih.gov/18035607/). Usado en §6.3,
  §12.1.
- **Kim, J. G. Y.; Wilson, A. C. (2016).** Loss of σI affects heat-shock response
  and virulence gene expression in *Bacillus anthracis*. *Microbiology* 162(4):
  564-574. DOI: [10.1099/mic.0.000236](https://doi.org/10.1099/mic.0.000236).
  PMID: [26744224](https://pubmed.ncbi.nlm.nih.gov/26744224/). Usado en §12.1.

### 14.4 Elementos accesorios del promotor

- **Estrem, S. T.; Gaal, T.; Ross, W.; Gourse, R. L. (1998).** Identification of
  an UP element consensus sequence for bacterial promoters. *PNAS* 95(17):
  9761-9766. DOI: [10.1073/pnas.95.17.9761](https://doi.org/10.1073/pnas.95.17.9761).
  PMID: [9707549](https://pubmed.ncbi.nlm.nih.gov/9707549/). Usado en §7.2.
- **Estrem, S. T.; Ross, W.; Gaal, T.; Chen, Z. W.; Niu, W.; Ebright, R. H.;
  Gourse, R. L. (1999).** Bacterial promoter architecture: subsite structure of
  UP elements and interactions with the carboxy-terminal domain of the RNA
  polymerase α subunit. *Genes & Development* 13(16): 2134-2147. DOI:
  [10.1101/gad.13.16.2134](https://doi.org/10.1101/gad.13.16.2134).
  PMID: [10465790](https://pubmed.ncbi.nlm.nih.gov/10465790/). Usado en §7.2.
- **Hawley, D. K.; McClure, W. R. (1983).** Compilation and analysis of
  *Escherichia coli* promoter DNA sequences. *Nucleic Acids Research* 11(8):
  2237-2255. DOI: [10.1093/nar/11.8.2237](https://doi.org/10.1093/nar/11.8.2237).
  PMID: [6344016](https://pubmed.ncbi.nlm.nih.gov/6344016/). Usado en §7.3.
- **Mitchell, J. E.; Zheng, D.; Busby, S. J. W.; Minchin, S. D. (2003).**
  Identification and analysis of "extended -10" promoters in *Escherichia coli*.
  *Nucleic Acids Research* 31(16): 4689-4695. DOI:
  [10.1093/nar/gkg694](https://doi.org/10.1093/nar/gkg694).
  PMID: [12907708](https://pubmed.ncbi.nlm.nih.gov/12907708/). Usado en §1, §7.1,
  §7.3, §9.4.
- **Shultzaberger, R. K.; Chen, Z.; Lewis, K. A.; Schneider, T. D. (2007).**
  Anatomy of *Escherichia coli* σ70 promoters. *Nucleic Acids Research* 35(3):
  771-788. DOI: [10.1093/nar/gkl956](https://doi.org/10.1093/nar/gkl956).
  PMID: [17189297](https://pubmed.ncbi.nlm.nih.gov/17189297/). Usado en §8, §9.5.
- **Travers, A. A. (1980).** Promoter sequence for stringent control of bacterial
  ribonucleic acid synthesis. *Journal of Bacteriology* 141(2): 973-976. DOI:
  [10.1128/jb.141.2.973-976.1980](https://doi.org/10.1128/jb.141.2.973-976.1980).
  PMID: [6154042](https://pubmed.ncbi.nlm.nih.gov/6154042/). Usado en §7.4.
- **Lamond, A. I.; Travers, A. A. (1985).** Stringent control of bacterial
  transcription. *Cell* 41(1): 6-8. DOI:
  [10.1016/0092-8674(85)90050-9](https://doi.org/10.1016/0092-8674(85)90050-9).
  PMID: [2581696](https://pubmed.ncbi.nlm.nih.gov/2581696/). Usado en §7.4.
- **Barker, M. M.; Gaal, T.; Josaitis, C. A.; Gourse, R. L. (2001).** Mechanism
  of regulation of transcription initiation by ppGpp. I. Effects of ppGpp on
  transcription initiation in vivo and in vitro. *Journal of Molecular Biology*
  305(4): 673-688. DOI:
  [10.1006/jmbi.2000.4327](https://doi.org/10.1006/jmbi.2000.4327).
  PMID: [11162084](https://pubmed.ncbi.nlm.nih.gov/11162084/). Usado en §7.4.
- **Wang, P. Z.; Doi, R. H. (1984).** Overlapping promoters transcribed by
  *Bacillus subtilis* σ55 and σ37 RNA polymerase holoenzymes during growth and
  stationary phases. *Journal of Biological Chemistry* 259(13): 8619-8625.
  PMID: [6330116](https://pubmed.ncbi.nlm.nih.gov/6330116/). Usado en §11.
- **Johnson, W. C.; Moran, C. P., Jr.; Losick, R. (1983).** Two RNA polymerase
  sigma factors from *Bacillus subtilis* discriminate between overlapping
  promoters for a developmentally regulated gene. *Nature* 302(5911): 800-804.
  DOI: [10.1038/302800a0](https://doi.org/10.1038/302800a0).
  PMID: [6405278](https://pubmed.ncbi.nlm.nih.gov/6405278/). Usado en §11.

### 14.5 Competencia, familias y límites de la PWM

- **Gruber, T. M.; Gross, C. A. (2003).** Multiple sigma subunits and the
  partitioning of bacterial transcription space. *Annual Review of Microbiology*
  57: 441-466. DOI:
  [10.1146/annurev.micro.57.030502.090913](https://doi.org/10.1146/annurev.micro.57.030502.090913).
  PMID: [14527287](https://pubmed.ncbi.nlm.nih.gov/14527287/). Usado en §9.6.
- **Österberg, S.; del Peso-Santos, T.; Shingler, V. (2011).** Regulation of
  alternative sigma factor use. *Annual Review of Microbiology* 65: 37-55. DOI:
  [10.1146/annurev.micro.112408.134219](https://doi.org/10.1146/annurev.micro.112408.134219).
  PMID: [21639785](https://pubmed.ncbi.nlm.nih.gov/21639785/). Usado en §9.6.
- **Feklístov, A.; Sharon, B. D.; Darst, S. A.; Gross, C. A. (2014).** Bacterial
  sigma factors: a historical, structural, and genomic perspective. *Annual
  Review of Microbiology* 68: 357-376. DOI:
  [10.1146/annurev-micro-092412-155737](https://doi.org/10.1146/annurev-micro-092412-155737).
  PMID: [25002089](https://pubmed.ncbi.nlm.nih.gov/25002089/). Usado en §9.6.
- **Lonetto, M.; Gribskov, M.; Gross, C. A. (1992).** The sigma 70 family:
  sequence conservation and evolutionary relationships. *Journal of
  Bacteriology* 174(12): 3843-3849. DOI:
  [10.1128/jb.174.12.3843-3849.1992](https://doi.org/10.1128/jb.174.12.3843-3849.1992).
  PMID: [1597408](https://pubmed.ncbi.nlm.nih.gov/1597408/). Usado en §1.
- **Siggers, T.; Gordân, R. (2014).** Protein-DNA binding: complexities and
  multi-protein codes. *Nucleic Acids Research* 42(4): 2099-2111. DOI:
  [10.1093/nar/gkt1112](https://doi.org/10.1093/nar/gkt1112).
  PMID: [24243859](https://pubmed.ncbi.nlm.nih.gov/24243859/). Usado en §9.3.
- **Slattery, M.; Zhou, T.; Yang, L.; Dantas Machado, A. C.; Gordân, R.; Rohs, R.
  (2014).** Absence of a simple code: how transcription factors read the genome.
  *Trends in Biochemical Sciences* 39(9): 381-399. DOI:
  [10.1016/j.tibs.2014.07.002](https://doi.org/10.1016/j.tibs.2014.07.002).
  PMID: [25129887](https://pubmed.ncbi.nlm.nih.gov/25129887/). Usado en §9.3.

---

## Apéndice A. Registro de verificación (2026-11-26)

### A.1 Método

1. **Metadatos (DOI + PMID):** PubMed E-utilities (`esearch` → `esummary`) para
   cada título; el DOI se tomó del campo `articleids` del registro MEDLINE.
2. **Resumen/abstract:** `efetch` (`rettype=abstract`, `retmode=text`) para las
   referencias cuyas afirmaciones concretas se citan en el texto (Cowing 1985,
   Ide 1999, Qiu & Helmann 2001, Mitchell 2003, Rhodius 2012, Typas & Hengge
   2006, Lane & Darst 2006, Campagne 2014, Price 2001, Petersohn 2001, Hecker
   2007, Gruber & Gross 2003, Österberg 2011, Bush & Dixon 2012). Los de Barrios
   et al. (1999) y Nonaka et al. (2006) se consultaron en Europe PMC y en la
   página del editor, respectivamente.
3. **Texto completo:** Bush & Dixon (2012) se leyó completo en PMC (PMC3429621)
   para el consenso `YTGGCACGrNNNTTGCW` y la descripción del mecanismo bEBP/ATP.

### A.2 Incidencias y límites de la verificación

- **Europe PMC REST** respondió lento o con timeouts en parte de la sesión; los
  PMID se obtuvieron por E-utilities, que fue estable.
- **Crossref** devolvió HTTP 429 (rate limit) al intentar resolver DOIs en lote;
  por eso la verificación DOI-PMID se hizo vía MEDLINE.
- **PMC HTML** devolvió reCAPTCHA en algunos intentos; solo se usó cuando
  respondió 200 (Bush & Dixon 2012).
- **Typas et al. (2007)** no tiene DOI recuperado en esta sesión: se cita por
  PMID. **Briat et al. (1985)** y **Wang & Doi (1984)** son artículos de *JBC*
  sin DOI en el registro MEDLINE consultado: se citan por PMID.
- No se citó ninguna fuente cuyo título, revista, volumen o año no se haya
  recuperado de MEDLINE/PMC o de la página del editor. **No hay citas de
  memoria.**

### A.3 Qué NO se verificó

- No se accedió al texto completo de la mayoría de los artículos: las
  afirmaciones se apoyan en el abstract y, cuando se indica, en el texto completo
  de Bush & Dixon (2012). Los detalles experimentales finos (p. ej., número exacto
  de promotores en cada subconjunto) pueden diferir de lo que el abstract resume.
- No se verificó la posición exacta del registro del consenso σ54 en cada
  compilación; el documento evita afirmar una numeración única.
- No se revisaron bases de datos de promotores (RegulonDB, DBTBS) para contar
  promotores por sigma; ese trabajo corresponde a
  [`SIGMA_DATABASES.md`](SIGMA_DATABASES.md).

---

## Ver también

- [`SIGMA_BIOLOGY.md`](SIGMA_BIOLOGY.md) — familias, dominios, ortología y censo
  de sigma por especie.
- [`SIGMA_PREDICTION.md`](SIGMA_PREDICTION.md) — métodos de predicción,
  benchmark y evaluación.
- [`SIGMA_DATABASES.md`](SIGMA_DATABASES.md) — bases de datos de promotores y
  regulones.
- [`../SIGMA_MODULE.md`](../SIGMA_MODULE.md) — módulo sigma del proyecto.
- [`../PWM_SUBTILIS.md`](../PWM_SUBTILIS.md) — matrices propias de *B. subtilis*
  y su cobertura.
- [`../SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md) — capa de contexto y evidencia
  trans.
- [`../SIGMA_GROUPS.md`](../SIGMA_GROUPS.md) — agrupación por familia y
  colisiones de nombre.
