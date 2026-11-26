# Contexto genómico y regulones para asignar o priorizar factores sigma

> **Documento de revisión bibliográfica para `dna_poc`.**
> Ruta: `docs/literature/SIGMA_CONTEXT.md`
> Fecha de corte y verificación: **2026-11-26**
> Autor: subagente `lit-sigma-context` (revisión bibliográfica; **no modifica código**)
> Estado: **v1.0** — cada referencia de la tabla final fue resuelta contra Crossref,
> PubMed/PMC o la página del editor. No se incluye ninguna referencia no verificada.
> Ver Anexo A.
>
> **Tesis en una frase:** el contexto genómico y los regulones sirven para
> **priorizar** hipótesis de factor sigma, pero **no son prueba** de que un
> promotor sea reconocido por un sigma. Su valor depende de que la fuente sea
> independiente del gold, de que la cobertura sea suficiente y de que la ganancia
> sobreviva a una ablación con procedencia controlada.

---

## Convenciones

- `[Rxx]`: identificador de la referencia en la tabla de §13. Cada identificador
  tiene PMID, DOI o URL verificado; la tabla indica cuál se comprobó.
- `[PMID:<número>]`: artículo indexado en PubMed; el PMID fue resuelto y el
  título comparado con la fuente.
- `[DOI:<registro>]`: DOI resuelto en Crossref; título, autores y año comparados.
- `[URL:<url>]`: recurso sin DOI (JMLR, PMLR); página resuelta y título comparado.
- Las **observaciones internas** del proyecto se marcan como tales y **no** se
  presentan como literatura. Provienen de
  [`docs/SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md), [`docs/SIGMA_GOLD.md`](../SIGMA_GOLD.md)
  y [`docs/SIGMA_AUDIT.md`](../SIGMA_AUDIT.md).
- «Sigma» se usa para el factor; «sigmulón» para el conjunto de genes/operones
  transcritos por RNAP con ese sigma; «regulón» para el conjunto regulado por un
  factor de transcripción (TF). La distinción es de §2.8 y no es cosmética.
- «Contexto» se usa en sentido amplio: todo lo que no es el motivo del promotor
  evaluado (operón, vecindad, regulón, filogenia, co-expresión, función, condición).
  La arquitectura del promotor se trata aparte porque es secuencia, no contexto
  genómico; se incluye solo para comparar y para evitar confundir canales.

Este documento no es una revisión sistemática PRISMA. Es una revisión narrativa
con trazabilidad de fuentes, orientada a decisiones de diseño y de evaluación.
La estrategia de búsqueda y sus límites están en el Anexo B.

---

## 0. Resumen ejecutivo

1. **Las bases de datos curan literatura, no realizan experimentos.** RegulonDB,
   EcoCyc, SubtiWiki, PRODORIC y DBTBS compilan evidencia experimental y
   computacional con distintos códigos y niveles de confianza. Un objeto «curado»
   no es una observación directa de unión de sigma; es una interpretación de la
   literatura. [R01] [R02] [R03] [R04] [R05]
2. **RegulonDB y EcoCyc no son fuentes independientes.** La regulación de
   *E. coli* K-12 en ambas bases proviene del mismo equipo de curación de la UNAM
   y comparte códigos de evidencia ECO. Usar una como gold y la otra como feature
   «independiente» es un error de procedencia. [R01] [R02]
3. **El objeto correcto es el promotor específico de sigma, no el TSS ni el gen.**
   Un mismo TSS puede tener promotores distintos para sigmas distintas, y un gen
   puede tener TSS internos y antisentido. EcoCyc modela explícitamente cada TSS
   con sigma como un promotor distinto. [R02] [R07] [R26]
4. **El contexto predice mejor el módulo funcional/condición que la identidad del
   sigma.** «Flagelo», «nitrógeno», «heat shock», «stationary» y «envelope» son
   módulos; el mapa módulo→sigma es específico de especie y está reconfigurado
   por duplicación, pérdida y sustitución de sigmas. [R09] [R10] [R14] [R20] [R25]
5. **La asociación por culpa («guilt by association») funciona por encima del azar
   pero está confundida por multifuncionalidad y por densidad de anotación.** Los
   controles de Gillis y Pavlidis muestran que gran parte de la señal se explica
   por cuántas funciones ya tiene el gen, no por la asociación específica.
   [R35] [R36] [R37]
6. **La cobertura del contexto es limitada y sesgada.** Huynen *et al.* reportan
   ~37 % de genes con conservación de orden, ~6 % con fusión, ~8 % con
   co-ocurrencia en operones y ~11 % con perfiles filogenéticos; la combinación
   llega a ~50 %, con solapamiento. La predicción de operones tiene sensibilidad
   de 30–50 % en *E. coli*. [R38] [R41]
7. **El contexto debe entrar como score/prior con procedencia explícita, no como
   «evidencia de sigma».** La conversión correcta es un log-odds con suavizado,
   contracción, corrección de prevalencia y calibración; un LLR crudo no es una
   probabilidad. [R60] [R61] [R62] [R63]
8. **La ablación debe comparar sobre el mismo subconjunto, con procedencia
   controlada y con controles negativos.** Una mejora que desaparece al quitar
   features derivadas del gold, o que también aparece con una feature aleatoria,
   es sesgo de anotación, no señal biológica. [R53] [R55] [R66] [R69] [R70]
9. **La circularidad tiene al menos ocho formas en este problema:** co-procedencia
   gold/feature, fuga por anotación funcional, fuga por motivo, fuga por condición,
   fuga por ortología, negativos muestreados de la etiqueta, ajuste de
   hiperparámetros en test y fuga temporal entre versiones de la base. [R53] [R54]
   [R57] [R58]
10. **En los datos del proyecto, el contexto no mejora K-12 fuera de muestra y sí
    mejora *B. subtilis*, pero la mejora comparte procedencia con el gold.** Esto
    es precisamente el patrón que exige un *provenance holdout* antes de afirmar
    que el contexto «captura biología». Observación interna, no literatura:
    [`docs/SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md) §5–§8.

### 0.1 Qué se puede y qué no se puede afirmar

| Afirmación | ¿Soportada por la literatura? | Forma correcta |
|---|---|---|
| «Los genes flagelares tardíos suelen ser σ28 en *E. coli*» | Sí, como tendencia del módulo | «El módulo flagelar tardío está enriquecido en promotores σ28; no todos sus TSS son σ28» [R15] [R16] |
| «Si el operón es flagelar, el TSS es σ28» | No | «La pertenencia al módulo eleva el prior de σ28; es una hipótesis, no una asignación» |
| «El regulón de σ54 es el regulón de nitrógeno» | No | «σ54 controla nitrógeno, pero también motilidad, virulencia, biofilm y envelope» [R20] |
| «La base dice σ32, entonces es σ32» | No | «La base lo cura como σ32 con tal evidencia; la etiqueta es una interpretación» [R01] [R02] |
| «El contexto mejora el top-1, entonces es causal» | No | «El contexto mejora la predicción en esta población y anotación; la causalidad requiere intervención» [R59] [R73] |
| «El modelo calibrado da 0,8, entonces hay 80 % de probabilidad de sigma» | Solo si está calibrado y la prevalencia es la de despliegue | «El score calibrado estima 0,8 bajo el prior y el dominio declarados» [R61] [R62] |

---

## 1. Qué significa «contexto» y qué no

### 1.1 Tres objetos que se confunden

En este problema hay tres objetos distintos:

1. **TSS**: la posición genómica donde empieza la transcripción. Es un dato
   experimental (5' RACE, CAGE, TSS-seq, RNA-seq de extremos) o una predicción.
2. **Promotor**: el segmento de DNA que recluta un holoenzima RNAP–sigma. Un TSS
   puede tener varios promotores solapados, cada uno con su sigma y su condición.
   EcoCyc separa explícitamente los promotores por sigma aunque compartan TSS.
   [R02] [R07]
3. **Asignación de sigma**: la hipótesis de qué sigma reconoce ese promotor en
   esa condición. Es un objeto probabilístico, no una propiedad intrínseca del
   gen ni del TSS.

El contexto genómico es útil para el objeto 3 **solo en la medida en que se
refiera al promotor concreto**. Un operón puede tener un TSS primario σ70 y TSS
internos σ38 o σ28; asignar el sigma del operón al gen entero es una
aproximación de granularidad gruesa. [R16] [R26]

### 1.2 Niveles de contexto

| Nivel | Ejemplo | Granularidad | Qué puede apoyar |
|---|---|---|---|
| Secuencia | −24/−12 de σ54, −10 extendida | nt | Reconocimiento por el sigma; es el canal más proximal |
| Gen | función EC/GO, dominio Pfam | gen | Módulo funcional; no el sigma concreto |
| Operón | posición del TSS, tamaño, líder | operón/TSS | Co-transcripción; prior débil sobre el TSS primario |
| Regulón | pertenencia a sigmulón/regulón | gen/operón | Co-regulación; alto riesgo de circularidad |
| Genoma | perfiles filogenéticos, sintenia | genoma/familia | Conservación del módulo; no la identidad del sigma |
| Condición | heat shock, stationary, nitrógeno | experimento/población | Estado fisiológico que activa el sigma |
| Población | co-expresión en compendio | gen/conjunto | Módulo; confundido por condición y crecimiento |

### 1.3 Definición operativa de feature contextual

Una **feature contextual** es una variable calculada para un TSS/operón/gen que
no usa la etiqueta de sigma de ese TSS en la partición de evaluación. Esta
definición excluye:

- features calculadas a partir del gold (p. ej., «está en el sigmulón σ38 de
  RegulonDB») cuando el gold es ese mismo sigmulón;
- features calculadas con un modelo entrenado sobre las mismas etiquetas;
- features que codifican la fuente de anotación (p. ej., «tiene número EC»)
  porque miden densidad de anotación, no biología.

La definición **no** excluye features derivadas de la misma base si se usan para
anotar un genoma nuevo y se evalúan contra otra evidencia. El problema no es la
base en sí, sino la **co-procedencia entre feature y etiqueta** en la evaluación.

### 1.4 Qué autoriza cada nivel

- La secuencia autoriza una hipótesis mecanística de reconocimiento; sigue
  necesitando prueba de unión/transcripción.
- La función/operón autoriza una hipótesis de módulo; no autoriza una identidad
  de sigma.
- El regulón autoriza una hipótesis de co-regulación; si viene de la misma base
  que el gold, no autoriza una evaluación independiente.
- La condición autoriza una hipótesis de activación; no dice qué promotor usa
  el sigma dentro de esa condición.
- La filogenia autoriza una hipótesis de conservación; no dice que el sigma sea
  el mismo ortólogo ni que el sitio esté conservado funcionalmente.

---

## 2. Fuentes de contexto: bases de datos

### 2.1 RegulonDB v12.0

RegulonDB es la base de referencia de la regulación de la transcripción en
*E. coli* K-12. La versión 12.0 describe la reconstrucción de su infraestructura,
la actualización de códigos de evidencia con la ontología ECO y la existencia de
cinco *datamarts*: gen, operón, regulón, **sigmulón** y unidades
GENSOR. [R01]

Dos puntos son críticos para este documento:

- **Tiene un objeto sigmulón**: se puede consultar el conjunto de genes/operones
  asociados a un sigma, no solo el regulón de un TF. [R01]
- **Tiene niveles de confianza weak/strong/confirmed** y distingue evidencia
  física (p. ej., footprinting) de predicción computacional. La propia base
  advierte que se pueden filtrar las interacciones por método para **evitar
  circularidad en benchmarks**. [R01]

La advertencia de RegulonDB es la mejor formulación del problema: si el benchmark
usa las mismas interacciones que la base usa para curar, el resultado es
circular. [R01]

### 2.2 EcoCyc

EcoCyc es la base de organismo modelo de *E. coli* K-12; recopila y resume datos
experimentales y se actualiza por curación manual y por adición programática de
resultados de alto rendimiento. Usa códigos de evidencia detallados (p. ej., IDA
para ensayo directo, IMP para fenotipo mutante) y modela cada TSS con sigma como
un promotor distinto, porque cada uno puede tener regulación diferente. [R02]

EcoCyc y RegulonDB comparten el equipo de curación de la UNAM y los códigos de
evidencia; la propia RegulonDB lo declara. [R01] Por tanto, **EcoCyc no es una
fuente independiente de RegulonDB** para evaluar un modelo entrenado con
RegulonDB.

### 2.3 SubtiWiki

SubtiWiki es la base de *Bacillus subtilis* 168. Contiene información curada de
genes/proteínas, interacciones, elementos reguladores, datos de expresión, rutas
metabólicas y regulones; asigna cada gen a una o más categorías funcionales y a
regulones, y su navegador de regulación incluye interacciones de sigma
(gris). [R03]

Para el problema: SubtiWiki es una fuente de contexto funcional y de red, pero
la asignación de un gen a un regulón es una anotación curada, no una medición de
unión de sigma en ese promotor. [R03]

### 2.4 PRODORIC

PRODORIC es una de las mayores colecciones de sitios de unión de TF procariotas.
La versión 2022 reanotó todo el conjunto y distingue **tres tipos de evidencia
experimental** por sitio: evidencia de expresión *in vivo*, evidencia física de
unión proteína–DNA y evidencia de variación del sitio de unión. También incorporó
efectores, tipo de regulación y 109 nuevos TF. [R04]

Esa taxonomía es directamente reutilizable para ordenar features por cercanía
causal: un sitio con evidencia física está más cerca de «sigma une» que uno con
evidencia de expresión; uno con variación del sitio está más cerca de «el sitio
es funcional» que uno solo predicho. [R04]

### 2.5 DBTBS

DBTBS es una base de referencia de la regulación en *B. subtilis*: TF
caracterizados experimentalmente, sus secuencias de reconocimiento y los genes
que regulan. La versión 2008 contiene 1.475 promotores anotados, 120 TF, 736
operones regulados y 463 operones validados experimentalmente; además incluye
conservación de secuencias upstream en 40 especies Gram-positivas. [R05]

DBTBS es valiosa para contexto comparado, pero es una base antigua (2008) y su
contenido no cubre la literatura posterior. Para *B. subtilis* conviene
combinarla con SubtiWiki y con la anotación de BSGatlas usada por el proyecto;
no son fuentes independientes entre sí. [R03] [R05]

### 2.6 RegPrecise 3.0

RegPrecise es un recurso de regulones reconstruidos por **genómica comparativa**:
781 *regulogs* de TF inferidos en más de 160 genomas, más colecciones de
riboswitches. Cada regulón se caracteriza por un TF, un motivo predicho y un
conjunto de genes/operones con sitios de unión. [R06]

RegPrecise es explícita sobre su límite: las estructuras de operón no se han
estudiado y los regulones propagados son preliminares; además, las evidencias
comparativas no son evidencia experimental. [R06] Es un buen ejemplo de contexto
**inferido** que puede priorizar, pero no probar.

### 2.7 Tabla comparativa

| Base | Alcance | Granularidad | Modelo de evidencia | Aporte para sigma | Riesgo principal | Ref. |
|---|---|---|---|---|---|---|
| RegulonDB v12 | *E. coli* K-12 | gen, operón, regulón, sigmulón, GENSOR | ECO; weak/strong/confirmed; física vs computacional | Sigmulón y promotores curados | Co-procedencia con EcoCyc; curaduría, no experimento | [R01] |
| EcoCyc | *E. coli* K-12 | gen, promotor por sigma, regulón | Códigos IDA/IMP/etc.; curación manual | Promotor específico de sigma | Misma procedencia que RegulonDB | [R02] |
| SubtiWiki | *B. subtilis* 168 | gen, operón, regulón, interacción | Curación manual + datos de expresión | Red y módulos; sigma en el navegador | Anotación curada, no unión por TSS | [R03] |
| PRODORIC | Procariotas (múltiples) | TF, sitio, promotor, operón | 3 tipos de evidencia experimental | Sitios con evidencia física vs expresión | Cobertura desigual entre organismos | [R04] |
| DBTBS | *B. subtilis* | promotor, TF, operón | Experimental curada | Promotores y operones validados | Base 2008; solapamiento con SubtiWiki | [R05] |
| RegPrecise 3.0 | >160 genomas | regulón/regulog, TF | Genómica comparativa | Regulones inferidos y motivos | Inferencia, no experimento; operones no estudiados | [R06] |

### 2.8 Regulón, sigmulón y operón no son sinónimos

- **Operón**: genes co-transcritos desde un mismo promotor. Es una unidad de
  transcripción; su predicción tiene incertidumbre de límites. [R41] [R42] [R43]
- **Regulón**: genes regulados por un mismo TF. Es una unidad de control; puede
  incluir operones distintos y no dice qué sigma usa cada promotor.
- **Sigmulón**: genes/operones transcritos por RNAP con un sigma dado. Es la
  unidad relevante para este problema, pero su definición depende de la evidencia
  y de la condición. [R01]
- **Regulón de sigma** se usa a veces como sinónimo de sigmulón; conviene
  evitarlo porque «regulón» ya está ocupado por el TF.

Un error frecuente es usar la pertenencia a un regulón de TF como si fuera
pertenencia a un sigmulón. Un gen puede estar en el regulón de CRP y ser
transcrito por σ70 o σ38 según la condición; la pertenencia a CRP no resuelve el
sigma. [R23]

### 2.9 Lo que las bases no son

- No son un conjunto negativo: la ausencia de un gen en un sigmulón no significa
  que no haya promotor de ese sigma; solo que la base no tiene la anotación.
- No son ground truth experimental: son compilaciones curadas con evidencia
  heterogénea. [R01] [R02] [R51]
- No son independientes entre sí: RegulonDB/EcoCyc comparten curación; DBTBS/
  SubtiWiki/BSGatlas se solapan; PRODORIC y RegPrecise reutilizan fuentes. [R01]
  [R03] [R05] [R06]
- No son completas: la anotación de regulación está sesgada hacia genes y
  condiciones bien estudiados. [R50] [R51]

---

## 3. Canales de evidencia contextual

### 3.1 Operón y vecindad génica

**Idea.** Si un TSS inicia un operón cuyos genes tienen una función común, el
TSS probablemente comparte el programa regulatorio de esa función. Si el operón
es flagelar, el prior de σ28 sube; si es de asimilación de nitrógeno, el prior de
σ54 sube.

**Qué mide realmente.** Co-transcripción y co-regulación funcional, no unión de
sigma. La posición del TSS en el operón (primero vs interno) y la longitud del
líder son proxies de si el TSS es el promotor primario o un TSS interno.

**Evidencia.**
- La organización en operones es real y predictible: Salgado *et al.* reportan
  ~88 % de acierto en pares adyacentes y ~75 % de unidades de transcripción en
  *E. coli*, con 630–700 operones estimados. [R44]
- La predicción de operones tiene sensibilidad limitada: Ermolaeva *et al.*
  reportan 30–50 % de sensibilidad en *E. coli*, aunque alta precisión en pares
  conservados. [R41]
- Dam *et al.* muestran que la accuracy depende de la distancia intergénica y de
  si se conocen operones del genoma: 90,2/93,7 % con información previa,
  84,6/83,3 % sin ella. [R42]
- Brouwer *et al.* advierten que las estimaciones de performance no son
  comparables entre métodos y que queda margen de mejora. [R43]
- La conservación de orden de genes cubre ~37 % de los genes en el benchmark de
  Huynen *et al.*; la combinación con otras señales llega a ~50 %. [R38]

**Lectura para sigma.** Es un canal **correlacional y de granularidad gruesa**.
Sirve para elevar el prior de un módulo; no para decidir el sigma de un TSS
interno. La incertidumbre de los límites del operón se propaga al prior.

### 3.2 Regulón y red reguladora

**Idea.** Si un gen pertenece a un sigmulón curado, su TSS es candidato de ese
sigma. Si un TF que controla el operón co-ocurre con un sigma, el prior sube.

**Qué mide realmente.** Pertenencia a una anotación de co-regulación. En el caso
del sigmulón, es casi la etiqueta que se quiere predecir; en el caso del regulón
de TF, es un contexto indirecto.

**Evidencia.**
- RegulonDB tiene un datamart de sigmulón y niveles de confianza; la propia base
  recomienda filtrar por método para evitar circularidad. [R01]
- SubtiWiki asigna genes a regulones y muestra interacciones de sigma. [R03]
- RegPrecise reconstruye regulones por comparación de sitios en genomas
  relacionados; son inferencias, no experimentos. [R06]
- La red flagelar de *E. coli* se ha mapeado con ChIP-seq/RNA-seq: FlhDC y FliA
  tienen regulones definidos, y FliA tiene 52 sitios de unión, 30 de ellos
  intragenéticos. [R16]
- El regulón de σE se predijo combinando modelo de promotor y genomas
  relacionados: 89 unidades de transcripción en *E. coli* y ocho genomas
  relacionados, con una parte conservada y otra variable. [R14]
- El regulón de σ32 se triplicó en genes validados combinando expresión y
  validación de promotores; el 12 % de los promotores σ32 también son iniciados
  por σ70. [R21]
- El regulón de σS tiene solapamiento extenso con otros regulones globales; Lrp
  afecta la selectividad σS/σ70 de muchos promotores. [R23]
- La integración probabilística de evidencias en redes funcionales [R34] y en
  STRING [R39] [R40] muestra que contexto, co-expresión y ortología se pueden
  combinar; también muestra que la transferencia entre organismos arrastra la
  anotación de origen.

**Lectura para sigma.** Es el canal de **mayor riesgo circular** si el sigmulón
viene de la misma base que el gold. Como contexto independiente (p. ej., un
regulón mapeado por ChIP en otra condición), es útil pero sigue siendo
correlacional: un sitio puede estar en el sigmulón y no ser funcional en la
condición evaluada.

### 3.3 Perfiles filogenéticos y conservación

**Idea.** Genes que funcionan juntos tienden a conservarse o perderse juntos; un
perfil de presencia/ausencia similar sugiere función relacionada. La conservación
de la vecindad (sintenia) y de motivos upstream sugiere que el módulo y su
regulación se conservan.

**Qué mide realmente.** Co-evolución funcional y conservación de módulo, no
identidad de sigma. Un módulo conservado puede cambiar de sigma en un linaje.

**Evidencia.**
- Pellegrini *et al.* introducen los perfiles filogenéticos y muestran que
  perfiles similares tienden a estar funcionalmente ligados. [R31]
- Kensche *et al.* revisan el método y advierten que solo los pares de mayor
  score son de alta confianza y que hay que corregir el sesgo filogenético y
  combinar con otras evidencias. [R33]
- Huynen *et al.* comparan fusión, orden de genes y perfiles: la conservación de
  orden tiene la mayor cobertura (~37 %), y las interacciones son más fuertes
  cuanto más estricto es el requisito de vecindad, con menos falsos positivos.
  [R38]
- La fusión de genes es otra señal de acoplamiento funcional: dos proteínas que
  se fusionan en otro organismo probablemente interactúan o comparten función.
  [R30]
- La no independencia filogenética es un problema estadístico conocido: los
  linajes no son réplicas independientes. [R68]

**Lectura para sigma.** Es un canal **correlacional y de baja cobertura**. Sirve
para sostener que un módulo es antiguo/conservado, no para asignar el sigma de un
promotor. La conservación de un motivo es más proximal, pero puede reflejar
conservación de la secuencia por otras razones (estructura, GC, sesgo mutacional).

### 3.4 Co-expresión

**Idea.** Genes con perfiles de expresión similares a través de condiciones
suelen compartir regulación. Si un TSS se co-expresa con targets conocidos de
σ32, el prior de σ32 sube.

**Qué mide realmente.** Correlación fenotípica condicionada por el compendio de
condiciones. Está confundida por crecimiento, disponibilidad de nutrientes,
estrés global y estructura de operón.

**Evidencia.**
- Eisen *et al.* muestran que el clustering de expresión agrupa genes de función
  similar y que la co-expresión con genes conocidos da pistas funcionales. [R32]
- Nicolas *et al.*, en *B. subtilis*, mapean 2.935 promotores y los agrupan en
  regulones controlados por sigma, explicando ~66 % de la varianza de actividad
  transcripcional; también muestran que una fracción importante de RNAs
  antisentido puede venir de iniciación espuria por sigmas alternativas. [R26]
- Weber *et al.* usan expresión genome-wide para definir el regulón σS y
  encuentran solapamiento con otros regulones globales. [R23]

**Lectura para sigma.** Es un canal **correlacional** y muy dependiente del
compendio. Un compendio con la misma condición que definió el regulón produce
fuga; un compendio independiente y balanceado es más valioso. La co-expresión
puede preceder a la asignación de sigma, pero no la prueba.

### 3.5 Contexto funcional (EC/GO/ruta)

**Idea.** La función del producto génico identifica el módulo fisiológico
(flagelo, nitrógeno, chaperonas, envelope, stationary) y el módulo tiene una
sigma característica.

**Qué mide realmente.** Función anotada, que a su vez puede ser inferida por
homología, dominio o texto. Es un canal **de anotación**, no de medición.

**Evidencia.**
- Los módulos funcionales tienen sigmas características: flagelo/σ28, nitrógeno/
  σ54, heat shock/σ32, stationary/σ38, envelope/ECF. [R12] [R15] [R17] [R21]
  [R23] [R24]
- Pero la asociación no es exclusiva: σ54 también controla motilidad, virulencia,
  biofilm y envelope [R20]; σ32 tiene un componente de membrana y regula
  reguladores globales [R21]; σS se solapa con σ70 [R23]; las ECF son diversas
  [R12] [R13].
- La anotación funcional tiene errores sustanciales en bases públicas: Schnoes
  *et al.* reportan 5–63 % de misanotación promedio en seis superfamilias, y
  >80 % en 10 de 37 familias. [R50]

**Lectura para sigma.** Es el canal más interpretable y el más tentador, pero
también el más expuesto a **sesgo de anotación**: si la función se infirió por
homología a un gen ya anotado como flagelar, y el sigma se curó por esa misma
lógica, la feature y la etiqueta comparten cadena de inferencia. [R50] [R52]

### 3.6 Arquitectura del promotor (secuencia, no contexto genómico)

**Idea.** El sigma reconoce determinantes de secuencia: −35/−10 para σ70, −24/
−12 para σ54, −10 extendida, longitud del espaciador, discriminador, UP element.

**Qué mide realmente.** Compatibilidad de secuencia con el modelo de promotor de
un sigma. Es el canal más proximal, pero no es suficiente: un motivo puede estar
presente y no ser funcional, o puede ser funcional con una variante no capturada
por el PWM.

**Evidencia.**
- Barrios *et al.* compilan 186 promotores σ54 y actualizan el consenso −24/−12;
  el consenso extendido aumenta la probabilidad de identificar promotores
  genuinos. [R19]
- σ54 requiere un activador bEBP y ATP para remodelar RNAP; la presencia del
  motivo no implica transcripción. [R18]
- Nonaka *et al.* muestran que σ70 inicia desde el 12 % de los promotores σ32,
  lo que ilustra el solapamiento de arquitecturas. [R21]
- La predicción de promotores σE requirió un modelo específico y validación
  comparativa; la variabilidad de secuencia y espaciado dificulta la
  predicción. [R14]

**Lectura para sigma.** Es un canal **mecanísticamente proximal** pero no
suficiente. Debe tratarse como likelihood de secuencia, no como contexto; y debe
calibrarse contra un nulo (no usar el score crudo como probabilidad). [R60]

### 3.7 Condición experimental

**Idea.** El sigma se activa en una condición: heat shock, stationary, limitación
de nitrógeno, estrés de envelope. Si la condición del TSS es conocida, el prior
del sigma activo sube.

**Qué mide realmente.** Estado fisiológico poblacional. No dice qué promotor usa
el sigma dentro de esa condición, ni si el TSS es el que se activa.

**Evidencia.**
- σS se induce al entrar en stationary y en muchas otras condiciones de estrés;
  su control es transcripcional, traduccional y por proteólisis. [R24]
- La respuesta a heat shock de *E. coli* está mediada por σ32 y tiene un bucle de
  retroalimentación por chaperonas. [R21] [R22]
- La limitación de nitrógeno controla ~100 genes de la respuesta Ntr. [R17]
- Las ECF responden a señales de envelope y suelen estar reguladas por
  anti-sigma. [R12] [R13]
- En *B. subtilis*, la respuesta general de estrés está gobernada por σB, no por
  un ortólogo directo de σ32 o σ38. [R25]

**Lectura para sigma.** Es un canal **casi causal** para la actividad del sigma
(intervenir la condición cambia la actividad), pero **correlacional** para el
promotor concreto. La condición es poblacional; la heterogeneidad celular y la
dinámica pueden hacer que un TSS no se use aunque el sigma esté activo.

### 3.8 Tabla maestra evidencia/feature/fuente

> **Cómo leer:** «Fuerza» es la fuerza esperada de la asociación con sigma, no la
> fuerza de la prueba. «Independencia» es respecto del gold de evaluación.
> «Causal/Corr.» distingue si la feature puede, en principio, intervenir en el
> mecanismo (causal) o solo asociarse (correlacional). Ninguna feature de esta
> tabla, por sí sola, prueba que un promotor sea de un sigma.

| # | Feature | Qué mide | Fuente típica | Evidencia | Fuerza | Causal/Corr. | Independencia | Riesgo circular |
|---|---|---|---|---|---|---|---|---|
| F01 | `operon_position` (primero/interno, distancia al gen) | Si el TSS es el promotor primario del operón | Anotación de operones; distancia intergénica | [R41] [R42] [R44] | Media | Correlacional (co-transcripción) | Baja si la anotación viene del gold | Alto |
| F02 | `operon_size` | Número de genes co-transcritos | Anotación de operones | [R42] [R43] | Débil | Correlacional (proxy de anotación) | Baja | Alto |
| F03 | `leader_length` / 5' UTR | Tipo de TSS (primario vs interno) | Coordenadas TSS + gen | [R44] | Débil-media | Correlacional | Media | Medio |
| F04 | `gene_neighborhood` (±5 kb) | Funciones de genes vecinos | Anotación funcional | [R29] [R38] | Media | Correlacional | Media | Medio |
| F05 | `synteny_conservation` | Conservación de vecindad entre genomas | Genómica comparada | [R29] [R33] [R38] | Media | Correlacional | Media | Medio |
| F06 | `regulon_membership` | Pertenencia a regulón/sigmulón curado | RegulonDB/SubtiWiki/RegPrecise | [R01] [R03] [R06] | Fuerte | Correlacional | **Baja** | **Muy alto** |
| F07 | `tf_overlap` | TFs que regulan el operón; co-ocurrencia TF–sigma | RegulonDB/DBTBS | [R01] [R05] | Media | Correlacional | Baja | Alto |
| F08 | `phylogenetic_profile` | Presencia/ausencia de ortólogos | OrthoFinder/eggNOG/OrthoDB | [R31] [R33] | Media-débil | Correlacional | Media | Medio |
| F09 | `module_conservation` | Conservación del módulo funcional | Genómica comparada | [R14] [R20] | Media | Correlacional | Media | Medio |
| F10 | `coexpression` | Correlación de expresión con targets conocidos | RNA-seq compendio | [R32] [R26] | Media | Correlacional | Media (si el compendio es independiente) | Alto si comparte condición |
| F11 | `condition_metadata` | Heat shock, stationary, nitrógeno, envelope | Metadatos experimentales | [R17] [R21] [R24] | Fuerte | Casi causal (activa el sigma) | Alta | Bajo-medio |
| F12 | `functional_annotation` | EC/GO/KEGG/UniProt keywords | Anotación automática/curada | [R50] [R51] | Media | Correlacional | Media | Alto |
| F13 | `domain_architecture` | Dominios Pfam del producto | Pfam/InterPro | [R50] | Débil | Correlacional | Alta | Bajo |
| F14 | `promoter_architecture` | −24/−12, −35/−10, extended −10, espaciador, UP | Secuencia | [R19] [R14] [R21] | Fuerte | Causal-mecanística (reconocimiento) | Alta | Bajo (si no se deriva del gold) |
| F15 | `motif_conservation` | Motivo conservado upstream de ortólogos | Alineamiento + PWM | [R14] [R19] | Media-fuerte | Causal-mecanística (sitio) | Media | Medio |
| F16 | `sigma_gene_inventory` | Qué sigmas codifica el genoma | HMM/Pfam + nombre | [R12] [R13] | Fuerte (límite) | Causal (qué es posible) | Alta | Bajo |
| F17 | `expression_tss` | Uso de TSS/promotor medido | TSS-seq/RNA-seq | [R26] | Fuerte | Observacional | Alta (si es independiente) | Bajo-medio |
| F18 | `chip_binding` | Unión física del sigma o TF | ChIP-seq/ChIP-exo | [R16] [R14] | Muy fuerte | Causal-mecanística | Alta (si es independiente) | Bajo-medio |
| F19 | `perturbation` | Knockout/depleción de sigma cambia el TSS | Mutante/CRISPRi | [R23] [R21] | Muy fuerte | Causal (intervención) | Alta | Bajo |
| F20 | `gold_derived` | Cualquier feature tomada del gold | Gold | — | — | **Tautológica** | **Nula** | **Prohibido en evaluación** |

### 3.9 Matriz de madurez

| Canal | Cobertura típica | Costo | Independencia | Uso recomendado |
|---|---|---|---|---|
| Secuencia (F14/F15) | Alta | Bajo | Alta | Likelihood o feature principal; calibrar contra nulo |
| Condición (F11) | Media | Bajo | Alta | Prior de activación; no asignación |
| Operón (F01–F03) | Alta en genomas modelos | Bajo | Baja-media | Prior débil; ablación obligatoria |
| Función (F12/F13) | Media | Bajo | Media | Prior de módulo; controlar densidad de anotación |
| Co-expresión (F10) | Media | Medio | Media | Contexto independiente si el compendio lo es |
| Filogenia (F08/F09) | Media-baja | Medio | Media | Conservación de módulo; corregir filogenia |
| Regulón (F06/F07) | Media | Bajo | Baja | **Solo si es independiente del gold** |
| ChIP/perturbación (F18/F19) | Baja | Alto | Alta | Validación causal de un subconjunto |

---

## 4. Contexto funcional por sigma

Esta sección resume las asociaciones clásicas y, sobre todo, sus excepciones. La
regla general es: **el módulo funcional está más conservado que el nombre del
sigma**.

### 4.1 Flagelo → σ28 (FliA)

En *E. coli* y *Salmonella*, la expresión flagelar es una cascada jerárquica. El
regulador maestro FlhDC activa genes tempranos y la expresión de FliA (σ28); los
genes tardíos de ensamblaje son transcritos principalmente por FliA. La actividad
de σ28 está controlada por el anti-sigma FlgM, que se exporta cuando el hook-basal
body está completo; así la transcripción tardía se acopla al ensamblaje. [R15]

Fitzgerald *et al.* mapean la red flagelar de *E. coli* con ChIP-seq y RNA-seq:
definen un regulón FlhDC más restringido de lo esperado y expanden el regulón
FliA; 30 de los 52 sitios de unión de FliA están dentro de genes. [R16]

**Consecuencia:** «flagelo → σ28» es válido para los genes tardíos (clase III),
no para todo el regulón flagelar. Los genes tempranos son σ70 (con FlhDC) y
algunos componentes pueden ser σ54 en otras bacterias. Un TSS intragenético de
FliA no se detecta mirando solo el promotor del operón. [R15] [R16] [R20]

### 4.2 Nitrógeno → σ54 (RpoN)

La limitación de nitrógeno activa la respuesta Ntr (~100 genes) mediante
reguladores cuya actividad depende de metabolitos. [R17] Los promotores σ54 usan
−24/−12 en lugar de −35/−10; Barrios *et al.* compilan 186 promotores y mejoran
el consenso. [R19] σ54 requiere activadores bEBP que hidrolizan ATP para
remodelar RNAP; sin activador, el motivo no basta. [R18]

Pero σ54 no es solo nitrógeno: análisis comparativos la implican en motilidad,
virulencia, biofilm, LPS, flagelos y envelope, con una distribución dependiente
del filo. [R20] En *B. subtilis*, la regulación de nitrógeno también involucra
GlnR/TnrA y no se reduce a σ54.

**Consecuencia:** «nitrógeno → σ54» es una pista de módulo, no una regla de
asignación. El regulón σ54 es más amplio y heterogéneo que el módulo nitrógeno.
[R20]

### 4.3 Heat shock → σ32 (RpoH)

En *E. coli*, la respuesta a heat shock está mediada por σ32 y un bucle de
retroalimentación por chaperonas. [R22] Nonaka *et al.* triplican los genes
validados como transcritos por σ32, muestran que el regulón codifica numerosos
reguladores globales, que σ70 inicia desde el 12 % de los promotores σ32 y que
una fracción importante del regulón reside en la membrana. [R21]

**Consecuencia:** «heat shock → σ32» es específico de *E. coli* y proteobacterias;
en *B. subtilis* la respuesta general de estrés está gobernada por σB. [R25] No
transferir la regla entre especies sin verificar el inventario de sigmas y la
literatura de la especie.

### 4.4 Stationary → σ38 (RpoS)

σS es el regulador maestro de la respuesta general de estrés en *E. coli* y
bacterias relacionadas; se induce al entrar en stationary y en muchas otras
condiciones, y su control es transcripcional, traduccional y por proteólisis.
[R24] Weber *et al.* estiman que hasta el 10 % de los genes están bajo control
directo o indirecto de σS, identifican un consenso y muestran solapamiento
extenso con otros regulones globales; Lrp afecta la selectividad σS/σ70. [R23]

**Consecuencia:** «stationary → σ38» es una tendencia, no una regla; σS se solapa
con σ70 y su actividad depende de la condición y de la disponibilidad. En
*B. subtilis* el equivalente funcional es σB, no σ38. [R25]

### 4.5 ECF / envelope → σE/σ24 y otras ECF

Las ECF son el grupo más numeroso de sigmas alternativas; suelen ser pequeñas,
divergentes y reguladas por anti-sigma. [R11] [R12] Staroń *et al.* clasifican la
familia en más de 40 grupos mayores, con apoyo en la arquitectura de dominios de
los anti-sigma, la conservación del contexto genómico y los motivos de promotor.
[R13]

El regulón de σE en *E. coli* se predijo con un modelo de promotor validado y
genómica comparada: 89 unidades de transcripción en *E. coli* y ocho genomas
relacionados, con una parte conservada (LPS, porinas) y otra variable
(asociada a patogénesis). [R14]

**Consecuencia:** «envelope → ECF» es una generalización; hay muchas ECF con
funciones distintas y no todas responden a envelope. Además, la nomenclatura es
tramposa: **σE de *E. coli* es la ECF σ24; σE de *B. subtilis* es un sigma de
esporulación (madre), no una ECF.** [R27] [R28] No transferir contexto de
envelope a σE de *B. subtilis*.

### 4.6 Tabla de asociaciones y excepciones

| Módulo | Sigma canónica en *E. coli* | Equivalente funcional en *B. subtilis* | Fuerza | Excepciones y límites | Evidencia |
|---|---|---|---|---|---|
| Flagelo tardío | σ28 (FliA) | SigD (σD) | Alta para clase III | Genes tempranos σ70; sitios intragenéticos; otras bacterias usan σ54 | [R15] [R16] [R20] |
| Nitrógeno | σ54 (RpoN) | SigL (σL) + GlnR/TnrA | Alta para Ntr | σ54 también motilidad, virulencia, biofilm, envelope | [R17] [R18] [R19] [R20] |
| Heat shock | σ32 (RpoH) | σB (estrés general) y otras | Alta en *E. coli* | Especie-específico; σ32 no es universal | [R21] [R22] [R25] |
| Stationary | σ38 (RpoS) | σB (SigB) | Media | Solapamiento con σ70; no todos los genes de stationary | [R23] [R24] [R25] |
| Envelope/ECF | σE/σ24 y otras ECF | σW/σX/σV/σY/σZ; σE es esporulación | Media | ECF diversas; colisión de nombres | [R12] [R13] [R14] [R27] [R28] |
| Esporulación | — | σH → σF/σE → σG/σK | Alta en *B. subtilis* | Cascada de desarrollo; no hay equivalente en *E. coli* | [R27] [R28] |

### 4.7 Por qué la asociación no es una regla

1. **Solapamiento:** el mismo promotor puede ser iniciado por más de un sigma
   (σ70 y σ32; σ70 y σS). [R21] [R23]
2. **Regulones amplios:** σ54 y σ32 controlan procesos que exceden el módulo
   clásico. [R20] [R21]
3. **Sitios intragenéricos/antisentido:** la señal puede estar dentro de genes o
   en transcripción espuria. [R16] [R26]
4. **Especificidad de especie:** el mapa módulo→sigma cambia; σE no es la misma
   σE en *E. coli* y *B. subtilis*. [R27]
5. **Anotación:** la función del gen puede ser inferida y estar mal anotada.
   [R50] [R52]
6. **Condición:** el módulo puede estar presente pero no activo; el sigma puede
   estar activo pero no unir ese promotor. [R24]

---

## 5. Predicción de regulones y sitios de unión

### 5.1 PWM y escaneo de motivos

El enfoque clásico construye una matriz de frecuencias (PWM/PSSM) a partir de
sitios conocidos y escanea el genoma con un umbral. MEME/FIMO son las
herramientas de referencia para descubrir y escanear motivos. [R45] [R46]

Ventajas: interpretable, rápido, funciona con pocos sitios. Límites: la PWM
asume independencia entre posiciones; el umbral es arbitrario; la tasa de falsos
positivos depende del tamaño y composición del genoma; y el score no es una
probabilidad. [R60]

Para σ54, el motivo −24/−12 es relativamente informativo; aun así, la presencia
no implica activación, porque se requiere bEBP. [R18] [R19] Para σ28 y σ32, la
degeneración y el solapamiento con σ70 hacen que el motivo solo sea un componente.

### 5.2 Genómica comparada y *footprinting* filogenético

La idea es que un sitio funcional está conservado en ortólogos, mientras que un
sitio espurio no. Rhodius *et al.* aplican un modelo de promotor σE validado y
predicen 89 unidades de transcripción en *E. coli* y ocho genomas relacionados;
la parte conservada es funcionalmente coherente y la variable es
organismo-específica. [R14]

Ventajas: reduce falsos positivos; permite transferir un regulón entre especies
cercanas. Límites: requiere ortología y alineamiento; los regulones pueden
reconfigurarse; un sitio conservado puede ser no funcional; la conservación de
secuencia no implica conservación de sigma. [R14] [R33]

### 5.3 Regulones integrativos

Los regulones mejor caracterizados combinan varias evidencias:

- **σE:** modelo de promotor + comparación de genomas. [R14]
- **σ32:** expresión genome-wide + validación de promotores + bibliotecas de
  sobreexpresión. [R21]
- **σS:** perfil de expresión + consenso + solapamiento con otros regulones. [R23]
- **Flagelar:** ChIP-seq + RNA-seq de FlhDC y FliA. [R16]
- **σ54:** compilación de promotores + comparación de genomas + activadores.
  [R19] [R20]
- **σB de *B. subtilis*:** microarrays + HMM de reconocimiento + mapeo de
  5'. [R25]
- **Transcriptoma condicional:** 2.935 promotores agrupados en regulones de
  sigma. [R26]

La lección común: **la integración mejora la cobertura y la confianza, pero cada
fuente aporta su propio sesgo**. Un regulón integrado no es una verdad
independiente; es un consenso de evidencias. [R14] [R21] [R23] [R26]

### 5.4 Herramientas

| Herramienta | Qué hace | Entrada | Evaluación | Límite para sigma | Ref. |
|---|---|---|---|---|---|
| MEME/FIMO | Descubre motivos y escanea secuencias | Sitios conocidos o genoma | Descubrimiento + scanning | Motivo ≠ sigma; umbral arbitrario | [R45] [R46] |
| SigmoID | Identifica promotores, TFBS y terminadores; visualiza contexto; ayuda a corregir anotación | Genoma + bases | Caso de *Pectobacterium*: >900 sitios, ~40 % del cromosoma | Usa RegPrecise/RegulonDB; circularidad si se evalúa con ellas | [R47] |
| iPromoter-2L | Predice promotor y tipo (σ24/σ28/σ32/σ38/σ54/σ70) | Secuencia | Validación cruzada | Solo secuencia; no contexto; clases desbalanceadas | [R48] |
| ProPr54 | Predice promotores y regulón σ54 con CNN/RNN híbrida | Secuencia | Reportada por los autores | Entrenado en promotores σ54 conocidos; validación independiente limitada | [R49] |

### 5.5 Qué mide cada evaluación

| Evaluación | Qué mide | Qué no mide |
|---|---|---|
| Validación cruzada en sitios conocidos | Capacidad de recuperar sitios del mismo conjunto | Generalización a promotores no curados |
| *Holdout* por promotor | Generalización a promotores no vistos | Generalización a otra especie/condición |
| *Holdout* por genoma | Transferencia entre especies | Transferencia a linajes distantes |
| *Holdout* temporal | Novedad respecto de la versión de la base | Causalidad |
| Validación experimental | Unión/actividad real | Coberturagenómica completa |

### 5.6 Límites de la predicción de regulones

- Los conjuntos positivos son pequeños para sigmas raras (σ54: 99 sitios en
  K-12, 11 en *B. subtilis*; σ32: 13 en *B. subtilis*). [Interno: `SIGMA_GOLD.md`]
- Las PWMs no capturan cooperatividad, geometría, nucleoides ni accesibilidad.
- La comparación entre especies requiere ortología correcta y alineamiento.
- Los benchmarks de descubrimiento de motivos dependen de gold standards
  curados; los sesgos del gold se propagan. [R57]
- Los ChIP-seq pueden contener motivos de TF no objetivo; en humanos se reportó
  hasta 45 % de peaks con motivos «zingers». El mismo tipo de confusión puede
  aparecer en bacterias. [R58]
- Un regulón predicho no es un sigmulón validado; y un sigmulón curado no es una
  prueba de unión en cada promotor. [R01] [R06]

---

## 6. Causal vs correlacional

### 6.1 Definiciones operativas

- **Causal**: una intervención sobre el sigma (knockout, depleción, mutación del
  motivo, cambio de condición controlado) cambia el uso del TSS o la
  transcripción del gen.
- **Mecanístico**: el sigma/RNAP une el promotor y lo usa; se prueba con unión
  física y transcripción dependiente de sigma.
- **Correlacional**: la feature y el sigma co-varían en la población de
  entrenamiento, sin intervención.
- **Tautológico**: la feature se deriva de la etiqueta; no hay hipótesis
  independiente.

### 6.2 DAG mínimo

```text
        A (anotación funcional / base de datos)
        │                 │
        ▼                 ▼
        F (feature)      Y (etiqueta de sigma)
        ▲                 ▲
        │                 │
        C (condición) ────┘
        ▲
        │
        S (secuencia del promotor) ──► Y
        │
        └──────────────► F (si la feature incluye motivos)
```

Caminos de confusión:

- `A → F` y `A → Y`: la anotación produce tanto la feature como la etiqueta.
- `C → F` y `C → Y`: la condición produce tanto el módulo como el sigma.
- `Y → F`: si la feature se deriva de la etiqueta (regulón/sigmulón del gold).
- `S → F` y `S → Y`: la secuencia produce tanto el motivo como la asignación.

Para estimar un efecto causal de `F` sobre `Y` habría que bloquear `A` y `C`,
o intervenir `F`. Un modelo observacional no lo hace. [R73] Por eso el contexto
puede **predecir** sin **explicar**.

### 6.3 Tabla de afirmaciones

| Afirmación | Tipo | Qué la sostendría | Qué no la sostiene |
|---|---|---|---|
| «El gen es flagelar, por tanto el TSS es σ28» | Correlacional | Unión de σ28 en ese TSS, mutación del motivo, transcripción dependiente | Anotación flagelar del gen |
| «El operón está en el sigmulón σ38, por tanto es σ38» | Tautológica si el sigmulón es el gold | Evidencia independiente de unión/actividad | La propia anotación del sigmulón |
| «La condición es stationary, por tanto el TSS es σ38» | Casi causal (actividad) | Medición de uso de TSS en stationary + dependencia de σS | El metadato de condición solo |
| «El promotor tiene −24/−12, por tanto es σ54» | Mecanístico parcial | Unión de σ54 + activación por bEBP | El motivo solo |
| «El motivo está conservado, por tanto es funcional» | Correlacional | Ensayo de actividad del sitio, mutagénesis | La conservación sola |
| «El knockout de σ32 baja la expresión, por tanto hay un promotor σ32» | Causal (si el diseño es limpio) | Knockout + TSS mapping + complementación | La expresión global sola |
| «El modelo acierta más con contexto, por tanto el contexto causa el sigma» | No secuencial | Ablación, intervención, *provenance holdout* | El acierto predictivo solo |

### 6.4 Jerarquía de evidencia para sigma

1. **Unión física y transcripción dependiente de sigma** en el mismo promotor:
   EMSA, footprinting, ChIP-seq/ChIP-exo, transcripción in vitro, mutagénesis del
   motivo, complementación. Es el estándar causal.
2. **Perturbación**: knockout/depleción de sigma con cambio de TSS/actividad.
3. **Genética del promotor**: mutación del sitio que abole la actividad y
   revierte con la secuencia canónica.
4. **Perfil de uso de TSS en condiciones** con dependencia de sigma.
5. **Conservación de motivo y contexto** (hipótesis).
6. **Anotación funcional/regulón** (hipótesis débil, alto riesgo circular).

### 6.5 Por qué el contexto no es prueba de sigma

Un LLR alto significa «esta feature es más frecuente en σX que en el fondo», no
«este promotor es σX». La inferencia va del efecto (módulo) a la causa (sigma) y
es una **inferencia inversa**: muchas causas pueden producir el mismo módulo.
[R59] Además, la feature y la etiqueta pueden compartir la misma causa
(anotación), lo que produce asociación sin causalidad. [R73]

La frase correcta es:

> «El contexto eleva el prior de σX para este TSS; es consistente con σX, pero
> no demuestra que σX lo reconozca.»

La frase incorrecta es:

> «El contexto demuestra que es σX.»

---

## 7. De contexto a log-odds calibrado

### 7.1 Notación y objetivo

Sea `Y ∈ {σ1, …, σK}` el sigma de un promotor (multietiqueta en general) y
`X = (x1, …, xM)` un vector de features contextuales. Se busca una cantidad
interpretable y calibrada:

```text
p(s | x) = P(Y = s | X = x)
```

No se busca un score arbitrario que ordene candidatos. La diferencia importa:
un ranking puede ser útil aunque no esté calibrado; una probabilidad mal
calibrada se propaga a cualquier decisión posterior (umbral, abstención,
priorización experimental).

### 7.2 LLR por evidencia

Para una feature categórica `e` con valores `v`, un sigma `s` y un conjunto de
entrenamiento:

```text
P(v | s) = (n_{s,v} + α) / (N_s + α·(|V| + 1))
P(v)     = (n_v   + α) / (N   + α·(|V| + 1))        α = 1 (Dirichlet)

LLR_e(s, v) = log P(v | s) − log P(v)
```

`N_s` es el número de sitios de entrenamiento con sigma `s`; `n_v` es el número
de sitios con valor `v` en el fondo. El fondo puede ser toda la población de
TSS o un conjunto negativo emparejado; la elección cambia el LLR y debe
declararse.

Propiedades:

- Un sitio multietiqueta cuenta para cada sigma verdadero en el numerador y una
  vez en el fondo.
- Si `N_s = 0`, la tabla de ese sigma queda vacía; no se inventa un LLR.
- `α = 1` y `min_support = 5` son constantes declaradas, no ajustadas al test.
- El LLR es una asociación; solo es un likelihood ratio si el entrenamiento es
  representativo del despliegue y las features se miden igual. [R64]

### 7.3 Combinación

La forma naive Bayes es:

```text
score(s | x) = b_s + Σ_e λ_e · LLR_e(s, x_e)
```

donde:

- `b_s = logit(π_s)` si `π_s` es un prior explícito; en otro caso es un sesgo
  aprendido.
- `λ_e` es un factor de contracción/temperatura por evidencia; `λ_e = 0` apaga
  la evidencia y reproduce el baseline.
- `λ_e` se estima en un *inner loop* de entrenamiento, nunca en test.

Si los LLR son likelihood ratios y `b_s` es un prior, `score` es log-odds
posterior. Si los LLR se estiman con las mismas etiquetas que se evalúan, `score`
es un **score discriminativo**; llamarlo «posterior» es incorrecto.

### 7.4 Naive Bayes y dependencia

La suma asume independencia condicional entre evidencias. En este problema las
features están correlacionadas: `operon_position`, `leader` y `operon_size`
miden aspectos del mismo TSS; `functional_annotation` y `gene_neighborhood`
comparten la misma anotación; `regulon_membership` y `tf_overlap` comparten la
misma base.

Consecuencia: la suma sobrecuenta la evidencia. Hay cuatro salidas:

1. **Dropear** features redundantes (medir correlación; V de Cramér, información
   mutua).
2. **Agrupar** en categorías conjuntas (p. ej., módulo × posición) con
   suavizado.
3. **Modelo regularizado** (regresión logística con L2/group lasso) que estima
   los coeficientes conjuntamente.
4. **Modelo jerárquico** con efectos por especie/sigma y *partial pooling*.

La opción 1 es la más barata y la que menos sobreajusta con 2.000–3.000 sitios.
La opción 2 es útil cuando hay una hipótesis mecanística clara. Las opciones 3–4
requieren más datos o más supuestos. [R64]

### 7.5 Prior, score y evidencia

Hay tres objetos que no deben confundirse:

| Objeto | Definición | Ejemplo | ¿Es una probabilidad? |
|---|---|---|---|
| Prior `π_s` | Frecuencia esperada de `s` antes de ver la feature | Censo de sigmas del genoma; prevalencia de la condición | No; es una distribución |
| LLR `LLR_e(s,v)` | Cuánto cambia la odds de `s` al ver `v` | `P(flagelar|σ28)/P(flagelar)` | No; es un ratio |
| Score `score(s|x)` | Suma de prior y LLR ponderados | `logit π_s + Σ λ_e LLR_e` | No, hasta calibrar |
| Posterior `p(s|x)` | Probabilidad calibrada bajo el modelo | `softmax(score/T)` | Sí, si está calibrado |

Un prior **no** se estima con el gold de evaluación. Si se estima con el gold,
es una frecuencia de entrenamiento, no un prior. Un LLR estimado con el gold es
una asociación, no una likelihood. Un score sin calibrar no es una probabilidad.

### 7.6 Contracción y suavizado

Con pocos positivos, los LLR extremos son inestables. Opciones:

```text
LLR_shrunk = (n / (n + k)) · LLR
```

o un prior jerárquico:

```text
β_{e,s} ~ Normal(0, τ_e²)
```

donde `τ_e` se estima con todos los sigmas. El efecto es encoger los LLR de
features con poco soporte hacia 0, que es la decisión conservadora: no mover el
ranking sin evidencia.

Reglas prácticas:

- `min_support = 5` por categoría; las categorías con menos se colapsan a
  «otras» o se neutralizan.
- `α = 1` es un prior débil; no elimina el sesgo de categorías raras.
- Nunca usar la categoría «missing» como señal; `missing → LLR = 0`. Si se
  modela como categoría, el clasificador aprende de qué falta (cobertura de
  anotación), no de lo que se midió. Observación interna: en K-12, modelar
  `none` como categoría producía una subida aparente de +1,1 pp que desaparecía
  al neutralizarla; [`docs/SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md) §2.1.

### 7.7 Calibración

Después de obtener `score`, se ajusta un mapa a probabilidad en un conjunto de
calibración (o *inner fold*):

- **Platt scaling** (binario): `p = 1 / (1 + exp(−(a·score + b)))`.
- **Temperature scaling** (multiclase): `p_s = softmax(score_s / T)`.
- **Isotonic regression**: no paramétrica; requiere más datos y puede
  sobreajustar.

`a`, `b` y `T` se ajustan **solo con datos de calibración**, nunca con el test.
La calibración no cambia el ranking; cambia la confianza. [R63] [R71]

Métricas:

| Métrica | Qué mide | Uso |
|---|---|---|
| Log loss | Calidad probabilística conjunta | Proper score; penaliza confianza errónea |
| Brier | Error cuadrático de probabilidad | Proper score; interpretable |
| ECE | Diferencia media entre confianza y acierto | Diagnóstico de calibración; depende de bins |
| Reliability diagram | Curva confianza–acierto | Inspección visual |
| Top-1 / top-k | Discriminación | Útil para priorización, no para probabilidad |
| Cobertura | Fracción con features suficientes | Debe reportarse junto con la métrica |
| Decisión neta | Utilidad a un umbral | Cuando hay costo asimétrico [R72] |

Los *proper scoring rules* son la referencia para comparar modelos
probabilísticos; el ECE es útil pero no es proper y puede manipularse con bins.
[R62]

### 7.8 Corrección de prevalencia

El gold está enriquecido en sitios de sigma; el genoma de despliegue tiene
muchos más TSS sin sigma anotado. Si el entrenamiento es caso-control con
prevalencia `π_train` y el despliegue tiene `π_deploy`, la corrección de
intercepto es:

```text
logit(p_deploy) = logit(p_train)
                + log(π_deploy / (1 − π_deploy))
                − log(π_train  / (1 − π_train))
```

Para multiclase, se suma `log(π_deploy_s / π_train_s)` al score de cada clase.
Esta es la corrección de King y Zeng para datos de eventos raros. [R61]

Sin corrección, un modelo entrenado en un gold balanceado produce
probabilidades infladas en un genoma real. La corrección asume que los LLR son
transportables; si no lo son, hay que recalibrar en el dominio de despliegue.

### 7.9 Ejemplo numérico ilustrativo

> **Aviso:** estos números son un ejemplo aritmético, no datos del proyecto ni
> estimaciones publicadas. Sirven para mostrar la sensibilidad al prior y al LLR.

Supuesto: prevalencia de σ28 `π = 0,04`; contexto «módulo flagelar» con
`P(flagelar|σ28) = 0,55` y `P(flagelar|no σ28) = 0,11`, luego
`LLR = log(0,55/0,11) = 1,609`; PWM con `LLR = 2,303`.

| Escenario | LLR contexto | LLR PWM | Odds posterior | p(σ28) |
|---|---:|---:|---:|---:|
| Solo prior | 0 | 0 | 0,042 | 0,040 |
| Prior + contexto | 1,609 | 0 | 0,209 | 0,172 |
| Prior + contexto + PWM | 1,609 | 2,303 | 0,480 | 0,676 |
| Contexto sobreestimado (LLR real 0,405) + PWM | 0,405 | 2,303 | 0,625 | 0,385 |
| Dos features flagelares correlacionadas, suma naive 3,218 | 3,218 | 0 | 1,042 | 0,510 |
| Las mismas dos, LLR conjunto 1,900 | 1,900 | 0 | 0,279 | 0,218 |

Lecturas:

- El contexto solo sube σ28 de 4 % a 17 %; no lo «prueba».
- Un PWM fuerte lo sube a 68 %; la combinación es la que decide.
- Si el LLR contextual está inflado por anotación (0,405 en vez de 1,609), la
  probabilidad cae de 68 % a 39 %: **la calibración del contexto importa tanto
  como su signo**.
- Dos features correlacionadas sumadas naive dan 51 %; el LLR conjunto da 22 %.
  La independencia mal asumida produce sobreconfianza.

### 7.10 Reglas de implementación

1. Estimar los LLR solo en entrenamiento, con *cross-fitting* por fold.
2. Suavizar con Dirichlet; colapsar categorías raras; documentar `α` y
   `min_support`.
3. `missing → LLR = 0`; reportar cobertura por feature.
4. Medir correlación entre features y dropear/agrupar las redundantes.
5. Estimar `λ_e` en *inner CV*; `λ_e = 0` debe reproducir el baseline exacto.
6. Calibrar en un fold separado; nunca en test.
7. Corregir prevalencia si el gold es caso-control.
8. Reportar log loss, Brier, ECE, top-1/top-k y cobertura.
9. Reportar por sigma y por estrato de cobertura (no solo el promedio).
10. Abstenerse si la cobertura es baja o las evidencias se contradicen.
11. Guardar la tabla de LLR, la versión de la fuente y la fecha; auditoría.
12. No llamar «probabilidad» a un score sin calibración; no llamar «prior» a una
    frecuencia del gold.

---

## 8. Ablación

### 8.1 Qué pregunta responde cada diseño

| Diseño | Qué varía | Qué prueba | Riesgo |
|---|---|---|---|
| Baseline | Nada | Piso de desempeño | Baseline débil infla la ganancia |
| *Add-one* | +1 feature | Señal marginal | Cobertura distinta |
| *Leave-one-out* | −1 feature | Redundancia | Cobertura distinta |
| Grupal | Bloque de features | Canal completo | Redundancia interna |
| Condicional | Feature dado el resto | Independencia condicional | Colinealidad |
| Permutación | Feature barajada | Nulo empírico | Preserva marginal, no estructura |
| Control negativo | Feature aleatoria o de anotación | Sesgo | Elección del control |
| *Provenance holdout* | Quitar features de la misma fuente | Circularidad | Reduce datos |
| Cross-species | Entrenar/testear en especies distintas | Generalización | Dominio |
| Temporal | Versión antigua/nueva de la base | Novedad | Disponibilidad |

### 8.2 Protocolo mínimo

```text
1. Definir unidad (TSS/operón) y folds (por coordenada, no por sitio).
2. Baseline: solo secuencia (PWM), sin contexto.
3. Para cada fold:
   a. Ajustar LLR con train.
   b. Elegir λ en inner CV.
   c. Predecir test.
4. Agregar predicciones out-of-fold.
5. Repetir para cada configuración (baseline, +feature, −feature, all,
   context-only, random).
6. Comparar sobre el mismo subconjunto (cobertura).
7. Tests pareados + IC bootstrap.
8. Provenance holdout + label permutation + controles negativos.
```

### 8.3 Cobertura y comparación justa

Un error frecuente es comparar el modelo con contexto sobre todos los TSS
contra el baseline sobre todos los TSS, cuando la feature solo existe en el 27 %
de los TSS. La ganancia puede deberse a que ese 27 % es «más fácil», no a la
feature.

Regla: **comparar siempre sobre el subconjunto donde la feature está
disponible**, y reportar el baseline en ese mismo subconjunto. Si el baseline ya
es mejor en ese subconjunto, la feature no aporta.

### 8.4 Permutación y controles negativos

- **Permutación dentro de estratos**: barajar la feature dentro de genoma,
  sigma o bin de cobertura, para preservar la marginal y destruir la asociación.
  Repetir 1.000 veces; el Δ observado debe estar fuera de la distribución nula.
- **Control negativo 1**: feature aleatoria con la misma cardinalidad.
- **Control negativo 2**: longitud del gen, número de anotaciones, número de
  dominios, «tiene EC», «tiene nombre corto». Si estos controles dan una
  ganancia similar, el modelo está midiendo densidad de anotación.
- **Control negativo 3**: etiquetas barajadas. Si el modelo con contexto
  «aprende» con etiquetas barajadas, hay fuga.
- **Control de multifuncionalidad**: predecir sigma solo con el número de
  funciones anotadas del gen; si el contexto no supera ese baseline, la señal
  es multifuncionalidad. [R35] [R36]

### 8.5 Tests y métricas

| Comparación | Test | Reportar |
|---|---|---|
| Top-1 pareado | McNemar | Δ pp, IC 95 %, p |
| Top-3 pareado | McNemar por TSS | Δ pp, IC 95 % |
| Δ top-1 entre modelos | Bootstrap por TSS/genoma | IC 95 %, tamaño de efecto |
| Log loss / Brier | Bootstrap pareado | Δ, IC 95 % |
| Calibración | ECE + reliability | Curva por sigma |
| Múltiples features | Holm/BH | Ajuste por comparaciones |

McNemar es apropiado para clasificadores pareados sobre el mismo test; el
bootstrap por unidad preserva la dependencia. [R65] [R67] No usar un test de
proporciones independientes cuando los modelos se evalúan sobre los mismos TSS.

### 8.6 Poder estadístico

- K-12: 3.314 filas; 1 pp ≈ 33 TSS.
- *B. subtilis*: 2.226 filas; 1 pp ≈ 22 TSS.
- σ54 K-12: 99 positivos; 1 pp ≈ 1 sitio; el IC es enorme.
- σ32 *B. subtilis*: 13 positivos; cualquier conclusión por sigma es frágil.

Por eso hay que reportar IC y no solo el punto. Una mejora de +0,1 pp en K-12
son ~3 TSS: no es una señal robusta.

### 8.7 Ablación no es causalidad

Quitar una feature cambia el modelo, no la biología. Un Δ grande al quitar
`operon_position` significa que el modelo usaba esa feature, no que la posición
del operón cause el sigma. Para causalidad hace falta intervención o un diseño
cuasi-experimental. [R73]

### 8.8 Lectura de resultados

| Patrón | Interpretación | Acción |
|---|---|---|
| Δ top-1 > 0, IC excluye 0, sobrevive provenance | Señal contextual | Mantener, con prior explícito |
| Δ top-1 ≈ 0, Δ log loss > 0 | Solo calibración | Usar en modo confianza, no en top-1 |
| Δ desaparece en provenance holdout | Sesgo de anotación | No usar para evaluar; solo anotar |
| Δ también con feature aleatoria | Sesgo de modelo/cobertura | Revisar diseño |
| Δ solo en una especie | Específico de anotación/especie | No generalizar |
| Δ solo en genes muy anotados | Multifuncionalidad | Controlar densidad |
| Δ negativo consistente | Contexto mal especificado | Apagar (`λ = 0`) |

### 8.9 Ejemplo interno (no literatura)

Números de [`docs/SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md) §5–§6:

| Conjunto | Baseline top-1 | Contexto top-1 | Δ | Lectura |
|---|---:|---:|---:|---|
| K-12 | 56,98 % | 56,98 % | 0,00 pp | El selector elige escala 0 en los 5 folds |
| *B. subtilis* | 29,69 % | 43,40 % | +13,70 pp [11,86; 15,68] | Ganancia fuerte en la población anotada |
| *B. subtilis* con tabla curada | 29,69 % | 46,63 % | +16,94 pp | La tabla comparte procedencia con el gold |

Ablación K-12 a escala 1,0: `operon_position` −0,66 pp, `leader` −0,09,
`operon_size` −0,28, `pathway_ec` +0,06, `conservation` −0,13, todas −0,78.
Ninguna supera el ruido. Ablación *B. subtilis*: `operon_position` +8,18,
`leader` +3,91, `operon_size` +2,16, todas +13,70; sin `operon_position` +6,87.

Lectura: en K-12 no hay señal out-of-fold; en *B. subtilis* sí, pero la
procedencia compartida y la subida con la tabla curada obligan a un
*provenance holdout* antes de llamarla biología. Es exactamente el patrón de
§9.

---

## 9. Circularidad

### 9.1 Definición

Hay circularidad cuando la etiqueta de evaluación y una feature comparten un
ancestro causal o procedural: la misma base, el mismo paper, la misma cadena de
anotación por homología, la misma condición, o la misma versión temporal. El
modelo puede predecir la **anotación** en lugar de la biología.

No es solo «usar el gold como feature». Es cualquier camino por el que la
etiqueta influya en la feature, o por el que ambas dependan de una tercera causa
común (la anotación). [R53] [R54] [R55]

### 9.2 Los diez bucles de este problema

| # | Bucle | Ejemplo | Detección | Mitigación |
|---|---|---|---|---|
| C1 | Co-procedencia de bases | Gold RegulonDB, feature EcoCyc | Auditoría de procedencia; RegulonDB declara que comparten curación | Usar una fuente externa o *holdout* |
| C2 | Sigmulón como feature y gold | Feature «está en el sigmulón σ38» y gold «σ38» | La feature se deriva de la etiqueta | Prohibir en evaluación; usar otra evidencia |
| C3 | Anotación funcional transferida | «Flagelar» inferido por homología a un gen ya flagelar | Trazar la cadena de anotación | Anotación independiente o experimental |
| C4 | Motivo descubierto en el gold | MEME sobre los mismos sitios que se evalúan | Descubrir en train, evaluar en test | *Holdout* por promotor; motivo de otra fuente |
| C5 | PWM de la misma base | PWM de RegulonDB para evaluar RegulonDB | La PWM y el gold comparten sitios | PWM de otra fuente o *leave-one-site-out* |
| C6 | Condición compartida | Compendio de co-expresión con el experimento que definió el regulón | Mismos experimentos en train/test | Compendio independiente; *holdout* por estudio |
| C7 | Ortología con anotación transferida | Función de ortólogo transferida del genoma gold | Trazar la transferencia | Ortología + anotación independiente |
| C8 | Negativos muestreados de la etiqueta | Negativos = genes «no están en el sigmulón» | El negativo es la etiqueta invertida | Negativos por ausencia de evidencia, no por etiqueta |
| C9 | Hiperparámetros en test | Elegir la escala en el barrido de test | El barrido usa predicciones out-of-fold | *Nested CV*; *inner loop* |
| C10 | Fuga temporal | Entrenar y testear con la misma versión de la base | No hay *time split* | Versión antigua → versión nueva [R08] |

### 9.3 Pruebas de detección

| Prueba | Cómo se hace | Qué demuestra |
|---|---|---|
| Auditoría de procedencia | Tabla feature → fuente → versión → paper → ¿cita el gold? | Identifica C1–C7 |
| *Same-source ablation* | Quitar todas las features de la fuente del gold y re-evaluar | Si la ganancia desaparece, era procedencia |
| Permutación de etiquetas | Barajar sigma dentro de estrato; re-entrenar | Si hay ganancia, hay fuga |
| Control negativo | Feature aleatoria / longitud / número de anotaciones | Si hay ganancia, es sesgo |
| Cross-species | Entrenar K-12, testear *B. subtilis* (y viceversa) | Si la ganancia desaparece, es específica de anotación |
| Temporal | Base antigua → base nueva | Si la ganancia desaparece, era novedad de curaduría |
| *Annotation-blind holdout* | Evaluar solo en TSS con evidencia experimental independiente | Reduce C3/C4/C5 |
| Concentración | Ver si la ganancia se concentra en genes muy anotados | Detecta multifuncionalidad [R35] |

### 9.4 Jerarquía de mitigación

1. **Fuente independiente**: la feature no comparte base, paper ni cadena de
   anotación con el gold. Es la única mitigación fuerte.
2. ***Provenance holdout***: se quitan las features de la misma fuente y se
   re-evalúa. Es la prueba mínima aceptable.
3. **Cross-species / temporal**: muestra si la señal es específica de la
   anotación o del linaje.
4. **Permutación / controles negativos**: descarta fugas gruesas y sesgos de
   modelo.
5. **Documentación**: si no se puede mitigar, se declara y se baja la claims.

### 9.5 Ejemplo interno

- **K-12**: RegulonDB y EcoCyc comparten equipo de curación. Una feature de
  EcoCyc no es independiente del gold de RegulonDB. [R01] [R02]
- ***B. subtilis***: el gold y las tablas `genes.tsv`/`gt_operons.tsv` vienen de
  BSGatlas/DBTBS. La ganancia de +13,70 pp y la subida a +16,94 pp con la tabla
  curada son consistentes con dependencia de anotación. [`docs/SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md)
  §6, §8.
- **Bug del 44,79 %**: la corrida usaba PWMs de *E. coli* sobre *B. subtilis*;
  no era un resultado de la PWM propia. Es un recordatorio de que la procedencia
  hay que verificarla en el código, no solo en el paper. [`docs/SIGMA_CONTEXT.md`](../SIGMA_CONTEXT.md)
  §7.

### 9.6 Checklist de circularidad

- [ ] ¿La feature viene de la misma base que el gold?
- [ ] ¿La feature cita el mismo paper que la etiqueta?
- [ ] ¿La feature se calculó usando la etiqueta?
- [ ] ¿La PWM/motivo se descubrió en los sitios de evaluación?
- [ ] ¿El compendio de expresión comparte experimentos con el gold?
- [ ] ¿La anotación funcional se transfirió por homología desde el gold?
- [ ] ¿Los negativos se definieron por ausencia en el gold?
- [ ] ¿La escala/pesos se eligieron en test?
- [ ] ¿Train y test usan la misma versión de la base?
- [ ] ¿La ganancia se concentra en genes muy anotados?

Si alguna respuesta es «sí», la evaluación no es independiente. Se puede usar
para anotar, no para afirmar que el contexto captura biología.

---

## 10. Cuándo el contexto ayuda y cuándo es solo sesgo de anotación

### 10.1 Señales de que ayuda

- **Fuente independiente** del gold (otro paper, otra base, otra modalidad).
- **Cobertura alta** y *missingness* no informativo.
- **Granularidad correcta**: la feature se mide en el TSS/operón, no solo en el
  gen.
- **Proximidad mecanística**: condición, motivo conservado, unión física.
- **Estabilidad entre especies/linajes** o, si es específica, declarada como
  específica.
- **Ganancia sobrevive** a cobertura, *provenance holdout*, cross-species,
  permutación y controles negativos.
- **Mejora proper scores** (log loss/Brier), no solo top-1.
- **Efecto dosis-respuesta**: más evidencia independiente → más confianza, sin
  saltos arbitrarios.

### 10.2 Señales de sesgo de anotación

- **Misma fuente** que el gold, o fuentes que comparten curación (RegulonDB/
  EcoCyc).
- **La feature codifica densidad de anotación**: número de funciones, número de
  dominios, longitud del gen, «tiene EC».
- **Ganancia concentrada** en genes muy estudiados.
- **Feature aleatoria o de control** produce una ganancia similar.
- **La ganancia desaparece** al quitar features de la misma fuente.
- **La ganancia desaparece** en otra especie o en una versión temporal distinta.
- **La feature es la etiqueta disfrazada**: sigmulón, regulón, «está en la base».
- **El modelo aprende la prevalencia** de la base, no la biología: predice σ70
  para casi todo porque es la clase mayoritaria.
- **Solo mejora ECE**, no top-1: puede ser calibración, no señal.

### 10.3 Tabla de decisión

| Condición | Interpretación | Acción |
|---|---|---|
| Fuente independiente + cobertura > 50 % + sobrevive provenance | Contexto útil | Prior con peso calibrado |
| Fuente independiente + cobertura < 20 % | Poca evidencia | Peso 0 o abstención |
| Misma fuente que el gold | Circular | No usar en evaluación; solo anotar |
| Ganancia solo en calibración | Útil para confianza | Modo confianza, no top-1 |
| Ganancia solo en una especie | Específica | Declarar; no generalizar |
| Ganancia solo en genes muy anotados | Multifuncionalidad | Controlar densidad o descartar |
| Control negativo también gana | Sesgo | Revisar features/diseño |
| Permutación de etiquetas no da 0 | Fuga | Bloquear y re-evaluar |

### 10.4 Umbrales operativos propuestos

No son leyes; son política de proyecto para evitar afirmaciones fuertes:

- **Cobertura mínima** para usar una feature en top-1: 20 % de los TSS de
  evaluación; por debajo, peso 0 o abstención.
- **Ganancia mínima** para mantener una feature: Δ top-1 ≥ 1 pp con IC 95 % que
  excluya 0, **y** mejora de log loss, **y** supervivencia a *provenance
  holdout*.
- **Si Δ top-1 < 1 pp pero log loss mejora**: usar en modo confianza, no en
  top-1.
- **Si la ganancia desaparece en cross-species**: declarar «específica de la
  anotación/especie» y no usarla como regla general.
- **Si la feature es de la misma fuente que el gold**: no usarla para evaluar;
  usarla solo para anotar un genoma nuevo y evaluar contra otra evidencia.

### 10.5 Cómo elegir un prior sin circularidad

1. **Prevalencia**: usar un censo independiente de TSS/sigmas, no la frecuencia
   del gold.
2. **Módulo→sigma**: usar literatura experimental de la especie, no la
   anotación que se quiere evaluar.
3. **Condición**: usar metadatos experimentales independientes del gold.
4. **Conservación**: usar alineamientos y ortología, no la etiqueta transferida.
5. **Documentar** la fuente y la fecha; si no es independiente, llamarlo
   «frecuencia de entrenamiento», no «prior».

---

## 11. Límites

### 11.1 Límites de datos

1. **No hay validación experimental en el proyecto**: el gold es RegulonDB/DBTBS,
   no un ensayo de unión. [`docs/SIGMA_GOLD.md`](../SIGMA_GOLD.md)
2. **El gold no tiene negativos**: la ausencia de un TSS no es un negativo; es
   ausencia de anotación. [R01] [R02]
3. **La anotación es incompleta y sesgada** hacia genes/condiciones estudiados.
   [R50] [R51] [R56]
4. **Las bases no son independientes**: RegulonDB/EcoCyc comparten curación;
   DBTBS/SubtiWiki/BSGatlas se solapan. [R01] [R03] [R05]
5. **Las versiones cambian**: una evaluación con la versión 12 no es comparable
   con una versión anterior. [R08]

### 11.2 Límites biológicos

6. **Un TSS puede tener varios promotores/sigmas**; la asignación es
   multietiqueta y dependiente de condición. [R02] [R07] [R21] [R23]
7. **El módulo→sigma no es universal**: σE de *E. coli* ≠ σE de *B. subtilis*;
   σ32/σ38 no tienen el mismo papel en todas las especies. [R25] [R27]
8. **Los regulones se reconfiguran**: la parte conservada y la variable pueden
   diferir entre genomas relacionados. [R14] [R20]
9. **Hay TSS intragenéricos y antisentido** que el contexto de operón no
   captura. [R16] [R26]
10. **La condición es poblacional**: un sigma activo no garantiza que ese TSS se
    use en todas las células. [R24]

### 11.3 Límites de método

11. **La predicción de operones tiene sensibilidad limitada** (30–50 % en
    *E. coli*) y límites inciertos. [R41] [R42] [R43]
12. **La anotación funcional tiene errores** (5–63 % en superfamilias
    enzimáticas). [R50] [R52]
13. **La filogenia no son réplicas independientes**; los tests de conservación
    pueden ser anti-conservadores. [R63]
14. **La co-expresión está confundida** por crecimiento, condición y estructura
    de operón. [R32] [R26]
15. **Los LLR son inestables con pocos positivos**; σ raras tienen IC enormes.
16. **La calibración no es transferible** entre especies/condiciones sin
    recalibración. [R61] [R63]

### 11.4 Vacíos de la literatura

- No hay un benchmark sistemático de features de contexto para asignación de
  sigma entre especies.
- No hay un marco consensuado de priors calibrados para sigma.
- No hay una forma estándar de *provenance holdout* en RegulonDB/EcoCyc.
- Los regulones predichos se evalúan a menudo con las mismas bases que los
  inspiraron; faltan evaluaciones con *holdout* por paper/estudio.
- Las ECF y las sigmas raras tienen muy pocos datos experimentales.
- No hay una métrica estándar de «anotación independiente» para genomas nuevos.

---

## 12. Recomendaciones

### 12.1 P0 — antes de usar contexto

1. **Tabla de procedencia** por feature: fuente, versión, fecha, paper, ¿cita el
   gold?, ¿se deriva de la etiqueta?
2. **Baseline fuerte** (PWM/secuencia) y **mismo subconjunto** de evaluación.
3. **Cobertura y missingness** por feature; `missing → LLR = 0`.
4. ***Nested CV***: elegir `λ`/escala en *inner loop*; nunca en test. [R69]
5. **Reportar** top-1/top-3, log loss, Brier, ECE, cobertura, IC bootstrap y
   McNemar. [R62] [R65] [R67]
6. **Controles negativos** y **permutación de etiquetas**. [R66] [R70]

### 12.2 P1 — para que la conclusión sea defendible

7. ***Provenance holdout***: quitar features de la fuente del gold y re-evaluar.
8. **Cross-species** (K-12 ↔ *B. subtilis*) y **temporal** (versión antigua →
   nueva). [R08]
9. **Ablación condicional/grupal**, no solo *add-one*.
10. **Calibración** Platt/temperature + corrección de prevalencia. [R61] [R63]
11. **Reportar context-only y oracle** (módulo verdadero) como cotas.
12. **Abstención** por cobertura/conflicto; no forzar una etiqueta.

### 12.3 P2 — para generalizar

13. **Modelo jerárquico multi-especie** con *partial pooling*; efectos de sigma
    y especie.
14. **Integrar datos independientes**: TSS-seq, ChIP-seq, reporteros,
    mutagénesis.
15. **Validación experimental** de un subconjunto: es la única forma de pasar de
    correlación a causalidad.
16. **Documentar el lenguaje**: «prioriza», «es consistente con», no «prueba».

### 12.4 Plantilla de reporte

```text
Feature | Fuente | Versión | Granularidad | Cobertura | Independencia
        | Dirección esperada | Δ top-1 [IC] | Δ log loss | Δ ECE
        | Δ provenance holdout | Δ cross-species | Veredicto
```

Ejemplo de veredicto:

```text
operon_position | BSGatlas/DBTBS | 2024 | TSS/operón | 100 %
  | Baja (comparte gold) | σ70/σ38 | +8,18 pp [..]
  | +0,02 | −0,01 | +0,00 | −0,05 | No usar para evaluar; solo anotar
```

### 12.5 Guía de lenguaje

| En vez de… | Decir… |
|---|---|
| «El contexto prueba que es σ28» | «El contexto eleva el prior de σ28» |
| «El gen es σ54» | «El gen está anotado en el módulo de nitrógeno; su TSS es candidato σ54» |
| «La base dice σ32» | «La base lo cura como σ32 con evidencia X» |
| «El modelo predice σ38» | «El modelo prioriza σ38; la probabilidad calibrada es p» |
| «El regulón de σ54 es el de nitrógeno» | «σ54 controla nitrógeno y otros módulos» |
| «La accuracy es 90 %» | «Top-1 90 % en esta población; cobertura X; IC Y» |

### 12.6 Protocolo mínimo para este proyecto

1. Congelar una versión de gold y de features; registrar fecha.
2. Marcar cada feature con nivel de independencia (A/B/C/D):
   - A: experimental independiente (ChIP, perturbación, TSS-seq).
   - B: computacional independiente (conservación, co-expresión de otro
     compendio).
   - C: anotación derivada (operón, función, ruta).
   - D: derivada del gold (sigmulón, regulón de la misma base).
3. Evaluar A+B por separado; C solo como prior; D prohibido en evaluación.
4. Correr *provenance holdout* y cross-species.
5. Si el contexto no supera el baseline en A+B, no presentarlo como mejora.
6. Reportar el resultado negativo de K-12 como evidencia de que el contexto no
   es universal.
7. No afirmar causalidad sin validación experimental.

### 12.7 Frase final

**El contexto genómico es una herramienta de priorización, no un prueba de
sigma.** Ayuda cuando es independiente, cubre, es proximal y sobrevive a la
ablación; sesga cuando comparte procedencia con el gold, codifica densidad de
anotación o se ajusta en test. En los datos del proyecto, K-12 no muestra señal
contextual out-of-fold y *B. subtilis* sí, pero con procedencia compartida: la
conclusión honesta es que el contexto **puede priorizar** y que **todavía no
está probado** que capture el mecanismo.

---

## 13. Referencias verificadas

> **Nota:** todas las entradas fueron resueltas el 2026-11-26. La columna
> «Verificación» indica el método. Los años de cita son los del volumen de la
> revista (PubMed `pubdate`), no necesariamente el *online first* de Crossref.
> Ver Anexo A para las discrepancias.

### A. Bases de datos y semántica

**R01** — Salgado H, Gama-Castro S, Lara P, *et al.* (2024). RegulonDB v12.0:
a comprehensive resource of transcriptional regulation in *E. coli* K-12.
*Nucleic Acids Research* 52(D1):D255–D264. PMID: 37971353.
DOI: [10.1093/nar/gkad1072](https://doi.org/10.1093/nar/gkad1072).
Verificación: Crossref + PubMed + PMC10767902.

**R02** — Keseler IM, Gama-Castro S, Mackie A, *et al.* (2021). The EcoCyc
Database in 2021. *Frontiers in Microbiology* 12:711077. PMID: 34394059.
DOI: [10.3389/fmicb.2021.711077](https://doi.org/10.3389/fmicb.2021.711077).
Verificación: Crossref + PubMed + PMC8357350.

**R03** — Pedreira T, Elfmann C, Stülke J (2022). The current state of
SubtiWiki, the database for the model organism *Bacillus subtilis*.
*Nucleic Acids Research* 50(D1):D875–D882. PMID: 34664671.
DOI: [10.1093/nar/gkab943](https://doi.org/10.1093/nar/gkab943).
Verificación: Crossref + PubMed + PMC8728116.

**R04** — Dudek C-A, Jahn D (2022). PRODORIC: state-of-the-art database of
prokaryotic gene regulation. *Nucleic Acids Research* 50(D1):D295–D302.
PMID: 34850133.
DOI: [10.1093/nar/gkab1110](https://doi.org/10.1093/nar/gkab1110).
Verificación: Crossref + PubMed + PMC8728284.

**R05** — Sierro N, Makita Y, de Hoon M, *et al.* (2008). DBTBS: a database of
transcriptional regulation in *Bacillus subtilis* containing upstream intergenic
conservation information. *Nucleic Acids Research* 36(Database issue):D93–D96.
PMID: 17962296.
DOI: [10.1093/nar/gkm910](https://doi.org/10.1093/nar/gkm910).
Verificación: Crossref + PubMed + PMC2247474.

**R06** — Novichkov PS, Kazakov AE, Ravcheev DA, *et al.* (2013). RegPrecise
3.0 — a resource for genome-scale exploration of transcriptional regulation in
bacteria. *BMC Genomics* 14:745. PMID: 24175918.
DOI: [10.1186/1471-2164-14-745](https://doi.org/10.1186/1471-2164-14-745).
Verificación: Crossref + PubMed + PMC3840689.

**R07** — Mejía-Almonte C, Busby SJW, Wade JT, *et al.* (2020). Redefining
fundamental concepts of transcription initiation in bacteria. *Nature Reviews
Genetics* 21(11):699–714. PMID: 32665585.
DOI: [10.1038/s41576-020-0254-8](https://doi.org/10.1038/s41576-020-0254-8).
Verificación: Crossref + PubMed.

**R08** — Gama-Castro S, Salgado H, Santos-Zavaleta A, *et al.* (2016).
RegulonDB version 9.0: high-level integration of gene regulation, coexpression,
motif clustering and beyond. *Nucleic Acids Research* 44(D1):D133–D143.
PMID: 26527724.
DOI: [10.1093/nar/gkv1156](https://doi.org/10.1093/nar/gkv1156).
Verificación: Crossref + PubMed.

### B. Biología de sigma y regulones

**R09** — Gruber TM, Gross CA (2003). Multiple sigma subunits and the
partitioning of bacterial transcription space. *Annual Review of Microbiology*
57:441–466. PMID: 14527287.
DOI: [10.1146/annurev.micro.57.030502.090913](https://doi.org/10.1146/annurev.micro.57.030502.090913).

**R10** — Feklistov A, Sharon BD, Darst SA, *et al.* (2014). Bacterial sigma
factors: a historical, structural, and genomic perspective. *Annual Review of
Microbiology* 68:357–376. PMID: 25002089.
DOI: [10.1146/annurev-micro-092412-155737](https://doi.org/10.1146/annurev-micro-092412-155737).

**R11** — Paget MS (2015). Bacterial sigma factors and anti-sigma factors:
structure, function and distribution. *Biomolecules* 5(3):1245–1265.
PMID: 26131973.
DOI: [10.3390/biom5031245](https://doi.org/10.3390/biom5031245).

**R12** — Helmann JD (2002). The extracytoplasmic function (ECF) sigma factors.
*Advances in Microbial Physiology* 46:47–110. PMID: 12073657.
DOI: [10.1016/S0065-2911(02)46002-X](https://doi.org/10.1016/S0065-2911(02)46002-X).

**R13** — Staroń A, Sofia HJ, Dietrich S, *et al.* (2009). The third pillar of
bacterial signal transduction: classification of the extracytoplasmic function
(ECF) sigma factor protein family. *Molecular Microbiology* 74(3):557–581.
PMID: 19737356.
DOI: [10.1111/j.1365-2958.2009.06870.x](https://doi.org/10.1111/j.1365-2958.2009.06870.x).

**R14** — Rhodius VA, Suh WC, Nonaka G, *et al.* (2006). Conserved and variable
functions of the σE stress response in related genomes. *PLoS Biology*
4(1):e2. PMID: 16336047.
DOI: [10.1371/journal.pbio.0040002](https://doi.org/10.1371/journal.pbio.0040002).

**R15** — Chilcott GS, Hughes KT (2000). Coupling of flagellar gene expression
to flagellar assembly in *Salmonella enterica* serovar Typhimurium and
*Escherichia coli*. *Microbiology and Molecular Biology Reviews* 64(4):694–708.
PMID: 11104815.
DOI: [10.1128/MMBR.64.4.694-708.2000](https://doi.org/10.1128/MMBR.64.4.694-708.2000).

**R16** — Fitzgerald DM, Bonocora RP, Wade JT (2014). Comprehensive mapping of
the *Escherichia coli* flagellar regulatory network. *PLoS Genetics*
10(10):e1004649. PMID: 25275371.
DOI: [10.1371/journal.pgen.1004649](https://doi.org/10.1371/journal.pgen.1004649).

**R17** — Reitzer L (2003). Nitrogen assimilation and global regulation in
*Escherichia coli*. *Annual Review of Microbiology* 57:155–176. PMID: 12730324.
DOI: [10.1146/annurev.micro.57.030502.090820](https://doi.org/10.1146/annurev.micro.57.030502.090820).

**R18** — Bush M, Dixon R (2012). The role of bacterial enhancer binding
proteins as specialized activators of σ54-dependent transcription.
*Microbiology and Molecular Biology Reviews* 76(3):497–529. PMID: 22933558.
DOI: [10.1128/MMBR.00006-12](https://doi.org/10.1128/MMBR.00006-12).

**R19** — Barrios H, Valderrama B, Morett E (1999). Compilation and analysis of
σ54-dependent promoter sequences. *Nucleic Acids Research* 27(22):4305–4313.
PMID: 10536136.
DOI: [10.1093/nar/27.22.4305](https://doi.org/10.1093/nar/27.22.4305).

**R20** — Francke C, Groot Kormelink T, Hagemeijer Y, *et al.* (2011).
Comparative analyses imply that the enigmatic sigma factor 54 is a central
controller of the bacterial exterior. *BMC Genomics* 12:385. PMID: 21806785.
DOI: [10.1186/1471-2164-12-385](https://doi.org/10.1186/1471-2164-12-385).

**R21** — Nonaka G, Blankschien M, Herman C, *et al.* (2006). Regulon and
promoter analysis of the *E. coli* heat-shock factor, σ32, reveals a multifaceted
cellular response to heat stress. *Genes & Development* 20(13):1776–1789.
PMID: 16818608.
DOI: [10.1101/gad.1428206](https://doi.org/10.1101/gad.1428206).

**R22** — Guisbert E, Yura T, Rhodius VA, *et al.* (2008). Convergence of
molecular, modeling, and systems approaches for an understanding of the
*Escherichia coli* heat shock response. *Microbiology and Molecular Biology
Reviews* 72(3):545–554. PMID: 18772288.
DOI: [10.1128/MMBR.00007-08](https://doi.org/10.1128/MMBR.00007-08).

**R23** — Weber H, Polen T, Heuveling J, *et al.* (2005). Genome-wide analysis
of the general stress response network in *Escherichia coli*: σS-dependent
genes, promoters, and sigma factor selectivity. *Journal of Bacteriology*
187(5):1591–1603. PMID: 15716429.
DOI: [10.1128/JB.187.5.1591-1603.2005](https://doi.org/10.1128/JB.187.5.1591-1603.2005).

**R24** — Hengge-Aronis R (2002). Signal transduction and regulatory mechanisms
involved in control of the σS (RpoS) subunit of RNA polymerase. *Microbiology
and Molecular Biology Reviews* 66(3):373–395. PMID: 12208995.
DOI: [10.1128/MMBR.66.3.373-395.2002](https://doi.org/10.1128/MMBR.66.3.373-395.2002).

**R25** — Price CW, Fawcett P, Cérémonie H, *et al.* (2001). Genome-wide
analysis of the general stress response in *Bacillus subtilis*. *Molecular
Microbiology* 41(4):757–774. PMID: 11532142.
DOI: [10.1046/j.1365-2958.2001.02534.x](https://doi.org/10.1046/j.1365-2958.2001.02534.x).

**R26** — Nicolas P, Mäder U, Dervyn E, *et al.* (2012). Condition-dependent
transcriptome reveals high-level regulatory architecture in *Bacillus subtilis*.
*Science* 335(6072):1103–1106. PMID: 22383849.
DOI: [10.1126/science.1206848](https://doi.org/10.1126/science.1206848).

**R27** — Stragier P, Losick R (1996). Molecular genetics of sporulation in
*Bacillus subtilis*. *Annual Review of Genetics* 30:297–341. PMID: 8982457.
DOI: [10.1146/annurev.genet.30.1.297](https://doi.org/10.1146/annurev.genet.30.1.297).

**R28** — Errington J (2003). Regulation of endospore formation in *Bacillus
subtilis*. *Nature Reviews Microbiology* 1(2):117–126. PMID: 15035041.
DOI: [10.1038/nrmicro750](https://doi.org/10.1038/nrmicro750).

### C. Métodos de contexto genómico

**R29** — Overbeek R, Fonstein M, D'Souza M, *et al.* (1999). The use of gene
clusters to infer functional coupling. *Proceedings of the National Academy of
Sciences USA* 96(6):2896–2901. PMID: 10077608.
DOI: [10.1073/pnas.96.6.2896](https://doi.org/10.1073/pnas.96.6.2896).

**R30** — Marcotte EM, Pellegrini M, Ng HL, *et al.* (1999). Detecting protein
function and protein-protein interactions from genome sequences. *Science*
285(5428):751–753. PMID: 10427000.
DOI: [10.1126/science.285.5428.751](https://doi.org/10.1126/science.285.5428.751).

**R31** — Pellegrini M, Marcotte EM, Thompson MJ, *et al.* (1999). Assigning
protein functions by comparative genome analysis: protein phylogenetic profiles.
*Proceedings of the National Academy of Sciences USA* 96(8):4285–4288.
PMID: 10200254.
DOI: [10.1073/pnas.96.8.4285](https://doi.org/10.1073/pnas.96.8.4285).

**R32** — Eisen MB, Spellman PT, Brown PO, *et al.* (1998). Cluster analysis
and display of genome-wide expression patterns. *Proceedings of the National
Academy of Sciences USA* 95(25):14863–14868. PMID: 9843981.
DOI: [10.1073/pnas.95.25.14863](https://doi.org/10.1073/pnas.95.25.14863).

**R33** — Kensche PR, van Noort V, Dutilh BE, *et al.* (2008). Practical and
theoretical advances in predicting the function of a protein by its phylogenetic
distribution. *Journal of the Royal Society Interface* 5(19):151–170.
PMID: 17535793.
DOI: [10.1098/rsif.2007.1047](https://doi.org/10.1098/rsif.2007.1047).

**R34** — Lee I, Date SV, Adai AT, *et al.* (2004). A probabilistic functional
network of yeast genes. *Science* 306(5701):1555–1558. PMID: 15567862.
DOI: [10.1126/science.1099511](https://doi.org/10.1126/science.1099511).

**R35** — Gillis J, Pavlidis P (2011). The impact of multifunctional genes on
"guilt by association" analysis. *PLoS ONE* 6(2):e17258. PMID: 21364756.
DOI: [10.1371/journal.pone.0017258](https://doi.org/10.1371/journal.pone.0017258).

**R36** — Gillis J, Pavlidis P (2012). "Guilt by association" is the exception
rather than the rule in gene networks. *PLoS Computational Biology*
8(3):e1002444. PMID: 22479173.
DOI: [10.1371/journal.pcbi.1002444](https://doi.org/10.1371/journal.pcbi.1002444).

**R37** — Pavlidis P, Gillis J (2012). Progress and challenges in the
computational prediction of gene function using networks. *F1000Research*
1:14. PMID: 23936626.
DOI: [10.12688/f1000research.1-14.v1](https://doi.org/10.12688/f1000research.1-14.v1).

**R38** — Huynen M, Snel B, Lathe W, *et al.* (2000). Predicting protein
function by genomic context: quantitative evaluation and qualitative inferences.
*Genome Research* 10(8):1204–1210. PMID: 10958638.
DOI: [10.1101/gr.10.8.1204](https://doi.org/10.1101/gr.10.8.1204).

**R39** — von Mering C, Jensen LJ, Snel B, *et al.* (2005). STRING: known and
predicted protein-protein associations, integrated and transferred across
organisms. *Nucleic Acids Research* 33(Database issue):D433–D437.
PMID: 15608232.
DOI: [10.1093/nar/gki005](https://doi.org/10.1093/nar/gki005).

**R40** — Szklarczyk D, Kirsch R, Koutrouli M, *et al.* (2023). The STRING
database in 2023: protein–protein association networks and functional enrichment
analyses for any sequenced genome of interest. *Nucleic Acids Research*
51(D1):D638–D646. PMID: 36370105.
DOI: [10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000).

**R41** — Ermolaeva MD, White O, Salzberg SL (2001). Prediction of operons in
microbial genomes. *Nucleic Acids Research* 29(5):1216–1221. PMID: 11222772.
DOI: [10.1093/nar/29.5.1216](https://doi.org/10.1093/nar/29.5.1216).

**R42** — Dam P, Olman V, Harris K, *et al.* (2007). Operon prediction using
both genome-specific and general genomic information. *Nucleic Acids Research*
35(1):288–298. PMID: 17170009.
DOI: [10.1093/nar/gkl1018](https://doi.org/10.1093/nar/gkl1018).

**R43** — Brouwer RWW, Kuipers OP, van Hijum SAFT (2008). The relative value of
operon predictions. *Briefings in Bioinformatics* 9(5):367–375. PMID: 18420711.
DOI: [10.1093/bib/bbn019](https://doi.org/10.1093/bib/bbn019).

**R44** — Salgado H, Moreno-Hagelsieb G, Smith TF, *et al.* (2000). Operons in
*Escherichia coli*: genomic analyses and predictions. *Proceedings of the
National Academy of Sciences USA* 97(12):6652–6657. PMID: 10823905.
DOI: [10.1073/pnas.110147297](https://doi.org/10.1073/pnas.110147297).

### D. Predicción de promotores y regulones

**R45** — Bailey TL, Boden M, Buske FA, *et al.* (2009). MEME SUITE: tools for
motif discovery and searching. *Nucleic Acids Research* 37(Web Server
issue):W202–W208. PMID: 19458158.
DOI: [10.1093/nar/gkp335](https://doi.org/10.1093/nar/gkp335).

**R46** — Grant CE, Bailey TL, Noble WS (2011). FIMO: scanning for occurrences
of a given motif. *Bioinformatics* 27(7):1017–1018. PMID: 21330290.
DOI: [10.1093/bioinformatics/btr064](https://doi.org/10.1093/bioinformatics/btr064).

**R47** — Nikolaichik Y, Damienikan AU (2016). SigmoID: a user-friendly tool
for improving bacterial genome annotation through analysis of transcription
control signals. *PeerJ* 4:e2056. PMID: 27257541.
DOI: [10.7717/peerj.2056](https://doi.org/10.7717/peerj.2056).

**R48** — Liu B, Yang F, Huang D-S, *et al.* (2018). iPromoter-2L: a two-layer
predictor for identifying promoters and their types by multi-window-based
PseKNC. *Bioinformatics* 34(1):33–40. PMID: 28968797.
DOI: [10.1093/bioinformatics/btx579](https://doi.org/10.1093/bioinformatics/btx579).

**R49** — Achterberg T, de Jong A (2025). ProPr54 web server: predicting σ54
promoters and regulon with a hybrid convolutional and recurrent deep neural
network. *NAR Genomics and Bioinformatics* 7(1):lqae188. PMID: 39781509.
DOI: [10.1093/nargab/lqae188](https://doi.org/10.1093/nargab/lqae188).

### E. Anotación, sesgo y circularidad

**R50** — Schnoes AM, Brown SD, Dodevski I, *et al.* (2009). Annotation error in
public databases: misannotation of molecular function in enzyme superfamilies.
*PLoS Computational Biology* 5(12):e1000605. PMID: 20011109.
DOI: [10.1371/journal.pcbi.1000605](https://doi.org/10.1371/journal.pcbi.1000605).

**R51** — Salzberg SL (2019). Next-generation genome annotation: we still
struggle to get it right. *Genome Biology* 20:92. PMID: 31097009.
DOI: [10.1186/s13059-019-1715-2](https://doi.org/10.1186/s13059-019-1715-2).

**R52** — Gilks WR, Audit B, De Angelis D, *et al.* (2002). Modeling the
percolation of annotation errors in a database of protein sequences.
*Bioinformatics* 18(12):1641–1649. PMID: 12490449.
DOI: [10.1093/bioinformatics/18.12.1641](https://doi.org/10.1093/bioinformatics/18.12.1641).

**R53** — Kapoor S, Narayanan A (2023). Leakage and the reproducibility crisis
in machine-learning-based science. *Patterns* 4(9):100804. PMID: 37720327.
DOI: [10.1016/j.patter.2023.100804](https://doi.org/10.1016/j.patter.2023.100804).

**R54** — Kaufman S, Rosset S, Perlich C, *et al.* (2012). Leakage in data
mining: formulation, detection, and avoidance. *ACM Transactions on Knowledge
Discovery from Data* 6(4):1–21. DOI:
[10.1145/2382577.2382579](https://doi.org/10.1145/2382577.2382579).
Verificación: Crossref (sin PMID).

**R55** — Whalen S, Schreiber J, Noble WS, *et al.* (2022). Navigating the
pitfalls of applying machine learning in genomics. *Nature Reviews Genetics*
23(3):169–181. PMID: 34837041.
DOI: [10.1038/s41576-021-00434-9](https://doi.org/10.1038/s41576-021-00434-9).

**R56** — Zhou N, Jiang Y, Bergquist T, *et al.* (2019). The CAFA challenge
reports improved protein function prediction and new functional annotations for
hundreds of genes through experimental screens. *Genome Biology* 20:244.
PMID: 31744546.
DOI: [10.1186/s13059-019-1835-8](https://doi.org/10.1186/s13059-019-1835-8).

**R57** — Tompa M, Li N, Bailey TL, *et al.* (2005). Assessing computational
tools for the discovery of transcription factor binding sites. *Nature
Biotechnology* 23(1):137–144. PMID: 15637633.
DOI: [10.1038/nbt1053](https://doi.org/10.1038/nbt1053).

**R58** — Worsley Hunt R, Wasserman WW (2014). Non-targeted transcription
factors motifs are a systemic component of ChIP-seq datasets. *Genome Biology*
15(7):412. PMID: 25070602.
DOI: [10.1186/s13059-014-0412-4](https://doi.org/10.1186/s13059-014-0412-4).

**R59** — Poldrack RA (2006). Can cognitive processes be inferred from
neuroimaging data? *Trends in Cognitive Sciences* 10(2):59–63. PMID: 16406760.
DOI: [10.1016/j.tics.2005.12.004](https://doi.org/10.1016/j.tics.2005.12.004).

### F. Calibración, evaluación y causalidad

**R60** — Karlin S, Altschul SF (1990). Methods for assessing the statistical
significance of molecular sequence features by using general scoring schemes.
*Proceedings of the National Academy of Sciences USA* 87(6):2264–2268.
PMID: 2315319.
DOI: [10.1073/pnas.87.6.2264](https://doi.org/10.1073/pnas.87.6.2264).

**R61** — King G, Zeng L (2001). Logistic regression in rare events data.
*Political Analysis* 9(2):137–163.
DOI: [10.1093/oxfordjournals.pan.a004868](https://doi.org/10.1093/oxfordjournals.pan.a004868).
Verificación: Crossref (sin PMID).

**R62** — Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical Association*
102(477):359–378.
DOI: [10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437).
Verificación: Crossref (sin PMID).

**R63** — Niculescu-Mizil A, Caruana R (2005). Predicting good probabilities with
supervised learning. *Proceedings of the 22nd International Conference on Machine
Learning (ICML '05)*:625–632.
DOI: [10.1145/1102351.1102430](https://doi.org/10.1145/1102351.1102430).
Verificación: Crossref (sin PMID).

**R64** — Domingos P, Pazzani M (1997). On the optimality of the simple Bayesian
classifier under zero-one loss. *Machine Learning* 29(2–3):103–130.
DOI: [10.1023/A:1007413511361](https://doi.org/10.1023/A:1007413511361).
Verificación: Crossref (sin PMID).

**R65** — Dietterich TG (1998). Approximate statistical tests for comparing
supervised classification learning algorithms. *Neural Computation*
10(7):1895–1923. PMID: 9744903.
DOI: [10.1162/089976698300017197](https://doi.org/10.1162/089976698300017197).

**R66** — Breiman L (2001). Random forests. *Machine Learning* 45(1):5–32.
DOI: [10.1023/A:1010933404324](https://doi.org/10.1023/A:1010933404324).
Verificación: Crossref (sin PMID).

**R67** — Efron B (1979). Bootstrap methods: another look at the jackknife.
*The Annals of Statistics* 7(1):1–26.
DOI: [10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552).
Verificación: Crossref (sin PMID).

**R68** — Felsenstein J (1985). Phylogenies and the comparative method. *The
American Naturalist* 125(1):1–15.
DOI: [10.1086/284325](https://doi.org/10.1086/284325).
Verificación: Crossref (sin PMID).

**R69** — Cawley GC, Talbot NLC (2010). On over-fitting in model selection and
subsequent selection bias in performance evaluation. *Journal of Machine
Learning Research* 11:2079–2107.
URL: [https://jmlr.org/papers/v11/cawley10a.html](https://jmlr.org/papers/v11/cawley10a.html).
Verificación: URL resuelta; título comparado.

**R70** — Fisher A, Rudin C, Dominici F (2019). All models are wrong, but many
are useful: learning a variable's importance by studying an entire class of
prediction models simultaneously. *Journal of Machine Learning Research*
20(177):1–81.
URL: [https://jmlr.org/papers/v20/18-760.html](https://jmlr.org/papers/v20/18-760.html).
Verificación: URL resuelta; título comparado.

**R71** — Guo C, Pleiss G, Sun Y, *et al.* (2017). On calibration of modern
neural networks. *Proceedings of the 34th International Conference on Machine
Learning (PMLR 70)*:1321–1330.
URL: [https://proceedings.mlr.press/v70/guo17a.html](https://proceedings.mlr.press/v70/guo17a.html).
Verificación: URL resuelta; título comparado.

**R72** — Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method
for evaluating prediction models. *Medical Decision Making* 26(6):565–574.
PMID: 17099194.
DOI: [10.1177/0272989X06295361](https://doi.org/10.1177/0272989X06295361).

**R73** — Pearl J (2009). *Causality: Models, Reasoning, and Inference*
(2nd ed.). Cambridge University Press.
DOI: [10.1017/CBO9780511803161](https://doi.org/10.1017/CBO9780511803161).
Verificación: Crossref (libro, sin PMID).

---

## Anexo A. Registro de verificación

### A.1 Método

- **Crossref REST** `https://api.crossref.org/works/{DOI}` para cada DOI:
  título, autores, revista, volumen, páginas y año.
- **NCBI E-utilities** `esearch`/`esummary`/`efetch` para PMIDs y abstracts;
  `pubdate` para el año de cita.
- **PMC** para el texto completo de las seis bases (R01–R06), con los PMCIDs:
  PMC10767902, PMC8357350, PMC8728116, PMC8728284, PMC2247474, PMC3840689.
- **JMLR/PMLR** para R69–R71: página resuelta y `<title>` comparado.
- Fecha de verificación: **2026-11-26**.

### A.2 Correcciones y erratas detectadas

| Caso | Lo que se encontró | Resolución |
|---|---|---|
| DOI candidato `10.1371/journal.pone.0033066` | Correspondía a «Composition and Evolution of the Vertebrate and Mammalian Selenoproteomes», no a GBA | Se reemplazó por `10.1371/journal.pone.0017258` (Gillis & Pavlidis 2011) |
| DOI candidato `10.1093/nar/gks1103` | Correspondía a UUCD, no a RegPrecise | Se reemplazó por `10.1186/1471-2164-14-745` |
| Año de R14 (Rhodius) | Crossref `issued` = 2005; PubMed/volumen = 2006 | Se cita 2006 (PLoS Biol 4:e2) |
| Año de R42 (Dam) | Crossref `issued` = 2006; PubMed/volumen = 2007 | Se cita 2007 (NAR 35:288–298) |
| Año de R40 (STRING) | Crossref `issued` = 2022; PubMed/volumen = 2023 | Se cita 2023 (NAR 51:D638–D646) |
| Año de R55 (Whalen) | Crossref `issued` = 2021; PubMed/volumen = 2022 | Se cita 2022 (Nat Rev Genet 23:169–181) |
| R01, R03, R04 | Crossref `issued` un año antes del volumen | Se citan 2024, 2022, 2022 |
| R05 (DBTBS) | Online 2007; volumen 2008 | Se cita 2008 (NAR 36:D93–D96) |
| R49 (ProPr54) | Abstract de PubMed mal formado/truncado | Solo se usa la afirmación del título: predice promotores y regulón σ54 con CNN/RNN |
| R02 (EcoCyc) | Abstract de PubMed truncado | Las afirmaciones se limitan al texto completo de PMC |
| R12 (Helmann) | Crossref no devuelve autores del capítulo | Se usa PMID 12073657 y el título; autor verificado en PubMed |
| R54, R61–R64, R66–R68 | Sin PMID (ACM, JASA, ICML, Machine Learning, Annals, Am Nat) | DOI verificado en Crossref; se declara «sin PMID» |

### A.3 Referencias no usadas o descartadas

- Se descartó un DOI de «RegPrecise» que resolvía a otra base (UUCD).
- Se descartó un DOI de «GBA» que resolvía a selenoproteomas.
- No se incluyeron preprints ni resúmenes de congreso.
- No se incluyeron referencias no verificables por DOI/PMID/URL.

---

## Anexo B. Estrategia de búsqueda y límites

### B.1 Consultas principales

- PubMed: `sigma factor regulon`, `sigma54 regulon`, `ECF sigma factor`,
  `flagellar sigma 28`, `RpoS regulon`, `RpoH regulon`, `guilt by association
  gene function`, `genomic context operon prediction`, `phylogenetic profiling
  function`, `operon prediction`, `annotation bias`, `data leakage genomics`,
  `calibration log odds`, `proper scoring rules`.
- Crossref: búsquedas bibliográficas por título para bases y métodos.
- PMC: texto completo de RegulonDB, EcoCyc, SubtiWiki, PRODORIC, DBTBS y
  RegPrecise para verificar semántica de evidencia y sigmulón.

### B.2 Criterios

- Inclusión: artículos primarios, revisiones y bases de datos revisadas por pares;
  inglés; sin límite de fecha; prioridad a trabajos canónicos y a revisiones
  recientes.
- Exclusión: preprints no revisados, resúmenes de congreso, referencias sin
  DOI/PMID/URL verificable, afirmaciones no comprobables.
- No es PRISMA: la selección es narrativa y está sesgada hacia *E. coli*,
  *B. subtilis* y métodos con código/datos públicos.

### B.3 Límites

- Solo inglés.
- No se hizo búsqueda sistemática en Embase/Scopus.
- La literatura de contexto genómico es enorme; se priorizaron los trabajos
  fundacionales y las revisiones que evalúan rendimiento.
- La ausencia de una referencia no implica ausencia de evidencia.
- Las bases cambian; las cifras de versión son una foto de 2026-11-26.

---

## Anexo C. Matriz de afirmaciones clave

| # | Afirmación | Cita | Tipo |
|---|---|---|---|
| 1 | RegulonDB tiene un datamart de sigmulón | [R01] | Base de datos |
| 2 | RegulonDB y EcoCyc comparten equipo de curación | [R01] | Base de datos |
| 3 | EcoCyc modela cada TSS con sigma como promotor distinto | [R02] [R07] | Base de datos |
| 4 | SubtiWiki asigna genes a regulones | [R03] | Base de datos |
| 5 | PRODORIC distingue 3 tipos de evidencia | [R04] | Base de datos |
| 6 | DBTBS tiene 1.475 promotores y 463 operones validados | [R05] | Base de datos |
| 7 | RegPrecise reconstruye regulones por genómica comparativa | [R06] | Base de datos |
| 8 | Un TSS puede tener varios promotores/sigmas | [R02] [R07] | Concepto |
| 9 | σ28 controla genes flagelares tardíos; FlgM lo inhibe | [R15] | Estudio primario |
| 10 | FliA tiene 52 sitios, 30 intragenéricos | [R16] | Estudio primario |
| 11 | σ54 usa −24/−12; 186 promotores compilados | [R19] | Estudio primario |
| 12 | σ54 requiere bEBP/ATP | [R18] | Revisión |
| 13 | σ54 controla nitrógeno, motilidad, virulencia, biofilm, envelope | [R20] | Estudio comparativo |
| 14 | σ32 regulón; σ70 inicia 12 % de sus promotores | [R21] | Estudio primario |
| 15 | σS controla hasta 10 % de genes; solapamiento con otros regulones | [R23] | Estudio primario |
| 16 | σB gobierna la respuesta general de estrés en *B. subtilis* | [R25] | Estudio primario |
| 17 | 2.935 promotores agrupados en regulones de sigma | [R26] | Estudio primario |
| 18 | ECF: >40 grupos mayores; anti-sigma | [R12] [R13] | Revisión |
| 19 | σE: 89 unidades de transcripción en genomas relacionados | [R14] | Estudio primario |
| 20 | σE de *B. subtilis* es de esporulación, no ECF | [R27] [R28] | Revisión |
| 21 | La multifuncionalidad confunde GBA | [R35] | Estudio metodológico |
| 22 | GBA es la excepción, no la regla | [R36] | Estudio metodológico |
| 23 | El contexto genómico cubre ~37–50 % de genes | [R38] | Benchmark |
| 24 | La predicción de operones tiene sensibilidad 30–50 % | [R41] | Estudio primario |
| 25 | La accuracy de operones depende de información previa | [R42] | Estudio primario |
| 26 | La misanotación en bases públicas es 5–63 % | [R50] | Estudio primario |
| 27 | Los errores de anotación percolan | [R52] | Modelado |
| 28 | El leakage tiene 8 tipos y afecta 294 papers | [R53] | Revisión |
| 29 | ChIP-seq puede tener motivos no objetivo | [R58] | Estudio primario |
| 30 | Naive Bayes puede ser óptimo para clasificar pero no para calibrar | [R64] | Teoría |
| 31 | La corrección de prevalencia es necesaria en eventos raros | [R61] | Método |
| 32 | Los proper scoring rules son la referencia | [R62] | Método |
| 33 | La calibración se ajusta en datos separados | [R63] [R71] | Método |
| 34 | La importancia por permutación requiere cuidado | [R66] [R70] | Método |
| 35 | El *nested CV* evita el sesgo de selección | [R69] | Método |
| 36 | La filogenia introduce no independencia | [R68] | Método |
| 37 | La inferencia inversa no prueba la causa | [R59] | Concepto |
| 38 | El contexto puede predecir sin explicar | [R73] | Concepto |

---

## Anexo D. Glosario mínimo

- **Sigma (σ):** subunidad de RNAP que confiere especificidad de promotor.
- **Sigmulón:** conjunto de genes/operones transcritos por RNAP con un sigma.
- **Regulón:** conjunto de genes regulados por un TF; no es un sigmulón.
- **Operón:** genes co-transcritos desde un mismo promotor.
- **TSS:** sitio de inicio de la transcripción (1-based en este proyecto).
- **Promotor:** segmento de DNA que recluta un holoenzima RNAP–sigma; puede haber
  varios por TSS.
- **PWM/PSSM:** matriz de frecuencias/pesos de un motivo.
- **LLR:** log-likelihood ratio; en este documento, a menudo una asociación
  estimada, no una likelihood.
- **Prior:** distribución previa explícita; no es una frecuencia del gold.
- **Posterior:** probabilidad calibrada bajo el modelo y el prior.
- **Calibración:** acuerdo entre probabilidad predicha y frecuencia observada.
- **ECE:** expected calibration error; diagnóstico, no proper score.
- **Brier / log loss:** proper scoring rules para probabilidades.
- **GBA:** *guilt by association*; predecir una propiedad por asociación con
  vecinos/partners.
- **Perfil filogenético:** vector de presencia/ausencia de ortólogos.
- **Sintenia:** conservación del orden de genes.
- **Procedencia:** fuente, versión, paper y cadena de anotación de una feature.
- **Leakage:** fuga de información del test al entrenamiento o a la feature.
- **Circularidad:** la etiqueta y la feature comparten un ancestro causal o
  procedural.
- **Ablación:** quitar/añadir features para medir su contribución.
- **Nested CV:** CV externa para evaluar, interna para elegir hiperparámetros.
- **PU learning:** aprendizaje con positivos y no etiquetados (el gold no tiene
  negativos).
- **ECF:** *extracytoplasmic function*; grupo 4 de la familia σ70.
- **bEBP:** *bacterial enhancer binding protein*; activador AAA+ de σ54.
- **Anti-sigma:** proteína que secuestra un sigma.

---

*Fin del documento. No se realizó commit ni se modificó código. La única
escritura es este archivo de revisión.*

