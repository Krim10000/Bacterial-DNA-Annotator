# Factores sigma bacterianos: familias, arquitectura de promotores, ortología y distribución

> **Documento de revisión bibliográfica para `dna_poc`.**
> Ruta: `docs/literature/SIGMA_BIOLOGY.md`
> Fecha de corte y verificación: **2026-11-26**
> Autor: subagente `lit-sigma-biology` (revisión bibliográfica; **no modifica código**)
> Estado: **v1.0** — todas las referencias de la tabla final fueron resueltas una a una
> contra Europe PMC, NCBI E-utilities, Crossref, OpenAlex, InterPro o el PDF del
> editor. No se incluye ninguna referencia no verificada. Ver Anexo A.

---

## Convenciones

- `[PMID:<número>]`: artículo indexado en PubMed; el PMID fue resuelto y el título
  comparado con la fuente.
- `[DOI:<registro>]`: DOI resuelto; título, autores y año comparados.
- `[InterPro:PFxxxxx]`: entrada de InterPro/Pfam resuelta vía API; el nombre y la
  definición se citan literalmente en el Anexo B.
- Las **observaciones internas** del proyecto se marcan explícitamente como tales y
  no se presentan como literatura.

Este documento no es una revisión sistemática PRISMA. Es una revisión narrativa con
trazabilidad de fuentes, orientada a decisiones de diseño. La estrategia de búsqueda
y sus límites están en el Anexo C.

---

## 0. Resumen ejecutivo

1. **Hay dos familias estructuralmente no relacionadas**: la familia σ70 (que incluye
   a las sigma alternativas y a las ECF) y la familia σ54, sin homología de secuencia
   detectable con σ70. [PMID:26131973] [PMID:12073657] [PMID:11763973]
2. **La familia σ70 se subdivide en cuatro grupos** (1: housekeeping; 2: primarias no
   esenciales; 3: alternativas; 4: ECF), aunque el número de grupos y el estatus de
   las ECF varían entre autores. [PMID:26131973] [PMID:12540296] [PMID:1597408]
3. **La arquitectura de dominios explica la arquitectura del promotor**: σ2 reconoce
   la caja –10 y une el core; σ3 reconoce la –10 extendida y une el core; σ4 reconoce
   la caja –35; σ1.1 actúa como compuerta. [PMID:26131973] [PMID:11931761]
   [PMID:22136875]
4. **Cada familia tiene una geometría de promotor distinta**: σ70 usa –35/–10 (con
   –10 extendida, UP y discriminador); σ54 usa –24/–12 y requiere activadores
   dependientes de ATP; las ECF usan motivos –35/–10 específicos de grupo.
   [PMID:18521075] [PMID:10536136] [PMID:22933558] [PMID:33398323]
5. **El número de sigma por genoma varía de 1 a >100**, y depende de la biología del
   linaje y de la anotación. *E. coli* K-12 tiene 7; *B. subtilis* 168 tiene 17–18;
   *M. tuberculosis* H37Rv 13; *P. aeruginosa* PAO1 24; *S. coelicolor* A3(2) 65–66;
   *M. genitalium* 1. [PMID:38709457] [PMID:14993308] [PMID:17064287] [PMID:18070067]
   [PMID:30371884] [PMID:11763973]
6. **El nombre nativo no es un grupo ortólogo**: RpoD/SigA/HrdB son ortólogos
   housekeeping; FliA/SigD son ortólogos flagelares; RpoN/SigL son ortólogos σ54;
   pero σE, σF, σH y σB designan funciones distintas en *E. coli* y *B. subtilis*.
   [PMID:30371884] [PMID:11763973] [PMID:11121773] [PMID:16336047]
7. **La transferencia automática de nombres es un problema documentado**, no una
   hipótesis: la nomenclatura de sigma se describió como confusa y necesitada de
   revisión; la anotación funcional automática tiene tasas de error de 5–63 % en
   superfamilias enzimáticas de bases públicas; y los modelos HMM de dominio tienen
   falsos negativos declarados por sus propios curadores.
   [PMID:11763973] [PMID:20011109] [InterPro:PF04545]
8. **Implicación de diseño**: no debe existir un panel fijo de seis sigma. Debe
   existir un **inventario por especie** con nombre nativo, familia, subgrupo,
   arquitectura de dominios, ortólogo, contexto genómico, evidencia y confianza; y un
   **modelo de promotor por familia/subgrupo**, no un escaneo σ70 por defecto.
   Ver §11.

---

## 1. Las dos familias: σ70 y σ54

### 1.1 Qué las define

La clasificación de primer nivel no se basa en el nombre del gen ni en el producto
anotado, sino en la **homología de secuencia y la arquitectura de dominios**.
Paget (2015) lo formula explícitamente: las sigma se clasifican en dos familias
distintas según su homología con dos sigma de *E. coli*: σ70, responsable de la mayor
parte de la transcripción durante el crecimiento, y σ54 (σN), estructuralmente no
relacionada, que dirige la transcripción en respuesta a señales ambientales y requiere
proteínas activadoras y ATP para abrir el DNA. [PMID:26131973]

Mittenhuber (2002), en un censo de 31 genomas completamente secuenciados, es igual de
explícito: “No sequence similarity can be detected between the σ70- and the σ54-
family”. [PMID:11763973] La partición del espacio de transcripción entre sigma
—cuántas sigma compiten por el core y cómo se reparte el programa— fue revisada por
Gruber & Gross (2003). [PMID:14527287]

### 1.2 Tabla comparativa de familias

| Familia | Relación con σ70 | Tamaño típico | Dominios | Geometría de promotor | Activador | Ejemplos |
|---|---|---|---:|---|---|---|
| **σ70** (grupos 1–3) | — | ~40–90 kDa | σ1.1 (solo grupo 1), σ2, σ3, σ4, linkers | –35 / –10; –10 extendida; UP; discriminador | Generalmente no requiere activador para unión; sí para iniciación en algunos promotores | RpoD/σ70, RpoS/σ38, RpoH/σ32, FliA/σ28 |
| **σ70 grupo 4 (ECF)** | Sí, dentro de σ70 | ~20 kDa | σ2 y σ4 (a veces σ3); sin σ1.1 | –35 / –10 específicos de grupo; espaciador variable | Frecuentemente anti-sigma y activación por liberación | RpoE/σE, FecI, SigW, SigX, AlgU |
| **σ54** | No; familia separada | ~50–60 kDa | AID, DBD, CBD (plegamiento distinto) | –24 / –12; sin –35/–10 canónicas | Sí: bEBP + ATP | RpoN/σ54, SigL |

**Citas de la tabla:** σ70/σ54 como familias separadas [PMID:26131973]
[PMID:11763973]; tamaños y dominios [PMID:26131973] [PMID:11931761]; ECF como grupo 4
[PMID:26131973] [PMID:12073657]; σ54 con –24/–12 y bEBP/ATP [PMID:10536136]
[PMID:22933558]; estructura de σ54 [PMID:26293966].

### 1.3 Lo que no hay que concluir

- **“σ54 es una sigma alternativa más”**: no lo es en el mismo sentido. Su mecanismo
  de apertura del promotor es dependiente de activador y ATP, y su estructura no es
  homóloga a σ70. [PMID:22933558] [PMID:26293966]
- **“ECF es una familia aparte”**: la literatura la trata mayoritariamente como un
  subgrupo/grupo 4 dentro de σ70, no como una tercera familia homóloga.
  [PMID:26131973] [PMID:12073657] [PMID:19737356]
- **“El nombre del producto basta para asignar familia”**: falso; ver §7 y §8.

---

## 2. La familia σ70: grupos, dominios y arquitectura

### 2.1 Los cuatro grupos

La subdivisión más citada es la de cuatro grupos. [PMID:26131973] [PMID:12540296]
Lonetto et al. (1992) establecieron las bases de la comparación de secuencia y las
regiones conservadas. [PMID:1597408] La perspectiva histórica, estructural y genómica
fue actualizada por Feklistov et al. (2014) [PMID:25002089], y la comparación clásica
de estructura y función se remonta a Helmann & Chamberlin (1988) [PMID:3052291].
Wösten (1998) revisó 14 clases de sigma, lo que muestra que la granularidad de la
clasificación ha cambiado con el tiempo. [PMID:9818380]

| Grupo | Nombre funcional | ¿Esencial? | Dominios | Promotor | Ejemplos | Citas |
|---|---|---|---|---|---|---|
| **1** | Primaria / housekeeping | Sí | σ1.1, σ2, σ3, σ4 | –35/–10 canónicas; UP; discriminador | RpoD (σ70) *E. coli*; SigA *B. subtilis*; HrdB *S. coelicolor* | [PMID:26131973] [PMID:12540296] [PMID:30371884] |
| **2** | Primaria no esencial / estrés | No | σ2, σ3, σ4; sin σ1.1 | Similar a grupo 1; –13/–14 C en σS | RpoS (σ38/σS) *E. coli*; SigB (σB) *B. subtilis* | [PMID:26131973] [PMID:11251833] [PMID:21639793] |
| **3** | Alternativas | No | σ2, σ3, σ4 (variable) | –10/–35 divergentes; –10 extendida en algunos | RpoH (σ32), FliA (σ28), σE/σK/σF/σG de esporulación | [PMID:26131973] [PMID:19400790] [PMID:11121773] |
| **4** | ECF | No | σ2 y σ4; sin σ1.1; a veces sin σ3 | –35/–10 específicos de grupo | RpoE (σE) *E. coli*; FecI; SigW/SigX/SigM *B. subtilis* | [PMID:26131973] [PMID:12073657] [PMID:16336047] |

**Advertencia de clasificación.** El límite entre grupo 2 y grupo 3, y entre grupo 3
y grupo 4, no es universal. Distintas revisiones usan distintos criterios (presencia
de σ1.1, esencialidad, arquitectura de dominios, filogenia). Esto se discute en §9.

### 2.2 Arquitectura de dominios de σ70

Las sigma de la familia σ70 tienen hasta cuatro regiones conservadas que pliegan en
dominios separados, conectados por linkers flexibles:

| Región/dominio | Contacto principal | Función | Citas |
|---|---|---|---|
| **σ1.1** (σR1.1) | Canal activo de RNAP; DNA | Compuerta; impide unión no productiva de σ libre al promotor; se desplaza al formarse el complejo abierto | [PMID:26131973] |
| **σ1.2–2.4 (σ2)** | –10 (cadena no molde) y β′ coiled-coil | Reconocimiento de –10, melting, unión al core; es la región más conservada | [PMID:26131973] [PMID:22136875] |
| **σ3.0–3.2 (σ3)** | DNA justo antes de –10; core | Reconocimiento de –10 extendida; unión al core; el linker σ3.2 forma el “σ finger” | [PMID:26131973] [PMID:11931761] |
| **σ4.1–4.2 (σ4)** | –35 (HTH) y β flap | Reconocimiento de –35; interfaz con activadores | [PMID:26131973] [PMID:11931761] |
| **NCR** (región no conservada) | Variable | Tamaño y composición variable; implicada en unión al core y escape en σ70 | [PMID:26131973] |

La estructura de los dominios σ2, σ3 y σ4 de σA de *Thermus aquaticus*, y de σ4 unido
a DNA de –35, fue resuelta por Campbell et al. (2002); los dominios están conectados
por linkers flexibles y las superficies expuestas son importantes para unir RNAP.
[PMID:11931761] La holoenzima σ70 de *E. coli* fue resuelta por Murakami (2013).
[PMID:23389035] La base estructural del reconocimiento de –10 por σ2 fue resuelta por
Feklistov & Darst (2011). [PMID:22136875] La estructura del complejo de iniciación
con burbuja completa fue resuelta por Zhang et al. (2012). [PMID:23086998]

### 2.3 Arquitectura canónica del promotor σ70

Los elementos mejor caracterizados en *E. coli* son:

| Elemento | Posición | Consenso/característica | Citas |
|---|---|---|---|
| Caja –35 | ~–35 | `TTGACA` (consenso clásico) | [PMID:26131973] [PMID:18521075] |
| Caja –10 | ~–10 | `TATAAT`; A–11 y T–7 se voltean y se entierran en bolsillos de σ2 | [PMID:26131973] [PMID:22136875] |
| Espaciador –35/–10 | entre ambas | Longitud óptima ~17 pb; su longitud afecta la iniciación | [PMID:19139410] [PMID:18521075] |
| –10 extendida | –15/–14 | `TG`; reconocida por σ3; puede compensar la ausencia de –35 | [PMID:26131973] [PMID:11931761] |
| UP element | aguas arriba de –35 | Rico en A/T; unido por el dominio αCTD de RNAP | [PMID:18521075] [DOI:10.1038/nrmicro.2016.103] |
| Discriminador | –6/+1 | `GGG` óptimo; interacciona con σ1.2; afecta estabilidad del complejo abierto | [PMID:26131973] |

**Matiz importante:** estos elementos son de σ70 de *E. coli* y de grupos 1–2. No deben
extrapolarse sin más a grupos 3, ECF o σ54. La revisión de Haugen et al. (2008)
describe la plasticidad de la arquitectura del promotor y cómo factores que no unen
DNA (ppGpp, DksA, factores de transcripción) modulan la iniciación.
[PMID:18521075] Browning & Busby (2016) actualizan el panorama de regulación local y
global. [DOI:10.1038/nrmicro.2016.103]

---

## 3. La familia σ54

### 3.1 Propiedades

σ54 (σN; gen `rpoN`) es una familia separada. Sus promotores tienen posiciones
conservadas alrededor de –24 y –12 en lugar de las cajas –35 y –10 típicas.
[PMID:10536136] La RNAP–σ54 forma un complejo cerrado incompetente para la
iniciación; la transición a complejo abierto requiere un activador (bEBP,
*Bacterial Enhancer Binding Protein*) y hidrólisis de ATP. [PMID:22933558]
La estructura de RNAP–σ54 mostró que σ54 es realmente distinta de σ70 y explicó por
qué es incompetente para iniciar por sí sola. [PMID:26293966]

| Propiedad | σ54 | Citas |
|---|---|---|
| Promotor | –24 (GG) / –12 (GC); sin –35/–10 canónicas | [PMID:10536136] |
| Activador | bEBP con dominio AAA+; ATP | [PMID:22933558] |
| Estructura | AID, DBD, CBD; plegamiento distinto de σ70 | [PMID:26293966] |
| Regulación | Típicamente nitrogeno, flagelos, pili, metabolismo | [PMID:22933558] |
| Homología con σ70 | No detectable | [PMID:11763973] [PMID:26131973] |

### 3.2 Consecuencia para el escaneo de promotores

Un escáner que busque –35/–10 no encontrará promotores σ54; y un escáner que busque
–24/–12 con geometría σ70 generará falsos positivos. El módulo `sigma` del proyecto ya
declara que σ54 se escanea con geometría –24/–12, no –35/–10 (observación interna:
`modules/sigma/module.py`). La literatura respalda esa separación.
[PMID:10536136] [PMID:22933558]

**Vacío:** la compilación de Barrios et al. (1999) contiene 186 promotores –24/–12;
sigue siendo una de las compilaciones experimentales más grandes, pero está sesgada a
organismos modelo y a promotores bien estudiados. [PMID:10536136]

---

## 4. La familia ECF (grupo 4 de σ70)

### 4.1 Qué son

Las ECF (*extracytoplasmic function*) son el grupo más numeroso de sigma alternativas
y el tercer mecanismo fundamental de transducción de señales bacterianas, después de
los sistemas de uno y dos componentes. [PMID:19737356] Helmann (2002) las definió como
un subgrupo de σ70 que responde a estímulos extracitoplasmáticos. [PMID:12073657]

Son proteínas pequeñas (~20 kDa) con σ2 y σ4 (a veces σ3), sin σ1.1. Su gen suele
estar co-transcrito con un gen de anti-sigma adyacente, y la liberación de la sigma
ocurre por proteólisis regulada, *partner switching* o sensado directo por la
anti-sigma. [PMID:26131973] [PMID:19737356]

### 4.2 Clasificación y censo

| Métrica | Staroń et al. 2009 | Casas-Pastor et al. 2021 | Interpretación |
|---|---:|---:|---|
| ECF por genoma (media) | ~6 | ~10 | El censo creció con la base de datos y la reclasificación |
| Grupos filogenéticos ECF | >40 | 157 | Expansión ~50× respecto de clasificaciones previas |
| Grupos nuevos | — | 22 | Nuevos grupos sin proteínas previamente clasificadas |
| Criterio | Filogenia + contexto genómico + anti-sigma + motivos de promotor | Filogenia + vecindad génica + mecanismo de regulación + motivos | La clasificación no depende del nombre del producto |

**Citas:** [PMID:19737356] [PMID:33398323].

### 4.3 Promotores ECF

Los promotores ECF tienen cajas –35 y –10, pero los motivos son **específicos de
grupo** y el espaciador es variable. Casas-Pastor et al. (2021) usan la clasificación
en grupos y subgrupos para predecir motivos comunes de promotor e inferir mecanismos
regulatorios. [PMID:33398323] Staroń et al. (2009) ya habían señalado que la
clasificación se apoya en motivos diana potenciales. [PMID:19737356]

Para σE de *E. coli*, Rhodius et al. (2006) desarrollaron y validaron un modelo de
predicción de promotores, y predijeron 89 unidades de transcripción únicas en *E. coli*
K-12 y ocho genomas relacionados; el artículo destaca que la variabilidad de secuencias
y el espaciado de motivos dificultan la predicción. [PMID:16336047]

**Consecuencia:** no existe un “promotor ECF” único. Un modelo PWM por grupo ECF es
más defendible que un consenso global. [PMID:33398323] [PMID:19737356]

### 4.4 Expansión y ortología

Las expansiones de ECF son reales y linaje-específicas. En el censo de 2002,
*P. aeruginosa* tenía 19 ECF, *M. loti* 18 y *C. crescentus* 14; el texto señala que
muchas de estas expansiones parecen haber evolucionado por duplicación génica.
[PMID:11763973] La filogenia de ECF y la transferencia horizontal complican la
asignación 1:1 de ortólogos. [PMID:33575573] [PMID:33398323]

---

## 5. Sigma alternativas: función, nombre y arquitectura de promotor

### 5.1 Tabla de referencia rápida

La siguiente tabla resume las sigma alternativas más citadas. **No es un panel**: es
una lista de casos bien estudiados, con sus nombres nativos y sus colisiones.

| Sigma | Nombre(s) nativo(s) | Familia/grupo | Función | Arquitectura de promotor | Organismo modelo | Citas |
|---|---|---|---|---|---|---|
| **σ70** | RpoD; σD; SigA; HrdB | σ70 grupo 1 | Housekeeping | –35/–10; UP; discriminador | *E. coli*, *B. subtilis*, *S. coelicolor* | [PMID:38709457] [PMID:30371884] |
| **σ38/σS** | RpoS; σS | σ70 grupo 2 | Estrés general, fase estacionaria | Similar a σ70; C en –13/–14 | *E. coli* | [PMID:26131973] [PMID:11251833] [PMID:21639793] |
| **σ32/σH** | RpoH; σH | σ70 grupo 3 | Heat shock | –35/–10 divergentes; no canónicas | *E. coli* | [PMID:38709457] [PMID:26131973] |
| **σ28/σF** | FliA; σF; σ28 | σ70 grupo 3 | Flagelos, quimiotaxis | –10 larga `GCCGATAA` (ext –10 + core –10) | *E. coli* | [PMID:19400790] [PMID:38709457] |
| **σE (ECF)** | RpoE; σE; AlgU en *Pseudomonas* | σ70 grupo 4 | Estrés de envoltura | –35/–10 ECF; modelo específico | *E. coli*, *P. aeruginosa* | [PMID:16336047] [PMID:18070067] |
| **σFecI** | FecI; σFecI | σ70 grupo 4 | Transporte de citrato férrico | –35/–10 ECF | *E. coli* | [PMID:38709457] [PMID:11763973] |
| **σ54/σN** | RpoN; σN; NtrA; SigL | σ54 | Nitrógeno, flagelos, pili | –24/–12; bEBP + ATP | *E. coli*, *B. subtilis* | [PMID:10536136] [PMID:22933558] [PMID:11763973] |
| **σB** | SigB; σB (histórico σ37) | σ70 grupo 2/3 | Estrés general | Promotores tipo σB | *B. subtilis* | [PMID:7708009] [PMID:11763973] |
| **σE (esporulación)** | SigE; σE | σ70 grupo 3 | Esporulación, célula madre | Promotores de esporulación | *B. subtilis* | [PMID:11121773] [PMID:7708009] |
| **σF (esporulación)** | SigF; σF | σ70 grupo 3 | Esporulación, forespora | Promotores de esporulación | *B. subtilis* | [PMID:11121773] [PMID:7708009] |
| **σG** | SigG; σG | σ70 grupo 3 | Esporulación, forespora tardía | Promotores de esporulación | *B. subtilis* | [PMID:11121773] [PMID:7708009] |
| **σK** | SigK; σK | σ70 grupo 3 | Esporulación, célula madre tardía | Promotores de esporulación | *B. subtilis* | [PMID:11121773] [PMID:7708009] |
| **σH (B. subtilis)** | SigH; σH | σ70 grupo 3 | Transición/competencia | Promotores de transición | *B. subtilis* | [PMID:11121773] |
| **σI** | SigI; σI (YkoZ) | σ70 grupo 3/4 (discutido) | Heat shock; ECF-like | Promotores SigI | *B. subtilis* | [PMID:11157964] [PMID:17185538] |
| **σW/σX/σM/σV/σY/σZ/YlaC** | SigW, SigX, SigM, SigV, SigY, SigZ, YlaC | σ70 grupo 4 (ECF) | Estrés de pared, sales, oxidativo, etc. | –35/–10 ECF específicos | *B. subtilis* | [PMID:14993308] [PMID:11763973] |

### 5.2 La cascada de esporulación de *B. subtilis*

La esporulación de *B. subtilis* es el ejemplo canónico de sigma alternativas
organizadas en cascada. Haldenwang (1995) describió 10 sigma identificadas en ese
momento y los controles postraduccionales: σE y σK se sintetizan como precursores
inactivos que se activan por proteólisis; σB, σF y σG son inhibidas por anti-sigma
factores que las secuestran. [PMID:7708009] Kroos & Yu (2000) revisaron la regulación
de la actividad de sigma durante el desarrollo, incluida la cascada σH → σF → σE → σG
→ σK. [PMID:11121773]

**Lección de diseño:** la regulación por anti-sigma y por proteólisis significa que la
presencia del gen no implica actividad; y que el contexto genómico (operón
sigma/anti-sigma) es evidencia funcional. [PMID:26131973] [PMID:11121773]

### 5.3 Colisiones de nombre

| Nombre | Significado en *E. coli* | Significado en *B. subtilis* | ¿Son ortólogos? | Citas |
|---|---|---|---|---|
| **σE** | RpoE, ECF, estrés de envoltura | SigE, esporulación (célula madre) | **No** | [PMID:16336047] [PMID:11121773] |
| **σF** | FliA, flagelar (σ28) | SigF, esporulación (forespora) | **No** | [PMID:19400790] [PMID:11121773] |
| **σH** | RpoH, heat shock (σ32) | SigH, transición/competencia | **No** | [PMID:38709457] [PMID:11121773] |
| **σB** | No existe (RpoS es el análogo funcional) | SigB, estrés general | **No 1:1** | [PMID:11763973] [PMID:7708009] |
| **σD** | RpoD (housekeeping) | SigD, flagelar (ortólogo de FliA) | **No** | [PMID:38709457] [PMID:11763973] |
| **σA** | No se usa (RpoD) | SigA, housekeeping | **Sí, con RpoD** | [PMID:30371884] [PMID:11763973] |

**Cita clave:** Mittenhuber (2002) afirma que “la nomenclatura actual para los genes
que codifican sigma es confusa y debería revisarse”, y documenta que el ortólogo de
FliA en *B. subtilis* se llama σD. [PMID:11763973]

---

## 6. ¿Cuántos sigma tiene cada especie?

### 6.1 Conteos modernos para especies modelo

| Especie | Conteo reportado | Composición (si se reporta) | Fuente | Caveat |
|---|---:|---|---|---|
| *E. coli* K-12 | **7** | 4 σ70 + 1 σ54 + 2 ECF (RpoD, RpoS, RpoH, FliA; RpoN; RpoE, FecI) | [PMID:38709457] [PMID:18070067] [PMID:11763973] | Consistente entre fuentes |
| *B. subtilis* 168 | **17–18** | 18 en el censo de 2002 (10 σ70 + 1 σ54 + 7 ECF, incluye SigI y YlaC); “al menos 17” en Yoshimura 2004 | [PMID:14993308] [PMID:11763973] | “Al menos 17”; el 18 depende de la inclusión |
| *M. tuberculosis* H37Rv | **13** | 3 σ70 + 0 σ54 + 10 ECF | [PMID:17064287] [PMID:11763973] | 13 “putativos” |
| *P. aeruginosa* PAO1 | **24** | 4 σ70 + 1 σ54 + 19 ECF | [PMID:18070067] [PMID:11763973] | 19 ECF es expansión real |
| *S. coelicolor* A3(2) | **65–66** | HrdB + 65 alternativas; 66 total | [PMID:30371884] | El mismo artículo dice “66” y “otros 65” |
| *M. genitalium* | **1** | 1 σ70 (housekeeping) | [PMID:30371884] [PMID:11763973] [PMID:7569993] | Genoma mínimo, 470 genes |
| *S. aureus* | **3** | 3 σ70; 0 σ54; 0 ECF | [PMID:11763973] | Censo 2002; puede haber cambiado |
| *C. crescentus* | **17** | 2 σ70 + 1 σ54 + 14 ECF | [PMID:11763973] | Expansión ECF |
| *Synechocystis* 6803 | **8** | 6 σ70 + 0 σ54 + 2 ECF | [PMID:11763973] | Censo 2002; revisar con anotación actual |
| *D. radiodurans* | **2** | 1 σ70-like no ortóloga + 1 ECF; 0 σ54 | [PMID:11763973] | Excepción a “todos tienen housekeeping” |
| *Sorangium cellulosum* | **109** | No reportada en la fuente | [PMID:26131973] | Ejemplo extremo de expansión |

### 6.2 Censo histórico de 31 genomas (2002)

Mittenhuber (2002) comparó 215 sigma de 31 genomas completamente secuenciados y
reportó, entre otros: *B. halodurans* 19 (8 σ70 + 1 σ54 + 10 ECF), *C. acetobutylicum*
14 (10 + 1 + 3), *M. loti* 23 (3 + 2 + 18), *P. aeruginosa* 24 (4 + 1 + 19),
*C. crescentus* 17 (2 + 1 + 14). [PMID:11763973]

El mismo censo concluyó que, con la excepción de *D. radiodurans*, todas las especies
poseen un factor sigma housekeeping; que las especies de vida libre poseen más sigma
alternativas que patógenos o simbiontes; y que las especies versátiles tienen más ECF.
[PMID:11763973]

**Genomas de referencia:** *E. coli* K-12 [PMID:9278503], *B. subtilis* 168
[PMID:9384377], *M. tuberculosis* H37Rv [PMID:9634230] y *S. coelicolor* A3(2)
[PMID:12000953].

**Advertencia:** estos números son de 2002, con los genomas y la anotación de entonces.
No deben usarse como verdad actual sin reanotar. La propia literatura posterior muestra
cambios: *S. coelicolor* se reporta como “65” o “66” según cómo se cuente.
[PMID:30371884]

### 6.3 Por qué el número varía

1. **Biología real:** expansiones de ECF en actinobacterias, mixobacterias y
   proteobacterias; pérdida de sigma en genomas pequeños y simbiontes.
   [PMID:11763973] [PMID:24739623] [PMID:33398323]
2. **Anotación:** genes divididos, fusionados, truncados o mal llamados.
   [PMID:20011109] [PMID:27342282]
3. **Criterio de conteo:** ¿se incluyen pseudogenes? ¿se incluyen factores sin
   dominio completo? ¿se cuenta el housekeeping? [PMID:30371884]
4. **Versión de la base de datos:** los grupos ECF pasaron de >40 a 157.
   [PMID:19737356] [PMID:33398323]

Chávez et al. (2020) cuantificaron la representación de familias de sigma en genomas
bacterianos y arqueanos; en los genomas bacterianos analizados, el sigma housekeeping
estuvo presente en todos. [PMID:33077635]

---

## 7. Ortología vs nombres nativos

### 7.1 El principio

Un **nombre nativo** es una etiqueta histórica, no una hipótesis de ortología. La
ortología se define por la historia genealógica (especiación) y se infiere por
filogenia, no por similitud de nombre ni por similitud global de secuencia.
[PMID:31727128] [PMID:34597405] [PMID:36350662]

En sigma esto es especialmente delicado porque:

- hay duplicaciones y pérdidas diferenciales;
- hay transferencia horizontal;
- hay expansiones linaje-específicas de ECF;
- la misma función puede ser realizada por sigma no ortólogas;
- el mismo nombre puede designar funciones distintas.

[PMID:33575573] [PMID:33398323] [PMID:11763973]

### 7.2 Mapa de ortología por clase funcional

| Clase funcional | *E. coli* K-12 | *B. subtilis* 168 | *M. tuberculosis* H37Rv | *S. coelicolor* A3(2) | *P. aeruginosa* PAO1 | Confianza de ortología | Citas |
|---|---|---|---|---|---|---|---|
| Housekeeping (grupo 1) | RpoD (σ70/σD) | SigA (σA) | SigA | HrdB (SCO5820) | RpoD | **Alta** (1:1, esencial) | [PMID:30371884] [PMID:11763973] |
| σ54 | RpoN (σN) | SigL | Ausente | No verificado | RpoN (NtrA) | **Alta** cuando está presente | [PMID:11763973] [PMID:18070067] |
| Flagelar (grupo 3) | FliA (σF/σ28) | SigD (σD) | No verificado | No verificado | RpoF/FliA | **Alta** | [PMID:11763973] [PMID:19400790] |
| Estrés general (grupo 2) | RpoS (σS/σ38) | SigB (σB) | No verificado | No verificado | RpoS | **Baja/analogía funcional**, no 1:1 | [PMID:11763973] [PMID:21639793] |
| Heat shock (grupo 3) | RpoH (σH/σ32) | Sin ortólogo RpoH; SigI es ECF-like | No verificado | No verificado | RpoH | **RpoH restringido a proteobacterias** | [PMID:11763973] [PMID:11157964] |
| ECF de envoltura | RpoE (σE) | SigW/SigX/SigM/SigV/SigY/SigZ/YlaC | 10 ECF | SigE/SigR (ECF) | AlgU (σ22) | **Solo a nivel de grupo ECF**, no 1:1 | [PMID:16336047] [PMID:14993308] [PMID:18070067] [PMID:33398323] |

**Notas:**

- “Alta” significa que la literatura consultada sostiene la equivalencia funcional y
  la presencia del ortólogo; no significa que se haya verificado filogenéticamente en
  este documento.
- Las celdas “no verificado” son deliberadas: no se rellenan por inferencia.
- La ausencia de RpoH en *B. subtilis* no es un vacío de anotación: el censo de 2002
  indica que la ocurrencia de RpoH de *E. coli* está restringida a proteobacterias.
  [PMID:11763973]
- La presencia de σ54 es variable; su ausencia en *M. tuberculosis* y *S. aureus* está
  reportada en el censo de 2002. [PMID:11763973]

### 7.3 Colisiones y alias que el sistema debe manejar

| Alias | Equivalencia | Organismo | Fuente |
|---|---|---|---|
| RpoD / σ70 / σD / SigA / HrdB | Housekeeping | *E. coli* / *B. subtilis* / *S. coelicolor* | [PMID:30371884] [PMID:38709457] |
| RpoS / σS / σ38 | Estrés general | *E. coli* | [PMID:26131973] |
| RpoH / σH / σ32 | Heat shock | *E. coli* | [PMID:38709457] [PMID:26131973] |
| FliA / σF / σ28 | Flagelar | *E. coli* | [PMID:19400790] [PMID:38709457] |
| RpoN / σN / σ54 / NtrA | σ54 | *E. coli*, *P. aeruginosa* | [PMID:18070067] [PMID:10536136] |
| SigL | σ54 | *B. subtilis* | [PMID:11763973] |
| SigD | FliA/σ28 ortólogo | *B. subtilis* | [PMID:11763973] |
| AlgU / RpoE / σ22 | ECF de envoltura | *P. aeruginosa* | [PMID:18070067] |
| σB / σ37 | Estrés general | *B. subtilis* | [PMID:11763973] [PMID:7708009] |

**Regla:** el alias no es identidad. “RpoE” en *E. coli* y “RpoE” en *P. aeruginosa*
pueden ser ortólogos, pero “σE” en *B. subtilis* no lo es. [PMID:16336047]
[PMID:11121773] [PMID:18070067]

---

## 8. Problemas de anotación automática y transferencia de nombres

### 8.1 Lo que hacen los pipelines

Las pipelines de anotación procariota predicen genes y asignan productos por
combinación de HMMs, homología y reglas. PGAP (NCBI Prokaryotic Genome Annotation
Pipeline) es el ejemplo mejor documentado. [PMID:27342282] RefSeq añade curaduría a esas anotaciones automáticas. [PMID:29112715]

El problema no es que la homología sea inútil; es que **la homología global no
resuelve la ortología fina** y que **el nombre del producto se propaga**. En
superfamilias enzimáticas se midió una misanotación media de 5–63 % en GenBank NR y
TrEMBL, con sobrepredicción como categoría dominante. [PMID:20011109]

### 8.2 Modos de fallo específicos de sigma

| Modo de fallo | Evidencia | Consecuencia | Mitigación |
|---|---|---|---|
| Transferir el nombre sin ortología | Nomenclatura “confusa”; σE/σF/σH con significados distintos | Familia/función equivocada | Nombre nativo + orthogroup + familia por dominio |
| Producto genérico “sigma factor” | Nomenclatura confusa [PMID:11763973]; anotación automática [PMID:27342282] | Sin familia asignable | HMM + arquitectura de dominios |
| Hit parcial de un solo dominio | σ70_r4 tiene umbral alto y falsos negativos declarados [InterPro:PF04545] | Fragmento llamado sigma completo | Cobertura de dominios y longitud |
| Modelo HMM incorrecto en el panel | PF04997 es `RNA_pol_Rpb1_1`, no sigma [InterPro:PF04997] | Falsos positivos | Auditar el panel |
| Duplicaciones y parálogos | Expansiones ECF por duplicación [PMID:11763973] | Colapso de copias o conteo inflado | copy_number + orthogroup |
| Genes truncados/pseudogenes | Anotación automática [PMID:27342282] | Dominio incompleto | Flags de completitud |
| Promotores solapados y cross-talk | Competencia y solapamiento funcional [PMID:10982868] [PMID:16892065] | Asignación única incorrecta | Modelo probabilístico + margen |
| Default σ70 | El escaneo clásico asume –35/–10 [PMID:18521075]; el sistema previo usaba `assumed_default` (interno) | Se pierden alternativas | PWM por familia |

**Observación interna (no literatura):** la auditoría del proyecto
(`docs/SIGMA_AUDIT.md`) documenta que el sistema previo asignaba `sigma70` por defecto
cuando encontraba un par –35/–10, y que un producto génico con “rpoE” en el nombre
podía clasificarse como sigmaE aunque no fuera un factor sigma. Estos son ejemplos
internos del mismo modo de fallo, no evidencia bibliográfica.

### 8.3 El caso de los modelos Pfam

Las definiciones de dominio provienen de Pfam 2021. [PMID:33125078] La verificación de
las entradas de InterPro/Pfam usadas por el proyecto arroja un hallazgo relevante:

| Pfam | Nombre verificado | ¿Es un dominio sigma? | Comentario |
|---|---|---|---|
| PF00140 | Sigma-70 factor, region 1.2 | Sí, grupos 1–2 | Marcador de grupos 1–2; no cubre ECF |
| PF04542 | Sigma-70 region 2 | Sí, núcleo σ70 | Reconocimiento de –10 y unión al core |
| PF04539 | Sigma-70 region 3 | Sí, grupos 1–3 | –10 extendida y unión al core |
| PF04545 | Sigma-70, region 4 | Sí, –35 | **Umbral artificialmente alto; muchos falsos negativos** (texto de Pfam) |
| PF04997 | RNA polymerase Rpb1, domain 1 | **No** | Es el dominio clamp de β′; no es sigma |
| PF07638 | ECF sigma factor | Sí, ECF | Específico de ECF |
| PF00309 | Sigma-54 factor, AID | Sí, σ54 | Dominio de interacción con activador |
| PF04552 | Sigma-54, DNA binding domain | Sí, σ54 | Unión a –24/–12 |
| PF04963 | Sigma-54 factor, core binding domain | Sí, σ54 | Unión al core |
| PF00158 | Sigma-54 interaction domain | **No es sigma** | Es el dominio de bEBPs; no confundir |
| PF12791 | Anti-sigma factor N-terminus | Anti-sigma | Contexto, no factor sigma |

**Citas:** [InterPro:PF00140] [InterPro:PF04542] [InterPro:PF04539] [InterPro:PF04545]
[InterPro:PF04997] [InterPro:PF07638] [InterPro:PF00309] [InterPro:PF04552]
[InterPro:PF04963] [InterPro:PF00158] [InterPro:PF12791].

**Consecuencia:** un panel que use PF04997 como modelo sigma, o que use solo
PF00140/PF04542/PF04539/PF04545, no cubre ECF ni σ54 de forma adecuada. Es un problema
de panel, no de biología. Ver §11.

### 8.4 Lo que la literatura no resuelve

No identificamos una tabla universal de “nombre → familia → orthogroup” para sigma. Las
revisiones comparan especies, pero no publican un mapeo curado y completo de alias.
[PMID:11763973] [PMID:33398323] Esto obliga a construir el mapeo localmente y a
declarar su procedencia.

---

## 9. Discrepancias entre fuentes

| Tema | Fuente A | Fuente B | Lectura |
|---|---|---|---|
| Número de grupos σ70 | 4 grupos (1–4) [PMID:26131973] [PMID:12540296] | 3 grupos + ECF como subfamilia separada [PMID:1597408] [PMID:11763973] | No es contradicción de datos: es granularidad y estatus de ECF |
| Clases de sigma | 14 clases [PMID:9818380] | 4 grupos σ70 + σ54 [PMID:26131973] | “Clase” en Wösten no equivale a “grupo” actual |
| σ de *B. subtilis* | 18 [PMID:11763973] | “Al menos 17” [PMID:14993308] | Depende de incluir SigI/YlaC y del criterio de conteo |
| ECF de *B. subtilis* | 7 (SigM, V, W, X, Y, Z, YlaC) [PMID:11763973] [PMID:14993308] | SigI descrita como heat shock de grupo IV/ECF-like [PMID:11157964] [PMID:17185538] | SigI es ECF-like de grupo IV; el límite ECF/grupo 3 es difuso |
| σ de *S. coelicolor* | 66 total [PMID:30371884] | “Otros 65” [PMID:30371884] | Mismo artículo: 66 incluye HrdB; 65 son las alternativas |
| ECF por genoma | ~6 [PMID:19737356] | ~10 [PMID:33398323] | Crecimiento de la base y reclasificación, no necesariamente contradicción biológica |
| Grupos ECF | >40 [PMID:19737356] | 157 [PMID:33398323] | Expansión ~50× de la clasificación |
| Housekeeping en todas las bacterias | “Todas excepto *D. radiodurans*” [PMID:11763973] | Anotaciones modernas podrían discrepar | Afirmación de 2002; no verificar como universal actual |
| RpoH | Restringida a proteobacterias [PMID:11763973] | *B. subtilis* tiene respuesta heat shock por SigI [PMID:11157964] | No es el mismo ortólogo: analogía funcional, no ortología |
| σE | ECF en *E. coli* [PMID:16336047] | Esporulación en *B. subtilis* [PMID:11121773] | Colisión de nombre, no discrepancia biológica |
| Número de sigma y patogenicidad | Especies de vida libre tienen más sigma [PMID:11763973] | Genomas pequeños/simbiontes pierden sistemas regulatorios [PMID:24739623] | Consistente: la reducción génica reduce repertorio |

**Regla para el documento del proyecto:** cuando dos fuentes discrepan, se reporta el
rango y la fuente, no se elige un número sin declarar el criterio.

---

## 10. Vacíos de literatura

Los siguientes vacíos se identificaron en la búsqueda del 2026-11-26. No son
afirmaciones de inexistencia absoluta; son ausencias en la literatura recuperada con
la estrategia del Anexo C.

1. **No identificamos un benchmark público, cuantitativo y multi-especie de exactitud de
   anotación de sigma** (familia, orthogroup y nombre). La literatura compara conteos
   y filogenias, pero no reporta precisión/recall de la asignación de nombres en
   genomas completos. [PMID:11763973] [PMID:33398323]
2. **Los modelos de promotor están sesgados a *E. coli* y *B. subtilis*.** La
   compilación de promotores σ54 más citada es de 1999 y anterior a la genómica
   masiva; el modelo de σE se validó en *E. coli* y genomas relacionados.
   [PMID:10536136] [PMID:16336047]
3. **Los motivos de promotor ECF son mayoritariamente predichos, no validados
   experimentalmente**, y dependen de la clasificación en grupos.
   [PMID:33398323] [PMID:19737356]
4. **El pareo sigma/anti-sigma y la estructura de operones están incompletos** en
   muchos genomas; la evidencia suele ser vecindad génica, no experimentos.
   [PMID:19737356] [PMID:26131973]
5. **No identificamos un mapeo sistemático de la ausencia de σ54 y de la pérdida de
   sigma en linajes reducidos**; los censos son puntuales. [PMID:11763973] [PMID:24739623]
6. **La función de muchas ECF y de muchos parálogos de expansión es desconocida**;
   *S. coelicolor* es el ejemplo extremo. [PMID:30371884] [PMID:33398323]
7. **No hay nomenclatura unificada**; los alias y las colisiones siguen sin
   resolución comunitaria. [PMID:11763973]
8. **Gran parte de la evidencia cuantitativa sobre competición y cross-talk entre
   sigma proviene de *E. coli***; su generalización a otras especies es incierta.
   [PMID:10982868] [PMID:16892065]
9. **Los conteos no son comparables entre versiones de anotación**; no hay una
   práctica estándar para reportar “sigma por genoma” con versión y criterio.
   [PMID:30371884] [PMID:11763973]
10. **La relación entre repertorio de sigma y estilo de vida está descrita
    cualitativamente**, no modelada de forma predictiva. [PMID:11763973]
    [PMID:24739623]

---

## 11. Implicaciones para nuestro diseño

### 11.1 Principio rector

**No existe un panel fijo de seis sigma. Existe un inventario por especie.**

La evidencia para rechazar un panel fijo es directa:

- el número varía de 1 a 109 entre especies [PMID:30371884] [PMID:26131973];
- la composición cambia: hay especies sin σ54, sin RpoH, con 19 ECF o con 1 sola
  sigma [PMID:11763973] [PMID:18070067];
- los nombres colisionan [PMID:11763973] [PMID:11121773] [PMID:16336047];
- las expansiones y pérdidas son linaje-específicas [PMID:33398323]
  [PMID:33575573];
- la anotación automática introduce errores y los conteos no son directamente
  comparables [PMID:20011109] [PMID:27342282].

Por tanto, el sistema debe **descubrir** el repertorio de cada genoma y **declarar**
lo que no puede asignar, en lugar de forzar seis categorías.

### 11.2 Modelo de datos propuesto

| Campo | Tipo | Obligatorio | Propósito |
|---|---|---|---|
| `genome_id` | string | Sí | Identificador de la especie/cepa |
| `protein_id` / `locus_tag` | string | Sí | Trazabilidad al gen |
| `native_name` | string | No | Nombre tal como está anotado; puede ser vacío |
| `aliases` | list[str] | No | RpoD/SigA/σ70, FliA/SigD/σ28, etc. |
| `family` | enum | Sí | `sigma70`, `sigma54`, `ecf`, `unknown` |
| `subgroup` | enum | No | `group1`, `group2`, `group3`, `group4_ecf`, `sigma54`, `unknown` |
| `orthogroup_id` | string | No | Grupo ortólogo (OrthoFinder/eggNOG/OrthoDB) |
| `domain_architecture` | list[str] | Sí | `r1.2`, `r2`, `r3`, `r4`, `sigma54_AID`, `sigma54_DBD`, `sigma54_CBD` |
| `domain_coverage` | float | Sí | Fracción de dominios esperados presentes |
| `length_aa` | int | Sí | Detecta fragmentos |
| `completeness` | enum | Sí | `full`, `partial_N`, `partial_C`, `fragment` |
| `copy_number` | int | Sí | Distingue copia única de expansión |
| `genomic_context` | object | No | Operón, anti-sigma adyacente, dirección |
| `promoter_model_ids` | list[str] | No | PWMs aplicables |
| `evidence` | list[str] | Sí | `hmm`, `phylogeny`, `name`, `context`, `experiment` |
| `confidence` | enum | Sí | `high`, `medium`, `low` |
| `flags` | list[str] | No | `name_conflict`, `partial_domain`, `paralog_expansion`, `no_ortholog` |

**Regla:** `family` se asigna por dominio, no por nombre; `subgroup` se asigna por
arquitectura + filogenia; `orthogroup_id` se asigna por ortología; `native_name` se
conserva siempre.

### 11.3 Pipeline recomendado

1. **Predicción de proteínas** del genoma.
2. **Escaneo HMM** con panel curado (HMMER [PMID:22039361]):
   - σ70: PF00140, PF04542, PF04539, PF04545;
   - ECF: PF07638;
   - σ54: PF00309, PF04552, PF04963;
   - anti-sigma (contexto): PF12791 y familias específicas.
3. **Filtro de arquitectura**: exigir σ2 + σ4 para σ70/ECF; exigir AID + DBD + CBD
   para σ54; marcar cobertura parcial como `partial_domain`.
4. **Ortología**: OrthoFinder o eggNOG-mapper; asignar `orthogroup_id`; comparar con
   un set de referencia curado. [PMID:31727128] [PMID:34597405] [PMID:36350662]
5. **Reconciliación de nombres**: mantener `native_name`; poblar `aliases` solo con
   evidencia; nunca sobrescribir el nombre nativo.
6. **Modelo de promotor**: asignar PWM por familia/subgrupo; σ54 con geometría –24/–12;
   ECF por grupo; no usar σ70 por defecto. [PMID:10536136] [PMID:33398323]
7. **Integración de evidencia**: HMM + arquitectura + orthogroup + contexto +
   promotor; calcular confianza.
8. **Salida por especie**: lista de factores, ausencias declaradas, expansiones
   marcadas, fragmentos marcados.
9. **QC**: comparar el conteo con el rango de literatura de la especie (§11.5);
   marcar outliers para revisión.
10. **Evaluación**: conjuntos curados (RegulonDB para *E. coli*, SubtiWiki para
    *B. subtilis*), métricas por familia y por orthogroup.

### 11.4 Verificación del panel HMM

El panel actual del proyecto incluye PF00140, PF04542, PF04539, PF04545 y PF04997
(observación interna: `scripts/regulation/12_promoters_sigma.py` y
`data/regulation/*/promoters_meta.json`). La verificación de InterPro muestra:

- **PF04997 no es un dominio sigma**: es `RNA polymerase Rpb1, domain 1`
  (`RNA_pol_Rpb1_1`), el dominio clamp de β′. [InterPro:PF04997]
- **PF04545 tiene umbral alto y falsos negativos** declarados por Pfam.
  [InterPro:PF04545]
- **PF00140 (región 1.2) es un marcador de los grupos 1–2 y no cubre ECF**, así que
  un panel basado solo en él subrepresenta ECF. [InterPro:PF00140]
- **σ54 no está cubierta** por los cinco modelos anteriores; requiere PF00309,
  PF04552 y PF04963. [InterPro:PF00309] [InterPro:PF04552] [InterPro:PF04963]
- **PF00158 es el dominio de interacción de bEBPs**, no de sigma; usarlo como marcador
  sigma produciría falsos positivos. [InterPro:PF00158]

**Recomendación:** retirar PF04997 y añadir PF07638 (ECF) y los tres modelos de σ54;
mantener PF04545 con corrección por longitud/cobertura; documentar la versión de Pfam.

### 11.5 Expectativas por especie (para QC, no para forzar)

| Especie | Rango esperado | Fuente | Uso |
|---|---:|---|---|
| *E. coli* K-12 | 7 | [PMID:38709457] [PMID:18070067] | Ancla |
| *B. subtilis* 168 | 17–18 | [PMID:14993308] [PMID:11763973] | Rango |
| *M. tuberculosis* H37Rv | 13 | [PMID:17064287] | Ancla |
| *P. aeruginosa* PAO1 | 24 | [PMID:18070067] | Ancla |
| *S. coelicolor* A3(2) | 65–66 | [PMID:30371884] | Rango |
| *M. genitalium* | 1 | [PMID:30371884] [PMID:11763973] | Ancla |
| *S. aureus* | 3 | [PMID:11763973] | Solo referencia histórica; no usar para QC actual |
| *C. crescentus* | 17 | [PMID:11763973] | Solo referencia histórica; no usar para QC actual |
| *Synechocystis* 6803 | 8 | [PMID:11763973] | Solo referencia histórica; no usar para QC actual |
| *D. radiodurans* | 2 | [PMID:11763973] | Solo referencia histórica; no usar para QC actual |
| *Sorangium cellulosum* | 109 | [PMID:26131973] | Caso extremo |

**Cómo usarlo:** si el inventario de una especie cae fuera del rango, no se corrige
automáticamente; se marca `outlier` y se revisa la causa (anotación, panel, biología).
No usar estos números como filtros duros.

### 11.6 Reglas de decisión

1. **El nombre no es familia.** Un producto “sigma factor” sin dominio no es una
   sigma asignable. [PMID:11763973] [PMID:20011109]
2. **El nombre no es orthogroup.** RpoD/SigA/HrdB son ortólogos; σE de *E. coli* y
   σE de *B. subtilis* no. [PMID:30371884] [PMID:11121773]
3. **No colapsar parálogos.** Las expansiones ECF son reales; `copy_number` debe ser
   explícito. [PMID:33398323] [PMID:11763973]
4. **No forzar un promotor.** Sin PWM de la familia, la asignación de promotor es
   `unknown`, no `sigma70`. [PMID:18521075] [PMID:33398323]
5. **Declarar la ausencia.** “No se encontró σ54” es un resultado; no un fallo.
   [PMID:11763973]
6. **Separar familia, subgrupo, orthogroup y función.** Son cuatro capas distintas.
7. **Reportar evidencia y confianza.** Toda llamada debe poder rastrearse a HMM,
   arquitectura, contexto o experimento.
8. **No extrapolar promotores de *E. coli* a otras especies** sin validación.
   [PMID:16336047] [PMID:33398323]
9. **Usar `unknown` como salida válida.** El sistema debe preferir la abstención a la
   falsa precisión. [PMID:20011109]
10. **Versionar el panel y la anotación.** Los conteos dependen de ambos.
    [PMID:30371884] [PMID:11763973]

### 11.7 Evaluación propuesta

| Nivel | Qué se evalúa | Conjunto | Métrica |
|---|---|---|---|
| Familia | σ70 vs σ54 vs ECF vs unknown | Curado manual + InterPro | Precisión/recall por familia |
| Subgrupo | Grupo 1/2/3/4 | Curado manual | Matriz de confusión |
| Orthogroup | Consistencia con OrthoFinder/eggNOG | Genomas de referencia | Fracción de asignaciones estables |
| Copia | copy_number vs literatura | Especies de §11.5 | Error absoluto |
| Promotor | Asignación de sigma a promotor | RegulonDB/SubtiWiki/literatura | Precisión/recall; top-k |
| Abstención | `unknown` correcto | Casos ambiguos | Tasa de abstención calibrada |

**Baseline obligatorio:** un clasificador por nombre/producto. Si el sistema basado en
dominios + ortología no supera al baseline por nombre, el diseño no está aportando.

### 11.8 Lo que el sistema no debe afirmar

- “Especie X tiene N sigma” sin versión de anotación y criterio de conteo.
- “Este promotor es de σX” solo por secuencia, sin PWM y sin margen.
- “SigB es el ortólogo de RpoS” (es una analogía funcional, no una ortología 1:1).
- “No hay sigma” cuando el panel no cubre esa familia.
- “El número de ECF predice virulencia” (no hay modelo predictivo validado).

---

## 12. Referencias verificadas

> Todas las referencias listadas fueron resueltas contra al menos una fuente de
> metadatos (Europe PMC, NCBI E-utilities, Crossref, OpenAlex) y, cuando fue posible,
> se leyó el abstract o el texto completo. La columna “Verificación” indica la fuente
> consultada. No se incluyen referencias no verificadas.

| ID | Referencia | PMID | DOI/URL | Verificación |
|---|---|---|---|---|
| R01 | Paget MS. Bacterial sigma factors and anti-sigma factors: structure, function and distribution. *Biomolecules*. 2015;5(3):1245–1265. | 26131973 | 10.3390/biom5031245 | Europe PMC + PubMed abstract |
| R02 | Feklistov A, Sharon BD, Darst SA, Gross CA. Bacterial sigma factors: a historical, structural, and genomic perspective. *Annu Rev Microbiol*. 2014;68:357–376. | 25002089 | 10.1146/annurev-micro-092412-155737 | Europe PMC |
| R03 | Gruber TM, Gross CA. Multiple sigma subunits and the partitioning of bacterial transcription space. *Annu Rev Microbiol*. 2003;57:441–466. | 14527287 | 10.1146/annurev.micro.57.030502.090913 | Europe PMC |
| R04 | Lonetto M, Gribskov M, Gross CA. The sigma 70 family: sequence conservation and evolutionary relationships. *J Bacteriol*. 1992;174(12):3843–3849. | 1597408 | 10.1128/jb.174.12.3843-3849.1992 | Europe PMC |
| R05 | Helmann JD, Chamberlin MJ. Structure and function of bacterial sigma factors. *Annu Rev Biochem*. 1988;57:839–872. | 3052291 | 10.1146/annurev.bi.57.070188.004203 | PubMed esummary |
| R06 | Wösten MM. Eubacterial sigma-factors. *FEMS Microbiol Rev*. 1998;22(3):127–150. | 9818380 | 10.1111/j.1574-6976.1998.tb00364.x | Europe PMC |
| R07 | Helmann JD. The extracytoplasmic function (ECF) sigma factors. *Adv Microb Physiol*. 2002;46:47–110. | 12073657 | 10.1016/S0065-2911(02)46002-X | Europe PMC |
| R08 | Mittenhuber G. An inventory of genes encoding RNA polymerase sigma factors in 31 completely sequenced eubacterial genomes. *J Mol Microbiol Biotechnol*. 2002;4(1):77–91. | 11763973 | http://www.caister.com/backlist/jmmb/v/v4/07.pdf | PubMed esummary + PDF completo |
| R09 | Staroń A, Sofia HJ, Dietrich S, Ulrich LE, Liesegang H, Mascher T. The third pillar of bacterial signal transduction: classification of the extracytoplasmic function (ECF) sigma factor protein family. *Mol Microbiol*. 2009;74(3):557–581. | 19737356 | 10.1111/j.1365-2958.2009.06870.x | Europe PMC + abstract |
| R10 | Casas-Pastor D, Müller RR, Jaenicke S, Brinkrolf K, Becker A, Buttner MJ, Gross CA, Mascher T, Goesmann A, Fritz G. Expansion and re-classification of the extracytoplasmic function (ECF) σ factor family. *Nucleic Acids Res*. 2021;49(2):986–1005. | 33398323 | 10.1093/nar/gkaa1229 | Europe PMC + abstract |
| R11 | Pinto D, da Fonseca RR. Evolution of the extracytoplasmic function σ factor protein family. *NAR Genom Bioinform*. 2020;2(1):lqz026. | 33575573 | 10.1093/nargab/lqz026 | Europe PMC |
| R12 | Paget MS, Helmann JD. The sigma70 family of sigma factors. *Genome Biol*. 2003;4(1):203. | 12540296 | 10.1186/gb-2003-4-1-203 | PubMed abstract |
| R13 | Browning DF, Busby SJW. Local and global regulation of transcription initiation in bacteria. *Nat Rev Microbiol*. 2016;14(10):638–650. | — | 10.1038/nrmicro.2016.103 | Crossref |
| R14 | Haugen SP, Ross W, Gourse RL. Advances in bacterial promoter recognition and its control by factors that do not bind DNA. *Nat Rev Microbiol*. 2008;6(7):507–519. | 18521075 | 10.1038/nrmicro1912 | PubMed esummary |
| R15 | Feklistov A, Darst SA. Structural basis for promoter –10 element recognition by the bacterial RNA polymerase σ subunit. *Cell*. 2011;147(6):1257–1269. | 22136875 | 10.1016/j.cell.2011.10.041 | PubMed esummary |
| R16 | Zhang Y, Feng Y, Chatterjee S, Tuske S, Ho MX, Arnold E, Ebright RH. Structural basis of transcription initiation. *Science*. 2012;338(6110):1076–1080. | 23086998 | 10.1126/science.1227786 | PubMed esummary |
| R17 | Murakami KS. X-ray crystal structure of Escherichia coli RNA polymerase σ70 holoenzyme. *J Biol Chem*. 2013;288(13):9126–9134. | 23389035 | 10.1074/jbc.M112.430900 | PubMed esummary |
| R18 | Campbell EA, Muzzin O, Chlenov M, Sun JL, Olson CA, Weinman O, Trester-Zedlitz ML, Darst SA. Structure of the bacterial RNA polymerase promoter specificity sigma subunit. *Mol Cell*. 2002;9(3):527–539. | 11931761 | 10.1016/S1097-2765(02)00470-7 | PubMed abstract |
| R19 | Hook-Barnard IG, Hinton DM. The promoter spacer influences transcription initiation via σ70 region 1.1 of Escherichia coli RNA polymerase. *Proc Natl Acad Sci USA*. 2009;106(3):737–742. | 19139410 | 10.1073/pnas.0808133106 | PubMed esummary |
| R20 | Barrios H, Valderrama B, Morett E. Compilation and analysis of sigma(54)-dependent promoter sequences. *Nucleic Acids Res*. 1999;27(22):4305–4313. | 10536136 | 10.1093/nar/27.22.4305 | PubMed abstract |
| R21 | Bush M, Dixon R. The role of bacterial enhancer binding proteins as specialized activators of σ54-dependent transcription. *Microbiol Mol Biol Rev*. 2012;76(3):497–529. | 22933558 | 10.1128/MMBR.00006-12 | PubMed esummary |
| R22 | Yang Y, et al. Structures of the RNA polymerase–σ54 reveal new and conserved regulatory strategies. *Science*. 2015;349(6250):882–885. | 26293966 | 10.1126/science.aab1478 | OpenAlex |
| R23 | Koo BM, Rhodius VA, Campbell EA, Gross CA. Mutational analysis of Escherichia coli σ28 and its target promoters reveals recognition of a composite –10 region, comprised of an “extended –10” motif and a core –10 element. *Mol Microbiol*. 2009;72(4):830–843. | 19400790 | 10.1111/j.1365-2958.2009.06691.x | PubMed abstract |
| R24 | Rhodius VA, Suh WC, Nonaka G, West J, Gross CA. Conserved and variable functions of the σE stress response in related genomes. *PLoS Biol*. 2006;4(1):e2. | 16336047 | 10.1371/journal.pbio.0040002 | PubMed abstract |
| R25 | Becker G, Hengge-Aronis R. What makes an Escherichia coli promoter σS dependent? Role of the –13/–14 nucleotide promoter positions and region 2.5 of σS. *Mol Microbiol*. 2001;39(5):1153–1165. | 11251833 | 10.1111/j.1365-2958.2001.02313.x | PubMed esummary |
| R26 | Maeda H, Fujita N, Ishihama A. Competition among seven Escherichia coli σ subunits: relative binding affinities to the core RNA polymerase. *Nucleic Acids Res*. 2000;28(18):3497–3503. | 10982868 | 10.1093/nar/28.18.3497 | OpenAlex |
| R27 | Wade JT, Castro Roa D, Grainger DC, Hurd D, Busby SJ, Struhl K, Nudler E. Extensive functional overlap between σ factors in Escherichia coli. *Nat Struct Mol Biol*. 2006;13(9):806–814. | 16892065 | 10.1038/nsmb1130 | PubMed esummary |
| R28 | Blattner FR, et al. The complete genome sequence of Escherichia coli K-12. *Science*. 1997;277(5331):1453–1462. | 9278503 | 10.1126/science.277.5331.1453 | PubMed abstract |
| R29 | Kunst F, et al. The complete genome sequence of the gram-positive bacterium Bacillus subtilis. *Nature*. 1997;390(6657):249–256. | 9384377 | 10.1038/36786 | PubMed abstract |
| R30 | Bentley SD, et al. Complete genome sequence of the model actinomycete Streptomyces coelicolor A3(2). *Nature*. 2002;417(6885):141–147. | 12000953 | 10.1038/417141a | PubMed abstract |
| R31 | Cole ST, et al. Deciphering the biology of Mycobacterium tuberculosis from the complete genome sequence. *Nature*. 1998;393(6685):537–544. | 9634230 | 10.1038/31159 | PubMed abstract |
| R32 | Rodrigue S, Provvedi R, Jacques PE, Gaudreau L, Manganelli R. The sigma factors of Mycobacterium tuberculosis. *FEMS Microbiol Rev*. 2006;30(6):926–941. | 17064287 | 10.1111/j.1574-6976.2006.00040.x | PubMed abstract |
| R33 | Yoshimura M, Asai K, Sadaie Y, Yoshikawa H. Interaction of Bacillus subtilis extracytoplasmic function (ECF) sigma factors with the N-terminal regions of their potential anti-sigma factors. *Microbiology*. 2004;150(Pt 3):591–599. | 14993308 | 10.1099/mic.0.26712-0 | PubMed abstract |
| R34 | Potvin E, Sanschagrin F, Levesque RC. Sigma factors in Pseudomonas aeruginosa. *FEMS Microbiol Rev*. 2008;32(1):38–55. | 18070067 | 10.1111/j.1574-6976.2007.00092.x | PubMed abstract |
| R35 | Kędzierska B, Stodolna A, Bryszkowska K, Dylewski M, Potrykus K. A simple and unified protocol to purify all seven Escherichia coli RNA polymerase sigma factors. *J Appl Genet*. 2024;65(3):615–625. | 38709457 | 10.1007/s13353-024-00870-3 | PubMed abstract |
| R36 | Haldenwang WG. The sigma factors of Bacillus subtilis. *Microbiol Rev*. 1995;59(1):1–30. | 7708009 | 10.1128/mr.59.1.1-30.1995 | PubMed abstract |
| R37 | Kroos L, Yu YT. Regulation of sigma factor activity during Bacillus subtilis development. *Curr Opin Microbiol*. 2000;3(6):553–560. | 11121773 | 10.1016/S1369-5274(00)00140-5 | PubMed esummary |
| R38 | Giovannoni SJ, Cameron Thrash J, Temperton B. Implications of streamlining theory for microbial ecology. *ISME J*. 2014;8(8):1553–1565. | 24739623 | 10.1038/ismej.2014.60 | PubMed abstract |
| R39 | Battesti A, Majdalani N, Gottesman S. The RpoS-mediated general stress response in Escherichia coli. *Annu Rev Microbiol*. 2011;65:189–213. | 21639793 | 10.1146/annurev-micro-090110-102946 | PubMed esummary |
| R40 | Šmídová K, Ziková A, Pospíšil J, Schwarz M, Bobek J, Vohradsky J. DNA mapping and kinetic modeling of the HrdB regulon in Streptomyces coelicolor. *Nucleic Acids Res*. 2019;47(2):621–633. | 30371884 | 10.1093/nar/gky1018 | PubMed abstract + texto completo PMC |
| R41 | Fraser CM, et al. The minimal gene complement of Mycoplasma genitalium. *Science*. 1995;270(5235):397–403. | 7569993 | 10.1126/science.270.5235.397 | PubMed abstract |
| R42 | Schnoes AM, Brown SD, Dodevski I, Babbitt PC. Annotation error in public databases: misannotation of molecular function in enzyme superfamilies. *PLoS Comput Biol*. 2009;5(12):e1000605. | 20011109 | 10.1371/journal.pcbi.1000605 | Europe PMC + abstract |
| R43 | Tatusova T, DiCuccio M, Badretdin A, Chetvernin V, Nawrocki EP, Zaslavsky L, Lomsadze A, Pruitt KD, Borodovsky M, Ostell J. NCBI prokaryotic genome annotation pipeline. *Nucleic Acids Res*. 2016;44(14):6614–6624. | 27342282 | 10.1093/nar/gkw569 | PubMed esummary |
| R44 | Haft DH, et al. RefSeq: an update on prokaryotic genome annotation and curation. *Nucleic Acids Res*. 2018;46(D1):D851–D860. | 29112715 | 10.1093/nar/gkx1068 | PubMed esummary |
| R45 | Emms DM, Kelly S. OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biol*. 2019;20(1):238. | 31727128 | 10.1186/s13059-019-1832-y | OpenAlex |
| R46 | Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol*. 2021;38(12):5825–5829. | 34597405 | 10.1093/molbev/msab293 | OpenAlex |
| R47 | Kuznetsov D, Tegenfeldt F, Manni M, Seppey M, Berkeley M, Kriventseva EV, Zdobnov EM. OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity. *Nucleic Acids Res*. 2023;51(D1):D445–D451. | 36350662 | 10.1093/nar/gkac998 | Crossref + OpenAlex |
| R48 | Mistry J, Chuguransky S, Williams L, Qureshi M, Salazar GA, Sonnhammer ELL, Tosatto SCE, Paladin L, Raj S, Richardson LJ, Finn RD, Bateman A. Pfam: the protein families database in 2021. *Nucleic Acids Res*. 2021;49(D1):D412–D419. | 33125078 | 10.1093/nar/gkaa913 | Crossref + OpenAlex |
| R49 | Eddy SR. Accelerated profile HMM searches. *PLoS Comput Biol*. 2011;7(10):e1002195. | 22039361 | 10.1371/journal.pcbi.1002195 | Crossref + OpenAlex |
| R50 | Chávez J, Devos DP, Merino E. Complementary tendencies in the use of regulatory elements (transcription factors, sigma factors, and riboswitches) in bacteria and archaea. *J Bacteriol*. 2020 Dec 18;203(2):e00413-20. | 33077635 | 10.1128/jb.00413-20 | OpenAlex |
| R51 | Zuber U, Drzewiecki K, Hecker M. Putative sigma factor SigI (YkoZ) of Bacillus subtilis is induced by heat shock. *J Bacteriol*. 2001;183(4):1472–1475. | 11157964 | 10.1128/JB.183.4.1472-1475.2001 | PubMed abstract |
| R52 | Asai K, Ootsuji T, Obata K, Matsumoto T, Fujita Y, Sadaie Y. Regulatory role of RsgI in sigI expression in Bacillus subtilis. *Microbiology*. 2007;153(Pt 1):92–101. | 17185538 | 10.1099/mic.0.29239-0 | PubMed abstract |

**Nota sobre R08:** *J Mol Microbiol Biotechnol* 2002;4(1):77–91 tiene PMID
11763973 y no tiene DOI registrado en las fuentes consultadas; se verificaron el
esummary de PubMed y el PDF completo. Es la fuente primaria del censo comparativo de
31 genomas y de la afirmación sobre la nomenclatura confusa. PubMed registra un
erratum en 4(5):514.

**Nota sobre R13:** se cita por DOI verificado en Crossref; el PMID no se usó en este
documento.

---

## Anexo A. Registro de verificación

**Fecha:** 2026-11-26 (America/Santiago).

**Herramientas:**

- Europe PMC REST (`/search`, `resultType=lite` y `core`).
- NCBI E-utilities (`esearch`, `esummary`, `efetch`, `db=pubmed` y `db=pmc`).
- Crossref REST (`/works/{DOI}`, `/works?query.title=...`).
- OpenAlex (`/works/{DOI}`, filtros por DOI).
- InterPro/Pfam API (`/api/entry/pfam/{PF...}`).

**Controles aplicados:**

1. Cada PMID fue resuelto y el título devuelto se comparó con la cita.
2. Cada DOI fue resuelto y el título/año se comparó.
3. Cuando el abstract no contenía la afirmación, se buscó el texto completo
   (PMC/PDF) o se restringió la afirmación a lo que el abstract sostiene.
4. Se descartaron candidatos con PMID incorrecto. Ejemplos de correcciones:
   - Helmann & Chamberlin 1988: el PMID 3052287 corresponde a otro artículo;
     el correcto es **3052291**.
   - Wösten 1998: el PMID 9812375 corresponde a otro artículo; el correcto es
     **9818380**.
   - Schnoes 2009: el PMID 19390634 corresponde a otro artículo; el correcto es
     **20011109**.
   - Koo 2009: el PMID 19400810 corresponde a otro artículo; el correcto es
     **19400790**.
   - Yang 2015: el PMID 26229118 corresponde a otro artículo; el correcto es
     **26293966**.
   - OrthoFinder: el DOI 10.1186/s13059-019-1812-2 corresponde a MetaCell; el
     correcto es **10.1186/s13059-019-1832-y**.
5. No se incluyeron preprints, tesis ni páginas web no institucionales.
6. Las entradas de Pfam se citan por su accession verificado; el nombre y la
   descripción se transcriben de la entrada de InterPro.

**Cobertura:** 52 referencias primarias/revisiones + 11 entradas de InterPro/Pfam.
**Referencias no verificadas incluidas:** 0.

---

## Anexo B. Entradas InterPro/Pfam verificadas

| Accession | Nombre verificado | Tipo | Definición relevante | URL |
|---|---|---|---|---|
| PF00140 | Sigma-70 factor, region 1.2 | family | Región 1.2; marcador de grupos 1–2, no cubre ECF | https://www.ebi.ac.uk/interpro/entry/pfam/PF00140/ |
| PF04542 | Sigma-70 region 2 | domain | Región más conservada; –10 y unión al core | https://www.ebi.ac.uk/interpro/entry/pfam/PF04542/ |
| PF04539 | Sigma-70 region 3 | family | Hélice 3; –10 extendida y unión al core | https://www.ebi.ac.uk/interpro/entry/pfam/PF04539/ |
| PF04545 | Sigma-70, region 4 | domain | HTH de unión a –35; umbral alto, falsos negativos | https://www.ebi.ac.uk/interpro/entry/pfam/PF04545/ |
| PF04997 | RNA polymerase Rpb1, domain 1 | domain | Dominio clamp de β′; **no es sigma** | https://www.ebi.ac.uk/interpro/entry/pfam/PF04997/ |
| PF07638 | ECF sigma factor | family | Sigma ECF | https://www.ebi.ac.uk/interpro/entry/pfam/PF07638/ |
| PF00309 | Sigma-54 factor, AID | family | Dominio de interacción con activador | https://www.ebi.ac.uk/interpro/entry/pfam/PF00309/ |
| PF04552 | Sigma-54, DNA binding domain | domain | Unión a –24/–12 | https://www.ebi.ac.uk/interpro/entry/pfam/PF04552/ |
| PF04963 | Sigma-54 factor, core binding domain | domain | Unión al core | https://www.ebi.ac.uk/interpro/entry/pfam/PF04963/ |
| PF00158 | Sigma-54 interaction domain | domain | Dominio de bEBPs; no es sigma | https://www.ebi.ac.uk/interpro/entry/pfam/PF00158/ |
| PF12791 | Anti-sigma factor N-terminus | domain | Anti-sigma (contexto) | https://www.ebi.ac.uk/interpro/entry/pfam/PF12791/ |

**Verificación:** consultas directas a `https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}/`
y a la búsqueda `?search=sigma-54`, `?search=ECF`, `?search=anti-sigma`, el
2026-11-26.

---

## Anexo C. Estrategia de búsqueda y límites

**Bases:** PubMed/Europe PMC, Crossref, OpenAlex, InterPro/Pfam.

**Términos principales:**

- `bacterial sigma factors review`
- `sigma70 family ECF sigma54`
- `sigma factor promoter architecture -35 -10`
- `sigma54 -24 -12 promoter`
- `ECF sigma factor classification`
- `sigma factor orthologs nomenclature`
- `Bacillus subtilis sigma factors`
- `Streptomyces coelicolor sigma factors`
- `Mycobacterium tuberculosis sigma factors`
- `Pseudomonas aeruginosa sigma factors`
- `Escherichia coli seven sigma factors`
- `sigma factor annotation error`
- `sigma factor misannotation`
- `Pfam sigma70 domain`

**Criterios de inclusión:**

- Artículos revisados por pares, en inglés.
- Revisiones o estudios primarios con datos de familia, estructura, promotor,
  distribución u ortología.
- Bases de datos curadas (Pfam/InterPro) para definiciones de dominios.

**Criterios de exclusión:**

- Preprints no revisados (no citados).
- Tesis y páginas web no institucionales (no citadas).
- Artículos sin metadata resoluble.
- Afirmaciones que solo podían sostenerse con texto de pago no accesible.

**Límites:**

- No se accedió a todos los textos completos; en esos casos la afirmación se limita
  al abstract o a la información verificable.
- Los censos de 2002 reflejan la anotación de entonces; se marcan como históricos.
- La búsqueda no es sistemática; los vacíos del §10 son ausencias en la literatura
  recuperada, no prueba de inexistencia.
- No se hizo filogenia propia; las afirmaciones de ortología se apoyan en las fuentes
  citadas y se marcan como “alta”, “baja” o “no verificada”.

---

## Anexo D. Matriz de afirmaciones clave

| # | Afirmación | Cita | Tipo |
|---|---|---|---|
| 1 | σ70 y σ54 son familias estructuralmente no relacionadas | [PMID:26131973] [PMID:11763973] | Revisión + censo |
| 2 | La familia σ70 tiene cuatro grupos | [PMID:26131973] [PMID:12540296] | Revisión |
| 3 | σ2 reconoce –10 y une el core | [PMID:26131973] [PMID:22136875] | Revisión + estructura |
| 4 | σ4 reconoce –35 mediante HTH | [PMID:11931761] [PMID:26131973] | Estructura + revisión |
| 5 | σ3 reconoce la –10 extendida | [PMID:11931761] [PMID:26131973] | Estructura + revisión |
| 6 | σ1.1 actúa como compuerta | [PMID:26131973] | Revisión |
| 7 | El promotor σ70 canónico usa –35/–10, UP y discriminador | [PMID:18521075] [PMID:26131973] | Revisión |
| 8 | σ54 usa –24/–12 y requiere bEBP + ATP | [PMID:10536136] [PMID:22933558] | Estudio + revisión |
| 9 | σ28 reconoce una –10 larga `GCCGATAA` | [PMID:19400790] | Estudio primario |
| 10 | σS depende de –13/–14 | [PMID:11251833] | Estudio primario |
| 11 | σE de *E. coli* es ECF y tiene modelo de promotor validado | [PMID:16336047] | Estudio primario |
| 12 | Las ECF son el grupo más numeroso de sigma alternativas | [PMID:19737356] | Revisión |
| 13 | El censo ECF pasó de ~6 a ~10 por genoma y de >40 a 157 grupos | [PMID:19737356] [PMID:33398323] | Revisiones |
| 14 | *E. coli* K-12 tiene 7 sigma | [PMID:38709457] [PMID:18070067] [PMID:11763973] | Estudio + revisión + censo |
| 15 | *B. subtilis* tiene 17–18 sigma | [PMID:14993308] [PMID:11763973] | Estudio + censo |
| 16 | *M. tuberculosis* tiene 13 sigma | [PMID:17064287] [PMID:11763973] | Revisión + censo |
| 17 | *P. aeruginosa* tiene 24 sigma, 19 ECF | [PMID:18070067] [PMID:11763973] | Revisión + censo |
| 18 | *S. coelicolor* tiene 66 sigma (65 alternativas + HrdB) | [PMID:30371884] | Estudio primario |
| 19 | *M. genitalium* tiene 1 sigma | [PMID:30371884] [PMID:11763973] | Estudio + censo |
| 20 | RpoD/SigA/HrdB son ortólogos housekeeping | [PMID:30371884] [PMID:11763973] | Estudio + censo |
| 21 | FliA/SigD son ortólogos flagelares | [PMID:11763973] | Censo |
| 22 | RpoN/SigL son ortólogos σ54 | [PMID:11763973] | Censo |
| 23 | RpoH está restringida a proteobacterias en el censo de 2002 | [PMID:11763973] | Censo |
| 24 | σE de *E. coli* y σE de *B. subtilis* no son ortólogas | [PMID:16336047] [PMID:11121773] | Estudio + revisión |
| 25 | La nomenclatura de sigma es confusa y requiere revisión | [PMID:11763973] | Censo |
| 26 | La anotación automática tiene errores de 5–63 % en superfamilias | [PMID:20011109] | Estudio |
| 27 | PGAP y RefSeq anotan por HMM/homología y curan parcialmente | [PMID:27342282] [PMID:29112715] | Artículos de pipeline |
| 28 | Pfam PF04545 tiene umbral alto y falsos negativos | [InterPro:PF04545] | Base de datos |
| 29 | Pfam PF04997 no es un dominio sigma | [InterPro:PF04997] | Base de datos |
| 30 | OrthoFinder/eggNOG/OrthoDB infieren ortología, no nombres | [PMID:31727128] [PMID:34597405] [PMID:36350662] | Métodos |

---

## Anexo E. Glosario mínimo

- **Sigma (σ):** subunidad disociable de RNAP que confiere especificidad de promotor.
- **Familia σ70:** familia homóloga a σ70/RpoD; incluye grupos 1–4.
- **ECF:** *extracytoplasmic function*; grupo 4 de σ70, típicamente pequeñas y
  reguladas por anti-sigma.
- **σ54:** familia separada, dependiente de activadores bEBP y ATP.
- **Grupo 1:** sigma housekeeping, esencial.
- **Grupo 2:** sigma primaria no esencial (RpoS; SigB según clasificación).
- **Grupo 3:** sigma alternativas (RpoH, FliA, esporulación).
- **Anti-sigma:** proteína que secuestra una sigma e impide su unión a RNAP.
- **bEBP:** *bacterial enhancer binding protein*; activador AAA+ de σ54.
- **Orthogroup:** conjunto de genes que descienden de un gen ancestral común por
  especiación; se infiere por filogenia.
- **Parálogo:** gen duplicado; no es un ortólogo 1:1.
- **PWM:** *position weight matrix*; modelo de motivo.
- **HMM:** *hidden Markov model*; modelo de perfil para dominios proteicos.

---

*Fin del documento. No se realizó commit ni se modificó código.*
