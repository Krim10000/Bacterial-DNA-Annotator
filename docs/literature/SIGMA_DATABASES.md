# Bases de datos curadas de TSS, promotores, sigma y regulación bacteriana

> **Estado:** revisión documental, corte **2026-11-26**. No se modificó código
> ni se hizo commit. Este documento audita fuentes; no re-anota el gold ni
> cambia `data/sigma/gold/`.
>
> **Pregunta:** ¿qué base de datos puede aportar pares **TSS–sigma** con
> evidencia y licencia utilizables como gold para el catálogo de 400?
>
> **Respuesta corta:** solo cuatro recursos generales contienen TSS y sigma en
> la misma entidad: **RegulonDB**, **DBTBS**, **BSGatlas** y **EcoCyc** (a los
> que se suma **Pro54DB**, específico de σ54). De ellos, solo **RegulonDB** y
> **DBTBS** son anclas primarias; **BSGatlas** deriva de DBTBS y de un estudio
> de tiling array, y **EcoCyc** comparte el equipo de curación de RegulonDB. El resto no sirve como gold σ–TSS: o tiene TSS sin sigma (PPD,
> MtbRegList), o sigma sin TSS (PRODORIC), o solo TF binding (CollecTF,
> RegPrecise, CoryneRegNet, SubtiWiki, P2TF), o es predicción/anotación
> (CDBProm, MicroScope, P2TF).

---

## 0. Convenciones y definiciones

| Término | Definición operativa en este documento |
|---|---|
| **TSS** | Coordenada 1-based del primer nucleótido transcrito, con hebra y contig. No es el centro del promotor ni el sitio de unión de un TF. |
| **Promotor** | Región que recluta la RNAP holoenzima; puede tener −35/−10, −10 extendida o −24/−12. Un promotor puede tener más de un TSS. |
| **Asignación σ** | Factor sigma que la fuente atribuye a ese TSS/promotor. Es una etiqueta curada o inferida, no una consecuencia automática de la secuencia. |
| **TF binding** | Ocupación de DNA por un regulador secuencia-específico. **No implica sigma** y **no se deriva de sigma**. |
| **Regulón** | Conjunto de genes/operones co-regulados por un TF (o por un factor sigma, pero eso es otra capa). |
| **Gold σ–TSS** | Fila con genoma/contig, posición TSS, hebra, sigma nativo, evidencia y PMID. Sin esos seis campos no es gold utilizable. |

Códigos usados en las tablas:

- TSS / σ: `✅` explícito y descargable; `◐` indirecto o incompleto; `❌` ausente.
- Evidencia: `E` experimental curada; `H` alto rendimiento; `P` predicción;
  `M` mixta (experimental + predicción/herencia).
- Uso como gold: `A` ancla primaria; `B` verificación/cross-check (no
  independiente); `C` solo contexto, nunca etiqueta σ; `D` no usar.

---

## 1. Tabla comparativa: recursos con TSS, promotor o sigma

| Base | Especies | TSS | σ | Arquitectura | Condición | Evidencia | Versión / fecha | Licencia | Gold |
|---|---|---|---|---|---|---|---|---|---|
| **RegulonDB** | *E. coli* K-12 MG1655 | ✅ coord. + evidencia | ✅ por promotor (nativo) | ✅ cajas −10/−35 | ✅ parcial (growth conditions) | `E`+`H`; confianza S/W/C | v12.0 (2024-01); **v14.5** (2026-02-26) | Datos v12.0: **CC BY-NC-SA 4.0**; artículo CC BY-NC 4.0 | **A** |
| **DBTBS** | *B. subtilis* 168 | ✅ “Start position” (declarada) | ✅ SigA/SigB/SigK/… | ✅ −35/−10 | ❌ | `E` (literatura) | Release 5 (2007-09); flat file 2008-10-21 | **Sin licencia de datos explícita**; artículo CC BY-NC 2.0 UK | **A** (congelado) |
| **BSGatlas** | *B. subtilis* 168 | ✅ 3.390 TSS unificados | ✅ heredado, con `Resolution` | ◐ (mapa de TSS, no cajas) | ❌ | `M` (curado + tiling array + propagación) | v1.1 (2020-12; paper 2021-02) | **CC BY 4.0** (GFF3 y artículo) | **B** |
| **EcoCyc** | *E. coli* K-12 MG1655 | ✅ promotores/TSS | ✅ | ✅ | ✅ nutrient conditions/efectores | `E`+`H` | EcoCyc 2025; RegulonDB v14.5 cita EcoCyc 29.5 | **BioCyc**: licencia obligatoria; datos EcoCyc gratuitos con licencia; comercial por suscripción | **B** |
| **PRODORIC** | 27 organismos | ◐ solo referencia al TSS | ✅ matrices −35/−10 | ✅ | ✅ efectores y tipo de regulación | `E` para TFBS; σ derivado de TSS | Release 2021 (paper 2022-01) | **CC BY-NC 4.0** (datos, explícito) | **C** |
| **PPD** | 63 especies / 74 cepas | ✅ TSS | ❌ | ◐ solo secuencia | ❌ | `E` (curada de papers) | v1.0; última actualización 2020-09-22 | Gratuita, **© LinGroup**, sin licencia abierta | **C** |
| **MtbRegList** | *M. tuberculosis* | ✅ TSS experimental | ❌ | ✅ motivos | ❌ | `M` (predicho + caracterizado + TSS) | 2005; sitio inaccesible (403/WAF) | Sin licencia | **D** |
| **Pro54DB** (extra) | 43 especies | ◐ según entrada | ✅ σ54 implícito | ✅ −24/−12 | ❌ | `E` (210 promotores) | v1.0; 2016-08-31 | Sin licencia explícita | **C** (solo σ54) |

## 2. Tabla comparativa: recursos de regulación/TF (sin TSS–sigma)

| Base | Especies | Qué dato | TF binding | Condición | Evidencia | Versión / fecha | Licencia | Gold |
|---|---|---|---|---|---|---|---|---|
| **CollecTF** | Bacteria; >100 especies | TFBS experimentales | ✅ | ❌ | `E` | 2014; update 2016 | Artículos CC BY 3.0/4.0; datos sin licencia explícita | **C** |
| **RegPrecise** | >160 genomas | regulones reconstruidos | ✅ (predichos) | ✅ efectores | `P` (genómica comparativa) | v3.2; paper 2013-11 | Artículo CC BY 2.0; datos sin licencia explícita | **D** |
| **CoryneRegNet 7** | *Corynebacterium* (+4 modelo) | redes TF–gen | ✅ | ❌ | `M` (experimental + predicho) | 7.0; 2020-05 | Artículo CC BY 4.0; datos sin licencia explícita | **C** |
| **SubtiWiki** | *B. subtilis* 168 | gen/proteína/regulación/regulón | ✅ (relaciones gen–gen) | ✅ expresión | `E` (curaduría manual) | v5; paper 2025-01 | Wiki CC BY; artículo CC BY-NC 4.0; API v5 sin aviso específico | **C** |
| **P2TF** | 1.987 genomas + 43 metagenomas | proteínas TF predichas | ❌ (no sitios) | ❌ | `P` (arquitectura de dominios) | 2012-11; URL muerta | Artículo CC BY 2.0; datos sin licencia explícita | **D** |
| **MicroScope** | >11.800 genomas | anotación genómica comparativa | ❌ | ◐ fenotipos/transcriptómica | `M` (automática + experta) | 2020-01 | Artículo CC BY 4.0; términos de datos de la plataforma no verificados | **D** |
| **CDBProm** | >6.000 genomas | promotores predichos | ❌ | ❌ | `P` (XGBoost) | 2024-02 | Artículo CC BY 4.0; descargas sin licencia declarada | **D** |

> **Regla de separación.** Las columnas “σ” y “TF binding” son ortogonales.
> Un TFBS no etiqueta sigma; una matriz sigma no identifica el TF. Varias
> bases (DBTBS, PRODORIC, P2TF) mezclan “sigma factors” dentro de su lista de
> “transcription factors”; contarlas como TFs o usarlas para inferir sigma es
> un error de categoría.

---

## 3. Qué base sirve para sigma–TSS y cuál no

| Base | ¿TSS? | ¿σ? | ¿Ambos en la misma fila? | ¿Independiente del gold local? | Veredicto |
|---|---|---|---|---|---|
| RegulonDB | ✅ | ✅ | ✅ | Sí (ancla) | **A** |
| DBTBS | ✅ | ✅ | ✅ | Sí (ancla, pero 2008) | **A** |
| BSGatlas | ✅ | ✅ | ✅ | No: deriva de DBTBS/BsubCyc/Nicolas | **B** |
| EcoCyc | ✅ | ✅ | ✅ | No: misma curación UNAM que RegulonDB | **B** |
| PRODORIC | ◐ | ✅ | ❌ | Parcial | **C** |
| PPD | ✅ | ❌ | ❌ | Parcial | **C** |
| MtbRegList | ✅ | ❌ | ❌ | No verificable (caída) | **D** |
| Pro54DB | ◐ | ✅ (σ54) | ◐ (solo σ54; coord. no verificada) | Parcial | **C** |
| CollecTF | ❌ | ❌ | ❌ | — | **C** |
| RegPrecise | ❌ | ❌ | ❌ | — | **D** |
| CoryneRegNet | ❌ | ❌ | ❌ | — | **C** |
| SubtiWiki | ❌ | ❌ | ❌ | — | **C** |
| P2TF | ❌ | ❌ | ❌ | — | **D** |
| MicroScope | ❌ | ❌ | ❌ | — | **D** |
| CDBProm | ❌ (predicho) | ❌ | ❌ | — | **D** |

**Lectura para el catálogo de 400:** solo **RegulonDB** y **DBTBS** permiten
construir un gold σ–TSS defendible. **BSGatlas** y **EcoCyc** sirven para
verificar consistencia, no para declarar independencia. **PRODORIC** y **PPD**
cubren piezas sueltas (matriz sigma / TSS sin sigma) y no se pueden unir sin
crear una etiqueta nueva. El resto es contexto TF o predicción.

---

## 4. Fichas por base

### 4.1 RegulonDB — *E. coli* K-12

- **Qué tiene.** Promotores con `posTSS`, hebra, sigma nativo, cajas −10/−35,
  gen/operón, evidencia y PMIDs; además TFBS/TFRS, unidades de transcripción,
  terminadores y condiciones de crecimiento. v12.0 reporta **4.050
  promotores**, **238 TFs** y **6.110 interacciones regulatorias**.
- **Especies.** *E. coli* K-12 MG1655; es una base de una sola especie.
- **Evidencia.** Curaduría de literatura de bajo rendimiento + colecciones HT
  (16 datasets de TSS, 22 de ChIP de sigma, 6 de RpoB). Cada promotor trae
  `confidenceLevel` S/W/C y evidencia; no todo es experimental de la misma
  calidad.
- **Versión.** Paper v12.0 (2024-01); release actual **v14.5 / v14.5.1**
  (2026-02-26), sincronizada con EcoCyc 29.5.
- **Licencia.** Datos v12.0.0 en Zenodo: **CC BY-NC-SA 4.0**; artículo
  CC BY-NC 4.0; el código tiene licencia propia en GitHub.
- **Limitación para gold.** (i) Una sola especie; (ii) S/W/C no es una
  probabilidad; (iii) la colección HT puede incluir ruido; (iv) el sigma se
  separa por promotor cuando dos sigma comparten TSS, así que la clave
  TSS+sigma es la correcta, no TSS sola; (v) **comparte curación con EcoCyc**.
- **Refs.** [R1], [R2], [R3].

### 4.2 DBTBS — *B. subtilis* 168

- **Qué tiene.** Entradas de promotor con “Start position”, dirección, factor
  de unión (incluye sigma, p. ej. `SigK`), cajas −35/−10, gen/operón,
  referencias y año; además TFBS, operones y terminadores. Release 5:
  **1.475 promotores**, **120 TFs** (incluye sigma), **736 operones**,
  **463 terminadores**, 947 publicaciones.
- **Especies.** *B. subtilis* 168.
- **Evidencia.** Curaduría de literatura; el paper lo describe como
  “experimentally characterized”. No hay niveles de confianza.
- **Versión.** Release 5 (2007-09); flat file `dbtbs.xml` actualizado
  **2008-10-21**. BSGatlas (2021) lo describe como proyecto discontinuado.
- **Licencia.** **No hay licencia de datos explícita** en el sitio; el artículo
  es CC BY-NC 2.0 UK. La descarga exige nombre, email, institución y propósito
  (formulario), por lo que no es automatizable ni redistribuible sin permiso.
- **Limitación para gold.** (i) Congelado desde 2008; (ii) sin confianza por
  fila; (iii) el número de entradas de promotor (1.475) no equivale a TSS
  únicos (el scrape local recuperó 865 TSS de 1.160 páginas de operón);
  (iv) licencia ausente.
- **Refs.** [R4], [R5].

### 4.3 BSGatlas — *B. subtilis* 168

- **Qué tiene.** Anotación integrada: **3.390 TSS**, 2.566 TTS, 2.266
  operones y 5.335 UTRs. Cada TSS lleva `Sigma`, `Resolution` y `Based on`
  (DBTBS, BsubCyc, Nicolas et al., SubtiWiki, etc.).
- **Especies.** *B. subtilis* 168 (RefSeq ASM904v1).
- **Evidencia.** **Mixta**: TSS curados de DBTBS/BsubCyc, TSS de tiling array
  de Nicolas et al. (predicción con resolución limitada) y unificación. Solo
  **12,4 %** de los TSS están compartidos entre las fuentes de origen; el sigma
  se toma de la entrada de menor resolución del grupo unificado. De los sitios
  sigma de DBTBS con posición relativa al TSS, se infirieron **644 TSS de alta
  resolución**.
- **Versión.** v1.1 (Zenodo 2020-12-04; paper 2021-02).
- **Licencia.** **CC BY 4.0** para el GFF3 (Zenodo) y el artículo.
- **Limitación para gold.** (i) Es una integración, no una medición
  independiente; (ii) el sigma puede venir de Nicolas et al. (tiling array),
  no de DBTBS; (iii) no trae cajas −10/−35; (iv) `Resolution` no es una
  probabilidad.
- **Refs.** [R6], [R7].

### 4.4 EcoCyc — *E. coli* K-12

- **Qué tiene.** Promotores, TSS, unidades de transcripción, interacciones
  regulatorias, terminadores, condiciones nutricionales y efectores. Desde la
  revisión de definiciones, cada TSS con distinto sigma es un promotor
  distinto, aunque la coordenada sea la misma; un promotor puede tener varios
  TSS en ~5 bp.
- **Especies.** *E. coli* K-12 MG1655 (con propagación a cepas de *E. coli*).
- **Evidencia.** Curaduría manual + HT; tabla 2021→24.5: 197 promotores nuevos
  y 228 modificados.
- **Versión.** EcoCyc 2025 (EcoSal Plus 2025); RegulonDB v14.5 reporta
  sincronización con EcoCyc 29.5.
- **Licencia.** **BioCyc**: el acceso a archivos de datos requiere licencia y
  suscripción; excepción: los archivos de **EcoCyc son gratuitos, pero con
  licencia**; uso comercial, suscripción. El artículo de 2021 es CC BY 4.0; el
  de 2025 es de suscripción (verificar términos antes de reutilizar figuras o
  texto).
- **Limitación para gold.** (i) **No es independiente de RegulonDB**: el paper
  de RegulonDB dice que la curación de regulación de ambas bases viene del
  mismo equipo de la UNAM y que v12.0 está sincronizada con EcoCyc 27.0;
  (ii) licencia restrictiva para redistribución; (iii) una sola especie.
- **Refs.** [R10], [R11], [R12].

### 4.5 PRODORIC — múltiples bacterias

- **Qué tiene.** TFBS experimentales, matrices PWM, efectores, tipo de
  regulación y evidencia; **27 organismos**, **362 TFs** y **3.238 TFBS** en
  el release 2021. Los sigma se guardan como una entrada con −35 y −10
  separadas por un espaciador variable.
- **Especies.** 27 organismos (mayoría modelos bacterianos); no es exhaustivo.
- **Evidencia.** TFBS re-anotados con evidencia experimental; el paper advierte
  que **la mayoría de los sitios sigma se definen a partir de TSS determinados
  por primer extension y pueden no tener evidencia experimental propia**.
- **Versión.** Release 2021; paper 2022-01.
- **Licencia.** **CC BY-NC 4.0** para los datos (explícito en el paper).
- **Limitación para gold.** (i) No es un catálogo TSS-level: el TSS es una
  referencia para definir el sitio sigma, no una coordenada descargable con
  evidencia por fila; (ii) no cubre las 400 especies; (iii) la matriz sigma no
  dice qué TSS existe.
- **Refs.** [R13].

### 4.6 PPD — 63 especies

- **Qué tiene.** **129.148 promotores** experimentales en **63 especies / 74
  cepas**, con posición TSS, hebra, nombre, gen regulado, locus tag, secuencia
  y referencias; visor genómico, búsqueda y descarga CSV.
- **Especies.** 63 especies procariotas (bacterias y algunas arqueas); incluye
  *E. coli*, *B. subtilis*, *P. aeruginosa* PAO1 (solo 10 promotores en la
  estadística del sitio), *S. aureus* MW2 (2.821) y *C. glutamicum* (3.581).
  **No** incluye *M. tuberculosis* ni *S. coelicolor*.
- **Evidencia.** Curaduría manual de publicaciones; el sitio se describe como
  “experimentally verified”.
- **Versión.** v1.0; última actualización **2020-09-22**.
- **Licencia.** Acceso libre y sin registro, pero **© LinGroup, all rights
  reserved**; no hay licencia abierta verificada.
- **Limitación para gold.** (i) **No tiene sigma**: el CSV descargable trae
  `id, Type, PromoterName, SpeciesName, GeneName, GeneStart, GeneEnd,
  GeneProduct, TSSPosition, Strand, Locus_tag, PromoterSeq`; la búsqueda
  avanzada solo ofrece nombre, secuencia, TSS, producto y locus tag; (ii)
  cobertura muy desigual por especie; (iii) sin condiciones; (iv) licencia no
  abierta.
- **Refs.** [R15], [R16].

### 4.7 MtbRegList — *M. tuberculosis*

- **Qué tiene.** Motivos regulatorios predichos y caracterizados,
  cross-referenciados con TFs, y **TSS identificados experimentalmente**;
  búsquedas genómicas flexibles.
- **Especies.** *M. tuberculosis*.
- **Evidencia.** Mixta: motivos predichos + caracterizados + TSS
  experimentales.
- **Versión.** 2005; sin actualización posterior verificada.
- **Licencia.** Sin licencia de datos; artículo de suscripción.
- **Limitación para gold.** (i) **No se encontró campo sigma** y el abstract
  no lo menciona; al no poder acceder al sitio, la ausencia no es verificable
  al 100 %; (ii) el sitio original (`USherbrooke.ca/vers/MtbRegList`) ya no
  responde (403/WAF, 2026-11-26); (iii) es una base de 2005, no apta para un
  gold actual.
- **Refs.** [R19].

### 4.8 Pro54DB — σ54 (extra)

- **Qué tiene.** **210 promotores σ54** experimentales, 297 genes regulados,
  **43 especies**, 133 publicaciones; promotor, gen y referencias.
- **Especies.** 43 especies bacterianas.
- **Evidencia.** Experimental curada de literatura.
- **Versión.** v1.0; última actualización **2016-08-31**.
- **Licencia.** Sin licencia explícita; © LinGroup.
- **Limitación para gold.** (i) Solo σ54; (ii) pequeño y antiguo; (iii) no se
  verificó que cada entrada tenga coordenada TSS genómica descargable;
  (iv) útil como set específico de σ54, no como gold multiespecie.
- **Refs.** [R24].

### 4.9 CollecTF — TFBS en Bacteria

- **Qué tiene.** **>9.750 TFBS experimentales**, >390 instancias de TF, >240
  TFs y >100 especies; integra evidencia experimental, efectos regulatorios y
  motivos; exporta a RefSeq.
- **Especies.** Bacteria (amplio).
- **Evidencia.** Solo experimental; los sitios solo in silico no se reportan.
- **Versión.** 2014; actualización 2016; portal actualmente inestable
  (`collectf.org` sirve un shell JS; el dominio antiguo de UMBC no responde).
- **Licencia.** Artículos CC BY 3.0 / CC BY 4.0; **sin licencia de datos
  explícita**.
- **Limitación para gold.** (i) **No tiene TSS ni sigma**; (ii) TFBS ≠
  promotor; (iii) portal inestable.
- **Refs.** [R21], [R22].

### 4.10 RegPrecise — regulones predichos

- **Qué tiene.** **781 regulogs de TF** en **>160 genomas**, con regulones
  reconstruidos por genómica comparativa, efectores, familias de TF y
  propagación a genomas cercanos.
- **Especies.** Múltiples taxones bacterianos.
- **Evidencia.** **Predicción** (genómica comparativa) a partir de regulones
  de referencia curados; la propagación es automática.
- **Versión.** v3.2; paper 2013-11.
- **Licencia.** Artículo CC BY 2.0; datos sin licencia explícita.
- **Limitación para gold.** (i) **No tiene TSS ni sigma**; (ii) es predicción,
  no medición; (iii) no usar como gold de σ–TSS, solo como contexto de
  regulones.
- **Refs.** [R23].

### 4.11 CoryneRegNet 7 — corinebacterias

- **Qué tiene.** TFBS y redes TF–gen para *Corynebacterium*; incluye una capa
  “experimental” para *C. glutamicum* 13032, *B. subtilis* 168, *E. coli*
  K-12 y *M. tuberculosis* H37Rv, y una capa “predicha” para el género.
- **Especies.** Género *Corynebacterium* (+4 organismos modelo).
- **Evidencia.** Mixta: experimental en modelos, predicha en el género.
- **Versión.** 7.0 (2020-05).
- **Licencia.** Artículo CC BY 4.0; datos sin licencia explícita.
- **Limitación para gold.** **No tiene TSS ni sigma**; la capa predicha no es
  gold.
- **Refs.** [R20].

### 4.12 SubtiWiki — *B. subtilis*

- **Qué tiene.** Gen/proteína, categorías funcionales, interacciones,
  regulación gen–gen, regulones, expresión, estructuras y datos de fitness;
  API REST documentada (Swagger) con endpoints de gen, operón, regulación y
  regulón.
- **Especies.** *B. subtilis* 168.
- **Evidencia.** Curaduría manual de literatura; no es una base de mediciones
  TSS.
- **Versión.** v5 (paper 2025-01); API `/v5/api`.
- **Licencia.** Wiki histórica: CC BY (“Creative Commons Attribution”); el
  registro institucional indica CC BY 4.0; el artículo es CC BY-NC 4.0; la API
  v5 no muestra aviso de licencia específico.
- **Limitación para gold.** **No ofrece un catálogo TSS-level ni un campo
  sigma por promotor**: sus campos de regulación son `regulator_gene_id`,
  `regulated_gene_id`, `mechanism`, `mode`, `regulon_id`, `operon_id`. Es
  contexto TF, no sigma.
- **Refs.** [R8], [R9].

### 4.13 P2TF — inventario de TFs predichos

- **Qué tiene.** **372.877 TFs** de **1.987 genomas** completamente
  secuenciados y 43 metagenomas; clasificación por arquitectura de dominios,
  contexto génico y familia; incluye sigma factors como subconjunto.
- **Especies.** Procariontes (amplio).
- **Evidencia.** **Predicción** por dominios; no hay sitios de unión ni
  promotores.
- **Versión.** 2012-11; el dominio `p2tf.org` ya no sirve la base (parking/
  redirección, 2026-11-26).
- **Licencia.** Artículo CC BY 2.0; datos sin licencia explícita; URL muerta.
- **Limitación para gold.** (i) No tiene TSS ni sigma por promotor; (ii) es
  predicción de proteínas; (iii) no usar como gold. Solo inventario de
  candidatos sigma/TF.
- **Refs.** [R14].

### 4.14 MicroScope — anotación genómica

- **Qué tiene.** Anotación automática y experta de genomas microbianos,
  comparativa, metabolismo, pangenomas y datos post-genómicos; **>11.800
  genomas**. No ofrece un catálogo curado de promotores, TSS ni sigma.
- **Especies.** >11.800 genomas procariotas.
- **Evidencia.** `M`: pipeline automático + curación experta de funciones
  génicas; no es curación de TSS.
- **Versión.** 2020-01 (paper); plataforma en desarrollo continuo.
- **Licencia.** Artículo CC BY 4.0; términos de datos de la plataforma no
  verificados.
- **Limitación para gold.** No sirve para σ–TSS; sirve para contexto génico,
  pangenoma y funciones.
- **Refs.** [R18].

### 4.15 CDBProm — promotores predichos

- **Qué tiene.** **>24 millones de promotores predichos** (el sitio dice >12
  millones; discrepancia a verificar) en **>6.000 genomas**; XGBoost entrenado
  con 15.654 promotores curados; asociación promotor–gen.
- **Especies.** >6.000 genomas bacterianos.
- **Evidencia.** **Predicción**; sin sigma; sin validación experimental propia
  por promotor.
- **Versión.** 2024-02.
- **Licencia.** Artículo CC BY 4.0; descargas sin licencia declarada.
- **Limitación para gold.** (i) Es predicción, no medición; (ii) se entrenó y
  validó con promotores de RegulonDB (3.861), lo que crea circularidad directa
  si se usa como gold; (iii) no tiene sigma; (iv) no usar.
- **Refs.** [R17].

---

## 5. Riesgos de circularidad

### 5.1 Linaje compartido

```text
RegulonDB  <---- misma curación UNAM ---->  EcoCyc
     \                                        /
      \--> gold K-12 (3.314 pares) <---------/
                (sincronización v12.0 <-> EcoCyc 27.0; v14.5 <-> 29.5)

DBTBS ─┐
BsubCyc├──> BSGatlas v1.1 ──> gt_tss.tsv ──> gold "DBTBS" B. subtilis (2.226)
Nicolas┘
```

- **RegulonDB y EcoCyc no son dos evidencias independientes.** El paper de
  RegulonDB afirma que la curación de regulación de ambas bases viene del mismo
  equipo de la UNAM y que la versión 12.0 está sincronizada con EcoCyc 27.0.
  Usar una para entrenar y la otra para evaluar no elimina la fuga de
  curación.
- **BSGatlas no es una fuente primaria.** Integra DBTBS, BsubCyc, SubtiWiki y
  Nicolas et al.; solo 12,4 % de los TSS están compartidos entre fuentes. El
  sigma se hereda de la entrada de menor resolución. Evaluar en BSGatlas
  después de calibrar con DBTBS mide, en parte, la misma información.
- **DBTBS y BsubCyc están discontinuados** (lo dice BSGatlas 2021); cualquier
  “actualización” de B. subtilis pasa por BSGatlas, que es derivado.

### 5.2 Predicción promovida a anotación

- **CDBProm** es un clasificador; sus 24 M de promotores son predicciones. Se
  entrenó con promotores curados (incluida validación en 3.861 promotores de
  RegulonDB). Usarlo como gold sería evaluar el modelo contra su propia
  distribución.
- **RegPrecise** propaga regulones por genómica comparativa; la capa
  “propagated” es automática.
- **P2TF** infiere TFs por arquitectura de dominios; no mide unión.
- **MicroScope** anota funciones, no TSS.
- **BSGatlas**: 2.294 de sus 3.390 TSS (67,7 %) solo citan “Nicolas et al.”
  (tiling array); no son curaduría de literatura primaria.

### 5.3 Fuga en la evaluación

- **Misma PMID / mismo operón**: un TSS y su vecino pueden caer en train y
  test si se parte por posición. Hay que partir por operón/locus y por PMID.
- **Misma PWM**: si la PWM se estima con los mismos TSS que se usan para
  evaluar, el test no es ciego. Estimar PWMs solo en el fold de train.
- **Homología**: ortólogos entre especies comparten promotores; un split
  aleatorio por TSS no mide generalización.
- **Multi-sigma**: 124 TSS de K-12 tienen más de un sigma; una métrica top-1
  penaliza por diseño. Reportar multi-etiqueta.
- **Doble conteo**: el mismo TSS puede aparecer en RegulonDB, EcoCyc, BSGatlas
  y PPD con etiquetas distintas. Deduplicar por `(genoma, contig, pos, strand,
  sigma)` y conservar la procedencia.

### 5.4 Auditación local: el gold “DBTBS” es BSGatlas

| Archivo local | Filas | Procedencia real | Problema |
|---|---:|---|---|
| `data/dbtbs/dbtbs_tss.tsv` | 865 | DBTBS release 5 (scrape) | Congelado 2008; sin licencia; no es el universo DBTBS |
| `data/dbtbs/gt_tss.tsv` | 3.390 | **BSGatlas v1.1 GFF** (`Sigma`, `Resolution`, `Based on`) | No es DBTBS; el sigma puede venir de Nicolas et al. |
| gold B. subtilis (`SIGMA_GOLD.md`) | 2.226 | subconjunto A/B/D/I/L de `gt_tss.tsv` | **1.777 (79,8 %) solo citan Nicolas et al.**; 379 (17,0 %) mencionan DBTBS; 356 BsubCyc |
| gold K-12 | 3.314 | RegulonDB `PromoterSet.tsv` | 928 S, 2.307 W, 73 C; no todo es evidencia fuerte |

Conteo de procedencia del subconjunto B. subtilis que hoy se usa como gold:

```text
2.226 filas A/B/D/I/L   (categorías solapadas; no suman 100 %)
├── 1.777  "Based on: Nicolas et al." (tiling array)  79,8 %
├──   379  mencionan DBTBS                            17,0 %
├──   356  mencionan BsubCyc                          16,0 %
└──    19  sin "Based on"                              0,9 %
```

Si se exige soporte **DBTBS o BsubCyc** y `Resolution: 0`, el subconjunto
estricto baja a **430 TSS** (A 374, B 36, D 13, L 6, I 1). Es un gold mucho
más defendible que 2.226, aunque más pequeño y limitado a *B. subtilis*. Los
2.226 siguen siendo útiles como **“silver”** para cobertura y para comparar con
trabajos previos, pero no como evidencia curada.

> **No es un bug de código:** es una diferencia de procedencia entre
> `SIGMA_GOLD.md` (etiqueta “DBTBS”) y el archivo que realmente consume
> (`gt_tss.tsv`, BSGatlas). Corregir la etiqueta en la documentación futura,
> no en este documento.

---

## 6. Recomendación de combinación para el catálogo de 400

### 6.1 Principio

El catálogo de 400 no puede tener “gold” para 400 genomas: la literatura
revisada solo ofrece **TSS–sigma experimental para dos especies**
(*E. coli* K-12 y *B. subtilis* 168). Todo lo demás es:

- TSS sin sigma (PPD, MtbRegList);
- sigma sin TSS (PRODORIC);
- TF binding (CollecTF, RegPrecise, CoryneRegNet, SubtiWiki);
- predicción/anotación (CDBProm, MicroScope, P2TF, BSGatlas-Nicolas).

Por tanto, la combinación correcta es **2 anclas + 0 gold adicional + 398
candidatos con abstención**, no “7 especies con gold parcial”.

### 6.2 Asignación por niveles

| Nivel | Fuente | Uso | Filas / alcance |
|---|---|---|---|
| **A1** | RegulonDB v14.5, promotores con `confidenceLevel=S` | calibración y test K-12 | 928 filas (amplio: 3.314) |
| **A2** | DBTBS release 5 crudo (865 TSS); BSGatlas v1.1 solo como mapper de coordenadas | calibración y test B. subtilis | 865 crudos; 430 estrictos (DBTBS/BsubCyc + `Resolution: 0`) |
| **B1** | EcoCyc 2025 | cross-check de consistencia, no independiente | ~E. coli |
| **B2** | PRODORIC | verificar arquitectura −35/−10 y prior de familia | 27 organismos, 3.238 TFBS |
| **B3** | PPD | definir ventana TSS, **nunca** etiqueta sigma | 63 especies |
| **B4** | Pro54DB | validación externa específica σ54 (si las coordenadas se confirman) | 210 promotores, 43 especies |
| **C** | CollecTF, RegPrecise, CoryneRegNet, SubtiWiki, P2TF | features de TF/contexto; controles negativos | variable |
| **D** | CDBProm, MicroScope, P2TF, RegPrecise, BSGatlas-Nicolas | **excluir de gold** | — |

> **No confundir A2 con un solo conjunto.** Los 865 TSS son el universo DBTBS
> crudo (literatura, 2008). Los 430 son el subconjunto del gold BSGatlas con
> soporte DBTBS/BsubCyc y resolución 0, mapeado a la misma assembly. El primero
> es la fuente; el segundo, la intersección de alta confianza.

### 6.3 Reglas de construcción

1. **Unidad de gold:** `(genoma, contig, pos, strand, sigma)` + `source`,
   `source_version`, `evidence_type`, `pmid`, `resolution`, `derived_from`.
2. **Prioridad de fuente:** RegulonDB > DBTBS (crudo) > BSGatlas (solo si
   `Based on` incluye DBTBS/BsubCyc y `Resolution: 0`) > EcoCyc (solo
   cross-check). El subconjunto de 430 no es independencia de proceso: BsubCyc
   también es curaduría de literatura y puede solaparse con DBTBS; es
   independencia de la evidencia débil de tiling array, no de la curación.
3. **Deduplicar** por la clave anterior; conservar todos los sigma de un mismo
   TSS; no colapsar multi-sigma.
4. **Excluir** TSS predichos, propagados o con `Based on` puramente
   Nicolas/CDBProm/RegPrecise.
5. **Partir** por operón/locus y por PMID; nunca el mismo paper en train y
   test. Añadir un split temporal por versión del gold.
6. **No evaluar** con una PWM estimada en el mismo conjunto; estimar PWMs solo
   en train.
7. **No usar** BioCyc/EcoCyc como test independiente de RegulonDB; ni BSGatlas
   como test independiente de DBTBS.
8. **Para los 398 genomas sin gold:** salida `candidate` / `no_model` /
   `sin_sospecha`; reportar cobertura y abstención, no accuracy.
9. **Para las 7 especies “parciales”** (Mtb, PAO1, *S. aureus*, *S.
   coelicolor*, *C. glutamicum*): no hay gold σ–TSS en estas bases. PPD tiene
   TSS sin sigma (PAO1 solo 10 promotores; *S. aureus* 2.821; *C. glutamicum*
   3.581); PRODORIC tiene matrices sin TSS; MtbRegList no tiene campo sigma
  verificable y está caída. Reportarlas como “sin gold”, no como “parcial”.
10. **Licencias:** CC BY-NC permite uso no comercial con atribución; CC
    BY-NC-SA añade *share-alike*. DBTBS, PPD, Pro54DB, RegPrecise, CollecTF y
    CoryneRegNet no tienen licencia de datos explícita. Para redistribuir el
    catálogo, pedir permiso o separar la capa de datos de la capa de
    anotaciones derivadas.

### 6.4 Qué reportar

- Cobertura del inventario sigma (2.314 candidatos del plan) y de modelos PWM.
- Fracción de TSS con evidencia en los 2 anclas.
- Métricas supervisadas **solo** en K-12 y *B. subtilis*: top-1, top-3,
  macro-F1, MCC, precisión/recall por sigma, multi-etiqueta.
- Para los 398: ablación, abstención y estabilidad, no “accuracy global”.
- Un anexo de procedencia por fila (fuente, versión, evidencia, PMID).

---

## 7. Verificación de DOI/PMID/URL y licencias

Corte: **2026-11-26**. Método: PubMed E-utilities para PMID/DOI, Crossref/
OpenAlex para resolución, Zenodo para licencia de datasets, y acceso directo a
los sitios para versión/estado.

| Recurso | DOI | PMID | URL / estado | Licencia verificada |
|---|---|---|---|---|
| RegulonDB v12.0 | `10.1093/nar/gkad1072` | 37971353 | `regulondb.ccg.unam.mx` (200) | Datos v12.0: CC BY-NC-SA 4.0 (Zenodo 10.5281/zenodo.8346680) |
| RegulonDB v14.5 | — | — | `github.com/regulondbunam/regulondb-releases/releases/tag/1.4.5` (200) | Verificar licencia del release |
| DBTBS 2008 | `10.1093/nar/gkm910` | 17962296 | `dbtbs.hgc.jp` (200; descarga con registro) | Sin licencia de datos; artículo CC BY-NC 2.0 UK |
| BSGatlas | `10.1099/mgen.0.000524` | 33539279 | `rth.dk/resources/bsgatlas/` (200, v1.1) | CC BY 4.0 (GFF3 Zenodo 10.5281/zenodo.4305869) |
| SubtiWiki 2025 | `10.1093/nar/gkae957` | 39441067 | `subtiwiki.uni-goettingen.de/v5` (200; SPA/API) | Wiki CC BY; artículo CC BY-NC 4.0; API sin aviso |
| EcoCyc 2025 | `10.1128/ecosalplus.esp-0019-2024` | 40304522 | `ecocyc.org` (200) | BioCyc: licencia + suscripción; EcoCyc data free con licencia |
| EcoCyc 2021 | `10.3389/fmicb.2021.711077` | 34394059 | PMC8357350 | Artículo CC BY 4.0 |
| PRODORIC 2022 | `10.1093/nar/gkab1110` | 34850133 | `prodoric.de` (200) | Datos CC BY-NC 4.0 |
| P2TF 2012 | `10.1186/1471-2164-13-628` | 23153078 | `p2tf.org` (**muerta/parking**, 2026-11-26) | Artículo CC BY 2.0; datos sin licencia |
| PPD 2021 | `10.1016/j.jmb.2021.166860` | 33539888 | `lin-group.cn/database/ppd/` (200, v1.0 2020-09-22) | Gratuita; © LinGroup; sin licencia abierta |
| CDBProm 2024 | `10.1093/nargab/lqae018` | 38385146 | `aw.iimas.unam.mx/cdbprom/` (200) | Artículo CC BY 4.0; descargas sin licencia declarada |
| MicroScope 2020 | `10.1093/nar/gkz926` | 31647104 | `mage.genoscope.cns.fr/microscope` (200; Turnstile) | Artículo CC BY 4.0 |
| MtbRegList 2005 | `10.1093/bioinformatics/bti321` | 15722376 | `USherbrooke.ca/vers/MtbRegList` (**403/WAF**) | Sin licencia |
| CoryneRegNet 7 | `10.1038/s41597-020-0484-9` | 32393779 | `exbio.wzw.tum.de/coryneregnet/` (200, 7.0) | Artículo CC BY 4.0 |
| CollecTF 2014 | `10.1093/nar/gkt1123` | 24234444 | `collectf.org` (shell JS; inestable) | Artículo CC BY 3.0; datos sin licencia |
| CollecTF 2016 | `10.1093/database/baw055` | 27114493 | PMC4843526 | Artículo CC BY 4.0 |
| RegPrecise 3.0 | `10.1186/1471-2164-14-745` | 24175918 | `regprecise.lbl.gov` (200, v3.2) | Artículo CC BY 2.0; datos sin licencia |
| Pro54DB 2017 | `10.1093/bioinformatics/btw630` | 28171531 | `lin-group.cn/database/pro54db/` (200, v1.0 2016-08-31) | Sin licencia explícita |

### 7.1 Erratas y discrepancias detectadas

1. **“Gold DBTBS” = BSGatlas.** El archivo local `gt_tss.tsv` es el GFF de
   BSGatlas v1.1, no la tabla de DBTBS. 1.777 de 2.226 filas del gold solo
   citan Nicolas et al. (tiling array).
2. **CDBProm: 24 M vs 12 M.** El paper dice “over 24 million promoters”; el
   sitio dice “over 12 million”. Verificar la cifra al descargar.
3. **DBTBS: 1.475 promotores ≠ 1.475 TSS.** El scrape local recuperó 865 TSS
   de 1.160 páginas de operón; el número de entradas de promotor no es el
   número de TSS únicos.
4. **EcoCyc no es independiente.** La sincronización con RegulonDB y el equipo
   de curación compartido están declarados en el paper de RegulonDB v12.0.
5. **PPD no tiene sigma.** Verificado en el CSV descargable (columnas `id` a
   `PromoterSeq`, sin campo sigma) y en la búsqueda avanzada (nombre,
   secuencia, TSS, producto, locus tag). No sirve para asignación σ.
6. **PRODORIC no es TSS-level.** Guarda sitios sigma como matrices −35/−10;
   el TSS es la referencia experimental de origen, no una coordenada del
   catálogo.
7. **P2TF está muerta.** El dominio ya no sirve la base (2026-11-26); citar el
   paper, no la URL.
8. **MtbRegList está caída.** El sitio devuelve 403/WAF; no es una fuente
   operativa para el catálogo.

---

## 8. Referencias

- **[R1]** Salgado H, Gama-Castro S, Lara P, *et al.* RegulonDB v12.0: a
  comprehensive resource of transcriptional regulation in *E. coli* K-12.
  *Nucleic Acids Res.* 2024;52(D1):D255–D264. DOI:
  [10.1093/nar/gkad1072](https://doi.org/10.1093/nar/gkad1072). PMID:
  [37971353](https://pubmed.ncbi.nlm.nih.gov/37971353/). PMCID: PMC10767902.
- **[R2]** RegulonDB releases. GitHub, tag 1.4.5 (RegulonDB v14.5/v14.5.1,
  EcoCyc 29.5), 2026-02-26.
  <https://github.com/regulondbunam/regulondb-releases/releases/tag/1.4.5>.
- **[R3]** RegulonDB v12.0.0. Zenodo dataset.
  [10.5281/zenodo.8346680](https://doi.org/10.5281/zenodo.8346680).
- **[R4]** Sierro N, Makita Y, de Hoon MJL, Nakai K. DBTBS: a database of
  transcriptional regulation in *Bacillus subtilis* containing upstream
  intergenic conservation information. *Nucleic Acids Res.*
  2008;36(Database issue):D93–D96. DOI:
  [10.1093/nar/gkm910](https://doi.org/10.1093/nar/gkm910). PMID:
  [17962296](https://pubmed.ncbi.nlm.nih.gov/17962296/). PMCID: PMC2247474.
- **[R5]** DBTBS release 5. <https://dbtbs.hgc.jp/> (flat file actualizado
  2008-10-21; descarga con registro).
- **[R6]** Geissler AS, Anthon C, Alkan F, *et al.* BSGatlas: a unified
  *Bacillus subtilis* genome and transcriptome annotation atlas with enhanced
  information access. *Microb Genom.* 2021;7(2):000524. DOI:
  [10.1099/mgen.0.000524](https://doi.org/10.1099/mgen.0.000524). PMID:
  [33539279](https://pubmed.ncbi.nlm.nih.gov/33539279/). PMCID: PMC8208703.
- **[R7]** BSGatlas GFF3 v1.1. Zenodo.
  [10.5281/zenodo.4305869](https://doi.org/10.5281/zenodo.4305869).
- **[R8]** Elfmann C, *et al.* A new framework for SubtiWiki, the database for
  the model organism *Bacillus subtilis*. *Nucleic Acids Res.*
  2025;53(D1):D864–D870. DOI:
  [10.1093/nar/gkae957](https://doi.org/10.1093/nar/gkae957). PMID:
  [39441067](https://pubmed.ncbi.nlm.nih.gov/39441067/). PMCID: PMC11701700.
- **[R9]** SubtiWiki v5 API (Swagger).
  <https://subtiwiki.uni-goettingen.de/v5/api/swagger.json>.
- **[R10]** Keseler IM, *et al.* The EcoCyc database (2025). *EcoSal Plus.*
  2025. DOI:
  [10.1128/ecosalplus.esp-0019-2024](https://doi.org/10.1128/ecosalplus.esp-0019-2024).
  PMID: [40304522](https://pubmed.ncbi.nlm.nih.gov/40304522/).
- **[R11]** Keseler IM, *et al.* The EcoCyc Database in 2021. *Front
  Microbiol.* 2021;12:711077. DOI:
  [10.3389/fmicb.2021.711077](https://doi.org/10.3389/fmicb.2021.711077).
  PMID: [34394059](https://pubmed.ncbi.nlm.nih.gov/34394059/). PMCID:
  PMC8357350.
- **[R12]** BioCyc downloads and licensing.
  <https://biocyc.org/download.shtml>.
- **[R13]** Dudek C-A, *et al.* PRODORIC: state-of-the-art database of
  prokaryotic gene regulation. *Nucleic Acids Res.* 2022;50(D1):D295–D302.
  DOI: [10.1093/nar/gkab1110](https://doi.org/10.1093/nar/gkab1110). PMID:
  [34850133](https://pubmed.ncbi.nlm.nih.gov/34850133/). PMCID: PMC8728284.
- **[R14]** Ortet P, De Luca G, Whitworth DE, Barakat M. P2TF: a comprehensive
  resource for analysis of prokaryotic transcription factors. *BMC Genomics.*
  2012;13:628. DOI:
  [10.1186/1471-2164-13-628](https://doi.org/10.1186/1471-2164-13-628). PMID:
  [23153078](https://pubmed.ncbi.nlm.nih.gov/23153078/). PMCID: PMC3532121.
- **[R15]** Su W, Liu M-L, Yang Y-H, *et al.* PPD: a manually curated database
  for experimentally verified prokaryotic promoters. *J Mol Biol.*
  2021;433(11):166860. DOI:
  [10.1016/j.jmb.2021.166860](https://doi.org/10.1016/j.jmb.2021.166860).
  PMID: [33539888](https://pubmed.ncbi.nlm.nih.gov/33539888/).
- **[R16]** PPD v1.0. <http://lin-group.cn/database/ppd/> (última
  actualización 2020-09-22).
- **[R17]** Martinez GS, *et al.* CDBProm: the Comprehensive Directory of
  Bacterial Promoters. *NAR Genom Bioinform.* 2024;6(1):lqae018. DOI:
  [10.1093/nargab/lqae018](https://doi.org/10.1093/nargab/lqae018). PMID:
  [38385146](https://pubmed.ncbi.nlm.nih.gov/38385146/). PMCID: PMC10880602.
- **[R18]** Vallenet D, *et al.* MicroScope: an integrated platform for the
  annotation and exploration of microbial gene functions through genomic,
  pangenomic and metabolic comparative analysis. *Nucleic Acids Res.*
  2020;48(D1):D579–D589. DOI:
  [10.1093/nar/gkz926](https://doi.org/10.1093/nar/gkz926). PMID:
  [31647104](https://pubmed.ncbi.nlm.nih.gov/31647104/). PMCID: PMC7145621.
- **[R19]** Jacques P-É, *et al.* MtbRegList, a database dedicated to the
  analysis of transcriptional regulation in *Mycobacterium tuberculosis*.
  *Bioinformatics.* 2005;21(10):2563–2565. DOI:
  [10.1093/bioinformatics/bti321](https://doi.org/10.1093/bioinformatics/bti321).
  PMID: [15722376](https://pubmed.ncbi.nlm.nih.gov/15722376/).
- **[R20]** Parise M, *et al.* CoryneRegNet 7, the reference database and
  analysis platform for corynebacterial gene regulatory networks. *Sci Data.*
  2020;7:142. DOI:
  [10.1038/s41597-020-0484-9](https://doi.org/10.1038/s41597-020-0484-9).
  PMID: [32393779](https://pubmed.ncbi.nlm.nih.gov/32393779/). PMCID:
  PMC7214426.
- **[R21]** Kılıç S, White ER, Sagitova DM, Cornish JP, Erill I. CollecTF: a
  database of experimentally validated transcription factor-binding sites in
  Bacteria. *Nucleic Acids Res.* 2014;42(D1):D156–D160. DOI:
  [10.1093/nar/gkt1123](https://doi.org/10.1093/nar/gkt1123). PMID:
  [24234444](https://pubmed.ncbi.nlm.nih.gov/24234444/). PMCID: PMC3965012.
- **[R22]** Kılıç S, *et al.* From data repositories to submission portals:
  rethinking the role of domain-specific databases in CollecTF. *Database
  (Oxford).* 2016;2016:baw055. DOI:
  [10.1093/database/baw055](https://doi.org/10.1093/database/baw055). PMID:
  [27114493](https://pubmed.ncbi.nlm.nih.gov/27114493/). PMCID: PMC4843526.
- **[R23]** Novichkov PS, *et al.* RegPrecise 3.0 — a resource for genome-scale
  exploration of transcriptional regulation in bacteria. *BMC Genomics.*
  2013;14:745. DOI:
  [10.1186/1471-2164-14-745](https://doi.org/10.1186/1471-2164-14-745). PMID:
  [24175918](https://pubmed.ncbi.nlm.nih.gov/24175918/). PMCID: PMC3840689.
- **[R24]** Liang Z-Y, *et al.* Pro54DB: a database for experimentally
  verified sigma-54 promoters. *Bioinformatics.* 2017;33(3):467. DOI:
  [10.1093/bioinformatics/btw630](https://doi.org/10.1093/bioinformatics/btw630).
  PMID: [28171531](https://pubmed.ncbi.nlm.nih.gov/28171531/).

---

## 9. Enlaces internos

- [`../SIGMA_GOLD.md`](../SIGMA_GOLD.md) — composición del gold actual.
- [`../SIGMA_400_PLAN.md`](../SIGMA_400_PLAN.md) — plan de extensión a 400.
- [`../TSS_LAYER.md`](../TSS_LAYER.md) — procedencia de TSS/TTS locales.
- [`SIGMA_BIOLOGY.md`](SIGMA_BIOLOGY.md) — biología de familias sigma.
- [`SIGMA_PREDICTION.md`](SIGMA_PREDICTION.md) — métodos y evaluación.
- [`GENERALIZATION.md`](GENERALIZATION.md) — generalización entre especies.
