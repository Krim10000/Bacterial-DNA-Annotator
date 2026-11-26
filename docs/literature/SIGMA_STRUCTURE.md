# Estructuras, dominios y bolsillos de factores sigma: qué informan y qué no

> **Documento de revisión bibliográfica y estructural para `dna_poc`.**
> Ruta: `docs/literature/SIGMA_STRUCTURE.md`
> Fecha de corte y verificación: **2026-11-26**
> Autor: subagente `lit-sigma-structure` (revisión; **no modifica código**)
> Estado: **v1.0** — cada DOI fue resuelto en Crossref, cada PMID en NCBI
> E-utilities, cada entrada PDB en la API de RCSB y cada entrada Pfam en la API
> de InterPro. No se incluye ninguna referencia no verificada. Ver Anexo A.
> Alcance: estructura experimental y predicha de σ; **no** es un catálogo de
> promotores ni una revisión de regulación.

---

## Convenciones

- `[Rxx]`: referencia de la tabla de §11; el PMID y el DOI están en esa tabla y
  fueron resueltos uno a uno.
- `[PDB:xxxx]`: entrada del Protein Data Bank consultada por API (título,
  método, resolución, ligandos y cita primaria). La tabla de §2.5 es la fuente
  de las afirmaciones estructurales de este documento.
- `[PFxxxxx]`: entrada Pfam/InterPro consultada por API; el nombre y la
  definición se transcriben de la entrada.
- **Observación interna**: afirmación sobre el repositorio, verificada leyendo
  código o datos; se marca explícitamente y no se presenta como literatura.

Este documento no es una revisión sistemática PRISMA. Es una revisión narrativa
con trazabilidad de fuentes, orientada a una decisión de diseño: **qué evidencia
estructural puede entrar al anotador de sigma y cuál no**. La estrategia de
búsqueda y sus límites están en el Anexo B.

---

## 0. Resumen ejecutivo

1. **La estructura responde preguntas de familia, arquitectura e interfaz; no
   responde preguntas de promotor.** Una σ puede estar perfectamente modelada y
   no decir nada sobre qué secuencia de DNA reconoce ni en qué condición se usa.
2. **Hay dos familias estructuralmente distintas**: la familia σ70 (grupos 1–4,
   con las ECF como grupo 4) y la familia σ54, con arquitectura, mecanismo y
   promotor diferentes. [R01] [R02] [R10] [R11] [R12]
3. **El mapa de dominios es el puente entre estructura y función**: σ1.1
   (compuerta), σ2 (caja –10 y unión al core), σ3 (–10 extendida y core), σ4
   (–35) en σ70; AID/CBD/DBD en σ54. Las ECF conservan principalmente σ2 y σ4
   unidas por un linker no conservado. [R01] [R02] [R13] [R14] [R27] [R28]
4. **La cobertura PDB es asimétrica**: σ70 y σ54 están bien cubiertas en
   complejos con RNAP; las sigma alternativas suelen estar como dominios
   aislados o en complejos inhibidores con su anti-σ. No encontramos una
   estructura de alta resolución de una σ completa, libre y en conformación
   competente para unir DNA.
5. **Los “bolsillos” de una σ son, casi siempre, bolsillos de reconocimiento de
   base o interfaces proteína–proteína**: los aromáticos de σ2 que atrapan A(–11)
   y T(–7), el surco mayor de –35 en σ4, el surco mayor de –24 en σ54, y las
   superficies de contacto con RNAP y anti-σ. [R14] [R15] [R13] [R24]
6. **En las entradas PDB revisadas, los ligandos no poliméricos son iones
   (Mg²⁺, Zn²⁺, K⁺, Cd²⁺) y ADP; el ADP y el Zn²⁺ están unidos a dominios
   anti-σ o a la RNAP, no a un dominio σ como ligando funcional.** La biología
   de ligandos de estos sistemas pertenece a los sensores y anti-σ, no a la σ.
   [R38] [R39] [R42] [R45] [R47] [R48]
7. **Detectar bolsillos (fpocket, P2Rank) predice “ligandabilidad”, no
   función.** En una σ, un bolsillo grande suele coincidir con una interfaz o
   con una cavidad sin significado biológico; nunca con una firma de promotor.
   [R59] [R60]
8. **AlphaFold amplía la cobertura, pero no cambia la lógica**: pLDDT mide
   confianza de plegamiento, no función; las regiones desordenadas (σ1.1,
   linkers) salen con pLDDT bajo; la conformación libre de σ70 es compacta y
   distinta de la conformación unida a RNAP; los ligandos y cofactores no están
   en el modelo. [R34] [R35] [R53] [R55] [R57] [R58]
9. **Lo que sí se puede inferir de estructura**: familia (σ70 vs σ54),
   arquitectura y límites de dominios, integridad (truncamientos reales vs
   truncamientos del modelo), superficies de interfaz con RNAP y con anti-σ, y
   la posición del dominio de unión a DNA. §3.
10. **Lo que no se puede inferir**: el PWM o promotor específico, qué σ se usa
    en una condición, el nombre nativo/ortólogo, la identidad del anti-σ, la
    fuerza del promotor y la regulación por ligandos. §4.
11. **Regla de diseño**: la estructura entra como evidencia **ortogonal de
    familia, dominio, integridad e interfaz**; nunca como evidencia de
    promotor. La especificidad de promotor es un problema de secuencia, grupo
    filogenético y PWM, no de geometría de bolsillo. §9.

---

## 1. Marco: dos familias, cuatro grupos y un mapa de dominios

### 1.1 La familia σ70 y sus grupos

La familia σ70 se define por homología de secuencia y por cuatro regiones
conservadas (1 a 4), con la mayor conservación en las regiones 2 y 4, que unen
RNAP y reconocen –10 y –35 respectivamente. [R04] [R05] [R06] [R10] La familia se
divide en cuatro grupos: 1 (housekeeping), 2 (primarias no esenciales), 3
(alternativas) y 4 (ECF, *extracytoplasmic function*). [R10] [R07] La partición
del espacio de transcripción entre σ housekeeping y alternativas es el marco
clásico para interpretar la competencia por el core. [R03] Las ECF son el grupo
más numeroso y divergente de sigma alternativas: la clasificación de referencia identificó más
de 40 grupos mayores, y la revisión posterior, sobre toda la base NCBI, unos
157 grupos filogenéticos con una media de ~10 ECF por genoma. [R08] [R09]

Que las ECF sean “grupo 4” no significa que sean un dominio menor: son el
sistema de transducción de señal más abundante después de los sistemas de uno y
dos componentes, y su clasificación predice el motivo de promotor **a nivel de
grupo filogenético**, no a partir de la estructura. [R08] [R09]

### 1.2 La familia σ54

σ54 es la otra familia. Sus promotores tienen posiciones conservadas alrededor
de –24 y –12 en lugar de –35 y –10. [R11] A diferencia de σ70, la holoenzima
σ54–RNAP forma un complejo transcripcionalmente silencioso y requiere
activadores dependientes de ATP (bEBP, *bacterial enhancer binding proteins*),
que se ensamblan como anillos hexaméricos y contactan σ54 a través del motivo
GAFTGA. [R12] [R33] La estructura del dominio C-terminal de σ54 unido al
elemento –24 muestra que el motivo conservado RpoN box se inserta en el surco
mayor del DNA. [R24] La comparación de las estructuras de RNAP–σ54 con las de
RNAP–σ70 muestra dominios funcionalmente análogos pero una disposición e
interacciones distintas con RNAP, lo que explica el mecanismo diferente. [R21]

**Consecuencia operativa:** σ54 no se puede anotar con el panel de dominios ni
con la geometría de promotor de σ70. El escaneo de –35/–10 no aplica; el de
–24/–12 sí. Esta distinción ya está incorporada como advertencia en el módulo
`modules/sigma/module.py` (observación interna).

### 1.3 Mapa de dominios y qué se puede leer en cada uno

**Tabla 1. Regiones/dominios de σ y su lectura estructural**

| Región / dominio | Pfam | Función canónica | Contacto estructural | Qué se puede leer en una estructura | Qué NO se puede concluir |
|---|---|---|---|---|---|
| σ1.1 (N-terminal) | sin Pfam dedicado | Compuerta; en la holoenzima bloquea el canal de DNA; autoinhibe la σ libre | Canal de DNA de RNAP; contactos internos en la conformación compacta | Presencia/posición; desorden o compactación; no es un dominio de unión a DNA | Que un pLDDT bajo implique pérdida de función |
| σ1.2 | PF00140 | Marcador de grupos 1–2; función estructural | — | Presencia/ausencia del dominio | Que su ausencia cambie el promotor |
| σ2.1–2.4 | PF04542 | σ2.1–2.2: unión al core (β′ *coiled-coil*); σ2.3–2.4: reconocimiento de –10 y *base flipping* | β′ clamp; ssDNA de –10 | Los bolsillos que atrapan A(–11) y T(–7); la interfaz con RNAP | El PWM completo; la identidad del promotor |
| σ3.0–3.2 | PF04539 | σ3.0: –10 extendida; σ3.1–3.2: unión al core | RNAP; DNA | Presencia del módulo de –10 extendida; interfaz con core | Qué promotores usan –10 extendida *in vivo* |
| σ4.1–4.2 | PF04545 | HTH que reconoce –35; superficie de unión a β-flap/activadores | Surco mayor de –35; β-flap | El surco de lectura de –35 y la geometría HTH | El PWM de –35; la fuerza del promotor |
| Linker σ2/σ4 | — | En ECF es no conservado y se inserta en el canal activo de RNAP | Canal activo de RNAP | Posición del linker; no su función | Función reguladora a partir de la longitud |
| σ54 AID | PF00309 | Interacción con el activador (bEBP); inhibe la iniciación antes del activador | bEBP (GAFTGA) | La interfaz σ54–activador | El promotor –24/–12 |
| σ54 CBD | PF04963 | Unión al core de RNAP | β/β′ | La interfaz con el core; el “punto de fractura” conformacional | El promotor |
| σ54 DBD | PF04552 | Reconocimiento de –24 (RpoN box) | Surco mayor de –24 | El surco de lectura de –24 | El elemento –12 completo |
| Anti-σ N-terminal | PF12791 | Dominio de unión a σ en anti-σ tipo RsgI | σ | Contexto de familia de anti-σ | Que sea universal a todos los anti-σ |

Fuentes: definiciones Pfam/InterPro verificadas (PF00140, PF04542, PF04539,
PF04545, PF00309, PF04963, PF04552, PF12791) [R62]; arquitectura y función
[R01] [R02] [R13] [R14]; ECF y linker [R27] [R28]; σ54 [R12] [R21] [R24]
[R25].

**Nota sobre Pfam:** la entrada PF04545 (región 4) declara explícitamente que su
umbral se fijó alto para evitar solapamientos con otras familias HTH y que, por
tanto, **tiene falsos negativos**. [PF04545] Un “no hit” de PF04545 no es
evidencia de que la proteína no tenga dominio de unión a –35; es una limitación
del modelo. Esta advertencia es relevante para cualquier anotador que use el
panel HMM como verdad de familia.

### 1.4 Comparación de familias

**Tabla 2. σ70 (grupos 1–3), ECF (grupo 4) y σ54**

| Propiedad | σ70 grupos 1–3 | ECF (grupo 4) | σ54 |
|---|---|---|---|
| Elementos de promotor | –35/–10 (con –10 extendida, UP, discriminador) | –35/–10 específicos de grupo; el espaciador y un T-tract proximal a –35 modulan la selectividad | –24/–12 |
| Arquitectura | σ1.1–σ4, cuatro regiones | σ2 y σ4 como módulos principales, con región 3 no conservada y linker no conservado | AID–CBD–DBD |
| Unión al core | σ2.2, σ3.1, σ4 | σ2/σ4 y linker en el canal activo | CBD |
| Activador | Normalmente no requiere ATPasa | Normalmente no requiere ATPasa | bEBP AAA+ dependiente de ATP |
| Anti-σ | Frecuente (Rsd, AsiA, FlgM, RseA…) | Frecuente y a menudo transmembrana (RseA, RsiW, MucA, ChrR…) | Menos caracterizado |
| Cobertura estructural | Alta (holoenzima, open complex, dominios) | Media (complejos de iniciación y anti-σ) | Media-alta (dominios y complejos) |
| Qué informa la estructura | Familia, interfaz, mecanismo | Grupo estructural; interfaz con anti-σ | Familia, interfaz con activador/core |

Fuentes: [R01] [R02] [R08] [R09] [R11] [R12] [R21] [R27] [R28] [R32].

---

## 2. Qué estructuras existen y qué muestra cada una

### 2.1 σ70: dominios sueltos, holoenzima y complejos de iniciación

El repertorio estructural de σ70 es el más completo. La estructura de los
dominios σ2, σ3 y σ4 de σA de *Thermus aquaticus*, y del dominio σ4 unido al
elemento –35, estableció que σ4 usa un HTH para leer –35 y que las superficies
expuestas de cada dominio son importantes para unir RNAP. [R13] La estructura
de σ2 unido a ssDNA de –10 mostró el mecanismo de *base flipping*: A(–11) y
T(–7) se extruyen de la pila y se entierran en bolsillos del dominio σ2. [R14]
Las estructuras de complejos de iniciación de *T. thermophilus* y *E. coli*
añadieron el discriminador, el *unwinding* y la organización del DNA en el
complejo abierto. [R15] [R18] [R19] [R20]

La holoenzima de *T. thermophilus* a 2.6 Å [R16], la holoenzima de *E. coli* a
3.7 Å [R17] y las estructuras cryo-EM del complejo abierto [R20] muestran la
disposición de σ sobre el core y, en el caso de *E. coli*, la posición de σ1.1
en el canal de DNA, donde bloquea la entrada al sitio activo. [R17]

**Lectura:** la estructura de σ70 permite mapear con precisión las superficies
de unión a RNAP, la posición de σ1.1 y los determinantes de –10/–35. No permite
predecir el promotor completo de una σ70 concreta, porque el PWM depende de la
identidad de los residuos que contactan las bases y de la geometría del
espaciador, no solo del pliegue.

### 2.2 σ54: una familia con otra geometría

Las estructuras de σ54 cubren el dominio C-terminal unido a –24 (NMR) [R24],
el dominio de unión al core (CBD) [R25], un fragmento de σN unido a promotor
(3.4 Å) [R26], la holoenzima RNAP–σ54 (3.8 Å) [R21] y los complejos cerrado,
intermedio y abierto por cryo-EM. [R22] [R23] El conjunto muestra que σ54 usa
un dominio HTH para el surco mayor de –24, que su disposición sobre RNAP difiere
de la de σ70 y que la apertura del DNA es un proceso activo que requiere el
activador. [R21] [R22] [R23] [R24]

**Lectura:** la estructura confirma familia y mecanismo (activador, apertura,
interfaz), y localiza el dominio de unión a –24. No convierte a σ54 en una σ70
“con otro promotor”: el escaneo y el modelo de promotor deben ser distintos.

### 2.3 ECF y sigma alternativas: dominios, complejos de iniciación y anti-σ

Las ECF son estructuralmente σ70-like en sus dominios σ2 y σ4, pero con
diferencias funcionales importantes. La estructura σE–RseA mostró que los dos
dominios mayores de σE son casi idénticos a los de otras σ70 y que RseA inhibe
ocluciendo los dos determinantes primarios de unión al core. [R36] Las
estructuras de iniciación con σE de *E. coli* y σH de *M. tuberculosis*
mostraron que las ECF usan un linker σ2/σ4 no conservado que se inserta en el
canal activo y sale por el canal de salida de RNA; el mecanismo de
reconocimiento y apertura es distinto del de las σ primarias. [R27] [R28] La
estructura del dominio σ4 de SigW unido a –35 y la comparación de SigR R4 con
σE muestran el modo de lectura de –35 en ECF y que un solo residuo (Met188 en
SigR frente a Arg171 en σE) puede cambiar la preferencia por una base. [R30] [R31] El espaciador y
un T-tract proximal a –35 también modulan la selectividad entre σM, σW, σX y σV.
[R32]

Las sigma alternativas de la familia σ70 no ECF están peor cubiertas. σ28
(flagelar) tiene estructura completa en complejo con FlgM, pero en una
conformación inactiva. [R37] σF (esporulación) está en complejo con SpoIIAB,
con ~80 % de la σ desordenada. [R38] Las σI de *Clostridium thermocellum* se
resolvieron en complejos abiertos y muestran un modo de reconocimiento distinto
del de otras σ70, con diferencias entre parálogos. [R29]

**Lectura:** las ECF y las alternativas se pueden reconocer por arquitectura y
por familia, y sus interfaces con anti-σ son estructuralmente visibles. Pero la
cobertura es desigual y muchas veces el modelo es un dominio o un complejo
inhibidor, no la σ libre competente para DNA.

### 2.4 Anti-σ y reguladores: la estructura informa regulación, no promotor

Los complejos σ–anti-σ revelan mecanismos distintos: RseA ocluye los
determinantes de unión al core de σE [R36]; FlgM entierra la hélice de
reconocimiento de –35 de σ28 y empaqueta los dominios [R37]; SpoIIAB, una
quinasa de la superfamilia GHKL, une ADP/ATP y deja desordenada la mayor parte
de σF [R38] [R39]; AsiA remodela σ4 y secuestra su superficie de DNA [R40];
Rsd se une a σ4 y bloquea tanto la unión al core como la unión a –35 [R41].

En ECF, ChrR aporta un módulo de dos dominios: un dominio anti-σ N-terminal que
coordina Zn²⁺ y un dominio C-terminal tipo cupin que coordina otro Zn²⁺ y
responde al oxígeno singlete. [R42] CnrY usa una pequeña “perilla” hidrofóbica
como determinante de la interfaz con CnrH. [R43] SigK–RskA muestra un puente
disulfuro en la región de σK que interactúa con –35, con un papel redox en la
estabilidad del complejo. [R44] RsiW, un anti-σ de la familia ZAS con motivo
CHCC, une Zn²⁺ y se une a un surco de SigW; la oxidación del motivo no disocia
el complejo. [R45] MucA se une a AlgU y MucB lo reconoce en el periplasma, con
un disulfuro en MucB que afecta su oligomerización. [R46]

RsrA es el ejemplo canónico de anti-σ con Zn²⁺: el Zn²⁺ no es necesario para
unir σR pero sí para la actividad anti-σ funcional y para modular la reactividad
de los tioles; la oxidación libera Zn²⁺, colapsa RsrA y libera σR. [R47] [R48]
En el lado del sensor, RseB tiene un dominio que se parece al *barrel* de
proteínas que unen el ancla lipídica de lipoproteínas, y DegS se activa por
unión de péptidos a su dominio PDZ. [R49] [R51] La proteólisis de RseA y de
FecR por RseP conecta la señal de envoltura con la liberación de la σ. [R50]
[R52]

**Lectura:** la estructura de un anti-σ o de un sensor informa **cómo se
regula** la σ (secuestro, proteólisis, redox, Zn²⁺, nucleótido), no **qué
promotor** reconoce. Es evidencia de regulación, no de especificidad.

### 2.5 Tabla estructura / familia / evidencia

**Tabla 3. Entradas PDB verificadas y qué evidencia aporta cada una**

| PDB | Molécula / complejo | Familia / grupo | Método / Å | DNA / ligando | Evidencia que aporta | Qué NO muestra |
|---|---|---|---|---|---|---|
| 1KU7 | σA R4 de *T. aquaticus* + DNA –35 | σ70 G1 | X-ray 2.4 | DNA | HTH de σ4 en el surco mayor de –35 | PWM completo; σ usada en una condición |
| 3UGO | σ70 R2 de *E. coli* + ssDNA –10 | σ70 G1 | X-ray 2.1 | ssDNA; K⁺ | Bolsillos de A(–11) y T(–7); *base flipping* | –35; promotor completo |
| 4YG2 | RNAP–σ70 de *E. coli* | σ70 G1 | X-ray 3.7 | Mg²⁺, Zn²⁺ | σ1.1 en el canal de DNA; disposición de dominios sobre el core | DNA; promotor |
| 4YLN | Complejo de iniciación de *E. coli* (burbuja completa, espaciador 17 bp) | σ70 G1 | X-ray 5.5 | DNA/RNA; Mg²⁺, Zn²⁺ | Contacto σ–RNA naciente; comparación con espaciador 16 bp (4YLO) | Detalle atómico; secuencia específica |
| 6CA0 | Complejo abierto σ70 de *E. coli* | σ70 G1 | cryo-EM 5.75 | DNA; Mg²⁺, Zn²⁺ | Arquitectura del complejo abierto | Cadenas laterales; PWM |
| 1IW7 | Holoenzima de *T. thermophilus* | σ70 G1 | X-ray 2.6 | Mg²⁺, Pb²⁺ | Interfaz σ–core a alta resolución | DNA; promotor |
| 4XLN | Iniciación de *T. aquaticus* | σ70 G1 | X-ray 4.0 | DNA/RNA; Mg²⁺, Zn²⁺ | Burbuja y contactos σ–core | Especificidad de secuencia |
| 2O8K | Dominio RpoN de σ54 + DNA –24 | σ54 | NMR | DNA | RpoN box en el surco mayor de –24 | –12; holoenzima |
| 2K9L | Dominio CBD de σ54 | σ54 | NMR | — | Pliegue de unión al core; punto de fractura conformacional | DNA; promotor |
| 5UI5 | σN de *A. aeolicus* + promotor | σ54 | X-ray 3.4 | DNA | Contactos σN–promotor (–24) | Activación por bEBP; interfaz en la holoenzima |
| 6GH5 | Complejo abierto σ54–RNAP | σ54 | cryo-EM 3.4 | DNA | Carga y apertura del DNA; estabilización de la burbuja | Otras σ; promotor |
| 1OR7 | σE de *E. coli* + RseA (dominio citoplásmico) | ECF | X-ray 2.0 | — | RseA ocluye los determinantes de unión al core | DNA; promotor |
| 5FGM | SigR R4 de *S. coelicolor* | ECF | X-ray 2.6 | — | Pliegue de R4 de ECF; base de la comparación con σE | DNA; PWM |
| 6JHE | SigW R4 de *B. subtilis* + DNA –35 | ECF | X-ray 3.1 | DNA | Lectura de –35 por R4 de ECF; exclusión mutua con RsiW | PWM; comportamiento *in vivo* |
| 5OR5 | σE R2 modificada + ssDNA (GTAAAA) | ECF | NMR | ssDNA | Contactos R2–DNA en una σE modificada (estructura no publicada) | Especificidad fisiológica |
| 5ZX2 | RNAP–σH de *M. tuberculosis* + promotor | ECF | X-ray 2.8 | DNA/RNA; Zn²⁺ | Iniciación por una ECF; interfaz σH–core | Otras ECF; promotor |
| 6JCY | RNAP–σH/E quimérica de *M. tuberculosis* | ECF | X-ray 3.1 | DNA/RNA; Mg²⁺, Zn²⁺ | Linker σ2/σ4 en el canal activo | σH/E nativa; promotor |
| 8I23 | Complejo abierto σI1 de *C. thermocellum* | σ70 (σI) | cryo-EM 3.0 | DNA | Modo de reconocimiento distinto de σI; diferencias entre parálogos | Otras σI; PWM |
| 1RP3 | σ28 + FlgM | σ70 (alternativa) | X-ray 2.3 | — | FlgM entierra la hélice de –35; σ intacta pero inactiva | DNA; promotor |
| 1L0O | SpoIIAB + σF | σ70 (alternativa) | X-ray 2.9 | ADP, Mg²⁺ | Anti-σ quinasa GHKL; σF mayormente desordenada | DNA; promotor |
| 1TLH | AsiA + σ70 R4 | σ70 G1 | NMR | — | AsiA remodela σ4 y secuestra su superficie de DNA | DNA; promotor |
| 2P7V | Rsd + σ70 R4 | σ70 G1 | X-ray 2.6 | Mg²⁺ | Rsd bloquea unión al core y a –35 | DNA; promotor |
| 2Q1Z | SigE + ChrR | ECF | X-ray 2.4 | Zn²⁺ | Módulo anti-σ con Zn²⁺; base de la familia | DNA; promotor |
| 4CXF | CnrH + CnrY | ECF | X-ray 1.75 | — | Interfaz mínima anti-σ (perilla hidrofóbica) | DNA; promotor |
| 4NQW | SigK + RskA | ECF | X-ray 2.4 | Cd²⁺ | Disulfuro redox en la región –35 de σK | DNA; promotor |
| 5WUQ | SigW + RsiW | ECF | X-ray 2.8 | Zn²⁺ | Surco de unión de un anti-σ ZAS; Zn²⁺ estructural | DNA; promotor |
| 6IN7 | AlgU + MucA | ECF | X-ray 1.96 | — | Interfaz σ–anti-σ conservada; ausencia de disulfuro en AlgU | DNA; promotor |
| 5FRF | RsrA (Zn²⁺, reducida) | anti-σ (ZAS) | NMR | Zn²⁺ | Sitio de Zn²⁺ y núcleo hidrofóbico del sensor | Complejo con σR; promotor |
| 1SOT / 1SOZ | DegS ± péptido activador | regulador | X-ray 2.3 / 2.4 | péptido | Activación por unión de péptido al PDZ | σ; promotor |

**Cómo leer esta tabla.** La columna “evidencia” es deliberadamente
conservadora: describe lo que la estructura *muestra*, no lo que se infiere. Por
ejemplo, 1KU7 muestra cómo σ4 lee –35, pero no permite escribir el PWM de una σ
concreta; 2Q1Z muestra un Zn²⁺ en ChrR, pero no autoriza a decir que la σ tiene
un ligando.

---

## 3. Qué se puede inferir de una estructura

### 3.1 Familia: σ70 vs σ54, y ECF dentro de σ70

La estructura permite distinguir los dos pliegues y las dos arquitecturas:
σ70 tiene dominios σ2/σ3/σ4 sobre un core y una σ1.1 N-terminal; σ54 tiene
AID/CBD/DBD y una disposición distinta sobre RNAP. [R21] [R24] [R25] Las ECF se
reconocen dentro de σ70 por la presencia de σ2 y σ4 con una región 3 no conservada
y un linker no conservado. [R27] [R28]

**Pero la asignación de familia no debe hacerse solo por estructura.** Un HTH
de σ4 es un pliegue muy común; la similitud estructural de un dominio aislado
con σ4 no convierte a una proteína en σ. La asignación robusta combina:
(a) perfiles HMM de la familia, (b) arquitectura completa (orden y número de
dominios), (c) contexto génico y (d) estructura como confirmación. [R62] [R63]

### 3.2 Arquitectura y límites de dominios

La estructura (experimental o AlphaFold) permite localizar los límites de los
dominios y comprobar si la arquitectura es coherente con la familia. Es útil
para:

- detectar dominios ausentes o duplicados;
- distinguir un dominio mal modelado de un dominio realmente ausente;
- localizar el dominio de unión a DNA (σ4 en σ70; DBD en σ54);
- ver si una ECF conserva σ2 y σ4 o si es una variante con arquitectura
  reducida.

La estructura **no** sustituye al HMM: un modelo puede estar incompleto aunque
el gen esté completo, y un HMM puede fallar aunque la estructura sea clara. La
combinación de ambos es más informativa que cualquiera por separado.

### 3.3 Integridad: truncamientos, pseudogenes y errores de anotación

Un modelo estructural puede sugerir que una proteína está truncada, pero hay
que separar tres casos:

1. **Gen truncado real** (codón de parada interno, frameshift, deleción): el
   modelo y la secuencia coinciden en el truncamiento.
2. **Modelo incompleto** (región no modelada, desorden, límite del dominio):
   la secuencia codifica una proteína completa, pero el modelo no la cubre.
3. **Error de anotación** (inicio de gen mal asignado, fusión/fisión): la
   secuencia traducida no es la proteína real. La misanotación funcional en
   bases públicas es un problema documentado, con tasas del 5–63 % en
   superfamilias enzimáticas. [R66]

**Regla:** un truncamiento estructural solo se convierte en “pseudogen” cuando
lo confirma la secuencia nucleotídica (marco, parada, promotor, evidencias de
transcripción). La estructura sola no alcanza.

### 3.4 Interfaz con RNAP

La estructura localiza las superficies de unión al core:

- σ2.1–2.2 (β′ *coiled-coil*/clamp);
- σ3.1–3.2 (core);
- σ4 (β-flap y superficies de activadores);
- σ1.1 (canal de DNA; en la holoenzima de *E. coli* bloquea la entrada);
- σ54 AID/CBD (interfaz con bEBP y core).

Fuentes: [R13] [R16] [R17] [R21] [R27] [R28]. La estructura también muestra
que las ECF usan un linker σ2/σ4 que se inserta en el canal activo, un detalle
que no se deduce de la secuencia. [R27]

**Qué se puede afirmar:** “esta variante conserva los residuos/superficies de
unión al core” o “esta variante tiene una deleción en la interfaz”. **Qué no se
puede afirmar:** “esta σ es la que se usa en esta condición”.

### 3.5 Interfaz con anti-σ

La estructura permite clasificar el mecanismo de inhibición, que es información
de regulación:

| Mecanismo | Ejemplo | Qué muestra la estructura | Fuente |
|---|---|---|---|
| Oclusión de determinantes de unión al core | RseA–σE | RseA tapa las dos superficies primarias de unión a RNAP | [R36] |
| Enterramiento de la hélice de –35 | FlgM–σ28 | FlgM estabiliza una conformación compacta e inactiva | [R37] |
| Oclusión de la superficie de RNAP + desorden | SpoIIAB–σF | GHKL; ADP; 80 % de σF desordenada | [R38] [R39] |
| Remodelación de σ4 | AsiA–σ70 R4 | σ4 pierde la hélice C-terminal y cambia de conformación | [R40] |
| Bloqueo de σ4 (core y –35) | Rsd–σ70 R4 | Rsd interfiere con las dos funciones de σ4 | [R41] |
| Zn²⁺ estructural/redox | ChrR, RsrA, RsiW | Sitios de Zn²⁺ y cambios redox | [R42] [R45] [R47] [R48] |
| Disulfuro en la región –35 | SigK–RskA | Puente disulfuro redox en σK | [R44] |
| Perilla hidrofóbica | CnrY–CnrH | Interfaz mínima de una clase de anti-σ | [R43] |

**Qué se puede afirmar:** “la superficie X es una interfaz plausible con un
anti-σ” o “la variante conserva el sitio de Zn²⁺”. **Qué no se puede afirmar:**
“el anti-σ es RsiW” o “esta σ está activa en esta condición”.

### 3.6 Conformación libre vs unida: el modelo no es la verdad conformacional

La σ libre no es simplemente la σ unida a RNAP sin RNAP. En σ70, σ1.1
autoinhibe la unión a DNA estabilizando una conformación compacta incompatible
con DNA; [R34] estudios de *cross-linking* masivo en células confirman una
conformación compacta de σ70 apo, distinta de la extendida unida a RNAP. [R35]
En la holoenzima, σ1.1 está en el canal de DNA. [R17] La estructura de σ28 con
FlgM muestra otra conformación compacta e inactiva. [R37]

**Consecuencia:** un modelo AlphaFold de una σ libre puede parecerse a la
conformación compacta o a la extendida, y ninguna de las dos es “la” estructura
funcional en todos los contextos. La posición relativa de dominios en un modelo
apo no debe usarse para inferir interfaces in vivo sin evidencia adicional.

---

## 4. Qué NO se puede inferir

### 4.1 El promotor específico o el PWM

La estructura muestra *cómo* se lee una base (bolsillo, contacto, geometría),
pero no *qué* secuencia se lee. El PWM requiere:

- la identidad de los residuos que contactan bases (no solo el pliegue);
- la geometría del espaciador y de los elementos vecinos;
- datos de unión o de expresión para calibrar;
- en muchos casos, comparación con parálogos y mutantes.

Ejemplos que lo ilustran: un solo residuo de SigR R4 (Met188 vs Arg171 de σE)
cambia la preferencia por la base –31′; [R31] la comparación de dos σI parálogas
revela determinantes de especificidad que no se ven en una sola estructura;
[R29] el T-tract del espaciador modula la selectividad ECF. [R32] La estructura
es necesaria para interpretar estos datos, pero no los sustituye.

### 4.2 Qué σ se usa en una condición

La estructura no dice qué σ está activa, qué anti-σ la secuestra, qué proteasa
la libera o qué factor global modula su disponibilidad. Eso
es biología de regulación: anti-σ, proteólisis, redox, disponibilidad de core,
competencia entre σ. [R01] [R36] [R37] [R45] [R48] La estructura de un complejo
σ–anti-σ puede sugerir un mecanismo, pero no la condición fisiológica.

### 4.3 El nombre nativo y la ortología

La estructura no distingue “SigW” de “SigX” si el pliegue es el mismo. Los
nombres nativos no coinciden necesariamente con los grupos filogenéticos ni con
la ortología. [R08] [R09] La asignación de nombre requiere
ortología, contexto génico y, cuando exista, grupo filogenético; la estructura
solo puede confirmar familia/arquitectura. [R08] [R09]

### 4.4 La identidad del anti-σ y la regulación por ligandos

La estructura de una σ sola no dice qué anti-σ la regula. La interfaz se puede
predecir por co-plegamiento (AF3) o por homología, pero es una hipótesis. La
estructura de un anti-σ sí puede revelar un sitio de Zn²⁺, un bolsillo de
nucleótido o una interfaz con un sensor; eso informa regulación, no promotor.
[R38] [R42] [R47] [R49]

### 4.5 La fuerza del promotor y la transferencia entre especies

La estructura no estima la afinidad de un promotor concreto ni la fuerza de la
iniciación. Tampoco garantiza transferencia entre especies: un pliegue
conservado puede convivir con motivos de promotor distintos. La clasificación
filogenética de ECF predice motivos por grupo, pero esa predicción es de
secuencia, no de estructura. [R08] [R09]

### 4.6 Tabla de afirmaciones y evidencia necesaria

**Tabla 4. Qué afirmación autoriza cada capa de evidencia**

| Afirmación | ¿Basta la estructura sola? | Evidencia necesaria |
|---|---|---|
| “Pertenece a la familia σ70” | No | HMM + arquitectura + estructura de confirmación |
| “Es una ECF del grupo X” | No | Filogenia/grupo + arquitectura + contexto |
| “Es σ54” | Casi | Pliegue/arquitectura + HMM de σ54 + –24/–12 |
| “El dominio σ4 está completo” | Sí, con reservas | Modelo + secuencia + HMM |
| “Es un pseudogen” | No | Secuencia nucleotídica + marco + evidencia de transcripción |
| “Une RNAP” | No | Interfaz estructural + homología + experimento |
| “La regula el anti-σ X” | No | Co-plegamiento + genómica + experimento |
| “Reconoce el promotor P” | No | PWM/experimento; la estructura solo ayuda a interpretar |
| “Se usa en la condición C” | No | Expresión/regulón/regulación |
| “Tiene un ligando L” | No | Estructura del complejo + bioquímica; no docking |

---

## 5. Bolsillos y ligandos: la pregunta central

### 5.1 Qué es un “bolsillo” en una σ

En las estructuras de σ, los bolsillos y cavidades que aparecen se pueden
clasificar en cuatro tipos:

1. **Bolsillos de reconocimiento de base.** En σ2, A(–11) y T(–7) se extruyen
   de la pila y se entierran en bolsillos del dominio; [R14] en σ4, el surco
   mayor de –35; [R13] [R30] en σ54, el surco mayor de –24. [R24] Son bolsillos
   de DNA, no de ligando.
2. **Interfaces proteína–proteína.** Superficies de unión al core (σ2.2, σ3.1,
   σ4) y a anti-σ (RseA, FlgM, Rsd, AsiA, RsiW, MucA…). [R13] [R36] [R37]
   [R41] [R45] [R46]
3. **Sitios de metal.** Zn²⁺ en anti-σ (ChrR, RsrA, RsiW), Cd²⁺ en RskA,
   Mg²⁺ en la RNAP. [R38] [R42] [R44] [R45] [R47]
4. **Sitios de nucleótido.** ADP/ATP en SpoIIAB, una quinasa anti-σ. [R38]
   [R39]

**Ninguno de los cuatro tipos codifica por sí mismo la especificidad de
promotor.** El tipo 1 informa el mecanismo de lectura; los tipos 2–4 informan
unión y regulación.

### 5.2 Census de ligandos en las entradas revisadas

En la Tabla 3, los ligandos no poliméricos observados son:

- **Mg²⁺**: 4YG2, 4YLN, 6CA0, 4XLN, 1IW7, 2P7V, 1L0O;
- **Zn²⁺**: 4YG2, 4YLN, 6CA0, 4XLN, 5ZX2, 6JCY, 2Q1Z, 5WUQ, 5FRF;
- **K⁺**: 3UGO;
- **ADP**: 1L0O (SpoIIAB);
- **Cd²⁺**: 4NQW (RskA);
- **Pb²⁺**: 1IW7 (cristalización).

No encontramos, en este conjunto, una entrada con un compuesto *drug-like*
unido a un dominio σ. El ADP está unido a SpoIIAB y el Zn²⁺ a dominios anti-σ o
a la RNAP. Esto es un resultado acotado al conjunto revisado, no una prueba de
que no exista ningún ligando de σ en la literatura; lo que sí se puede afirmar
es que **la biología de ligandos de estos sistemas está en los anti-σ y
sensores, no en el dominio σ**.

### 5.3 ¿Informan función sigma o solo regulación/unión?

**Informan mecanismo y unión; no informan especificidad de promotor.**

- Un bolsillo de σ2 que atrapa A(–11) dice que esa σ usa *base flipping* en
  –10; no dice el resto del PWM.
- Un surco de σ4 que lee –35 dice que hay un HTH funcional; no dice qué bases
  prefiere, salvo que se conozcan los residuos que las contactan.
- Un sitio de Zn²⁺ en un anti-σ dice que hay un sensor redox; no dice qué
  promotor se activa.
- Un sitio de ADP en SpoIIAB dice que hay una quinasa; no dice qué σ se libera
  en qué condición.

La razón es simple: la especificidad de promotor es una propiedad de la
**secuencia de la superficie de contacto** y de la geometría del DNA, no de la
existencia de una cavidad. Dos σ con el mismo pliegue y bolsillos equivalentes
pueden tener PWM distintos; dos parálogas con estructuras casi idénticas pueden
diferir en especificidad. [R29] [R31]

### 5.4 Dónde sí hay ligandos biológicos: anti-σ y sensores

**Tabla 5. Ligandos/cofactores en sistemas σ–anti-σ y su lectura**

| Sistema | Ligando/cofactor | Dónde está | Qué informa | Qué no informa |
|---|---|---|---|---|
| ChrR–σE (*R. sphaeroides*) | Zn²⁺ (ASD y cupin) | Anti-σ | Sensor redox; módulo de dos dominios | Promotor de σE |
| RsrA–σR (*S. coelicolor*) | Zn²⁺ | Anti-σ | Sensor de disulfuro; liberación de σR | Promotor de σR |
| RsiW–SigW (*B. subtilis*) | Zn²⁺ (CHCC) | Anti-σ (ZAS) | Zn²⁺ estructural; la oxidación no disocia | Promotor de SigW |
| SpoIIAB–σF (*Bacillus*) | ADP/ATP, Mg²⁺ | Anti-σ quinasa | Ciclo de fosforilación; afinidad | Promotor de σF |
| DegS–RseA | Péptido de OMP | Sensor PDZ | Activación proteolítica | Promotor de σE |
| RseB–RseA | Lípidos/lipoproteínas (bolsillo tipo *barrel*, sugerido por la estructura) | Sensor periplásmico | Señal de envoltura | Promotor de σE |
| FecR–σFecI | Citrato férrico (sistema Fec) | Regulador/sensor | Señal de hierro; proteólisis por RseP | Promotor de σFecI |

Fuentes: [R42] [R45] [R47] [R48] [R51] [R49] [R52]. En FecR, la evidencia
citada es la regulación dependiente de citrato férrico y la proteólisis por
RseP, no una estructura del complejo con el ligando; se marca como evidencia de
regulación, no estructural.

### 5.5 Regla para el anotador

1. Un bolsillo que solapa la superficie de unión a DNA (σ2/σ4/DBD) se etiqueta
   como **“superficie de unión a DNA”**, nunca como “promotor”.
2. Un bolsillo que solapa una interfaz conocida con RNAP o anti-σ se etiqueta
   como **“hipótesis de interfaz”**, con confianza baja.
3. Un bolsillo que no solapa nada conocido se etiqueta como **“cavidad sin
   interpretación”** y no se usa.
4. Un ligando detectado por homología o *docking* no se anota como regulación
   sin estructura del complejo o bioquímica.
5. El *score* de P2Rank o el volumen de fpocket **no** son features de función
   ni de familia.

---

## 6. AlphaFold y PDB: qué aportan y cuáles son sus límites

### 6.1 Cobertura

AlphaFold DB contiene más de 214 millones de estructuras predichas y cubre
proteomas completos, incluidas bacterias; [R55] la primera versión ya incluía
21 proteomas modelo. [R54] El PDB experimental contiene del orden de 10⁵ estructuras
de proteínas únicas, muchas menos que secuencias, y está sesgado hacia proteínas
cristalizables, organismos modelo y fragmentos. [R53] [R61] Para σ, esto
significa que la mayoría de las especies tendrá un modelo AlphaFold, pero solo
unas pocas tendrán estructuras experimentales de su σ.

### 6.2 pLDDT y PAE: qué miden

- **pLDDT** es una confianza por residuo en la estructura local; no es una
  medida de función ni de conservación. [R53]
- **PAE** (*predicted aligned error*) mide la confianza en la posición relativa
  de dos residuos; es útil para juzgar la orientación entre dominios. [R53]
- Regiones con pLDDT muy bajo suelen coincidir con regiones intrínsecamente
  desordenadas; interpretarlas como “estructura” es un error. [R58]

En σ, σ1.1, los linkers y las colas N/C-terminales suelen salir con baja
confianza. Eso **no** significa que sean prescindibles: σ1.1 es funcional
(compuerta/autoinhibición) y los linkers determinan la disposición de dominios.
[R17] [R27] [R34]

### 6.3 Conformación dependiente de contexto

AlphaFold predice una conformación “promedio” condicionada por el MSA y la
plantilla; no modela el equilibrio entre conformaciones. En σ70, la forma libre
es compacta y la unida a RNAP es extendida; [R34] [R35] en σ28, FlgM estabiliza
una forma compacta. [R37] Por eso, la orientación de dominios en un modelo apo
no debe usarse como evidencia de interfaz in vivo.

Terwilliger et al. compararon predicciones AlphaFold con mapas
cristalográficos y encontraron casos de alta confianza que diferían en la
orientación global de dominios y en la conformación de cadenas laterales; su
recomendación es tratar las predicciones como hipótesis excepcionalmente útiles,
no como estructuras experimentales. [R57]

### 6.4 Ligandos, cofactores y complejos

AlphaFold2 no modela ligandos, cofactores ni modificaciones covalentes; [R57]
AlphaFold3 sí predice complejos con proteínas, ácidos nucleicos, iones y
moléculas pequeñas, con mejoras notables frente a *docking* especializado, pero
sigue siendo una predicción, no una medida. [R56] Para σ:

- un modelo de σ sola no contiene Zn²⁺, Mg²⁺ ni DNA;
- un modelo de σ–anti-σ (AF3) puede sugerir la interfaz, pero hay que validarla;
- un modelo de σ–DNA puede sugerir contactos, pero no sustituye un PWM
  experimental.

### 6.5 Búsqueda estructural (Foldseek) y sus límites

Foldseek permite buscar vecinos estructurales de forma rápida y sensible.
[R63] Es útil para detectar miembros divergentes de una familia o para
confirmar que un pliegue es σ70-like. Pero:

- un HTH de σ4 es un pliegue común; un hit a σ4 no convierte a la proteína en
  σ;
- un hit a un dominio aislado no informa la arquitectura completa;
- la similitud estructural no es ortología ni nombre nativo.

**Uso recomendado:** confirmación de familia/pliegue con arquitectura completa
y HMM; nunca asignación de nombre o promotor por *hit* estructural.

### 6.6 Sesgos del PDB

- **Fragmentos:** muchas entradas son dominios aislados (1KU7, 3UGO, 5FGM,
  6JHE, 2O8K, 2K9L), no la proteína completa.
- **Complejos inhibidores:** las σ alternativas suelen estar en complejos con
  anti-σ (1RP3, 1L0O, 1OR7, 5WUQ, 6IN7), no en la forma libre.
- **Cristalización:** iones y aditivos (Pb²⁺, Cd²⁺, K⁺) pueden aparecer sin
  significado biológico.
- **Resolución:** 4–6 Å en complejos grandes limita el detalle de cadenas
  laterales; no todos los contactos son visibles.
- **Organismos modelo:** *E. coli*, *T. aquaticus*, *B. subtilis*,
  *M. tuberculosis* dominan; la mayoría de las especies no tiene estructura.

**Tabla 6. Límite → consecuencia práctica → mitigación**

| Límite | Consecuencia | Mitigación |
|---|---|---|
| pLDDT bajo en σ1.1/linkers | No es “pérdida de función” | Anotar desorden; no penalizar |
| Modelo apo compacto | Orientación de dominios no es la de la holoenzima | Usar estructura de complejo o PAE con cautela |
| Sin ligandos en AF2 | No inferir Zn²⁺/ADP/DNA | Buscar en PDB/AlphaFill o AF3; validar |
| AF3 para complejos grandes | RNAP–σ (~450 kDa): accuracy no validada | Usar estructuras experimentales cuando existan |
| Foldseek con dominios aislados | Falsos positivos de familia | Exigir arquitectura completa y HMM |
| PDB sesgado a modelos | Especies no modelo sin estructura | AlphaFold + HMM + contexto; declarar incertidumbre |
| Baja resolución | Contactos laterales inciertos | No afirmar contactos específicos |

---

## 7. Features estructurales para el anotador

Esta sección traduce la evidencia a decisiones. Se prioriza por costo/beneficio
y por riesgo de sobreafirmar.

### 7.1 P0 — usar (alto valor, bajo riesgo)

**F1. Arquitectura de dominios (HMMER/Pfam).**
Escanear el panel σ70/σ54/ECF/anti-σ y reportar, por proteína: dominio, Pfam,
inicio, fin, e-value, bitscore y cobertura. Es la feature más informativa y la
más barata. [R62] Advertencia: PF04545 tiene falsos negativos declarados;
[PF04545] no convertir un “no hit” en “no hay dominio”.

**F2. Asignación de familia/grupo.**
Combinar HMM + arquitectura + grupo filogenético (cuando exista) para asignar
familia (σ70/σ54) y, dentro de σ70, grupo/ECF. La estructura entra solo como
confirmación. [R08] [R09] [R10]

**F3. Integridad estructural.**
Medir cobertura del modelo sobre la secuencia y presencia de los dominios
clave: σ2 y σ4 en σ70; CBD y DBD en σ54. Distinguir “modelo parcial” de “gen
truncado”. Si el gen está completo y el modelo no, es limitación del modelo; si
ambos están truncados, es candidato a pseudogen (requiere secuencia).

**F4. Mapa de interfaces.**
Anotar en la secuencia las regiones de interfaz conocidas: σ2.1–2.2 (core),
σ3.1–3.2 (core), σ4 (β-flap/activadores), σ1.1 (canal), AID/CBD (σ54),
linker σ2/σ4 (ECF). Reportar si una variante tiene deleciones o sustituciones
en estas regiones. [R13] [R17] [R27] [R28]

**F5. Calidad del modelo por dominio.**
pLDDT medio por dominio y PAE entre dominios. Usar umbrales solo como política
interna (p. ej., >70 para “modelo utilizable”), declarados como tales, no como
verdad biológica. [R53] [R57]

**F6. Desorden.**
Marcar regiones de pLDDT bajo como “desordenadas/no modeladas”, no como
“ausentes”. Especial atención a σ1.1 y linkers. [R58] [R34]

### 7.2 P1 — usar con cautela (hipótesis, no anotación)

**F7. Vecino estructural (Foldseek).**
Solo para confirmar pliegue/familia en casos divergentes; exigir arquitectura
completa y coherencia con HMM. No usar para nombre nativo ni promotor. [R63]

**F8. Bolsillo que solapa interfaz.**
Si un bolsillo P2Rank/fpocket solapa residuos de una interfaz conocida, marcar
“posible interfaz” con confianza baja; si no solapa nada, descartar. [R59] [R60]

**F9. Co-plegamiento σ–anti-σ (AF3).**
Generar hipótesis de pareja y de interfaz; validar con contexto génico
(operón, genes vecinos) y, si es posible, con datos experimentales. No anotar
“regulado por X” solo por un modelo. [R56] [R57]

### 7.3 P2 — no usar (riesgo alto, valor nulo)

- Volumen, *score* o *druggability* de un bolsillo como evidencia de función.
- *Docking* de moléculas pequeñas como evidencia de regulación.
- *Docking* de DNA o de un promotor como evidencia de especificidad.
- pLDDT como confianza funcional.
- Similitud estructural a un dominio aislado como familia o nombre.
- Presencia de un ligando en un homólogo como evidencia de que la consulta lo
  une.
- Interfaz predicha como interacción demostrada.

### 7.4 Esquema propuesto (JSON)

```json
{
  "structure_evidence": {
    "source": "PDB|AFDB|ESMFold",
    "id": "AF-EXAMPLE-F1",
    "model_confidence": {"mean_plddt": 82.1, "pae_domains": 4.2},
    "family_call": "sigma70",
    "group_call": "ECF",
    "domains": [
      {"name": "sigma70_r2", "pfam": "PF04542", "start": 45, "end": 130,
       "plddt": 88.0, "coverage": 0.97}
    ],
    "completeness": {"sigma2": true, "sigma4": true, "sigma54_cbd": null},
    "interfaces": {
      "rnap_core": ["sigma2.1-2.2", "sigma3.1-3.2", "sigma4"],
      "anti_sigma": "not_evaluated"
    },
    "pockets": [
      {"rank": 1, "score": 12.3, "overlap": "sigma4_dna_surface",
       "interpretation": "dna_binding_surface", "use": "none"}
    ],
    "verdict": "family_support",
    "caveats": ["sigma1.1_low_plddt", "no_ligand_in_model"]
  }
}
```

**Campos obligatorios:** `source`, `family_call`, `domains`, `completeness`,
`verdict`, `caveats`. **Campos prohibidos como evidencia de promotor:**
`pockets`, `ligands`, `docking`, `plddt` (como función), `foldseek_hit` (como
nombre).

### 7.5 Reglas de decisión

1. Si HMM y estructura coinciden en familia → `family_call` con confianza alta.
2. Si HMM no da hit pero la estructura y la arquitectura son σ70-like →
   `candidate_divergent`, revisión manual; **no** asignar nombre.
3. Si falta σ2 o σ4 en el modelo pero el gen está completo → limitación del
   modelo; no penalizar la proteína.
4. Si el gen tiene parada/frameshift y el modelo está truncado en el mismo
   punto → `pseudogene_candidate`; requiere secuencia y evidencia adicional.
5. Si un bolsillo solapa σ2/σ4/DBD → `dna_binding_surface`; no es promotor.
6. Si un bolsillo solapa una interfaz conocida → `interface_hypothesis`,
   confianza baja.
7. Si un bolsillo no solapa nada → `uninterpreted_cavity`; no se reporta.
8. Si no hay estructura → `not_available`; el resto del anotador sigue
   funcionando sin ella.

**Observación interna.** `scripts/sigma/105_build_pocket_dictionary.py` y
`scripts/sigma/111_foldseek_clusters.py` ya declaran explícitamente que no
asignan función; esta revisión respalda esa decisión y sugiere mantener los
pockets como evidencia de interfaz, nunca como feature de promotor.

---

## 8. Qué descartar como evidencia de promotor

**Tabla 7. Artefactos y evidencias que no autorizan una afirmación de promotor**

| Señal | Por qué no sirve | Qué sí autoriza |
|---|---|---|
| Volumen/*score* de bolsillo | Mide cavidad, no especificidad | Nada; a lo sumo, interfaz |
| *Docking* de DNA | No modela flexibilidad ni agua; puede dar poses plausibles sin ser ciertas | Hipótesis estructural para diseñar mutantes |
| *Docking* de molécula pequeña | No se conoce un ligando fisiológico de σ en el conjunto revisado | Nada |
| Similitud a σ4/HTH | Pliegue promiscuo; muchas proteínas lo tienen | Confirmar un dominio HTH, no familia |
| pLDDT alto | Confianza de plegamiento local | Calidad del modelo |
| pLDDT bajo | Desorden o incertidumbre | Región flexible/no modelada |
| Interfaz con RNAP | Informa unión, no promotor | Mecanismo de iniciación |
| Interfaz con anti-σ | Informa regulación | Mecanismo de secuestro |
| Presencia de Zn²⁺/ADP | Informa sensor/quinasa | Regulación |
| *Foldseek hit* a una σ conocida | Similitud de pliegue, no ortología | Confirmar familia/pliegue |
| Modelo AF de σ–DNA | Predicción, no medida | Hipótesis de contactos |
| Expresión de la σ | No dice qué promotores usa | Contexto fisiológico |

**Regla de oro:** si la afirmación es “esta σ reconoce el promotor P”, la
evidencia debe ser una PWM, un experimento de unión/transcripción o una
predicción calibrada con datos; la estructura puede ayudar a interpretar el
resultado, pero no lo genera.

---

## 9. Identidad/familia vs especificidad de promotor

### 9.1 Dos preguntas distintas

| Pregunta | Capa de evidencia | Métrica/artefacto | Afirmación permitida | Afirmación prohibida |
|---|---|---|---|---|
| ¿Es una σ? | HMM + arquitectura + estructura | Dominios, pliegue, cobertura | “Proteína de la familia σ70/σ54” | “Es la σ X” |
| ¿Qué familia/grupo? | HMM + filogenia + arquitectura | Grupo, contexto génico | “ECF del grupo Y” | “Reconoce el promotor Z” |
| ¿Qué ortólogo? | Ortología + contexto | Árbol, sintenía | “Ortólogo de SigW” | “Es SigW” |
| ¿Qué promotor? | PWM/experimento | PWM, unión, transcripción | “Promotor predicho –35/–10” | “Promotor demostrado” |
| ¿En qué condición? | Expresión/regulón/regulación | RNA-seq, ChIP, mutantes | “Hipótesis de uso en C” | “Se usa en C” |
| ¿Qué ligando? | Estructura del complejo + bioquímica | PDB, ITC, MS | “Sitio de Zn²⁺ en el anti-σ” | “La σ une Zn²⁺” |

### 9.2 Por qué la estructura no puede separar identidad de promotor

La identidad de familia es una propiedad de **pliegue + arquitectura +
homología**. La especificidad de promotor es una propiedad de **la superficie de
contacto con el DNA y su geometría**, que puede variar sin cambiar el pliegue.
Ejemplos:

- Dos σI parálogas tienen el mismo pliegue y distinta especificidad; la
  comparación de estructuras revela los determinantes, pero no se pueden
  predecir sin comparar. [R29]
- SigR y σE comparten pliegue; un solo residuo de R4 cambia la preferencia por
  –31′. [R31]
- Las ECF de distintos grupos comparten σ2/σ4 pero tienen motivos de promotor
  distintos; la clasificación filogenética, no la estructura, predice el motivo.
  [R08] [R09]
- σ28 y σ70 comparten pliegue, pero σ28 reconoce una –10 compuesta (extendida +
  core) con residuos específicos. [R65]
- σS y σ70 reconocen secuencias muy similares; la selectividad depende de
  posiciones –13/–14 y de la región 2.5. [R64]

### 9.3 Implicación para el score del anotador

La estructura **no debe sumar puntos a la evidencia de promotor**. Puede:

- habilitar/deshabilitar un modelo de promotor (p. ej., si es σ54, usar
  –24/–12);
- aportar contexto de familia para elegir el panel de PWM;
- marcar una variante como “dominio de unión a DNA incompleto”, lo que baja la
  confianza de cualquier predicción de promotor;
- aportar una hipótesis de interfaz para explicar por qué una σ no se usa.

Lo que no puede es convertir una predicción débil de promotor en fuerte por el
hecho de que el modelo sea bonito.

---

## 10. Vacíos y límites de esta revisión

1. **No hay estructura de una σ completa, libre y competente para DNA.** Las
   estructuras completas conocidas son complejos inhibidores (σ28/FlgM,
   σE/RseA) o conformaciones compactas (σ70 apo). Esto limita lo que se puede
   decir sobre la conformación activa.
2. **Cobertura desigual de alternativas.** σ28, σF, σE, σH, σK, σW, σI y σ54
   tienen estructuras; muchas otras σ alternativas solo tienen modelos.
3. **No hay benchmark sistemático de detectores de bolsillos en σ.** fpocket y
   P2Rank están validados para sitios de unión de ligandos; en σ, donde la
   mayoría de las cavidades son interfaces o surcos de DNA, su valor es
   heurístico y no está calibrado. [R59] [R60]
4. **La biología de ligandos de σ está poco explorada.** La mayoría de los
   “ligandos” observados son iones o nucleótidos unidos a anti-σ/RNAP; no se
   puede descartar que existan ligandos de σ no caracterizados, pero no hay
   evidencia estructural en el conjunto revisado.
5. **AlphaFold3 para complejos grandes.** La RNAP–σ es un complejo de ~450 kDa;
   la evaluación publicada de AF3 se centra en complejos de tamaño moderado y su
   accuracy en un complejo de este tamaño no está establecida en esta revisión.
   [R56]
6. **Motivos de promotor de muchas ECF son inferidos, no medidos.** La
   clasificación filogenética predice motivos, pero la validación experimental
   es limitada. [R08] [R09]
7. **Sesgo de organismo modelo.** Las generalizaciones estructurales provienen
   de unas pocas especies; la extrapolación a genomas no modelo debe declararse
   como hipótesis.

**Recomendación operativa.** Usar la estructura solo como evidencia de
familia/arquitectura/integridad/interfaz. En el anotador: (1) HMM + arquitectura
como base; (2) estructura para confirmar familia y detectar truncamientos;
(3) pLDDT/PAE para calidad y desorden; (4) bolsillos solo como hipótesis de
interfaz; (5) nunca derivar PWM, promotor ni condición de la estructura. Si hay
que elegir una sola feature estructural, elegir la arquitectura de dominios; si
hay que elegir una sola prohibición, prohibir el bolsillo como evidencia de
promotor.

---

## 11. Referencias verificadas

> Todas las referencias fueron resueltas contra Crossref (DOI) y NCBI
> E-utilities (PMID), y el título devuelto se comparó con la cita. Las entradas
> PDB se resolvieron contra la API de RCSB; las Pfam, contra la API de InterPro.
> No se incluye ninguna referencia no verificada. Ver Anexo A.

| ID | Referencia | PMID | DOI |
|---|---|---|---|
| R01 | Paget M. *Bacterial Sigma Factors and Anti-Sigma Factors: Structure, Function and Distribution*. Biomolecules. 2015;5(3):1245-1265. | 26131973 | 10.3390/biom5031245 |
| R02 | Feklistov A, Sharon BD, Darst SA, et al. *Bacterial Sigma Factors: A Historical, Structural, and Genomic Perspective*. Annual Review of Microbiology. 2014;68(1):357-376. | 25002089 | 10.1146/annurev-micro-092412-155737 |
| R03 | Gruber TM, Gross CA. *Multiple Sigma Subunits and the Partitioning of Bacterial Transcription Space*. Annual Review of Microbiology. 2003;57(1):441-466. | 14527287 | 10.1146/annurev.micro.57.030502.090913 |
| R04 | Lonetto M, Gribskov M, Gross CA. *The sigma 70 family: sequence conservation and evolutionary relationships*. Journal of Bacteriology. 1992;174(12):3843-3849. | 1597408 | 10.1128/jb.174.12.3843-3849.1992 |
| R05 | Helmann JD, Chamberlin MJ. *Structure and function of bacterial sigma factors*. Annual Review of Biochemistry. 1988;57(1):839-872. | 3052291 | 10.1146/annurev.bi.57.070188.004203 |
| R06 | Wösten MMSM. *Eubacterial sigma-factors*. FEMS Microbiology Reviews. 1998;22(3):127-150. | 9818380 | 10.1111/j.1574-6976.1998.tb00364.x |
| R07 | Helmann JD. *The extracytoplasmic function (ECF) sigma factors*. Advances in Microbial Physiology. 2002;46:47-110. | 12073657 | 10.1016/S0065-2911(02)46002-X |
| R08 | Staroń A, Sofia HJ, Dietrich S, et al. *The third pillar of bacterial signal transduction: classification of the extracytoplasmic function (ECF) σ factor protein family*. Molecular Microbiology. 2009;74(3):557-581. | 19737356 | 10.1111/j.1365-2958.2009.06870.x |
| R09 | Casas-Pastor D, Müller RR, Jaenicke S, et al. *Expansion and re-classification of the extracytoplasmic function (ECF) σ factor family*. Nucleic Acids Research. 2021;49(2):986-1005. | 33398323 | 10.1093/nar/gkaa1229 |
| R10 | Paget MS, Helmann JD. *The sigma70 family of sigma factors*. Genome Biology. 2003;4(1):203. | 12540296 | 10.1186/gb-2003-4-1-203 |
| R11 | Barrios H, Valderrama B, Morett E. *Compilation and analysis of σ54-dependent promoter sequences*. Nucleic Acids Research. 1999;27(22):4305-4313. | 10536136 | 10.1093/nar/27.22.4305 |
| R12 | Bush M, Dixon R. *The Role of Bacterial Enhancer Binding Proteins as Specialized Activators of σ54-Dependent Transcription*. Microbiology and Molecular Biology Reviews. 2012;76(3):497-529. | 22933558 | 10.1128/mmbr.00006-12 |
| R13 | Campbell EA, Muzzin O, Chlenov M, et al. *Structure of the Bacterial RNA Polymerase Promoter Specificity σ Subunit*. Molecular Cell. 2002;9(3):527-539. | 11931761 | 10.1016/s1097-2765(02)00470-7 |
| R14 | Feklistov A, Darst SA. *Structural Basis for Promoter −10 Element Recognition by the Bacterial RNA Polymerase σ Subunit*. Cell. 2011;147(6):1257-1269. | 22136875 | 10.1016/j.cell.2011.10.041 |
| R15 | Zhang Y, Feng Y, Chatterjee S, et al. *Structural Basis of Transcription Initiation*. Science. 2012;338(6110):1076-1080. | 23086998 | 10.1126/science.1227786 |
| R16 | Vassylyev DG, Sekine Si, Laptenko O, et al. *Crystal structure of a bacterial RNA polymerase holoenzyme at 2.6 Å resolution*. Nature. 2002;417(6890):712-719. | 12000971 | 10.1038/nature752 |
| R17 | Murakami KS. *X-ray Crystal Structure of Escherichia coli RNA Polymerase σ70 Holoenzyme*. Journal of Biological Chemistry. 2013;288(13):9126-9134. | 23389035 | 10.1074/jbc.m112.430900 |
| R18 | Zuo Y, Steitz TA. *Crystal Structures of the E. coli Transcription Initiation Complexes with a Complete Bubble*. Molecular Cell. 2015;58(3):534-540. | 25866247 | 10.1016/j.molcel.2015.03.010 |
| R19 | Bae B, Feklistov A, Lass-Napiorkowska A, et al. *Structure of a bacterial RNA polymerase holoenzyme open promoter complex*. eLife. 2015;4. | 26349032 | 10.7554/elife.08504 |
| R20 | Narayanan A, Vago FS, Li K, et al. *Cryo-EM structure of Escherichia coli σ70 RNA polymerase and promoter DNA complex revealed a role of σ non-conserved region during the open complex formation*. Journal of Biological Chemistry. 2018;293(19):7367-7375. | 29581236 | 10.1074/jbc.ra118.002161 |
| R21 | Yang Y, Darbari VC, Zhang N, et al. *Structures of the RNA polymerase-σ54 reveal new and conserved regulatory strategies*. Science. 2015;349(6250):882-885. | 26293966 | 10.1126/science.aab1478 |
| R22 | Glyde R, Ye F, Darbari VC, et al. *Structures of RNA Polymerase Closed and Intermediate Complexes Reveal Mechanisms of DNA Opening and Transcription Initiation*. Molecular Cell. 2017;67(1):106-116.e4. | 28579332 | 10.1016/j.molcel.2017.05.010 |
| R23 | Glyde R, Ye F, Jovanovic M, et al. *Structures of Bacterial RNA Polymerase Complexes Reveal the Mechanism of DNA Loading and Transcription Initiation*. Molecular Cell. 2018;70(6):1111-1120.e3. | 29932903 | 10.1016/j.molcel.2018.05.021 |
| R24 | Doucleff M, Pelton JG, Lee PS, et al. *Structural Basis of DNA Recognition by the Alternative Sigma-factor, σ54*. Journal of Molecular Biology. 2007;369(4):1070-1078. | 17481658 | 10.1016/j.jmb.2007.04.019 |
| R25 | Hong E, Doucleff M, Wemmer DE. *Structure of the RNA Polymerase Core-Binding Domain of σ54 Reveals a Likely Conformational Fracture Point*. Journal of Molecular Biology. 2009;390(1):70-82. | 19426742 | 10.1016/j.jmb.2009.04.070 |
| R26 | Campbell EA, Kamath S, Rajashankar KR, et al. *Crystal structure of Aquifex aeolicus σN bound to promoter DNA and the structure of σN-holoenzyme*. Proceedings of the National Academy of Sciences. 2017;114(10). | 28223493 | 10.1073/pnas.1619464114 |
| R27 | Fang C, Li L, Shen L, et al. *Structures and mechanism of transcription initiation by bacterial ECF factors*. Nucleic Acids Research. 2019;47(13):7094-7104. | 31131408 | 10.1093/nar/gkz470 |
| R28 | Li L, Fang C, Zhuang N, et al. *Structural basis for transcription initiation by bacterial ECF σ factors*. Nature Communications. 2019;10(1). | 30858373 | 10.1038/s41467-019-09096-y |
| R29 | Li J, Zhang H, Li D, et al. *Structure of the transcription open complex of distinct σI factors*. Nature Communications. 2023;14(1). | 37833284 | 10.1038/s41467-023-41796-4 |
| R30 | Kwon E, Devkota SR, Pathak D, et al. *Structural analysis of the recognition of the -35 promoter element by SigW from Bacillus subtilis*. PLOS ONE. 2019;14(8):e0221666. | 31461489 | 10.1371/journal.pone.0221666 |
| R31 | Kim KY, Park JK, Park S. *In Streptomyces coelicolor SigR, methionine at the -35 element interacting region 4 confers the -31′-adenine base selectivity*. Biochemical and Biophysical Research Communications. 2016;470(2):257-262. | 26775842 | 10.1016/j.bbrc.2016.01.075 |
| R32 | Gaballa A, Guariglia-Oropeza V, Dürr F, et al. *Modulation of extracytoplasmic function (ECF) sigma factor promoter selectivity by spacer region sequence*. Nucleic Acids Research. 2017;46(1):134-145. | 29069433 | 10.1093/nar/gkx953 |
| R33 | Sallai L, Tucker PA. *Crystal structure of the central and C-terminal domain of the σ54-activator ZraR*. Journal of Structural Biology. 2005;151(2):160-170. | 16005641 | 10.1016/j.jsb.2005.05.006 |
| R34 | Schwartz EC, Shekhtman A, Dutta K, et al. *A Full-Length Group 1 Bacterial Sigma Factor Adopts a Compact Structure Incompatible with DNA Binding*. Chemistry & Biology. 2008;15(10):1091-1103. | 18940669 | 10.1016/j.chembiol.2008.09.008 |
| R35 | Joron K, Zamel J, Kalisman N, et al. *Evidence for a compact σ70 conformation in vitro and in vivo*. iScience. 2024;27(6):110140. | 38957792 | 10.1016/j.isci.2024.110140 |
| R36 | Campbell EA, Tupy JL, Gruber TM, et al. *Crystal Structure of Escherichia coli σE with the Cytoplasmic Domain of Its Anti-σ RseA*. Molecular Cell. 2003;11(4):1067-1078. | 12718891 | 10.1016/s1097-2765(03)00148-5 |
| R37 | Sorenson MK, Ray SS, Darst SA. *Crystal Structure of the Flagellar σ/Anti-σ Complex σ28/FlgM Reveals an Intact σ Factor in an Inactive Conformation*. Molecular Cell. 2004;14(1):127-138. | 15068809 | 10.1016/s1097-2765(04)00150-9 |
| R38 | Campbell EA, Masuda S, Sun JL, et al. *Crystal Structure of the Bacillus stearothermophilus Anti-σ Factor SpoIIAB with the Sporulation σ Factor σF*. Cell. 2002;108(6):795-807. | 11955433 | 10.1016/s0092-8674(02)00662-1 |
| R39 | Masuda S, Murakami KS, Wang S, et al. *Crystal Structures of the ADP and ATP Bound Forms of the Bacillus Anti-σ Factor SpoIIAB in Complex with the Anti-anti-σ SpoIIAA*. Journal of Molecular Biology. 2004;340(5):941-956. | 15236958 | 10.1016/j.jmb.2004.05.040 |
| R40 | Lambert LJ, Wei Y, Schirf V, et al. *T4 AsiA blocks DNA recognition by remodeling σ70 region 4*. The EMBO Journal. 2004;23(15):2952-2962. | 15257291 | 10.1038/sj.emboj.7600312 |
| R41 | Patikoglou GA, Westblade LF, Campbell EA, et al. *Crystal Structure of the Escherichia coli Regulator of σ70, Rsd, in Complex with σ70 Domain 4*. Journal of Molecular Biology. 2007;372(3):649-659. | 17681541 | 10.1016/j.jmb.2007.06.081 |
| R42 | Campbell EA, Greenwell R, Anthony JR, et al. *A Conserved Structural Module Regulates Transcriptional Responses to Diverse Stress Signals in Bacteria*. Molecular Cell. 2007;27(5):793-805. | 17803943 | 10.1016/j.molcel.2007.07.009 |
| R43 | Maillard AP, Girard E, Ziani W, et al. *The Crystal Structure of the Anti-σ Factor CnrY in Complex with the σ Factor CnrH Shows a New Structural Class of Anti-σ Factors Targeting Extracytoplasmic Function σ Factors*. Journal of Molecular Biology. 2014;426(12):2313-2327. | 24727125 | 10.1016/j.jmb.2014.04.003 |
| R44 | Shukla J, Gupta R, Thakur KG, et al. *Structural basis for the redox sensitivity of the Mycobacterium tuberculosis SigK–RskA σ–anti-σ complex*. Acta Crystallographica Section D Biological Crystallography. 2014;70(4):1026-1036. | 24699647 | 10.1107/s1399004714000121 |
| R45 | Devkota SR, Kwon E, Ha SC, et al. *Structural insights into the regulation of Bacillus subtilis SigW activity by anti-sigma RsiW*. PLOS ONE. 2017;12(3):e0174284. | 28319136 | 10.1371/journal.pone.0174284 |
| R46 | Li S, Lou X, Xu Y, et al. *Structural basis for the recognition of MucA by MucB and AlgU in Pseudomonas aeruginosa*. The FEBS Journal. 2019;286(24):4982-4994. | 31297938 | 10.1111/febs.14995 |
| R47 | Li W, Bottrill AR, Bibb MJ, et al. *The Role of Zinc in the Disulphide Stress-regulated Anti-sigma Factor RsrA from Streptomyces coelicolor*. Journal of Molecular Biology. 2003;333(2):461-472. | 14529630 | 10.1016/j.jmb.2003.08.038 |
| R48 | Rajasekar KV, Zdanowski K, Yan J, et al. *The anti-sigma factor RsrA responds to oxidative stress by reburying its hydrophobic core*. Nature Communications. 2016;7(1). | 27432510 | 10.1038/ncomms12194 |
| R49 | Wollmann P, Zeth K. *The Structure of RseB: A Sensor in Periplasmic Stress Response of E. coli*. Journal of Molecular Biology. 2007;372(4):927-941. | 17692869 | 10.1016/j.jmb.2007.06.039 |
| R50 | Lima S, Guo MS, Chaba R, et al. *Dual Molecular Signals Mediate the Bacterial Response to Outer-Membrane Stress*. Science. 2013;340(6134):837-841. | 23687042 | 10.1126/science.1235358 |
| R51 | Wilken C, Kitzing K, Kurzbauer R, et al. *Crystal Structure of the DegS Stress Sensor*. Cell. 2004;117(4):483-494. | 15137941 | 10.1016/s0092-8674(04)00454-4 |
| R52 | Yokoyama T, Niinae T, Tsumagari K, et al. *The Escherichia coli S2P intramembrane protease RseP regulates ferric citrate uptake by cleaving the sigma factor regulator FecR*. Journal of Biological Chemistry. 2021;296:100673. | 33865858 | 10.1016/j.jbc.2021.100673 |
| R53 | Jumper J, Evans R, Pritzel A, et al. *Highly accurate protein structure prediction with AlphaFold*. Nature. 2021;596(7873):583-589. | 34265844 | 10.1038/s41586-021-03819-2 |
| R54 | Varadi M, Anyango S, Deshpande M, et al. *AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models*. Nucleic Acids Research. 2021;50(D1):D439-D444. | 34791371 | 10.1093/nar/gkab1061 |
| R55 | Varadi M, Bertoni D, Magana P, et al. *AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences*. Nucleic Acids Research. 2024;52(D1):D368-D375. | 37933859 | 10.1093/nar/gkad1011 |
| R56 | Abramson J, Adler J, Dunger J, et al. *Accurate structure prediction of biomolecular interactions with AlphaFold 3*. Nature. 2024;630(8016):493-500. | 38718835 | 10.1038/s41586-024-07487-w |
| R57 | Terwilliger TC, Liebschner D, Croll TI, et al. *AlphaFold predictions are valuable hypotheses and accelerate but do not replace experimental structure determination*. Nature Methods. 2024;21(1):110-116. | 38036854 | 10.1038/s41592-023-02087-4 |
| R58 | Ruff KM, Pappu RV. *AlphaFold and Implications for Intrinsically Disordered Proteins*. Journal of Molecular Biology. 2021;433(20):167208. | 34418423 | 10.1016/j.jmb.2021.167208 |
| R59 | Le Guilloux V, Schmidtke P, Tuffery P. *Fpocket: An open source platform for ligand pocket detection*. BMC Bioinformatics. 2009;10(1). | 19486540 | 10.1186/1471-2105-10-168 |
| R60 | Krivák R, Hoksza D. *P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure*. Journal of Cheminformatics. 2018;10(1). | 30109435 | 10.1186/s13321-018-0285-8 |
| R61 | Berman HM. *The Protein Data Bank*. Nucleic Acids Research. 2000;28(1):235-242. | 10592235 | 10.1093/nar/28.1.235 |
| R62 | Paysan-Lafosse T, Blum M, Chuguransky S, et al. *InterPro in 2022*. Nucleic Acids Research. 2022;51(D1):D418-D427. | 36350672 | 10.1093/nar/gkac993 |
| R63 | van Kempen M, Kim SS, Tumescheit C, et al. *Fast and accurate protein structure search with Foldseek*. Nature Biotechnology. 2023;42(2):243-246. | 37156916 | 10.1038/s41587-023-01773-0 |
| R64 | Becker G, Hengge-Aronis R. *What makes an Escherichia coli promoter sigmaS dependent? Role of the -13/-14 nucleotide promoter positions and region 2.5 of sigmaS*. Molecular Microbiology. 2001;39(5):1153-1165. | 11251833 | 10.1046/j.1365-2958.2001.02313.x |
| R65 | Koo B, Rhodius VA, Campbell EA, et al. *Mutational analysis of Escherichia coli σ28 and its target promoters reveals recognition of a composite −10 region, comprised of an ‘extended −10’ motif and a core −10 element*. Molecular Microbiology. 2009;72(4):830-843. | 19400790 | 10.1111/j.1365-2958.2009.06691.x |
| R66 | Schnoes AM, Brown SD, Dodevski I, et al. *Annotation Error in Public Databases: Misannotation of Molecular Function in Enzyme Superfamilies*. PLoS Computational Biology. 2009;5(12):e1000605. | 20011109 | 10.1371/journal.pcbi.1000605 |

**Notas de verificación:**

- **R06 (Wösten 1998):** el DOI de Wiley y el de Elsevier resuelven al mismo
  artículo; se usa el de Wiley, verificado en la tabla. El PMID 9818380
  corresponde al título “Eubacterial sigma-factors”.
- **R19 (Bae 2015):** el DOI correcto es `10.7554/eLife.08504`; Crossref
  devuelve sub-DOI de *decision letter* y *author response*, que se
  descartaron.
- **R57 (Terwilliger et al.):** publicado en *Nature Methods* 21(1):110–116
  (2024); Crossref lo registra *online* en 2023. Se usa 2024 en el texto.
- **R64 (Becker & Hengge-Aronis 2001):** el PMID se recuperó por título; el
  DOI resuelve correctamente en Crossref.
- **R44 (Shukla et al. 2014):** el título de Crossref pega “the” y
  “Mycobacterium”; se corrigió el espacio en la tabla.

---

## Anexo A. Registro de verificación

**Fecha:** 2026-11-26 (America/Santiago).

**Herramientas:**

- Crossref REST (`/works/{DOI}`, `query.bibliographic`).
- NCBI E-utilities (`esearch`, `esummary`, `efetch`; `db=pubmed`).
- RCSB PDB REST (`/rest/v1/core/entry/{ID}`) y búsqueda `search.rcsb.org`.
- InterPro/Pfam API (`/api/entry/pfam/{PF...}`).
- Europe PMC (búsqueda inicial; se usó para localizar candidatos, no como
  única fuente de verdad).

**Controles aplicados:**

1. Cada DOI fue resuelto y el título/año se comparó con la cita.
2. Cada PMID fue resuelto y el título se comparó; se descartaron PMID
   incorrectos.
3. Cada PDB fue consultado por ID; se verificaron título, método, resolución,
   ligandos y cita primaria.
4. Cada Pfam fue consultado por accession; se transcribió el nombre y la
   definición.
5. Cuando el abstract no contenía la afirmación, se restringió la afirmación a
   lo que el abstract sostiene; no se infirió de títulos.
6. No se incluyeron preprints, tesis ni páginas no institucionales.

**Correcciones y descartes documentados:**

- `1SIG` es un fragmento de σ70 de *E. coli* (1996), no el complejo σ4–DNA de
  Campbell 2002; se usa `1KU7` para el complejo σ4–DNA de –35.
- `5BY3` y `5OQ0` se probaron como posibles σ54–RNAP y resultaron ser otras
  proteínas; se descartaron.
- `5OR5` figura en el PDB como “not published” (estructura NMR de una región 2
  de σE modificada). Se usa solo como evidencia de contactos R2–DNA; no se cita
  como artículo ni se usa para afirmar especificidad fisiológica.
- `4YG2` se depositó en 2015 pero la publicación es Murakami 2013 [R17].
- `6GH5` y `8I23` son cryo-EM; su resolución (3.4 y 3.0 Å) no permite afirmar
  contactos de cadenas laterales con la misma confianza que una estructura de
  1.8 Å.
- `4NQW` contiene Cd²⁺, un ion de cristalización; no se interpreta como ligando
  fisiológico.
- `1IW7` contiene Pb²⁺, un ion de cristalización; no se interpreta como ligando
  fisiológico.

**Cobertura:** 66 referencias (46 primarias, 20 revisiones/recursos) + 30
entradas PDB (29 filas; 1SOT y 1SOZ comparten fila) + 8 entradas Pfam. **Referencias no verificadas incluidas:** 0.

---

## Anexo B. Estrategia de búsqueda y límites

**Bases:** PubMed/NCBI, Crossref, Europe PMC, RCSB PDB, InterPro/Pfam.

**Términos principales:**

- `sigma factor structure`, `sigma70 structure`, `sigma54 structure`,
  `ECF sigma structure`, `RNA polymerase sigma holoenzyme`,
  `sigma factor anti-sigma complex`, `sigma region 1.1`, `sigma region 4`,
  `-10 element recognition`, `-24 promoter sigma54`, `RpoN box`,
  `anti-sigma zinc`, `ZAS anti-sigma`, `RsrA`, `RseA`, `RsiW`, `ChrR`,
  `MucA`, `SpoIIAB`, `FlgM`, `AsiA`, `Rsd`.
- `AlphaFold`, `AlphaFold Protein Structure Database`, `AlphaFold 3`,
  `pLDDT`, `PAE`, `intrinsically disordered`.
- `fpocket`, `P2Rank`, `ligand binding site prediction`,
  `protein structure search Foldseek`.
- `sigma factor promoter specificity`, `ECF promoter motif`,
  `sigma54 -24/-12`, `sigma28 extended -10`, `sigmaS -13/-14`.

**Criterios de inclusión:** artículo revisado por pares con DOI/PMID
verificable; estructura PDB con cita primaria; entrada Pfam/InterPro con
accession.

**Criterios de exclusión:** preprints no publicados, tesis, páginas web no
institucionales, referencias sin DOI/PMID verificable, estructuras sin cita
primaria.

**Límites:** la búsqueda no fue sistemática (no PRISMA); se priorizaron
estructuras y revisiones representativas, no un censo exhaustivo. La
literatura sobre ligandos de σ es escasa; la ausencia de evidencia no es
evidencia de ausencia.

---

## Anexo C. Matriz de afirmaciones clave

| Afirmación | Tipo de evidencia | Fuente | Fuerza |
|---|---|---|---|
| σ70 y σ54 son familias distintas | Revisión + estructuras | [R01] [R02] [R21] [R24] | Alta |
| σ70 tiene grupos 1–4; ECF es grupo 4 | Revisión | [R10] [R08] [R09] | Alta |
| σ2 reconoce –10 con *base flipping* | Estructura + bioquímica | [R14] [R15] | Alta |
| σ4 reconoce –35 por HTH | Estructura + genética | [R13] [R30] [R31] | Alta |
| σ54 reconoce –24 por RpoN box | Estructura NMR | [R24] | Alta |
| σ54 requiere bEBP ATPasa | Revisión + estructura | [R12] [R33] | Alta |
| ECF usan linker σ2/σ4 en el canal activo | Estructura cryo-EM | [R27] [R28] | Alta |
| La σ libre es compacta y autoinhibida | Estructura + cross-linking | [R34] [R35] | Alta |
| RseA ocluye unión al core | Estructura | [R36] | Alta |
| FlgM entierra la hélice –35 | Estructura | [R37] | Alta |
| SpoIIAB es GHKL y une ADP/ATP | Estructura | [R38] [R39] | Alta |
| RsrA usa Zn²⁺ para sensar redox | Bioquímica + estructura | [R47] [R48] | Alta |
| RsiW une Zn²⁺ pero la oxidación no disocia | Estructura | [R45] | Alta |
| Un solo residuo de R4 cambia –31′ | Estructura + mutantes | [R31] | Alta |
| σI parálogas difieren en especificidad | Estructuras comparadas | [R29] | Alta |
| AlphaFold no modela ligandos | Revisión | [R57] | Alta |
| pLDDT bajo ≈ desorden | Revisión | [R58] | Alta |
| fpocket/P2Rank predicen ligandabilidad | Método | [R59] [R60] | Alta |
| No hay ligando drug-like de σ en el PDB revisado | Censo PDB | Tabla 3 | Media (acotada) |
| La estructura no predice el PWM | Inferencia de la revisión | §4, §9 | Alta |
| La estructura no dice qué σ se usa in vivo | Inferencia de la revisión | §4.2 | Alta |

---

## Anexo D. Glosario mínimo

- **σ70:** familia principal de factores sigma; cuatro grupos; incluye ECF.
- **σ54:** familia distinta; promotores –24/–12; requiere activador ATPasa.
- **ECF:** *extracytoplasmic function*; grupo 4 de σ70; el grupo más numeroso
  de sigma alternativas.
- **Región 1.1:** dominio N-terminal de σ70; compuerta/autoinhibición.
- **Región 2:** dominio más conservado; –10 y unión al core.
- **Región 3:** –10 extendida y unión al core.
- **Región 4:** HTH de unión a –35.
- **AID/CBD/DBD:** dominios de σ54 (activador, core, DNA).
- **RpoN box:** motivo conservado de σ54 que lee –24.
- **Anti-σ:** proteína que secuestra/inhibe una σ.
- **ZAS:** familia de anti-σ con Zn²⁺ (p. ej., RsrA, RsiW).
- **bEBP:** activador dependiente de ATP de σ54.
- **pLDDT:** confianza por residuo de AlphaFold.
- **PAE:** error alineado predicho entre residuos.
- **PWM:** matriz de pesos de posición; modelo de motivo.
