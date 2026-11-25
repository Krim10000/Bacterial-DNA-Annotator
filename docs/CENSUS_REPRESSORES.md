# Censo de represores: caracterizar todos, no revisar uno

> **Cambio de enfoque (2026-11-24).** El objetivo no es resolver un represor
> aislado ni cristalizar proteínas. El objetivo es **caracterizar todos los
> represores anotados** de un genoma (y de los 400 genomas del proyecto) con un
> censo exhaustivo, determinista y con confianza explícita.
>
> Este documento reemplaza la lógica de "piloto de 20 casos" como **alcance
> final**. El piloto sigue existiendo, pero solo como **calibración** del censo.
>
> Marco conceptual: [`FRAMEWORK_INFERENCIA_ACOTADA.md`](FRAMEWORK_INFERENCIA_ACOTADA.md).
> Brechas técnicas: [`FRAMEWORK_GAP_ANALYSIS.md`](FRAMEWORK_GAP_ANALYSIS.md).
> Aprendizaje en K-12: [`K12_ATLAS_SUBCASOS.md`](K12_ATLAS_SUBCASOS.md).
> Backtest contra RegulonDB: [`K12_BACKTEST.md`](K12_BACKTEST.md).
> Flujo de trabajo revisado: [`FLUJO_TRABAJO_REVISADO.md`](FLUJO_TRABAJO_REVISADO.md).

---

## 0. El cambio de enfoque

### 0.1 Lo que se abandona

- Cristalizar proteínas nuevas.
- Depender de AF3/RF2NA para cada represor.
- Revisar un TF por vez.
- Exigir un ΔΔG medido para poder anotar.
- Truncar la salida a "top 3–5" sin mostrar el resto.

### 0.2 Lo que se adopta

Un **censo**:

```text
para cada genoma:
  para cada represor anotado:
    para cada región reguladora del genoma:
      puntuar el sitio con el motivo de su familia
      filtrar por conservación evolutiva
      filtrar por contexto del operón
      asignar confianza
    para cada operón diana:
      proponer efectores por química del regulón
```

El resultado no es "la diana de LacI". Es una tabla completa:

```text
(genoma, represor, operador candidato, operón diana, efector candidato,
 score, evidencia, confianza, estado)
```

### 0.3 Por qué esto es coherente con la evolución

La evolución no optimizó un represor para una secuencia única en abstracto.
Optimizó un **sistema regulador** dentro de un genoma concreto:

- un motivo de unión con una geometría y una degeneración concretas;
- un promotor con una arquitectura concreta;
- un operón con una función metabólica concreta;
- una vecindad génica que se conserva o se pierde junta.

Por eso el censo no busca "la secuencia que encaja", sino **la combinación
motivo + conservación + contexto + función** que se repite en genomas
emparentados. Eso es leer la evolución, no adivinarla.

> **Matiz:** evolución optimiza, pero también produce redundancia, *cross-talk*,
> sitios crípticos y deriva. El censo debe separar señal funcional de ruido; la
> conservación es el mejor filtro disponible, no una garantía.

---

## 1. Universo del censo

### 1.1 Qué se cuenta como represor

No los 1.555 productos que dicen "repressor" sin más. Hay que separar:

| Clase | Ejemplo | ¿Entra al censo de operadores? |
|---|---|---|
| Represor secuencia-específico | LacI, TetR, GntR, MarR, PurR | **Sí** |
| Represor de familia con DBD | LexA, ArgR, HrcA, Rex | **Sí** |
| Proteína arquitectónica | H-NS, Fis, IHF, HU | **No** como operador discreto; se anota aparte |
| Antirepresor fágico | KilAC, AntA/AntB | **No**; no une ADN como represor |
| Factor sigma / activador | RpoS, CRP | **No** en este censo; otra capa |

**Consecuencia:** de los 1.555 "repressor" hay que quedarse con los que tienen
un DBD y un operador discreto. El resto no se pierde: se clasifica y se reporta
como "no aplica a censo de operadores".

### 1.2 Qué se cuenta como región reguladora

No todo el genoma, pero sí **todas** las regiones plausibles:

- promotores de todos los operones (`promoters.tsv`);
- regiones intergénicas con TSS o caja −10/−35;
- región upstream del gen del represor (divergón);
- región downstream si el operón es divergente;
- no se incluyen CDS ni regiones codificantes, salvo que haya un promotor
  interno anotado.

**El espacio es finito y enumerable.** Para un genoma de 5 Mb:

- ~1,5 × 10⁶ ventanas intergénicas por hebra;
- ~10⁷ posiciones si se incluye todo el genoma;
- con 400 genomas, ~6 × 10⁸ ventanas.

Eso es grande, pero no infinito. Se puede recorrer.

---

## 2. Motor del censo: motivo + conservación + contexto

### 2.1 Capa 1 — Biblioteca de motivos por familia

No se deriva un PWM por estructura para cada proteína. Se usa una biblioteca de
motivos por familia, que es más estable y ya existe para muchas familias.

| Fuente | Qué aporta | Uso |
|---|---|---|
| RegulonDB | 82 PWM de *E. coli* | Semilla y calibración |
| CollecTF | Sitios curados de bacterias | Validación y ampliación |
| DBTBS | *Bacillus subtilis* | Validación |
| RegPrecise | Regulones y efectores | Contexto y motivos |
| PRODORIC | Reguladores procariotas | Ampliación |
| Literatura de familia | Motivos canónicos | Prior de familia |

Para familias sin PWM curado: **descubrimiento comparativo**. Se toman las
regiones upstream de ortólogos en los 400 genomas, se alinean y se busca un
motivo conservado con MEME/DREME. Eso es exactamente "leer la evolución".

**Salida:** `motif_library.tsv` con `family`, `motif_id`, `pwm`, `width`,
`n_sites`, `source`, `quality`, `genomes_supported`.

### 2.2 Capa 2 — Escaneo exhaustivo

Para cada genoma y cada motivo de las familias presentes:

- escanear **todas** las regiones reguladoras, ambas hebras;
- calcular score PWM, p-valor y FDR;
- guardar **todos** los hits que pasen el umbral, no solo el top-1;
- registrar coordenadas, hebra, secuencia, distancia al TSS y operón.

No se muestrea. No se elige un subconjunto "representativo". Se enumera.

**Salida:** `motif_hits.tsv` con todos los hits y su p-valor.

### 2.3 Capa 3 — Conservación evolutiva

Un hit de PWM es barato; un hit conservado en ortólogos es caro de explicar por
azar. Para cada gen/operón:

1. Construir grupos de ortólogos con MMseqs2/DIAMOND sobre las proteínas de los
   400 genomas.
2. Alinear las regiones promotoras de los ortólogos.
3. Verificar si el hit aparece en la misma posición y orientación.
4. Calcular:
   - `conservation_score`: fracción de ortólogos con hit alineado;
   - `synteny_score`: conservación de la vecindad génica;
   - `phylogenetic_support`: número de linajes independientes.

**Salida:** `conserved_sites.tsv` con el hit y su soporte filogenético.

### 2.4 Capa 4 — Contexto genómico

El sitio no vive solo. Se puntúa:

- **operón diana**: genes corriente abajo y su función;
- **divergón**: orientación del gen del represor respecto del operón;
- **proximidad**: distancia al gen del represor;
- **coherencia funcional**: familia del represor vs función del operón;
- **arquitectura**: tamaño del operón, líder, atenuador, terminador.

**Salida:** `site_context.tsv`.

### 2.5 Capa 5 — Efectores candidatos

Sin docking obligatorio y sin cristalografía:

- **GBA d=0**: metabolitos que participan directamente en reacciones del
  operón regulado; excluir cofactores ubicuos (ATP, NAD, H₂O, H⁺);
- **prior de familia**: efector conocido de la familia (LacI/GalR → azúcares;
  TetR → tetraciclina; Fur → metales; MarR → salicilato);
- **coherencia de ruta**: mismo módulo metabólico;
- **conservación del regulón**: si el operón y el TF se conservan juntos.

**Salida:** `effector_candidates.tsv` con top 10 y `status=hypothesis`.

### 2.6 Capa 6 — Integración

\[
S = w_m \cdot \text{motif} + w_c \cdot \text{conservation} +
w_x \cdot \text{context} + w_f \cdot \text{function} + w_p \cdot \text{proximity}
\]

Los pesos se calibran en K-12 y se reporta su sensibilidad. La salida incluye
**todos** los sitios que pasan el umbral, ordenados, con confianza.

---

## 3. Cómo se ve la salida

### 3.1 Tabla principal

```text
genome
tf_locus_tag
tf_name
tf_family
tf_class                 # sequence_specific | architectural | antirepressor | other
operator_contig
operator_start
operator_end
operator_strand
operator_sequence
motif_id
motif_score
motif_pvalue
motif_fdr
conservation_score
synteny_score
target_operon
target_genes
target_function
distance_to_tf
orientation              # divergon | same | opposite
effector_candidates      # lista ordenada, con ChEBI
confidence               # high | medium | low
status                   # hypothesis | corroborated | documented | excluded
evidence                 # motivo, conservación, contexto, GBA, estructura
```

### 3.2 Lo que el censo **sí** entrega

- Todos los represores anotados, no una selección.
- Todos los operadores candidatos por represor, no solo el top-1.
- Todos los efectores candidatos por operón, no solo uno.
- Una medida de conservación y contexto para cada fila.
- Una tabla comparable entre genomas.

### 3.3 Lo que el censo **no** entrega

- Una verdad experimental. `status=documented` solo con PMID o estructura del
  complejo.
- Una única diana por represor. Muchos represores tienen varios operadores.
- Una promesa de que el top-1 es el funcional *in vivo*.
- Una dirección inductor/correpresor sin dato experimental.

---

## 4. Determinismo y trazabilidad

El censo es determinista en el sentido computacional:

- mismas entradas + mismos parámetros → mismas salidas;
- versión fijada de genomas, motivos, ortólogos y umbrales;
- hash SHA-256 por genoma y por tabla;
- sin red durante la corrida;
- sin muestreo aleatorio; los controles negativos usan semillas fijas.

No es determinista en el sentido biológico:

- la asignación de función sigue siendo una inferencia;
- la conservación no prueba función;
- el ranking depende de pesos y umbrales.

Por eso la salida lleva `confidence`, `evidence` y `status`.

---

## 5. Límites que hay que declarar

1. **No todos los "represores" son secuencia-específicos.** H-NS, Fis, IHF y
   otros son arquitectónicos; no tienen un operador discreto. Se excluyen del
   censo de operadores y se anotan en otra capa.
2. **Los motivos son degenerados.** Un PWM captura una familia, no una única
   secuencia. Muchos hits son near-sites no funcionales; la conservación ayuda,
   pero no elimina el ruido.
3. **Cross-talk.** Un TF puede unir sitios de otra familia; un operón puede
   tener sitios de varios TF. El censo debe permitir múltiples TF por operón.
4. **Familias huérfanas.** Sin PWM curado ni suficientes ortólogos, el
   descubrimiento comparativo puede fallar. Se reporta `no_motif`.
5. **Transferencia horizontal.** Un TF reciente puede no tener sitios
   conservados; su confianza será baja.
6. **Regulones superpuestos.** Un operón puede estar regulado por varios TF;
   no se fuerza una asignación única.
7. **Validación experimental.** El censo prioriza; no sustituye EMSA/ITC.

---

## 6. Plan de ejecución

### 6.1 Pasada 1 — Cobertura (días)

- Clasificar los 1.555 "represores" en secuencia-específicos / arquitectónicos /
  antirepresores / otros.
- Construir la biblioteca de motivos por familia.
- Escanear todos los promotores de los 400 genomas con los motivos de familia.
- Guardar todos los hits con p-valor y FDR.

**Entregable:** `census_hits.tsv` (todas las filas, todos los genomas).

### 6.2 Pasada 2 — Evolución (semanas)

- Construir grupos de ortólogos.
- Alinear promotores.
- Calcular conservación y sintenia.
- Integrar contexto y función.
- Rankear y asignar confianza.

**Entregable:** `census_ranked.tsv`, `census_report.md`.

### 6.3 Pasada 3 — Efectores (semanas)

- Aplicar GBA d=0 a cada operón diana.
- Añadir prior de familia.
- Excluir cofactores ubicuos.
- Reportar top 10 con `status=hypothesis`.

**Entregable:** `census_effectors.tsv`.

### 6.4 Pasada 4 — Calibración (1–2 semanas)

- K-12 contra RegulonDB: 3.444 sitios.
- *Held-out*: CollecTF/DBTBS/RegPrecise.
- Controles negativos: motivo barajado, PWM de otra familia, región no promotora.
- Métricas: precisión, recall, FDR, enriquecimiento por conservación.

**Entregable:** `census_benchmark.md` con números reales.

### 6.5 Pasada 5 — Opcional, estructura

Solo para los casos de confianza baja o familias sin motivo:

- estructura PDB/AFDB existente;
- DeepPBS/NA-MPNN si está disponible;
- docking solo como apoyo, no como verdad.

No se cristaliza nada. No se predice estructura para todos.

---

## 7. Qué se necesita para empezar

| Necesidad | Estado | Acción |
|---|---|---|
| Genomas y anotaciones | Existe (400) | Usar `data/regulation/` |
| Operones/promotores | Existe | Usar `operons.tsv`, `promoters.tsv` |
| PWM de familia | Parcial (82 de *E. coli*) | Ampliar con CollecTF/DBTBS/RegPrecise |
| Ortólogos | Parcial (39.407 pares de sintenia) | Construir grupos con MMseqs2/DIAMOND |
| Alineamiento de promotores | No existe | Implementar |
| GBA d=0 | Existe (prototipo) | Integrar al censo |
| Esquema de salida | No existe | Definir y validar |
| Tests | No hay para el pipeline 18 | Añadir |

---

## 8. Criterios de aceptación

- [ ] Todos los represores anotados aparecen en la tabla, incluso si su
      resultado es `no_motif` o `excluded`.
- [ ] Todos los hits pasan por un filtro de FDR, no por un top-k arbitrario.
- [ ] Cada fila tiene conservación, contexto, confianza y procedencia.
- [ ] El censo corre en los 400 genomas sin red y con hash reproducible.
- [ ] La calibración en K-12 reporta precisión y recall, no solo recall.
- [ ] Los controles negativos están incluidos.
- [ ] Ningún `status=documented` sin evidencia experimental.
- [ ] El informe declara los límites y los fallos.

---

## 9. Respuesta en una frase

**No hay que cristalizar ni revisar un represor a la vez: hay que enumerar
todos los sitios posibles de todos los represores, puntuarlos con el motivo de
su familia, filtrarlos por conservación evolutiva y contexto génico, y publicar
el censo completo con confianza y procedencia. La estructura existente ayuda,
pero no es el motor; el motor es la evolución leída como conservación.**
