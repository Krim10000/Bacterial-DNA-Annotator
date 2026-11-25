# Backtest K-12: qué calza con lo que existe y qué no

> **Objetivo:** correr el censo completo sobre *E. coli* K-12 y comparar cada
> predicción contra RegulonDB, para responder tres preguntas:
>
> 1. **¿Qué recupera?** Sitios, operones y efectores ya conocidos.
> 2. **¿Qué falla?** Falsos negativos y conflictos, con causa.
> 3. **¿Qué propone que no está?** Candidatos nuevos, con evidencia y sin
>    llamarlos descubrimientos.
>
> **No es entrenamiento.** El pipeline se congela antes de mirar la comparación.
> **No es el accuracy final.** K-12 es el backtest; los 400 genomas son la
> generalización.

---

## 0. Qué es y qué no es

| Es | No es |
|---|---|
| Validación retrospectiva contra verdad curada | Entrenar con las etiquetas de K-12 |
| Análisis de errores por TF, familia y subcaso | Un número global de accuracy |
| Detección de candidatos nuevos | Prueba de que los candidatos son reales |
| Calibración de umbrales y confianza | Garantía de que funciona en otros filos |

La regla es: **primero se corre, después se compara, después se ajusta.** Si se
ajustan los umbrales mirando la comparación y se reporta el mismo K-12, el
resultado es circular. Para ajustar, se usa *leave-one-family-out* o una
partición de TF; para reportar, se usa el conjunto completo solo una vez.

---

## 1. Insumos

### 1.1 Verdad curada de K-12

| Dato | Valor | Archivo |
|---|---:|---|
| TF curados | 243 | `data/regulondb/TFSet.tsv` |
| Interacciones TF–promotor | 5.785 filas | `data/regulondb/TF-RISet.tsv` |
| Con secuencia de sitio | 4.733 (206 TF) | `TF-RISet.tsv`, columna `10)tfrsSeq` |
| Con secuencia + S/C + evidencia EXP | 3.438 (184 TF) | Filtro de confianza y evidencia |
| Represor + S/C | 2.119 (148 TF) | Filtro `riFunction=repressor` |
| Con coordenadas | 5.785 | `tfrsLeft`, `tfrsRight`, `strand` |
| Con operón diana | 5.785 | `targetTuOrGene` |
| Con función | 5.785 | `riFunction` (repressor/activator/dual) |
| Con confianza | 5.785 | `confidenceLevel` (S/C/W…) |
| TF con efector | 93 | `TFSet.tsv` |
| Represores decisivos | 36 | `data/allostery/` |
| Operones curados | 827 (≥2 genes) | `data/regulondb/gt/gt_operons.tsv` |

### 1.2 Predicción

- Genoma K-12, genes, operones y promotores de `data/regulation/ecoli_k12/`.
- 82 PWM de `data/regulondb/gt/tf_pwms.tsv`.
- Biblioteca de motivos por familia (a construir).
- 400 genomas para conservación, con **leave-one-genome-out** (K-12 excluido del
  cálculo de conservación de sus propios sitios).

### 1.3 Aclaración de conteos

K-12 tiene **243 TF curados**, no 1.555. Los 1.555 son etiquetas de "repressor"
en los 400 genomas. El backtest de K-12 usa 243 TF, de los cuales ~148 tienen
sitio con confianza S/C y ~126 tienen sitio a ≤100 pb del primer gen. No mezclar
los dos universos.

---

## 2. Qué significa "calzar"

Hay cuatro niveles. Un resultado puede calzar en uno y fallar en otro.

| Nivel | Criterio | Qué mide |
|---|---|---|
| **Sitio** | Solapamiento ≥1 pb con el sitio curado, misma hebra | ¿Encontró la región correcta? |
| **Operón** | `target_operon` predicho = `targetTuOrGene` curado | ¿Asignó el gen/operón correcto? |
| **TF** | El sitio conocido está entre los candidatos del TF | ¿Recupera la diana de ese TF? |
| **Efector** | Mismo ChEBI/InChIKey, o mismo compuesto base | ¿Propuso la molécula correcta? |
| **Dirección** | Inductor/correpresor coincide | Solo para los 36 decisivos |

### 2.1 Matching de sitios

- **Exacto:** mismo intervalo (±0 pb) y misma hebra.
- **Solapante:** ≥1 pb de solapamiento y misma hebra.
- **Cercano:** centro a ≤5 pb, misma hebra.
- **Mismo promotor:** mismo `promoterID`/`targetTuOrGene`, aunque el sitio esté
  desplazado.

Reportar los cuatro, no solo el exacto. Un sitio desplazado 2 pb puede ser el
mismo sitio funcional; un sitio en otro promotor no.

### 2.2 Matching de efectores

- Mismo ChEBI o mismo InChIKey (bloque esqueleto).
- Si no hay ChEBI, mismo nombre normalizado.
- No comparar por cadena SMILES: se pierden formas protonadas.
- Reportar `exact`, `same_compound`, `same_pathway`, `different`.

### 2.3 Matching de dirección

Solo para los 36 decisivos. Reportar:

- accuracy;
- balanced accuracy (27 inductor / 9 correpresor);
- recall por clase;
- MCC.

Nunca reportar solo accuracy: predecir siempre "inductor" acierta 75 %.

---

## 3. Pipeline del backtest

```text
K-12 genome + genes + operones + promotores
        │
        ├─ 1. Universo: 243 TF → clasificar → quedarse con represores
        │
        ├─ 2. Motivo: PWM de TF/familia; fallback a familia
        │
        ├─ 3. Escaneo exhaustivo: todos los promotores, ambas hebras
        │       → todos los hits con p-valor y FDR
        │
        ├─ 4. Conservación: ortólogos de los 400, leave-one-genome-out
        │
        ├─ 5. Contexto: operón, orientación, distancia, función
        │
        ├─ 6. Efector: GBA d=0 + prior de familia
        │
        └─ 7. Ranking + confianza + procedencia
```

### 3.1 Universo

- 243 TF de RegulonDB.
- Clasificar: `sequence_specific`, `architectural`, `antirepressor`, `other`.
- Los arquitectónicos y antirepresores se reportan como `not_applicable` al
  censo de operadores; no se cuentan como falsos negativos.
- Los TF sin PWM se reportan como `no_motif`; no se inventa un motivo.

### 3.2 Motivo

- Asignar PWM por nombre de TF; si no hay, por familia.
- Si la familia tiene varias PWM, usar todas y quedarse con la mejor por sitio.
- Si no hay PWM de familia, `no_motif`.
- Registrar `motif_source` y `motif_quality`.

### 3.3 Escaneo

- Regiones: promotores de operones, intergénicas con TSS, upstream del gen del
  TF, regiones divergentes.
- Ambas hebras.
- Todos los hits con FDR ≤ umbral; no top-k.
- Guardar coordenadas, hebra, secuencia, score, p-valor, FDR.

### 3.4 Conservación

- Grupos de ortólogos de los 400 genomas.
- Alineamiento de promotores.
- Para K-12, calcular conservación **excluyendo K-12** (leave-one-genome-out).
- `conservation_score`: fracción de ortólogos con hit alineado.
- `synteny_score`: conservación de la vecindad.
- Si no hay ortólogos suficientes: `conservation=NA`, no cero.

### 3.5 Contexto

- Operón diana, genes, función.
- Distancia al TSS y al gen del TF.
- Orientación: divergón, mismo sentido, opuesto.
- Coherencia funcional TF–operón.
- Arquitectura: tamaño del operón, líder, terminador.

### 3.6 Efector

- GBA d=0 sobre el regulón del TF.
- Prior de familia.
- Excluir cofactores ubicuos.
- Top 10 con `status=hypothesis`.
- No usar docking para todos; solo para casos concretos si hace falta.

---

## 4. Categorías de comparación

Cada fila del reporte cae en una de estas clases:

| Clase | Definición | Interpretación |
|---|---|---|
| `recovered` | Conocido y predicho, coincide | El pipeline calza |
| `known_missed` | Conocido, no predicho | Falso negativo; analizar causa |
| `novel` | Predicho, no conocido | Candidato nuevo o falso positivo |
| `conflict` | Conocido y predicho, pero distinto | Modelo y base discrepan |
| `unresolved` | Ni conocido ni predicho | Sin información |
| `not_applicable` | Arquitectónico/antirepresor | Fuera del censo de operadores |
| `no_motif` | TF sin PWM ni familia | Dato faltante, no fallo del modelo |

### 4.1 Causas de `known_missed`

No todas son errores del modelo. Clasificar:

- `no_pwm`: no hay PWM para el TF/familia.
- `motif_degenerate`: el PWM no captura el sitio.
- `outside_window`: el sitio está fuera de la región escaneada.
- `operon_boundary`: el operón predicho está mal cortado.
- `threshold`: el hit existe pero queda bajo el umbral.
- `conservation_filter`: el hit existe pero no se conserva.
- `context_filter`: el hit existe pero el contexto lo penaliza.
- `annotation_gap`: RegulonDB tiene el sitio pero no la secuencia.
- `bug`: error de código.

### 4.2 Causas de `novel`

- `same_operon_new_site`: sitio nuevo en un operón ya conocido.
- `new_operon`: operón no reportado para ese TF.
- `new_tf`: TF sin sitios reportados.
- `artifact`: hit espurio, repetitivo o de baja complejidad.
- `cross_talk`: sitio de otra familia.

Los `novel` no se publican como hallazgo sin evidencia: conservación, contexto,
motivo de calidad y, si se puede, validación experimental.

---

## 5. Salidas

### 5.1 Archivos

| Archivo | Contenido |
|---|---|
| `data/census/k12/k12_predictions.tsv` | Todos los hits rankeados |
| `data/census/k12/k12_comparison.tsv` | Una fila por (TF, sitio conocido/predicho) con clase y causa |
| `data/census/k12/k12_confusion.json` | Conteos por clase, familia y subcaso |
| `data/census/k12/k12_report.md` | Qué calza, qué no, qué es nuevo |
| `data/census/k12/k12_metrics.json` | Precisión, recall, F1, calibración |

### 5.2 Tabla de ejemplo

| TF | Sitio conocido | Sitio predicho | Match | Efector conocido | Efector predicho | Clase |
|---|---|---|---|---|---|---|
| LacI | O1 | O1 | exacto | alolactosa | alolactosa | `recovered` |
| PurR | PurR box | PurR box | solapante | hipoxantina | hipoxantina | `recovered` |
| TrpR | trpO | — | — | triptófano | triptófano | `known_missed` |
| MarR | marO | marO | exacto | salicilato | salicilato | `recovered` |
| H-NS | — | — | — | — | — | `not_applicable` |

La tabla real tendrá miles de filas; el reporte debe agregarlas por familia y
por subcaso.

---

## 6. Métricas

### 6.1 Sitios

- recall a nivel sitio (solapante y exacto);
- precisión a nivel sitio;
- F1;
- recall@k por TF (k=1, 3, 5, 10);
- fracción de TF con al menos un sitio recuperado.

### 6.2 Operones

- recall de operón diana;
- precisión de operón diana;
- F1;
- exactitud de frontera (si el operón predicho coincide).

### 6.3 Efectores

- recall@1, @3, @5, @10;
- precisión@k;
- enriquecimiento sobre GBA solo;
- cobertura: cuántos TF tienen al menos un efector propuesto.

### 6.4 Dirección

- accuracy;
- balanced accuracy;
- MCC;
- recall por clase (inductor/correpresor).

### 6.5 Por subcaso

Todo lo anterior desagregado por:

- familia;
- clase (`simple`, `multi_operator`, `dual`, `metal`, `xenobiotic`, `orphan`);
- confianza;
- filo (solo cuando se pase a los 400).

**Nunca reportar solo el agregado.** Un buen global puede esconder que todas las
familias grandes fallan.

---

## 7. Qué esperar de la primera corrida

### 7.1 Lo que probablemente calce

- TF con PWM fuerte y sitio canónico: LacI, PurR, TrpR, LexA, Fur, MarR, ArgR.
- Operones bien delimitados: lac, trp, pur, SOS.
- Efectores de familia con química clara: LacI/GalR, TetR, Fur.

### 7.2 Lo que probablemente falle

- TF sin PWM o con PWM degenerado.
- Familias con un solo caso en K-12.
- Sitios fuera de −250 pb.
- Operones largos con fronteras mal cortadas (F1 7+ ≈ 0,44).
- TF arquitectónicos (H-NS, Fis, IHF): no aplica.
- Efectores no metabólicos (MarR, MprA): GBA ciego.
- Dirección: 36 casos, 9 correpresores; no esperar milagros.

### 7.3 Lo que probablemente aparezca como `novel`

- Sitios en operones ya conocidos, desplazados o adicionales.
- Operones no reportados para un TF.
- Sitios de TF sin anotación en RegulonDB.
- Falsos positivos por motivos degenerados.

Los `novel` son la parte más interesante, pero también la más propensa a ruido.
Hay que mirarlos con conservación y contexto, no con entusiasmo.

---

## 8. Plan de implementación

| Fase | Qué se hace | Entregable |
|---|---|---|
| **B0** | Arreglar el bug de escaneo y crear el harness | `scripts/census/00_backtest_k12.py` |
| **B1** | Correr con PWM actuales, sin conservación | `k12_predictions.tsv` v0 |
| **B2** | Comparar contra RegulonDB | `k12_comparison.tsv`, `k12_confusion.json` |
| **B3** | Añadir conservación leave-one-genome-out | v1 con `conservation_score` |
| **B4** | Añadir contexto y efector GBA | v2 con `effector_candidates` |
| **B5** | Reporte y revisión manual de `novel` | `k12_report.md` |
| **B6** | Ajuste por LOFO, no por K-12 completo | umbrales calibrados |

### 8.1 Regla de oro

No ajustar el pipeline con las mismas filas que se reportan. Si se ajusta un
umbral, se hace en una partición por familia y se reporta en la otra.

---

## 9. Criterios de aceptación

- [ ] Los 243 TF aparecen en la tabla, incluso si son `no_motif` o
      `not_applicable`.
- [ ] Todos los sitios conocidos con secuencia se comparan.
- [ ] Cada fila tiene clase, causa, evidencia y confianza.
- [ ] El reporte separa `recovered`, `known_missed`, `novel` y `conflict`.
- [ ] Las métricas se reportan por familia y por subcaso.
- [ ] Los `novel` tienen un nivel de evidencia, no una etiqueta de hallazgo.
- [ ] El backtest es reproducible: hash de entradas y salidas.
- [ ] Ningún `confirmed` sin experimento.

---

## 10. Respuesta en una frase

**Correr el censo completo sobre K-12 y compararlo contra RegulonDB no es para
“ver si acierta” en abstracto: es para separar lo que el modelo recupera, lo que
la base no tiene, lo que el modelo no puede ver por falta de datos, y los
candidatos nuevos que merecen una segunda mirada.**
