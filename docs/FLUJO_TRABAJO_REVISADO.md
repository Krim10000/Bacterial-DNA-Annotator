# Flujo de trabajo revisado

> **De “un represor → una diana” a “todos los represores → todos los sitios,
> con K-12 como backtest y la evolución como filtro”.**
>
> Este documento reemplaza el flujo del pipeline maestro propuesto. Conserva la
> idea de filtrar por capas, pero corrige la falsa precisión: nada de
> “promotor >9 bits = ON”, “spacer 19–20 = MerR” ni “choque 3D = represor”.
>
> Documentos relacionados:
> [`CENSUS_REPRESSORES.md`](CENSUS_REPRESSORES.md),
> [`K12_BACKTEST.md`](K12_BACKTEST.md),
> [`K12_ATLAS_SUBCASOS.md`](K12_ATLAS_SUBCASOS.md),
> [`SIGMA_AUDIT.md`](SIGMA_AUDIT.md),
> [`REVIEW_PIPELINE_MAESTRO.md`](REVIEW_PIPELINE_MAESTRO.md).

---

## 0. Principios

1. **K-12 primero, como backtest.** Es el único genoma con verdad curada
   suficiente. No se usa para entrenar y reportar sobre lo mismo.
2. **Censo, no caso único.** Todos los represores anotados, todos los sitios
   posibles, todos los genomas.
3. **Secuencia y contexto primero; estructura después.** El 3D es caro y
   ruidoso; se reserva para el top N.
4. **Sin cristalizar.** Se usan PDB/AFDB existentes. La estructura es evidencia
   auxiliar, no requisito.
5. **Nada se descarta en silencio.** Cada fila descartada tiene `reason` y queda
   en `unresolved`/`discarded`.
6. **Hipótesis con confianza.** `status=hypothesis` por defecto;
   `documented` solo con experimento o estructura del complejo.

---

## 1. Inputs

| Input | Fuente | Uso |
|---|---|---|
| Genomas | 400 + K-12 | Censo y conservación |
| Genes/proteínas | `data/regulation/`, `data/annotations/` | Universo de TF |
| Operones | `operons.tsv` | Contexto y dianas |
| Promotores | `promoters.tsv` | Ventanas de escaneo |
| PWMs de TF | `tf_pwms.tsv` (82) | Motivos |
| PWMs de sigma | `promoter_pwms.tsv` (σ24/28/32/38/54/70) | Asignación de sigma |
| RegulonDB | `TF-RISet.tsv`, `TFSet.tsv`, `gt_promoters.tsv` | Verdad K-12 |
| Ortólogos/sintenía | 400 genomas, 39.407 pares | Conservación |
| Estructuras | PDB/AFDB existentes | Evidencia Tier 2 |
| RNA-seq | opcional | Prior de expresión, no verdad |

---

## 2. Fase 0 — Backtest K-12

**Objetivo:** saber qué recupera el pipeline antes de escalar.

```text
K-12 → predicción → comparación contra RegulonDB → reporte
```

### 2.1 Universo

- 243 TF curados.
- 148 represores con sitio S/C.
- 2.119 sitios represor S/C; 1.869 con PWM asignable.
- 82 PWM de TF; 62/148 represores tienen PWM.

### 2.2 Comparación

Cuatro niveles:

| Nivel | Criterio |
|---|---|
| Sitio | solapamiento ≥1 pb, misma hebra |
| Operón | mismo `targetTuOrGene` |
| TF | sitio conocido entre los candidatos del TF |
| Efector | mismo ChEBI/InChIKey |

### 2.3 Clases

| Clase | Significado |
|---|---|
| `recovered` | conocido y predicho, coincide |
| `known_missed` | conocido, no predicho |
| `conflict` | conocido y predicho, distinto |
| `novel` | predicho, no conocido |
| `unresolved` | ni conocido ni predicho |
| `not_applicable` | arquitectónico/antirepresor |
| `no_motif` | sin PWM/familia |

### 2.4 Regla anti-circularidad

- El pipeline se congela antes de comparar.
- Los umbrales se ajustan en una partición por familia (`leave-one-family-out`).
- El reporte final se hace sobre el conjunto completo una vez.

**Entregable:** `data/census/k12/k12_comparison.tsv`, `k12_report.md`.

---

## 3. Fase 1 — Censo Tier 0: secuencia y contexto

**Objetivo:** enumerar todos los candidatos posibles sin 3D.

### 3.1 Universo

Todos los TF anotados, clasificados:

- `sequence_specific`: entra al censo de operadores;
- `architectural` (H-NS, Fis, IHF, HU): `not_applicable`;
- `antirepressor`: `not_applicable`;
- `sigma_factor`: otra capa;
- `other`: `unresolved`.

### 3.2 Biblioteca de motivos

- PWM de TF/familia (RegulonDB, CollecTF, DBTBS, RegPrecise, PRODORIC).
- Descubrimiento comparativo para familias sin PWM: ortólogos + MEME/DREME.
- Si no hay motivo: `no_motif`, no se inventa.

### 3.3 Escaneo

- Todas las regiones reguladoras: promotores, intergénicas, divergón.
- Ambas hebras.
- Todos los hits con p-valor y FDR; no top-k.
- Coordenadas, hebra, secuencia, score.

### 3.4 Sigma

La asignación de sigma es una capa aparte, no un default. Ver
[`SIGMA_AUDIT.md`](SIGMA_AUDIT.md) y el módulo ejecutable
[`SIGMA_MODULE.md`](SIGMA_MODULE.md).

- σ70: −35/−10, espaciador 15–19, −10 extendida, UP.
- σ54: −24/−12, espaciador 10–11, bEBP, IHF.
- σ32/σ38/σ28/σ24: matrices propias.
- Si el margen entre sigmas es bajo: `sigma_status=ambiguous`.
- Si no hay señal: `unresolved`, no `sigma70`.

### 3.5 Contexto

- Operón diana, genes, función.
- Distancia al TSS y al gen del TF.
- Orientación: divergón, mismo sentido, opuesto.
- Coherencia funcional TF–operón.

**Entregable:** `census_hits.tsv`, `census_context.tsv`.

---

## 4. Fase 2 — Tier 1: evolución

**Objetivo:** separar sitios funcionales de *near-sites*.

1. Grupos de ortólogos de los 400 genomas.
2. Alineamiento de promotores.
3. Conservación del hit en posición y orientación.
4. Sintenía de la vecindad.
5. `conservation_score`, `synteny_score`, `phylogenetic_support`.

Para K-12, conservación **leave-one-genome-out**: se calcula con los otros 399,
no con K-12 mismo.

**Entregable:** `census_ranked.tsv` con conservación.

---

## 5. Fase 3 — Tier 2: estructura y RNAP (top N)

**Objetivo:** añadir evidencia geométrica solo donde hace falta.

- Solo top N por genoma (p. ej. 20–50), no todos.
- PDB/AFDB existentes; sin cristalizar.
- Complejo TF–ADN; oligómero si hay plantilla.
- Scores, no etiquetas:
  - `interface_score`
  - `steric_clash_score`
  - `rnap_contact_score`
  - `dna_bend_angle`
  - `mechanism_probability`

La salida es un vector de mecanismos, no `REPRESOR ESTÉRICO`. Si no se puede
distinguir: `mechanism=unknown`.

**Entregable:** `census_structure.tsv` solo para el top N.

---

## 6. Fase 4 — Tier 3: rutas de rescate

Solo si Tier 0–2 no explican:

| Ruta | Qué busca | Requisito |
|---|---|---|
| DNA looping | dos sitios, fase helicoidal, oligómero | geometría + distancia |
| Anti-activación | sitio solapado con activador | mecanismo explícito |
| Sigma alternativa | σ54/σ32/σ38/σ28/σ24 | PWM + arquitectura |
| Riboswitch/atenuador | UTR 5′, Rfam, terminador | Infernal + secuencia |
| Trampa de escape RNAP | sitio que bloquea la salida | evidencia de ocupación |

Cada rescate es un módulo opcional con `evidence`. No es un comodín para forzar
una explicación.

**Entregable:** `census_rescue.tsv`, con `status=hypothesis`.

---

## 7. Fase 5 — Integración y salida

### 7.1 Esquema

```text
genome
tf_locus_tag
tf_name
tf_family
tf_class
operator_contig
operator_start
operator_end
operator_strand
operator_sequence
motif_id
motif_score
motif_pvalue
motif_fdr
sigma
sigma_score
sigma_margin
sigma_status
conservation_score
synteny_score
target_operon
target_genes
target_function
distance_to_tf
orientation
effector_candidates
mechanism
confidence
status
evidence
provenance
```

### 7.2 Reglas

- `status=documented` solo con PMID o estructura del complejo.
- `effector_direction=unknown` por defecto.
- `confidence=high` exige motivo + conservación + contexto, no solo PWM.
- Toda fila se rastrea a un archivo y una versión.
- Hash SHA-256 por genoma y por tabla.

---

## 8. Gates y métricas

### 8.1 Gate K-12

| Métrica | Umbral propuesto |
|---|---|
| TF con al menos un sitio recuperado | ≥70 % de los que tienen PWM |
| Recall de sitio (solapante) | ≥0,60 |
| Precisión de sitio | ≥0,50 |
| Recall de operón diana | ≥0,50 |
| Cobertura de efector | ≥0,50 en los 43 con efector |
| Sigma: recall por sigma | reportar; no agregar |

Los umbrales son propuestos; se congelan antes de correr.

### 8.2 Gate de generalización

- 3 genomas *held-out* con datos de CollecTF/DBTBS/RegPrecise.
- Caída máxima de 20 puntos vs K-12.
- Reporte por familia y por subcaso, no solo global.

### 8.3 Controles negativos

- Motivo barajado.
- PWM de otra familia.
- Región no promotora.
- Ligando *decoy*.
- Sigma equivocado.

---

## 9. Orden de implementación

| Fase | Tarea | Entregable |
|---|---|---|
| P0 | Corregir bug de escaneo y devolver coordenadas | `18_genome_pipeline.py` fix |
| P0 | Backtest K-12 con PWM y operones existentes | `k12_comparison.tsv` |
| P0 | Escáner de sigma por arquitectura | `sigma_promoters.tsv` |
| P1 | Biblioteca de motivos por familia | `motif_library.tsv` |
| P1 | Ortólogos + conservación | `census_ranked.tsv` |
| P2 | 3D top N | `census_structure.tsv` |
| P2 | Rutas de rescate | `census_rescue.tsv` |
| P3 | RNA-seq y lógica | `logic_report.md` |

---

## 10. Qué NO hacer

- No cristalizar proteínas.
- No etiquetar `sigma70` por defecto.
- No definir operones por 45 pb.
- No tratar el 3D como juez binario.
- No inferir lógica AND/OR sin perturbaciones.
- No usar Fisher como anti-alucinación.
- No llamar “descubrimiento” a un `novel` sin conservación y contexto.
- No escalar a 400 sin pasar el gate K-12.

---

## 11. Frase final

**El flujo revisado es: K-12 para calibrar, secuencia y contexto para enumerar,
evolución para filtrar, estructura solo para el top N, y rescates solo cuando
hacen falta. Todo con `evidence`, `confidence`, `status` y sin cristalizar.**
