# Plan de trabajo — anotación de ADN y completado de rutas metabólicas

**Proyecto:** DNA_TOKEN_BAC (`dna_poc`)
**Objetivo:** encontrar, en el ADN no anotado (o anotado de forma imprecisa), los
genes que faltan para completar las rutas metabólicas de 400 bacterias.

---

## 0. Resumen ejecutivo

Ya existe un motor que genera y prioriza ORFs en regiones no anotadas, validado
con 99% de recall y 0% de falsos positivos. Lo que falta es el tramo final:
**agrupar los candidatos de los 400 genomas, evaluarlos barato y cerrar las
brechas metabólicas** (candidato → enzima faltante → evidencia).

El plan tiene 6 fases. La Fase 1 se puede lanzar hoy: corre el motor sobre los
400 genomas (≈6 h), deduplica los candidatos y produce el catálogo global.

---

## 1. Visión del sistema completo

```
ADN (400 genomas)
 │
 ├─ Anotación estructural (ya existe)
 │    pyrodigal (CDS) · Infernal (rRNA) · tRNAscan (tRNA)
 │    ORCA/DoriC (oriC) · transposones · defensa · BGCs
 │
 ├─ Anotación funcional (ya existe)
 │    DIAMOND vs Swiss-Prot → functional.tsv (51,6% con hit)
 │    (falta) Pfam/HMM → dominios
 │
 ├─ Red metabólica (ya existe)
 │    proteína → UniProt → Rhea → reacciones → brechas
 │    (nuevo) filtro de conservación por filo → brechas reales
 │
 ├─ ADN NO anotado → ORFs (ya existe, 09_ghost)
 │    6 marcos · ambigüedad IUPAC/N · RBS · potencial codificante
 │
 ├─ Deduplicación entre 400 genomas (FALTA)
 │    diamond cluster 95% → clusters + nº de genomas
 │
 ├─ Evaluación barata (parcial)
 │    DB de proteínas conocidas · Pfam · "sentido biológico" en CPU
 │
 ├─ Cierre de brechas (FALTA el empalme)
 │    candidato ↔ enzima que produce el metabolito faltante
 │
 └─ Estructura solo para el top (parcial)
      ESMFold → Foldseek (TEM-1 validado al 99,6%)
```

---

## 2. Inventario: lo que ya tenemos

### Datos
| recurso | ruta | estado |
|---|---|---|
| Genomas | `data/raw` + `data/raw_excluded` | 400 FASTA |
| Anotación ab initio | `data/annotations_ab_initio/<GCF>/` | 400, 1.337.745 CDS |
| Anotación funcional | `functional.tsv` por genoma | 400, 690.157 hits |
| Redes metabólicas | `data/metabolism/*_network.json` | 400 |
| Brechas / trazas | `*_network_trace.json` | 400 |
| Tokenizador + modelo ADN | `data/tokenizer`, checkpoints | 16k BPE, 97M params |
| DB Swiss-Prot | `data/db/swissprot.dmnd` | 281 MB |
| DB estructuras | `data/db/foldseek_swissprot` | 3 GB |
| DB proteomas | `data/db/proteomes.dmnd` | 492 MB (400 genomas) |

### Herramientas
`pyrodigal` · `Infernal/cmsearch` · `tRNAscan-SE` · `ORCA` · `DoriC` ·
`DIAMOND 2.1.11` · `HMMER` · `Foldseek` · `torch CPU` · `fair-esm` ·
`sentencepiece`.

### Motor nuevo (`scripts/ghost/`, 15 módulos)
| módulo | función |
|---|---|
| `seqio.py`, `orfs.py` | regiones, 6 marcos, ambigüedad |
| `scoring.py` | codones, RBS, composición, ranking |
| `pipeline.py`, `run.py` | orquestación y CLI |
| `protein_lm.py` | ESM-2 8M en CPU |
| `validate.py` | control positivo/negativo |
| `gap_demand.py`, `biomass.py`, `gap_conservation.py` | priorización de brechas |
| `gap_link.py`, `gap_targets.py` | enzimas objetivo y cierre |
| `structure.py` | ESMFold + Foldseek |
| `report.py` | reporte Markdown |

---

## 3. El hueco: "completar cadenas metabólicas"

Se ataca por dos lados que hay que empalmar.

### Lado A — Candidatos (ADN no anotado)
Tres niveles, de más barato a más caro:

| tier | qué es | tamaño | cómo se obtiene |
|---|---|---|---|
| **T1** | ORFs en regiones intergénicas | ~250.000 en 400 genomas | ya implementado |
| **T2** | proteínas anotadas como hypothetical/uncharacterized | **53.596** | ya existen en `functional.tsv`; solo hay que evaluarlas |
| **T3** | anotación imprecisa: genes truncados, marco/start alternativo | a medir | re-enumerar ORFs en la ventana del gen ±300 nt |

T2 es el hallazgo más rentable: no hay que leer ADN, ya hay 53.596 cadenas y
varias pueden ser la enzima faltante mal nombrada.

### Lado B — Objetivos (enzimas faltantes)
| paso | estado |
|---|---|
| brechas crudas | 26.308 "genuine" (ruido) |
| filtro de conservación por filo | **11/43** y **6/31** objetivos reales |
| enzimas que podrían cerrarlas | `gap_targets` (Rhea → Swiss-Prot) |

### El empalme que falta
`candidato ↔ enzima objetivo`. Hoy `gap_targets` solo mira el mejor hit
Swiss-Prot de cada ORF y da 0 cierres. Hay que ampliar la evidencia:
homología directa, dominio Pfam, similitud de embedding y (al final) estructura.

---

## 4. Arquitectura propuesta

```
scripts/ghost/ghost/
  # ya existe
  seqio.py  orfs.py  scoring.py  pipeline.py  run.py
  protein_lm.py  validate.py  report.py  structure.py
  gap_demand.py  biomass.py  gap_conservation.py  gap_link.py  gap_targets.py

  # NUEVO
  dedup.py        # [HECHO] diamond cluster 95% sobre los ORFs de los 400 genomas
  tiers.py        # T2 (hypothetical) y T3 (re-enumerar ventanas imprecisas)
  eval_known.py   # DIAMOND vs Swiss-Prot + hmmscan vs Pfam (descarga Pfam)
  bio_sense.py    # "sentido biológico" en CPU: ESM-2 8M o n-gramas AA
  gap_closure.py  # empalme candidato ↔ enzima objetivo, con confianza
  global_report.py# catálogo y reporte de los 400
```

### `dedup.py`
- `diamond cluster --approx-id 90` sobre el FASTA global de candidatos.
- Salida: `cluster_id`, representante, nº de genomas, nº de copias.
- Reduce ~250k ORFs a unos pocos miles de clusters → abarata todo lo demás.

### `tiers.py`
- T2: lee `functional.tsv`, toma las 53.596 hypothetical; las une al catálogo.
- T3: para genes con `partial` o sin RBS, re-enumera la ventana ±300 nt y
  genera cadenas alternativas (extensión, truncamiento, marco desplazado).

### `eval_known.py`
- DIAMOND vs Swiss-Prot para todos los clusters.
- `hmmscan` vs Pfam-A (descarga ~2 GB) para asignar familia/dominio.
- Salida por cluster: hit, dominio, e-value, cobertura.

### `bio_sense.py`
Dos opciones, ambas en CPU:
1. **ESM-2 8M** (recomendada): 0,3 s/péptido; con dedup son ~30–50 min para
   todos los clusters. AUC 0,987 separando proteínas de secuencias barajadas.
2. **Sustituto sin ESM**: modelo de n-gramas de aminoácidos entrenado con las
   1,34 M proteínas del proyecto + composición + Pfam. Más débil, cero
   dependencias.

> Aclaración: ESM-2 8M **no requiere GPU**. Lo que requiere GPU es ESMFold
> (plegamiento). El costo en CPU se controla deduplicando antes.

### `gap_closure.py`
Para cada brecha priorizada y su conjunto de enzimas objetivo:
- ¿algún cluster tiene hit directo a una enzima objetivo? (homología)
- ¿algún cluster tiene el dominio Pfam de esa familia? (dominio)
- ¿algún cluster es vecino en el embedding de una enzima objetivo? (remoto)
- (fase 5) ¿la estructura pliega como esa familia? (Foldseek)

Salida: `gap → candidato → evidencia → confianza (alta/media/baja)`.

---

## 5. Fases, entregables y criterios de éxito

### Fase 1 — Catálogo global (HECHO para 310 de vida libre)
- Correr T1 en los 400 genomas (`run_batch.sh`).
- Deduplicar con `diamond cluster`.
- **Entregable:** `out/global/candidate_clusters.tsv` (representante, genomas,
  hit, rank).
- **Éxito:** ≥90% de los genomas procesados y clusters estables.

### Fase 2 — T2 y T3 (1–2 sesiones) — T2 HECHO · T3 IMPLEMENTADO (lote completo pendiente)
- **T2 hecho:** 53.650 candidatas, 1.164 clusters multi-genoma
  (`out/tier2/CONSERVED.md`).
- Generar cadenas T3 en ventanas imprecisas.
- **Entregable:** catálogo unificado T1+T2+T3.
- **Éxito:** T2 aporta candidatos con rank comparable a T1.

### Fase 3 — Evaluación barata (1 sesión + cómputo)
- Descargar Pfam-A, correr `hmmscan` sobre clusters.
- Puntuar todos los clusters con ESM-2 8M (o el sustituto).
- **Entregable:** `cluster_scores.tsv` con hit, dominio, ESM, conservación.
- **Éxito:** separación clara entre representantes con dominio y sin él.

### Fase 4 — Cierre de brechas (EN CURSO)
- Modelo de brechas refinado (8.655) y primer cierre con evidencia de dominio.
- `gap_closure` sobre las brechas priorizadas de los 400.
- **Entregable:** `gap_closures.tsv` y reporte por genoma.
- **Éxito:** al menos un cierre con evidencia media/alta en el corpus; si no,
  cuantificar cuántas brechas quedan sin explicación homológica.

### Fase 5 — Estructura y validación (1 semana)
- Desiertos de anotación integrados (CDS, RNA, oriC, móviles, defensa).
- Plegar el top de cada genoma (ESMFold local o API cuando vuelva).
- Foldseek contra Swiss-Prot; revisar el top 50 a mano.
- **Entregable:** tabla final "enzima faltante → candidato → estructura".
- **Éxito:** precisión revisada ≥50% en el top 20.

### Fase 6 — Integración (3–4 sesiones) — EN CURSO
- **GUI Streamlit** para subir FASTA y marcar anotaciones (`webapp/`).
- Reporte global + vista en `app_dashboard.py`.
- **Entregable:** dashboard con brechas, candidatos y evidencia.

---

## 6. Decisiones que necesito de ti

1. **"Sentido biológico":** ¿usamos ESM-2 8M en CPU (ya validado) o prefieres
   solo Pfam + n-gramas, sin ESM?
2. **Pfam:** ¿autorizas descargar Pfam-A (~2 GB) para asignar dominios?
3. **Dedup:** umbral **95%** (elegido con el experimento de `dedup_eval`).
4. **Lote 400:** ¿lo lanzo ya en segundo plano?
5. **ESMFold:** ¿esperamos a la API o montamos la versión local?

---

## 7. Riesgos

| riesgo | mitigación |
|---|---|
| Brechas de Rhea incompletas | filtro de conservación por filo; reportar cobertura |
| ORFs cortos con hits espurios | exigir cobertura, filo y longitud; Pfam |
| Costo CPU de ESM | deduplicar antes; solo representantes |
| Anotaciones heterogéneas entre genomas | normalizar por versión y genoma |
| ESMFold caído | Foldseek + homología remota por embedding como puente |

---

## 7.1 Etiquetado con procedencia

Toda salida lleva `label`, `label_source`, `label_evidence`,
`label_confidence`, `label_interpretation`, `label_pipeline` y `label_date`.
Ver `LABELS.md`. Ninguna etiqueta existe sin fuente.

## 8. Primer sprint (concreto)

1. Lanzar `run_batch.sh` sobre los 400 (≈6 h, reanudable).
2. Mientras corre: implementar `dedup.py` y `tiers.py`.
3. Al terminar: deduplicar y correr `eval_known.py` (Swiss-Prot).
4. Entregar `candidate_clusters.tsv` + top 100 con evidencia.

**Criterio de aceptación del sprint:** catálogo global deduplicado, con hit
conocido o dominio para cada cluster, y top 100 revisable.
