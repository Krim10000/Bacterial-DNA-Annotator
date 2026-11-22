# Diccionario de etiquetas y procedencia

Regla del proyecto: **ninguna etiqueta existe sin fuente**. Toda fila de salida
lleva siete campos de procedencia.

| campo | significado |
|---|---|
| `label` | etiqueta de la fila (ver tabla) |
| `label_source` | herramienta/base/versión que la origina |
| `label_evidence` | valores concretos que la sustentan |
| `label_confidence` | `high` / `medium` / `low` |
| `label_interpretation` | qué significa en lenguaje biológico |
| `label_pipeline` | versión del motor (`09_ghost 0.1.0`) |
| `label_date` | fecha UTC de generación |

## Etiquetas

| etiqueta | significado | fuente típica | confianza |
|---|---|---|---|
| `annotated_known` | gen con función específica | DIAMOND/Swiss-Prot, Pfam | alta |
| `predicted_hypothetical` | gen predicho sin función | pyrodigal + Swiss-Prot | media |
| `conserved_hypothetical` | "hypothetical" en ≥2 genomas | Tier 2 + dedup 95% | media |
| `candidate_novel` | ORF no anotado con evidencia de codificación | 09_ghost (codones/RBS/ESM) | media |
| `candidate_conserved` | ORF no anotado en ≥2 genomas | 09_ghost + dedup 95% | alta |
| `desert_region` | región ≥3 kb sin señales | GFF + ORCA + móviles + defensa + ARN | media/baja |
| `desert_candidate` | desierto con ORFs conservados | desierto + dedup | alta |
| `gap_biosynthetic` | brecha enzimática, no transporte, con demanda | Rhea + Swiss-Prot + filo | media/alta |
| `gap_closure_candidate` | candidato que podría cerrar la brecha | Pfam → EC → Rhea | media |
| `gap_closure_confirmed` | cierre validado (dirección verificada) | revisión curada | alta |
| `discarded_artifact` | brecha descartada (ión/transporte/ROS/sin demanda) | refinador | — |

## Niveles de confianza

- **high**: dos o más evidencias independientes coinciden (p. ej. dominio Pfam +
  reacción Rhea + conservación).
- **medium**: una evidencia fuerte (dominio o conservación) sin la otra.
- **low**: solo indicio débil (p. ej. región vacía sin ORFs).

## Ejemplos

```json
{
  "label": "desert_candidate",
  "label_source": "GFF (pyrodigal/Infernal/tRNAscan) + ORCA/DoriC + móviles + defensa + ARN + 09_ghost + conservación",
  "label_evidence": "len=3448; gc=0.424; orfs=5; conserved=1; max_rank=0.436",
  "label_confidence": "high",
  "label_interpretation": "Hueco de anotación con ORFs presentes en ≥2 genomas: candidato a gen/operón no anotado.",
  "label_pipeline": "09_ghost 0.1.0",
  "label_date": "2026-11-19T15:00:00Z"
}
```

```json
{
  "label": "gap_closure_candidate",
  "label_source": "Pfam-A → EC → Rhea + conservación",
  "label_evidence": "pfam=Glyoxalase; EC=4.4.1.5; RHEA:19069; genoma=GCF_000199675.1; ng=3",
  "label_confidence": "medium",
  "label_interpretation": "Glioxalasa conservada en un genoma sin producción de (R)-S-lactoilglutationa; dirección de Rhea reversible, verificar.",
  "label_pipeline": "09_ghost 0.1.0",
  "label_date": "2026-11-19T15:00:00Z"
}
```

## Estados de decisión

La etiqueta dice *qué es*; el estado dice *cómo se decidió* y *qué se intentó*.
Ver `PIPELINE_LOGIC.md` y `ghost/status.py`:

- `validated`, `consensus`, `lineage_specific`, `anomalous`, `unresolved`, `novel`.
- Campos: `status`, `status_validated_by`, `status_reason`, `status_attempted`,
  `status_confidence`, `status_evidence`.
