# Cinco filtros estructurales

Objetivo: **reducir candidatos** antes de docking/validación. No cierran brechas.

## F1 — Foldseek

- Busca homología estructural en PDB100/AFDB.
- Fuente: `Foldseek`.
- Pasa si hay hits.
- No confirma función.

## F2 — EnzyMM/PyJess

- Busca plantillas catalíticas de M-CSA.
- Fuente: `EnzyMM/PyJess + M-CSA`.
- Pasa si hay hit.
- `no_template` no es ausencia.

## F3 — P2Rank

- Predice bolsillos.
- Pasa si hay bolsillo.
- Se cruza con residuos catalíticos: `pocket_active_hits`.
- Bolsillo ≠ sitio activo.

## F4 — AlphaFill

- Ligandos/cofactores transplantados.
- Pasa si hay ligando con clash ≤20.
- No identifica sustrato fisiológico.

## F5 — Sitio activo local

- Residuos `Active site` conservados del mapa UniProt/estructura.
- Pasa si hay al menos uno.
- Un residuo no demuestra reacción.

## Salidas

`GENOMA.structural_filters.tsv`:

- `f1_foldseek`, `f2_enzymm`, `f3_pocket`, `f4_alphafill`, `f5_active_site`
- `n_pass`
- `priority`: `high`, `review`, `low`
- `source`
- `limitation`

## Regla de prioridad

- `high`: F2 y F5.
- `review`: ≥2 filtros.
- `low`: 0–1.

No se elimina automáticamente un candidato.

## GUI

En el cierre ESMFold:

**Ejecutar 5 filtros**

Muestra tabla, prioridad y créditos.

## Pipeline

```bash
python -m webapp.pipeline genoma.fna \
  --outdir results/ \
  --esmfold "neuros nap" \
  --esmfold-gaps 3 \
  --structure-evidence \
  --structural-annotation \
  --structural-filters
```

## Coste

- Foldseek: $0, local/API.
- EnzyMM: $0, local, ~30–40 s/PDB.
- P2Rank: $0, local Java, ~30 s/PDB.
- AlphaFill: $0, API.
- Sitio activo local: $0, instantáneo.
- No GPU.
- No créditos Neurosnap.

## Límites

- Filtros computacionales, no experimentales.
- `no_template`, `no_pocket` y `no_ligand` no son negativos.
- La procedencia se guarda en cada fila.

## Benchmark

- 10 estructuras: **213,8 s** (~3,6 min).
- Promedio: **21,4 s/estructura**.
- 647.588 proteínas: **~3.850 h** en 1 núcleo; **~120 h** en 32.
- **No se lanza el barrido completo.**
- Detalle: `docs/FIVE_STRUCTURAL_FILTERS_BENCHMARK.md`.

## Hipótesis revisables

Los 5 filtros alimentan:

- `GENOMA.structural_hypotheses.tsv/.json`;
- pestaña GUI **Hipótesis estructural**;
- decisiones `review_decisions.tsv`.

Detalle: `docs/STRUCTURAL_HYPOTHESES_INTEGRATION.md`.
