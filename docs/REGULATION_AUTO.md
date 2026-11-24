# Regulación automática por reglas

`webapp/regulation_layer.py` construye, para un FASTA nuevo:

```text
<genoma>.operons_auto.tsv
<genoma>.promoters_auto.tsv
```

## Operón

- mismo contig y misma hebra;
- CDS ordenados 5'→3';
- hueco ≤ `--operon-max-gap` (default 300 pb).

Etiqueta:

```text
evidence = rule_based
status   = hypothesis_not_experimental
confidence = low/medium
```

No es cotranscripción.

## Promotor

- ventana fija de 50 pb upstream del primer gen;
- consenso sigma70 `-35 TTGACA` y `-10 TATAAT`;
- espaciador 15–19 pb;
- score heurístico; `box_found=yes` solo con coincidencias mínimas.

No es un TSS. No modela sigma alternativo, UP-element, `-10` extendida,
reguladores ni promotores internos.

## Límites

- no usa RegulonDB/DBTBS ni red;
- no hay TSS/TTS;
- no demuestra transcripción;
- `medium` sigue siendo una hipótesis.

## Estado en los 400

Los 400 tienen operones/promotores precalculados. El pipeline nuevo ahora genera
una versión por reglas; no sustituye los datos precalculados ni el benchmark
experimental.
