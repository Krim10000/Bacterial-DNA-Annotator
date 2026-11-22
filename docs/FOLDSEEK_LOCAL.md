# Foldseek local por lotes

## Por qué

La API de Foldseek devuelve HTTP 429 con lotes grandes. La base local elimina esa dependencia.

Base disponible:

```text
data/db/foldseek_swissprot
```

- ~2,4 GB.
- 590.183 estructuras.
- Sin API.

## Consulta por lotes

```bash
foldseek createdb lote_pdb/ lote_db --threads 8

foldseek search lote_db data/db/foldseek_swissprot \
  resultados_db tmp --threads 8

foldseek convertalis lote_db data/db/foldseek_swissprot \
  resultados_db lote.m8 \
  --format-output "query,target,fident,alnlen,evalue,bits"
```

## Uso en los filtros

El resultado `.m8` se pasa como `foldseek_hits` a `run_filters`:

```python
run_filters(pdb, genome, uniprot, foldseek_hits=hits)
```

Así F1 no consulta la API y F2–F5 usan sus cachés.

## Coste observado

- 1.000 PDB de consulta: ~0,21 GB.
- Base: ~2,4 GB.
- `createdb`: segundos.
- `search` contra 590k: ~15 min.
- Filtros con cache: ~7 min.

## Regla

- No parallelizar contra la API.
- Paralelizar F2–F4.
- F1 local por lotes.
- No guardar 647k PDB a la vez.

## Pipeline

```bash
python -m webapp.pipeline genoma.fna \
  --outdir results/ \
  --structure-manifest manifest.tsv \
  --foldseek-local-db data/db/foldseek_swissprot \
  --structural-filters
```

El pipeline:

1. crea una DB de consulta del lote;
2. busca contra la base local;
3. pasa los hits a `run_filters`;
4. escribe `*.structural_hypotheses.tsv/.json`.

No consulta la API.
