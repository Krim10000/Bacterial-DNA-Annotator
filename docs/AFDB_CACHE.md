# Caché AFDB/PAE

## Qué hay

- `data/metabolism/targets/afdb25k/*.pdb`: modelos AFDB locales.
- `data/metabolism/targets/afdb25k/*.pae.json`: PAE descargado.
- `data/metabolism/targets/afdb_catalog.tsv`: catálogo por accesión.

## Descargador

```bash
python scripts/structure/fetch_afdb_cache.py --limit 25 --pae --workers 3
```

- reúne accesiones únicas de `data/regulation/*/genes.tsv`;
- no descarga dos veces la misma accesión;
- reanudable;
- catálogo con estado, hash, tamaño y URL;
- PAE opcional con `--pae`.
- HTTP 404 se marca `not_found` y no se reintenta en cada lote; HTTP 429 sí.
- Si AFDB responde HTTP 429, esperar y bajar `--workers`/subir `--delay`:
  ```bash
  python scripts/structure/fetch_afdb_cache.py --pae --workers 1 --delay 1.5
  ```
  429 no es un 404: se reintenta; el catálogo lo conserva como error temporal.

## Estrategia de espacio

- una accesión UniProt = una entrada;
- deduplicar por hash de secuencia cuando dos accesiones tengan la misma proteína;
- comprimir PDB con gzip/zstd;
- descargar bajo demanda, no empaquetar todo;
- PAE solo para modelos que se van a evaluar.

## Estimación

- 111.569 accesiones únicas;
- 1.192 modelos iniciales;
- ~0,205 MB por PDB;
- faltantes ~110.377 → ~22,6 GB en PDB sin comprimir;
- con compresión razonable: ~6–8 GB;
- PAE añade espacio aparte.

## Límites

- AFDB no cubre proteínas sin UniProt;
- un modelo AFDB monomérico no es estructura cuaternaria;
- PAE no es evidencia experimental.
