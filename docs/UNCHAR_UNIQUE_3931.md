# Deduplicación de 25.058 clusters

## Hallazgo

Los 25.058 clusters de las 54.031 proteínas "uncharacterized" contienen solo:

| métrica | valor |
|---|---:|
| Clusters | 25.058 |
| Accesiones UniProt únicas | **3.931** |
| Con AFDB descargable | 1.192 |
| Sin AFDB (404) | 2.739 |

## Implicación

No hay que correr 25.000 estructuras.

```text
25.058 clusters
  → 3.931 accesiones únicas
  → 1.192 modelos AFDB
  → filtros una vez por accesión
  → propagar a clusters
```

Esto reduce el trabajo estructural de 25k a ~1,2k.

## Cobertura

- 1.192 accesiones tienen modelo AFDB.
- 2.739 no tienen modelo descargable (404).
- Para esas no se puede aplicar F1–F4 sin plegar.

## Regla

- Los filtros se ejecutan por **accesión única**.
- El resultado se propaga a todos los clusters/genomas con esa accesión.
- No se repite F1–F4 por copia genómica.

## Archivos

- `UNCHAR_25K_MANIFEST.tsv`
- `UNCHAR_UNIQUE_3931.tsv`
- `UNCHAR_AFDB_PRESENT.tsv`
