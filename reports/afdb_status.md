# Estado del catálogo AFDB/PAE

Generado: `2026-11-25T16:33:14Z`

Script: `scripts/structure/afdb_status.py` · modo solo lectura.
Catálogo: `data/metabolism/targets/afdb_catalog.tsv` · caché: `data/metabolism/targets/afdb25k`.

Este informe **no descarga**, **no consulta la red** y **no modifica** el catálogo ni los archivos de la caché. Solo lee rutas, tamaños y estados locales.

## Resumen

| métrica | valor |
|---|---:|
| universo de accesiones | 111.569 |
| filas del catálogo | 2.001 |
| filas `ok` | 1.976 |
| filas `rate_limited` (429) | 0 |
| filas `not_found` (404) | 25 |
| filas `error` (otros) | 0 |
| PDB locales | 1.976 |
| PAE locales | 1.764 |
| PDB locales sin fila de catálogo | 0 |
| PDB locales sin PAE | 212 |
| pendientes de PDB | 109.568 |
| cobertura PDB sobre el universo | 1.7711% |
| cobertura contable (PDB + 404) | 1.7935% |

## Inventario del catálogo

| estado | filas | % del catálogo | ¿reintentable? | ¿permanente? | criterio |
|---|---:|---:|:---:|:---:|---|
| `ok` | 1.976 | 98.75% | no | no | fila `ok`/`already_present` con archivo local |
| `rate_limited` | 0 | 0.00% | sí | no | 429 o `Retry-After` agotado |
| `not_found` | 25 | 1.25% | no | sí | 404; no se vuelve a pedir |
| `error` | 0 | 0.00% | sí | no | 5xx, timeout, red u otro 4xx; se reintenta en otro lote |
| **total** | **2.001** | **100.00%** | | | |

> El catálogo guarda **15** filas con `status=error` en bruto. La auditoría las reclasifica por el código HTTP del campo `error`: **0** son 429 (reintentables) y **15** son 404 (permanentes). El catálogo **no** se reescribe; la clasificación es efectiva para el plan.

## Cobertura sobre el universo

| medida | n | % del universo |
|---|---:|---:|
| modelos PDB presentes en disco | 1.976 | 1.7711% |
| filas de catálogo `ok` con archivo | 1.976 | 1.7711% |
| PDB locales + 404 permanentes | 2.001 | 1.7935% |
| pendientes de PDB | 109.568 | 98.2065% |

## Pendientes

| concepto | n |
|---|---:|
| nunca intentadas (sin fila de catálogo y sin PDB) | 109.568 |
| reintentables con fila (429/error, sin PDB) | 0 |
| permanentes 404 (no se reintentan) | 25 |
| **total pendiente de PDB** | **109.568** |

## Caché local

| archivo | n | bytes | cobertura |
|---|---:|---:|---:|
| PDB | 1.976 | 480.415.381 | 1.7711% del universo |
| PAE | 1.764 | 887.985.444 | 89.2713% de los PDB locales |
| PAE huérfanos (sin PDB) | 0 | | |

## Plan de lotes recomendado

Sin descargar nada: solo se calcula el próximo ritmo. Supuestos del ETA: 2 solicitudes por accesión nueva (metadatos + PDB), 0.5 s de latencia media por solicitud, factor de reintento 1.0000 por el 0.00% de 429 observado, y 30 s de pausa entre lotes. No incluye esperas largas de `Retry-After`: el ETA es una estimación optimista.

| paso | workers | delay (s) | limit | PAE | ETA del lote | condición |
|---|---:|---:|---:|:---:|---:|---|
| `produccion` | 2 | 0.5 | 100 | no | 1.7 min | Sin 429 observados en el catálogo. |

- Ritmo sostenido (`produccion`, 2 workers, delay 0.5, limit 100): **39.6 h** para los 109.568 PDB pendientes (~1096 lotes).
- PAE de los 212 PDB locales que no lo tienen: hacerlo **después** del PDB y solo si esos modelos se van a evaluar; con el mismo ritmo son ~8.8 min.

## Siguiente lote recomendado

```bash
python scripts/structure/fetch_afdb_cache.py --limit 100 --workers 2 --delay 0.5 --retries 3
```

Accesiones en el orden natural del descargador (alfabético):

| # | accesión |
|---:|---|
| 1 | `A0A7U4E834` |
| 2 | `A0A7U9P668` |
| 3 | `A0A7V4WV16` |
| 4 | `A0A7W3N5X5` |
| 5 | `A0A7W3RCJ3` |
| 6 | `A0A7W3RDA3` |
| 7 | `A0A829NF98` |
| 8 | `A0A834R821` |
| 9 | `A0A848M4Z0` |
| 10 | `A0A861BA10` |

Hoy las primeras pendientes son las que quedaron en 429; el descargador las reintenta porque solo los 404 son permanentes.

## Criterio de avance

- Después del lote, ejecutar `python scripts/structure/afdb_status.py --dry-run` y contar los nuevos `ok`/`rate_limited`.
- Si >= 9/10 quedan `ok`, pasar al paso `rampa`.
- Si <= 8/10, mantener 1 worker / 2,0 s y repetir 10; **no acelerar**.
- No usar `--pae` durante la recuperación de PDB: duplica solicitudes.

## Artefacto TSV

`reports/afdb_status.tsv` tiene una fila por accesión con evidencia local o de catálogo (2.001 filas), ordenada por estado (reintentables primero). Columnas: `accession`, `status`, `catalog_status`, `local_pdb`, `local_pae`, `pending`, `retryable`, `permanent`, `in_catalog`, `in_universe`, `http_code`, `pdb_path`, `pae_path`, `error`.

Las accesiones nunca intentadas no se listan una por una: se cuentan en el resumen y el siguiente lote se toma del orden determinista del descargador.

## Reproducir

```bash
.venv/bin/python scripts/structure/afdb_status.py --dry-run
.venv/bin/python scripts/structure/afdb_status.py
.venv/bin/python -m pytest tests/test_afdb_status.py -q
```
