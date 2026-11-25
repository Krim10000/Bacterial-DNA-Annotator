# Estado del catálogo AFDB/PAE

Generado: `2026-11-25T12:52:13Z`

Script: `scripts/structure/afdb_status.py` · modo solo lectura.
Catálogo: `data/metabolism/targets/afdb_catalog.tsv` · caché: `data/metabolism/targets/afdb25k`.

Este informe **no descarga**, **no consulta la red** y **no modifica** el catálogo ni los archivos de la caché. Solo lee rutas, tamaños y estados locales.

## Resumen

| métrica | valor |
|---|---:|
| universo de accesiones | 111.569 |
| filas del catálogo | 620 |
| filas `ok` | 478 |
| filas `rate_limited` (429) | 127 |
| filas `not_found` (404) | 15 |
| filas `error` (otros) | 0 |
| PDB locales | 1.670 |
| PAE locales | 478 |
| PDB locales sin fila de catálogo | 1.188 |
| PDB locales sin PAE | 1.192 |
| pendientes de PDB | 109.884 |
| cobertura PDB sobre el universo | 1.4968% |
| cobertura contable (PDB + 404) | 1.5103% |

## Inventario del catálogo

| estado | filas | % del catálogo | ¿reintentable? | ¿permanente? | criterio |
|---|---:|---:|:---:|:---:|---|
| `ok` | 478 | 77.10% | no | no | fila `ok`/`already_present` con archivo local |
| `rate_limited` | 127 | 20.48% | sí | no | 429 o `Retry-After` agotado |
| `not_found` | 15 | 2.42% | no | sí | 404; no se vuelve a pedir |
| `error` | 0 | 0.00% | sí | no | 5xx, timeout, red u otro 4xx; se reintenta en otro lote |
| **total** | **620** | **100.00%** | | | |

> El catálogo guarda **142** filas con `status=error` en bruto. La auditoría las reclasifica por el código HTTP del campo `error`: **127** son 429 (reintentables) y **15** son 404 (permanentes). El catálogo **no** se reescribe; la clasificación es efectiva para el plan.

## Cobertura sobre el universo

| medida | n | % del universo |
|---|---:|---:|
| modelos PDB presentes en disco | 1.670 | 1.4968% |
| filas de catálogo `ok` con archivo | 478 | 0.4284% |
| PDB locales + 404 permanentes | 1.685 | 1.5103% |
| pendientes de PDB | 109.884 | 98.4897% |

## Pendientes

| concepto | n |
|---|---:|
| nunca intentadas (sin fila de catálogo y sin PDB) | 109.761 |
| reintentables con fila (429/error, sin PDB) | 123 |
| permanentes 404 (no se reintentan) | 15 |
| **total pendiente de PDB** | **109.884** |

## Caché local

| archivo | n | bytes | cobertura |
|---|---:|---:|---:|
| PDB | 1.670 | 388.697.119 | 1.4968% del universo |
| PAE | 478 | 346.346.286 | 28.6228% de los PDB locales |
| PAE huérfanos (sin PDB) | 0 | | |

## Plan de lotes recomendado

Sin descargar nada: solo se calcula el próximo ritmo. Supuestos del ETA: 2 solicitudes por accesión nueva (metadatos + PDB), 0.5 s de latencia media por solicitud, factor de reintento 1.2048 por el 20.48% de 429 observado, y 30 s de pausa entre lotes. No incluye esperas largas de `Retry-After`: el ETA es una estimación optimista.

| paso | workers | delay (s) | limit | PAE | ETA del lote | condición |
|---|---:|---:|---:|:---:|---:|---|
| `sonda` | 1 | 2 | 10 | no | 60 s | Siempre. Mide si la API se recuperó de los 429. |
| `rampa` | 1 | 1.5 | 25 | no | 2.0 min | Solo si la sonda deja <= 1/10 rate_limited. |
| `estable` | 2 | 1 | 50 | no | 1.5 min | Ritmo sostenido si el 429 baja a <= 5%. |
| `produccion` | 2 | 0.75 | 100 | no | 2.5 min | Solo si el 429 es <= 2% en el paso anterior. |

- Ritmo sostenido (`estable`, 2 workers, delay 1, limit 50): **3.1 días** para los 109.884 PDB pendientes (~2198 lotes).
- Si el 429 cae a <= 2%, el paso `produccion` baja el total a **2.3 días**.
- Si el 429 persiste, quedarse en `sonda`: **11.5 días**. A ese ritmo la API no admite cerrar la deuda completa; hay que esperar.
- PAE de los 1.192 PDB locales que no lo tienen: hacerlo **después** del PDB y solo si esos modelos se van a evaluar; con el mismo ritmo son ~48.1 min.

## Siguiente lote recomendado

```bash
python scripts/structure/fetch_afdb_cache.py --limit 10 --workers 1 --delay 2 --retries 3
```

Accesiones en el orden natural del descargador (alfabético):

| # | accesión |
|---:|---|
| 1 | `A0A0H3JQ59` |
| 2 | `A0A0H3JR16` |
| 3 | `A0A0H3JRU9` |
| 4 | `A0A0H3JT43` |
| 5 | `A0A0H3JT80` |
| 6 | `A0A0H3JTK0` |
| 7 | `A0A0H3JTL0` |
| 8 | `A0A0H3JU78` |
| 9 | `A0A0H3JUG6` |
| 10 | `A0A0K8P6T7` |

Hoy las primeras pendientes son las que quedaron en 429; el descargador las reintenta porque solo los 404 son permanentes.

## Criterio de avance

- Después del lote, ejecutar `python scripts/structure/afdb_status.py --dry-run` y contar los nuevos `ok`/`rate_limited`.
- Si >= 9/10 quedan `ok`, pasar al paso `rampa`.
- Si <= 8/10, mantener 1 worker / 2,0 s y repetir 10; **no acelerar**.
- No usar `--pae` durante la recuperación de PDB: duplica solicitudes.

## Artefacto TSV

`reports/afdb_status.tsv` tiene una fila por accesión con evidencia local o de catálogo (1.808 filas), ordenada por estado (reintentables primero). Columnas: `accession`, `status`, `catalog_status`, `local_pdb`, `local_pae`, `pending`, `retryable`, `permanent`, `in_catalog`, `in_universe`, `http_code`, `pdb_path`, `pae_path`, `error`.

Las accesiones nunca intentadas no se listan una por una: se cuentan en el resumen y el siguiente lote se toma del orden determinista del descargador.

## Reproducir

```bash
.venv/bin/python scripts/structure/afdb_status.py --dry-run
.venv/bin/python scripts/structure/afdb_status.py
.venv/bin/python -m pytest tests/test_afdb_status.py -q
```
