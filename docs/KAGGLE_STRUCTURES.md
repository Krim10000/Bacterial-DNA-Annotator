# Estructuras ESMFold de Kaggle (`data/Kaggler`)

## Qué contiene

ZIP con:

```text
<lote>/metrics.tsv
<lote>/pdb/<uniprot>.pdb
```

`metrics.tsv` trae:

```text
uniprot status length plddt ptm seconds error
```

Son estructuras **predichas por ESMFold**, no experimentales. Los ZIP se van
añadiendo con el tiempo en `data/Kaggler`; la ingesta está pensada para
reejecutarse sin duplicar ni perder lo anterior.

## Ingesta incremental

```bash
# 1) Ver qué entraría y auditar sin escribir nada
.venv/bin/python scripts/structure/ingest_kaggle_zips.py --dry-run --report --prune-missing

# 2) Ingesta real + auditoría
.venv/bin/python scripts/structure/ingest_kaggle_zips.py --report --prune-missing

# 3) Reporte TSV para revisar aparte (también válido con --dry-run)
.venv/bin/python scripts/structure/ingest_kaggle_zips.py \
    --dry-run --report --prune-missing \
    --report-tsv reports/kaggle_ingest_report.tsv
```

La CLI recorre `data/Kaggler/*.zip` ordenados por nombre. Para cada ZIP imprime
un resumen por lote:

```text
set_006.zip: nuevos=79 duplicados=0 mejorados=0 omitidos=0
```

y al final:

```text
total catálogo=1644 nuevos=79 duplicados=1565 mejorados=0 omitidos=0 errores=0
reporte: zips_nuevos=1 pdb_sin_catalogo=0 zip_incompletos=0 \
accesiones_duplicadas=0 pdb_corruptos=0
prune-missing: pdb_ausentes=0 marcados=0 (no se borró nada)
```

- `nuevos`: accesiones que no estaban en el catálogo (o que estaban en
  `status=missing` y vuelven a tener PDB).
- `duplicados`: PDB de una accesión ya catalogada que **no** reemplazan la
  copia existente porque su pLDDT es menor o igual. Reejecutar la ingesta
  completa produce este contador, pero no cambia el catálogo.
- `mejorados`: la copia nueva tenía mayor pLDDT y reemplazó a la anterior.
- `omitidos`: PDB no ingeridos por `status` de métricas fallido, PDB truncado
  o ZIP incompleto.
- `errores`: ZIP ilegible.

## Catálogo

Ruta: `data/metabolism/targets/esmfold_catalog.tsv`.

| columna | contenido |
|---|---|
| `accession` | UniProt/identificador de la accesión |
| `status` | `ok`, o `missing` si `--prune-missing` detectó el PDB local ausente |
| `plddt` | pLDDT medio tal como lo reporta `metrics.tsv` |
| `ptm` | pTM reportado por `metrics.tsv` |
| `length` | longitud reportada por `metrics.tsv` |
| `source_zip` | ZIP del lote que aportó la copia conservada |
| `source_member` | ruta del PDB dentro del ZIP |
| `sha256` | SHA-256 de los bytes del PDB (no del ZIP) |
| `bytes` | tamaño del PDB |
| `duplicate_of` | ZIP de otra copia desplazada: la reemplazada cuando hubo mejora, o la descartada por peor pLDDT; vacío si no hubo duplicado |
| `error` | motivo de un estado no utilizable |
| `ingested_at` | fecha de la primera ingesta; se conserva en las reejecuciones |

Una accesión aparece **una sola vez**. Si dos lotes traen la misma accesión:

1. se conserva el PDB con mayor pLDDT;
2. la copia descartada se reporta como `duplicate_accession` (nunca se borra
   del ZIP de origen);
3. `duplicate_of` guarda el ZIP desplazado cuando hubo reemplazo.

`--report` cuenta accesiones presentes en **más de un lote**. Releer el mismo
ZIP no cuenta como duplicado entre lotes (es idempotencia).

## ZIP incompletos o corruptos

La CLI es estricta por defecto. Un ZIP se reporta y **no se ingiere** si:

- no tiene `metrics.tsv`;
- no tiene PDB;
- un PDB está truncado: vacío, con bytes NUL, sin registros
  `ATOM`/`HETATM`, con líneas `ATOM` incompletas o sin registro `END` final;
- un PDB no tiene fila en `metrics.tsv`.

El reporte incluye `incomplete_zip`, `corrupt_pdb` y/o `missing_metrics`, y el
proceso termina con código 1. Se descarta el ZIP completo para no dejar una
ingesta a medias: un PDB truncado bloquea el resto del lote hasta corregirlo.
Un `status` de métricas distinto de `done`/`ok`/`already_present` solo omite
ese PDB (`metrics_status`); no marca el ZIP como incompleto.

La función `ingest()` conserva `require_complete=False` por compatibilidad con
llamadas históricas que ingerían lotes sin `metrics.tsv` (por ejemplo, el
contrato de `tests/test_structure_catalog.py`). La CLI pasa
`require_complete=True`. `--allow-incomplete` reproduce el modo permisivo y
**no** debe usarse para reingesta nueva.

## Auditoría y `--prune-missing`

`--report` (o `--prune-missing`, que lo implica) detecta, sin borrar nada:

| `kind` | qué detecta |
|---|---|
| `new_zip` | ZIP sin filas previas en el catálogo |
| `incomplete_zip` | ZIP sin `metrics.tsv`, sin PDB o descartado por PDB truncado |
| `corrupt_pdb` | PDB vacío, truncado o ilegible dentro del ZIP |
| `missing_metrics` | PDB sin fila en `metrics.tsv` |
| `metrics_status` | fila de métricas con `status` fallido |
| `duplicate_accession` | misma accesión en más de un lote |
| `duplicate_hash` | accesiones distintas con PDB de bytes idénticos |
| `orphan_pdb` | PDB en `esmfold_kaggle/` sin fila en el catálogo |
| `missing_file` | fila del catálogo cuyo PDB local no existe (solo con `--prune-missing`) |

`--prune-missing` **no poda ni borra**: reporta el faltante y, fuera de
`--dry-run`, cambia `status` a `missing` y conserva `plddt`, `ptm`, `length`,
`sha256`, `source_zip` y `source_member`. Si más adelante reaparece el PDB en
un ZIP, la fila se restaura. Los PDB huérfanos tampoco se borran.

`--report-tsv RUTA` escribe el reporte con las columnas:

```text
kind  accession  source_zip  path  sha256  bytes  related  detail
```

Si se indica explícitamente, el TSV se escribe incluso con `--dry-run`; es un
diagnóstico solicitado. El catálogo y los PDB nunca se modifican en seco.

## Idempotencia

- filas ordenadas por accesión;
- escritura atómica (`*.tmp` + `os.replace`);
- `ingested_at` original conservado;
- un PDB ya catalogado con pLDDT mayor o igual no se reescribe;
- correr dos veces con las mismas entradas produce el mismo catálogo;
- repetir el dry-run produce los mismos conteos.

## Estado real (2026-11-25)

`data/Kaggler` tenía 16 ZIP, 1.644 PDB y 1.565 accesiones ya catalogadas.
El dry-run real de `--report --prune-missing` dio:

| métrica | valor |
|---|---:|
| ZIP | 16 |
| catálogo en disco | 1.565 |
| catálogo tras ingesta simulada | 1.644 |
| nuevos | 79 (`set_006.zip`) |
| duplicados (ya catalogados) | 1.565 |
| mejorados | 0 |
| omitidos | 0 |
| ZIP nuevos | 1 |
| ZIP incompletos | 0 |
| PDB corruptos | 0 |
| PDB sin catálogo | 0 |
| accesiones duplicadas entre lotes | 0 |
| PDB locales ausentes | 0 |

El dry-run no escribió el catálogo ni los PDB: en disco siguen 1.565 filas y
1.565 PDB. Los números cambian con cada ZIP nuevo; el dry-run es la fuente de
verdad.

## Uso en la GUI

La ficha del gen busca, en orden:

1. `pdb` declarado;
2. AFDB en `afdb25k/`;
3. ESMFold en `esmfold_kaggle/`.

Si hay PDB, muestra descarga y visor 3D. Si no hay PAE, lo dice.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_kaggle_refresh.py tests/test_ingest_kaggle.py -q
```

Cubren ZIP nuevo, duplicado entre lotes, hash duplicado, PDB truncado, ZIP sin
`metrics.tsv`, PDB sin fila de métricas, `status` fallido, auditoría,
`--prune-missing` y el dry-run de la CLI. Usan `tmp_path`; no tocan
`data/metabolism/targets`.

## Límites

- ESMFold no da PAE.
- pLDDT/pTM bajos no significan que la proteína no exista.
- Una estructura predicha no confirma función, interfaz ni ensamblaje.
- No es estructura cuaternaria.
- El detector de truncamiento es heurístico (registros `ATOM`, ancho mínimo y
  `END` final); un PDB mutilado que conserve `END` podría pasar. El CRC del
  ZIP detecta archivos dañados, no contenido semánticamente incompleto.
- `duplicate_of` registra solo el último ZIP desplazado; el historial completo
  está en el reporte (`duplicate_accession`).
- El modo estricto descarta el ZIP completo: un solo PDB malo bloquea el lote.
  Es deliberado para no dejar ingestas parciales.
- `--prune-missing` comprueba el PDB local de `esmfold_kaggle/`; no consulta
  fuentes remotas ni arregla el PDB.
- La detección de duplicados entre lotes se apoya en `source_zip` y
  `source_member`; filas externas con `source_member` vacío pueden producir
  falsos duplicados si el mismo ZIP aporta el mismo PDB con otra ruta.
- `sha256` es del PDB descomprimido, no del ZIP.
