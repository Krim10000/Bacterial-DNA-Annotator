# Catálogo unificado de estructuras AFDB + ESMFold

## Qué resuelve

Hasta ahora había dos catálogos con esquemas distintos:

- `data/metabolism/targets/afdb_catalog.tsv` (descarga AFDB);
- `data/metabolism/targets/esmfold_catalog.tsv` (ingesta de ZIP Kaggle).

Este frente agrega un catálogo único, deduplicado y con procedencia, sin
borrar ni mover los PDB originales. La ingesta de ZIP se reutiliza mediante
`scripts/structure/ingest_kaggle_zips.py`.

## Archivos

| archivo | papel |
|---|---|
| `scripts/structure/ingest_kaggle_zips.py` | ingesta de ZIP ESMFold; expone `read_tsv_rows`/`write_tsv_rows` reutilizados por el unificador |
| `scripts/structure/unify_structure_catalog.py` | construye catálogo unificado + reportes; `--dry-run`; sin red |
| `tests/test_structure_catalog.py` | catálogos sintéticos y ZIP sintético |
| `docs/STRUCTURE_CATALOG.md` | este documento |

## Entradas y salidas

Entradas por defecto:

```text
data/metabolism/targets/afdb_catalog.tsv
data/metabolism/targets/esmfold_catalog.tsv
data/metabolism/targets/afdb25k/*.pdb
data/metabolism/targets/esmfold_kaggle/*.pdb
```

Salidas por defecto:

```text
data/metabolism/targets/structure_catalog.tsv
reports/structure_catalog_report.tsv
reports/structure_catalog_report.md
```

Ninguna salida se escribe con `--dry-run`.

## Esquema del catálogo

Una fila por accesión, ordenada por `accession`:

| columna | tipo | significado |
|---|---|---|
| `accession` | texto | accesión UniProt; clave única del catálogo |
| `source` | `AFDB` \| `ESMFold` | fuente del modelo conservado |
| `plddt` | decimal 0–1 o vacío | ESMFold: `metrics.tsv`; AFDB: media de B-factor de CA del primer modelo / 100 |
| `ptm` | decimal 0–1 o vacío | solo ESMFold; AFDB no lo publica en el PDB |
| `length` | entero o vacío | ESMFold: `metrics.tsv`; AFDB: residuos contados en el PDB |
| `path` | ruta POSIX | relativa a la raíz del repo si el archivo está dentro; absoluta si no |
| `sha256` | hex o vacío | SHA-256 del archivo PDB |
| `bytes` | entero | tamaño del PDB |
| `status` | `ok` \| `missing` \| `error` \| `not_found` | salud del archivo/fila; no indica rol de deduplicación |
| `duplicate_of` | accesión o vacío | accesión canónica cuando otra accesión tiene el mismo `sha256` |
| `updated` | ISO-8601 o vacío | `updated` de AFDB o `ingested_at` de ESMFold; nunca la hora de generación |

`duplicate_of` vacío no significa “único”: una fila puede ser la única de su
accesión y aun así compartir hash con otra; en ese caso la canónica tiene
`duplicate_of` vacío y las demás apuntan a ella.

## Reglas de deduplicación

1. **Una fila por accesión.** Si el mismo catálogo de origen repite una
   accesión, se conserva la fila con archivo utilizable, mejor estado y mayor
   pLDDT; las demás se reportan como `duplicate_accession`.
2. **Misma accesión en AFDB y ESMFold.** Se prefiere AFDB si su archivo
   existe y su estado es `ok`. Si AFDB está en error o le falta el PDB, se
   conserva ESMFold. Si ninguna copia es utilizable, se conserva la fila
   AFDB para no perder la procedencia del error.
3. **Mismo `sha256` en accesiones distintas.** Se elige una canónica con la
   clave: archivo utilizable > estado > fuente AFDB > mayor pLDDT > accesión
   lexicográfica menor. Las demás conservan su fila y archivo, con
   `duplicate_of` apuntando a la canónica.
4. **Nunca se borra ni se mueve un archivo.** Los PDB huérfanos y los
   duplicados siguen en disco.
5. **Sin red.** El unificador solo lee catálogos, PDB y tamaño/hash local.

La preferencia AFDB es de procedencia, no una comparación de calidad entre
pLDDT de AFDB y ESMFold: las escalas se normalizan a 0–1, pero provienen de
sistemas distintos.

## Reportes

`reports/structure_catalog_report.tsv` tiene las columnas:

```text
kind  accession  source  path  sha256  bytes  related  detail
```

`kind` puede ser:

| kind | significado |
|---|---|
| `duplicate_accession` | accesión repetida (entre fuentes o dentro de una fuente); `related` indica el ganador |
| `duplicate_hash` | mismo `sha256` en otra accesión; `related` es la accesión canónica |
| `orphan_pdb` | PDB en disco no registrado en el catálogo de esa fuente, o registrado sin ruta |
| `path_mismatch` | PDB en disco cuya accesión está registrada con otra ruta |
| `missing_file` | el catálogo declara una ruta que no existe |

El Markdown resume conteos, muestra hasta `--report-limit` filas por sección
(200 por defecto) y documenta los límites. El TSV no se trunca.

## Uso

Primero, siempre en seco:

```bash
.venv/bin/python scripts/structure/unify_structure_catalog.py --dry-run
```

Salida esperada (formato):

```text
[dry-run] AFDB=... ESMFold=... unificado=... (...) dups_acc=... dups_hash=... huérfanos=... path_mismatch=... faltantes=...
dry-run: no se escribió ningún archivo
```

Para automatización:

```bash
.venv/bin/python scripts/structure/unify_structure_catalog.py --dry-run --json
```

Escritura real (solo cuando el dry-run sea aceptable):

```bash
.venv/bin/python scripts/structure/unify_structure_catalog.py
```

Opciones útiles:

```text
--no-orphans         no escanear directorios de PDB
--no-hash-orphans    no calcular SHA-256 de huérfanos
--report-limit N     filas por sección en el Markdown
--allow-empty        escribir catálogo vacío si no hay entradas
--json               resumen JSON por stdout
```

## Estado de referencia (dry-run)

Corrida del **2026-11-24** sobre los datos locales, sin escribir:

| métrica | valor |
|---|---:|
| filas AFDB de entrada | 620 |
| filas ESMFold de entrada | 1.565 |
| accesiones unificadas | 2.184 |
| fuente AFDB | 620 |
| fuente ESMFold | 1.564 |
| estado `ok` | 2.042 |
| estado `error` (AFDB 404/otros) | 142 |
| `missing` / `not_found` | 0 / 0 |
| duplicados por accesión | 1 |
| duplicados por hash | 0 |
| PDB huérfanos | 1.192 |
| rutas no coincidentes | 0 |
| archivos faltantes | 0 |

La única accesión en ambos catálogos es `A0A089QKZ7`; se conservó AFDB.
Los 1.192 huérfanos son PDB de `afdb25k/` que el catálogo AFDB actual no
registra (incluye 4 accesiones con fila en error y sin ruta). No se
borraron: solo se reportan.

Estos números cambian cuando se actualizan los catálogos; el dry-run es la
fuente de verdad.

## Idempotencia

- filas ordenadas por accesión y reportes por `(kind, accession, source, path)`;
- no se usa la hora actual en el catálogo ni en el Markdown;
- escritura atómica (`*.tmp` + `os.replace`);
- correr dos veces con las mismas entradas produce los mismos bytes;
- `ingest_kaggle_zips.py` conserva el `ingested_at` original y no reescribe
  un PDB ya catalogado con peor pLDDT.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_structure_catalog.py -q
.venv/bin/python -m pytest tests/test_structure_catalog.py tests/test_ingest_kaggle.py -q
```

Las pruebas usan `tmp_path`, catálogos sintéticos y un ZIP sintético. No
tocan `data/metabolism/targets`.

## Límites

- `plddt` AFDB es una media de B-factor de CA del primer modelo; no sustituye
  el pLDDT por residuo ni el valor oficial de AlphaFold DB.
- AFDB no aporta pTM en el PDB: `ptm` queda vacío.
- `length` AFDB se cuenta desde el PDB, no desde la secuencia UniProt; puede
  diferir del largo de la proteína completa si el modelo está truncado.
- El catálogo AFDB puede traer columnas aditivas (`sequence_sha256`,
  `duplicate_of`, `source`, `compressed`); el unificador las ignora y usa
  las columnas clásicas (`pdb`, `sha256`, `pdb_bytes`, `updated`).
- La deduplicación por hash detecta bytes idénticos, no la misma proteína
  predicha con coordenadas distintas ni duplicados por secuencia.
- Solo se inspeccionan `*.pdb` en los directorios configurados; no hay mmCIF,
  ensamblados biológicos, ligandos ni validación geométrica.
- No se valida que el PDB sea biológicamente correcto ni que la accesión
  corresponda a la secuencia modelada.
- El reporte de huérfanos puede ser grande (1.192 filas en la corrida de
  referencia); el Markdown lo trunca, el TSV no.
- El catálogo no reemplaza la evidencia experimental ni justifica por sí solo
  anotaciones funcionales.
- No se publican ni se suben los PDB; la política de datos sigue en
  `docs/CATALOG_PUBLISHING.md` y `DATA_NOTICE.md`.
