# Caché AFDB/PAE

## Qué hay

- `data/metabolism/targets/afdb25k/*.pdb` o `*.pdb.gz`: modelos AFDB locales.
- `data/metabolism/targets/afdb25k/*.pae.json` o `*.pae.json.gz`: PAE.
- `data/metabolism/targets/afdb_catalog.tsv`: catálogo por accesión.

El catálogo anterior sigue siendo legible: las columnas nuevas se añaden al
final y los lectores viejos pueden ignorarlas.

## Descargador

```bash
# Piloto conservador: como máximo 10 accesiones nuevas por lote
python scripts/structure/fetch_afdb_cache.py --limit 10 --pae \
    --workers 1 --delay 1.5

# Ver el plan sin red y sin escribir nada
python scripts/structure/fetch_afdb_cache.py --limit 10 --dry-run

# Archivos nuevos comprimidos (no toca los ya existentes)
python scripts/structure/fetch_afdb_cache.py --limit 10 --pae --compress
```

Opciones principales:

| Opción | Efecto |
|---|---|
| `--limit N` | máximo de accesiones **pendientes** en el lote; `0` = todas |
| `--pae` | descarga también el PAE JSON |
| `--compress` | guarda descargas nuevas como `.pdb.gz` / `.pae.json.gz` |
| `--workers N` | hilos; para 429 bajar a `1` |
| `--delay S` | espera antes de cada solicitud de red por worker (def. `0.5`) |
| `--retries N` | intentos totales por accesión (def. `3`) |
| `--timeout S` | timeout de metadatos; los archivos usan `3 * timeout` |
| `--backoff-base S` | base del backoff exponencial (def. `2`) |
| `--backoff-max S` | tope del backoff (def. `120`) |
| `--retry-after-max S` | tope para `Retry-After` (def. `120`) |
| `--flush-every N` | guarda el catálogo cada N accesiones (def. `500`; el anterior era `5`) |
| `--dry-run` | imprime pendientes; no descarga ni escribe |

## Catálogo

| Columna | Significado |
|---|---|
| `accession` | accesión UniProt única (clave del catálogo) |
| `status` | `ok`, `already_present`, `not_found`, `error`, `rate_limited` |
| `pdb`, `pae` | ruta local (`.pdb`/`.pdb.gz`, `.pae.json`/`.pae.json.gz`) |
| `sha256` | SHA-256 del **contenido descomprimido** del PDB |
| `pdb_bytes`, `pae_bytes` | bytes **en disco** (comprimidos si aplica) |
| `pdb_url`, `pae_url` | URLs usadas, para reanudar sin pedir metadatos |
| `error` | último error legible |
| `updated` | marca de tiempo local |
| `sequence_sha256` | SHA-256 de la secuencia normalizada (sin espacios, A-Z) |
| `duplicate_of` | accesión canónica si otra accesión tiene la misma secuencia |
| `source` | `alphafold_db` si vino de la API; `local_cache` si ya estaba |
| `compressed` | `true` si el PDB o el PAE local están en gzip |

La deduplicación por accesión la garantiza la clave del catálogo. La
deduplicación por secuencia agrupa filas con el mismo `sequence_sha256`:

- la secuencia se toma de `uniprotSequence` de la API;
- para archivos locales sin metadatos se reconstruye desde `SEQRES`;
- el canónico es la accesión con modelo local (`ok`/`already_present`) y, en
  empate, la menor lexicográficamente;
- las demás filas quedan con `duplicate_of=<canónica>`. No se borra ningún
  archivo: el marcador es informativo.

## Compresión y lectura transparente

- `--compress` solo afecta a descargas nuevas. Un `.pdb`/`.pae.json` existente
  se reutiliza tal cual y **nunca se reescribe ni se borra**.
- El módulo expone lectura transparente para `.gz` y archivos planos:

```python
from scripts.structure.fetch_afdb_cache import (
    cached_path, read_cached_bytes, read_cached_text,
)

pdb = cached_path(targets, "P12345", "pdb")   # .pdb o .pdb.gz
texto = read_cached_text(pdb)
pae = cached_path(targets, "P12345", "pae")   # .pae.json o .pae.json.gz
```

- `compressed=false` puede convivir con `compressed=true` en el mismo lote:
  cada fila describe sus propios archivos.
- **Límite actual de consumidores**: `webapp/gene_card.py`,
  `webapp/structural_b.py` y `webapp/pae_utils.py` buscan nombres sin
  comprimir. Mientras no se les añada un resolvedor, dejar `--compress`
  apagado para la caché que consume la GUI, o usar los helpers de arriba en
  código nuevo.

## Reanudación

- Antes de descargar se revisan los archivos locales (planos o `.gz`).
- Si el PDB existe y no se pide PAE, no hay ninguna petición de red.
- Si solo falta el PAE y el catálogo ya tiene `pae_url`, se descarga solo el
  PAE (sin volver a pedir metadatos ni el PDB).
- El lote se ordena por accesión; `--limit` corta solo pendientes, por lo que
  repetir el mismo comando avanza de forma determinista.
- Al arrancar, una pasada local sin red completa `sha256`, `sequence_sha256`,
  rutas y `compressed` de los PDB ya presentes, aunque el catálogo estuviera
  incompleto.

## Auditoría de estado (solo lectura)

Antes de cada lote conviene medir la deuda sin descargar:

```bash
# Ver métricas, plan y siguiente lote; no escribe nada
python scripts/structure/afdb_status.py --dry-run

# Generar reports/afdb_status.tsv y reports/afdb_status.md
python scripts/structure/afdb_status.py
```

El script **no hace red** y **no modifica** el catálogo ni la caché: solo lee
el catálogo, los archivos locales y `data/regulation/*/genes.tsv`. Clasifica
cada fila del catálogo en un estado efectivo:

| estado | significado | ¿se reintenta? |
|---|---|---|
| `ok` | `ok`/`already_present` con PDB local | no hace falta |
| `rate_limited` | 429 o `Retry-After` agotado | **sí** (temporal) |
| `not_found` | 404 en la API o en un archivo | **no** (permanente) |
| `error` | 5xx, timeout, red u otro 4xx | sí |
| `local_only` | PDB en disco sin fila de catálogo | no; cuenta como cobertura |

El catálogo puede conservar `status=error` genérico. La auditoría mira el
código HTTP del campo `error` y separa 429 de 404 **sin reescribir el
catálogo**; esa clasificación es la que usa el plan.

Cobertura: el universo son las accesiones únicas de
`data/regulation/*/genes.tsv` (111.569). El informe distingue cobertura de
catálogo (`ok` con archivo), cobertura en disco (cualquier PDB local) y
pendientes = universo − PDB local − 404 permanentes. También cuenta PAE
presentes, PAE que faltan para los PDB locales y PDB huérfanos.

Plan de lotes: si hay 429, empieza con `--limit 10 --workers 1 --delay 2.0` y
sube solo si el lote queda limpio (≥ 9/10 `ok`). El `.md` incluye
workers/delay/limit por paso, ETA aproximada y las accesiones del siguiente
lote en el orden determinista del descargador. Durante la recuperación de PDB
no usar `--pae`: duplica las solicitudes.

Si el catálogo tiene filas con `"HTTP Error 404"`, el informe avisa de una
discrepancia: `fetch_afdb_cache.is_permanent()` solo reconoce `"HTTP 404"`, así
que el descargador las sigue contando como pendientes. La auditoría las marca
`not_found` (permanentes) y recomienda un `--limit` algo mayor para cerrarlas en
la misma pasada. La corrección de fondo es ampliar `is_permanent()` a
`"HTTP Error 404"`; la auditoría no modifica el descargador.

## 404, 429 y errores

- **404** en la API o en un archivo: `status=not_found`, `error` contiene
  `HTTP 404` o `HTTP Error 404` (el descargador solo reconoce la primera
  forma; la auditoría reconoce ambas). Es permanente: no se reintenta ni se
  vuelve a pedir en ejecuciones posteriores (salvo que se borre/edite esa fila
  del catálogo).
- **429**: se respeta `Retry-After` (segundos o fecha HTTP) con tope
  `--retry-after-max`; sin cabecera se usa backoff exponencial con jitter.
  Si se agotan los intentos queda `status=rate_limited` (temporal) y el
  siguiente lote lo vuelve a intentar. Un 429 **no** es un 404.
- **5xx / red / timeout**: backoff exponencial; si se agotan, `status=error`
  (temporal) y se reintenta en otro lote.
- **Otros 4xx** (400/401/403): `status=error` sin reintento inmediato.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_afdb_cache.py -q
.venv/bin/python -m pytest tests/test_afdb_status.py -q
```

Las pruebas son sintéticas y no usan red (`urlopen` se sustituye por un doble).
Las de la caché cubren compresión/descompresión, catálogo y deduplicación, 404
permanente, 429 con `Retry-After` y agotamiento, reanudación (PDB y PAE),
límite por lote y `--dry-run`. Las de la auditoría cubren un catálogo
sintético con `ok`/429/404/error, cobertura sobre el universo, pendientes,
siguiente lote, plan/ETA, archivos `.gz`, TSV/MD y la garantía de que ni el
catálogo ni la caché cambian.

## Estrategia de espacio

- una accesión UniProt = una entrada;
- deduplicar por hash de secuencia cuando dos accesiones tengan la misma
  proteína (marcador `duplicate_of`, sin borrar archivos);
- comprimir PDB con gzip (`--compress`): PDB y PAE ocupan menos;
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

- AFDB no cubre proteínas sin UniProt.
- un modelo AFDB monomérico no es estructura cuaternaria.
- PAE no es evidencia experimental.
- No ejecutar descargas masivas durante pruebas: usar `--limit 10`,
  `--workers 1` y `--delay 1.5`; `--dry-run` para inspeccionar.
