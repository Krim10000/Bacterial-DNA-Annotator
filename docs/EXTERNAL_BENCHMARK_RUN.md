# Benchmark externo CPU: Bakta, Prokka y DFAST

**Runner:** `benchmarks/external/external_benchmark.py`  
**Gold local:** `benchmarks/external/make_gold.py`  
**Receta CPU:** `benchmarks/external/Dockerfile` y `benchmarks/external/docker-compose.yml`  
**Alternativa local:** `benchmarks/external/prepare_tools.sh`  
**Pruebas:** `tests/test_external_benchmark_real.py` y `tests/test_external_benchmark.py`

> **Estado honesto en este entorno (2026-11-25):** `bakta`, `prokka` y `dfast`
> **no están instalados** (`shutil.which` no los encuentra) y **no hay Bakta DB**
> (`BAKTA_DB` sin definir; `~/.bakta_db` y `/data/bakta_db` no existen).
> Por eso el informe real generado deja las tres filas en `not_available` y
> **no contiene métricas inventadas**. La imagen Docker **no se construyó ni se ejecutó** y las tres herramientas **no se ejecutaron**. Este documento
> describe la receta reproducible, el gold pequeño y los comandos exactos para
> ejecutar cuando el entorno tenga las dependencias.

## Qué mide y qué no

El runner recibe un FASTA y una referencia local en GFF3 o GenBank/GBFF.
Compara **CDS, tRNA y rRNA** de cada herramienta contra la referencia:

- por `seqid` (con normalización opcional de accesiones NCBI), tipo, hebra y
  coordenadas 1-based inclusivas;
- emparejamiento 1:1 por **solapamiento recíproco** (0,5 por defecto);
- TP, FP, FN, precisión, recall y F1 por categoría y globales;
- coincidencia exacta de coordenadas y diferencias de borde (`Δstart`, `Δend`);
- diferencias concretas: falsos negativos, falsos positivos, pares con hebra
  distinta y pares con bordes distintos;
- atributos `product`, `gene` y `locus_tag` en pares emparejados;
- porcentaje de CDS hipotéticas cuando la anotación trae `product`;
- tiempo de pared y RSS máximo del proceso hijo cuando se ejecuta la
  herramienta. El RSS es `ru_maxrss` de `RUSAGE_CHILDREN` (máximo acumulado de
  hijos, no un perfil instantáneo).

No evalúa función, dominios, rutas metabólicas, expresión ni calidad biológica
global. La referencia es el gold de esta comparación: hereda sus criterios y sus
errores. Tampoco declara que una herramienta sea superior; para eso hacen falta
un gold independiente y ejecuciones pareadas.

## Estados

| estado | significado | ¿hay métricas? |
|---|---|---|
| `ok` | la anotación se cargó y se comparó | sí |
| `not_available` | falta el binario o falta una Bakta DB válida; no se ejecutó nada | **no** |
| `not_run` | el binario existe, pero no se pasó `--run` | **no** |
| `failed` | la herramienta falló, no dejó salida legible o no se pudo parsear | **no** |
| `timeout` | se agotó `--timeout` | **no** |
| `not_comparable` | no hay contigs en común o la salida no tiene features | **no** |
| `invalid_input` | falta un archivo o no se puede parsear | **no** |

Una herramienta ausente **no es un cero ni un fracaso del método**. Su fila
queda como `not_available` con la razón exacta (por ejemplo,
`binario 'bakta' no encontrado en PATH (shutil.which)`) y sin categorías,
agregados ni F1 inventados.

Bakta tiene una regla adicional: si el binario existe pero la **Bakta DB** no
está o no parece válida, `--run` deja Bakta en `not_available` y **no ejecuta**
el binario. No se descarga la base automáticamente.

## Entradas

- **FASTA** del genoma, plano o `.gz`. Se usa para los contigs y para normalizar
  seqids.
- **Referencia** en GFF3 (`.gff`, `.gff3`, `.gtf`, con o sin `.gz`) o
  GenBank/GBFF (`.gbff`, `.gbk`, `.gb`, `.genbank`, con o sin `.gz`).
  El formato se autodetecta; se puede forzar con `--reference-format`.
- Opcionalmente, **GFF/GBFF precalculados** por herramienta con
  `--query-gff TOOL=RUTA`. Esto permite comparar salidas producidas en otro
  host o contenedor sin tener el binario instalado; la fuente queda como
  `precomputed`.

Los CDS marcados como `pseudo`/`pseudogene` se excluyen por defecto
(`--keep-pseudo` para conservarlos). En GFF3, varias líneas con el mismo `ID`
(por ejemplo un CDS con frameshift programado) se fusionan en un feature con el
intervalo envolvente, igual que GenBank. En GenBank, las ubicaciones compuestas
también se representan por su intervalo envolvente; se registra esta limitación.

### Normalización de seqid

Una referencia NCBI puede usar `NC_000913.3` y el FASTA una cabecera
`>NC_000913.3 Escherichia coli...`. El runner construye alias a partir de:

- la cabecera completa y su primer token;
- formatos tipo `gi|...|ref|NC_000913.3|`;
- accesiones con y sin versión (`NC_000913.3` y `NC_000913`).

Se puede desactivar con `--no-seqid-aliasing`. Los seqids no encontrados y los
alias ambiguos se reportan en `sequence_resolution` y en `warnings`.

## Gold pequeño e independiente desde un GFF3/GBFF local

`make_gold.py` extrae una ventana de un ensamblado y su anotación local (por
ejemplo NCBI RefSeq) y la reindexa a 1. El gold **no se deriva de Bakta,
Prokka ni DFAST** y no descarga nada. Es pequeño para poder ejecutar el
benchmark en CPU y comparar exactamente las mismas coordenadas.

Ejemplo real ya generado en este repositorio
(`Chlorobaculum tepidum TLS`, ventana de 60 000 bp, 55 CDS + 1 tRNA):

```bash
.venv/bin/python benchmarks/external/make_gold.py \
  --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
  --annotation data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
  --out-dir benchmarks/external/gold/chlorobaculum_tepidum_60kb \
  --name chlorobaculum_tepidum_60kb \
  --max-bp 60000 --flank 250 --min-features 20 --force
```

Salidas del gold:

| archivo | contenido |
|---|---|
| `genome.fna` | ventana FASTA sin comprimir; coordenadas relativas a la ventana |
| `reference.gff3` | CDS/tRNA/rRNA completos de la ventana, 1-based inclusivo |
| `gold_manifest.json` | procedencia, SHA256, selección, conteos y comando `--run` |
| `run_benchmark.sh` | ejecuta el benchmark con `--run` para ese genoma |

El manifiesto incluye:

- SHA256 y tamaño del FASTA y de la anotación de origen;
- SHA256 y tamaño de `genome.fna` y `reference.gff3`;
- `seqid`, ventana (`window_start`/`window_end`), bp extraídos y features
  omitidos por coordenadas o seqid;
- conteos CDS/tRNA/rRNA y limitaciones;
- el comando exacto con `--run`, el script y las rutas de JSON/MD/evidencia.

Se puede usar un GBFF local **sin FASTA externo**: en ese caso las secuencias se
extraen del propio GBFF con Biopython. Para un GFF3 hace falta `--genome`.

> **Límite:** el gold hereda los criterios y errores de la anotación de origen;
> no es un gold curado a mano. Las coordenadas del gold están reindexadas y no
> coinciden con el ensamblado completo.

## Bakta DB: detección y validación local

El runner **no descarga** la base. La busca en este orden:

1. `--bakta-db RUTA`;
2. variable de entorno `BAKTA_DB`;
3. `~/.bakta_db`;
4. `/data/bakta_db`.

En cada ruta comprueba que sea un directorio y busca marcadores conocidos:
`version.json`, `db_info.json`, `bakta.db`, `amrfinderplus`, `uniprot`
(también acepta una base anidada en `<ruta>/db`). Si existe `version.json` o
`db_info.json`, registra la versión y el SHA256 de ese JSON pequeño. **No
recorre ni hashea la base completa**, que puede ser enorme.

Estados de la base en el JSON (`environment.bakta_db` y
`validation.bakta_db`):

| estado | significado | ¿se usa con `--run`? |
|---|---|---|
| `ok` | directorio con marcadores y versión legible | sí |
| `incomplete` | directorio sin marcadores conocidos | **no** |
| `invalid` | la ruta existe pero no es un directorio | **no** |
| `not_found` | no se encontró ninguna ruta candidata | **no** |

Comprobación local sin ejecutar Bakta:

```bash
export BAKTA_DB=/ruta/a/bakta_db
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools bakta --run --threads 1
```

Si la base no es válida, Bakta queda `not_available` con la razón y no se
ejecuta. La regla se aplica al ejecutable canónico `bakta`; un envoltorio con
otro nombre puede gestionar su propia base y el runner deja que la ejecución
reporte el fallo real. Prokka y DFAST no requieren una base obligatoria en el
runner, aunque sus dependencias pueden necesitar configuración adicional según
la versión.

## Uso local (sin Docker)

Sin `--run` el runner solo detecta binarios y valida la base. Si no están
instalados, las tres filas quedan `not_available`:

```bash
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools all \
  --json reports/external_benchmark/benchmark.json \
  --md reports/external_benchmark/benchmark.md
```

Con `--run`, la base Bakta se toma de `BAKTA_DB` o de `--bakta-db`:

```bash
export BAKTA_DB=/ruta/a/bakta_db   # opcional si no se usa Bakta
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools all \
  --run \
  --bakta-db /ruta/a/bakta_db \
  --threads 1 \
  --timeout 3600 \
  --save-gffs reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence \
  --json reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.json \
  --md reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.md
```

El gold ya trae el comando exacto y un script equivalente:

```bash
cat benchmarks/external/gold/chlorobaculum_tepidum_60kb/gold_manifest.json
bash benchmarks/external/gold/chlorobaculum_tepidum_60kb/run_benchmark.sh
```

`run_benchmark.sh` usa `PROJECT_ROOT` o el `git` del repositorio para localizar
la raíz, `.venv/bin/python` si existe, `--threads 1`, `--run`, evidencia y
JSON/MD bajo `reports/external_benchmark/<gold>/`. Se le pueden añadir
argumentos extra al final; por ejemplo `--timeout 7200`.

Para comparar una salida ya producida (sin binario ni `--run`):

```bash
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome genome.fna.gz \
  --reference reference.gbff.gz \
  --tools prokka \
  --query-gff prokka=/ruta/a/prokka.gff \
  --json reports/external_benchmark/prokka_compare.json \
  --md reports/external_benchmark/prokka_compare.md
```

El runner acepta `--extra-arg TOOL=ARG` repetible para opciones propias de una
versión concreta de la herramienta. Por ejemplo:
`--extra-arg dfast=--config --extra-arg dfast=/ruta/config.py`.

## Receta Docker CPU

### Construcción

Desde la raíz del repositorio:

```bash
docker build \
  -f benchmarks/external/Dockerfile \
  -t dna-token-bac/benchmark-external:cpu \
  .
```

La imagen parte de `mambaorg/micromamba:1.5.10`, instala Python 3.11 +
Biopython para el runner y crea **un entorno separado por herramienta**:

| herramienta | versión fijada por defecto | entorno |
|---|---:|---|
| Bakta | 1.12.2 | `bakta` |
| Prokka | 1.15.6 | `prokka` |
| DFAST | 1.6.2 | `dfast` |

Son versiones de Bioconda consultadas el 2026-11-24. Se pueden sobreescribir:

```bash
docker build \
  -f benchmarks/external/Dockerfile \
  -t dna-token-bac/benchmark-external:cpu \
  --build-arg BAKTA_VERSION=<version> \
  --build-arg PROKKA_VERSION=<version> \
  --build-arg DFAST_VERSION=<version> \
  .
```

No se garantiza que una versión concreta siga disponible en el canal; hay que
verificarla antes de fijarla. **El build requiere red** para descargar la imagen
base y los paquetes de conda-forge/bioconda. La imagen no incluye bases de datos
ni las descarga durante el build.

### Ejecución

```bash
docker run --rm \
  -v "$PWD:/work" \
  -v "/ruta/host/bakta_db:/data/bakta_db:ro" \
  dna-token-bac/benchmark-external:cpu \
  python /opt/benchmark/external/external_benchmark.py \
    --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
    --reference /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
    --tools all --run --bakta-db /data/bakta_db \
    --json /work/reports/external_benchmark/benchmark.json \
    --md /work/reports/external_benchmark/benchmark.md
```

Con Docker Compose:

```bash
BAKTA_DB_DIR=/ruta/host/bakta_db \
docker compose -f benchmarks/external/docker-compose.yml run --rm benchmark \
  --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools prokka,dfast --run \
  --json /work/reports/external_benchmark/benchmark.json \
  --md /work/reports/external_benchmark/benchmark.md
```

El servicio monta el repositorio en `/work` y `BAKTA_DB_DIR` en
`/data/bakta_db:ro`. Si no se ejecuta Bakta, se puede omitir el montaje y la
variable. Para comprobar la imagen sin ejecutar herramientas:

```bash
docker run --rm dna-token-bac/benchmark-external:cpu \
  python /opt/benchmark/external/external_benchmark.py --help
docker run --rm dna-token-bac/benchmark-external:cpu \
  sh -lc 'which bakta; which prokka; which dfast'
```

La imagen es **CPU-only**: no instala drivers propietarios ni dependencias de
aceleración. El runner fuerza 1 hilo/CPU por herramienta para que los tiempos
sean comparables.

### Alternativa sin Docker

`prepare_tools.sh` crea tres entornos micromamba locales con las mismas
versiones y envoltorios en `benchmarks/external/.tools/bin`:

```bash
bash benchmarks/external/prepare_tools.sh
export PATH="$PWD/benchmarks/external/.tools/bin:$PATH"
python benchmarks/external/external_benchmark.py --help
```

Requiere `micromamba` y red para resolver paquetes. No descarga bases de datos.
El directorio `.tools/` es de ejecución local y no debe versionarse.

## Validación de rutas, versiones, SHA256 y tiempos

Antes de ejecutar nada, el runner valida:

- que el FASTA y la referencia existan, sean archivos y se puedan parsear;
- el formato real de la referencia (`gff3`/`genbank`) y sus conteos;
- que cada binario detectado sea un archivo ejecutable;
- la Bakta DB (estado, marcadores y versión si está disponible);
- los parámetros (`--threads`, `--timeout`, umbral, máximo de diferencias).

Durante/después de una ejecución registra:

- comando real, `returncode`, `stdout`/`stderr` (cola), timeout;
- `started_at`/`finished_at` y `runtime_seconds` por herramienta;
- `max_rss_kb`/`max_rss_mib` de `RUSAGE_CHILDREN`;
- tamaño, mtime y SHA256 de la anotación de salida;
- versión de la herramienta (`<tool> --version`) y versión parseada;
- `validation.checks`, `validation.inputs`, `validation.versions` y
  `validation.timing` en el JSON.

El JSON y el Markdown se escriben de forma atómica en:

```text
reports/external_benchmark/benchmark.json
reports/external_benchmark/benchmark.md
```

o bajo `reports/external_benchmark/<gold>/` cuando se usa el script del gold.

## Métricas

Para cada categoría:

- `reference_n`, `query_n`;
- `TP`, `FP`, `FN` según emparejamiento 1:1;
- `precision = TP / (TP + FP)` (si el denominador es 0, `null`);
- `recall = TP / (TP + FN)` (si el denominador es 0, `null`);
- `F1 = 2·TP / (2·TP + FP + FN)` (si el denominador es 0, `null`);
- `exact_matches`: misma `seqid`, tipo, hebra, `start` y `end`;
- `boundary_shifted`: pares emparejados con algún borde distinto;
- `mean_reciprocal_overlap`;
- `mean_abs_start_delta`, `mean_abs_end_delta`;
- atributos `product`, `gene`, `locus_tag` de los pares;
- `hypothetical` de referencia y consulta.

El agregado global usa TP/FP/FN sumados (micro) y también reporta `macro_f1`
(promedio de los F1 por categoría). `F1 = 0` cuando hay FN pero no TP; la
precisión queda `null` si no hay predicciones.

`product` se compara de dos formas: literal normalizada y Jaccard de tokens.
Ninguna de las dos es accuracy funcional: dos herramientas pueden describir la
misma proteína con nombres distintos.

`% hypothetical` se calcula solo sobre CDS con `product` no vacío. Se consideran
hipotéticas las cadenas que contienen `hypothetical`, `uncharacterized`,
`unknown protein`, `protein of unknown function`, `DUF...` o `UPF...`.
`putative` a secas **no** cuenta como hipotética.

## Salidas

`--json` produce un documento con:

- `schema`, `generated_at`, `command`, `scope`;
- `parameters`: umbral, `run`, herramientas, timeout, hilos, pseudogenes,
  normalización de seqid y base Bakta;
- `environment`: Python, plataforma, CPUs, detección/versión de cada
  herramienta y `bakta_db`;
- `inputs`: FASTA y referencia con **SHA256**, tamaño, mtime, conteos y contigs;
- `validation`: checks, huellas de entrada, Bakta DB, versiones y tiempos;
- `results` por herramienta: estado, razón, fuente (`executed`/`precomputed`),
  versión, ejecución, anotación comparada y comparación completa;
- `not_available`: lista de herramientas seleccionadas sin binario/base;
- `warnings`, `limitations` y `reproduction` (incluido `run_command`).

Cada comparación incluye `categories` (CDS/tRNA/rRNA), `aggregate`,
`attributes`, `hypothetical`, `differences` y `sequence_resolution`. Las listas
de diferencias se truncan a `--max-differences` elementos por lista; los conteos
son completos y `differences.truncated` lo indica.

`--md` produce un informe legible con entorno, Bakta DB, entradas, checks de
validación, tabla resumen, tablas por herramienta, diferencias y limitaciones.

`--save-gffs DIR` copia los GFF/GBFF de salida y `command.txt`, `stdout.txt`,
`stderr.txt` por herramienta para auditar la ejecución.

## Reproducibilidad

1. Fijar las versiones de las tres herramientas (build args o
   `prepare_tools.sh`).
2. Registrar la digestión/ID de la imagen y la versión de micromamba.
3. Conservar el `--json`: incluye SHA256 del FASTA, la referencia, la salida de
   cada herramienta y la Bakta DB (`version.json`/`db_info.json`).
4. Ejecutar con `--threads 1` y un `--timeout` explícito.
5. Usar el mismo FASTA, la misma referencia y el mismo umbral para todas las
   herramientas.
6. Guardar los GFF/GBFF con `--save-gffs`.
7. Repetir por genoma con el `run_benchmark.sh` del gold correspondiente.

El runner intenta registrar `<tool> --version` al usar `--run` (se puede
desactivar con `--no-probe-versions`). Un fallo al obtener la versión no invalida
la comparación, pero queda visible.

## Limitaciones

1. **No hay resultados reales de herramientas en este repositorio**: no se
   construyó la imagen ni se ejecutaron Bakta/Prokka/DFAST. Cualquier informe
   generado aquí debe mostrar `not_available`, no métricas inventadas.
2. **La referencia es el gold**: si es NCBI, hereda sus criterios y errores.
3. **Emparejamiento voraz 1:1**: no es un alineamiento global ni un óptimo
   bipartito. Cambiar el umbral puede cambiar TP/FP/FN.
4. **Coordenadas**: se comparan `start`/`end` 1-based inclusivos. Genes con
   intrones o ubicaciones compuestas de GenBank usan el intervalo envolvente.
   En el gold, las coordenadas están reindexadas a la ventana.
5. **Nombres**: `product`, `gene` y `locus_tag` dependen de la base y la
   convención de cada herramienta; la coincidencia literal no mide función.
6. **Categorías**: no se comparan ncRNA, tmRNA, CRISPR, pseudogenes (si se
   excluyen), función, dominios ni rutas.
7. **Bases de datos**: Bakta necesita una Bakta DB local; Prokka/DFAST pueden
   necesitar bases o configuración adicional. La imagen no las incluye ni las
   descarga. Revisar sus licencias y términos.
8. **Red en el build**: construir la imagen requiere red. El benchmark en sí no
   descarga datos, pero una herramienta de terceros podría intentar acceder a
   servicios externos; `--run` es una decisión explícita.
9. **Sin conclusión de superioridad**: faltan un gold curado a mano,
   ejecuciones pareadas y análisis estadístico.
10. **Versiones**: los pines por defecto son una foto de Bioconda; si el canal
    deja de servirlos, el build falla y hay que actualizarlos explícitamente.

## Comando exacto y bloqueo actual

Comando que ejecutaría el benchmark completo con el gold real:

```bash
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools all \
  --run \
  --threads 1 \
  --timeout 3600 \
  --json reports/external_benchmark/benchmark.json \
  --md reports/external_benchmark/benchmark.md
```

**Bloqueo:** faltan los binarios `bakta`, `prokka` y `dfast` y falta una Bakta
DB local. Sin ellos no hay ejecución ni métricas; el comando anterior deja las
tres filas en `not_available` y escribe el informe en
`reports/external_benchmark/`.

Para desbloquear, elegir una vía:

```bash
# Opción A: imagen CPU (requiere Docker + red para el build)
docker build -f benchmarks/external/Dockerfile \
  -t dna-token-bac/benchmark-external:cpu .

# Opción B: entornos micromamba locales (requiere micromamba + red)
bash benchmarks/external/prepare_tools.sh
export PATH="$PWD/benchmarks/external/.tools/bin:$PATH"

# En ambas: conseguir una Bakta DB local y exportarla (no se descarga aquí)
export BAKTA_DB=/ruta/a/bakta_db
```

Después, comprobar `which bakta prokka dfast` y repetir el comando exacto.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_external_benchmark_real.py -q
.venv/bin/python -m pytest tests/test_external_benchmark.py -q
```

Las pruebas usan GFF3, GenBank y FASTA sintéticos, binarios falsos de Python,
directorios de Bakta DB sintéticos y binarios ausentes. **No usan red, no
descargan imágenes/bases y no ejecutan Bakta/Prokka/DFAST reales.** También
validan estáticamente la receta Docker, el compose, el script de preparación,
el gold y este documento.
