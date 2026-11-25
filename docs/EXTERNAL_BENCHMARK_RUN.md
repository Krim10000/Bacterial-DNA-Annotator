# Benchmark externo CPU: Bakta, Prokka y DFAST

**Runner:** `benchmarks/external/external_benchmark.py`  
**Receta CPU:** `benchmarks/external/Dockerfile` y `benchmarks/external/docker-compose.yml`  
**Alternativa local:** `benchmarks/external/prepare_tools.sh`  
**Pruebas:** `tests/test_external_benchmark.py`

> **Estado honesto:** en este repositorio **no se construyó ni se ejecutó** la
> imagen y **no se ejecutaron las tres herramientas**. El entorno de desarrollo
> no tiene esos binarios. Este documento describe una receta reproducible y un
> runner; no reporta métricas de herramientas que no se corrieron.

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
| `not_available` | el binario no está en `PATH`; no se ejecutó nada | **no** |
| `not_run` | el binario existe, pero no se pasó `--run` | **no** |
| `failed` | la herramienta falló, no dejó salida legible o no se pudo parsear | **no** |
| `timeout` | se agotó `--timeout` | **no** |
| `not_comparable` | no hay contigs en común o la salida no tiene features | **no** |
| `invalid_input` | falta un archivo o no se puede parsear | **no** |

Una herramienta ausente **no es un cero ni un fracaso del método**. Su fila
queda como `not_available` con la razón exacta (por ejemplo,
`binario 'bakta' no encontrado en PATH (shutil.which)`) y sin categorías,
agregados ni F1 inventados.

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

## Uso local (sin Docker)

```bash
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
  --reference data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
  --tools all \
  --json reports/external_benchmark_run/benchmark.json \
  --md reports/external_benchmark_run/benchmark.md
```

Sin `--run` el runner solo detecta binarios. Si no están instalados, las tres
filas quedan `not_available`.

Para ejecutar las herramientas instaladas:

```bash
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome genome.fna.gz \
  --reference reference.gff3 \
  --tools bakta,prokka,dfast \
  --run \
  --bakta-db /ruta/a/bakta_db \
  --threads 1 \
  --timeout 3600 \
  --save-gffs reports/external_benchmark_run/evidence \
  --json reports/external_benchmark_run/benchmark.json \
  --md reports/external_benchmark_run/benchmark.md
```

Bakta necesita una base local; el runner no la descarga. Prokka y DFAST pueden
necesitar bases o dependencias adicionales según la versión.

Para comparar una salida ya producida (sin binario ni `--run`):

```bash
.venv/bin/python benchmarks/external/external_benchmark.py \
  --genome genome.fna.gz \
  --reference reference.gbff.gz \
  --tools prokka \
  --query-gff prokka=/ruta/a/prokka.gff \
  --json reports/prokka_compare.json \
  --md reports/prokka_compare.md
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
    --genome /work/genome.fna.gz \
    --reference /work/reference.gff.gz \
    --tools all --run --bakta-db /data/bakta_db \
    --json /work/reports/external_benchmark_run/benchmark.json \
    --md /work/reports/external_benchmark_run/benchmark.md
```

Con Docker Compose:

```bash
BAKTA_DB_DIR=/ruta/host/bakta_db \
docker compose -f benchmarks/external/docker-compose.yml run --rm benchmark \
  --genome /work/genome.fna.gz \
  --reference /work/reference.gff.gz \
  --tools prokka,dfast --run \
  --json /work/reports/external_benchmark_run/benchmark.json \
  --md /work/reports/external_benchmark_run/benchmark.md
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
  normalización de seqid, base Bakta;
- `environment`: Python, plataforma, CPUs, detección y versión de cada
  herramienta;
- `inputs`: FASTA y referencia con **SHA256**, conteos y contigs;
- `results` por herramienta: estado, razón, fuente (`executed`/`precomputed`),
  versión, ejecución, anotación comparada y comparación completa;
- `not_available`: lista de herramientas seleccionadas sin binario;
- `warnings`, `limitations` y `reproduction`.

Cada comparación incluye `categories` (CDS/tRNA/rRNA), `aggregate`,
`attributes`, `hypothetical`, `differences` y `sequence_resolution`. Las listas
de diferencias se truncan a `--max-differences` elementos por lista; los conteos
son completos y `differences.truncated` lo indica.

`--md` produce un informe legible con entorno, entradas, tabla resumen, tablas
por herramienta, diferencias y limitaciones.

`--save-gffs DIR` copia los GFF/GBFF de salida y `command.txt`, `stdout.txt`,
`stderr.txt` por herramienta para auditar la ejecución.

## Reproducibilidad

1. Fijar las versiones de las tres herramientas (build args o
   `prepare_tools.sh`).
2. Registrar la digestión/ID de la imagen y la versión de micromamba.
3. Conservar el `--json`: incluye SHA256 del FASTA, la referencia y cada
   anotación comparada.
4. Ejecutar con `--threads 1` y un `--timeout` explícito.
5. Usar el mismo FASTA, la misma referencia y el mismo umbral para todas las
   herramientas.
6. Guardar los GFF/GBFF con `--save-gffs`.

El runner intenta registrar `<tool> --version` al usar `--run` (se puede
desactivar con `--no-probe-versions`). Un fallo al obtener la versión no invalida
la comparación, pero queda visible.

## Limitaciones

1. **No hay resultados reales en este repositorio**: no se construyó la imagen
   ni se ejecutaron las herramientas. Cualquier informe generado aquí debe
   mostrar `not_available`, no métricas inventadas.
2. **La referencia es el gold**: si es NCBI, hereda sus criterios y errores.
3. **Emparejamiento voraz 1:1**: no es un alineamiento global ni un óptimo
   bipartito. Cambiar el umbral puede cambiar TP/FP/FN.
4. **Coordenadas**: se comparan `start`/`end` 1-based inclusivos. Genes con
   intrones o ubicaciones compuestas de GenBank usan el intervalo envolvente.
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
9. **Sin conclusión de superioridad**: faltan un gold independiente,
   ejecuciones pareadas y análisis estadístico.
10. **Versiones**: los pines por defecto son una foto de Bioconda; si el canal
    deja de servirlos, el build falla y hay que actualizarlos explícitamente.

## Qué falta para ejecutarlo de verdad

- Construir la imagen (requiere Docker y red) o preparar los entornos locales
  con `prepare_tools.sh`.
- Obtener una **Bakta DB** local y montarla en `/data/bakta_db`; verificar sus
  términos de uso.
- Verificar dependencias/bases de Prokka y DFAST para la versión fijada.
- Elegir un genoma y una referencia gold adecuados (GFF3 o GBFF), con el mismo
  ensamblado y coordenadas.
- Ejecutar con `--run`, `--threads 1`, timeout explícito y `--save-gffs`.
- Revisar el JSON: si una herramienta falla o falta, su estado debe ser
  `failed`/`not_available`, nunca un cero silencioso.
- Para afirmar superioridad, repetir en varios genomas y con un gold
  independiente; este runner no lo sustituye.

## Pruebas

```bash
.venv/bin/python -m pytest tests/test_external_benchmark.py -q
```

Las pruebas usan GFF3 y GenBank sintéticos, FASTA sintético, un binario falso de
Python y binarios ausentes. **No usan red, no descargan imágenes y no ejecutan
Bakta/Prokka/DFAST.** También validan estáticamente la receta Docker, el compose,
el script de preparación y este documento.
