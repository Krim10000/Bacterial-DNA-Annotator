# Entorno Docker para el benchmark externo (CPU)

**Archivo de imagen:** `docker/benchmark.Dockerfile`  
**Script objetivo:** `tools/benchmark_external.py`  
**Herramientas:** Bakta, Prokka y DFAST en CPU.

> **Estado honesto:** esta imagen **no se construyó ni se ejecutó** en este
> repositorio. El Dockerfile es una receta reproducible, no una afirmación de
> que el build o las herramientas funcionen en todas las plataformas. Las
> versiones quedan sin pinear por defecto.

## Qué incluye

- Base `mambaorg/micromamba:latest` (Linux, CPU).
- Python 3.11 y `pytest` en el entorno base.
- Un entorno micromamba separado por herramienta:
  - `bakta` desde `conda-forge`/`bioconda`;
  - `prokka` desde `conda-forge`/`bioconda`;
  - `dfast` desde `conda-forge`/`bioconda`.
- Envoltorios en `/usr/local/bin/{bakta,prokka,dfast}` que ejecutan el binario
  real dentro de su entorno. Así `shutil.which` de `benchmark_external.py`
  encuentra las tres herramientas.
- `tools/benchmark_external.py` copiado en
  `/opt/benchmark/tools/benchmark_external.py`.
- Sin bases de datos incluidas y **sin descargas de bases durante el build**.
- Sin GPU, CUDA ni drivers propietarios.

## Requisitos

- Docker con BuildKit y capacidad de descargar la imagen base y paquetes de
  `conda-forge`/`bioconda`. **El build requiere red**; el script de benchmark no
  descarga datos.
- Espacio en disco suficiente para la imagen y los entornos conda. No se fija
  una cifra porque depende de las versiones que resuelva el solver.
- CPU. No se asume GPU.
- Un GFF de referencia local (por ejemplo el GFF NCBI de `data/annotations/`) y
  el FASTA del genoma. No se incluyen en la imagen.
- Para Bakta, una base local ``Bakta DB`` obtenida por el usuario. La imagen no
  la descarga ni la redistribuye.

## Construcción

Desde la raíz del proyecto:

```bash
docker build \
  -f docker/benchmark.Dockerfile \
  -t dna-benchmark:cpu \
  .
```

El Dockerfile acepta versiones opcionales. Si se dejan vacías, el solver elige
la última versión disponible; **no se garantiza una versión concreta**:

```bash
docker build \
  -f docker/benchmark.Dockerfile \
  -t dna-benchmark:cpu \
  --build-arg BAKTA_VERSION=<version> \
  --build-arg PROKKA_VERSION=<version> \
  --build-arg DFAST_VERSION=<version> \
  .
```

Antes de usar versiones fijas, comprobar que existan y resuelvan en los canales
actuales. El repositorio no fija ni inventa números de versión.

## Ejecución

Montar el repositorio en `/work` para que el script escriba los reportes en
`reports/`. El contenedor se ejecuta como root por defecto; si se usa
`--user "$(id -u):$(id -g)"`, verificar que el usuario pueda escribir en
`reports/` y en el directorio temporal.

### Prokka

```bash
docker run --rm \
  -v "$PWD:/work" \
  dna-benchmark:cpu \
  python /work/tools/benchmark_external.py \
    --genome /work/data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
    --reference-gff /work/data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
    --tool prokka --run \
    --md /work/reports/external_benchmark_prokka.md \
    --json /work/reports/external_benchmark_prokka.json
```

### DFAST

```bash
docker run --rm \
  -v "$PWD:/work" \
  dna-benchmark:cpu \
  python /work/tools/benchmark_external.py \
    --genome /work/data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
    --reference-gff /work/data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
    --tool dfast --run \
    --md /work/reports/external_benchmark_dfast.md \
    --json /work/reports/external_benchmark_dfast.json
```

### Bakta (requiere base local)

La base no está en la imagen. Se monta en `/data/bakta_db` y se pasa con
`--bakta-db`:

```bash
docker run --rm \
  -v "$PWD:/work" \
  -v "/ruta/host/bakta_db:/data/bakta_db:ro" \
  dna-benchmark:cpu \
  python /work/tools/benchmark_external.py \
    --genome /work/data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
    --reference-gff /work/data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
    --tool bakta --run --bakta-db /data/bakta_db \
    --md /work/reports/external_benchmark_bakta.md \
    --json /work/reports/external_benchmark_bakta.json
```

### Las tres a la vez

```bash
docker run --rm \
  -v "$PWD:/work" \
  -v "/ruta/host/bakta_db:/data/bakta_db:ro" \
  dna-benchmark:cpu \
  python /work/tools/benchmark_external.py \
    --genome /work/data/raw/GCF_000006985.1_ASM698v1_genomic.fna.gz \
    --reference-gff /work/data/annotations/GCF_000006985.1_ASM698v1.gff.gz \
    --tool all --run --bakta-db /data/bakta_db \
    --md /work/reports/external_benchmark.md \
    --json /work/reports/external_benchmark.json
```

Si Bakta falla por la base, Prokka/DFAST siguen su curso en la misma llamada;
el reporte conserva el error de Bakta y las métricas de las otras dos.

Para comprobar la imagen sin ejecutar herramientas:

```bash
docker run --rm dna-benchmark:cpu \
  python /opt/benchmark/tools/benchmark_external.py --help
docker run --rm dna-benchmark:cpu \
  sh -lc 'which bakta; which prokka; which dfast'
```

## Qué mide el benchmark dentro del contenedor

`tools/benchmark_external.py` compara el GFF de referencia con el GFF que
produzca cada herramienta:

- CDS, tRNA y rRNA por `seqid`, coordenadas 1-based inclusive y hebra;
- emparejamiento 1:1 por solapamiento recíproco (0,5 por defecto);
- coincidencia exacta de coordenadas;
- atributos `product`, `gene` y `locus_tag`;
- `% hypothetical` cuando el GFF trae `product`;
- tiempo de pared y RSS máximo de la ejecución externa.

No evalúa función, dominios, rutas ni expresa que una herramienta sea superior.
La referencia GFF es el gold de la comparación; sus criterios y errores se
heredan.

## Licencias y datos

- **Bakta:** GPL-3.0 según `docs/BENCHMARK.md`. La base Bakta tiene sus propios
  términos; el usuario debe obtenerla y respetarlos.
- **Prokka:** GPL-3.0 según `docs/BENCHMARK.md`. Sus bases y dependencias tienen
  licencias propias.
- **DFAST:** `docs/BENCHMARK.md` lo registra como "sin licencia declarada".
  Confirmar la licencia y los términos de sus bases en el proyecto upstream
  antes de redistribuir la imagen o sus datos.
- Las dependencias de `conda-forge`/`bioconda` conservan sus licencias.
- La imagen no redistribuye bases de datos ni genomas.
- Construir y distribuir la imagen es responsabilidad de quien la use; revisar
  las licencias de cada paquete y base.

## Límites

1. **No hay build ejecutado aquí.** No se puede afirmar que el solver resuelva,
   que los paquetes existan o que las herramientas corran en todas las
   arquitecturas.
2. **Versiones no pineadas.** Los `ARG` vacíos usan la última versión que
   resuelva el canal; para reproducibilidad estricta hay que fijar versiones y
   verificar que existan.
3. **Bases de datos fuera de la imagen.** Bakta necesita una base local;
   Prokka/DFAST pueden requerir bases o dependencias adicionales según la
   versión. Este Dockerfile no las descarga.
4. **El build requiere red.** La ejecución del benchmark no descarga datos,
   pero una herramienta de terceros podría intentar usar servicios externos;
   `--run` es una decisión explícita.
5. **CPU solamente.** No hay soporte GPU/CUDA.
6. **Sin gold funcional.** Las métricas son de features contra un GFF de
   referencia; no de función, dominios ni rutas.
7. **No hay comparación pareada ni conclusiones de superioridad.** Para eso
   faltan un gold independiente y ejecuciones controladas.

## Pruebas del entorno declarado

```bash
.venv/bin/python -m pytest tests/test_benchmark_env.py -q
```

Estas pruebas solo validan que el Dockerfile y este documento existen y
contienen los comandos/versiones declarados. No construyen la imagen, no usan
red y no ejecutan Bakta, Prokka ni DFAST.
