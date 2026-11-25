# Estado del benchmark externo CPU (Bakta / Prokka / DFAST)

**Generado:** 2026-11-25T16:09Z  
**Estado global:** `partial` — Prokka con métricas reales; Bakta `not_available`; DFAST `failed`.  
**Gold:** `benchmarks/external/gold/chlorobaculum_tepidum_60kb/` (55 CDS + 1 tRNA; 60 000 bp).  
**Ejecución:** Docker CPU, `--run --threads 1 --timeout 1800`, sin GPU y sin descargar bases grandes.

> No se inventaron métricas. Una herramienta ausente o fallida no es un cero:
> queda con su estado y su razón, sin F1.

## Resultado real

| herramienta | estado | ¿métricas? | detalle |
|---|---|---|---|
| Bakta 1.12.2 | `not_available` | no | falta la Bakta DB; no se ejecutó |
| Prokka 1.15.6 | `ok` | sí | F1 global **0.9643** (54 TP / 2 FP / 2 FN) en 35.2 s |
| DFAST 1.6.2 | `failed` | no | returncode=1 por falta de su base de referencia; no produjo anotación |

### Prokka (ejecución real)

| métrica | valor |
|---|---:|
| features de referencia | 56 (55 CDS + 1 tRNA) |
| TP / FP / FN | 54 / 2 / 2 |
| precisión / recall / F1 | 0.9643 / 0.9643 / **0.9643** |
| F1 macro | 0.9818 |
| coincidencias exactas | 45 |
| pares con borde distinto | 9 |
| F1 CDS / tRNA / rRNA | 0.9636 / 1.0000 / — |
| runtime / RSS | 35.168 s / 126.4 MiB |
| SHA256 del GFF | `bd7a273938250c0807b85058b9650fff81d9d578a17d14f30ad453e319c28e0e` |

El GFF comparado está en
`reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence/prokka/benchmark.gff`.

### Bakta — bloqueo

```
bakta: not_available
razón: Bakta DB no utilizable (no se encontró Bakta DB; definir BAKTA_DB o pasar
--bakta-db (el benchmark no la descarga)); no se ejecuta Bakta y no se reportan métricas
```

En el host no existen `~/.bakta_db` ni `/data/bakta_db`, y `BAKTA_DB` no estaba definida.
La imagen no incluye la base y **no se descargó**.

### DFAST — bloqueo

DFAST sí arrancó, pero abortó antes de anotar:

```
rRNA model file not found:
/opt/conda/envs/dfast/opt/dfast-1.6.2/db/rrna/Rfam_15.1_bac.hmm.
Run 'dfast_file_downloader.py --rrna rfam'.
```

`returncode=1`; no hay GFF/GBFF y por tanto no hay métricas DFAST.
La imagen solo trae `db/hmm_search` (208 KB), no el DB root completo de DFAST.

## Entorno y comando exacto

- Host: Docker **29.1.3**; `which bakta/prokka/dfast` no encuentra binarios en el host.
- Imagen: `dna-token-bac/benchmark-external:cpu`
  (`sha256:d039458ad3a5f78b2928b85e6322dfc1d219c4d77ae2782eafc64996d4b64a2f`, 8.58 GB).
- Versiones dentro de la imagen: Bakta 1.12.2, Prokka 1.15.6, DFAST 1.6.2.
- Build usado:

```bash
DOCKER_CONFIG="$PWD/.scratch/docker-config" \
docker build --progress=plain \
  -f benchmarks/external/Dockerfile \
  -t dna-token-bac/benchmark-external:cpu .
```

`DOCKER_CONFIG` apunta a un directorio escribible del workspace porque `~/.docker`
es de solo lectura en este sandbox. En un host normal basta `docker build ...`.

- Ejecución usada:

```bash
docker run --rm -v "$PWD:/work" -w /work \
  dna-token-bac/benchmark-external:cpu \
  --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools all --run --threads 1 --timeout 1800 \
  --save-gffs /work/reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence \
  --json /work/reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.json \
  --md /work/reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.md
```

Salida del runner:

```
bakta: not_available F1=None fuente=None
prokka: ok F1=0.9643 fuente=executed
dfast: failed F1=None fuente=None
```

## Correcciones aplicadas al harness (con test)

1. **`benchmarks/external/external_benchmark.py`**: DFAST no recibía `--force`.
   El runner crea el directorio de salida antes de invocar la herramienta y
   DFAST abortaba con `Output directory ... already exists`. Añadido `--force`.
   Test: `tests/test_external_benchmark.py::test_build_command_is_cpu_single_thread`.
2. **`benchmarks/external/Dockerfile`**: el `ENTRYPOINT ["python", ...]` fallaba
   con `exec: "python": executable file not found in $PATH` porque al sustituir
   el entrypoint de `mambaorg/micromamba` no se ejecuta su script de activación.
   Añadido `ENV PATH=/opt/conda/bin:$PATH` al final; `docker run ... --help`
   vuelve a funcionar. Test:
   `tests/test_external_benchmark.py::test_dockerfile_exists_and_is_cpu_only`.

Pruebas ejecutadas:

```bash
.venv/bin/python -m pytest tests/test_external_benchmark.py tests/test_external_benchmark_real.py
# 49 passed, 1 warning
```

## Qué falta para métricas de las tres herramientas

1. **Bakta DB local** (varios GB, no descargada): montarla y pasarla al runner.

```bash
export BAKTA_DB=/ruta/host/bakta_db
export DFAST_DB=/ruta/host/dfast_db
docker run --rm -w /work \
  -v "$PWD:/work" \
  -v "$BAKTA_DB:/data/bakta_db:ro" \
  -v "$DFAST_DB:/data/dfast_db:ro" \
  dna-token-bac/benchmark-external:cpu \
  --genome /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/genome.fna \
  --reference /work/benchmarks/external/gold/chlorobaculum_tepidum_60kb/reference.gff3 \
  --tools bakta,prokka,dfast --run --threads 1 --bakta-db /data/bakta_db \
  --extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db
```

2. **DFAST DB root** ya poblado (Rfam rRNA, CDD, NCBIfam, etc.; no descargado):

```bash
# montar el DB root en /data/dfast_db y pasarlo con --dbroot
--extra-arg dfast=--dbroot --extra-arg dfast=/data/dfast_db
```

Sin esas dos bases, el benchmark completo de tres herramientas no puede producir
métricas reales. Prokka ya quedó medido; Bakta y DFAST siguen bloqueados por datos,
no por el código del runner.

## Artefactos

- `reports/external_benchmark/status.json` — resumen máquina-legible.
- `reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.json` — informe completo real.
- `reports/external_benchmark/chlorobaculum_tepidum_60kb/benchmark.md` — informe legible.
- `reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence/prokka/benchmark.gff` — salida Prokka.
- `reports/external_benchmark/chlorobaculum_tepidum_60kb/evidence/dfast/stdout.txt` — error DFAST.
- `reports/external_benchmark/benchmark.json` y `.md` de nivel superior son
  copia del informe real de esta ejecución; el original del gold está bajo
  `chlorobaculum_tepidum_60kb/`.

**Sin commit** (según la instrucción). No se tocaron `app.py`, `browser.py` ni
`data/reference_genomes.json`.
