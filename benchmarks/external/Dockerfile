# syntax=docker/dockerfile:1.7
#
# Imagen CPU reproducible para benchmarks/external/external_benchmark.py.
# Incluye Bakta, Prokka y DFAST en entornos separados. No contiene bases de
# datos y no descarga bases durante el build. La base Bakta se monta en /data.
#
# Construcción:
#   docker build -f benchmarks/external/Dockerfile -t dna-token-bac/benchmark-external:cpu .
#
# Versiones fijadas por defecto (Bioconda, consultadas 2026-11-24). Se pueden
# sobreescribir con --build-arg si el canal deja de servirlas.

FROM mambaorg/micromamba:1.5.10

ARG BAKTA_VERSION=1.12.2
ARG PROKKA_VERSION=1.15.6
ARG DFAST_VERSION=1.6.2
ARG PYTHON_VERSION=3.11

LABEL org.opencontainers.image.title="DNA_TOKEN_BAC external benchmark (CPU)" \
      org.opencontainers.image.description="CPU-only image with Bakta, Prokka and DFAST for a reproducible GFF/GBFF benchmark; databases are not included." \
      benchmark.tools="bakta,prokka,dfast" \
      benchmark.cpu_only="true"

USER root

ENV PYTHONUNBUFFERED=1 \
    MAMBA_ROOT_PREFIX=/opt/conda \
    BAKTA_DB=/data/bakta_db

# Python del runner + Biopython para leer referencias GenBank/GBFF.
RUN micromamba install -y -n base -c conda-forge \
        "python=${PYTHON_VERSION}" \
        "biopython>=1.81" \
    && micromamba clean --all --yes

# Un entorno por herramienta: Bakta, Prokka y DFAST tienen dependencias que
# pueden entrar en conflicto si se instalan juntas.
RUN set -eux; \
    micromamba create -y -n bakta -c conda-forge -c bioconda "bakta=${BAKTA_VERSION}"; \
    micromamba clean --all --yes

RUN set -eux; \
    micromamba create -y -n prokka -c conda-forge -c bioconda "prokka=${PROKKA_VERSION}"; \
    micromamba clean --all --yes

RUN set -eux; \
    micromamba create -y -n dfast -c conda-forge -c bioconda "dfast=${DFAST_VERSION}"; \
    micromamba clean --all --yes

# Envoltorios visibles en PATH: el runner detecta los binarios con
# shutil.which y cada uno se ejecuta dentro de su entorno micromamba.
RUN set -eux; \
    MAMBA_BIN="$(command -v micromamba)"; \
    printf '%s\n' '#!/bin/sh' "exec ${MAMBA_BIN} run -n bakta bakta \"\$@\"" > /usr/local/bin/bakta; \
    printf '%s\n' '#!/bin/sh' "exec ${MAMBA_BIN} run -n prokka prokka \"\$@\"" > /usr/local/bin/prokka; \
    printf '%s\n' '#!/bin/sh' "exec ${MAMBA_BIN} run -n dfast dfast \"\$@\"" > /usr/local/bin/dfast; \
    chmod 0755 /usr/local/bin/bakta /usr/local/bin/prokka /usr/local/bin/dfast

# El ENTRYPOINT final invoca `python` directamente. Al sustituir el entrypoint de
# mambaorg/micromamba ya no se ejecuta su script de activación, así que hay que
# exponer el Python del entorno base en PATH.
ENV PATH=/opt/conda/bin:$PATH

# El runner también queda dentro de la imagen para poder ejecutar --help sin
# montar el repositorio.
COPY benchmarks/external/external_benchmark.py /opt/benchmark/external/external_benchmark.py
COPY docs/EXTERNAL_BENCHMARK_RUN.md /opt/benchmark/docs/EXTERNAL_BENCHMARK_RUN.md

WORKDIR /work
VOLUME ["/data"]

# Sin bases de datos incluidas. Montar /data/bakta_db y usar --bakta-db (o la
# variable BAKTA_DB) cuando se ejecute Bakta. No se descarga nada en el build.
ENTRYPOINT ["python", "/opt/benchmark/external/external_benchmark.py"]
CMD ["--help"]
