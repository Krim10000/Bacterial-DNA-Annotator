# syntax=docker/dockerfile:1
#
# Imagen CPU para ejecutar tools/benchmark_external.py contra Bakta, Prokka y
# DFAST. No se construye ni se ejecuta en este repositorio: es una receta.
#
# - No descarga bases de datos (Bakta DB, Prokka DB, DFAST DB) durante el build.
# - Cada herramienta vive en su propio entorno micromamba para evitar conflictos
#   de dependencias.
# - Los binarios se exponen en /usr/local/bin para que shutil.which los vea.
# - Las versiones quedan sin pinear por defecto; se pueden fijar con build-args.
#
# Construcción:
#   docker build -f docker/benchmark.Dockerfile -t dna-benchmark:cpu .
#
# Si se quieren versiones concretas:
#   docker build -f docker/benchmark.Dockerfile -t dna-benchmark:cpu \
#     --build-arg BAKTA_VERSION=<version> \
#     --build-arg PROKKA_VERSION=<version> \
#     --build-arg DFAST_VERSION=<version> .
#
# No se garantiza que una versión concreta exista o resuelva en el canal actual;
# el usuario debe comprobarlo en Bioconda/conda-forge.

FROM mambaorg/micromamba:latest

# Vacío = última versión que resuelva el solver. No se inventan versiones.
ARG BAKTA_VERSION=""
ARG PROKKA_VERSION=""
ARG DFAST_VERSION=""
ARG PYTHON_VERSION="3.11"

LABEL org.opencontainers.image.title="DNA_TOKEN_BAC external benchmark (CPU)" \
      org.opencontainers.image.description="CPU image with Bakta, Prokka and DFAST for tools/benchmark_external.py; databases are not included." \
      benchmark.tools="bakta,prokka,dfast"

USER root

ENV PYTHONUNBUFFERED=1 \
    MAMBA_ROOT_PREFIX=/opt/conda

# Python para el script de benchmark y pytest para verificaciones opcionales.
RUN micromamba install -y -n base -c conda-forge \
        "python=${PYTHON_VERSION}" \
        pytest \
    && micromamba clean --all --yes

# Entornos separados: Bakta, Prokka y DFAST tienen dependencias que pueden
# entrar en conflicto en un único entorno.
RUN set -eux; \
    if [ -n "${BAKTA_VERSION}" ]; then BAKTA_SPEC="bakta=${BAKTA_VERSION}"; else BAKTA_SPEC="bakta"; fi; \
    micromamba create -y -n bakta -c conda-forge -c bioconda "${BAKTA_SPEC}"; \
    micromamba clean --all --yes

RUN set -eux; \
    if [ -n "${PROKKA_VERSION}" ]; then PROKKA_SPEC="prokka=${PROKKA_VERSION}"; else PROKKA_SPEC="prokka"; fi; \
    micromamba create -y -n prokka -c conda-forge -c bioconda "${PROKKA_SPEC}"; \
    micromamba clean --all --yes

RUN set -eux; \
    if [ -n "${DFAST_VERSION}" ]; then DFAST_SPEC="dfast=${DFAST_VERSION}"; else DFAST_SPEC="dfast"; fi; \
    micromamba create -y -n dfast -c conda-forge -c bioconda "${DFAST_SPEC}"; \
    micromamba clean --all --yes

# Envoltorios: benchmark_external.py detecta los comandos con shutil.which.
# Cada envoltorio ejecuta el binario real dentro de su entorno micromamba.
RUN set -eux; \
    printf '%s\n' '#!/bin/sh' 'exec /usr/local/bin/micromamba run -n bakta bakta "$@"' > /usr/local/bin/bakta; \
    printf '%s\n' '#!/bin/sh' 'exec /usr/local/bin/micromamba run -n prokka prokka "$@"' > /usr/local/bin/prokka; \
    printf '%s\n' '#!/bin/sh' 'exec /usr/local/bin/micromamba run -n dfast dfast "$@"' > /usr/local/bin/dfast; \
    chmod 0755 /usr/local/bin/bakta /usr/local/bin/prokka /usr/local/bin/dfast

# Script de benchmark. El repositorio se monta normalmente en /work; esta copia
# permite ejecutar --help incluso sin montar el proyecto.
COPY tools/benchmark_external.py /opt/benchmark/tools/benchmark_external.py
COPY docs/EXTERNAL_BENCHMARK.md /opt/benchmark/docs/EXTERNAL_BENCHMARK.md

WORKDIR /work
VOLUME ["/data"]

# La imagen no incluye bases de datos. Montar /data y pasar --bakta-db cuando
# se ejecute Bakta. No se descarga nada en el build.
CMD ["python", "/opt/benchmark/tools/benchmark_external.py", "--help"]
