#!/usr/bin/env bash
# Build de la receta Bioconda. conda-build lo ejecuta con cwd = raíz del
# código fuente (source.path), no dentro del directorio de la receta.
set -euo pipefail

PYTHON_BIN="${PYTHON:-python}"

"${PYTHON_BIN}" -m pip install . --no-deps --no-build-isolation -vv

# webapp/annotator.py importa el paquete `ghost` desde scripts/ghost, pero
# pyproject.toml solo empaqueta `webapp*` y `modules*`. Copiamos ese paquete
# a site-packages para que `dna-annotate` y `python -m webapp.pipeline`
# funcionen fuera del checkout. Es un parche local hasta que el empaquetado
# de pyproject.toml incluya `ghost`.
if [[ -n "${SP_DIR:-}" ]]; then
    SITE_DIR="${SP_DIR}"
else
    SITE_DIR="$("${PYTHON_BIN}" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
fi

if [[ ! -d scripts/ghost/ghost ]]; then
    echo "ERROR: falta scripts/ghost/ghost; no se puede instalar el motor ghost" >&2
    exit 1
fi

mkdir -p "${SITE_DIR}/ghost"
find scripts/ghost/ghost -maxdepth 1 -type f -name '*.py' \
    -exec cp -p {} "${SITE_DIR}/ghost/" \;
