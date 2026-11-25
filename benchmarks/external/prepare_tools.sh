#!/usr/bin/env bash
# Prepara Bakta, Prokka y DFAST en entornos micromamba locales de CPU.
#
# - No descarga bases de datos (Bakta DB, referencias DFAST, etc.).
# - Las versiones están fijadas y se pueden sobreescribir por variable de entorno.
# - Crea envoltorios en $TARGET/bin; añadir ese directorio al PATH.
#
# Uso:
#   bash benchmarks/external/prepare_tools.sh
#   export PATH="$PWD/benchmarks/external/.tools/bin:$PATH"
#   python benchmarks/external/external_benchmark.py --help
#
# Requiere micromamba y red para resolver los paquetes de conda-forge/bioconda.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TARGET="${TARGET:-$ROOT/benchmarks/external/.tools}"
MAMBA="${MAMBA:-micromamba}"

BAKTA_VERSION="${BAKTA_VERSION:-1.12.2}"
PROKKA_VERSION="${PROKKA_VERSION:-1.15.6}"
DFAST_VERSION="${DFAST_VERSION:-1.6.2}"

if ! command -v "$MAMBA" >/dev/null 2>&1; then
    echo "error: no se encontró '$MAMBA' en PATH" >&2
    echo "instalar micromamba o definir MAMBA=/ruta/a/micromamba" >&2
    exit 2
fi

mkdir -p "$TARGET/envs" "$TARGET/bin"

create_env() {
    local tool="$1"
    local version="$2"
    echo "[prepare] creando entorno $tool=$version en $TARGET/envs/$tool"
    "$MAMBA" create -y -p "$TARGET/envs/$tool" \
        -c conda-forge -c bioconda "${tool}=${version}"
    "$MAMBA" clean --all --yes
}

create_env bakta "$BAKTA_VERSION"
create_env prokka "$PROKKA_VERSION"
create_env dfast "$DFAST_VERSION"

for tool in bakta prokka dfast; do
    cat > "$TARGET/bin/$tool" <<EOF
#!/bin/sh
exec "$MAMBA" run -p "$TARGET/envs/$tool" $tool "\$@"
EOF
    chmod 0755 "$TARGET/bin/$tool"
done

echo "[prepare] listo. Añadir al PATH:"
echo "  export PATH=\"$TARGET/bin:\$PATH\""
echo "[prepare] no se descargó ninguna base de datos."
