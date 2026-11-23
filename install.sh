#!/usr/bin/env bash
# install.sh — instalador sencillo para Ubuntu (y derivados).
#
#   bash install.sh                 # entorno + lanzador de menú
#   bash install.sh --with-data     # + Swiss-Prot y Rhea (~120 MB)
#   bash install.sh --with-tools    # + tRNAscan-SE, barrnap, HMMER, DIAMOND (sudo)
#   bash install.sh --full          # + todas las dependencias del pipeline
#   bash install.sh --dry-run       # muestra lo que haría, sin tocar nada
#   bash install.sh --uninstall     # quita el lanzador del menú
#
# No usa root para Python: crea un entorno virtual dentro del proyecto. Solo
# pide sudo si se piden las herramientas del sistema con --with-tools.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
VENV="$ROOT/.venv"
PY="$VENV/bin/python"
APP_DIR="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
DESKTOP="$APP_DIR/anotador-adn.desktop"
ICON="$ROOT/assets/icon.svg"
BIN_DIR="$HOME/.local/bin"
LOG="$ROOT/logs/install.log"

DRY=0
WITH_DATA=0
WITH_TOOLS=0
FULL=0
UNINSTALL=0
NO_DESKTOP=0

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY=1 ;;
    --with-data) WITH_DATA=1 ;;
    --with-tools) WITH_TOOLS=1 ;;
    --full) FULL=1 ;;
    --uninstall) UNINSTALL=1 ;;
    --no-desktop) NO_DESKTOP=1 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "Opción desconocida: $arg (usa --help)"; exit 2 ;;
  esac
done

if [[ "${LANG:-}" == es* ]]; then
  T_TITLE="Anotador de ADN bacteriano — instalación"
  T_ROOT="Proyecto"
  T_DRY="MODO PRUEBA: no se modifica nada."
  T_PY="Comprobando Python…"
  T_PY_OK="Python encontrado"
  T_PY_BAD="Se necesita Python 3.10 o superior."
  T_VENV="Creando entorno virtual en .venv (solo la primera vez)…"
  T_PIP="Instalando dependencias de Python…"
  T_DATA="Descargando Swiss-Prot y Rhea (~120 MB)…"
  T_TOOLS="Instalando herramientas del sistema (pedirá sudo)…"
  T_DESK="Creando lanzador en el menú de aplicaciones…"
  T_DONE="Listo. Abre «Anotador de ADN» en Aplicaciones o ejecuta:"
  T_UNINST="Quitando el lanzador del menú…"
  T_MISS="Herramientas no instaladas (la app las desactivará):"
else
  T_TITLE="Bacterial DNA annotator — installation"
  T_ROOT="Project"
  T_DRY="DRY RUN: nothing will be changed."
  T_PY="Checking Python…"
  T_PY_OK="Python found:"
  T_PY_BAD="Python 3.10 or newer is required."
  T_VENV="Creating virtual environment in .venv (first run only)…"
  T_PIP="Installing Python dependencies…"
  T_DATA="Downloading Swiss-Prot and Rhea (~120 MB)…"
  T_TOOLS="Installing system tools (sudo will be requested)…"
  T_DESK="Creating the application-menu launcher…"
  T_DONE="Done. Open “DNA annotator” in Applications, or run:"
  T_UNINST="Removing the menu launcher…"
  T_MISS="Tools not installed (the app will disable them):"
fi

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }
run() {
  if [[ $DRY -eq 1 ]]; then printf '  [dry-run] %s\n' "$*"; else "$@"; fi
}

if [[ $UNINSTALL -eq 1 ]]; then
  say "$T_UNINST"
  run rm -f "$DESKTOP"
  run rm -f "$BIN_DIR/dna-annotator"
  if command -v update-desktop-database >/dev/null 2>&1 && [[ $DRY -eq 0 ]]; then
    update-desktop-database "$APP_DIR" >/dev/null 2>&1 || true
  fi
  echo "OK"
  exit 0
fi

say "$T_TITLE"
echo "$T_ROOT: $ROOT"
[[ $DRY -eq 1 ]] && echo "$T_DRY"

# ------------------------------------------------------------------ Python
say "$T_PY"
if ! command -v python3 >/dev/null 2>&1; then
  echo "$T_PY_BAD"; exit 1
fi
python3 - <<'PY' || { echo "$T_PY_BAD"; exit 1; }
import sys
raise SystemExit(0 if sys.version_info >= (3, 10) else 1)
PY
echo "$T_PY_OK $(python3 --version 2>&1)"

# ------------------------------------------------------------- entorno
if [[ ! -x "$PY" ]]; then
  say "$T_VENV"
  run python3 -m venv "$VENV"
  run "$PY" -m pip install --quiet --upgrade pip
fi
if [[ ! -x "$PY" && $DRY -eq 0 ]]; then
  echo "ERROR: no se pudo crear el entorno virtual."; exit 1
fi

say "$T_PIP"
REQ="$ROOT/requirements-annotator.txt"
[[ $FULL -eq 1 ]] && REQ="$ROOT/requirements.txt"
run "$PY" -m pip install --quiet -r "$REQ"

# ------------------------------------------------------------------ datos
if [[ $WITH_DATA -eq 1 ]]; then
  say "$T_DATA"
  run "$PY" -m webapp.setup_data --data-root "$ROOT/data" || true
fi

# ----------------------------------------------------------- herramientas
if [[ $WITH_TOOLS -eq 1 ]]; then
  say "$T_TOOLS"
  MISSING=()
  for pkg in trnascan-se infernal hmmer diamond-aligner barrnap ncbi-blast+; do
    if dpkg -s "$pkg" >/dev/null 2>&1; then continue; fi
    if [[ $DRY -eq 1 ]]; then
      echo "  [dry-run] sudo apt-get install -y $pkg"
    elif sudo apt-get install -y "$pkg"; then
      :
    else
      MISSING+=("$pkg")
    fi
  done
  if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "$T_MISS ${MISSING[*]}"
  fi
fi

# ------------------------------------------------------------ lanzador
if [[ $NO_DESKTOP -eq 0 ]]; then
  say "$T_DESK"
  if [[ $DRY -eq 0 ]]; then
    mkdir -p "$APP_DIR" "$BIN_DIR" "$ROOT/logs"
    chmod +x "$ROOT/lanzar.sh" 2>/dev/null || true
    cat > "$DESKTOP" <<EOF
[Desktop Entry]
Type=Application
Version=1.0
Name=Anotador de ADN bacteriano
Name[en]=Bacterial DNA Annotator
Name[zh_CN]=细菌 DNA 注释器
Comment=Anota genomas bacterianos y busca brechas metabólicas
Comment[en]=Annotate bacterial genomes and find metabolic gaps
Comment[zh_CN]=注释细菌基因组并寻找代谢缺口
Exec=$ROOT/lanzar.sh
Path=$ROOT
Icon=$ICON
Terminal=false
Categories=Science;Education;Biology;
Keywords=DNA;genome;bioinformatics;annotation;
StartupNotify=true
EOF
    chmod +x "$DESKTOP"
    # comando de terminal opcional
    ln -sf "$ROOT/lanzar.sh" "$BIN_DIR/dna-annotator"
    if command -v update-desktop-database >/dev/null 2>&1; then
      update-desktop-database "$APP_DIR" >/dev/null 2>&1 || true
    fi
  else
    echo "  [dry-run] escribir $DESKTOP y $BIN_DIR/dna-annotator"
  fi
fi

say "$T_DONE"
echo "  $ROOT/lanzar.sh"
echo "  http://localhost:8501"
echo
STAR_URL="https://github.com/Krim10000/Bacterial-DNA-Annotator"
echo "⭐ Si te sirvió, dale una estrella: $STAR_URL"
