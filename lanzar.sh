#!/usr/bin/env bash
# lanzar.sh — abre la interfaz web del proyecto con doble clic.
#
# Resuelve los tropiezos habituales de un primer arranque:
#   · si no hay entorno virtual, lo crea e instala requirements.txt
#   · si el puerto 8501 está ocupado, busca el siguiente libre
#   · si la app ya está corriendo, solo abre el navegador
#   · si algo falla, muestra el motivo y NO cierra la ventana
#
# Uso:  doble clic, o `bash lanzar.sh`
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT" || { echo "No encuentro el directorio del proyecto"; read -rp "Enter para cerrar…"; exit 1; }

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }
pause() { read -rp "Pulsa Enter para cerrar…" _ || true; }

open_browser() {
  local url=$1
  for cmd in xdg-open open sensible-browser x-www-browser; do
    command -v "$cmd" >/dev/null 2>&1 && { "$cmd" "$url" >/dev/null 2>&1 & return 0; }
  done
  say "Abre esta dirección en tu navegador: $url"
}

say "Anotador de ADN bacteriano"
echo "Proyecto: $ROOT"

# ---------------------------------------------------------------- intérprete
PY="$ROOT/.venv/bin/python"
REQ="$ROOT/requirements-annotator.txt"
[ -f "$REQ" ] || REQ="$ROOT/requirements.txt"
if [ ! -x "$PY" ]; then
  echo "No hay entorno virtual; creándolo (solo la primera vez)…"
  if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: no hay python3 instalado."; pause; exit 1
  fi
  python3 -m venv .venv || { echo "ERROR al crear el entorno virtual."; pause; exit 1; }
  PY="$ROOT/.venv/bin/python"
  "$PY" -m pip install --quiet --upgrade pip || true
  "$PY" -m pip install --quiet -r "$REQ" || {
    echo "ERROR instalando dependencias (¿sin red?)."; pause; exit 1; }
fi

if ! "$PY" -c "import streamlit" 2>/dev/null; then
  echo "Falta Streamlit; instalando dependencias…"
  "$PY" -m pip install --quiet -r "$REQ" || { echo "ERROR instalando."; pause; exit 1; }
fi

# -------------------------------------------------------------- puerto libre
PORT=8501
if curl -s -m 2 -o /dev/null "http://127.0.0.1:${PORT}/" 2>/dev/null; then
  say "La interfaz ya está corriendo en http://127.0.0.1:${PORT}"
  open_browser "http://127.0.0.1:${PORT}"
  exit 0
fi
for p in 8501 8502 8503 8504 8505 8506 8507 8508; do
  if ! (echo > "/dev/tcp/127.0.0.1/$p") >/dev/null 2>&1; then PORT=$p; break; fi
done

URL="http://127.0.0.1:${PORT}"
LOG="$ROOT/logs/gui.log"
mkdir -p logs

say "Arrancando en ${URL}"
"$PY" -m streamlit run webapp/app.py \
  --server.headless true --server.port "$PORT" \
  --browser.gatherUsageStats false > "$LOG" 2>&1 &
GUI_PID=$!

# esperar a que responda (máximo 40 s)
for _ in $(seq 1 40); do
  sleep 1
  if curl -s -m 1 -o /dev/null "$URL" 2>/dev/null; then
    open_browser "$URL"
    say "Listo: ${URL}   (cierra esta ventana para detener la app)"
    wait "$GUI_PID"
    exit 0
  fi
  if ! kill -0 "$GUI_PID" 2>/dev/null; then break; fi
done

echo "ERROR: la interfaz no respondió. Últimas líneas del registro:"
tail -20 "$LOG" 2>/dev/null
pause
kill "$GUI_PID" 2>/dev/null
exit 1
