#!/usr/bin/env bash
# modal_guard.sh — control de GPU en Modal (límite: 10 GPU simultáneas por workspace).
#
# Regla del proyecto: todo agente que use Modal debe (1) avisarlo y (2) cerrar lo
# que abrió al terminar. Este script verifica y cierra.
#
#   bash scripts/allostery/modal_guard.sh status     # ¿qué está activo?
#   bash scripts/allostery/modal_guard.sh stop-all   # cerrar todo lo no parado
#
# Salida de `status`: código 0 si no hay nada activo, 1 si hay apps o contenedores.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
PY="python"
[ -x ".venv/bin/python" ] && PY=".venv/bin/python"

cmd="${1:-status}"

case "$cmd" in
  status)
    "$PY" - <<'EOF'
import json, subprocess, sys

def run(*args):
    return subprocess.run(["modal", *args], capture_output=True, text=True).stdout

apps = json.loads(run("app", "list", "--json") or "[]")
active = [a for a in apps if a.get("state") not in ("stopped", "ephemeral_stopped")]
print(f"apps: {len(apps)} total · {len(active)} activas")
for a in active:
    print(f"  {a['app_id']}  {a['state']:12s}  {a.get('description','')}")

try:
    cont = run("container", "list")
    n = sum(1 for line in cont.splitlines() if line.startswith("│") and "Container ID" not in line
            and line.strip("│ ").strip())
    print(f"contenedores activos (aprox.): {n}")
except Exception as e:
    print("contenedores: no consultable", e)

sys.exit(1 if active else 0)
EOF
    ;;
  stop-all)
    "$PY" - <<'EOF'
import json, subprocess

out = subprocess.run(["modal", "app", "list", "--json"], capture_output=True, text=True).stdout
apps = json.loads(out or "[]")
active = [a for a in apps if a.get("state") not in ("stopped", "ephemeral_stopped")]
if not active:
    print("nada que cerrar")
    raise SystemExit(0)
for a in active:
    print(f"cerrando {a['app_id']} ({a.get('description','')})…")
    subprocess.run(["modal", "app", "stop", "-y", a["app_id"]], capture_output=True)
print(f"cerradas {len(active)} apps")
EOF
    ;;
  *)
    echo "uso: $0 {status|stop-all}" >&2
    exit 2
    ;;
esac
