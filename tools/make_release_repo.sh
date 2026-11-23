#!/usr/bin/env bash
# make_release_repo.sh — construye un repositorio limpio para publicar v1.0.
#
#   bash tools/make_release_repo.sh [destino] [--force]
#
# Copia solo código, documentación, pruebas y datos de demostración pequeños.
# Usa una lista explícita: no hay forma de que se cuele data/, un modelo, un
# entorno virtual o un plan interno. Deja un commit inicial en la rama v1.0.
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${1:-$SRC/dist/DNA_TOKEN_BAC_v1.0}"
FORCE=0
[[ "${2:-}" == "--force" ]] && FORCE=1

if [[ -e "$TARGET" && -n "$(ls -A "$TARGET" 2>/dev/null)" && $FORCE -eq 0 ]]; then
  echo "ERROR: $TARGET existe y no está vacío. Usa --force o elige otro destino."
  exit 2
fi

echo "Origen : $SRC"
echo "Destino: $TARGET"
mkdir -p "$TARGET"

FILES=(
  README.md README.es.md README.zh.md LICENSE NOTICE THIRD_PARTY_LICENSES.md DATA_NOTICE.md
  AUTHORS.md CHANGELOG.md CITATION.cff CONTRIBUTING.md CODE_OF_CONDUCT.md
  SECURITY.md VERSION MANIFEST.in pyproject.toml setup.py
  requirements.txt requirements-annotator.txt
  install.sh lanzar.sh run_pipeline.sh config.yaml .gitignore
)
DIRS=(webapp modules scripts docs tools tests assets LICENSES reports .github .streamlit)

for f in "${FILES[@]}"; do
  [[ -f "$SRC/$f" ]] && cp -p "$SRC/$f" "$TARGET/$f"
done

for d in "${DIRS[@]}"; do
  [[ -d "$SRC/$d" ]] || continue
  mkdir -p "$TARGET/$d"
  rsync -a \
    --exclude='__pycache__' --exclude='*.pyc' --exclude='.pytest_cache' \
    --exclude='scripts_orginales' --exclude='_review' \
    --exclude='mocks' --exclude='AGENT_*' --exclude='ACTIVE_WORK.md' \
    --exclude='WORK_PLAN.md' --exclude='MENU_PROPOSALS.md' \
    "$SRC/$d/" "$TARGET/$d/"
done

# config.yaml no debe llevar rutas absolutas del autor
if [[ -f "$TARGET/config.yaml" ]]; then
  sed -i 's|^project_root:.*|project_root: .|' "$TARGET/config.yaml"
fi

# Datos mínimos: solo el README que explica por qué no hay datos.
mkdir -p "$TARGET/data"
cp -p "$SRC/data/README.md" "$TARGET/data/README.md"

# Nada de basura de Python.
find "$TARGET" -type d -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null || true
find "$TARGET" -name '*.pyc' -delete 2>/dev/null || true

cd "$TARGET"
if [[ ! -d .git ]]; then
  git init -q -b v1.0
fi
git add -A
git -c user.name="Camilo Rojas Bravo" \
    -c user.email="Camilo.Rojas@usach.cl" \
    commit -q -m "v1.0.0: anotador bacteriano con GUI, i18n ES/EN/ZH y navegador clicable" \
  || echo "Sin cambios que commitear."

echo
echo "Archivos versionados: $(git ls-files | wc -l)"
echo "Tamaño del árbol   : $(du -sh --exclude=.git . | cut -f1)"
echo "Tamaño de .git     : $(du -sh .git | cut -f1)"
echo
echo "Archivos >50 MB (debe estar vacío):"
git ls-files -z | while IFS= read -r -d '' f; do
  sz=$(du -m "$f" | cut -f1)
  [[ $sz -gt 50 ]] && echo "$sz MB $f"
done
echo
echo "Listo. Para publicar:"
echo "  cd $TARGET"
echo "  git remote add origin git@github.com:<usuario>/<repo>.git"
echo "  git push -u origin v1.0"
