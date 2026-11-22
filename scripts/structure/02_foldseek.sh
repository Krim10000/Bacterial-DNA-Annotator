#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
ROOT=data/metabolism/targets
mkdir -p "$ROOT/foldseek_hits" /tmp/fs_tmp
if [ ! -d "$ROOT/structures/esmfold" ] || ! ls "$ROOT"/structures/esmfold/*.pdb >/dev/null 2>&1; then
  echo "No hay PDB ESMFold todavía. Ejecuta primero la app Modal."; exit 1
fi
if [ ! -f /tmp/afdb_db ]; then
  foldseek createdb "$ROOT/structures/afdb" /tmp/afdb_db --threads 4
fi
for pdb in "$ROOT"/structures/esmfold/*.pdb; do
  name=$(basename "$pdb" .pdb)
  foldseek easy-search "$pdb" /tmp/afdb_db "$ROOT/foldseek_hits/$name.m8" /tmp/fs_tmp \
    --threads 4 --format-output "query,target,fident,alnlen,evalue,bits"
done
echo "Foldseek listo en $ROOT/foldseek_hits"
