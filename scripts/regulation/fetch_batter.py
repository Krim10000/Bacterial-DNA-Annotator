#!/usr/bin/env python3
"""fetch_batter.py — BATTER-TPE, detector alternativo de terminadores.

Qué es
------
BATTER (BActeria Transcript Three prime End Recognizer) predice los extremos 3'
de transcritos primarios. BATTER-TPE es su módulo de aprendizaje: un Longformer
de 8 capas con cabeza CRF que puntúa cada posición. Es **independiente de Rho**
en el sentido de que no asume el mecanismo: aprende del contexto, no de un
patrón fijo de tallo + poli-U. Eso es justo lo que le falta a nuestro detector
de horquillas, que solo ve estructuras intrínsecas.

Fuente: https://github.com/lulab/BATTER  (Zenodo: 10.5281/zenodo.18863501)
Preprint: bioRxiv 2023.10.02.560326

Licencia: **el repositorio no declara licencia**. Sin licencia, los derechos
quedan reservados: no se puede redistribuir ni incorporar al repositorio. Por
eso este script **descarga a `data/terminators/batter/`** (ignorado por git) en
la máquina de quien lo use, igual que se hace con Pfam o DBTBS. Si el autor
publica una licencia, este comentario debe actualizarse.

Compatibilidad
--------------
El modelo se entrenó con transformers 4.18 y torch antiguo. Con las versiones
actuales fallan dos cosas, que este script corrige sobre la copia descargada:

1. `load_state_dict` estricto: transformers ≥4.31 añade el buffer
   `position_ids`, que no está en el checkpoint -> `strict=False`.
2. `.view(...)` sobre tensores no contiguos en la decodificación CRF ->
   `.reshape(...)`.

Ambos parches están verificados: el modelo carga y produce predicciones.

Uso
---
    python scripts/regulation/fetch_batter.py           # descarga + parchea
    python scripts/regulation/fetch_batter.py --check   # solo verifica

Nota CPU: 32 MB de pesos, ~30 kb/min en un núcleo sin GPU. Un genoma de 4,6 Mb
tarda ~2,5 h; para evaluar se usa una región, no el genoma completo.
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data" / "terminators" / "batter"
RAW = "https://raw.githubusercontent.com/lulab/BATTER/main/"
FILES = [
    "model/batter.mdl.pt", "model/scoring.json", "model/fpr01.pkl",
    "model/TNF-by-clusters.txt", "model/TPE-FPR-by-clusters.txt",
    "model/RUT-FPR-by-clusters.txt",
    "config/64.8.256.4.json",
    "scripts/batter-tpe", "scripts/model.py", "scripts/crf.py",
    "scripts/dataset.py", "scripts/constants.py",
]
UA = {"User-Agent": "Mozilla/5.0"}


def download(rel: str) -> None:
    path = DEST / rel
    if path.exists() and path.stat().st_size > 100:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(RAW + rel, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as fh, open(path, "wb") as out:
        out.write(fh.read())
    print(f"  {rel:38s} {path.stat().st_size / 1e6:6.1f} MB")


def patch() -> None:
    """Compatibilidad con transformers/torch actuales (ver docstring)."""
    tpe = DEST / "scripts" / "batter-tpe"
    s = tpe.read_text()
    if "strict=False" not in s:
        s = s.replace("tagger.load_state_dict(state_dict)",
                      "tagger.load_state_dict(state_dict, strict=False)")
        tpe.write_text(s)
        print("  parche: load_state_dict(strict=False)")

    crf = DEST / "scripts" / "crf.py"
    s = crf.read_text()
    if ".view(batch_size" in s:
        import re
        s = re.sub(r"\.view\(batch_size", ".reshape(batch_size", s)
        crf.write_text(s)
        print("  parche: .view -> .reshape en la decodificación CRF")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    missing = [f for f in FILES if not (DEST / f).exists()]
    if a.check:
        print(f"BATTER en {DEST}: "
              f"{'completo' if not missing else 'faltan ' + str(missing)}")
        return 0 if not missing else 1

    print(f"descargando BATTER en {DEST.relative_to(ROOT)} …")
    for rel in FILES:
        download(rel)
    patch()
    print("listo. Ejemplo:")
    print(f"  cd {DEST.relative_to(ROOT)} && ../../../.venv/bin/python "
          "scripts/batter-tpe -f <genoma.fna> -d cpu -o salida.bed -w 500 -s 500")
    print("licencia: no declarada -> no se versiona (data/ está ignorado)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
