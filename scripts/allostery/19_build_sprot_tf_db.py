#!/usr/bin/env python3
"""Construye una DB DIAMOND de TFs de Swiss-Prot para identificar homólogos.

El FASTA completo de Swiss-Prot (`data/db/uniprot_sprot.fasta`, 288 MB) trae ~570.000
entradas. Para este pipeline solo interesan las que son factores de transcripción o
represores: se filtran por descripción y se construye una DB pequeña.

Salida:
    data/allostery/work/sprot_tf.faa    proteínas filtradas
    data/allostery/work/sprot_tf.dmnd   DB DIAMOND

Uso: python scripts/allostery/19_build_sprot_tf_db.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "db" / "uniprot_sprot.fasta"
WORK = ROOT / "data" / "allostery" / "work"
FAA = WORK / "sprot_tf.faa"
DMND = WORK / "sprot_tf.dmnd"

KEY = re.compile(r"(transcription factor|transcriptional repressor|"
                 r"transcriptional regulator|DNA-binding repressor|HTH-type)", re.I)


def main():
    if not SRC.exists():
        raise SystemExit(f"falta {SRC}")
    WORK.mkdir(parents=True, exist_ok=True)
    keep, buf, header, n_in, n_out = [], [], None, 0, 0
    with open(SRC, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if header and KEY.search(header):
                    keep.append("".join(buf))
                n_in += 1
                header = line[1:].strip()
                buf = [line]
            else:
                buf.append(line)
        if header and KEY.search(header):
            keep.append("".join(buf))
    FAA.write_text("".join(keep))
    n_out = len(keep)
    print(f"entradas Swiss-Prot: {n_in} · TFs conservados: {n_out}", file=sys.stderr)
    print(f"tamaño filtrado: {FAA.stat().st_size/1e6:.1f} MB", file=sys.stderr)
    r = subprocess.run(["diamond", "makedb", "--in", str(FAA), "-d", str(DMND), "--quiet"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(r.stderr[-800:])
    print(f"DB: {DMND} ({DMND.stat().st_size/1e6:.1f} MB)", file=sys.stderr)


if __name__ == "__main__":
    main()
