#!/usr/bin/env python3
"""Campaña: correr el pipeline de represores sobre N genomas y agregar resultados.

Entrada: archivos .fna/.fna.gz (por defecto, los más pequeños de data/raw).
Salida:  data/allostery/campaign/
           <genoma>/represores.tsv, summary.json
           agregado.tsv    una fila por genoma
           informe.md      tabla + familias más frecuentes

Uso: python scripts/allostery/22_campaign.py --n 6
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "data" / "allostery" / "work"
CAMP = ROOT / "data" / "allostery" / "campaign"
PIPE = ROOT / "scripts" / "allostery" / "18_genome_pipeline.py"
PY = str(ROOT / ".venv" / "bin" / "python")


def pick(n: int) -> list[Path]:
    files = sorted((ROOT / "data" / "raw").glob("*.fna.gz"), key=lambda p: p.stat().st_size)
    return files[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--glob", default="")
    args = ap.parse_args()
    CAMP.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    files = sorted((ROOT / "data" / "raw").glob(args.glob)) if args.glob else pick(args.n)

    rows = []
    for i, f in enumerate(files, 1):
        name = f.name.replace(".fna.gz", "")
        fa = WORK / f"{name}.fna"
        with gzip.open(f, "rb") as src, open(fa, "wb") as dst:
            shutil.copyfileobj(src, dst)
        out = CAMP / name
        print(f"[{i}/{len(files)}] {name}…", file=sys.stderr, flush=True)
        r = subprocess.run([PY, str(PIPE), str(fa), "--out", str(out),
                            "--find-repressors"], capture_output=True, text=True)
        if r.returncode != 0:
            print("  falló:", r.stderr[-300:], file=sys.stderr)
            continue
        s = json.loads((out / "summary.json").read_text())
        rows.append({"genoma": name, "bp": s["bp"], "genes": s["genes"],
                     "candidatos_tf": s["tf_candidatos"], "represores": s["reportados"],
                     "con_efector": s["con_efector_heredado"],
                     "tiempo_s": s["recursos"]["tiempos_s"]["total_s"],
                     "familias": json.dumps(s["familias"], ensure_ascii=False)})

    with open(CAMP / "agregado.tsv", "w") as fh:
        cols = ["genoma", "bp", "genes", "candidatos_tf", "represores", "con_efector",
                "tiempo_s", "familias"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    tot_fam = Counter()
    for r in rows:
        tot_fam.update(json.loads(r["familias"]))
    lines = ["# Campaña de represores en genomas desconocidos", "",
             f"Genomas procesados: **{len(rows)}**", "",
             "| genoma | Mb | genes | candidatos TF | represores | con efector | tiempo |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['genoma']} | {r['bp']/1e6:.2f} | {r['genes']} | "
                     f"{r['candidatos_tf']} | {r['represores']} | {r['con_efector']} | "
                     f"{r['tiempo_s']} s |")
    lines += ["", "## Familias más frecuentes", ""]
    for fam, n in tot_fam.most_common(12):
        lines.append(f"- {fam}: {n}")
    (CAMP / "informe.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"genomas": len(rows),
                      "represores_total": sum(r["represores"] for r in rows),
                      "tiempo_medio_s": round(sum(r["tiempo_s"] for r in rows) / max(1, len(rows)), 1)},
                     indent=1))


if __name__ == "__main__":
    main()
