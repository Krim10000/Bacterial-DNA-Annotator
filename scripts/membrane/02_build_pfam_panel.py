#!/usr/bin/env python3
"""Extrae un panel Pfam curado desde Pfam-A.hmm y lo deja presionado para hmmscan.

El Pfam-A.hmm completo (3,6 GB) no esta presionado en este checkout y hmmscan
requiere una base presionada. En vez de presionar 20k modelos, se extrae solo el
panel de membranas/transporte/secrecion definido en pfam_panel_curated.tsv.

Uso:
  python3 scripts/membrane/02_build_pfam_panel.py
"""
from __future__ import annotations
import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def read_panel(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                sys.exit(f"[panel] linea invalida: {line!r}")
            name, role, subgroup, tc, note = parts[:5]
            if role not in {"transporter", "sbp", "secretion", "membrane", "lipoprotein", "signal"}:
                sys.exit(f"[panel] rol invalido {role!r} en {name}")
            rows.append(dict(name=name, role=role, subgroup=subgroup, tc=tc, note=note))
    names = [r["name"] for r in rows]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        sys.exit(f"[panel] nombres duplicados: {sorted(dup)}")
    return rows


def extract(pfam: Path, wanted: set[str], out_hmm: Path) -> dict[str, str]:
    """Streaming de Pfam-A.hmm; escribe solo los registros pedidos."""
    found: dict[str, str] = {}
    buf: list[str] = []
    in_record = False
    name = ""
    acc = ""
    with open(pfam, errors="replace") as fh, open(out_hmm, "w") as out:
        for line in fh:
            if line.startswith("HMMER3/"):
                in_record = True
                buf = [line]
                name = acc = ""
                continue
            if not in_record:
                continue
            buf.append(line)
            if line.startswith("NAME "):
                name = line[5:].strip()
            elif line.startswith("ACC "):
                acc = line.split()[1]
            elif line.startswith("//"):
                if name in wanted:
                    out.write("".join(buf))
                    found[name] = acc
                in_record = False
                buf = []
    return found


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pfam", default=str(ROOT / ".cache/pfam/Pfam-A.hmm"))
    ap.add_argument("--curated", default=str(HERE / "pfam_panel_curated.tsv"))
    ap.add_argument("--outdir", default=str(ROOT / "data/membrane/pfam_panel"))
    args = ap.parse_args()

    pfam = Path(args.pfam)
    curated = read_panel(Path(args.curated))
    wanted = {r["name"] for r in curated}
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_hmm = outdir / "membrane_panel.hmm"

    if not pfam.exists():
        sys.exit(f"[panel] no existe {pfam}")
    print(f"[panel] modelos curados: {len(wanted)}; leyendo {pfam} ...", flush=True)
    found = extract(pfam, wanted, out_hmm)
    missing = sorted(wanted - set(found))
    if missing:
        sys.exit(f"[panel] modelos no encontrados en Pfam-A.hmm: {missing}")
    print(f"[panel] escritos {len(found)} modelos en {out_hmm}", flush=True)

    with open(outdir / "panel_models.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pfam_name", "pfam_acc", "role", "subgroup", "tc", "note"])
        for r in curated:
            w.writerow([r["name"], found[r["name"]], r["role"], r["subgroup"], r["tc"], r["note"]])

    print("[panel] hmmpress ...", flush=True)
    subprocess.run(["hmmpress", "-f", str(out_hmm)], check=True)
    print("[panel] OK", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
