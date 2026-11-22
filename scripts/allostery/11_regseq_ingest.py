#!/usr/bin/env python3
"""Capa 1 (paso 1) — ingesta y caracterización de las matrices de energía de Reg-Seq.

Qué hace:
  1. Lista las 1.000 matrices de `data/energy_matrices/` vía API de GitHub.
  2. Las descarga (caché local en data/allostery/work/regseq/matrices/).
  3. Parsea el formato `pos val_A val_C val_G val_T` (ΔΔG por base, en unidades k_BT
     según el README del repo).
  4. Caracteriza cada matriz: longitud, amplitud máxima, consenso (base de mínima energía)
     y el ΔΔG del consenso frente a la media de bases.
  5. Mapea promotor → TF usando RegulonDB (TFs que regulan ese gen), y marca ambigüedad.

NO calcula aún afinidades de sitios reales ni valida contra RegulonDB: eso es el paso 2.

Uso:
    python scripts/allostery/11_regseq_ingest.py [--limit N]
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "data" / "allostery" / "work" / "regseq"
MAT = WORK / "matrices"
RD = ROOT / "data" / "regulondb"
OUT = ROOT / "data" / "allostery"

API = "https://api.github.com/repos/RPGroup-PBoC/RegSeq/contents/data/energy_matrices"
RAW = "https://raw.githubusercontent.com/RPGroup-PBoC/RegSeq/master/data/energy_matrices/"


def get(url: str, timeout: int = 30) -> bytes:
    req = Request(url, headers={"User-Agent": "alloreg"})
    with urlopen(req, timeout=timeout) as fh:
        return fh.read()


def list_matrices() -> list[dict]:
    return json.loads(get(API))


def fetch_one(item: dict) -> tuple[str, int]:
    name = item["name"]
    dst = MAT / name
    if dst.exists() and dst.stat().st_size > 0:
        return name, dst.stat().st_size
    data = get(RAW + name, timeout=60)
    dst.write_bytes(data)
    return name, len(data)


def parse_matrix(path: Path) -> dict | None:
    """Formato: cabecera 'pos val_A val_C val_G val_T' y una fila por posición."""
    lines = [l for l in path.read_text(errors="ignore").splitlines() if l.strip()]
    if not lines:
        return None
    header = lines[0].split()
    if header[:1] != ["pos"] or len(header) < 5:
        return None
    pos, vals = [], []
    for line in lines[1:]:
        f = line.split()
        if len(f) < 5:
            continue
        try:
            pos.append(int(f[0]))
            vals.append([float(x) for x in f[1:5]])
        except ValueError:
            continue
    if not vals:
        return None
    n = len(vals)
    cons = []
    amp = 0.0
    for row in vals:
        i = min(range(4), key=lambda k: row[k])
        cons.append("ACGT"[i])
        amp = max(amp, max(row) - min(row))
    return {"n_pos": n, "from": pos[0], "to": pos[-1],
            "amplitud_max_kT": round(amp, 3),
            "consenso": "".join(cons)}


def load_regulondb_tfs():
    """gen -> TFs que lo regulan (para mapear promotor -> TF)."""
    tfs = collections.defaultdict(set)
    for r in csv.DictReader(open(RD / "NetworkRegulatorGene.tsv"), delimiter="\t"):
        g = (r["5)regulatedGeneName"] or "").strip()
        if g:
            tfs[g.lower()].add(r["2)regulatorName"])
    return tfs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    MAT.mkdir(parents=True, exist_ok=True)
    items = list_matrices()
    if args.limit:
        items = items[:args.limit]
    print(f"matrices listadas: {len(items)}", file=sys.stderr)

    sizes = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        for name, size in ex.map(fetch_one, items):
            sizes[name] = size
    total_kb = sum(sizes.values()) / 1024
    print(f"descargadas: {len(sizes)} · {total_kb:.0f} KB", file=sys.stderr)

    tfs = load_regulondb_tfs()
    rows = []
    for item in items:
        name = item["name"]
        m = parse_matrix(MAT / name)
        if not m:
            continue
        # nombre: <promotor>_<condicion>_<desde>_<hasta>_energymatrix
        stem = re.sub(r"_energymatrix$", "", name)
        parts = stem.split("_")
        promoter = parts[0]
        cond = "_".join(parts[1:-2]) if len(parts) > 3 else ""
        cand = sorted(tfs.get(promoter.lower(), ()))
        rows.append({"archivo": name, "promotor": promoter, "condicion": cond,
                     **m, "tfs_candidatos": cand,
                     "n_tfs": len(cand)})

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "regseq_matrices.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "archivo", "promotor", "condicion", "n_pos", "from", "to",
            "amplitud_max_kT", "consenso", "n_tfs", "tfs_candidatos"])
        w.writeheader()
        for r in rows:
            w.writerow({**r, "tfs_candidatos": "|".join(r["tfs_candidatos"])})

    prom = collections.Counter(r["promotor"] for r in rows)
    mapped = [r for r in rows if r["n_tfs"] == 1]
    amb = [r for r in rows if r["n_tfs"] > 1]
    none = [r for r in rows if r["n_tfs"] == 0]
    amps = sorted(r["amplitud_max_kT"] for r in rows)
    summary = {
        "n_matrices": len(rows),
        "n_promotores": len(prom),
        "tamano_total_kb": round(total_kb, 1),
        "amplitud_mediana_kT": amps[len(amps) // 2] if amps else None,
        "amplitud_p90_kT": amps[int(len(amps) * 0.9)] if amps else None,
        "amplitud_maxima_kT": amps[-1] if amps else None,
        "mapeo_tf_unico": len(mapped),
        "mapeo_ambiguo": len(amb),
        "sin_tf_en_regulondb": len(none),
        "top_promotores": prom.most_common(10),
    }
    (OUT / "regseq_ingest.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
