#!/usr/bin/env python3
"""02 — Ejecuta fastANI (pilotos x panel de referencia).

Salida cruda:   <WORK>/ani/ani_raw.tsv  (formato fastANI)
Salida tabular: data/taxonomy/ani_hits.tsv
  genome <TAB> ref <TAB> ani <TAB> orthologous_matches <TAB> total_fragments <TAB> af

`af` = fraccion de fragmentos de la query con match (proxy de cobertura), tal
como la calcula fastANI. No se transforma ni se suaviza ningun valor.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402


def read_manifest() -> dict[str, dict]:
    rows: dict[str, dict] = {}
    with open(L.OUT / "refpanel.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["genome"]] = r
    return rows


def main() -> int:
    L.ensure_dirs()
    panel = read_manifest()
    pilots = L.read_pilot_tsv()
    pilot_keys = {p["genome"] for p in pilots}

    qlist = L.Q_DIR / "queries.txt"
    rlist = L.Q_DIR / "refs.txt"
    with open(qlist, "w", encoding="utf-8") as fh:
        for p in pilots:
            fh.write(panel[p["genome"]]["fasta_path"] + "\n")
    with open(rlist, "w", encoding="utf-8") as fh:
        for g, r in sorted(panel.items()):
            if g in pilot_keys:  # no comparar un genoma consigo mismo
                continue
            fh.write(r["fasta_path"] + "\n")

    n_refs = sum(1 for g in panel if g not in pilot_keys)
    L.log(f"[ani] {len(pilots)} queries x {n_refs} referencias")

    raw = L.ANI_DIR / "ani_raw.tsv"
    cmd = [
        L.fastani_bin(),
        "--ql", str(qlist),
        "--rl", str(rlist),
        "-o", str(raw),
        "--minFraction", L.FASTANI_MIN_FRACTION,
        "-t", os.environ.get("TAXONOMY_THREADS", "8"),
    ]
    L.log("[ani] " + " ".join(cmd))
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        L.log(f"[ani] fastANI fallo con rc={proc.returncode}")
        return proc.returncode

    # fastANI imprime rutas absolutas; las devolvemos a claves de genoma
    path_to_key = {str(Path(r["fasta_path"]).resolve()): g for g, r in panel.items()}
    out_rows = []
    if raw.exists():
        with open(raw, encoding="utf-8") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                q, r, ani, matches, total = p[0], p[1], p[2], p[3], p[4]
                qk = path_to_key.get(str(Path(q).resolve()), Path(q).stem)
                rk = path_to_key.get(str(Path(r).resolve()), Path(r).stem)
                if qk in pilot_keys and rk == qk:
                    continue
                af = (int(matches) / int(total)) if int(total) else 0.0
                out_rows.append(
                    {
                        "genome": qk,
                        "ref": rk,
                        "ani": ani,
                        "orthologous_matches": matches,
                        "total_fragments": total,
                        "af": f"{af:.4f}",
                    }
                )
    out_rows.sort(key=lambda r: (r["genome"], -float(r["ani"])))
    with open(L.OUT / "ani_hits.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["genome", "ref", "ani", "orthologous_matches", "total_fragments", "af"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(out_rows)

    per_genome = {}
    for r in out_rows:
        per_genome.setdefault(r["genome"], []).append(r)
    for g in sorted(pilot_keys):
        hits = per_genome.get(g, [])
        if hits:
            L.log(f"[ani] {g}: {len(hits)} hits, mejor {hits[0]['ref']} {hits[0]['ani']}")
        else:
            L.log(f"[ani] {g}: sin hits con minFraction={L.FASTANI_MIN_FRACTION}")
    L.log(f"[ani] filas -> {L.OUT / 'ani_hits.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
