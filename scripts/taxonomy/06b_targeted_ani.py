#!/usr/bin/env python3
"""06b — Segunda pasada de ANI, dirigida, para tapar huecos del barrido por filo.

El primer barrido (06) compara cada genoma contra ≤5 referencias de su filo. Eso
deja sin ANI a genomas cuyo pariente cercano no entró en esa muestra. Esta pasada
solo trabaja con:

  * genomas que tienen un **conspecífico** en el panel (aunque ya tuvieran ANI: se
    sustituye por el valor contra su propia especie);
  * genomas **sin ningún hit** que tienen referencias en su **misma familia**.

Referencias = todas las del mismo `species_taxid` (conspecíficas) + todas las de la
misma familia (sin tope). Coste pequeño porque cada grupo tiene pocas queries.

Los hits se anexan al mismo `ani_hits_full.jsonl` con `scope_rank` real, y se
reconsolida el TSV (se conserva el mayor ANI por par genoma-referencia).
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_taxonomy as L  # noqa: E402

FULL = L.OUT / "full"


def main() -> int:
    meta = {r["genome"]: r for r in
            csv.DictReader(open(FULL / "genome_meta.tsv", encoding="utf-8"), delimiter="\t")}
    panel = {r["genome"]: r for r in
             csv.DictReader(open(FULL / "panel_full.tsv", encoding="utf-8"), delimiter="\t")}
    paths = {g: r["fasta_path"] for g, r in panel.items() if r.get("fasta_path")}
    path_to_key = {str(Path(v).resolve()): k for k, v in paths.items()}

    have = set()
    hits_path = FULL / "ani_hits_full.tsv"
    if hits_path.exists():
        with open(hits_path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                have.add(r["genome"])

    by_species = defaultdict(list)
    by_family = defaultdict(list)
    for g, r in panel.items():
        if r.get("species_taxid"):
            by_species[r["species_taxid"]].append(g)
        if r.get("rk_family"):
            by_family[r["rk_family"]].append(g)

    pilots = {r["genome"] for r in csv.DictReader(open(L.PILOT_TSV, encoding="utf-8"),
                                                 delimiter="\t")}
    targets: dict[str, set[str]] = defaultdict(set)  # family -> queries
    reasons: dict[str, str] = {}
    for g, r in meta.items():
        sp, fam = r.get("species_taxid") or "", r.get("family") or ""
        cons = [x for x in by_species.get(sp, []) if x != g]
        if cons:
            targets[fam or f"sp:{sp}"].add(g)
            reasons[g] = "conspecific"
        elif g not in have and fam and [x for x in by_family.get(fam, []) if x != g]:
            targets[fam].add(g)
            reasons[g] = "no_hit_family_refs"
        elif g in pilots and fam and [x for x in by_family.get(fam, []) if x != g]:
            # los 10 pilotos se recalculan contra TODA su familia para que la capa
            # completa no pierda el pariente cercano que sí ve la capa piloto
            targets[fam].add(g)
            reasons[g] = "pilot_family_full"

    L.log(f"[06b] genomas objetivo: {len(reasons)} (conspecíficos "
          f"{sum(1 for v in reasons.values() if v == 'conspecific')}) en "
          f"{len(targets)} grupos")

    appended = 0
    with open(FULL / "ani_hits_full.jsonl", "a", encoding="utf-8") as hf:
        for fam, queries in sorted(targets.items()):
            refs: set[str] = set()
            for q in queries:
                sp = meta[q].get("species_taxid") or ""
                refs.update(x for x in by_species.get(sp, []) if x != q)
            fam_key = fam if not fam.startswith("sp:") else ""
            if fam_key:
                refs.update(x for x in by_family.get(fam_key, []) if x not in queries)
            refs = sorted(r for r in refs if r in paths)
            if not refs:
                continue
            qfile = L.Q_DIR / "tgt_q.txt"
            rfile = L.Q_DIR / "tgt_r.txt"
            qfile.write_text("\n".join(paths[g] for g in queries) + "\n", encoding="utf-8")
            rfile.write_text("\n".join(paths[r] for r in refs) + "\n", encoding="utf-8")
            out = L.ANI_DIR / "tgt_out.tsv"
            out.unlink(missing_ok=True)
            t = time.time()
            proc = subprocess.run(
                [L.fastani_bin(), "--ql", str(qfile), "--rl", str(rfile),
                 "-o", str(out), "--minFraction", L.FASTANI_MIN_FRACTION,
                 "-t", os.environ.get("TAXONOMY_THREADS", "8")],
                capture_output=True, text=True)
            n = 0
            if proc.returncode == 0 and out.exists():
                with open(out, encoding="utf-8") as fh:
                    for line in fh:
                        p = line.rstrip("\n").split("\t")
                        if len(p) < 5:
                            continue
                        qk = path_to_key.get(str(Path(p[0]).resolve()), Path(p[0]).stem)
                        rk = path_to_key.get(str(Path(p[1]).resolve()), Path(p[1]).stem)
                        if qk == rk:
                            continue
                        af = int(p[3]) / int(p[4]) if int(p[4]) else 0.0
                        hf.write(json.dumps({"genome": qk, "ref": rk, "ani": p[2],
                                             "orthologous_matches": p[3],
                                             "total_fragments": p[4],
                                             "af": round(af, 4),
                                             "scope_rank": "family" if fam_key else "species",
                                             "scope_name": fam_key or fam}) + "\n")
                        n += 1
            appended += n
            L.log(f"[06b] {fam[:45]:45s} q={len(queries)} r={len(refs)} hits={n} "
                  f"{round(time.time()-t,1)}s")

    # reconsolidación
    best: dict[tuple[str, str], dict] = {}
    for line in (FULL / "ani_hits_full.jsonl").read_text(encoding="utf-8").splitlines():
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        key = (rec["genome"], rec["ref"])
        if key not in best or float(rec["ani"]) > float(best[key]["ani"]):
            best[key] = rec
    rows = sorted(best.values(), key=lambda r: (r["genome"], -float(r["ani"])))
    with open(hits_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "ref", "ani", "orthologous_matches",
                                           "total_fragments", "af", "scope_rank",
                                           "scope_name"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    L.log(f"[06b] {appended} hits nuevos | total consolidado: {len(rows)} pares")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
