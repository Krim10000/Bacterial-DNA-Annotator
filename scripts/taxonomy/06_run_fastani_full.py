#!/usr/bin/env python3
"""06 — ANI de los 400 genomas contra el panel, agrupado por filo.

Producto cartesiano 400 x 695 (278.000 pares) es inviable. Estrategia elegida
tras medir el coste real en esta máquina:

* Cada genoma se compara contra **todo el panel de su filo** (no contra el filo
  entero de la base ni contra todo el panel). El filo es la unidad porque:
  - cada genoma del panel pertenece a un solo filo -> se esboza una vez (695 en
    total), no una vez por grupo;
  - 55 grupos en vez de ~140, así el coste de arranque no domina;
  - el vecino más cercano se busca en el filo completo, sin recortes arbitrarios.
* Los 7 genomas sin filo en NCBI caen a **clase**; si tampoco tienen clase,
  quedan sin ANI (se declaran).
* Se conserva el alcance usado por genoma (`scope_rank`/`scope_name`) y el número
  de referencias del grupo.

Ejecución **resumible**: cada grupo anexa sus hits a `full/ani_hits_full.jsonl` y su
estado a `full/ani_groups.jsonl`. Si se corta, relanzar el script retoma donde
quedó. Al final consolida `full/ani_hits_full.tsv` y `full/ani_groups.json`.

No descarga nada: reutiliza binario, FASTA y GTDB del piloto.
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
PROGRESS = FULL / "ani_progress.log"
# Tope de referencias por grupo. El panel base es uno-por-especie, así que en la
# mayoría de filos no hay un pariente cercano real; el tope limita el coste y queda
# registrado por genoma (`n_refs`). Se priorizan conspecíficas, misma familia y
# material tipo.
MAX_REFS = int(os.environ.get("TAXONOMY_MAX_REFS", "5"))


def evenly(rows: list[str], limit: int) -> list[str]:
    rows = sorted(rows)
    if limit <= 0:
        return []
    if len(rows) <= limit:
        return rows
    if limit == 1:
        return rows[:1]
    step = (len(rows) - 1) / (limit - 1)
    return [rows[int(round(i * step))] for i in range(limit)]


def pick_refs(refs: list[str], queries: list[str], panel: dict[str, dict],
              cap: int = MAX_REFS) -> list[str]:
    """Conspecíficas > misma familia > material tipo > resto, hasta `cap`."""
    q_species = {panel.get(q, {}).get("species_taxid") or "" for q in queries}
    q_families = {panel.get(q, {}).get("rk_family") or "" for q in queries}
    cons = [r for r in refs
            if (panel.get(r, {}).get("species_taxid") or "") in q_species
            and (panel.get(r, {}).get("species_taxid") or "")]
    fam = [r for r in refs if r not in cons
           and (panel.get(r, {}).get("rk_family") or "") in q_families]
    types = [r for r in refs if r not in cons and r not in fam
             and (panel.get(r, {}).get("relation_to_type_material") or "na") != "na"]
    rest = [r for r in refs if r not in cons and r not in fam and r not in types]

    picked: list[str] = []
    for tier in (cons, fam, types, rest):
        room = cap - len(picked)
        if room <= 0:
            break
        picked += evenly(tier, room)
    return sorted(set(picked))


def main() -> int:
    FULL.mkdir(parents=True, exist_ok=True)
    meta = list(csv.DictReader(open(FULL / "genome_meta.tsv", encoding="utf-8"),
                               delimiter="\t"))
    panel = {r["genome"]: r for r in
             csv.DictReader(open(FULL / "panel_full.tsv", encoding="utf-8"), delimiter="\t")}
    paths = {g: r["fasta_path"] for g, r in panel.items() if r.get("fasta_path")}
    path_to_key = {str(Path(v).resolve()): k for k, v in paths.items()}

    # alcance por genoma
    scope_of: dict[str, tuple[str, str]] = {}
    for r in meta:
        if r.get("phylum"):
            scope_of[r["genome"]] = ("phylum", r["phylum"])
        elif r.get("class"):
            scope_of[r["genome"]] = ("class", r["class"])
    refs_by_scope: dict[tuple[str, str], list[str]] = defaultdict(list)
    for g, pr in panel.items():
        ph, cl = pr.get("rk_phylum") or "", pr.get("rk_class") or ""
        if ph:
            refs_by_scope[("phylum", ph)].append(g)
        elif cl:
            refs_by_scope[("class", cl)].append(g)

    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for r in meta:
        sc = scope_of.get(r["genome"])
        if sc:
            groups[sc].append(r["genome"])

    L.log(f"[ani-full] queries con alcance: {sum(len(v) for v in groups.values())}"
          f"/{len(meta)} | grupos: {len(groups)} | panel: {len(panel)}")

    groups_jsonl = FULL / "ani_groups.jsonl"
    hits_jsonl = FULL / "ani_hits_full.jsonl"
    done: set[tuple[str, str]] = set()
    if groups_jsonl.exists():
        for line in groups_jsonl.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("rc") == 0:
                done.add((rec["rank"], rec["name"]))
    # un grupo con hits ya escritos tampoco se repite (protege cortes a mitad)
    if hits_jsonl.exists():
        for line in hits_jsonl.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("scope_rank") and rec.get("scope_name"):
                done.add((rec["scope_rank"], rec["scope_name"]))
    L.log(f"[ani-full] grupos ya completados: {len(done)}")

    t0 = time.time()
    with open(groups_jsonl, "a", encoding="utf-8") as glog:
        for i, ((rank, name), queries) in enumerate(sorted(groups.items()), 1):
            if (rank, name) in done:
                continue
            # los propios queries son referencias validas entre si: el par
            # genoma-genoma se descarta despues (qk == rk), no antes
            all_refs = [r for r in refs_by_scope.get((rank, name), []) if r in paths]
            refs = pick_refs(all_refs, queries, panel)
            if not refs:
                refs = all_refs
            if not refs:
                glog.write(json.dumps({"rank": rank, "name": name,
                                       "n_queries": len(queries), "n_refs": 0,
                                       "rc": 0, "hits": 0, "seconds": 0.0}) + "\n")
                glog.flush()
                continue
            qfile = L.Q_DIR / "full_q.txt"
            rfile = L.Q_DIR / "full_r.txt"
            qfile.write_text("\n".join(paths[g] for g in queries) + "\n", encoding="utf-8")
            rfile.write_text("\n".join(paths[r] for r in refs) + "\n", encoding="utf-8")
            outfile = L.ANI_DIR / "full_out.tsv"
            outfile.unlink(missing_ok=True)

            cmd = [L.fastani_bin(), "--ql", str(qfile), "--rl", str(rfile),
                   "-o", str(outfile), "--minFraction", L.FASTANI_MIN_FRACTION,
                   "-t", os.environ.get("TAXONOMY_THREADS", "8")]
            start = time.time()
            proc = subprocess.run(cmd, capture_output=True, text=True)
            elapsed = round(time.time() - start, 1)
            n_hits = 0
            if proc.returncode == 0 and outfile.exists():
                with open(outfile, encoding="utf-8") as fh, \
                        open(hits_jsonl, "a", encoding="utf-8") as hf:
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
                                             "scope_rank": rank,
                                             "scope_name": name}) + "\n")
                        n_hits += 1
            rec = {"rank": rank, "name": name, "n_queries": len(queries),
                   "n_refs": len(refs), "rc": proc.returncode, "hits": n_hits,
                   "seconds": elapsed}
            if proc.returncode != 0:
                rec["error"] = (proc.stderr or "")[-400:]
                L.log(f"[ani-full] FALLO {rank}:{name} rc={proc.returncode}")
            glog.write(json.dumps(rec) + "\n")
            glog.flush()
            msg = (f"[ani-full] {i}/{len(groups)} {rank}:{name} q={len(queries)} "
                   f"r={len(refs)} hits={n_hits} {elapsed}s")
            L.log(msg)
            with open(PROGRESS, "a", encoding="utf-8") as pf:
                pf.write(msg + "\n")

    # ------------- consolidación (idempotente) -------------
    best: dict[tuple[str, str], dict] = {}
    if hits_jsonl.exists():
        for line in hits_jsonl.read_text(encoding="utf-8").splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            key = (rec["genome"], rec["ref"])
            if key not in best or float(rec["ani"]) > float(best[key]["ani"]):
                best[key] = rec
    rows = list(best.values())
    rows.sort(key=lambda r: (r["genome"], -float(r["ani"])))
    with open(FULL / "ani_hits_full.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "ref", "ani", "orthologous_matches",
                                           "total_fragments", "af", "scope_rank",
                                           "scope_name"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    glog_rows = [json.loads(l) for l in groups_jsonl.read_text(encoding="utf-8").splitlines()
                 if l.strip()]
    summary = {
        "n_groups": len(glog_rows),
        "n_groups_failed": sum(1 for g in glog_rows if g.get("rc") != 0),
        "n_pairs": sum(g["n_queries"] * g["n_refs"] for g in glog_rows),
        "n_hits": len(rows),
        "n_genomes_with_hits": len({r["genome"] for r in rows}),
        "elapsed_s": round(time.time() - t0, 1),
        "failures": [g for g in glog_rows if g.get("rc") != 0],
        "groups": glog_rows,
    }
    L.write_json(FULL / "ani_groups.json", summary)
    L.log(f"[ani-full] consolidado: {summary['n_hits']} hits / {summary['n_pairs']} "
          f"pares / {summary['n_groups']} grupos / {summary['n_groups_failed']} fallos")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
