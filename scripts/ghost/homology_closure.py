#!/usr/bin/env python3
"""homology_closure.py — cierre por homología, independiente de Pfam.

El cierre por dominio (A1) y este comparten el punto de partida (las enzimas
productoras de Rhea) pero **no la evidencia**: aquí no se usan dominios Pfam,
sino alineamiento directo DIAMOND de las proteínas del genoma contra las
secuencias completas de esas enzimas. Si ambos métodos coinciden, la brecha
pasa a `validated` — dos caminos independientes, que es lo que el vocabulario
del proyecto exige para esa etiqueta.

Umbrales (declarados, no ajustados al set):
  pident >= 30%   identidad de la alineación
  qcov   >= 70%   fracción de la proteína consultada alineada
  evalue <= 1e-5

Uso
---
    python scripts/ghost/homology_closure.py --build-db
    python scripts/ghost/homology_closure.py --workers 8
    python scripts/ghost/homology_closure.py --merge
"""
from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "data" / "metabolism"
ANNOT = ROOT / "data" / "annotations_ab_initio"
WORKLIST = META / "domain_worklist_top{top}.tsv"
OUT = META / "homology_closure_top{top}.tsv"
FINAL = META / "closure_final_top{top}.tsv"
DOMAIN_VALID = META / "closure_validation_top{top}.tsv"
SPROT = ROOT / "data" / "db" / "uniprot_sprot.fasta"
DB = ROOT / ".cache" / "producers"
MIN_PIDENT, MIN_QCOV, MAX_EVALUE = 30.0, 0.70, 1e-5


def build_db(top: int) -> None:
    """Extrae del Swiss-Prot solo las enzimas productoras y hace el índice."""
    acc: set[str] = set()
    for r in csv.DictReader(open(Path(str(WORKLIST).format(top=top))),
                            delimiter="\t"):
        acc.update(a for a in r["producer_accs"].split(";") if a)
    faa = ROOT / ".cache" / "producers.faa"
    keep, n = False, 0
    with open(faa, "w") as out:
        for line in open(SPROT):
            if line.startswith(">"):
                parts = line[1:].split("|")
                keep = len(parts) >= 2 and parts[1] in acc
                n += keep
            if keep:
                out.write(line)
    subprocess.run(["diamond", "makedb", "--in", str(faa), "-d", str(DB),
                    "--quiet"], check=True)
    print(f"DB: {n}/{len(acc)} secuencias -> {DB}.dmnd")


def _scan_one(job: tuple) -> tuple[str, list[dict]]:
    genome, specs = job
    faa = ANNOT / genome / "proteins.faa"
    if not faa.exists():
        return genome, [{"genome": genome, "gap": r["gap"], "acc": "",
                         "qseqid": "", "pident": "", "qcov": "", "evalue": "",
                         "status": "no_proteins"} for r in specs]
    with tempfile.TemporaryDirectory(dir=ROOT / ".cache") as td:
        out = Path(td) / "hits.tsv"
        subprocess.run(
            ["diamond", "blastp", "--query", str(faa), "--db", str(DB),
             "--out", str(out), "--outfmt", "6", "qseqid", "sseqid", "pident",
             "length", "qlen", "evalue", "--evalue", str(MAX_EVALUE),
             "--max-target-seqs", "5", "--threads", "1", "--quiet"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        best: dict[str, dict] = {}
        for line in open(out):
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            acc = p[1].split("|")[1] if "|" in p[1] else p[1]
            pident, alen, qlen = float(p[2]), int(p[3]), int(p[4])
            qcov = alen / qlen if qlen else 0.0
            if pident < MIN_PIDENT or qcov < MIN_QCOV:
                continue
            prev = best.get(acc)
            if prev is None or float(p[5]) < float(prev["evalue"]):
                best[acc] = {"pident": pident, "qcov": qcov, "evalue": p[5],
                             "qseqid": p[0]}
    rows = []
    for r in specs:
        accs = [a for a in r["producer_accs"].split(";") if a]
        hits = [(a, best[a]) for a in accs if a in best]
        if hits:
            a, h = min(hits, key=lambda x: float(x[1]["evalue"]))
            rows.append({"genome": genome, "gap": r["gap"], "acc": a,
                         "qseqid": h["qseqid"],
                         "pident": f"{h['pident']:.1f}",
                         "qcov": f"{h['qcov']:.2f}", "evalue": h["evalue"],
                         "status": "homology_hit"})
        else:
            rows.append({"genome": genome, "gap": r["gap"], "acc": "",
                         "qseqid": "", "pident": "", "qcov": "", "evalue": "",
                         "status": "no_hit"})
    return genome, rows


def merge(top: int) -> None:
    """Cruza dominio x homología y emite el veredicto final en 2x2."""
    dom = {(r["genome"], r["gap"]): r for r in csv.DictReader(
        open(Path(str(DOMAIN_VALID).format(top=top))), delimiter="\t")}
    hom = {(r["genome"], r["gap"]): r for r in csv.DictReader(
        open(Path(str(OUT).format(top=top))), delimiter="\t")}
    out = []
    for key in sorted(set(dom) | set(hom)):
        d, h = dom.get(key), hom.get(key)
        # `d` solo existe si el método de dominio emitió candidato: eso ES
        # evidencia, sea `new`, `reproduced_*` o `disagreement_domain`.
        d_ok = d is not None
        h_ok = bool(h and h["status"] == "homology_hit")
        same_protein = bool(d_ok and h_ok
                            and d.get("new_query") == h.get("qseqid"))
        if d_ok and h_ok and d["verdict"] == "disagreement_domain":
            verdict, status = "both_but_conflict", "anomalous"
        elif d_ok and h_ok and same_protein:
            verdict, status = "validated_protein", "validated"
        elif d_ok and h_ok:
            verdict, status = "validated_genome", "validated"
        elif d_ok:
            verdict, status = "domain_only", "consensus"
        elif h_ok:
            verdict, status = "homology_only", "candidate"
        else:
            verdict, status = "none", "unresolved"
        out.append({
            "genome": key[0], "gap": key[1],
            "domain": (d or {}).get("new_domain", ""),
            "domain_verdict": (d or {}).get("verdict", ""),
            "homology_acc": (h or {}).get("acc", ""),
            "homology_pident": (h or {}).get("pident", ""),
            "homology_qcov": (h or {}).get("qcov", ""),
            "homology_protein": (h or {}).get("qseqid", ""),
            "same_protein": "yes" if verdict == "validated_protein" else "no",
            "verdict": verdict, "status": status,
        })
    path = Path(str(FINAL).format(top=top))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"\ncruce dominio x homología -> {path.relative_to(ROOT)}")
    for k, v in Counter(o["verdict"] for o in out).most_common():
        print(f"  {k:16s} {v}")
    print(f"\n  validated (dos métodos): "
          f"{sum(1 for o in out if o['status'] == 'validated')}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--build-db", action="store_true")
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()

    if a.build_db:
        build_db(a.top)
        return 0
    if a.merge:
        merge(a.top)
        return 0

    plan: dict[str, list[dict]] = defaultdict(list)
    for r in csv.DictReader(open(Path(str(WORKLIST).format(top=a.top))),
                            delimiter="\t"):
        if r["producer_accs"]:
            plan[r["genome"]].append(r)
    todo = sorted(plan)
    if a.limit:
        todo = todo[:a.limit]
    print(f"homología: {len(todo)} genomas | DB {DB}.dmnd | "
          f"{a.workers} workers", flush=True)

    fields = ["genome", "gap", "acc", "qseqid", "pident", "qcov", "evalue",
              "status"]
    fh = open(Path(str(OUT).format(top=a.top)), "w", newline="")
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    n = 0
    with mp.Pool(a.workers) as pool:
        for genome, rows in pool.imap(_scan_one,
                                      [(g, plan[g]) for g in todo], chunksize=1):
            w.writerows(rows)
            fh.flush()
            n += 1
            if n % 20 == 0 or n == len(todo):
                print(f"  [{n}/{len(todo)}] {genome}", flush=True)
    fh.close()

    rows = list(csv.DictReader(open(Path(str(OUT).format(top=a.top))),
                               delimiter="\t"))
    print(f"\nsalida: {len(rows)} filas")
    for k, v in Counter(r["status"] for r in rows).most_common():
        print(f"  {k:14s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
