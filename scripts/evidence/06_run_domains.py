#!/usr/bin/env python3
"""06_run_domains.py — hmmscan del panel Pfam contra las proteomas del lote.

Un `hmmscan --cut_ga` por genoma y fuente (`ncbi` y `abinitio`). Para no
escanear 65 modelos en cada genoma, se construye un **panel mínimo por
genoma** con solo las familias Pfam que ese genoma necesita (las esperadas por
sus brechas, de `work/gap_domains.tsv`), extraídas con `hmmfetch` del panel del
paso 04. El resultado es el mismo que el panel completo para esas familias.

Los aciertos por dominio quedan en `work/domtbl/<genome>.tsv.gz` y el resumen
en `work/domain_hits.tsv`:

    genome  source  protein_id  pfam_acc  pfam_name  i_evalue  score  qcov_dom

`--cut_ga` es el corte por gathering threshold de Pfam (no se ajusta a mano).
El emparejamiento con los dominios esperados de cada brecha se hace en
07_build_evidence.py.

Uso:
    python scripts/evidence/06_run_domains.py [--workers 7]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (ABINITIO, ANNOT, WORK, read_tsv, write_json, write_tsv)

PANEL = WORK / "panel_evidence.hmm"
PANEL_DIR = WORK / "panels"
DOMTBL = WORK / "domtbl"
DOMTBL.mkdir(parents=True, exist_ok=True)
PANEL_DIR.mkdir(parents=True, exist_ok=True)


def _decompress(src: Path, dst: Path) -> None:
    with gzip.open(src, "rt") as fh, open(dst, "w") as out:
        out.write(fh.read())


def build_panel(accs: set[str], name_by_acc: dict[str, str]) -> Path | None:
    """Panel mínimo (hmmfetch + hmmpress) cacheado por conjunto de familias."""
    names = sorted({name_by_acc[a] for a in accs if a in name_by_acc})
    if not names:
        return None
    key = hashlib.sha1("\n".join(names).encode()).hexdigest()[:12]
    dest = PANEL_DIR / key / "panel.hmm"
    if dest.exists() and Path(str(dest) + ".h3i").exists():
        return dest
    d = dest.parent
    d.mkdir(parents=True, exist_ok=True)
    keys = d / "keys.txt"
    keys.write_text("\n".join(names) + "\n")
    tmp = d / f"panel.{os.getpid()}.hmm"
    with open(tmp, "w") as out:
        subprocess.run(["hmmfetch", "-f", str(PANEL), str(keys)],
                       check=True, stdout=out)
    subprocess.run(["hmmpress", "-f", str(tmp)], check=True,
                   stdout=subprocess.DEVNULL)
    os.replace(tmp, dest)
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        Path(str(tmp) + ext).unlink(missing_ok=True)
    subprocess.run(["hmmpress", "-f", str(dest)], check=True,
                   stdout=subprocess.DEVNULL)
    return dest


def _scan_one(job: tuple[str, str | None]) -> tuple[str, list[dict]]:
    genome, panel = job
    raw_path = DOMTBL / f"{genome}.tsv.gz"
    if raw_path.exists() and raw_path.stat().st_size > 0:
        return genome, _parse(genome, raw_path)
    sources = {"ncbi": ANNOT / f"{genome}.faa.gz",
               "abinitio": ABINITIO / genome / "proteins.faa"}
    if panel is None:
        # sin dominios esperados: no hay nada que escanear, pero se marca el
        # genoma como escaneado para no confundir "sin base" con "ausente".
        with gzip.open(raw_path, "wt") as gz:
            gz.write("# sin dominios Pfam esperados para este genoma\n")
        return genome, []
    tmp_out = DOMTBL / f"{genome}.tmp.domtbl"
    with tempfile.TemporaryDirectory(dir=WORK) as td, \
            gzip.open(raw_path, "wt") as gz:
        for name, src in sources.items():
            if not src.exists():
                continue
            q = src
            if str(q).endswith(".gz"):
                q = Path(td) / f"{name}.faa"
                _decompress(src, q)
            subprocess.run(
                ["hmmscan", "--cut_ga", "--cpu", "1", "--noali",
                 "--domtblout", str(tmp_out), "-o", os.devnull,
                 str(panel), str(q)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            with open(tmp_out) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    p = line.split()
                    if len(p) < 22:
                        continue
                    gz.write("\t".join([name] + p[:22]) + "\n")
            tmp_out.unlink(missing_ok=True)
    return genome, _parse(genome, raw_path)


def _parse(genome: str, raw_path: Path) -> list[dict]:
    best: dict[tuple, dict] = {}
    with gzip.open(raw_path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 23:
                continue
            # columnas: 0 source | 1 target name | 2 target acc | 4 query name
            #           6 qlen | 13 i-Evalue | 14 score | 18-19 ali from-to
            source, qname = p[0], p[4]
            qlen = int(float(p[6]))
            pfam_name, pfam_acc = p[1], p[2]
            i_evalue = float(p[13])
            score = float(p[14])
            ali_from, ali_to = int(p[18]), int(p[19])
            qcov = (ali_to - ali_from + 1) / qlen if qlen else 0.0
            acc = pfam_acc.split(".")[0]
            rec = {"genome": genome, "source": source, "protein_id": qname,
                   "pfam_acc": acc, "pfam_name": pfam_name,
                   "i_evalue": i_evalue, "score": score, "qcov_dom": qcov}
            key = (source, qname, acc)
            prev = best.get(key)
            if prev is None or i_evalue < prev["i_evalue"]:
                best[key] = rec
    return list(best.values())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=min(7, os.cpu_count() or 1))
    a = ap.parse_args()
    if not PANEL.exists():
        raise SystemExit(f"falta el panel: {PANEL} (corre 04_build_domain_panel.py)")

    sel = read_tsv(WORK / "selected_gaps.tsv")
    gd = read_tsv(WORK / "gap_domains.tsv")
    name_by_acc = {r["requested_acc"]: r["name"]
                   for r in read_tsv(WORK / "panel_meta.tsv")}

    want: dict[str, set[str]] = defaultdict(set)
    for r in gd:
        want[r["genome"]] |= {d for d in r["domains"].split(";") if d}
    genomes = sorted({r["genome"] for r in sel})
    pending = [g for g in genomes
               if not (DOMTBL / f"{g}.tsv.gz").exists()
               or (DOMTBL / f"{g}.tsv.gz").stat().st_size == 0]
    print(f"genomas: {len(genomes)}; ya escaneados: {len(genomes) - len(pending)}; "
          f"pendientes: {len(pending)}")

    # paneles mínimos: uno por conjunto distinto de familias, en serie
    panels: dict[str, Path | None] = {}
    job_panels: dict[str, Path | None] = {}
    for g in pending:
        accs = want.get(g, set())
        key = tuple(sorted(accs))
        if key not in panels:
            panels[key] = build_panel(set(accs), name_by_acc)
        job_panels[g] = panels[key]
    print(f"paneles distintos: {sum(1 for v in panels.values() if v)}")

    rows: list[dict] = []
    done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_scan_one, (g, job_panels[g])): g for g in pending}
        for fut in as_completed(futs):
            _, part = fut.result()
            rows.extend(part)
            done += 1
            if done % 10 == 0 or done == len(futs):
                print(f"  {done}/{len(futs)} pendientes", file=sys.stderr)
    # incorpora lo ya escaneado en corridas previas
    for g in genomes:
        if g not in job_panels:
            rows.extend(_parse(g, DOMTBL / f"{g}.tsv.gz"))

    rows.sort(key=lambda r: (r["genome"], r["source"], r["protein_id"],
                             r["pfam_acc"]))
    write_tsv(WORK / "domain_hits.tsv", rows,
              ["genome", "source", "protein_id", "pfam_acc", "pfam_name",
               "i_evalue", "score", "qcov_dom"])
    stats = {
        "genomes": len(genomes),
        "rows": len(rows),
        "sources": dict(Counter(r["source"] for r in rows)),
        "proteins_with_domain": len({(r["genome"], r["source"], r["protein_id"])
                                     for r in rows}),
        "pfam_detected": len({r["pfam_acc"] for r in rows}),
        "panel": str(PANEL),
        "per_genome_panels": sum(1 for v in panels.values() if v),
        "cut": "cut_ga",
    }
    write_json(WORK / "domains_scan.json", stats)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
