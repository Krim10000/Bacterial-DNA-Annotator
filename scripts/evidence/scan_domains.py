#!/usr/bin/env python3
"""scan_domains.py -- evidencia de dominio para los 77 accionables de `gap_review.tsv`.

Lee los 17 `annotation_gap` + 60 `domain_present` y, para cada genoma implicado,
escanea su proteoma (pyrodigal, `data/annotations_ab_initio/<genome>/proteins.faa`)
contra los modelos Pfam que el revisor de brechas asoció a cada gap
(`data/metabolism/NO_HOMOLOG_DOMAINS.tsv`, columna `domains`), más los modelos
del productor que faltaban en el panel original cuando este solo traía dominios
de la enzima consumidora:

  - `S-adenosyl 3-(methylsulfanyl)propylamine` -> `AdoMet_dc` (SpeD, productor;
    el panel original solo traía `Spermine_synth`, que es la enzima consumidora).

Dos pasadas por genoma:
  - `ga`      : `hmmscan --cut_ga` (corte por gathering threshold de Pfam).
  - `relaxed` : `-E 1e-3 --domE 1e-3` (señal divergente, para validación).

Los modelos se toman del panel ya construido por el pipeline
(`data/db/pfam/panel_domain_closure.hmm`, 144 modelos) y de
`_work/AdoMet_dc.hmm`; **no** se vuelve a leer `Pfam-A.hmm` (3,6 GB) por genoma.
Cada panel por genoma se ensambla concatenando los modelos y se cachea en
`data/metabolism/reannotation/_work/panels/`.

Salida: `data/metabolism/reannotation/domain_hits.tsv` (un hit de dominio por
fila: ORF candidato, e-value, score y cobertura del modelo). Solo se escribe
dentro de `data/metabolism/reannotation/`.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT / "data" / "metabolism" / "review" / "gap_review.tsv"
NOHOM = ROOT / "data" / "metabolism" / "NO_HOMOLOG_DOMAINS.tsv"
OUTDIR = ROOT / "data" / "metabolism" / "reannotation"
WORK = OUTDIR / "_work"
MODELS = WORK / "models"
PANELS = WORK / "panels"
PANEL_SRC = ROOT / "data" / "db" / "pfam" / "panel_domain_closure.hmm"
ADOMET_DC = WORK / "AdoMet_dc.hmm"
OUT_HITS = OUTDIR / "domain_hits.tsv"

ACTIONABLE = ("annotation_gap", "domain_present")
EXTRA_DOMAINS = {
    "S-adenosyl 3-(methylsulfanyl)propylamine": ["AdoMet_dc"],
}

HIT_COLS = [
    "genome", "orf", "orf_len", "domain", "mode", "full_evalue", "full_score",
    "dom_ievalue", "dom_score", "hmm_len", "ali_from", "ali_to", "ali_cov",
]


def load_cases() -> list[dict]:
    return [r for r in csv.DictReader(open(REVIEW), delimiter="\t")
            if r["classification"] in ACTIONABLE]


def load_expected() -> dict[tuple[str, str], list[str]]:
    exp = {}
    for r in csv.DictReader(open(NOHOM), delimiter="\t"):
        exp[(r["genome"], r["gap"])] = [d for d in r["domains"].split(";") if d]
    return exp


def genome_plan(cases: list[dict], expected: dict) -> dict[str, set[str]]:
    plan: dict[str, set[str]] = {}
    for c in cases:
        key = (c["genome"], c["gap"])
        doms = set(expected.get(key, []))
        doms |= set(EXTRA_DOMAINS.get(c["gap"], []))
        if not doms:
            raise SystemExit(f"sin dominios esperados para {key}")
        plan.setdefault(c["genome"], set()).update(doms)
    return plan


def _index_hmm(path: Path) -> int:
    """Escribe un archivo HMM por modelo a partir de un HMM multi-modelo."""
    n = 0
    name = None
    buf: list[str] = []
    for line in open(path):
        if line.startswith("HMMER3/"):
            buf = [line]
            name = None
            continue
        buf.append(line)
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("//") and name:
            (MODELS / f"{name}.hmm").write_text("".join(buf))
            n += 1
    return n


def split_models() -> None:
    """Indexa el panel de 144 modelos y los modelos extra como HMM por modelo."""
    MODELS.mkdir(parents=True, exist_ok=True)
    if len(list(MODELS.glob("*.hmm"))) < 144:
        _index_hmm(PANEL_SRC)
    if not ADOMET_DC.exists():
        raise SystemExit(
            f"falta {ADOMET_DC}; generarlo una vez con:\n"
            f"  hmmfetch -f .cache/pfam/Pfam-A.hmm <(echo AdoMet_dc) > {ADOMET_DC}")
    _index_hmm(ADOMET_DC)
    extra = WORK / "extra_models.hmm"
    if extra.exists():
        _index_hmm(extra)


def panel_for(domains: set[str]) -> Path:
    """Panel Pfam mínimo (cacheado por conjunto de modelos)."""
    key = hashlib.sha1("|".join(sorted(domains)).encode()).hexdigest()[:12]
    d = PANELS / key
    hmm = d / "panel.hmm"
    if (d / "panel.hmm.h3m").exists():
        return hmm
    missing = [x for x in domains if not (MODELS / f"{x}.hmm").exists()]
    if missing:
        raise SystemExit(f"modelos ausentes: {sorted(missing)}")
    d.mkdir(parents=True, exist_ok=True)
    with open(hmm, "w") as out:
        for name in sorted(domains):
            out.write((MODELS / f"{name}.hmm").read_text())
    subprocess.run(["hmmpress", "-f", str(hmm)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return hmm


def parse_domtbl(path: Path, genome: str, mode: str) -> list[dict]:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 22:
                continue
            tlen = int(f[2])
            ali_from, ali_to = int(f[17]), int(f[18])
            rows.append({
                "genome": genome,
                "orf": f[3],
                "orf_len": int(f[5]),
                "domain": f[0],
                "mode": mode,
                "full_evalue": f[6],
                "full_score": f[7],
                "dom_ievalue": f[12],
                "dom_score": f[13],
                "hmm_len": tlen,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "ali_cov": f"{round((ali_to - ali_from + 1) / tlen, 3):.3f}",
            })
    return rows


def scan_one(args: tuple[str, list[str], int]) -> list[dict]:
    genome, domains, threads = args
    faa = ROOT / "data" / "annotations_ab_initio" / genome / "proteins.faa"
    if not faa.exists():
        raise SystemExit(f"falta proteoma: {faa}")
    panel = panel_for(set(domains))
    rows: list[dict] = []
    with tempfile.TemporaryDirectory(dir=WORK) as td:
        for mode, extra in (("ga", ["--cut_ga"]),
                            ("relaxed", ["-E", "1e-3", "--domE", "1e-3"])):
            out = Path(td) / f"{mode}.domtbl"
            cmd = ["hmmscan", "--noali", "--domtblout", str(out),
                   "--cpu", str(threads), *extra, str(panel), str(faa)]
            for attempt in (1, 2):
                proc = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                                      stderr=subprocess.PIPE)
                if proc.returncode == 0:
                    break
                if attempt == 2:
                    raise RuntimeError(
                        f"hmmscan falló {genome} {mode}: "
                        f"{proc.stderr.decode()[-800:]}")
            rows.extend(parse_domtbl(out, genome, mode))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="solo N genomas")
    ap.add_argument("--genome", default="", help="un genoma concreto")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if OUT_HITS.exists() and not args.force:
        print(f"ya existe {OUT_HITS} (usa --force)", file=sys.stderr)
        return 0
    OUTDIR.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    PANELS.mkdir(parents=True, exist_ok=True)
    split_models()

    cases = load_cases()
    expected = load_expected()
    plan = genome_plan(cases, expected)
    genomes = sorted(plan)
    if args.genome:
        genomes = [g for g in genomes if g == args.genome]
    if args.limit:
        genomes = genomes[: args.limit]
    print(f"casos={len(cases)} genomas={len(genomes)} "
          f"modelos={sorted(set().union(*plan.values()))}", file=sys.stderr)

    tasks = [(g, sorted(plan[g]), args.threads) for g in genomes]
    rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, chunk in enumerate(ex.map(scan_one, tasks), 1):
            rows.extend(chunk)
            print(f"  [{i}/{len(tasks)}] {tasks[i-1][0]}: {len(chunk)} hits",
                  file=sys.stderr)

    rows.sort(key=lambda r: (r["genome"], r["mode"], r["domain"],
                             float(r["dom_ievalue"])))
    with open(OUT_HITS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HIT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"escrito {OUT_HITS} ({len(rows)} hits)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
