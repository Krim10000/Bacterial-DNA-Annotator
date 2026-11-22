#!/usr/bin/env python3
"""scan_domain_closure.py — cierra brechas por dominio sobre el lote preparado.

Para cada genoma del lote `domain_worklist_topN.tsv`, escanea sus proteínas una
sola vez contra el panel de dominios (144 modelos Pfam extraídos con hmmfetch)
y emite un candidato cuando el dominio del productor aparece en una proteína
del genoma. Un `hmmscan` por genoma, no uno por brecha.

Criterio de cierre (el mismo del A1, `data/metabolism/DOMAIN_CLOSURE.md`):
  - el dominio debe ser **específico**: los genéricos (`FAD_binding_4`,
    `NAD_binding_1`, Rossmann, aminotransferasas...) no cierran por sí solos;
  - `--cut-ga` para el corte por gathering threshold de Pfam;
  - la dirección Rhea (`LR`/`RL`) sube la confianza a `high`; `UN`/`BI` la deja
    en `medium` — no se inventa direccionalidad;
  - los dominios permitidos los fija `build_domain_worklist.py` con la política
    `core` (en todos los productores) > `majority` (>=60%) > `any`. Elegir por
    e-value entre TODOS los dominios producía falsos positivos (p. ej. APH en
    vez de GHMP_kinases_N para homoserina quinasa).

Estados de cada fila:
  candidate        dominio específico presente, dirección direccional
  candidate_weak   dominio específico presente, reacción reversible/incierta
  generic_only     solo dominios genéricos
  no_hit           ninguno de los dominios del productor aparece

Uso
---
    python scripts/ghost/scan_domain_closure.py --top 1000 --threads 4
    python scripts/ghost/scan_domain_closure.py --top 1000 --limit 10 --force
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import multiprocessing as mp
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANNOT = ROOT / "data" / "annotations_ab_initio"
PANEL_RAW = ROOT / "data" / "db" / "pfam" / "panel_domain_closure.hmm"
PANEL_DIR = ROOT / ".cache" / "panels"
WORKLIST = ROOT / "data" / "metabolism" / "domain_worklist_top{top}.tsv"
OUT = ROOT / "data" / "metabolism" / "domain_closure_candidates_top{top}.tsv"

GENERIC = {
    "FAD_binding_4", "FAD_binding_2", "NAD_binding_1", "NAD_binding_2",
    "Rossmann-like", "Aldedh", "Aminotran_1_2", "Aminotran_3", "Aminotran_4",
    "Aminotran_5", "Metallophos", "Radical_SAM", "GAF_2", "AAA_26",
}


def load_plan(path: Path) -> dict[str, list[dict]]:
    """genoma -> filas del lote (solo las que tienen dominios)."""
    plan: dict[str, list[dict]] = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        doms = [d for d in r["pfam_domains"].split(";") if d]
        if not doms:
            continue
        r["_domains"] = set(doms)
        plan[r["genome"]].append(r)
    return plan


_LOCK = None


def _init_lock() -> None:
    """Lock por proceso para que dos workers no construyan el mismo panel."""
    global _LOCK
    _LOCK = mp.Lock()


def panel_for(domains: set[str]) -> Path | None:
    """Panel mínimo con solo los dominios que este genoma necesita.

    Escanear los 144 modelos en cada genoma tardaba ~1,3 min; la mayoría de los
    genomas necesita 3-5. El panel se cachea por conjunto de dominios.
    """
    key = hashlib.sha1("|".join(sorted(domains)).encode()).hexdigest()[:12]
    d = PANEL_DIR / key
    hmm = d / "panel.hmm"
    if (d / "panel.hmm.h3m").exists():
        return hmm
    if _LOCK is not None:
        _LOCK.acquire()
        try:
            return _panel_build(d, hmm, domains)
        finally:
            _LOCK.release()
    return _panel_build(d, hmm, domains)


def _panel_build(d: Path, hmm: Path, domains: set[str]) -> Path | None:
    if (d / "panel.hmm.h3m").exists():
        return hmm
    d.mkdir(parents=True, exist_ok=True)
    names = d / "names.txt"
    names.write_text("\n".join(sorted(domains)) + "\n")
    with open(hmm, "w") as fh:
        subprocess.run(["hmmfetch", "-f", str(PANEL_RAW), str(names)],
                       stdout=fh, check=True)
    if not hmm.stat().st_size:
        return None
    subprocess.run(["hmmpress", "-f", str(hmm)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return hmm


def scan_genome(faa: Path, panel: Path, threads: int) -> dict[str, list[dict]]:
    """dominio -> hits (query, evalue, score) del genoma."""
    with tempfile.TemporaryDirectory(dir=ROOT / ".cache") as td:
        out = Path(td) / "domtbl"
        subprocess.run(
            ["hmmscan", "--cut_ga", "--domtblout", str(out), "--noali",
             "--cpu", str(threads), str(panel), str(faa)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        hits: dict[str, list[dict]] = defaultdict(list)
        for line in open(out):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            hits[f[0]].append({"query": f[3], "evalue": float(f[6]),
                               "score": float(f[7])})
    return hits


def _scan_one(job: tuple) -> tuple[str, list[dict]]:
    """Escanea un genoma. Se ejecuta en un proceso worker."""
    genome, specs, threads = job
    faa = ANNOT / genome / "proteins.faa"
    if not faa.exists():
        return genome, [{"genome": genome, "gap": r["gap"], "chebi": r["chebi"],
                         "role": r["role"], "support_400": r["support_400"],
                         "domain_policy": r.get("domain_policy", ""),
                         "domain": "", "specific": "", "query": "", "evalue": "",
                         "score": "", "direction": r["direction"],
                         "rhea_master": r["rhea_master"], "status": "no_proteins"}
                        for r in specs]
    panel = panel_for({d for r in specs for d in r["_domains"]})
    hits = scan_genome(faa, panel, threads) if panel else {}
    rows = []
    for r in specs:
        found = [(d, h) for d in sorted(r["_domains"]) for h in hits.get(d, [])]
        specific = [(d, h) for d, h in found if d not in GENERIC]
        if specific:
            d, h = min(specific, key=lambda x: x[1]["evalue"])
            status = "candidate" if r["direction"] in ("LR", "RL") else "candidate_weak"
            row = {"domain": d, "specific": "yes", "query": h["query"],
                   "evalue": f"{h['evalue']:.2e}", "score": f"{h['score']:.1f}",
                   "status": status}
        elif found:
            d, h = min(found, key=lambda x: x[1]["evalue"])
            row = {"domain": d, "specific": "no", "query": h["query"],
                   "evalue": f"{h['evalue']:.2e}", "score": f"{h['score']:.1f}",
                   "status": "generic_only"}
        else:
            row = {"domain": "", "specific": "", "query": "", "evalue": "",
                   "score": "", "status": "no_hit"}
        row.update({k: r[k] for k in ("genome", "gap", "chebi", "role",
                                      "support_400", "direction", "rhea_master")})
        row["domain_policy"] = r.get("domain_policy", "")
        rows.append(row)
    return genome, rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=1,
                    help="hilos de hmmscan POR worker (default 1)")
    ap.add_argument("--workers", type=int, default=0,
                    help="procesos en paralelo (0 = núcleos/cpu, máx 8)")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true",
                    help="re-escanea genomas ya presentes en la salida")
    a = ap.parse_args()

    workers = a.workers or min(mp.cpu_count(), 8)
    plan = load_plan(Path(str(WORKLIST).format(top=a.top)))
    out_path = Path(str(OUT).format(top=a.top))

    done: set[str] = set()
    if out_path.exists() and not a.force:
        done = {r["genome"] for r in csv.DictReader(open(out_path), delimiter="	")}
    todo = [g for g in sorted(plan) if g not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"lote {a.top}: {len(plan)} genomas | ya hechos {len(done)} | "
          f"a escanear {len(todo)} | {workers} workers x {a.threads} hilo(s)",
          flush=True)

    fields = ["genome", "gap", "chebi", "role", "support_400", "domain_policy", "domain",
              "specific", "query", "evalue", "score", "direction",
              "rhea_master", "status"]
    write_header = not out_path.exists() or a.force or not done
    fh = open(out_path, "w" if write_header else "a", newline="")
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="	")
    if write_header:
        w.writeheader()

    jobs = [(g, plan[g], a.threads) for g in todo]
    n = 0
    with mp.Pool(workers, initializer=_init_lock) as pool:
        # imap conserva el orden de los genomas: salida determinista
        for genome, rows in pool.imap(_scan_one, jobs, chunksize=1):
            w.writerows(rows)
            fh.flush()
            n += 1
            if n % 10 == 0 or n == len(jobs):
                print(f"  [{n}/{len(jobs)}] {genome}", flush=True)
    fh.close()

    rows = list(csv.DictReader(open(out_path), delimiter="\t"))
    counts = defaultdict(int)
    for r in rows:
        counts[r["status"]] += 1
    print(f"\nsalida: {out_path.relative_to(ROOT)} ({len(rows)} filas)")
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {k:14s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
