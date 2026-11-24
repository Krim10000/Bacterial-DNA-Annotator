#!/usr/bin/env python3
"""Lote del modelo genómico K-12 sobre los 400 genomas del proyecto.

Para cada genoma con `genes.tsv` y operones, ejecuta
`scripts/complexes/build_candidates.py --apply-only` con el modelo
`model_k12.genomic.json` y convierte los `clusters_<genoma>.tsv` en
`model_<genoma>.complexes.tsv`, con las columnas que consume
`webapp/complexes.py`.

Reglas del lote
---------------
- sin red: se pasan `--string /dev/null` y `--synteny /dev/null`;
- solo señales genómicas: el modelo es el respaldo sin STRING de K-12;
- no se usa STRING ni sintenia de K-12 en otros genomas;
- cada fila es una **sospecha** (`complex_candidate`, `yellow`); no cierra
  brechas ni demuestra interacción física.

Los operones se buscan primero en `data/operon_finder/<genoma>_operons.tsv` y,
si no existen, en `data/regulation/<genoma>/operons.tsv` (el layout real de
los 400). El genoma de entrenamiento (`ecoli_k12`) se excluye por defecto.

Uso
---
    python tools/run_complex_model_400.py --limit 3
    python tools/run_complex_model_400.py --threads 8 --force
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "scripts" / "complexes" / "build_candidates.py"
MODEL_NAME = "model_k12.genomic.json"
TRAINING_GENOME = "ecoli_k12"
BUILD_TIMEOUT = 900
COLS = ["complex_id", "name", "contig", "start", "end", "subunits",
        "expected", "missing", "status", "evidence", "confidence",
        "source", "ecs", "flag", "limitation"]
LIMITATION = ("Sospecha por modelo genómico K-12; no demuestra interacción "
              "física ni estequiometría. No cierra brechas.")


@dataclass(frozen=True)
class Job:
    genome: str
    genes: Path
    operons: Path
    operon_source: str


@dataclass
class JobResult:
    genome: str
    operons: str = ""
    clusters: int = 0
    altas: int = 0
    omitidos: int = 0
    fallos: int = 0
    status: str = "ok"
    error: str = ""


@dataclass
class GeneCoord:
    contig: str
    start: int
    end: int
    name: str = ""


def discover_genomes(data_root: Path) -> list[Job]:
    """Genomas con `genes.tsv` y operones, ordenados por nombre.

    Se prefiere `operon_finder/<genoma>_operons.tsv`; el lote histórico usa
    `regulation/<genoma>/operons.tsv`. El genoma de entrenamiento se excluye.
    """
    reg = Path(data_root) / "regulation"
    if not reg.is_dir():
        return []
    jobs: list[Job] = []
    for gdir in sorted(reg.iterdir(), key=lambda p: p.name):
        if not gdir.is_dir():
            continue
        genome = gdir.name
        if genome == TRAINING_GENOME:
            continue
        genes = gdir / "genes.tsv"
        if not genes.is_file():
            continue
        op = Path(data_root) / "operon_finder" / f"{genome}_operons.tsv"
        source = "operon_finder"
        if not op.is_file():
            op = gdir / "operons.tsv"
            source = "regulation"
        if not op.is_file():
            continue
        jobs.append(Job(genome, genes, op, source))
    return jobs


def load_gene_coords(path: Path) -> dict[str, GeneCoord]:
    """Coordenadas por locus_tag (o gene_id), igual que `read_genes_tsv`."""
    out: dict[str, GeneCoord] = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tag = (r.get("locus_tag") or r.get("gene_id") or "").strip()
            if not tag:
                continue
            try:
                start = int(float(r.get("start") or 0))
                end = int(float(r.get("end") or 0))
            except ValueError:
                continue
            contig = (r.get("contig") or "").strip()
            if not contig or start <= 0 or end < start:
                continue
            out[tag] = GeneCoord(
                contig=contig, start=start, end=end,
                name=(r.get("gene_name") or r.get("gene") or "").strip(),
            )
    return out


def _count_rows(path: Path) -> int:
    if not path.is_file():
        return 0
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in csv.DictReader(fh, delimiter="\t"))


def convert_clusters(clusters: list[dict], coords: dict[str, GeneCoord],
                     model_name: str = MODEL_NAME) -> tuple[list[dict], int]:
    """Convierte clusters del modelo en filas compatibles con `complexes.py`.

    Un cluster se omite si no tiene al menos dos genes con coordenadas en un
    mismo contig. Si hay genes en varios contigs, se conserva el contig con más
    genes localizados para que `start`/`end` sigan siendo un intervalo válido.
    """
    rows: list[dict] = []
    omitted = 0
    for r in clusters:
        tags = [g.strip() for g in str(r.get("genes") or "").split(",") if g.strip()]
        located = [g for g in tags if g in coords]
        by_contig: dict[str, list[str]] = {}
        for tag in located:
            by_contig.setdefault(coords[tag].contig, []).append(tag)
        if not by_contig:
            omitted += 1
            continue
        contig = max(
            by_contig,
            key=lambda c: (len(by_contig[c]), -located.index(by_contig[c][0])),
        )
        keep = by_contig[contig]
        if len(keep) < 2:
            omitted += 1
            continue

        raw_names = [n.strip() for n in
                     str(r.get("gene_names") or "").split(",")]
        name_by_tag = dict(zip(tags, raw_names))
        names = [coords[t].name or name_by_tag.get(t) or t for t in keep]
        starts = [coords[t].start for t in keep]
        ends = [coords[t].end for t in keep]

        mean_score = r.get("mean_score", "")
        min_score = r.get("min_score", "")
        ev = [f"modelo K-12 ({model_name})"]
        if r.get("evidence"):
            ev.append(str(r["evidence"]))
        if mean_score != "":
            ev.append(f"score medio={mean_score}; score mínimo={min_score}")
        rows.append({
            "complex_id": r.get("cluster_id", ""),
            "name": "Complejo candidato K-12: " + ", ".join(names),
            "contig": contig,
            "start": min(starts),
            "end": max(ends),
            "subunits": ";".join(f"{n}:{t}" for n, t in zip(names, keep)),
            "expected": "",
            "missing": "",
            "status": "complex_candidate",
            "evidence": "; ".join(ev),
            "confidence": r.get("confidence", ""),
            "source": "modelo K-12 (build_candidates.py)",
            "ecs": "",
            "flag": "yellow",
            "limitation": LIMITATION,
        })
    return rows, omitted


def write_complexes(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=COLS,
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def build_command(job: Job, tmp_dir: Path, model: Path) -> list[str]:
    return [
        sys.executable, str(BUILD), "--apply-only",
        "--genome", job.genome,
        "--genes", str(job.genes),
        "--operons", str(job.operons),
        "--model", str(model),
        "--synteny", "/dev/null",
        "--string", "/dev/null",
        "--string-species", "0",
        "--out-dir", str(tmp_dir),
    ]


def _run_build(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True,
                          timeout=BUILD_TIMEOUT)


def process_job(job: Job, out_dir: Path, model: Path, force: bool) -> JobResult:
    """Ejecuta un genoma. Nunca lanza: los fallos se informan en el resumen."""
    out_path = out_dir / f"model_{job.genome}.complexes.tsv"
    if out_path.is_file() and not force:
        return JobResult(genome=job.genome, operons=job.operon_source,
                         altas=_count_rows(out_path), status="omitido",
                         error="salida existente; usar --force para rehacer")
    try:
        with tempfile.TemporaryDirectory(prefix=f"cx_model_{job.genome}_") as td:
            tmp_dir = Path(td)
            cmd = build_command(job, tmp_dir, model)
            res = _run_build(cmd)
            if res.returncode != 0:
                tail = (res.stderr or res.stdout or "").strip()[-800:]
                raise RuntimeError(f"build_candidates rc={res.returncode}: {tail}")
            cluster_path = tmp_dir / f"clusters_{job.genome}.tsv"
            if not cluster_path.is_file():
                raise RuntimeError(f"falta {cluster_path.name} tras build_candidates")
            with open(cluster_path, newline="", encoding="utf-8",
                      errors="replace") as fh:
                clusters = list(csv.DictReader(fh, delimiter="\t"))
            coords = load_gene_coords(job.genes)
            rows, omitted = convert_clusters(clusters, coords)
            write_complexes(rows, out_path)
        return JobResult(genome=job.genome, operons=job.operon_source,
                         clusters=len(clusters), altas=len(rows),
                         omitidos=omitted)
    except Exception as exc:  # noqa: BLE001 - el lote sigue con otros genomas
        tmp = out_path.with_name(out_path.name + ".tmp")
        if tmp.exists():
            tmp.unlink()
        return JobResult(genome=job.genome, operons=job.operon_source,
                         fallos=1, status="fallo",
                         error=f"{type(exc).__name__}: {exc}"[:800])


def _reports_dir(data_root: Path) -> Path:
    # data/ -> reports/ en la raíz; para un data-root temporal, junto a él.
    base = data_root.parent if data_root.name == "data" else data_root
    return base / "reports"


def write_summary(results: list[JobResult], reports: Path, model: Path) -> None:
    reports.mkdir(parents=True, exist_ok=True)
    cols = ["genome", "clusters", "altas", "fallos", "omitidos", "status",
            "error"]
    tsv = reports / "complex_model_400_summary.tsv"
    with open(tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow({
                "genome": r.genome,
                "clusters": r.clusters,
                "altas": r.altas,
                "fallos": r.fallos,
                "omitidos": r.omitidos,
                "status": r.status,
                "error": r.error.replace("\t", " ").replace("\n", " "),
            })
    summary = {
        "model": str(model),
        "genomas": len(results),
        "clusters": sum(r.clusters for r in results),
        "altas": sum(r.altas for r in results),
        "fallos": sum(r.fallos for r in results),
        "omitidos": sum(r.omitidos for r in results),
        "ok": sum(r.status == "ok" for r in results),
        "saltados": sum(r.status == "omitido" for r in results),
        "resultados": [
            {"genome": r.genome, "operons": r.operons,
             "clusters": r.clusters, "altas": r.altas,
             "fallos": r.fallos, "omitidos": r.omitidos,
             "status": r.status, "error": r.error}
            for r in results
        ],
    }
    (reports / "complex_model_400_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "resultados"},
                     indent=2, ensure_ascii=False))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default: <data-root>/complexes")
    ap.add_argument("--limit", type=int, default=0,
                    help="número de genomas (0 = todos)")
    ap.add_argument("--force", action="store_true",
                    help="rehacer salidas existentes")
    ap.add_argument("--threads", type=int,
                    default=max(1, min(8, os.cpu_count() or 1)),
                    help="genomas en paralelo (subprocesos)")
    a = ap.parse_args(argv)

    if not BUILD.is_file():
        print(f"falta {BUILD}", file=sys.stderr)
        return 2
    data_root = Path(a.data_root)
    out_dir = Path(a.out_dir) if a.out_dir else data_root / "complexes"
    model = data_root / "complexes" / MODEL_NAME
    if not model.is_file():
        alt = out_dir / MODEL_NAME
        model = alt if alt.is_file() else model
    if not model.is_file():
        print(f"falta el modelo genómico {model}", file=sys.stderr)
        return 2

    jobs = discover_genomes(data_root)
    if a.limit > 0:
        jobs = jobs[:a.limit]
    if not jobs:
        print("no hay genomas con genes.tsv y operones", file=sys.stderr)
        return 2
    out_dir.mkdir(parents=True, exist_ok=True)

    threads = max(1, int(a.threads))
    results: dict[str, JobResult] = {}
    print(f"lote: {len(jobs)} genomas, threads={threads}, modelo={model.name}",
          flush=True)

    if threads == 1:
        for i, job in enumerate(jobs, 1):
            r = process_job(job, out_dir, model, a.force)
            results[r.genome] = r
            if r.status == "fallo":
                print(f"  fallo {r.genome}: {r.error}", flush=True)
            if i % 25 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}", flush=True)
    else:
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futures = {ex.submit(process_job, j, out_dir, model, a.force): j
                       for j in jobs}
            done = 0
            for fut in as_completed(futures):
                r = fut.result()
                results[r.genome] = r
                done += 1
                if r.status == "fallo":
                    print(f"  fallo {r.genome}: {r.error}", flush=True)
                if done % 25 == 0 or done == len(jobs):
                    print(f"  {done}/{len(jobs)}", flush=True)

    ordered = [results[j.genome] for j in jobs]
    write_summary(ordered, _reports_dir(data_root), model)
    ok = sum(r.status in ("ok", "omitido") for r in ordered)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
