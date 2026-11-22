#!/usr/bin/env python3
"""03 - DIAMOND all-hits de las proteínas del genoma contra la base de transporte.

No se queda con el mejor hit por proteína: conserva todos los hits por consulta
(evalue <= 1e-3, sin límite de targets) para poder buscar, gap por gap, el acc que esté ligado a la
reacción Rhea exacta del metabolito.

Entrada:  data/metabolism/transport_t1/pilot_genomes.tsv
          data/annotations_ab_initio/<genome>/proteins.faa (solo lectura)
          data/db/transport.dmnd
Salida:   data/metabolism/transport_t1/diamond/<genome>.hits.tsv
          data/metabolism/transport_t1/diamond/run_stats.json
"""
from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
common = importlib.import_module("00_common")

OUTFMT = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "evalue",
    "bitscore",
    "qcovhsp",
    "stitle",
]


def diamond_version() -> str:
    try:
        out = subprocess.run(
            ["diamond", "version"], capture_output=True, text=True, check=True
        )
        return (out.stdout or out.stderr).strip().splitlines()[0]
    except Exception as exc:  # pragma: no cover - entorno sin DIAMOND
        return f"no disponible ({exc})"


def validate_hits(path: Path, query_ids: set) -> int:
    n = 0
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(OUTFMT):
                raise SystemExit(f"ERROR {path}:{lineno}: {len(fields)} columnas")
            if fields[0] not in query_ids:
                raise SystemExit(f"ERROR {path}:{lineno}: qseqid desconocido {fields[0]!r}")
            if not fields[1].startswith("sp|"):
                raise SystemExit(f"ERROR {path}:{lineno}: sseqid inesperado {fields[1]!r}")
            float(fields[2]), float(fields[6]), float(fields[7]), float(fields[8])
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force", action="store_true", help="recalcula aunque exista la salida")
    args = ap.parse_args()

    if not common.TRANSPORT_DMND.exists():
        raise SystemExit(f"ERROR: falta {common.TRANSPORT_DMND}")
    if subprocess.run(["which", "diamond"], capture_output=True).returncode != 0:
        raise SystemExit("ERROR: `diamond` no está en PATH")

    pilot = list(common.read_tsv(common.OUT / "pilot_genomes.tsv"))
    if len(pilot) != 10:
        raise SystemExit(f"ERROR: se esperaban 10 pilotos, hay {len(pilot)}")
    common.DIAMOND_DIR.mkdir(parents=True, exist_ok=True)

    version = diamond_version()
    stats = {
        "generated_by": "scripts/metabolism/03_run_diamond.py",
        "diamond": version,
        "database": "data/db/transport.dmnd",
        "query": "data/annotations_ab_initio/<genome>/proteins.faa",
        "parameters": {
            "outfmt": "6 " + " ".join(OUTFMT),
            "evalue": common.SEARCH_EVALUE,
            "max_target_seqs": "sin límite (0)",
            "threads": args.threads,
            "note": "all-hits sin límite por proteína; no solo el mejor hit",
        },
        "genomes": [],
    }

    for row in pilot:
        genome = row["genome"]
        faa = common.ANN / genome / "proteins.faa"
        if not faa.exists():
            raise SystemExit(f"ERROR: falta {faa}")
        query_ids = common.fasta_ids(faa)
        if not query_ids:
            raise SystemExit(f"ERROR: {faa} sin proteínas")
        if len(query_ids) != len(set(query_ids)):
            raise SystemExit(f"ERROR: ids duplicados en {faa}")

        out = common.DIAMOND_DIR / f"{genome}.hits.tsv"
        t0 = time.time()
        reused = False
        if out.exists() and out.stat().st_size > 0 and not args.force:
            reused = True
        else:
            tmp = out.with_suffix(".hits.tsv.tmp")
            cmd = [
                "diamond",
                "blastp",
                "--query",
                str(faa),
                "--db",
                str(common.TRANSPORT_DMND),
                "--out",
                str(tmp),
                "--outfmt",
                "6",
                *OUTFMT,
                "--evalue",
                str(common.SEARCH_EVALUE),
                "--max-target-seqs",
                str(common.MAX_TARGET_SEQS),
                "--threads",
                str(args.threads),
                "--quiet",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                if tmp.exists():
                    tmp.unlink()
                raise SystemExit(
                    f"ERROR DIAMOND {genome}:\n{proc.stdout}\n{proc.stderr}"
                )
            tmp.replace(out)

        n_hits = validate_hits(out, set(query_ids))
        elapsed = time.time() - t0
        stats["genomes"].append(
            {
                "genome": genome,
                "proteins": len(query_ids),
                "hits": n_hits,
                "hits_file": f"diamond/{genome}.hits.tsv",
                "reused": reused,
                "seconds": round(elapsed, 1),
            }
        )
        print(
            f"[03] {genome}: {len(query_ids)} proteínas, {n_hits} hits "
            f"({'reusado' if reused else f'{elapsed:.1f}s'})",
            flush=True,
        )

    total_hits = sum(g["hits"] for g in stats["genomes"])
    if total_hits == 0:
        raise SystemExit("ERROR: DIAMOND no produjo ningún hit")
    stats["totals"] = {
        "genomes": len(stats["genomes"]),
        "proteins": sum(g["proteins"] for g in stats["genomes"]),
        "hits": total_hits,
    }
    common.write_json(common.DIAMOND_DIR / "run_stats.json", stats)
    print(f"[03] total: {total_hits} hits en {len(stats['genomes'])} genomas")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
