#!/usr/bin/env python3
"""03_scan_pfam.py — hmmscan del panel curado contra los proteomas ab initio.

Usa los protein_id del GFF ab initio (mismos identificadores que `proteins.faa`),
de modo que las coordenadas de cada hit son directamente las del CDS ab initio.

Se corre sin `--cut_ga` y con `--domE 1e-3`: se guarda el mejor dominio por
(proteína, modelo), se marca `passes_ga` si supera el umbral GA de Pfam y esa
marca se usa luego para graduar la confianza (un hit débil solo puede apoyar,
nunca sostener por sí solo una llamada).

Salida: data/mobile_typing/_work/<genoma>/pfam_hits.tsv
        data/mobile_typing/_work/<genoma>/pfam_scan_stats.json
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANN = ROOT / "data" / "annotations_ab_initio"
OUT = ROOT / "data" / "mobile_typing"
WORK = OUT / "_work"
PANEL = WORK / "pfam_panel" / "panel.hmm"
MANIFEST = WORK / "pfam_panel" / "panel_manifest.tsv"


def load_manifest(path=None):
    meta = {}
    with (path or MANIFEST).open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["accession"]] = {
                "name": row["name"],
                "desc": row["description"],
                "categories": row["categories"],
                "ga_seq": float(row["ga_seq"] or 0),
                "ga_dom": float(row["ga_dom"] or 0),
            }
    return meta


def parse_domtbl(path: Path):
    """Mejor dominio por (proteína, modelo): menor i-Evalue, luego mayor score."""
    best = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            # domtblout: target=modelo, query=proteína
            f = line.split()
            if len(f) < 22:
                continue
            target_acc = f[1].split(".")[0]
            query = f[3]
            ievalue = float(f[12])
            score = float(f[13])
            dom_score = float(f[13])
            hmm_from, hmm_to = int(f[15]), int(f[16])
            ali_from, ali_to = int(f[17]), int(f[18])
            env_from, env_to = int(f[19]), int(f[20])
            key = (query, target_acc)
            cand = {
                "target_acc": target_acc,
                "query": query,
                "ievalue": ievalue,
                "score": score,
                "dom_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "env_from": env_from,
                "env_to": env_to,
            }
            prev = best.get(key)
            if prev is None or (ievalue, -score) < (prev["ievalue"], -prev["score"]):
                best[key] = cand
    return list(best.values())


def main() -> int:
    if not PANEL.exists():
        print("ERROR: falta el panel; ejecuta 02_build_pfam_panel.py", file=sys.stderr)
        return 1
    meta = load_manifest()
    pilots = [row["genome"] for row in csv.DictReader((OUT / "pilot_genomes.tsv").open(), delimiter="\t")]
    for genome in pilots:
        work = WORK / genome
        work.mkdir(parents=True, exist_ok=True)
        faa = ANN / genome / "proteins.faa"
        domtbl = work / "pfam.domtbl"
        cmd = [
            "hmmscan", "--cpu", "8", "--noali", "--domE", "1e-3",
            "--domtblout", str(domtbl), "-o", "/dev/null", str(PANEL), str(faa),
        ]
        print(f"[{genome}] hmmscan …", flush=True)
        subprocess.run(cmd, check=True)
        hits = parse_domtbl(domtbl)
        rows = []
        n_ga = 0
        for hit in hits:
            info = meta[hit["target_acc"]]
            passes = hit["dom_score"] >= info["ga_dom"] > 0
            n_ga += int(passes)
            rows.append(
                {
                    "protein_id": hit["query"],
                    "pfam_acc": hit["target_acc"],
                    "pfam_name": info["name"],
                    "pfam_desc": info["desc"],
                    "categories": info["categories"],
                    "passes_ga": "yes" if passes else "no",
                    "ievalue": f"{hit['ievalue']:.3g}",
                    "score": f"{hit['score']:.1f}",
                    "hmm_from": hit["hmm_from"],
                    "hmm_to": hit["hmm_to"],
                    "ali_from": hit["ali_from"],
                    "ali_to": hit["ali_to"],
                    "env_from": hit["env_from"],
                    "env_to": hit["env_to"],
                }
            )
        rows.sort(key=lambda r: (r["protein_id"], r["pfam_acc"]))
        with (work / "pfam_hits.tsv").open("w", newline="") as fh:
            fields = list(rows[0].keys()) if rows else [
                "protein_id", "pfam_acc", "pfam_name", "pfam_desc", "categories",
                "passes_ga", "ievalue", "score", "hmm_from", "hmm_to", "ali_from",
                "ali_to", "env_from", "env_to",
            ]
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        stats = {
            "proteins": sum(1 for line in faa.open() if line.startswith(">")),
            "hits": len(hits),
            "hits_passing_ga": n_ga,
            "models": len(meta),
        }
        (work / "pfam_scan_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
        print(f"  hits={stats['hits']} (GA={stats['hits_passing_ga']}) proteins={stats['proteins']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
