#!/usr/bin/env python3
"""Anota con el pipeline los genomas de referencia descargados.

Escribe `results/<id>/` y una vista mínima en `data/regulation/<id>/` para que
el navegador pueda abrirlos. No descarga nada.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from webapp.pipeline_adapter import run_in_memory  # noqa: E402


def _write_genes(res: dict, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cols = ["genome", "contig", "start", "end", "strand", "gene_id",
            "locus_tag", "gene_name", "product", "ec", "acc", "gene_type",
            "source", "confidence"]
    rows = []
    for f in res.get("features") or []:
        if f.get("type") not in ("CDS", "tRNA", "rRNA", "ncRNA"):
            continue
        rows.append({
            "genome": res.get("genome", ""), "contig": f.get("contig", ""),
            "start": f.get("start", ""), "end": f.get("end", ""),
            "strand": f.get("strand", ""), "gene_id": f.get("id", ""),
            "locus_tag": f.get("id", ""), "gene_name": "",
            "product": f.get("product", ""), "ec": "",
            "acc": "", "gene_type": f.get("type", ""),
            "source": "pipeline", "confidence": "hypothesis",
        })
    with (outdir / "genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--level", type=int, choices=(1, 2), default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    mp = ROOT / "data" / "reference_genomes.json"
    rows = json.loads(mp.read_text(encoding="utf-8"))
    todo = [r for r in rows if r.get("status") == "downloaded"
            and (a.level is None or r.get("level") == a.level)]
    if a.limit:
        todo = todo[:a.limit]
    print(f"a anotar: {len(todo)}")
    for i, rec in enumerate(todo, 1):
        rid = rec["id"]
        fna = ROOT / rec["genome_path"]
        out = ROOT / "results" / rid
        if out.joinpath(f"{rid}.summary.json").exists() and not a.force:
            print(f"  [{i}/{len(todo)}] {rid}: ya anotado")
            continue
        print(f"  [{i}/{len(todo)}] {rid}: pipeline...", flush=True)
        try:
            res = run_in_memory(f"{rid}.fna", fna.read_bytes(), out,
                                ROOT / "data",
                                {"functional": True, "metabolic": False,
                                 "ncrna": True, "membrane": True, "oric": True,
                                 "taxonomy": True, "tss_tts": True},
                                threads=4)
            _write_genes(res, ROOT / "data" / "regulation" / rid)
            rec.update({"status": "annotated",
                        "system_summary": str((out / f"{rid}.summary.json").relative_to(ROOT)),
                        "system_result": str(out.relative_to(ROOT))})
            print(f"    ok CDS={len([f for f in res.get('features',[]) if f.get('type')=='CDS'])}", flush=True)
        except Exception as exc:  # noqa: BLE001
            rec.update({"status": "annotate_error", "error": f"{type(exc).__name__}: {exc}"})
            print(f"    ERROR {exc}", flush=True)
        mp.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
