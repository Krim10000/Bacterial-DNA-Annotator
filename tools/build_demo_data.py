"""Genera la demo incluida en el repositorio.

Toma `webapp/samples/demo.fna` (150 kb de *Chlorobaculum tepidum*), lo anota
con pyrodigal y guarda las pistas mínimas del navegador en `webapp/demo/`.

La demo es **ilustrativa**: los genes son predicciones ab initio y los operones
se agrupan por proximidad y hebra, no por evidencia experimental. Cada fila lo
declara con `evidence=demo` y `confidence=low`.

Uso:
    python tools/build_demo_data.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "webapp"))

from annotator import annotate  # noqa: E402

SAMPLE = ROOT / "webapp" / "samples" / "demo.fna"
OUT = ROOT / "webapp" / "demo"
GAP = 120  # pb entre genes para considerarlos parte del mismo operón (demo)


def build_genes(res: dict) -> list[dict]:
    genes = []
    for i, f in enumerate(res["features"], 1):
        if f.get("type") != "CDS":
            continue
        genes.append({
            "genome": "demo",
            "contig": f["contig"],
            "start": int(f["start"]),
            "end": int(f["end"]),
            "strand": "1" if str(f.get("strand")) in ("1", "+") else "-1",
            "gene_id": f"demo_{i:04d}",
            "locus_tag": f"DEMO_{i:04d}",
            "gene_name": "",
            "product": f.get("product") or "hypothetical protein",
            "ec": "",
            "source": "pyrodigal (demo)",
            "status": "demo",
            "confidence": "low",
            "evidence": "ab_initio",
        })
    genes.sort(key=lambda g: (g["contig"], g["start"]))
    return genes


def build_operons(genes: list[dict]) -> list[dict]:
    operons = []
    n = 0
    for contig in sorted({g["contig"] for g in genes}):
        group = []
        for g in [x for x in genes if x["contig"] == contig]:
            if not group:
                group = [g]
                continue
            prev = group[-1]
            same_strand = prev["strand"] == g["strand"]
            gap = g["start"] - prev["end"] - 1
            if same_strand and 0 <= gap <= GAP:
                group.append(g)
            else:
                if len(group) > 1:
                    n += 1
                    operons.append(_operon(group, n))
                group = [g]
        if len(group) > 1:
            n += 1
            operons.append(_operon(group, n))
    return operons


def _operon(group: list[dict], n: int) -> dict:
    gaps = [group[i]["start"] - group[i - 1]["end"] - 1 for i in range(1, len(group))]
    return {
        "operon_id": f"demo_op_{n:03d}",
        "label": f"demo_op_{n:03d}",
        "contig": group[0]["contig"],
        "start": min(g["start"] for g in group),
        "end": max(g["end"] for g in group),
        "strand": group[0]["strand"],
        "n_genes": len(group),
        "genes": ",".join(g["locus_tag"] for g in group),
        "gene_names": "",
        "gene_ids": ",".join(g["gene_id"] for g in group),
        "span_bp": max(g["end"] for g in group) - min(g["start"] for g in group) + 1,
        "min_gap": min(gaps) if gaps else 0,
        "max_gap": max(gaps) if gaps else 0,
        "status": "demo",
        "evidence": "proximity+strand (demo)",
        "confidence": "low",
    }


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    if not SAMPLE.exists():
        print("Falta", SAMPLE)
        return 1
    res = annotate("demo.fna", SAMPLE.read_bytes(),
                   {"cds": True, "trna": False, "rrna": False, "ghost": False,
                    "esm": False, "deserts": False})
    genes = build_genes(res)
    operons = build_operons(genes)
    write_tsv(OUT / "genes.tsv", genes, list(genes[0]))
    write_tsv(OUT / "operons.tsv", operons, list(operons[0]) if operons else
              ["operon_id", "label", "contig", "start", "end", "strand",
               "n_genes", "genes", "gene_names", "gene_ids", "span_bp",
               "min_gap", "max_gap", "status", "evidence", "confidence"])
    write_tsv(OUT / "promoters.tsv", [], [
        "operon_id", "label", "contig", "strand", "box35", "mm35", "box10",
        "mm10", "spacer", "ext10", "score", "dist_10_to_start", "n_sites",
        "sigma", "sigma_status", "status", "evidence", "confidence"])
    write_tsv(OUT / "hairpins.tsv", [], [
        "genome", "contig", "start", "end", "strand", "class", "stem", "gc",
        "loop", "polyU", "dG", "gap", "left_gene", "right_gene", "class_fixed"])
    (OUT / "README.md").write_text(
        "# Demo incluida\n\n"
        "Datos generados con `python tools/build_demo_data.py` a partir de\n"
        "`webapp/samples/demo.fna` (150 kb de *Chlorobaculum tepidum*).\n\n"
        "- `genes.tsv`: predicción ab initio con pyrodigal.\n"
        "- `operons.tsv`: agrupación **heurística** por proximidad (≤120 pb) y\n"
        "  misma hebra. No es evidencia experimental; se marca `demo`/`low`.\n"
        "- `promoters.tsv` y `hairpins.tsv`: vacíos; la demo no inventa esas\n"
        "  pistas.\n\n"
        "La demo sirve para probar el navegador y las fichas al hacer clic. No\n"
        "se usa para las métricas del proyecto ni para conclusiones biológicas.\n",
        encoding="utf-8",
    )
    print(f"genes={len(genes)} operones={len(operons)} -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
