#!/usr/bin/env python3
"""Dossier de sigma por especie a partir de salidas existentes.

No reanota: lee `sigma_factors.tsv` y `sigma_promoters.tsv`, agrega evidencia
de rutas y produce un resumen interpretable por especie.
"""
from __future__ import annotations
import argparse, csv, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def _read(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _pathways(value: str) -> list[str]:
    raw = (value or "").strip()
    if not raw or raw in {"[]", "no", "none", "null", "-"}:
        return []
    try:
        data = json.loads(raw)
    except (TypeError, ValueError):
        return [raw]
    out = []
    for item in data if isinstance(data, list) else [data]:
        if isinstance(item, dict):
            out.append(str(item.get("pathway") or item.get("name") or item.get("id") or item))
        elif item:
            out.append(str(item))
    return [x for x in out if x and x not in {"no", "None", "null"}]


def literature_status(f: dict[str, str]) -> str:
    status = (f.get("status") or "").lower()
    q = (f.get("quality_flag") or "").lower()
    family = (f.get("sigma_type") or "").lower()
    if family in {"sigma_unresolved", "", "none"}:
        return "unknown_family"
    if status in {"validated", "consensus"}:
        return "known_ortholog"
    if status == "novel" or q in {"novel_two_domains", "partial_domain"}:
        return "lineage_specific_or_partial"
    return "candidate"


def build(genome: str, data_root: Path) -> dict:
    base = data_root / "sigma" / genome
    factors = _read(base / "sigma_factors.tsv")
    promoters = _read(base / "sigma_promoters.tsv")
    by_sigma: dict[str, Counter] = defaultdict(Counter)
    for p in promoters:
        sigma = (p.get("sigma") or "").strip()
        if not sigma:
            continue
        if (p.get("pathway_support") or "").lower() == "yes":
            for pathway in _pathways(p.get("pathway_hits") or ""):
                by_sigma[sigma][pathway] += 1
    out_factors = []
    for f in factors:
        sigma = (f.get("sigma_type") or "").strip()
        paths = [p for p, _ in by_sigma.get(sigma, Counter()).most_common(3)]
        out_factors.append({
            "sigma_id": f.get("gene_id") or f.get("locus_tag") or "",
            "locus_tag": f.get("locus_tag", ""),
            "native_name": f.get("gene_name", ""),
            "family": sigma,
            "orthogroup": f.get("ortholog_group", ""),
            "product": f.get("product", ""),
            "status": literature_status(f),
            "source_status": f.get("status", ""),
            "quality_flag": f.get("quality_flag", ""),
            "confidence": f.get("confidence", ""),
            "promoter_hits": f.get("promoter_hits", ""),
            "pathways": paths,
            "hypothesis": ("asociación contextual con " + "; ".join(paths)) if paths else "",
            "note": f.get("factor_note", ""),
        })
    counts = Counter(x["status"] for x in out_factors)
    return {
        "genome": genome,
        "n_factors": len(out_factors),
        "n_known": counts.get("known_ortholog", 0),
        "n_lineage_or_partial": counts.get("lineage_specific_or_partial", 0),
        "n_unknown_family": counts.get("unknown_family", 0),
        "n_other": counts.get("candidate", 0),
        "n_promoters": len(promoters),
        "n_promoters_with_pathway": sum(1 for p in promoters if (p.get("pathway_support") or "").lower() == "yes"),
        "factors": out_factors,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "reports" / "sigma_dossiers")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    genomes = sorted(p.parent.name for p in (a.data_root / "sigma").glob("*/sigma_factors.tsv"))
    if a.limit:
        genomes = genomes[:a.limit]
    a.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for genome in genomes:
        d = build(genome, a.data_root)
        (a.out / f"{genome}.json").write_text(json.dumps(d, indent=1, ensure_ascii=False), encoding="utf-8")
        rows.append(d)
    summary = {
        "n_genomes": len(rows),
        "n_factors": sum(r["n_factors"] for r in rows),
        "n_known": sum(r["n_known"] for r in rows),
        "n_lineage_or_partial": sum(r["n_lineage_or_partial"] for r in rows),
        "n_unknown_family": sum(r["n_unknown_family"] for r in rows),
        "n_promoters": sum(r["n_promoters"] for r in rows),
        "n_promoters_with_pathway": sum(r["n_promoters_with_pathway"] for r in rows),
        "genomes_without_factors": sum(1 for r in rows if r["n_factors"] == 0),
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")
    lines = ["# Dossier sigma por especie", "",
             f"- Genomas: {summary['n_genomes']}",
             f"- Factores: {summary['n_factors']}",
             f"- Con ortólogo conocido: {summary['n_known']}",
             f"- Linaje-específicos/parciales: {summary['n_lineage_or_partial']}",
             f"- Familia desconocida: {summary['n_unknown_family']}",
             f"- Genes sin factores: {summary['genomes_without_factors']}",
             f"- Promotores candidatos: {summary['n_promoters']}",
             f"- Promotores con ruta: {summary['n_promoters_with_pathway']}", ""]
    (a.out / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
