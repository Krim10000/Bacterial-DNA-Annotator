#!/usr/bin/env python3
"""18b_rebuild_rna_strand.py — reconstruye rna_annotations.json desde el GFF de NCBI.

Motivo
------
`data/rna_annotations/rna_annotations.json` guardaba solo `start`/`end` por
feature: sin hebra, sin tipo específico y sin producto. Para unidades de
transcripción eso es inservible (no se sabe si un tRNA corta el operón o está
en la hebra contraria) y para el resto del pipeline pierde el isotipo y el Rfam.

Este script lo regenera desde los GFF de NCBI, que sí traen hebra, `product`,
`anticodon`, `locus_tag` y `Dbxref` (Rfam). No vuelve a correr cmsearch ni
tRNAscan: solo re-lee lo que ya está descargado.

Salida (misma ruta, esquema retrocompatible)
-------------------------------------------
`data/rna_annotations/rna_annotations.json`:

    {"<genoma>": {"rRNA": [...], "tRNA": [...], "ncRNA": [...]}}

Cada feature conserva `start`/`end` (los únicos campos que usaban los
consumidores actuales) y añade `strand`, `type`, `product`, `gene`,
`locus_tag`, `seqid`, `anticodon`, `rfam`, `partial`, `source`.

`ncRNA` agrupa tmRNA, SRP_RNA, RNase_P_RNA, riboswitch y ncRNA; el campo `type`
conserva la clase original.

Además escribe `data/rna_annotations/meta.json` (recuentos, cobertura de hebra,
hash, fuente) y `data/rna_annotations/MIGRATION.md` (qué cambió frente al
fichero anterior).

Uso
---
    python scripts/annotation/18b_rebuild_rna_strand.py --dry-run
    python scripts/annotation/18b_rebuild_rna_strand.py
    python scripts/annotation/18b_rebuild_rna_strand.py --limit 5 --verbose
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GFF_DIR = ROOT / "data" / "annotations"
OUT = ROOT / "data" / "rna_annotations" / "rna_annotations.json"
META = ROOT / "data" / "rna_annotations" / "meta.json"
MIGRATION = ROOT / "data" / "rna_annotations" / "MIGRATION.md"

GROUPS = {"tRNA": "tRNA", "rRNA": "rRNA"}
NC_TYPES = {"ncRNA", "tmRNA", "SRP_RNA", "RNase_P_RNA", "riboswitch",
            "antisense_RNA", "sRNA", "misc_RNA"}
WANTED = set(GROUPS) | NC_TYPES


def parse_attrs(field: str) -> dict[str, str]:
    """GFF3 columna 9 -> dict. Solo parte en el primer '=' de cada par."""
    out: dict[str, str] = {}
    for part in field.strip().split(";"):
        if not part or "=" not in part:
            continue
        key, val = part.split("=", 1)
        out[key.strip()] = urllib.parse.unquote(val.strip())
    return out


def group_of(ftype: str) -> str | None:
    if ftype in GROUPS:
        return GROUPS[ftype]
    if ftype in NC_TYPES:
        return "ncRNA"
    return None


def parse_gff(path: Path) -> tuple[dict, dict]:
    """Devuelve (features por grupo, stats del genoma)."""
    genes: dict[str, dict] = {}
    rows: list[tuple[str, dict]] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = parts
            if ftype == "gene":
                a = parse_attrs(attrs)
                gid = a.get("ID", "")
                if gid:
                    genes[gid] = a
                continue
            if ftype not in WANTED:
                continue
            rows.append((ftype, {"seqid": seqid, "start": int(start), "end": int(end),
                                 "strand": strand, "attrs": parse_attrs(attrs)}))

    feats: dict[str, list[dict]] = {"rRNA": [], "tRNA": [], "ncRNA": []}
    stats = {"no_strand": 0, "no_product": 0}
    for ftype, r in rows:
        a = r["attrs"]
        parent = genes.get(a.get("Parent", ""), {})
        product = a.get("product") or a.get("Name") or parent.get("Name") or ftype
        if not a.get("product"):
            stats["no_product"] += 1
        if r["strand"] not in "+-":
            stats["no_strand"] += 1
        rec = {
            "start": r["start"], "end": r["end"],
            "strand": r["strand"] if r["strand"] in "+-" else "?",
            "type": ftype,
            "product": product,
            "gene": a.get("gene") or parent.get("gene") or "",
            "locus_tag": a.get("locus_tag") or parent.get("locus_tag") or "",
            "seqid": r["seqid"],
            "anticodon": a.get("anticodon", ""),
            "rfam": next((d.split(":")[1] for d in a.get("Dbxref", "").split(",")
                          if d.startswith("RFAM:")), ""),
            "partial": (a.get("partial") == "true" or a.get("pseudo") == "true"),
            "source": "NCBI_GFF",
        }
        feats[group_of(ftype)].append(rec)

    for group in feats:
        feats[group].sort(key=lambda x: (x["seqid"], x["start"], x["end"]))
    return feats, stats


def find_gffs(limit: int | None) -> list[Path]:
    paths = sorted(GFF_DIR.glob("*.gff.gz"))
    return paths[:limit] if limit else paths


def genome_id(path: Path) -> str:
    return path.name.replace(".gff.gz", "")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def old_counts(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    try:
        d = json.load(open(path))
    except Exception:
        return {}
    return {k: sum(len(v) for v in val.values()) for k, val in d.items()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--baseline", type=Path, default=None,
                    help="fichero previo con el que comparar la migración "
                         "(por defecto, el propio destino antes de escribir)")
    args = ap.parse_args()

    gffs = find_gffs(args.limit)
    if not gffs:
        print("no hay GFF en data/annotations/", file=sys.stderr)
        return 1

    baseline = args.baseline or OUT
    before = old_counts(baseline)
    baseline_sha = sha256(baseline) if baseline.exists() else ""
    result: dict[str, dict] = {}
    totals = {"rRNA": 0, "tRNA": 0, "ncRNA": 0}
    no_strand = no_product = 0
    for i, path in enumerate(gffs, 1):
        feats, st = parse_gff(path)
        gid = genome_id(path)
        result[gid] = feats
        for g, items in feats.items():
            totals[g] += len(items)
        no_strand += st["no_strand"]
        no_product += st["no_product"]
        if args.verbose:
            print(f"  [{i}/{len(gffs)}] {gid}: "
                  f"tRNA={len(feats['tRNA'])} rRNA={len(feats['rRNA'])} "
                  f"ncRNA={len(feats['ncRNA'])}")
        elif i % 100 == 0:
            print(f"  [{i}/{len(gffs)}] ...")

    n_feats = sum(totals.values())
    with_strand = sum(1 for v in result.values() for items in v.values()
                      for it in items if it["strand"] in "+-")
    genomes_with_rrna = sum(1 for v in result.values() if v["rRNA"])
    genomes_with_trna = sum(1 for v in result.values() if v["tRNA"])

    print(f"\ngenomas: {len(result)} | features: {n_feats}")
    print(f"  tRNA {totals['tRNA']} (en {genomes_with_trna} genomas) | "
          f"rRNA {totals['rRNA']} (en {genomes_with_rrna}) | "
          f"ncRNA {totals['ncRNA']}")
    print(f"  con hebra: {with_strand}/{n_feats} "
          f"({100 * with_strand / n_feats:.1f}%) | sin producto: {no_product}")

    if args.dry_run:
        print("dry-run: no se escribe nada")
        return 0

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")))
    digest = sha256(OUT)

    after = {k: sum(len(v) for v in val.values()) for k, val in result.items()}
    gained = sum(1 for k in after if k in before and after[k] > before[k])
    lost = sum(1 for k in after if k in before and after[k] < before[k])
    new_genomes = len([k for k in after if k not in before])

    meta = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "data/annotations/*.gff.gz (GFF3 de NCBI)",
        "script": "scripts/annotation/18b_rebuild_rna_strand.py",
        "n_genomes": len(result),
        "n_features": n_feats,
        "by_group": totals,
        "with_strand": with_strand,
        "strand_coverage": round(with_strand / n_feats, 6) if n_feats else 0,
        "without_product": no_product,
        "genomes_with_trna": genomes_with_trna,
        "genomes_with_rrna": genomes_with_rrna,
        "sha256": digest,
        "note": "derivado del GFF; no se re-ejecutó cmsearch ni tRNAscan-SE",
    }
    META.write_text(json.dumps(meta, indent=1))

    lines = [
        "# Migración de rna_annotations.json",
        "",
        f"Generado por `18b_rebuild_rna_strand.py` el {meta['generated_at']}.",
        "",
        "Antes: `{start, end}`. Después: `{start, end, strand, type, product, "
        "gene, locus_tag, seqid, anticodon, rfam, partial, source}`.",
        "",
        "| métrica | valor |",
        "|---|---|",
        f"| genomas | {len(result)} |",
        f"| features | {n_feats} |",
        f"| tRNA / rRNA / ncRNA | {totals['tRNA']} / {totals['rRNA']} / {totals['ncRNA']} |",
        f"| con hebra | {with_strand}/{n_feats} ({100 * with_strand / n_feats:.1f}%) |",
        f"| sin producto | {no_product} |",
        f"| genomas con conteo mayor | {gained} |",
        f"| genomas con conteo menor | {lost} |",
        f"| genomas nuevos | {new_genomes} |",
        f"| sha256 | `{digest[:16]}…` |",
        f"| baseline | `{baseline_sha[:16]}…` ({len(before)} genomas) |",
        "",
        "Retrocompatible: los consumidores que leen `start`/`end` siguen "
        "funcionando; los campos nuevos son aditivos.",
        "",
    ]
    MIGRATION.write_text("\n".join(lines))
    print(f"escrito: {OUT.relative_to(ROOT)}")
    print(f"escrito: {META.relative_to(ROOT)}")
    print(f"escrito: {MIGRATION.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
