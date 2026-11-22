#!/usr/bin/env python3
"""14_parse_gff.py — parsea los GFF3 de NCBI a un JSON de etiquetas biológicas.

Qué cambia respecto a la versión anterior
-----------------------------------------
1. **Conserva los atributos.** Antes solo guardaba `start`/`end`/`strand` y
   descartaba `locus_tag`, `gene`, `product`, EC y `protein_id`: las etiquetas
   no se podían unir a un gen ni a una reacción. Ahora se conservan todos, más
   `seqid` y `partial`.
2. **Extrae EC** de `Dbxref=EC:...` y de `product="... (EC x.y.z.w)"`.
3. **No finge promotores.** NCBI no emite features de promotor en sus GFF; la
   versión anterior buscaba `regulatory`/`promoter` y reportaba `0 promotores`,
   que se lee como "este genoma no tiene promotores". Ahora `promoter` va vacío
   con `_meta.promoter_reason` explicando por qué y a qué capa pertenece
   (L2: `scripts/regulation/12_promoters_sigma.py`).
4. **Idempotente y determinista.** Orden estable, sin red.

Esquema de salida (`data/parsed_annotations/<genoma>_labels.json`)
----------------------------------------------------------------
    {
      "CDS":    [{start, end, strand, seqid, locus_tag, gene, product,
                  ec: [...], protein_id, partial}],
      "tRNA":   [...], "rRNA": [...], "ncRNA": [...],
      "promoter": [],
      "_meta":  {source, counts, ec_coverage, promoter_reason, script}
    }

Las claves de primer nivel se mantienen para no romper consumidores que lean
`d["CDS"]`; `_meta` es aditivo.

Uso
---
    python scripts/annotation/14_parse_gff.py --limit 5 --verbose
    python scripts/annotation/14_parse_gff.py
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANN_DIR = ROOT / "data" / "annotations"
OUT_DIR = ROOT / "data" / "parsed_annotations"

CDS_KEYS = ("CDS",)
RNA_GROUPS = {"tRNA": "tRNA", "rRNA": "rRNA"}
RNA_NC = {"ncRNA", "tmRNA", "SRP_RNA", "RNase_P_RNA", "riboswitch",
          "antisense_RNA", "sRNA", "misc_RNA"}
EC_RE = re.compile(r"EC[:\s]*(\d+\.\d+\.\d+\.(?:\d+|-))", re.IGNORECASE)
PROMOTER_REASON = (
    "NCBI no emite features de promotor en sus GFF (no hay líneas "
    "'regulatory' con promotor). Un 0 aquí no significa 'sin promotores'. "
    "La capa L2 los predice en scripts/regulation/12_promoters_sigma.py."
)


def parse_attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in field.strip().split(";"):
        if not part or "=" not in part:
            continue
        key, val = part.split("=", 1)
        out[key.strip()] = urllib.parse.unquote(val.strip())
    return out


def extract_ec(attrs: dict[str, str]) -> list[str]:
    found: set[str] = set()
    for key in ("Dbxref", "db_xref", "ec_number", "EC_number", "product"):
        for m in EC_RE.finditer(attrs.get(key, "")):
            found.add(m.group(1))
    return sorted(found)


def gbff_ec_map(path: Path) -> dict[str, list[str]]:
    """locus_tag -> EC, desde los qualifiers /EC_number del GenBank.

    Parser de líneas, no Biopython: SeqIO sobre 400 GBFF de ~4 MB tardaba
    ~24 s/genoma (>2,5 h); esto es I/O puro y baja a ~0,2 s.

    El GFF de NCBI de este set no trae EC (0 coincidencias), pero el GBFF sí
    (~700 por genoma). Sin esto la etiqueta EC queda vacía por la fuente.
    """
    if not path.exists():
        return {}
    out: dict[str, list[str]] = {}
    tag, ecs = None, set()

    def flush() -> None:
        if tag and ecs:
            out[tag] = sorted(ecs)

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("     "):
                continue
            body = line[5:].lstrip()
            if body.startswith("/"):
                q = body[1:].strip()
                if q.startswith("locus_tag="):
                    tag = q.split("=", 1)[1].strip().strip('"')
                elif q.startswith("EC_number="):
                    ecs.add(q.split("=", 1)[1].strip().strip('"').rstrip(";"))
            else:
                key = body.split()[0] if body.split() else ""
                if key != "CDS" or not line[5:].strip():
                    flush()
                    tag, ecs = None, set()
    flush()
    return out


def group_of(ftype: str) -> str | None:
    if ftype in RNA_GROUPS:
        return RNA_GROUPS[ftype]
    if ftype in RNA_NC:
        return "ncRNA"
    return None


def rna_record(ftype: str, seqid: str, start: int, end: int, strand: str,
               a: dict[str, str]) -> dict:
    return {
        "start": start, "end": end,
        "strand": strand if strand in "+-" else "?",
        "seqid": seqid, "type": ftype,
        "locus_tag": a.get("locus_tag", ""),
        "gene": a.get("gene", ""),
        "product": a.get("product") or a.get("Name") or ftype,
        "anticodon": a.get("anticodon", ""),
        "partial": a.get("partial") == "true" or a.get("pseudo") == "true",
    }


def parse_gff(path: Path, ec_map: dict[str, list[str]] | None = None) -> dict:
    feats: dict[str, list] = {"CDS": [], "tRNA": [], "rRNA": [], "ncRNA": [],
                              "promoter": []}
    ec_map = ec_map or {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _src, ftype, start, end, _s, strand, _ph, raw = parts
            a = parse_attrs(raw)
            if ftype == "CDS":
                ec = extract_ec(a)
                ec_source = "gff" if ec else ""
                if not ec:
                    tag = a.get("locus_tag", "")
                    if tag in ec_map:
                        ec, ec_source = ec_map[tag], "gbff"
                feats["CDS"].append({
                    "start": int(start), "end": int(end),
                    "strand": strand if strand in "+-" else "?",
                    "seqid": seqid,
                    "locus_tag": a.get("locus_tag", ""),
                    "gene": a.get("gene", "") or a.get("Name", ""),
                    "product": a.get("product", ""),
                    "ec": ec,
                    "ec_source": ec_source,
                    "protein_id": a.get("protein_id", ""),
                    "partial": (a.get("partial") == "true"
                                or a.get("pseudo") == "true"),
                })
            elif group_of(ftype):
                feats[group_of(ftype)].append(
                    rna_record(ftype, seqid, int(start), int(end), strand, a))
    for key in ("CDS", "tRNA", "rRNA", "ncRNA"):
        feats[key].sort(key=lambda x: (x["seqid"], x["start"], x["end"]))
    feats["_meta"] = {
        "source": str(path.relative_to(ROOT)),
        "script": "scripts/annotation/14_parse_gff.py",
        "counts": {k: len(v) for k, v in feats.items() if isinstance(v, list)},
        "ec_coverage": round(sum(1 for c in feats["CDS"] if c["ec"])
                             / len(feats["CDS"]), 4) if feats["CDS"] else 0.0,
        "ec_note": "el GFF de NCBI no trae EC; se enriquece desde /EC_number "
                   "del GBFF cuando existe",
        "promoter_reason": PROMOTER_REASON,
    }
    return feats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--no-gbff", action="store_true",
                    help="no enriquece EC desde el GenBank")
    args = ap.parse_args()

    gffs = sorted(ANN_DIR.glob("*.gff.gz"))
    if args.limit:
        gffs = gffs[:args.limit]
    if not gffs:
        print(f"no hay GFF en {ANN_DIR}", file=sys.stderr)
        return 1

    if not args.dry_run:
        OUT_DIR.mkdir(parents=True, exist_ok=True)

    totals = {"CDS": 0, "tRNA": 0, "rRNA": 0, "ncRNA": 0}
    cds_with_ec = 0
    for i, gff in enumerate(gffs, 1):
        genome = gff.name.replace(".gff.gz", "")
        ec_map = {}
        if not args.no_gbff:
            gbff = ANN_DIR / f"{genome}_genomic.gbff.gz"
            ec_map = gbff_ec_map(gbff)
        feats = parse_gff(gff, ec_map)
        for k in totals:
            totals[k] += len(feats[k])
        cds_with_ec += sum(1 for c in feats["CDS"] if c["ec"])
        if args.verbose:
            c = feats["_meta"]["counts"]
            print(f"  [{i}/{len(gffs)}] {genome}: CDS={c['CDS']} "
                  f"tRNA={c['tRNA']} rRNA={c['rRNA']} ncRNA={c['ncRNA']} "
                  f"EC={feats['_meta']['ec_coverage']:.1%} (gbff {len(ec_map)})")
        elif i % 100 == 0:
            print(f"  [{i}/{len(gffs)}] ...")
        if not args.dry_run:
            (OUT_DIR / f"{genome}_labels.json").write_text(
                json.dumps(feats, sort_keys=True, separators=(",", ":")))

    n_cds = totals["CDS"]
    print(f"\ngenomas: {len(gffs)} | CDS {n_cds} | tRNA {totals['tRNA']} | "
          f"rRNA {totals['rRNA']} | ncRNA {totals['ncRNA']}")
    print(f"  CDS con EC: {cds_with_ec} "
          f"({100 * cds_with_ec / n_cds:.1f}%)" if n_cds else "")
    print("  promotores: no aplicable — NCBI no los emite (ver _meta.promoter_reason)")
    if args.dry_run:
        print("dry-run: no se escribió nada")
    else:
        print(f"escrito: {OUT_DIR.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
