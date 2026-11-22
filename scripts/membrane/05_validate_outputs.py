#!/usr/bin/env python3
"""Valida TSV y GFF3 de data/membrane contra el GFF ab initio (solo lectura).

Comprueba:
  - cabecera exacta del TSV y tipos/confianzas permitidos;
  - coordenadas 1-based dentro del contig y hebra +/-;
  - protein_id presente en proteins.faa;
  - GFF3 parseable (9 columnas, atributos, IDs unicos, sequence-region);
  - correspondencia 1:1 entre filas TSV y features GFF.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import re
import sys
import urllib.parse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ALLOWED_TYPES = {"tm", "signal_peptide", "lipoprotein", "transporter", "secretion"}
ALLOWED_CONF = {"high", "medium", "low"}
HEADER = ["genome", "contig", "start", "end", "strand", "protein_id", "type",
          "family", "substrate", "evidence", "confidence"]
ATTR_TAG = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def read_pilot(path: Path) -> list[str]:
    genomes = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        if line.startswith("genome\t"):
            continue
        genomes.append(line.split("\t")[0])
    return genomes


def gff_contigs(path: Path) -> dict[str, int]:
    contigs = {}
    for line in open(path):
        if line.startswith("##sequence-region"):
            p = line.split()
            contigs[p[1]] = int(p[3])
    return contigs


def parse_attrs(text: str):
    attrs = {}
    for part in text.split(";"):
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"atributo sin '=': {part!r}")
        k, v = part.split("=", 1)
        if not ATTR_TAG.match(k):
            raise ValueError(f"tag invalido: {k!r}")
        if k in attrs:
            raise ValueError(f"tag duplicado: {k}")
        attrs[k] = urllib.parse.unquote(v)
    return attrs


def parse_gff(path: Path):
    seqregions = {}
    features = []
    errors = []
    first = None
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if first is None and line.strip():
                first = line.strip()
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) != 4:
                    errors.append(f"GFF linea {ln}: sequence-region invalido")
                    continue
                seqregions[p[1]] = int(p[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                errors.append(f"GFF linea {ln}: {len(f)} columnas (se esperan 9)")
                continue
            try:
                start, end = int(f[3]), int(f[4])
                attrs = parse_attrs(f[8])
            except ValueError as exc:
                errors.append(f"GFF linea {ln}: {exc}")
                continue
            features.append(dict(seqid=f[0], source=f[1], type=f[2], start=start,
                                 end=end, score=f[5], strand=f[6], phase=f[7],
                                 attrs=attrs, line=ln))
    return first, seqregions, features, errors


def validate_genome(genome: str, tsv_dir: Path, gff_dir: Path, annot_root: Path,
                    allow_missing: bool = False) -> dict:
    errors: list[str] = []
    warnings: list[str] = []
    tsv = tsv_dir / f"{genome}.tsv"
    gff = gff_dir / f"{genome}.gff3"
    annot = annot_root / genome
    if not tsv.exists() or not gff.exists():
        miss = [str(p) for p in (tsv, gff) if not p.exists()]
        if allow_missing:
            return dict(ok=True, pending=True, rows=0, gff_features=0, errors=[],
                        warnings=[f"pendiente: faltan {miss}"])
        return dict(ok=False, pending=False, rows=0, gff_features=0,
                    errors=[f"falta {m}" for m in miss], warnings=[])

    contigs = gff_contigs(annot / "annotation.gff3")
    proteins = {line[1:].split()[0] for line in open(annot / "proteins.faa")
                if line.startswith(">")}

    rows = []
    with open(tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        try:
            header = next(rd)
        except StopIteration:
            return dict(ok=False, rows=0, gff_features=0, errors=["TSV vacio"], warnings=[])
        if header != HEADER:
            errors.append(f"cabecera TSV inesperada: {header}")
        for ln, f in enumerate(rd, 2):
            if len(f) != len(HEADER):
                errors.append(f"TSV linea {ln}: {len(f)} columnas")
                continue
            rows.append(dict(zip(HEADER, f), line=ln))
    for r in rows:
        where = f"{r['protein_id']} ({r['type']})"
        if r["genome"] != genome:
            errors.append(f"{where}: genome {r['genome']} != {genome}")
        if r["contig"] not in contigs:
            errors.append(f"{where}: contig desconocido {r['contig']}")
            continue
        try:
            start, end = int(r["start"]), int(r["end"])
        except ValueError:
            errors.append(f"{where}: coordenadas no enteras")
            continue
        if not (1 <= start <= end <= contigs[r["contig"]]):
            errors.append(f"{where}: coordenadas {start}-{end} fuera de 1-{contigs[r['contig']]}")
        if r["strand"] not in {"+", "-"}:
            errors.append(f"{where}: hebra invalida {r['strand']}")
        if r["protein_id"] not in proteins:
            errors.append(f"{where}: protein_id no esta en proteins.faa")
        if r["type"] not in ALLOWED_TYPES:
            errors.append(f"{where}: tipo invalido {r['type']}")
        if r["confidence"] not in ALLOWED_CONF:
            errors.append(f"{where}: confianza invalida {r['confidence']}")
        if not r["family"].strip():
            errors.append(f"{where}: family vacia")
        if not r["evidence"].strip():
            errors.append(f"{where}: evidence vacia")
        if any(ch in r["evidence"] for ch in "\t\n\r"):
            errors.append(f"{where}: evidence con tabulador/salto")

    dup = [k for k, v in Counter((r["protein_id"], r["type"]) for r in rows).items() if v > 1]
    if dup:
        errors.append(f"filas duplicadas (protein_id,type): {dup[:5]}")

    first, seqregions, feats, gff_errors = parse_gff(gff)
    errors.extend(gff_errors)
    if first != "##gff-version 3":
        errors.append(f"GFF no empieza con ##gff-version 3 (empieza {first!r})")
    if set(seqregions) != set(contigs):
        errors.append(f"GFF sequence-region {sorted(seqregions)} != contigs GFF {sorted(contigs)}")
    for c, ln in seqregions.items():
        if contigs.get(c) != ln:
            errors.append(f"GFF sequence-region {c}={ln} != {contigs.get(c)}")

    ids = Counter()
    gff_keys = set()
    for ft in feats:
        where = f"GFF linea {ft['line']}"
        if ft["type"] != "polypeptide":
            errors.append(f"{where}: tipo {ft['type']} (se esperaba polypeptide)")
        if ft["source"] != "membrane-builder":
            errors.append(f"{where}: source {ft['source']}")
        if ft["strand"] not in {"+", "-"}:
            errors.append(f"{where}: hebra {ft['strand']}")
        if ft["phase"] != ".":
            errors.append(f"{where}: phase {ft['phase']} (se esperaba '.')")
        if ft["seqid"] not in seqregions:
            errors.append(f"{where}: seqid desconocido {ft['seqid']}")
        elif not (1 <= ft["start"] <= ft["end"] <= seqregions[ft["seqid"]]):
            errors.append(f"{where}: coordenadas fuera del contig")
        for k in ("ID", "protein_id", "membrane_type", "family", "substrate",
                  "evidence", "confidence"):
            if k not in ft["attrs"]:
                errors.append(f"{where}: falta atributo {k}")
        rid = ft["attrs"].get("ID", "")
        ids[rid] += 1
        if ft["attrs"].get("membrane_type") not in ALLOWED_TYPES:
            errors.append(f"{where}: membrane_type invalido")
        if ft["attrs"].get("confidence") not in ALLOWED_CONF:
            errors.append(f"{where}: confidence invalida")
        gff_keys.add((ft["attrs"].get("protein_id"), ft["attrs"].get("membrane_type"),
                      ft["start"], ft["end"], ft["strand"]))
    repeated = [k for k, v in ids.items() if v > 1]
    if repeated:
        errors.append(f"GFF IDs repetidos: {repeated[:5]}")

    tsv_keys = {(r["protein_id"], r["type"], int(r["start"]), int(r["end"]), r["strand"])
                for r in rows if r["start"].lstrip("-").isdigit() and r["end"].lstrip("-").isdigit()}
    if tsv_keys != gff_keys:
        only_t = sorted(tsv_keys - gff_keys)[:3]
        only_g = sorted(gff_keys - tsv_keys)[:3]
        errors.append(f"TSV/GFF no coinciden; solo TSV={only_t}; solo GFF={only_g}")

    if not rows:
        warnings.append("sin filas")
    return dict(ok=not errors, pending=False, rows=len(rows), gff_features=len(feats),
                errors=errors, warnings=warnings)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", default=str(ROOT / "data/membrane/pilot_genomes.tsv"))
    ap.add_argument("--outdir", default=str(ROOT / "data/membrane"))
    ap.add_argument("--tsv-dir", default=None,
                    help="directorio de TSV (por defecto --outdir)")
    ap.add_argument("--gff-dir", default=None,
                    help="directorio de GFF3 (por defecto <outdir>/gff)")
    ap.add_argument("--annot-root", default=str(ROOT / "data/annotations_ab_initio"))
    ap.add_argument("--json", default=str(ROOT / "data/membrane/validation.json"))
    ap.add_argument("--sample", type=int, default=0,
                    help="validar solo una muestra aleatoria de N genomas (0=todos)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--allow-missing", action="store_true",
                    help="tratar genomas sin salida como pendientes, no como error")
    args = ap.parse_args()

    annot_root = Path(args.annot_root)
    tsv_dir = Path(args.tsv_dir) if args.tsv_dir else Path(args.outdir)
    gff_dir = Path(args.gff_dir) if args.gff_dir else Path(args.outdir) / "gff"
    all_genomes = read_pilot(Path(args.pilot))
    if args.sample and args.sample < len(all_genomes):
        genomes = random.Random(args.seed).sample(all_genomes, args.sample)
        genomes.sort()
    else:
        genomes = all_genomes

    results = {}
    total_errors = 0
    for g in genomes:
        res = validate_genome(g, tsv_dir, gff_dir, annot_root,
                              allow_missing=args.allow_missing)
        results[g] = res
        total_errors += len(res["errors"])
        if not args.quiet:
            status = "OK" if res["ok"] else f"FALLA ({len(res['errors'])})"
            print(f"[validate] {g}: {status}, filas={res['rows']}, gff={res['gff_features']}")
            for e in res["errors"][:8]:
                print(f"    - {e}")
        elif not res["ok"]:
            print(f"[validate] FALLA {g}: {res['errors'][:3]}")

    pending_n = sum(1 for r in results.values() if r.get("pending"))
    payload = dict(
        ok=total_errors == 0 and pending_n == 0,
        processed_ok=total_errors == 0,
        checked_at=datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        tsv_dir=str(tsv_dir),
        gff_dir=str(gff_dir),
        sample=args.sample if args.sample else "all",
        genomes=results,
        totals=dict(genomes=len(genomes), expected=len(all_genomes),
                    rows=sum(r["rows"] for r in results.values()),
                    gff_features=sum(r["gff_features"] for r in results.values()),
                    errors=total_errors, pending=pending_n,
                    genomes_ok=sum(1 for r in results.values() if r.get("ok") and not r.get("pending"))),
    )
    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.json, "w") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print(f"[validate] {payload['totals']} -> {args.json}")
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
