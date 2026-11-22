#!/usr/bin/env python3
"""07_validate.py — valida TSV, GFF3 y espaciadores antes de integrar.

Comprueba:
  - TSV: cabecera exacta, tipos, coordenadas 1-based dentro del contig, strand válido.
  - GFF3: 9 columnas, `##sequence-region` por contig, coordenadas dentro del
    contig, ID único, atributos con ID, tipos SO permitidos.
  - Correspondencia TSV ↔ GFF (mismo número de filas y mismas coordenadas).
  - Espaciadores CRISPR dentro de su array y del contig, longitud = end-start+1.
  - Hash de determinismo guardado en summary.json == ficheros actuales.

Salida: data/mobile_typing/validation.json; exit code 1 si hay errores.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths as P  # noqa: E402

OUT = P.BASE
WORK = P.WORK

COLUMNS = ["genome", "contig", "start", "end", "strand", "element_type",
           "subtype", "gene", "evidence", "confidence"]
ALLOWED_TYPES = {"crispr", "cas", "transposon", "prophage", "defense", "amr", "virulence"}
ALLOWED_GFF = {"repeat_region", "CDS", "mobile_genetic_element", "gene"}
ALLOWED_CONF = {"high", "medium", "low"}


def read_tsv(path: Path):
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames, list(reader)


def parse_gff(path: Path):
    features, regions = [], []
    for line in path.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            regions.append((parts[1], int(parts[3])))
        elif line.startswith("#") or not line:
            continue
        else:
            features.append(line.split("\t"))
    return regions, features


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", type=int, default=0,
                        help="valida solo cada N-ésimo genoma (0 = todos)")
    args, _ = parser.parse_known_args()
    pilots = read_tsv(P.GENOMES)[1]
    if args.sample and args.sample > 1:
        pilots = pilots[:: args.sample]
    summary = json.loads(P.SUMMARY.read_text())
    errors, warnings, gaps = [], [], []
    stats = {"genomes": 0, "tsv_rows": 0, "gff_features": 0, "spacers": 0}
    per_genome = {}

    for pilot in pilots:
        genome = pilot["genome"]
        stats["genomes"] += 1
        contigs = {r["contig"]: int(r["length"]) for r in read_tsv(WORK / genome / "contigs.tsv")[1]}
        tsv_path = P.TSV_DIR / f"{genome}.tsv"
        gff_path = P.GFF_DIR / f"{genome}.gff3"
        if not tsv_path.exists() or not gff_path.exists():
            gaps.append(genome)
            continue
        if not (WORK / genome / "contigs.tsv").exists():
            gaps.append(genome)
            continue
        header, rows = read_tsv(tsv_path)
        if header != COLUMNS:
            errors.append(f"{genome}: cabecera TSV {header} != {COLUMNS}")
        gff_regions, gff_features = parse_gff(gff_path)
        if [c for c, _ in gff_regions] != list(contigs):
            errors.append(f"{genome}: ##sequence-region {gff_regions} != contigs {list(contigs)}")

        gff_coords, gff_ids = [], []
        for f in gff_features:
            if len(f) != 9:
                errors.append(f"{genome}: línea GFF con {len(f)} columnas")
                continue
            seqid, _, ftype, start, end, score, strand, phase, attrs = f
            if seqid not in contigs:
                errors.append(f"{genome}: GFF seqid {seqid} desconocido")
                continue
            if ftype not in ALLOWED_GFF:
                errors.append(f"{genome}: tipo GFF no permitido {ftype}")
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                errors.append(f"{genome}: coordenadas GFF no numéricas {start}-{end}")
                continue
            if not (1 <= start_i <= end_i <= contigs[seqid]):
                errors.append(f"{genome}: GFF fuera de contig {seqid}:{start_i}-{end_i}")
            if strand not in ("+", "-", ".", "?"):
                errors.append(f"{genome}: strand GFF inválido {strand}")
            if "ID=" not in attrs:
                errors.append(f"{genome}: feature GFF sin ID: {attrs[:60]}")
            else:
                gff_ids.append(attrs.split("ID=", 1)[1].split(";", 1)[0])
            gff_coords.append((seqid, start_i, end_i, ftype))
            stats["gff_features"] += 1
        if len(gff_ids) != len(set(gff_ids)):
            errors.append(f"{genome}: IDs GFF duplicados")
        if len(gff_features) != len(rows):
            errors.append(f"{genome}: {len(gff_features)} features GFF != {len(rows)} filas TSV")

        tsv_coords = []
        by_type = Counter()
        for row in rows:
            stats["tsv_rows"] += 1
            if row["genome"] != genome:
                errors.append(f"{genome}: genome de fila {row['genome']}")
            if row["element_type"] not in ALLOWED_TYPES:
                errors.append(f"{genome}: element_type inválido {row['element_type']}")
            if row["confidence"] not in ALLOWED_CONF:
                errors.append(f"{genome}: confidence inválida {row['confidence']}")
            contig = row["contig"]
            if contig not in contigs:
                errors.append(f"{genome}: contig TSV {contig} desconocido")
                continue
            try:
                start, end = int(row["start"]), int(row["end"])
            except ValueError:
                errors.append(f"{genome}: start/end no numéricos {row['start']}-{row['end']}")
                continue
            if not (1 <= start <= end <= contigs[contig]):
                errors.append(f"{genome}: TSV fuera de contig {contig}:{start}-{end}")
            strand = row["strand"]
            if strand not in ("+", "-") and not (strand == "." and row["element_type"] == "crispr"):
                errors.append(f"{genome}: strand TSV inválido {strand} en {row['element_type']}")
            if not row["gene"]:
                errors.append(f"{genome}: gene vacío en {contig}:{start}")
            if not row["evidence"]:
                errors.append(f"{genome}: evidence vacía en {contig}:{start}")
            by_type[row["element_type"]] += 1
            tsv_coords.append((contig, start, end))
        dup_keys = Counter((r["contig"], int(r["start"]), int(r["end"]), r["element_type"]) for r in rows)
        dups = [k for k, v in dup_keys.items() if v > 1]
        if dups:
            errors.append(f"{genome}: {len(dups)} filas duplicadas, p.ej. {dups[:3]}")
        gff_plain = [(c, s, e) for c, s, e, _ in gff_coords]
        if Counter(tsv_coords) != Counter(gff_plain):
            errors.append(f"{genome}: coordenadas TSV y GFF no coinciden")

        # espaciadores
        sp_path = P.SPACER_DIR / f"{genome}.tsv"
        n_sp = 0
        if sp_path.exists():
            sp_header, sp_rows = read_tsv(sp_path)
            expected_sp = ["genome", "array_id", "contig", "spacer_index", "start", "end", "length", "sequence"]
            if sp_header != expected_sp:
                errors.append(f"{genome}: cabecera de espaciadores inesperada")
            arrays = {r["gene"]: r for r in rows if r["element_type"] == "crispr"}
            for sp in sp_rows:
                n_sp += 1
                stats["spacers"] += 1
                arr = arrays.get(sp["array_id"])
                if not arr:
                    errors.append(f"{genome}: espaciador de array desconocido {sp['array_id']}")
                    continue
                c = sp["contig"]
                s, e = int(sp["start"]), int(sp["end"])
                if c not in contigs or not (1 <= s <= e <= contigs[c]):
                    errors.append(f"{genome}: espaciador fuera de contig {c}:{s}-{e}")
                if not (int(arr["start"]) <= s <= e <= int(arr["end"])):
                    errors.append(f"{genome}: espaciador fuera de su array {sp['array_id']}")
                if int(sp["length"]) != e - s + 1 or len(sp["sequence"]) != int(sp["length"]):
                    errors.append(f"{genome}: longitud de espaciador inconsistente {sp['array_id']}:{sp['spacer_index']}")
        stats["spacers"] = stats["spacers"]
        per_genome[genome] = {
            "tsv_rows": len(rows), "gff_features": len(gff_features),
            "spacers": n_sp, "by_type": dict(sorted(by_type.items())),
        }

    # determinismo: hashes guardados
    for genome, hashes in summary.get("outputs_sha256", {}).items():
        for kind, stored in hashes.items():
            path = P.TSV_DIR / f"{genome}.tsv" if kind.startswith("tsv") else P.GFF_DIR / f"{genome}.gff3"
            if not path.exists():
                errors.append(f"{genome}: falta {path.name} para verificar hash")
                continue
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != stored:
                errors.append(f"{genome}: hash {kind} no coincide (¿salida modificada a mano?)")

    status = "fail" if errors else ("pass_with_gaps" if gaps else "pass")
    result = {
        "status": status,
        "errors": errors,
        "warnings": warnings,
        "gaps": gaps,
        "coverage": {
            "genomes_expected": len(pilots),
            "genomes_validated": stats["genomes"],
            "genomes_without_output": len(gaps),
        },
        "stats": stats,
        "per_genome": per_genome,
        "checks": {
            "columns": COLUMNS,
            "coords_within_contig": not any("fuera de contig" in e for e in errors),
            "gff_parsed": not any("GFF" in e for e in errors),
            "tsv_gff_consistent": not any("no coinciden" in e for e in errors),
            "deterministic_hashes": not any("hash" in e for e in errors),
        },
    }
    (P.VALIDATION).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if errors:
        print(f"FALLO: {len(errors)} errores", file=sys.stderr)
        for error in errors[:30]:
            print(f"  - {error}", file=sys.stderr)
        return 1
    print(f"{'OK' if not gaps else 'OK CON HUECOS'}: {stats['genomes']} genomas, "
          f"{stats['tsv_rows']} filas TSV, {stats['gff_features']} features GFF, "
          f"{stats['spacers']} espaciadores, {len(gaps)} sin salida")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
