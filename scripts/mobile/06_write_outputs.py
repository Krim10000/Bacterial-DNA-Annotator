#!/usr/bin/env python3
"""06_write_outputs.py — escribe TSV, GFF3, espaciadores CRISPR y summary.json.

Contrato de columnas del TSV (BACKLOG.md):
    genome, contig, start, end, strand, element_type, subtype, gene, evidence, confidence

Convenciones:
  - start/end son 1-based inclusivos, coordenadas del GFF ab initio.
  - strand es + o -; los arrays CRISPR (no codificantes) usan "." salvo que
    tengan un gen cas a ≤10 kb, en cuyo caso heredan su hebra.
  - El GFF3 usa tipos SO estándar: `repeat_region` (CRISPR), `CDS` (cas),
    `mobile_genetic_element` (IS/transposón/profago), `gene` (defensa/AMR/virulencia).
"""
from __future__ import annotations

import sys

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths as P  # noqa: E402

OUT = P.BASE
WORK = P.WORK

COLUMNS = ["genome", "contig", "start", "end", "strand", "element_type",
           "subtype", "gene", "evidence", "confidence"]
GFF_TYPE = {
    "crispr": "repeat_region",
    "cas": "CDS",
    "transposon": "mobile_genetic_element",
    "prophage": "mobile_genetic_element",
    "defense": "gene",
    "amr": "gene",
    "virulence": "gene",
}
VALID_TYPES = set(GFF_TYPE)


def esc(value: str) -> str:
    return quote(str(value), safe="")


def load_tsv(path: Path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def contig_table(genome):
    rows = load_tsv(WORK / genome / "contigs.tsv")
    return [(r["contig"], int(r["length"])) for r in rows]


def mobile_element_type(row) -> str:
    if row["element_type"] == "transposon":
        if row["subtype"].startswith("IS"):
            return f"insertion_sequence:{row['subtype']}"
        if row["subtype"].startswith("Tn"):
            return f"transposon:{row['subtype']}"
        if "conjugative" in row["subtype"]:
            return f"conjugative_transposon:{row['subtype']}"
        return f"transposon:{row['subtype']}"
    if row["element_type"] == "prophage":
        return "prophage"
    return row["element_type"]


def gff_attributes(row, feature_id: str) -> str:
    attrs = [f"ID={feature_id}", f"Name={esc(row['gene'])}"]
    if row["gene"] and row["gene"] != "region":
        attrs.append(f"gene={esc(row['gene'])}")
    if row["element_type"] == "transposon":
        attrs.append(f"mobile_element_type={esc(mobile_element_type(row))}")
    if row["element_type"] == "prophage":
        attrs.append(f"mobile_element_type={esc('prophage:' + row['subtype'])}")
    if row["element_type"] == "crispr":
        attrs.append("rpt_type=CRISPR")
        attrs.append("rpt_family=CRISPR")
    attrs.append(f"Note={esc(row['evidence'])}")
    attrs.append(f"confidence={esc(row['confidence'])}")
    return ";".join(attrs)


def main() -> int:
    pilots = load_tsv(P.GENOMES)
    summary = {
        "module": "mobile_typing",
        "pilot_genomes": [p["genome"] for p in pilots],
        "columns": COLUMNS,
        "coordinate_system": "GFF ab initio (data/annotations_ab_initio/<genome>/annotation.gff3), 1-based inclusive",
    }
    counts_by_genome, counts_by_type, counts_by_conf, counts_by_subtype = {}, Counter(), Counter(), Counter()
    outputs_hash = {}
    all_rows = []
    arrays_total = 0
    spacers_total = 0
    failures = {}
    processed = []
    coverage_rows = {}
    meta_rows = {}
    if P.FULL:
        for row in load_tsv(P.GENOMES):
            meta_rows[row["genome"]] = row

    for pilot in pilots:
        genome = pilot["genome"]
        elements_path = WORK / genome / "elements.tsv"
        if not elements_path.exists():
            failures[genome] = "elements.tsv ausente (01/05 no completaron)"
            continue
        contigs = contig_table(genome)
        order = {contig: i for i, (contig, _) in enumerate(contigs)}
        lengths = dict(contigs)
        rows = load_tsv(WORK / genome / "elements.tsv")
        for row in rows:
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
        rows.sort(key=lambda r: (order.get(r["contig"], 9999), r["start"], r["end"],
                                 r["element_type"], r["subtype"]))
        # TSV
        tsv_path = P.TSV_DIR / f"{genome}.tsv"
        with tsv_path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                                    extrasaction="ignore", quoting=csv.QUOTE_MINIMAL)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        # GFF3
        gff_path = P.GFF_DIR / f"{genome}.gff3"
        gff_path.parent.mkdir(parents=True, exist_ok=True)
        with gff_path.open("w") as fh:
            fh.write("##gff-version 3\n")
            fh.write(f"##genome-build {genome}\n")
            for contig, length in contigs:
                fh.write(f"##sequence-region {contig} 1 {length}\n")
            for idx, row in enumerate(rows, start=1):
                ftype = GFF_TYPE[row["element_type"]]
                feature_id = f"{genome}_{row['element_type']}_{idx:04d}"
                strand = row["strand"] if row["strand"] in ("+", "-", ".") else "."
                phase = "0" if row["element_type"] == "cas" else "."
                fh.write("\t".join([
                    row["contig"], "mobile_typing", ftype,
                    str(row["start"]), str(row["end"]), ".", strand, phase,
                    gff_attributes(row, feature_id),
                ]) + "\n")
        # espaciadores CRISPR
        sp_path = P.SPACER_DIR / f"{genome}.tsv"
        sp_path.parent.mkdir(parents=True, exist_ok=True)
        sp_src = WORK / genome / "crispr_spacers.tsv"
        if sp_src.exists():
            # solo espaciadores de arrays que están en elements.tsv (evita
            # carreras con el pool CRISPR que puede añadir arrays después)
            array_ids = {r["gene"] for r in rows if r["element_type"] == "crispr"}
            sp_rows = [r for r in load_tsv(sp_src) if r["array_id"] in array_ids]
        else:
            sp_rows = []
        with sp_path.open("w", newline="") as fh:
            fields = ["genome", "array_id", "contig", "spacer_index", "start", "end", "length", "sequence"]
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in sp_rows:
                writer.writerow({k: row[k] for k in fields})
        arrays_total += sum(1 for r in rows if r["element_type"] == "crispr")
        spacers_total += len(sp_rows)

        # conteos
        by_type = Counter(r["element_type"] for r in rows)
        by_conf = Counter(r["confidence"] for r in rows)
        by_subtype = Counter(f"{r['element_type']}:{r['subtype']}" for r in rows)
        counts_by_genome[genome] = {
            "rows": len(rows),
            "by_type": dict(sorted(by_type.items())),
            "by_confidence": dict(sorted(by_conf.items())),
            "crispr_spacers": len(sp_rows),
        }
        counts_by_type.update(by_type)
        counts_by_conf.update(by_conf)
        counts_by_subtype.update(by_subtype)
        outputs_hash[genome] = {
            "tsv_sha256": hashlib.sha256(tsv_path.read_bytes()).hexdigest(),
            "gff_sha256": hashlib.sha256(gff_path.read_bytes()).hexdigest(),
        }
        all_rows.extend(rows)
        processed.append(genome)
        coverage_rows[genome] = {
            "pfam_scanned": (WORK / genome / "pfam_hits.tsv").exists(),
            "crispr_scanned": (WORK / genome / "crispr_stats.json").exists(),
            "has_fasta": (meta_rows.get(genome, {}).get("has_fasta", "yes") == "yes"),
        }
        print(f"[{genome}] {len(rows)} filas, {len(sp_rows)} espaciadores")

    # ---- summary -----------------------------------------------------------
    demoted_all = []
    prior_stats = Counter()
    for pilot in pilots:
        genome = pilot["genome"]
        for row in load_tsv(WORK / genome / "demoted.tsv"):
            demoted_all.append(row)
        stats_path = WORK / genome / "extract_stats.json"
        if stats_path.exists():
            prior_stats.update(json.loads(stats_path.read_text()))
    reasons = Counter(r["reason"] for r in demoted_all)
    pfam_stats = Counter()
    pfam_models = 0
    full_panel_genomes = 0
    core_info = OUT / "_work" / "pfam_core" / "panel_info.json"
    if core_info.exists():
        pfam_models = json.loads(core_info.read_text()).get("models", 0)
    for pilot in pilots:
        path = WORK / pilot["genome"] / "pfam_scan_stats.json"
        if path.exists():
            stats = json.loads(path.read_text())
            for key, value in stats.items():
                if key == "models":
                    if value != pfam_models:
                        full_panel_genomes += 1
                else:
                    pfam_stats[key] += value
    prior_by_section = Counter()
    demoted_by_section = Counter()
    for pilot in pilots:
        genome = pilot["genome"]
        for row in load_tsv(WORK / genome / "candidates.tsv"):
            prior_by_section[row["source_section"]] += 1
        for row in load_tsv(WORK / genome / "demoted.tsv"):
            demoted_by_section[row["source_section"]] += 1

    if P.FULL:
        classify_failures = {}
        for cf in sorted((WORK / "state").glob("classify_failures*.json")):
            classify_failures.update(json.loads(cf.read_text()))
        failures.update({g: v for g, v in classify_failures.items() if g not in failures})
        summary["coverage"] = {
            "genomes_expected": len(pilots),
            "genomes_processed": len(processed),
            "genomes_failed": len(failures),
            "genomes_with_fasta": sum(1 for v in coverage_rows.values() if v["has_fasta"]),
            "genomes_without_fasta": sum(1 for v in coverage_rows.values() if not v["has_fasta"]),
            "pfam_scanned": sum(1 for v in coverage_rows.values() if v["pfam_scanned"]),
            "crispr_scanned": sum(1 for v in coverage_rows.values() if v["crispr_scanned"]),
            "failures": failures,
            "without_pfam": [g for g, v in coverage_rows.items() if not v["pfam_scanned"]],
            "without_crispr": [g for g, v in coverage_rows.items() if not v["crispr_scanned"]],
        }
    summary.update({
        "counts": {
            "rows_total": sum(counts_by_genome[g]["rows"] for g in counts_by_genome),
            "by_type": dict(sorted(counts_by_type.items())),
            "by_confidence": dict(sorted(counts_by_conf.items())),
            "by_genome": counts_by_genome,
            "by_subtype": dict(sorted(counts_by_subtype.items())),
            "crispr_arrays": arrays_total,
            "crispr_spacers": spacers_total,
        },
        "evidence": {
            "prior_items_total": prior_stats.get("prior_total", 0),
            "prior_exact_ncbi_cds": prior_stats.get("exact_ncbi", 0),
            "prior_mapped_to_ab_initio": prior_stats.get("mapped_abi", 0),
            "prior_ncbi_only": prior_stats.get("ncbi_only", 0),
            "prior_unplaced": prior_stats.get("unplaced", 0),
            "prior_out_of_contig": prior_stats.get("out_of_contig", 0),
            "prior_by_section": dict(sorted(prior_by_section.items())),
            "prior_demoted_or_reclassified": len(demoted_all),
            "prior_demoted_by_section": dict(sorted(demoted_by_section.items())),
            "demotion_reasons": dict(sorted(reasons.items())),
            "pfam": {
                "panel_models_core": pfam_models,
                "genomes_with_full_339_panel": full_panel_genomes,
                "protein_sequences_scanned": pfam_stats.get("proteins", 0),
                "hits": pfam_stats.get("hits", 0),
                "hits_passing_ga": pfam_stats.get("hits_passing_ga", 0),
            },
        },
        "outputs_sha256": outputs_hash,
    })
    (P.SUMMARY).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    if P.FULL:
        print("Cobertura:", {k: summary["coverage"][k] for k in
                             ("genomes_processed", "genomes_failed", "pfam_scanned",
                              "crispr_scanned", "genomes_without_fasta")})
    print(f"\nTotal: {summary['counts']['rows_total']} filas, "
          f"{summary['counts']['crispr_arrays']} arrays, "
          f"{summary['counts']['crispr_spacers']} espaciadores")
    print("Por tipo:", summary["counts"]["by_type"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
