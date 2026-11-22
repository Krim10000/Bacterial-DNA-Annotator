#!/usr/bin/env python3
"""01_extract_candidates.py — une la detección previa (JSON) con coordenadas reales.

Problema: `data/mobile_elements.json`, `data/defense_systems.json` y
`data/defense_virulence.json` guardan solo `{start, end, name}` y **no** traen
contig ni strand. Sus coordenadas provienen del GFF de NCBI RefSeq
(`data/annotations/<genoma>.gff.gz`), no del GFF ab initio.

Método (verificado, no inventado):
  1. Se indexan los CDS del GFF de NCBI por `(start, end)`; cada ítem previo
     coincide exactamente con un CDS de NCBI (100 % en los 10 pilotos, se
     comprueba y se registra).
  2. Ese CDS de NCBI se mapea al CDS del GFF ab initio por solapamiento
     recíproco ≥ 0,8 y misma hebra. Así la coordenada final vive en el sistema
     de coordenadas del proyecto (`data/annotations_ab_initio`), que es el que
     usan el resto de capas y el que tiene `proteins.faa` para HMMER.
  3. Si no hay CDS ab initio equivalente, se conserva la coordenada de NCBI y se
     marca `resolution=ncbi_only` (la coordenada sigue siendo válida y está
     dentro del contig). Si el ítem previo no se puede situar en ningún contig,
     se marca `unplaced` y se excluye de la capa tipada.

Salidas en `data/mobile_typing/_work/`:
  <genoma>/cds.tsv          todos los CDS ab initio + anotación NCBI mapeada
  <genoma>/candidates.tsv   ítems previos resueltos a contig/start/end/strand
  <genoma>/extract_stats.json
"""
from __future__ import annotations

import bisect
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths as P  # noqa: E402

ANN = ROOT / "data" / "annotations_ab_initio"
NCBI_GFF = ROOT / "data" / "annotations"
RAW = ROOT / "data" / "raw"
OUT = P.BASE
WORK = P.WORK

JSONS = {
    "mobile_elements": ROOT / "data" / "mobile_elements.json",
    "defense_systems": ROOT / "data" / "defense_systems.json",
    "defense_virulence": ROOT / "data" / "defense_virulence.json",
}
SECTIONS = {
    "mobile_elements": ["mobile"],
    "defense_systems": ["crispr", "prophage", "defense"],
    "defense_virulence": ["amr", "virulence"],
}


def parse_gff_attrs(text: str) -> dict[str, str]:
    attrs = {}
    for part in text.strip().rstrip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def sequence_regions(path: Path, gzipped: bool):
    opener = gzip.open if gzipped else open
    out = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                parts = line.split()
                out.append((parts[1], int(parts[3])))
            elif not line.startswith("#"):
                break
    return out


def parse_gff_cds(path: Path, gzipped: bool):
    """Devuelve (contigs, cds) del GFF; cds = dicts con attrs crudos."""
    opener = gzip.open if gzipped else open
    contigs = sequence_regions(path, gzipped)
    cds = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "region":
                # NCBI no siempre emite ##sequence-region para plásmidos: la
                # feature `region` sí trae el contig y su longitud.
                if parts[0] not in {c for c, _ in contigs}:
                    contigs.append((parts[0], int(parts[4])))
                continue
            if parts[2] != "CDS":
                continue
            attrs = parse_gff_attrs(parts[8])
            cds.append(
                {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "id": attrs.get("ID", ""),
                    "gene": attrs.get("gene", ""),
                    "locus_tag": attrs.get("locus_tag", ""),
                    "product": attrs.get("product", ""),
                    "protein_id": attrs.get("protein_id", ""),
                }
            )
    return contigs, cds


def reciprocal_overlap(a_start, a_end, b_start, b_end) -> float:
    overlap = min(a_end, b_end) - max(a_start, b_start) + 1
    if overlap <= 0:
        return 0.0
    return overlap / min(a_end - a_start + 1, b_end - b_start + 1)


class AbiIndex:
    """Índice de CDS ab initio por contig para mapear desde coordenadas NCBI."""

    def __init__(self, cds):
        self.by_contig = defaultdict(list)
        for rec in cds:
            self.by_contig[rec["contig"]].append(rec)
        for contig in self.by_contig:
            self.by_contig[contig].sort(key=lambda r: (r["start"], r["end"]))
        self.starts = {c: [r["start"] for r in recs] for c, recs in self.by_contig.items()}

    def best(self, contig, start, end, strand, min_frac=0.8):
        recs = self.by_contig.get(contig)
        if not recs:
            return None, 0.0
        starts = self.starts[contig]
        # ventana acotada: un CDS no puede solapar otro a >200 kb (el gen más
        # largo del corpus ronda 50 kb); evita el escaneo O(n²) del contig.
        lo = bisect.bisect_left(starts, start - 200_000)
        best, best_frac = None, 0.0
        for i in range(max(0, lo), len(recs)):
            rec = recs[i]
            if rec["start"] > end + 200_000:
                break
            if rec["strand"] != strand:
                continue
            frac = reciprocal_overlap(start, end, rec["start"], rec["end"])
            if frac > best_frac:
                best, best_frac = rec, frac
        if best is None or best_frac < min_frac:
            return None, best_frac
        return best, best_frac


def load_prior(genome: str):
    items = []
    cache = {}
    for name, path in JSONS.items():
        if name not in cache:
            cache[name] = json.loads(path.read_text())
        payload = cache[name].get(genome)
        if not payload:
            continue
        for section in SECTIONS[name]:
            if name == "mobile_elements":
                records = payload if isinstance(payload, list) else []
            else:
                records = payload.get(section, []) if isinstance(payload, dict) else []
            for idx, rec in enumerate(records):
                items.append(
                    {
                        "source_json": name,
                        "source_section": section,
                        "prior_index": idx,
                        "prior_name": rec.get("name", ""),
                        "start": int(rec["start"]),
                        "end": int(rec["end"]),
                    }
                )
    return items


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int, default=1, help="índice 1-based del shard")
    parser.add_argument("--shards", type=int, default=1, help="número total de shards")
    args, _ = parser.parse_known_args()
    pilots = P.genome_list()
    if args.shards > 1:
        pilots = pilots[args.shard - 1 :: args.shards]
    if not pilots:
        print(f"ERROR: {P.GENOMES} vacío", file=sys.stderr)
        return 1

    totals = defaultdict(int)
    for genome in pilots:
        print(f"[{genome}]")
        abi_dir = ANN / genome
        abi_contigs, abi_cds = parse_gff_cds(abi_dir / "annotation.gff3", gzipped=False)
        ncbi_contigs, ncbi_cds = parse_gff_cds(NCBI_GFF / f"{genome}.gff.gz", gzipped=True)
        lengths = dict(abi_contigs)
        ncbi_lengths = dict(ncbi_contigs)
        shared = set(lengths) & set(ncbi_lengths)
        bad_shared = [c for c in shared if lengths[c] != ncbi_lengths[c]]
        if bad_shared:
            print(f"  ERROR: longitudes distintas NCBI/ab initio en {genome}: {bad_shared}", file=sys.stderr)
        missing_ncbi = sorted(set(lengths) - set(ncbi_lengths))
        if missing_ncbi:
            print(f"  nota: {len(missing_ncbi)} contigs sin GFF NCBI (solo ab initio): {missing_ncbi[:4]}")

        # índice NCBI (start,end) -> CDS
        ncbi_by_coords = defaultdict(list)
        for rec in ncbi_cds:
            ncbi_by_coords[(rec["start"], rec["end"])].append(rec)

        abi_index = AbiIndex(abi_cds)
        abi_by_id = {rec["id"]: rec for rec in abi_cds}
        # anotación NCBI mapeada a cada CDS ab initio (la mejor por id)
        abi_annot = {}
        for rec in ncbi_cds:
            hit, frac = abi_index.best(rec["contig"], rec["start"], rec["end"], rec["strand"])
            if hit is None:
                continue
            prev = abi_annot.get(hit["id"])
            if prev is None or frac > prev["frac"]:
                abi_annot[hit["id"]] = {
                    "ncbi_id": rec["id"],
                    "gene": rec["gene"],
                    "locus_tag": rec["locus_tag"],
                    "product": rec["product"],
                    "protein_id": rec["protein_id"],
                    "ncbi_start": rec["start"],
                    "ncbi_end": rec["end"],
                    "frac": round(frac, 3),
                }

        # anotación funcional ab initio (Swiss-Prot) de respaldo
        func_path = abi_dir / "functional_annotation.json"
        func_hits = {}
        if func_path.exists():
            func_hits = json.loads(func_path.read_text()).get("hits", {})

        # cds.tsv: todos los CDS ab initio con la mejor anotación disponible
        work = WORK / genome
        work.mkdir(parents=True, exist_ok=True)
        with (work / "cds.tsv").open("w", newline="") as fh:
            fields = [
                "cds_id", "contig", "start", "end", "strand", "gene", "locus_tag",
                "product", "product_source", "protein_id", "ncbi_cds_id",
                "ncbi_start", "ncbi_end", "map_frac",
            ]
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for rec in abi_cds:
                ann = abi_annot.get(rec["id"])
                fh_it = func_hits.get(rec["id"], {})
                if ann and ann["product"]:
                    gene = ann["gene"] or fh_it.get("gene", "") or rec["gene"]
                    product = ann["product"]
                    source = "ncbi"
                elif fh_it:
                    gene = fh_it.get("gene", "") or rec["gene"]
                    product = fh_it.get("product", "") or rec["product"]
                    source = "ab_initio_swissprot"
                else:
                    gene = rec["gene"]
                    product = rec["product"]
                    source = "ab_initio"
                writer.writerow(
                    {
                        "cds_id": rec["id"],
                        "contig": rec["contig"],
                        "start": rec["start"],
                        "end": rec["end"],
                        "strand": rec["strand"],
                        "gene": gene,
                        "locus_tag": ann["locus_tag"] if ann else rec["locus_tag"],
                        "product": product,
                        "product_source": source,
                        "protein_id": (ann["protein_id"] if ann and ann["protein_id"] else rec["protein_id"]),
                        "ncbi_cds_id": ann["ncbi_id"] if ann else "",
                        "ncbi_start": ann["ncbi_start"] if ann else "",
                        "ncbi_end": ann["ncbi_end"] if ann else "",
                        "map_frac": ann["frac"] if ann else "",
                    }
                )

        with (work / "contigs.tsv").open("w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["order", "contig", "length"])
            for order, (contig, length) in enumerate(abi_contigs, start=1):
                writer.writerow([order, contig, length])

        # candidatos previos
        items = load_prior(genome)
        stats = defaultdict(int)
        cand_rows = []
        for item in items:
            stats["prior_total"] += 1
            matches = ncbi_by_coords.get((item["start"], item["end"]), [])
            stats["exact_ncbi"] += 1 if matches else 0
            chosen = None
            ambiguous = False
            if len(matches) == 1:
                chosen = matches[0]
            elif len(matches) > 1:
                named = [m for m in matches if m["product"].lower() == item["prior_name"].lower()]
                if len(named) == 1:
                    chosen = named[0]
                else:
                    ambiguous = True
                    chosen = matches[0]
            if chosen is None:
                # sin CDS exacto: buscar contención del nombre en CDS que solapen
                cands = []
                for rec in ncbi_cds:
                    if rec["start"] <= item["start"] <= rec["end"] and item["prior_name"].lower() == rec["product"].lower():
                        cands.append(rec)
                if len(cands) == 1:
                    chosen = cands[0]
                    stats["containment"] += 1
            if chosen is None:
                stats["unplaced"] += 1
                cand_rows.append(
                    {
                        "genome": genome, "source_json": item["source_json"],
                        "source_section": item["source_section"], "prior_index": item["prior_index"],
                        "prior_name": item["prior_name"], "resolution": "unplaced",
                        "contig": "unplaced", "start": item["start"], "end": item["end"],
                        "strand": ".", "cds_id": "", "gene": "", "locus_tag": "",
                        "product": item["prior_name"], "ncbi_cds_id": "", "map_frac": "",
                        "note": "no_cds_at_coords",
                    }
                )
                continue

            if chosen["contig"] not in lengths:
                stats["contig_unknown"] += 1
                cand_rows.append(
                    {
                        "genome": genome, "source_json": item["source_json"],
                        "source_section": item["source_section"], "prior_index": item["prior_index"],
                        "prior_name": item["prior_name"], "resolution": "unplaced",
                        "contig": "unplaced", "start": item["start"], "end": item["end"],
                        "strand": ".", "cds_id": "", "gene": "", "locus_tag": "",
                        "product": item["prior_name"], "ncbi_cds_id": chosen["id"], "map_frac": "",
                        "note": "contig_not_in_ab_initio",
                    }
                )
                continue
            if not (1 <= chosen["start"] <= chosen["end"] <= lengths[chosen["contig"]]):
                stats["out_of_contig"] += 1
                cand_rows.append(
                    {
                        "genome": genome, "source_json": item["source_json"],
                        "source_section": item["source_section"], "prior_index": item["prior_index"],
                        "prior_name": item["prior_name"], "resolution": "unplaced",
                        "contig": "unplaced", "start": item["start"], "end": item["end"],
                        "strand": ".", "cds_id": "", "gene": "", "locus_tag": "",
                        "product": item["prior_name"], "ncbi_cds_id": chosen["id"], "map_frac": "",
                        "note": "coords_outside_contig",
                    }
                )
                continue

            abi_hit, frac = abi_index.best(chosen["contig"], chosen["start"], chosen["end"], chosen["strand"])
            if abi_hit is not None and abi_hit["id"] in abi_by_id:
                stats["mapped_abi"] += 1
                resolution = "abi_exact" if frac >= 0.999 else "abi_overlap"
                contig = abi_hit["contig"]
                start, end, strand = abi_hit["start"], abi_hit["end"], abi_hit["strand"]
                cds_id = abi_hit["id"]
            else:
                stats["ncbi_only"] += 1
                resolution = "ncbi_only"
                contig = chosen["contig"]
                start, end, strand = chosen["start"], chosen["end"], chosen["strand"]
                cds_id = ""
            if ambiguous:
                stats["ambiguous_coords"] += 1
            cand_rows.append(
                {
                    "genome": genome, "source_json": item["source_json"],
                    "source_section": item["source_section"], "prior_index": item["prior_index"],
                    "prior_name": item["prior_name"], "resolution": resolution,
                    "contig": contig, "start": start, "end": end, "strand": strand,
                    "cds_id": cds_id, "gene": chosen["gene"], "locus_tag": chosen["locus_tag"],
                    "product": chosen["product"], "ncbi_cds_id": chosen["id"],
                    "map_frac": round(frac, 3) if abi_hit is not None else "",
                    "note": "ambiguous_coords" if ambiguous else "",
                }
            )

        fields = [
            "genome", "source_json", "source_section", "prior_index", "prior_name",
            "resolution", "contig", "start", "end", "strand", "cds_id", "gene",
            "locus_tag", "product", "ncbi_cds_id", "map_frac", "note",
        ]
        with (work / "candidates.tsv").open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(cand_rows)

        stats["cds_abi"] = len(abi_cds)
        stats["cds_ncbi"] = len(ncbi_cds)
        stats["cds_with_ncbi_product"] = len(abi_annot)
        stats["contigs"] = len(abi_contigs)
        (work / "extract_stats.json").write_text(json.dumps(dict(stats), indent=2, sort_keys=True) + "\n")
        for key, value in stats.items():
            totals[key] += value if key not in ("cds_abi", "cds_ncbi", "contigs", "cds_with_ncbi_product") else 0
        print(
            f"  CDS abi={stats['cds_abi']} ncbi={stats['cds_ncbi']} | previos={stats['prior_total']} "
            f"exact_ncbi={stats['exact_ncbi']} abi={stats['mapped_abi']} ncbi_only={stats['ncbi_only']} "
            f"unplaced={stats['unplaced']}"
        )

    print("\nTOTAL:", dict(totals))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
