#!/usr/bin/env python3
"""Deduplicación global de candidatos entre genomas.

Toma los `*_ghost_candidates.json` de N genomas, agrupa los péptidos por
identidad (95% por defecto, umbral elegido en `dedup_eval`) y produce:

  - un representante por cluster,
  - en cuántos genomas aparece,
  - cuántas copias hay (parálogos),
  - el mejor locus/rank de cada cluster.

Así las etapas caras (ESM-2, Pfam) se corren una vez por cluster y no una vez
por ORF.

Uso:
    python3 -m ghost.dedup --inputs "out/*_ghost_candidates.json" \
        --out out/global --identity 95
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Occurrence:
    genome: str
    contig: str
    start: int
    end: int
    strand: str
    rank: float
    aa: str
    aa_len: int

    @property
    def locus(self) -> str:
        return f"{self.contig}:{self.start}-{self.end}{self.strand}"


def load_candidates(paths: list[Path]) -> list[Occurrence]:
    out: list[Occurrence] = []
    for path in paths:
        try:
            doc = json.load(open(path))
        except Exception:
            continue
        for c in doc.get("candidates", []):
            aa = c.get("aa", "")
            if not aa:
                continue
            out.append(Occurrence(
                genome=c.get("genome", path.name.split("_ghost")[0]),
                contig=c["contig"], start=int(c["start"]), end=int(c["end"]),
                strand=c["strand"], rank=float(c.get("rank", 0.0)),
                aa=aa, aa_len=len(aa),
            ))
    return out


def write_unique_fasta(occs: list[Occurrence], path: Path) -> dict[str, str]:
    """Deduplica por secuencia exacta. Devuelve seq_id -> secuencia."""
    by_seq: dict[str, str] = {}
    for o in occs:
        if o.aa not in by_seq:
            by_seq[o.aa] = f"s{len(by_seq):07d}"
    with open(path, "w") as fh:
        for seq, sid in by_seq.items():
            fh.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return by_seq


def run_cluster(unique_faa: Path, out_tsv: Path, identity: float, threads: int, workdir: Path) -> None:
    db = workdir / "unique.dmnd"
    subprocess.run(
        ["diamond", "makedb", "--in", str(unique_faa), "-d", str(db.with_suffix(""))],
        check=True, capture_output=True,
    )
    subprocess.run(
        ["diamond", "cluster", "--db", str(db), "--out", str(out_tsv),
         "--approx-id", f"{identity:g}", "--threads", str(threads),
         "--memory-limit", "8G", "--quiet"],
        check=True,
    )


def parse_clusters(path: Path) -> dict[str, list[str]]:
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                clusters[p[0]].append(p[1])
    return clusters


def aggregate(
    clusters: dict[str, list[str]],
    occs: list[Occurrence],
    seq_to_id: dict[str, str],
) -> tuple[list[dict], list[dict]]:
    """Une clusters (ids de secuencia) con ocurrencias (loci)."""
    id_to_seq = {v: k for k, v in seq_to_id.items()}
    occ_by_seq: dict[str, list[Occurrence]] = defaultdict(list)
    for o in occs:
        occ_by_seq[o.aa].append(o)

    cluster_rows: list[dict] = []
    member_rows: list[dict] = []
    for rep_id, members in clusters.items():
        occurrences: list[Occurrence] = []
        for sid in members:
            occurrences.extend(occ_by_seq.get(id_to_seq.get(sid, ""), []))
        if not occurrences:
            continue
        genomes = sorted({o.genome for o in occurrences})
        best = max(occurrences, key=lambda o: (o.rank, o.aa_len))
        cid = f"c{len(cluster_rows):06d}"
        cluster_rows.append({
            "cluster": cid,
            "n_genomes": len(genomes),
            "n_copies": len(occurrences),
            "aa_len": len(id_to_seq.get(rep_id, "")),
            "best_rank": round(best.rank, 4),
            "best_locus": best.locus,
            "best_genome": best.genome,
            "rep_id": rep_id,
            "rep_aa": id_to_seq.get(rep_id, ""),
            "genomes": ",".join(genomes[:20]) + ("..." if len(genomes) > 20 else ""),
        })
        for o in occurrences:
            member_rows.append({
                "cluster": cid, "genome": o.genome, "locus": o.locus,
                "rank": round(o.rank, 4), "aa_len": o.aa_len,
            })
    cluster_rows.sort(key=lambda r: (r["n_genomes"], r["best_rank"]), reverse=True)
    return cluster_rows, member_rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="globs de *_ghost_candidates.json")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--identity", type=float, default=95.0)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    paths: list[Path] = []
    for pat in a.inputs:
        paths.extend(sorted(Path(p) for p in glob.glob(pat)))
    if not paths:
        print("ERROR: no se encontraron candidatos")
        return 2

    a.out.mkdir(parents=True, exist_ok=True)
    occs = load_candidates(paths)
    genomes = {o.genome for o in occs}
    print(f"candidatos: {len(occs):,} de {len(genomes)} genomas")

    unique_faa = a.out / "unique.faa"
    seq_to_id = write_unique_fasta(occs, unique_faa)
    print(f"secuencias únicas: {len(seq_to_id):,}")

    clusters_tsv = a.out / "clusters_raw.tsv"
    run_cluster(unique_faa, clusters_tsv, a.identity, a.threads, a.out)
    clusters = parse_clusters(clusters_tsv)
    print(f"clusters a {a.identity:g}%: {len(clusters):,}")

    cluster_rows, member_rows = aggregate(clusters, occs, seq_to_id)

    with open(a.out / "dedup_clusters.tsv", "w", newline="") as fh:
        cols = [k for k in cluster_rows[0] if k != "rep_aa"] if cluster_rows else ["cluster"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(cluster_rows)
    with open(a.out / "dedup_members.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(member_rows[0].keys()) if member_rows else ["cluster"], delimiter="\t")
        w.writeheader()
        w.writerows(member_rows)
    with open(a.out / "representatives.faa", "w") as fh:
        for r in cluster_rows:
            fh.write(f">{r['cluster']}|ng={r['n_genomes']}|rank={r['best_rank']}\n")
            aa = r["rep_aa"]
            for i in range(0, len(aa), 60):
                fh.write(aa[i:i + 60] + "\n")

    multi = sum(1 for r in cluster_rows if r["n_genomes"] >= 2)
    summary = {
        "genomes": len(genomes),
        "candidates": len(occs),
        "unique_sequences": len(seq_to_id),
        "clusters": len(cluster_rows),
        "multi_genome_clusters": multi,
        "identity": a.identity,
        "compression_vs_unique": round(1 - len(cluster_rows) / max(1, len(seq_to_id)), 4),
        "compression_vs_candidates": round(1 - len(cluster_rows) / max(1, len(occs)), 4),
    }
    (a.out / "dedup_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"  -> {a.out/'dedup_clusters.tsv'}")
    print(f"  -> {a.out/'representatives.faa'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
