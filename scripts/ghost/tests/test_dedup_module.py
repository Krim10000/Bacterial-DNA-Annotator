"""Pruebas del módulo de deduplicación global."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.dedup import (  # noqa: E402
    Occurrence,
    aggregate,
    load_candidates,
    write_unique_fasta,
)


def write_json(path: Path, genome: str, aas: list[str]) -> None:
    cands = [{
        "genome": genome, "contig": "chr", "start": i * 100 + 1,
        "end": i * 100 + 90, "strand": "+", "aa": aa, "rank": 0.5 + i * 0.01,
    } for i, aa in enumerate(aas)]
    path.write_text(json.dumps({"stats": {}, "candidates": cands}))


class TestDedupModule(unittest.TestCase):
    def test_load_candidates(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "GCF_A_ghost_candidates.json"
            write_json(p, "GCF_A", ["MKAVL", "MKAVL", "GGGG"])
            occs = load_candidates([p])
            self.assertEqual(len(occs), 3)
            self.assertEqual(occs[0].genome, "GCF_A")
            self.assertEqual(occs[0].locus, "chr:1-90+")

    def test_unique_fasta(self):
        occs = [
            Occurrence("g1", "c", 1, 90, "+", 0.5, "MKAVL", 5),
            Occurrence("g2", "c", 1, 90, "+", 0.6, "MKAVL", 5),
            Occurrence("g1", "c", 200, 290, "+", 0.4, "GGGG", 4),
        ]
        with tempfile.TemporaryDirectory() as td:
            fasta = Path(td) / "u.faa"
            mapping = write_unique_fasta(occs, fasta)
            self.assertEqual(len(mapping), 2)
            self.assertEqual(fasta.read_text().count(">"), 2)

    def test_aggregate_counts_genomes_and_copies(self):
        occs = [
            Occurrence("g1", "c", 1, 90, "+", 0.5, "MKAVL", 5),
            Occurrence("g2", "c", 1, 90, "+", 0.7, "MKAVL", 5),
            Occurrence("g1", "c", 300, 390, "+", 0.4, "GGGG", 4),
        ]
        seq_to_id = {"MKAVL": "s0000000", "GGGG": "s0000001"}
        clusters = {"s0000000": ["s0000000"], "s0000001": ["s0000001"]}
        rows, members = aggregate(clusters, occs, seq_to_id)
        top = [r for r in rows if r["n_genomes"] == 2]
        self.assertEqual(len(top), 1)
        self.assertEqual(top[0]["n_copies"], 2)
        self.assertEqual(top[0]["best_rank"], 0.7)
        self.assertEqual(len(members), 3)

    def test_aggregate_same_genome_copies(self):
        occs = [
            Occurrence("g1", "c", 1, 90, "+", 0.5, "MKAVL", 5),
            Occurrence("g1", "c", 200, 290, "+", 0.6, "MKAVL", 5),
        ]
        clusters = {"s0000000": ["s0000000"]}
        rows, _ = aggregate(clusters, occs, {"MKAVL": "s0000000"})
        self.assertEqual(rows[0]["n_genomes"], 1)
        self.assertEqual(rows[0]["n_copies"], 2)

    def test_empty_cluster_ignored(self):
        occs = [Occurrence("g1", "c", 1, 90, "+", 0.5, "MKAVL", 5)]
        rows, members = aggregate({"s9": ["s9"]}, occs, {"MKAVL": "s0000000"})
        self.assertEqual(rows, [])
        self.assertEqual(members, [])


if __name__ == "__main__":
    unittest.main()
