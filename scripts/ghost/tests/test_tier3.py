"""Pruebas de la extracción Tier 3 (modelos alternativos)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.tier3 import (  # noqa: E402
    alt_start_within,
    classify,
    gene_protein,
    is_imprecise,
)


def gene(**kw) -> dict:
    base = {
        "contig": "c", "start": 1, "end": 100, "strand": "+", "id": "g1",
        "conf": 100.0, "partial": "00", "start_type": "ATG",
    }
    base.update(kw)
    return base


class TestTier3(unittest.TestCase):
    def test_gene_protein(self):
        seq = "ATG" + "GCT" * 10 + "TAA"
        self.assertEqual(gene_protein(seq, gene(start=1, end=len(seq))), "M" + "A" * 10)

    def test_alt_start_within(self):
        seq = "ATG" + "GTG" + "GCT" * 20
        self.assertTrue(alt_start_within(seq, 1, 60, "+"))
        seq2 = "ATG" + "GCT" * 20
        self.assertFalse(alt_start_within(seq2, 1, 60, "+"))

    def test_is_imprecise(self):
        clean = "ATG" + "GCT" * 100
        self.assertFalse(is_imprecise(gene(end=303), clean))
        self.assertTrue(is_imprecise(gene(end=303, conf=79), clean))
        self.assertTrue(is_imprecise(gene(end=303, partial="01"), clean))
        self.assertTrue(is_imprecise(gene(end=120), clean[:120]))
        self.assertTrue(is_imprecise(gene(end=303, start_type="ATT"), clean))
        self.assertTrue(is_imprecise(gene(end=303), "ATG" + "GTG" + "GCT" * 99))

    def test_classify(self):
        g = gene(start=100, end=400, strand="+")
        self.assertEqual(classify(g, {"start": 50, "end": 400, "strand": "+"}), "extended_5p")
        self.assertEqual(classify(g, {"start": 130, "end": 400, "strand": "+"}), "truncated_5p")
        self.assertEqual(classify(g, {"start": 101, "end": 400, "strand": "+"}), "frameshift")
        self.assertEqual(classify(g, {"start": 50, "end": 400, "strand": "-"}), "antisense")
        gm = gene(start=100, end=400, strand="-")
        self.assertEqual(classify(gm, {"start": 100, "end": 450, "strand": "-"}), "extended_5p")


if __name__ == "__main__":
    unittest.main()
