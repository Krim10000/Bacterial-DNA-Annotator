"""Pruebas de la extracción Tier 2."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.tiers import (  # noqa: E402
    parse_functional,
    parse_gff_coords,
    read_selected_proteins,
)

FUNCTIONAL = (
    # qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
    # sstart, send, evalue, bitscore, stitle
    "p1\tsp|A|HYP\t40\t100\t0\t0\t1\t100\t1\t100\t1e-20\t80\t"
    "sp|A|HYP Hypothetical protein OS=X\n"
    "p2\tsp|B|REAL\t90\t200\t0\t0\t1\t200\t1\t200\t1e-100\t400\t"
    "sp|B|REAL DNA polymerase III OS=Y\n"
    "p3\tsp|C|DUF\t30\t80\t0\t0\t1\t80\t1\t80\t1e-5\t40\t"
    "sp|C|DUF Protein with DUF1234 domain OS=Z\n"
    "p4\tsp|D|UNC\t35\t150\t0\t0\t1\t150\t1\t150\t1e-9\t60\t"
    "sp|D|UNC Uncharacterized protein OS=W\n"
)

GFF = (
    "##gff-version 3\n"
    "chr\ttest\tCDS\t10\t309\t.\t+\t0\tID=p1;partial=00\n"
    "chr\ttest\tCDS\t400\t999\t.\t-\t0\tID=p2\n"
)


class TestTiers(unittest.TestCase):
    def test_parse_functional(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "functional.tsv"
            p.write_text(FUNCTIONAL)
            sel = parse_functional(p)
            self.assertIn("p1", sel)
            self.assertNotIn("p2", sel)     # enzima conocida
            self.assertIn("p3", sel)        # DUF
            self.assertIn("p4", sel)        # uncharacterized
            # la descripción no debe arrastrar columnas numéricas ni tabs
            self.assertNotIn("\t", sel["p1"])
            self.assertIn("Hypothetical", sel["p1"])

    def test_parse_gff_coords(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.gff3"
            p.write_text(GFF)
            coords = parse_gff_coords(p)
            self.assertEqual(coords["p1"], ("chr", 10, 309, "+"))
            self.assertEqual(coords["p2"], ("chr", 400, 999, "-"))

    def test_read_selected_proteins(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "p.faa"
            p.write_text(">p1\nMKAVL\n>p2\nGGGG\n")
            seqs = read_selected_proteins(p, {"p1"})
            self.assertEqual(seqs, {"p1": "MKAVL"})


if __name__ == "__main__":
    unittest.main()
