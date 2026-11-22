"""Pruebas unitarias del núcleo: secuencias, ORFs y scoring."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.orfs import (  # noqa: E402
    find_orfs_oriented,
    intergenic_regions,
    merge_intervals,
    padded_regions,
)
from ghost.scoring import CodonModel, protein_features, rank_score, rbs_score  # noqa: E402
from ghost.seqio import (  # noqa: E402
    CODON_TABLE,
    read_fasta,
    read_gff_features,
    revcomp,
    translate,
    translate_codon,
)


class TestSeqIO(unittest.TestCase):
    def test_codon_table(self):
        self.assertEqual(CODON_TABLE["ATG"], "M")
        self.assertEqual(CODON_TABLE["TAA"], "*")
        self.assertEqual(CODON_TABLE["TGA"], "*")
        self.assertEqual(CODON_TABLE["TGG"], "W")

    def test_translate_codon(self):
        self.assertEqual(translate_codon("ATG"), ("M", frozenset({"M"})))
        self.assertEqual(translate_codon("TAA"), ("*", frozenset({"*"})))
        # GCN codifica solo Ala
        self.assertEqual(translate_codon("GCN"), ("A", frozenset({"A"})))
        # TTN puede ser Phe o Leu
        aa, poss = translate_codon("TTN")
        self.assertEqual(aa, "X")
        self.assertEqual(poss, frozenset({"F", "L"}))

    def test_translate_to_stop(self):
        self.assertEqual(translate("ATGAAATAA", to_stop=True), "MK")
        self.assertEqual(translate("ATGTAA", to_stop=False), "M*")

    def test_revcomp(self):
        self.assertEqual(revcomp("ATGC"), "GCAT")
        self.assertEqual(revcomp(revcomp("ACGTRYSWKMBDHVN")), "ACGTRYSWKMBDHVN")

    def test_read_fasta(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.fna"
            p.write_text(">a desc\nACGT\nACGT\n>b\nTTTT\n")
            self.assertEqual(list(read_fasta(p)), [("a", "ACGTACGT"), ("b", "TTTT")])

    def test_read_gff(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.gff3"
            p.write_text(
                "##gff-version 3\n##sequence-region chr 1 100\n"
                "chr\ttest\tCDS\t10\t30\t.\t+\t0\tID=c1\n"
                "chr\ttest\tgene\t40\t60\t.\t-\t.\tID=g1\n"
            )
            lengths, feats = read_gff_features(p)
            self.assertEqual(lengths, {"chr": 100})
            self.assertEqual(feats, [("chr", 10, 30, "+", "CDS")])


class TestIntervals(unittest.TestCase):
    def test_merge(self):
        self.assertEqual(
            merge_intervals([(1, 10), (5, 20), (30, 40)]), [(1, 20), (30, 40)]
        )
        self.assertEqual(merge_intervals([(1, 10), (11, 20)]), [(1, 20)])

    def test_intergenic(self):
        self.assertEqual(
            intergenic_regions(100, [(10, 20)], min_len=5), [(1, 9), (21, 100)]
        )
        # ningún hueco suficientemente largo
        self.assertEqual(intergenic_regions(20, [(1, 20)], min_len=5), [])

    def test_padded(self):
        regs = padded_regions(1000, [(100, 200)], min_len=50, pad=30)
        self.assertEqual(regs, [(1, 129), (171, 1000)])
        # con pad, un hueco pequeño puede alcanzar el mínimo
        regs2 = padded_regions(1000, [(100, 200), (220, 300)], min_len=50, pad=30)
        self.assertIn((171, 249), regs2)


class TestORFs(unittest.TestCase):
    def test_frame_0(self):
        seq = "ATG" + "GCT" * 60 + "TAA"
        hits = find_orfs_oriented(seq, min_aa=30)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["start_idx"], 0)
        self.assertEqual(len(hits[0]["aa"]), 61)
        self.assertEqual(hits[0]["stop_codon"], "TAA")

    def test_other_frames(self):
        # ORF desplazado un nucleótido: antes del fix se perdía
        seq = "A" + "ATG" + "GCT" * 40 + "TAA"
        hits = find_orfs_oriented(seq, min_aa=20)
        self.assertTrue(any(h["start_idx"] == 1 for h in hits))

    def test_ambiguity_variants(self):
        seq = "ATG" + "TTN" + "GCT" * 30 + "TAA"
        hit = find_orfs_oriented(seq, min_aa=20)[0]
        self.assertGreaterEqual(hit["n_ambiguous"], 1)
        self.assertTrue(hit["variants"] or "X" in hit["aa"])

    def test_ambiguous_stop(self):
        seq = "ATG" + "GCT" * 5 + "TGR" + "GCT" * 40 + "TAA"
        hits = [h for h in find_orfs_oriented(seq, min_aa=5) if not h["partial"]]
        self.assertEqual(len(hits), 1)
        self.assertTrue(hits[0]["ambig_stop"])
        self.assertEqual(len(hits[0]["aa"]), 6)

    def test_min_length(self):
        seq = "ATG" + "GCT" * 3 + "TAA"
        self.assertEqual(find_orfs_oriented(seq, min_aa=10), [])


class TestScoring(unittest.TestCase):
    def test_codon_model(self):
        m = CodonModel()
        for _ in range(50):
            m.add_cds("ATG" + "GCT" * 40 + "TAA")
            m.add_background("ACGTACGTACGTACGT")
        m.fit()
        sc = m.score("ATG" + "GCT" * 40)
        self.assertIn("codon_llr", sc)
        self.assertGreater(sc["codon_llr"], m.score("ACGTACGTACGT")["codon_llr"])

    def test_rbs(self):
        seq = "AGGAGG" + "A" * 5 + "ATG"
        score, motif = rbs_score(seq, start_idx=11)
        self.assertGreater(score, 0.5)
        self.assertTrue(motif)

    def test_protein_features(self):
        f = protein_features("MKAVL")
        self.assertIn("gravy", f)
        self.assertGreater(f["hydro"], 0.5)

    def test_rank_range(self):
        c = {"codon_llr": 0.6, "aa_len": 120, "rbs": 0.7, "start_w": 1.0, "frac_bad": 0.1}
        r = rank_score(c)
        self.assertTrue(0.0 <= r <= 1.0)


if __name__ == "__main__":
    unittest.main()
