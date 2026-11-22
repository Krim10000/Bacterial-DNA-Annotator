"""Prueba end-to-end sobre un genoma sintético (sin datos del proyecto)."""
from __future__ import annotations

import csv
import gzip
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.pipeline import build_candidates, write_outputs  # noqa: E402


def make_genome(root: Path) -> tuple[str, int]:
    """Genoma con un CDS conocido, un espaciador y un ORF con codones ambiguos."""
    genome = "GCF_TEST.1_ASM1v1"
    (root / "raw").mkdir(parents=True)
    (root / "annotations_ab_initio" / genome).mkdir(parents=True)

    rng = random.Random(7)
    known = "ATG" + "GCT" * 100 + "TAA"
    spacer = "".join(rng.choice("ACGT") for _ in range(60))
    ambiguous = "ATG" + "GCN" * 6 + "TTN" * 4 + "GCT" * 40 + "TAA"
    contig = known + spacer + ambiguous

    with gzip.open(root / "raw" / f"{genome}_genomic.fna.gz", "wt") as fh:
        fh.write(">chr\n")
        for i in range(0, len(contig), 70):
            fh.write(contig[i:i + 70] + "\n")

    with open(root / "annotations_ab_initio" / genome / "annotation.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region chr 1 {len(contig)}\n")
        fh.write(f"chr\ttest\tCDS\t1\t{len(known)}\t.\t+\t0\tID=cds1\n")
    return genome, len(known)


class TestSyntheticEndToEnd(unittest.TestCase):
    def test_pipeline(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "data"
            genome, known_len = make_genome(root)
            res = build_candidates(
                root, genome, min_region_len=30, min_aa=20, max_variants=4
            )

            # el CDS conocido no debe reaparecer como candidato
            self.assertTrue(res.n_orfs >= 1)
            for c in res.candidates:
                self.assertGreater(c["start"], known_len,
                                   "se colaron fragmentos del CDS anotado")

            # el ORF ambiguo debe existir y estar marcado
            amb = [c for c in res.candidates if c["n_ambiguous"] > 0]
            self.assertTrue(amb, "no se generó el ORF con ambigüedad")
            self.assertTrue(any("X" in c["aa"] or c["variants"] > 0 for c in amb))

            # campos mínimos
            c = res.candidates[0]
            for key in ("codon_llr", "rbs", "rank", "ann_frac", "aa_len", "start_codon"):
                self.assertIn(key, c)
            self.assertTrue(0.0 <= c["rank"] <= 1.0)

            # salidas
            paths = write_outputs(res, Path(td) / "out")
            self.assertTrue(paths["tsv"].exists())
            self.assertTrue(paths["json"].exists())
            rows = list(csv.DictReader(open(paths["tsv"]), delimiter="\t"))
            self.assertEqual(len(rows), len(res.candidates))
            doc = json.load(open(paths["json"]))
            self.assertIn("stats", doc)
            self.assertIn("candidates", doc)


if __name__ == "__main__":
    unittest.main()
