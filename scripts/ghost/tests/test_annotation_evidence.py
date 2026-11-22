"""Pruebas de la evidencia de anotación del motor de genes fantasma."""
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "webapp"))

from family_logic import explain_families  # noqa: E402
from producer_homology import annotate as annotate_homologs  # noqa: E402
from transport_layer import annotate_gaps as annotate_transport  # noqa: E402
from pipeline import write_outputs  # noqa: E402


class TestFamilyLogic(unittest.TestCase):
    def test_present_and_sin_ruta(self):
        gaps = [
            {"gap": "A", "pathways": "map00260 Glycine, serine and threonine metabolism",
             "homolog_status": "homolog_strong", "transport_status": "no_transport"},
            {"gap": "B", "pathways": "sin ruta", "homolog_status": "no_homolog",
             "transport_status": "no_transport"},
        ]
        out = {r["family"]: r for r in explain_families(gaps)}
        self.assertIn("enzima presente no mapeada", out["map00260 Glycine, serine and threonine metabolism"]["conclusion"])
        self.assertIn("mapeo KEGG incompleto", out["sin ruta"]["conclusion"])


class TestProducerHomology(unittest.TestCase):
    def test_strong_and_not_checked(self):
        gaps = [
            {"gap": "A", "producers": ["RHEA:1"], "chebi": "C"},
            {"gap": "B", "producers": ["RHEA:9"], "chebi": "D"},
        ]
        r2a = {"RHEA:1": {"ACC1"}}
        hits = {"ACC1": {"acc": "ACC1", "pident": 50.0, "qcov": 0.9,
                         "bits": 100.0, "desc": "enzima"}}
        annotate_homologs(gaps, hits, r2a, {"ACC1"})
        self.assertEqual(gaps[0]["homolog_status"], "homolog_strong")
        self.assertEqual(gaps[0]["homolog_rhea"], "RHEA:1")
        self.assertEqual(gaps[1]["homolog_status"], "not_checked")


class TestTransportLayer(unittest.TestCase):
    def test_t1_and_no_transport(self):
        gaps = [{"gap": "A", "chebi": "C"}, {"gap": "B", "chebi": "D"}]
        tbc = {"C": {"RHEA:2"}}
        r2a = {"RHEA:2": {"TR"}}
        hits = {"TR": {"acc": "TR", "pident": 50.0, "qcov": 0.9,
                       "evalue": 1e-30, "bits": 100.0, "desc": "transportador"}}
        annotate_transport(gaps, hits, tbc, r2a)
        self.assertEqual(gaps[0]["transport_status"], "t1_bacterial")
        self.assertEqual(gaps[1]["transport_status"], "no_transport")


class TestPipelineOutputs(unittest.TestCase):
    def test_write_outputs(self):
        res = {"genome": "demo", "name": "demo.fna", "bp": 100, "gc": 0.5,
               "gff": "##gff-version 3\n",
               "proteins_fasta": ">p\nMK\n",
               "candidates_tsv": "contig\tstart\n",
               "report_md": "# Informe\n",
               "features": [], "deserts": [], "candidates": [],
               "coverage": 0.9, "annotated_bp": 90}
        with tempfile.TemporaryDirectory() as td:
            write_outputs(res, Path(td), None)
            for name in ("demo.gff3", "demo.proteins.faa", "demo.candidates.tsv",
                         "demo.report.md", "demo.summary.json"):
                self.assertTrue((Path(td) / name).exists(), name)


if __name__ == "__main__":
    unittest.main()
