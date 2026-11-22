"""Pruebas de los helpers del experimento de deduplicación."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.dedup_eval import metrics  # noqa: E402
from ghost.dedup_pairs import cluster_of, parse_allvall  # noqa: E402


class TestDedupHelpers(unittest.TestCase):
    def test_parse_allvall(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "a.tsv"
            p.write_text(
                "a\ta\t100\t100\t100\t100\n"        # self: se descarta
                "a\tb\t95\t80\t120\t120\n"          # cobertura 0.67: se descarta
                "a\tc\t90\t100\t100\t100\n"         # válido
                "c\ta\t92\t100\t100\t100\n"         # mismo par, mayor identidad
                "b\tc\t99\t100\t100\t100\n"         # válido
            )
            pairs = parse_allvall(p)
            self.assertNotIn(("a", "a"), pairs)
            self.assertNotIn(("a", "b"), pairs)
            self.assertEqual(pairs[("a", "c")], 92.0)
            self.assertEqual(pairs[("b", "c")], 99.0)

    def test_cluster_of(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "c.tsv"
            p.write_text("rep1\ta\nrep1\tb\nrep2\tc\n")
            cmap = cluster_of(p)
            self.assertEqual(cmap, {"a": "rep1", "b": "rep1", "c": "rep2"})

    def test_metrics_pure(self):
        clusters = {"r1": ["a", "b"], "r2": ["c"]}
        key = {"a": "k1", "b": "k1", "c": "k2"}
        m = metrics(clusters, key)
        self.assertEqual(m["sequences"], 3)
        self.assertEqual(m["clusters"], 2)
        self.assertEqual(m["mixed_clusters"], 0)
        self.assertEqual(m["pair_precision"], 1.0)
        self.assertEqual(m["pair_recall"], 1.0)

    def test_metrics_mixed(self):
        clusters = {"r1": ["a", "b", "c"]}
        key = {"a": "k1", "b": "k2", "c": "k3"}
        m = metrics(clusters, key)
        self.assertEqual(m["mixed_clusters"], 1)
        self.assertEqual(m["pair_precision"], 0.0)
        self.assertEqual(m["pair_recall"], 1.0)


if __name__ == "__main__":
    unittest.main()
