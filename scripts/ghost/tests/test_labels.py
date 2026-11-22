"""Pruebas del esquema de etiquetado."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.labels import LABELS, label_row  # noqa: E402


class TestLabels(unittest.TestCase):
    def test_valid(self):
        row = label_row({}, label="desert_candidate", source="s",
                        evidence="e", confidence="high", interpretation="i")
        for k in ("label", "label_source", "label_evidence", "label_confidence",
                  "label_interpretation", "label_pipeline", "label_date"):
            self.assertIn(k, row)
        self.assertEqual(row["label"], "desert_candidate")

    def test_unknown_label(self):
        with self.assertRaises(ValueError):
            label_row({}, label="inventada", source="s", evidence="e",
                      confidence="high", interpretation="i")

    def test_bad_confidence(self):
        with self.assertRaises(ValueError):
            label_row({}, label="desert_candidate", source="s", evidence="e",
                      confidence="alta", interpretation="i")

    def test_all_labels_documented(self):
        self.assertTrue(all(isinstance(v, str) and v for v in LABELS.values()))


if __name__ == "__main__":
    unittest.main()
