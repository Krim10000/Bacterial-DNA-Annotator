"""Pruebas de las regiones vacías."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ghost.empty_regions import gc_n  # noqa: E402


class TestEmptyRegions(unittest.TestCase):
    def test_gc_n(self):
        g, n = gc_n("GGCCNN")
        self.assertAlmostEqual(g, 4 / 6)
        self.assertAlmostEqual(n, 2 / 6)
        self.assertEqual(gc_n(""), (0.0, 0.0))
        self.assertEqual(gc_n("ATAT"), (0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
