"""Verifica que todas las CLIs del paquete arrancan y muestran ayuda."""
from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULES = [
    "run",
    "analyze",
    "validate",
    "report",
    "structure",
    "gap_link",
    "gap_targets",
    "gap_demand",
    "gap_conservation",
    "biomass",
]


class TestCLIs(unittest.TestCase):
    def test_help(self):
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")
        for mod in MODULES:
            with self.subTest(module=mod):
                res = subprocess.run(
                    [sys.executable, "-m", f"ghost.{mod}", "--help"],
                    capture_output=True, text=True, env=env, cwd=ROOT, timeout=120,
                )
                self.assertEqual(res.returncode, 0, f"{mod}: {res.stderr[-400:]}")
                self.assertIn("usage", (res.stdout + res.stderr).lower())


if __name__ == "__main__":
    unittest.main()
