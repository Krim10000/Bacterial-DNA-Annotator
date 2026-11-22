"""Pruebas del estado y los contadores del motor de genes fantasma."""
import unittest

from ghost.status import STATUSES, Decision, apply_status, decide


class TestStatus(unittest.TestCase):
    def test_valid(self):
        d = decide("validated", validated_by="ORCA + DoriC",
                   reason="coinciden", attempted="ORCA; GC skew")
        self.assertEqual(d.status, "validated")

    def test_invalid_status(self):
        with self.assertRaises(ValueError):
            decide("clarísimo")

    def test_invalid_confidence(self):
        with self.assertRaises(ValueError):
            decide("validated", confidence="alta")

    def test_apply(self):
        row = apply_status({}, decide("anomalous", reason="Δ grande",
                                      attempted="GC skew"))
        self.assertEqual(row["status"], "anomalous")
        self.assertEqual(row["status_reason"], "Δ grande")
        self.assertEqual(row["status_confidence"], "")

    def test_all_statuses_documented(self):
        for status in ("validated", "consensus", "lineage_specific",
                       "anomalous", "unresolved", "novel"):
            self.assertIn(status, STATUSES)


if __name__ == "__main__":
    unittest.main()
