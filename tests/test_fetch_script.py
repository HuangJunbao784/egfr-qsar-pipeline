"""Tests for the ChEMBL fetch helper (no network access required)."""
from __future__ import annotations

import importlib.util
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("fetch_chembl_egfr", ROOT / "scripts" / "fetch_chembl_egfr.py")
fetch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fetch)


class Pic50Tests(unittest.TestCase):
    def test_unit_conversion(self):
        self.assertAlmostEqual(fetch.to_pic50(1.0), 9.0)
        self.assertAlmostEqual(fetch.to_pic50(10.0), 8.0)
        self.assertAlmostEqual(fetch.to_pic50(1000.0), 6.0)

    def test_matches_formula(self):
        for value in (0.5, 3.2, 250.0, 10000.0):
            self.assertAlmostEqual(fetch.to_pic50(value), 9.0 - math.log10(value))

    def test_default_target_is_egfr(self):
        self.assertEqual(fetch.EGFR, "CHEMBL203")
        self.assertTrue(fetch.BASE.startswith("https://"))


if __name__ == "__main__":
    unittest.main()
