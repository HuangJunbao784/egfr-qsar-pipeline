"""Unit tests for the EGFR QSAR pipeline (standard library unittest only)."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import full_qsar_pipeline as fqp  # noqa: E402


class DescriptorTests(unittest.TestCase):
    def test_descriptor_list(self):
        self.assertEqual(len(fqp.DESCRIPTORS), 12)
        for name in ("MW", "logP", "TPSA", "HBA", "HBD", "RotBonds", "AroRings"):
            self.assertIn(name, fqp.DESCRIPTORS)

    def test_ethanol_molecular_weight(self):
        df = pd.DataFrame({"mol": [fqp.Chem.MolFromSmiles("CCO")]})
        fqp.compute_descriptors(df)
        self.assertAlmostEqual(df["MW"].iloc[0], 46.069, places=2)
        self.assertEqual(int(df["HBD"].iloc[0]), 1)
        self.assertEqual(int(df["HBA"].iloc[0]), 1)
        self.assertEqual(int(df["AroRings"].iloc[0]), 0)


class DataLoadingTests(unittest.TestCase):
    def _write_csv(self, text: str) -> Path:
        handle = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, encoding="utf-8")
        handle.write(text)
        handle.close()
        return Path(handle.name)

    def test_invalid_smiles_are_dropped(self):
        path = self._write_csv("smiles,pIC50\nCCO,5.0\nnot_a_smiles,6.0\nc1ccccc1,7.0\n")
        df = fqp.load_dataset(path)
        self.assertEqual(len(df), 2)
        self.assertTrue(df["mol"].notna().all())

    def test_missing_columns_raise(self):
        path = self._write_csv("smiles,activity\nCCO,5.0\n")
        with self.assertRaises(SystemExit):
            fqp.load_dataset(path)

    def test_resolve_data_accepts_explicit_path(self):
        path = self._write_csv("smiles,pIC50\nCCO,5.0\n")
        self.assertEqual(fqp.resolve_data(path), path)

    def test_resolve_data_finds_repository_copy(self):
        repo_data = Path(fqp.__file__).resolve().parent.parent / "data" / "egfr_pic50.csv"
        if not repo_data.exists():
            self.skipTest("bundled dataset not present")
        resolved = fqp.resolve_data(Path("definitely-missing.csv"))
        self.assertTrue(resolved.exists())
        self.assertGreater(resolved.stat().st_size, 0)


class ModelTests(unittest.TestCase):
    def test_fit_and_score_returns_expected_keys(self):
        rng = np.random.default_rng(0)
        X = rng.normal(size=(60, 3))
        y = 2.0 * X[:, 0] - X[:, 1] + rng.normal(scale=0.1, size=60)
        model = fqp.MODEL_FACTORIES["RandomForest"](seed=1)
        result = fqp.fit_and_score(model, X[:40], y[:40], X[40:], y[40:])
        for key in ("model", "train_r2", "test_r2", "rmse", "pred"):
            self.assertIn(key, result)
        self.assertTrue(np.isfinite(result["test_r2"]))
        self.assertEqual(len(result["pred"]), 20)

    def test_cross_validate_returns_mean_and_std(self):
        rng = np.random.default_rng(1)
        X = rng.normal(size=(50, 3))
        y = X[:, 0] + rng.normal(scale=0.1, size=50)
        mean, std = fqp.cross_validate(fqp.MODEL_FACTORIES["RandomForest"](seed=2), X, y, seed=2)
        self.assertTrue(np.isfinite(mean))
        self.assertGreaterEqual(std, 0.0)


if __name__ == "__main__":
    unittest.main()
