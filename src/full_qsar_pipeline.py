#!/usr/bin/env python3
"""EGFR QSAR pipeline.

ChEMBL 33 EGFR inhibitors -> 12 RDKit 2D descriptors -> Random Forest vs
Gradient Boosting -> 5-fold cross-validation -> permutation importance ->
one four-panel figure and two CSV tables.

Example
-------
    python full_qsar_pipeline.py --data ../data/egfr_pic50.csv --out ../results
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Descriptors, rdMolDescriptors
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_score, train_test_split

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore", category=DeprecationWarning)

DESCRIPTORS = {
    "MW": Descriptors.MolWt,
    "logP": Descriptors.MolLogP,
    "HBD": rdMolDescriptors.CalcNumHBD,
    "HBA": rdMolDescriptors.CalcNumHBA,
    "TPSA": rdMolDescriptors.CalcTPSA,
    "RotBonds": rdMolDescriptors.CalcNumRotatableBonds,
    "AroRings": rdMolDescriptors.CalcNumAromaticRings,
    "FracCSP3": rdMolDescriptors.CalcFractionCSP3,
    "HeavyAtoms": lambda mol: mol.GetNumHeavyAtoms(),
    "qed": QED.qed,
    "NumRings": rdMolDescriptors.CalcNumRings,
    "NOCount": Descriptors.NOCount,
}


def load_dataset(csv_path: Path) -> pd.DataFrame:
    """Read SMILES / pIC50, drop unparsable molecules, cache RDKit mol objects."""
    df = pd.read_csv(csv_path)
    missing = {"smiles", "pIC50"} - set(df.columns)
    if missing:
        raise SystemExit(f"{csv_path} is missing columns: {sorted(missing)}")

    df = df.dropna(subset=["smiles", "pIC50"]).reset_index(drop=True)
    mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    valid = np.array([mol is not None for mol in mols])
    if (~valid).any():
        print(f"  Dropped {int((~valid).sum())} unparsable SMILES")
    df = df.loc[valid].reset_index(drop=True)
    df["mol"] = [mol for mol, ok in zip(mols, valid) if ok]
    return df


def compute_descriptors(df: pd.DataFrame) -> list[str]:
    """Compute every descriptor on the cached mol objects (one parse per molecule)."""
    for name, func in DESCRIPTORS.items():
        df[name] = df["mol"].apply(func)
    return list(DESCRIPTORS)


def fit_and_score(model, X_train, y_train, X_test, y_test) -> dict:
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    return {
        "model": model,
        "train_r2": float(r2_score(y_train, model.predict(X_train))),
        "test_r2": float(r2_score(y_test, pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_test, pred))),
        "pred": pred,
    }


def build_figure(results, y_test, y_all, importances, names, seed, out_png: Path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    rf = results["RandomForest"]

    ax = axes[0, 0]
    ax.scatter(y_test, rf["pred"], c="steelblue", alpha=0.3, s=10)
    lo, hi = float(np.min(y_all)), float(np.max(y_all))
    ax.plot([lo, hi], [lo, hi], "r--", lw=1)
    ax.set_xlabel("True pIC50")
    ax.set_ylabel("Predicted pIC50")
    ax.set_title(f"Random Forest: Test R2 = {rf['test_r2']:.3f}")

    ax = axes[0, 1]
    order = np.argsort(importances)
    ax.barh(np.arange(len(names)), importances[order], color="steelblue")
    ax.set_yticks(np.arange(len(names)))
    ax.set_yticklabels([names[i] for i in order])
    ax.set_xlabel("Permutation importance (drop in R2)")
    ax.set_title("Feature Importance (RF, test set)")

    ax = axes[1, 0]
    residuals = y_test - rf["pred"]
    ax.hist(residuals, bins=50, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(0.0, color="red", linestyle="--")
    ax.set_xlabel("Residual (true - predicted)")
    ax.set_ylabel("Count")
    ax.set_title("Residual Distribution")

    ax = axes[1, 1]
    labels = list(results)
    xs = np.arange(len(labels))
    width = 0.35
    ax.bar(xs - width / 2, [results[k]["train_r2"] for k in labels], width,
           label="Train R2", color="lightblue", edgecolor="navy")
    ax.bar(xs + width / 2, [results[k]["test_r2"] for k in labels], width,
           label="Test R2", color="steelblue", edgecolor="navy")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylabel("R2")
    ax.set_title(f"Train vs Test (seed = {seed})")
    ax.legend()

    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", default="../data/egfr_pic50.csv", type=Path,
                        help="CSV with 'smiles' and 'pIC50' columns")
    parser.add_argument("--out", default="../results", type=Path,
                        help="output directory for figures and tables")
    parser.add_argument("--seed", default=42, type=int)
    parser.add_argument("--test-size", default=0.2, type=float)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("EGFR QSAR Pipeline")
    print("=" * 60)

    print("\n[1] Loading data ...")
    df = load_dataset(args.data)
    print(f"  Molecules: {len(df)}")

    print("\n[2] Computing descriptors ...")
    names = compute_descriptors(df)
    X = df[names].to_numpy(dtype=float)
    y = df["pIC50"].to_numpy(dtype=float)
    print(f"  Features: {X.shape[1]}  |  pIC50 range: {y.min():.2f} - {y.max():.2f}")

    print("\n[3] Train/test split and model fitting ...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.seed
    )
    models = {
        "RandomForest": RandomForestRegressor(
            n_estimators=200, random_state=args.seed, n_jobs=-1
        ),
        "GradientBoost": GradientBoostingRegressor(
            n_estimators=200, max_depth=5, random_state=args.seed
        ),
    }
    results = {}
    for name, model in models.items():
        results[name] = fit_and_score(model, X_train, y_train, X_test, y_test)
        r = results[name]
        print(f"  {name:14s} Train R2={r['train_r2']:.4f}  "
              f"Test R2={r['test_r2']:.4f}  RMSE={r['rmse']:.3f}")

    print("\n[4] 5-fold cross-validation ...")
    cv = KFold(n_splits=5, shuffle=True, random_state=args.seed)
    cv_scores = {}
    for name, model in models.items():
        scores = cross_val_score(model, X, y, cv=cv, scoring="r2", n_jobs=-1)
        cv_scores[name] = scores
        print(f"  {name:14s} CV R2 = {scores.mean():.4f} +/- {scores.std():.4f}")

    print("\n[5] Permutation importance (RandomForest, test set) ...")
    rf = results["RandomForest"]["model"]
    perm = permutation_importance(
        rf, X_test, y_test, n_repeats=10, random_state=args.seed, n_jobs=-1
    )
    importance = perm.importances_mean
    for i, idx in enumerate(np.argsort(importance)[::-1]):
        print(f"  {i + 1:2d}. {names[idx]:12s} {importance[idx]:+.4f}")

    print("\n[6] Writing figure and tables ...")
    fig_path = args.out / "qsar_full_analysis.png"
    build_figure(results, y_test, y, importance, names, args.seed, fig_path)

    importance_df = pd.DataFrame({
        "Descriptor": names,
        "PermutationImportance": importance,
        "ImpurityImportance": rf.feature_importances_,
    }).sort_values("PermutationImportance", ascending=False)
    importance_df.to_csv(args.out / "feature_importance.csv", index=False)

    pred_df = pd.DataFrame({
        "True_pIC50": y_test,
        "Predicted_pIC50": results["RandomForest"]["pred"],
        "Residual": y_test - results["RandomForest"]["pred"],
    })
    pred_df.to_csv(args.out / "predictions.csv", index=False)

    for path in (fig_path, args.out / "feature_importance.csv", args.out / "predictions.csv"):
        print(f"  Saved: {path}")

    best = max(results, key=lambda k: results[k]["test_r2"])
    print("\n" + "=" * 60)
    print(f"Best model : {best}")
    print(f"Test R2    : {results[best]['test_r2']:.3f}")
    print(f"5-fold CV  : {cv_scores[best].mean():.3f} +/- {cv_scores[best].std():.3f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
