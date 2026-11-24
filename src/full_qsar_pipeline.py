#!/usr/bin/env python3
"""EGFR QSAR pipeline.

6,896 ChEMBL EGFR inhibitors -> 2D descriptors and/or Morgan fingerprints ->
Random Forest vs Gradient Boosting -> random and Murcko-scaffold splits ->
5-fold CV -> permutation or impurity importance.

Examples
--------
    python full_qsar_pipeline.py --data ../data/egfr_pic50.csv --out ../results
    python full_qsar_pipeline.py --features both --split both
    egfr-qsar --features fp --split scaffold      # after pip install -e .
"""
from __future__ import annotations

import argparse
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold, cross_val_score

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

MODEL_FACTORIES = {
    "RandomForest": lambda seed: RandomForestRegressor(
        n_estimators=200, random_state=seed, n_jobs=-1
    ),
    "GradientBoost": lambda seed: GradientBoostingRegressor(
        n_estimators=200, max_depth=5, random_state=seed
    ),
}

FP_RADIUS = 2
FP_BITS = 2048


def resolve_data(path: Path) -> Path:
    """Accept an explicit path, the current directory, or the repository layout."""
    repo_data = Path(__file__).resolve().parent.parent / "data" / "egfr_pic50.csv"
    for candidate in (path, Path("data/egfr_pic50.csv"), repo_data):
        if candidate.exists():
            return candidate
    raise SystemExit(f"Data file not found: {path} (pass --data explicitly)")


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


def build_features(df: pd.DataFrame, kind: str):
    """Return (matrix, feature names) for the requested feature set."""
    desc_names = compute_descriptors(df)
    desc = df[desc_names].to_numpy(dtype=float)
    if kind == "desc":
        return desc, desc_names

    fingerprints = np.array([
        list(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, nBits=FP_BITS))
        for mol in df["mol"]
    ], dtype=np.float32)
    fp_names = [f"morgan_{i:04d}" for i in range(fingerprints.shape[1])]
    if kind == "fp":
        return fingerprints, fp_names
    return np.hstack([desc, fingerprints]), desc_names + fp_names


def make_split(df: pd.DataFrame, method: str, test_size: float, seed: int):
    """Return (train_idx, test_idx, group_labels) for a random or scaffold split.

    The scaffold split keeps every Murcko scaffold entirely on one side, so the
    test set measures generalisation to unseen chemotypes rather than to close
    analogues.
    """
    n = len(df)
    indices = np.arange(n)
    if method == "random":
        np.random.default_rng(seed).shuffle(indices)
        cut = int(round((1.0 - test_size) * n))
        return indices[:cut], indices[cut:], None

    scaffolds = np.array([
        MurckoScaffold.MurckoScaffoldSmiles(mol=mol) for mol in df["mol"]
    ])
    by_scaffold = defaultdict(list)
    for i, scaffold in enumerate(scaffolds):
        by_scaffold[scaffold].append(i)

    groups = sorted(by_scaffold)
    np.random.default_rng(seed).shuffle(groups)

    target = int(round(test_size * n))
    test_idx: list[int] = []
    for scaffold in groups:
        if len(test_idx) >= target:
            break
        test_idx.extend(by_scaffold[scaffold])

    test_idx = np.array(sorted(test_idx))
    train_idx = np.array(sorted(set(range(n)) - set(test_idx)))
    return train_idx, test_idx, scaffolds


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


def cross_validate(model, X, y, method: str, groups, seed: int):
    if method == "random":
        cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    else:
        cv = GroupKFold(n_splits=5)
    scores = cross_val_score(model, X, y, cv=cv, groups=groups, scoring="r2", n_jobs=-1)
    return float(scores.mean()), float(scores.std())


def run_seed_sweep(df, X, y, method: str, seeds, test_size: float, out_dir: Path) -> None:
    print(f"\n[seed sweep] split={method}")
    rows = []
    for seed in seeds:
        train_idx, test_idx, groups = make_split(df, method, test_size, seed)
        for name, factory in MODEL_FACTORIES.items():
            res = fit_and_score(factory(seed), X[train_idx], y[train_idx], X[test_idx], y[test_idx])
            cv_mean, cv_std = cross_validate(factory(seed), X, y, method, groups, seed)
            rows.append({
                "seed": seed, "split": method, "model": name,
                "train_r2": res["train_r2"], "test_r2": res["test_r2"],
                "cv_r2_mean": cv_mean, "cv_r2_std": cv_std, "rmse": res["rmse"],
            })
    per_seed = pd.DataFrame(rows)
    per_seed.to_csv(out_dir / "seed_sweep.csv", index=False)
    summary = per_seed.groupby("model").agg(
        test_r2_mean=("test_r2", "mean"), test_r2_std=("test_r2", "std"),
        cv_r2_mean=("cv_r2_mean", "mean"), cv_r2_std=("cv_r2_std", "mean"),
    ).round(4)
    print(summary.to_string())
    print(f"  Saved: {out_dir / 'seed_sweep.csv'}")


def build_figure(results, y_test, y_all, importances, names, title, out_png: Path):
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
    top = np.argsort(importances)[::-1][:15]
    ax.barh(np.arange(len(top)), importances[top][::-1], color="steelblue")
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels([names[i] for i in top[::-1]], fontsize=8)
    ax.set_xlabel("Importance")
    ax.set_title("Top-15 features")

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
    ax.set_title(title)
    ax.legend()

    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", default="../data/egfr_pic50.csv", type=Path)
    parser.add_argument("--out", default="../results", type=Path)
    parser.add_argument("--seed", default=42, type=int)
    parser.add_argument("--seeds", default=None,
                        help="comma-separated seeds for a sweep of the first split, e.g. 42,7,2026")
    parser.add_argument("--test-size", default=0.2, type=float)
    parser.add_argument("--features", default="desc", choices=["desc", "fp", "both"],
                        help="desc = 12 2D descriptors, fp = 2048-bit Morgan fingerprints")
    parser.add_argument("--split", default="random", choices=["random", "scaffold", "both"],
                        help="random split or Murcko-scaffold split (unseen chemotypes)")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    data_path = resolve_data(args.data)

    print("=" * 60)
    print("EGFR QSAR Pipeline")
    print("=" * 60)

    print("\n[1] Loading data ...")
    df = load_dataset(data_path)
    print(f"  Molecules: {len(df)}")

    print(f"\n[2] Building features ({args.features}) ...")
    X, names = build_features(df, args.features)
    y = df["pIC50"].to_numpy(dtype=float)
    print(f"  Features: {X.shape[1]}  |  pIC50 range: {y.min():.2f} - {y.max():.2f}")

    methods = ["random", "scaffold"] if args.split == "both" else [args.split]
    rows, per_split = [], {}
    for method in methods:
        train_idx, test_idx, groups = make_split(df, method, args.test_size, args.seed)
        print(f"\n[3] Split={method}: train {len(train_idx)} / test {len(test_idx)}")
        results = {}
        for name, factory in MODEL_FACTORIES.items():
            res = fit_and_score(factory(args.seed), X[train_idx], y[train_idx],
                                X[test_idx], y[test_idx])
            cv_mean, cv_std = cross_validate(factory(args.seed), X, y, method, groups, args.seed)
            results[name] = res
            rows.append({
                "split": method, "features": args.features, "model": name,
                "train_r2": res["train_r2"], "test_r2": res["test_r2"],
                "cv_r2_mean": cv_mean, "cv_r2_std": cv_std, "rmse": res["rmse"],
            })
            print(f"  {name:14s} Train={res['train_r2']:.3f}  Test={res['test_r2']:.3f}  "
                  f"CV={cv_mean:.3f} +/- {cv_std:.3f}")
        per_split[method] = (results, test_idx)

    comparison = pd.DataFrame(rows)
    comparison.to_csv(args.out / "split_comparison.csv", index=False)
    print("\n" + comparison.round(4).to_string(index=False))
    print(f"\n  Saved: {args.out / 'split_comparison.csv'}")

    primary = methods[0]
    results, test_idx = per_split[primary]
    y_test = y[test_idx]

    if args.features == "desc":
        print("\n[4] Permutation importance (RandomForest, test set) ...")
        rf = results["RandomForest"]["model"]
        importance = permutation_importance(
            rf, X[test_idx], y_test, n_repeats=10, random_state=args.seed, n_jobs=-1
        ).importances_mean
        for i, idx in enumerate(np.argsort(importance)[::-1]):
            print(f"  {i + 1:2d}. {names[idx]:12s} {importance[idx]:+.4f}")
    else:
        print("\n[4] Fingerprints present: reporting impurity importance (top 15).")
        importance = results["RandomForest"]["model"].feature_importances_
        for i, idx in enumerate(np.argsort(importance)[::-1][:15]):
            print(f"  {i + 1:2d}. {names[idx]:14s} {importance[idx]:.4f}")

    print("\n[5] Writing figure and tables ...")
    fig_path = args.out / "qsar_full_analysis.png"
    build_figure(results, y_test, y, importance, names,
                 f"{args.features} features, {primary} split", fig_path)

    importance_df = pd.DataFrame({
        "Descriptor": names, "Importance": importance,
    }).sort_values("Importance", ascending=False)
    importance_df.to_csv(args.out / "feature_importance.csv", index=False)

    pred_df = pd.DataFrame({
        "True_pIC50": y_test,
        "Predicted_pIC50": results["RandomForest"]["pred"],
        "Residual": y_test - results["RandomForest"]["pred"],
    })
    pred_df.to_csv(args.out / "predictions.csv", index=False)
    for path in (fig_path, args.out / "feature_importance.csv", args.out / "predictions.csv"):
        print(f"  Saved: {path}")

    if args.seeds:
        seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
        run_seed_sweep(df, X, y, primary, seeds, args.test_size, args.out)

    print("\n" + "=" * 60)
    print(f"Features: {args.features}  |  split: {args.split}")
    print("=" * 60)


if __name__ == "__main__":
    main()
