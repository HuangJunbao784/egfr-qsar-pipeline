"""
EGFR QSAR Pipeline — ChEMBL 33 Data
=====================================
从ChEMBL下载EGFR抑制剂 → RDKit特征 → 多模型对比 → 交叉验证 → 特征重要性 → 可视化

Author: HuangJunbao
Date: 2026-06
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, cross_val_score, KFold, GridSearchCV
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings("ignore")

# ============================================
# 1. 加载数据
# ============================================
print("=" * 60)
print("EGFR QSAR Pipeline")
print("=" * 60)

df = pd.read_csv("../data/egfr_pic50.csv")
print(f"\n[1] 数据量: {len(df)} molecules")

# ============================================
# 2. 计算分子描述符
# ============================================
print(f"\n[2] Computing molecular descriptors...")

# 2D描述符
df["MW"] = df["smiles"].apply(lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)))
df["logP"] = df["smiles"].apply(lambda s: Descriptors.MolLogP(Chem.MolFromSmiles(s)))
df["HBD"] = df["smiles"].apply(lambda s: Descriptors.NumHDonors(Chem.MolFromSmiles(s)))
df["HBA"] = df["smiles"].apply(lambda s: Descriptors.NumHAcceptors(Chem.MolFromSmiles(s)))
df["TPSA"] = df["smiles"].apply(lambda s: Descriptors.TPSA(Chem.MolFromSmiles(s)))
df["RotBonds"] = df["smiles"].apply(lambda s: Descriptors.NumRotatableBonds(Chem.MolFromSmiles(s)))
df["AroRings"] = df["smiles"].apply(lambda s: Descriptors.NumAromaticRings(Chem.MolFromSmiles(s)))
df["FracCSP3"] = df["smiles"].apply(lambda s: Descriptors.FractionCSP3(Chem.MolFromSmiles(s)))
df["HeavyAtoms"] = df["smiles"].apply(lambda s: Descriptors.HeavyAtomCount(Chem.MolFromSmiles(s)))
df["qed"] = df["smiles"].apply(lambda s: Descriptors.qed(Chem.MolFromSmiles(s)))
df["NumRings"] = df["smiles"].apply(lambda s: Descriptors.RingCount(Chem.MolFromSmiles(s)))
df["NOCount"] = df["smiles"].apply(lambda s: Descriptors.NOCount(Chem.MolFromSmiles(s)))

descriptor_names = ["MW", "logP", "HBD", "HBA", "TPSA", "RotBonds",
                    "AroRings", "FracCSP3", "HeavyAtoms", "qed", "NumRings", "NOCount"]

# 特征矩阵和目标值
X_desc = df[descriptor_names].values
y = df["pIC50"].values

print(f"  Descriptors: {X_desc.shape[1]} features")
print(f"   pIC50 range: {y.min():.1f} ~ {y.max():.1f}")

# ============================================
# 3. Train/Test Split + 模型对比
# ============================================
print(f"\n[3] Model comparison...")

X_train, X_test, y_train, y_test = train_test_split(
    X_desc, y, test_size=0.2, random_state=42
)

models = {
    "RandomForest": RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1),
    "GradientBoost": GradientBoostingRegressor(n_estimators=200, max_depth=5, random_state=42),
}

results = {}
for name, model in models.items():
    model.fit(X_train, y_train)
    train_r2 = model.score(X_train, y_train)
    test_r2 = model.score(X_test, y_test)
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    results[name] = {"model": model, "train_r2": train_r2, "test_r2": test_r2, "rmse": rmse}
    print(f"  {name:15s}: Train R2={train_r2:.4f}, Test R2={test_r2:.4f}, RMSE={rmse:.3f}")

# ============================================
# 4. 5折交叉验证
# ============================================
print(f"\n[4] 5-fold Cross-validation...")

cv = KFold(n_splits=5, shuffle=True, random_state=42)
for name, model in models.items():
    scores = cross_val_score(model, X_desc, y, cv=cv, scoring="r2", n_jobs=-1)
    print(f"  {name:15s}: CV R2 = {scores.mean():.4f} (+/- {scores.std():.4f})")

# ============================================
# 5. 特征重要性
# ============================================
print(f"\n[5] Feature importance (RandomForest):")

best_rf = models["RandomForest"]
importances = best_rf.feature_importances_
sorted_idx = np.argsort(importances)[::-1]

for i, idx in enumerate(sorted_idx):
    bar = "#" * int(importances[idx] * 100)
    print(f"  {i+1:2d}. {descriptor_names[idx]:12s}: {importances[idx]:.3f} {bar}")

# ============================================
# 6. 可视化
# ============================================
print(f"\n[6] Generating figures...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# (a) 预测 vs 真实
rf = results["RandomForest"]["model"]
y_pred_rf = rf.predict(X_test)
axes[0, 0].scatter(y_test, y_pred_rf, c="steelblue", alpha=0.3, s=10)
axes[0, 0].plot([y.min(), y.max()], [y.min(), y.max()], "r--", lw=1)
axes[0, 0].set_xlabel("True pIC50")
axes[0, 0].set_ylabel("Predicted pIC50")
axes[0, 0].set_title(f"Random Forest: Test R2={results['RandomForest']['test_r2']:.3f}")

# (b) 特征重要性棒图
axes[0, 1].barh(range(12), importances[sorted_idx][::-1][-12:], color="steelblue")
axes[0, 1].set_yticks(range(12))
axes[0, 1].set_yticklabels([descriptor_names[i] for i in sorted_idx[::-1]][-12:])
axes[0, 1].set_xlabel("Importance")
axes[0, 1].set_title("Feature Importance (RF)")

# (c) 残差分布
residuals = y_test - y_pred_rf
axes[1, 0].hist(residuals, bins=50, color="steelblue", edgecolor="white", alpha=0.8)
axes[1, 0].axvline(x=0, color="red", linestyle="--")
axes[1, 0].set_xlabel("Residual (True - Predicted)")
axes[1, 0].set_ylabel("Count")
axes[1, 0].set_title("Residual Distribution")

# (d) 模型对比
model_names = list(results.keys())
train_scores = [results[n]["train_r2"] for n in model_names]
test_scores = [results[n]["test_r2"] for n in model_names]

x = np.arange(len(model_names))
w = 0.35
axes[1, 1].bar(x - w/2, train_scores, w, label="Train R2", color="lightblue", edgecolor="navy")
axes[1, 1].bar(x + w/2, test_scores, w, label="Test R2", color="steelblue", edgecolor="navy")
axes[1, 1].set_xticks(x)
axes[1, 1].set_xticklabels(model_names)
axes[1, 1].set_ylabel("R2 Score")
axes[1, 1].set_title("Train vs Test Performance")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("../results/qsar_full_analysis.png", dpi=200, bbox_inches="tight")
plt.close()
print("   Saved: results/qsar_full_analysis.png")

# ============================================
# 7. 保存结果
# ============================================
print(f"\n[7] Saving results...")

# 特征重要性表
importance_df = pd.DataFrame({
    "Descriptor": descriptor_names,
    "Importance": importances
}).sort_values("Importance", ascending=False)
importance_df.to_csv("../results/feature_importance.csv", index=False)

# 预测结果
pred_df = pd.DataFrame({
    "True_pIC50": y_test,
    "Predicted_pIC50": y_pred_rf,
    "Residual": residuals
})
pred_df.to_csv("../results/predictions.csv", index=False)

print("   Saved: results/feature_importance.csv")
print("   Saved: results/predictions.csv")
print(f"\n{'='*60}")
print(f"Pipeline complete!")
print(f"  Best model: RandomForest")
print(f"  Test R2:    {results['RandomForest']['test_r2']:.3f}")
print(f"  5-fold CV:  {cross_val_score(best_rf, X_desc, y, cv=5, scoring='r2').mean():.3f}")
print(f"{'='*60}")
