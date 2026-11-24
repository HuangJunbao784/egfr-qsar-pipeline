# EGFR QSAR Pipeline

**ChEMBL 33 · RDKit 2D descriptors · Random Forest / Gradient Boosting · 5-fold cross-validation**

从 SMILES 出发，用 12 个二维分子描述符预测 EGFR 激酶抑制剂的 pIC50，比较两种集成模型，并给出交叉验证、置换重要性、残差与特征排序。

## 快速开始

```bash
pip install -r requirements.txt
python src/full_qsar_pipeline.py --data data/egfr_pic50.csv --out results
```

脚本会在 `results/` 生成一张四联图与两张表，全程无需联网。

## 结果（seed = 42，80 / 20 划分）

| 模型 | Train R² | Test R² | 5 折 CV R² | RMSE |
|---|---|---|---|---|
| Random Forest | 0.90 | **0.45** | 0.49 ± 0.03 | 1.08 |
| Gradient Boosting | 0.72 | 0.45 | 0.45 ± 0.02 | 1.08 |

### 特征重要性：两种口径

| 排名 | 置换重要性（测试集） | 不纯度重要性（RF 内置） |
|---|---|---|
| 1 | MW 0.263 | TPSA 0.122 |
| 2 | RotBonds 0.208 | logP 0.121 |
| 3 | HBA 0.176 | AroRings 0.118 |
| 4 | TPSA 0.139 | MW 0.112 |
| 5 | AroRings 0.127 | qed 0.106 |

两种口径的差异本身是信息：2D 描述符彼此高度相关（MW、HeavyAtoms、RotBonds 都在刻画分子大小与柔性），置换重要性会在相关特征之间分摊权重，因此把更多权重给了 MW 与 RotBonds；而不纯度重要性更平均地铺在 TPSA、logP、AroRings 上。两种口径共同支持的一点是：**极性表面积与芳香环数对活性有实质贡献**。

### 结论

- 12 个 2D 描述符可解释约 45% 的 pIC50 方差，属于 QSAR 中"可用于初筛排序"的水平；
- 训练 R² 0.90 与测试 R² 0.45 的差距说明模型存在过拟合，但 5 折 CV（0.49 ± 0.03）与测试集接近，说明结论不是单次划分的偶然；
- 要进一步提升，需要引入 3D 构象或蛋白-配体相互作用特征，而不是继续堆叠 2D 描述符。

## 数据

`data/egfr_pic50.csv` 共 6,896 条 EGFR 抑制剂记录，两列：`smiles`、`pIC50`（pIC50 = −log10 IC50 [M]）。数据取自 ChEMBL 33；如需重新拉取，可调用 ChEMBL API 按 target（EGFR, CHEMBL203）筛选活性记录后导出为同格式。

## 项目结构

```
.
├── data/egfr_pic50.csv          # 6,896 EGFR inhibitors
├── src/full_qsar_pipeline.py    # CLI: descriptors -> models -> CV -> figures
├── results/
│   ├── qsar_full_analysis.png   # 4-panel figure
│   ├── feature_importance.csv   # permutation + impurity importance
│   └── predictions.csv          # test-set predictions and residuals
├── requirements.txt
└── LICENSE
```

## 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--data` | `../data/egfr_pic50.csv` | 输入 CSV |
| `--out` | `../results` | 输出目录 |
| `--seed` | `42` | 划分与模型随机种子 |
| `--test-size` | `0.2` | 测试集比例 |

## English

A QSAR pipeline for EGFR kinase inhibitors: 6,896 ChEMBL 33 compounds, 12 RDKit 2D descriptors, Random Forest vs Gradient Boosting, 5-fold cross-validation, permutation importance and residual analysis. Run `python src/full_qsar_pipeline.py --data data/egfr_pic50.csv --out results`.

## License

MIT — see [LICENSE](LICENSE).
