# EGFR QSAR Pipeline

**ChEMBL 33 · RDKit 2D descriptors · Random Forest / Gradient Boosting · 5-fold cross-validation**

从 SMILES 出发，用 12 个二维分子描述符预测 EGFR 激酶抑制剂的 pIC50，比较两种集成模型，并给出交叉验证、置换重要性、残差与特征排序。

## 快速开始

```bash
pip install -r requirements.txt
python src/full_qsar_pipeline.py --data data/egfr_pic50.csv --out results
```

也可以安装为命令行工具：

```bash
pip install -e .
egfr-qsar --data data/egfr_pic50.csv --out results
egfr-qsar --seeds 42,7,2026 --out results   # 多种子评估，写入 results/seed_sweep.csv
```

脚本会在 `results/` 生成一张四联图与若干表格，训练与评估本身无需联网。

## 结果（seed = 42，80 / 20 划分）

| 模型 | Train R² | Test R² | 5 折 CV R² | RMSE |
|---|---|---|---|---|
| Random Forest | 0.90 | **0.45** | 0.49 ± 0.03 | 1.08 |
| Gradient Boosting | 0.72 | 0.45 | 0.45 ± 0.02 | 1.08 |

### 多种子稳定性（seed = 42 / 7 / 2026）

| 模型 | Test R² (mean ± sd) | CV R² (mean ± sd) | RMSE (mean) |
|---|---|---|---|
| Random Forest | 0.476 ± 0.029 | 0.484 ± 0.024 | 1.049 |
| Gradient Boosting | 0.459 ± 0.019 | 0.454 ± 0.021 | 1.066 |

单次划分的 Test R² 会在 0.45 附近浮动约 ±0.03，CV 的波动更小。逐 seed 明细见 `results/seed_sweep.csv`。

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

`data/egfr_pic50.csv` 共 6,896 条 EGFR 抑制剂记录，两列：`smiles`、`pIC50`（pIC50 = −log10 IC50 [M]）。数据取自 ChEMBL 33，由 `scripts/fetch_chembl_egfr.py` 从 ChEMBL API 重建：

```bash
python scripts/fetch_chembl_egfr.py --out data/egfr_pic50.csv
```

脚本按 target `CHEMBL203`（EGFR）拉取 IC50、只保留单位为 nM 且关系为等号、数值为正的记录，换算成 pIC50 后按 SMILES 取中位数去重。若用 conda 的 Python 遇到 `[ASN1: NOT_ENOUGH_DATA]` 握手错误，换用 python.org 安装的（OpenSSL 1.1.1）解释器即可。

## 项目结构

```
.
├── data/egfr_pic50.csv          # 6,896 EGFR inhibitors
├── src/full_qsar_pipeline.py    # CLI: descriptors -> models -> CV -> figures
├── scripts/fetch_chembl_egfr.py # rebuild the dataset from the ChEMBL API
├── tests/                       # unittest suite (11 cases)
├── results/
│   ├── qsar_full_analysis.png   # 4-panel figure
│   ├── feature_importance.csv   # permutation + impurity importance
│   ├── predictions.csv          # test-set predictions and residuals
│   └── seed_sweep.csv           # per-seed metrics (--seeds)
├── pyproject.toml               # pip install -e . -> egfr-qsar
├── requirements.txt
└── LICENSE
```

## 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--data` | `../data/egfr_pic50.csv` | 输入 CSV（找不到时会依次回退到 `data/` 与仓库根目录） |
| `--out` | `../results` | 输出目录 |
| `--seed` | `42` | 图表与表格所用的随机种子 |
| `--seeds` | 无 | 多种子评估，如 `42,7,2026`，结果写入 `seed_sweep.csv` |
| `--test-size` | `0.2` | 测试集比例 |

## 测试

```bash
python -m unittest discover -s tests -v
```

11 个用例覆盖描述符计算、无效 SMILES 过滤、数据路径解析、模型拟合与交叉验证、ChEMBL 单位换算；CI 会在每次 push 时跑一遍。

## English

A QSAR pipeline for EGFR kinase inhibitors: 6,896 ChEMBL 33 compounds, 12 RDKit 2D descriptors, Random Forest vs Gradient Boosting, 5-fold cross-validation, permutation importance and residual analysis. Run `python src/full_qsar_pipeline.py --data data/egfr_pic50.csv --out results`.

## License

MIT — see [LICENSE](LICENSE).
