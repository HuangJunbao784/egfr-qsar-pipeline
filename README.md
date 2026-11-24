# EGFR QSAR Pipeline

**ChEMBL 33 · 12 个 2D 描述符 + 2048 位 Morgan 指纹 · 随机划分 / 骨架划分 · Random Forest vs Gradient Boosting · 5 折交叉验证**

同一个数据集上比较三套特征（描述符、指纹、两者拼接）与两种划分（随机、Murcko 骨架），并把"随机划分会高估泛化"这件事量化出来。

## 快速开始

```bash
pip install -r requirements.txt
python src/full_qsar_pipeline.py                              # 默认：12 描述符 + 随机划分
python src/full_qsar_pipeline.py --features both --split both # 完整对比
```

或安装为命令行工具：

```bash
pip install -e .
egfr-qsar --features both --split both --out results
egfr-qsar --seeds 42,7,2026 --out results   # 多种子扫描，写入 results/seed_sweep.csv
```

## 结果（Random Forest，seed = 42，80 / 20）

| 特征集 | 随机划分 Test R² | Murcko 骨架划分 Test R² |
|---|---|---|
| 12 个描述符 | 0.451 | 0.289 |
| Morgan 指纹（2048） | 0.707 | 0.560 |
| 描述符 + 指纹 | **0.713** | **0.552** |

Gradient Boosting 在"描述符 + 指纹"上的同一组数字：随机 0.659 / 骨架 0.498。

### 三个结论

1. **指纹是主要杠杆**：同一模型、同一划分，只把 12 个描述符换成 2048 位 Morgan 指纹，随机划分 R² 从 0.451 升到 0.707。二维描述符能表达的信息量就是有限的。
2. **随机划分会高估泛化**：换成骨架划分（训练集与测试集不共享任何 Murcko 骨架，测的是"对新化学型"的外推），同一个模型掉到 0.552–0.560。随机划分下大量近似类似物同时落在训练与测试集，模型记住了骨架而不是规律。
3. **描述符单独跨骨架迁移能力弱**：只有描述符时骨架划分 R² 仅 0.289，与随机划分的 0.451 差距最大。拼接指纹后反而略低于纯指纹（0.552 vs 0.560）——冗余的大小类特征在分布偏移下会拖累模型。

### 数据噪声与天花板

6,896 行记录里只有 56 行是重复 SMILES；有 ≥3 次测量的化合物仅 18 个，而这些化合物的 pIC50 **中位极差达 1.95 个 log 单位**（10% 超过 4.8）。这说明数据集混合了不同 assay 类型与实验条件。测量误差本身有 1–2 个 log 时，随机划分 R² 的物理天花板就在 0.5–0.6 附近。

### 多种子（描述符 + 随机划分，seed = 42 / 7 / 2026）

| 模型 | Test R² (mean ± sd) | CV R² (mean ± sd) |
|---|---|---|
| Random Forest | 0.474 ± 0.010 | 0.484 ± 0.024 |
| Gradient Boosting | 0.436 ± 0.014 | 0.454 ± 0.021 |

### 特征重要性

- **含指纹时**用不纯度重要性（2048 维上做置换重要性不现实）：`morgan_1452` 最高，`logP`、`MW`、`qed` 仍在前 10。
- **仅描述符时**用置换重要性（测试集）：MW、RotBonds、HBA 居前，TPSA 与 AroRings 紧随其后。

## 数据

`data/egfr_pic50.csv` 共 6,896 条 EGFR 抑制剂记录，两列：`smiles`、`pIC50`（pIC50 = −log10 IC50 [M]）。数据取自 ChEMBL 33，可由 `scripts/fetch_chembl_egfr.py` 从 ChEMBL API 重建：

```bash
python scripts/fetch_chembl_egfr.py --out data/egfr_pic50.csv
```

脚本按 target `CHEMBL203`（EGFR）拉取 IC50，只保留单位为 nM、关系为等号、数值为正的记录，换算成 pIC50 后按 SMILES 取中位数去重。若用 conda 的 Python 遇到 `[ASN1: NOT_ENOUGH_DATA]` 握手错误，换用 python.org 安装的（OpenSSL 1.1.1）解释器即可。

## 项目结构

```
.
├── data/egfr_pic50.csv          # 6,896 EGFR inhibitors
├── src/full_qsar_pipeline.py    # CLI: features -> splits -> models -> figures
├── scripts/fetch_chembl_egfr.py # rebuild the dataset from the ChEMBL API
├── tests/                       # unittest suite (17 cases)
├── results/
│   ├── qsar_full_analysis.png   # 4-panel figure
│   ├── split_comparison.csv     # every feature x split combination
│   ├── feature_importance.csv   # permutation (desc) or impurity (fp)
│   ├── predictions.csv          # test-set predictions and residuals
│   └── seed_sweep.csv           # per-seed metrics (--seeds)
├── pyproject.toml               # pip install -e . -> egfr-qsar
├── requirements.txt
└── LICENSE
```

## 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--data` | `../data/egfr_pic50.csv` | 输入 CSV（找不到时回退到 `data/` 与仓库根目录） |
| `--out` | `../results` | 输出目录 |
| `--features` | `desc` | `desc` / `fp` / `both` |
| `--split` | `random` | `random` / `scaffold` / `both` |
| `--seed` | `42` | 图表与表格所用的随机种子 |
| `--seeds` | 无 | 多种子扫描（作用于第一种划分），如 `42,7,2026` |
| `--test-size` | `0.2` | 测试集比例 |

随机划分用 5 折 KFold；骨架划分用 5 折 GroupKFold（按 scaffold 分组），避免同一骨架跨折泄漏。

## 测试

```bash
python -m unittest discover -s tests -v
```

17 个用例覆盖描述符与指纹形状、无效 SMILES 过滤、数据路径解析、随机/骨架划分的无重叠与可复现性、模型拟合与交叉验证、ChEMBL 单位换算；CI 在每次 push 时跑一遍。

## 局限

- 数据混合了不同 assay 类型，尚未按 assay 清洗；这是当前 R² 的主要短板，比换模型更能提升骨架划分下的表现。
- 未使用 3D 构象或蛋白-配体相互作用特征。这类靶点相关特征通常比 2D 指纹更能跨骨架迁移。
- 骨架划分采用"随机分配 scaffold"（seed 固定），不是按字母序排列的极端划分；两种协议给出的数字可以差很多，报告时应写清协议。

## English

QSAR pipeline for EGFR kinase inhibitors: 6,896 ChEMBL 33 compounds, 12 RDKit 2D descriptors and/or 2048-bit Morgan fingerprints, Random Forest vs Gradient Boosting, evaluated under both a random split and a Murcko-scaffold split with grouped cross-validation. The scaffold split is the honest number: fingerprints reach R² ≈ 0.56 on unseen chemotypes versus ≈ 0.71 on a random split.

## License

MIT — see [LICENSE](LICENSE).
