# 简历 — 项目经历部分

## 教育背景

**吉林大学** | 药物化学硕士 | 2024 - 至今（预计2027毕业）
研究方向：计算生物学与AI药物化学
GPA: [你填]

**吉林大学** | [本科专业]学士 | [年份]

---

## 项目经历

### EGFR激酶抑制剂QSAR建模 | 个人项目
[GitHub: github.com/HuangJunbao784/egfr-qsar-pipeline]

- 从ChEMBL 33数据库获取6,896个EGFR抑制剂的生物活性数据（pIC50）
- 使用RDKit提取12个2D分子描述符（MW, logP, TPSA, 氢键供体/受体等）构建特征矩阵
- 对比Random Forest和Gradient Boosting两个模型，Random Forest达到测试集R² = 0.45, 5折交叉验证R² = 0.49
- 特征重要性分析显示TPSA、logP和芳香环数为Top 3特征，与EGFR激酶口袋的理化性质一致
- 完整pipeline自动化：数据清洗 → 特征提取 → 模型训练 → 交叉验证 → 可视化出图

### MD轨迹自动化分析工具 | 个人项目
[GitHub: github.com/HuangJunbao784/md-analysis-tool]

- 基于Python + MDTraj开发AMBER/GROMACS轨迹自动化分析脚本
- 实现RMSD/RMSF/Rg/PCA本质动力学/自由能景观/配体-蛋白氢键分析
- 6合1分析图自动生成，替代GROMACS命令行 + VMD手动操作流程
- 应用于500ns蛋白-配体复合物轨迹，识别出与配体形成氢键占有率>50%的关键残基

### [你的论文题目] | 硕士课题
[JCIM / JMC / 等目标期刊]

- [一句话：你用了什么方法，解决了什么问题]
- 整合分子对接、MD模拟、QM计算建立完整CADD分析管线
- [接收后补充发表信息]

---

## 技能

| 类别 | 具体技能 |
|------|---------|
| 计算化学 | 分子对接（AutoDock Vina / Glide）、MD模拟（GROMACS / AMBER）、MD轨迹分析、量子力学计算（Gaussian）、自由能计算 |
| 编程 & 数据 | Python（pandas, matplotlib, RDKit, MDTraj）、scikit-learn、PyTorch、Git |
| AI/ML | QSAR建模、随机森林、梯度提升、交叉验证、特征工程、过拟合诊断 |
| 工具 | Jupyter Notebook、VS Code、Git/GitHub |

---

## 放在简历上的 GitHub 链接

```
github.com/HuangJunbao784/egfr-qsar-pipeline
github.com/HuangJunbao784/md-analysis-tool
```

---

## 简历一句话总结（放在简历最上面）

> 药物化学硕士，精通分子对接、MD模拟和QM计算，自学Python和机器学习，独立完成6,896个EGFR抑制剂的QSAR建模（测试R²=0.45）和500ns MD轨迹自动化分析工具开发。代码在GitHub开源。
