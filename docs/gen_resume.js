const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        HeadingLevel, AlignmentType, BorderStyle, WidthType, ShadingType } = require("docx");

const BLUE = "1F4E79";
const GRAY = "666666";
const LIGHT_BLUE = "E8F0F8";

const border = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: border, bottom: border, left: border, right: border };
const cellMargins = { top: 60, bottom: 60, left: 100, right: 100 };

function section(title, children) {
    return [
        new Paragraph({
            spacing: { before: 260, after: 100 },
            border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: BLUE, space: 4 } },
            children: [new TextRun({ text: title, bold: true, size: 28, font: "Microsoft YaHei", color: BLUE })]
        }),
        ...children
    ];
}

function projectBlock(title, date, items) {
    const paras = [
        new Paragraph({
            spacing: { before: 140, after: 40 },
            children: [
                new TextRun({ text: title, bold: true, size: 24, font: "Microsoft YaHei" }),
                new TextRun({ text: "  |  " + date, size: 20, font: "Microsoft YaHei", color: GRAY }),
            ]
        })
    ];
    for (const item of items) {
        paras.push(new Paragraph({
            spacing: { before: 20, after: 20 },
            indent: { left: 360 },
            children: [new TextRun({ text: "• " + item, size: 21, font: "Microsoft YaHei" })]
        }));
    }
    return paras;
}

function skillRow(label, content) {
    return new TableRow({
        children: [
            new TableCell({
                borders, shading: { fill: LIGHT_BLUE, type: ShadingType.CLEAR }, margins: cellMargins,
                width: { size: 1800, type: WidthType.DXA },
                children: [new Paragraph({ children: [new TextRun({ text: label, bold: true, size: 20, font: "Microsoft YaHei" })] })]
            }),
            new TableCell({
                borders, margins: cellMargins,
                width: { size: 7200, type: WidthType.DXA },
                children: [new Paragraph({ children: [new TextRun({ text: content, size: 20, font: "Microsoft YaHei" })] })]
            })
        ]
    });
}

const doc = new Document({
    styles: {
        default: { document: { run: { font: "Microsoft YaHei", size: 21 } } }
    },
    sections: [{
        properties: {
            page: {
                size: { width: 11906, height: 16838 },
                margin: { top: 1000, right: 1300, bottom: 1000, left: 1300 }
            }
        },
        children: [
            // Header
            new Paragraph({ spacing: { after: 40 }, children: [new TextRun({ text: "黄君宝", bold: true, size: 44, font: "Microsoft YaHei", color: BLUE })] }),
            new Paragraph({ spacing: { after: 40 }, children: [new TextRun({ text: "吉林大学  |  药物化学硕士  |  计算生物学与AI药物化学方向", size: 22, font: "Microsoft YaHei", color: GRAY })] }),
            new Paragraph({ spacing: { after: 200 }, children: [new TextRun({ text: "15879344038@163.com  |  github.com/HuangJunbao784", size: 20, font: "Microsoft YaHei", color: GRAY })] }),

            // 教育背景
            ...section("教育背景", [
                new Paragraph({ spacing: { before: 60, after: 30 }, children: [
                    new TextRun({ text: "吉林大学", bold: true, size: 22, font: "Microsoft YaHei" }),
                    new TextRun({ text: "（985）— 药物化学硕士", size: 21, font: "Microsoft YaHei" })
                ]}),
                new Paragraph({ indent: { left: 360 }, spacing: { before: 20, after: 20 }, children: [
                    new TextRun({ text: "研究方向：计算生物学与AI药物化学    |    2024.09 - 2027.06（预计）", size: 20, font: "Microsoft YaHei", color: GRAY })
                ]}),
                new Paragraph({ spacing: { before: 60, after: 30 }, children: [
                    new TextRun({ text: "江西中医药大学", bold: true, size: 22, font: "Microsoft YaHei" }),
                    new TextRun({ text: " — 药学学士", size: 21, font: "Microsoft YaHei" })
                ]}),
                new Paragraph({ indent: { left: 360 }, spacing: { before: 20, after: 20 }, children: [
                    new TextRun({ text: "系统学习药学与化学核心课程，掌握常规药物分析、药理学、药剂学实验操作    |    2020.09 - 2024.06", size: 20, font: "Microsoft YaHei", color: GRAY })
                ]}),
            ]),

            // 项目经历
            ...section("项目经历", [
                ...projectBlock("硕士课题：格列本脲CYP2C9/CYP3A4选择性代谢的分子机制：MD与QM联合研究\nMolecular Mechanism of Isoform-Specific Regioselective Metabolism of Glibenclamide by CYP2C9 and CYP3A4: A Combined MD and QM Study", "2025.03 - 至今", [
                    "采用分子对接-500ns全原子MD模拟-MM/GBSA自由能计算-QM反应能垒计算的多尺度策略，系统解析两种CYP同工酶对格列本脲代谢选择性差异的分子机制",
                    "完成3个体系共1.5μs分子动力学模拟，通过NAC统计与距离监测定量分析反应可及性，明确CYP2C9优势代谢位点为环己基C4、CYP3A4为乙基桥C2'",
                    "残基水平能量分解鉴定关键作用残基：CYP2C9三苯丙氨酸簇（Phe100/114/476）、CYP3A4五苯丙氨酸簇（Phe108/213/215/220/304），揭示π-π堆积与疏水作用的稳定机制",
                    "基于Gaussian 16完成4个反应位点氢提取能垒计算，构建\"热力学稳定性→反应可及性→电子反应性\"完整证据链，首次提出四因素协同解释模型",
                    "论文已完成，目标投稿JCIM / PCCP等计算化学领域期刊"
                ]),
                ...projectBlock("EGFR激酶抑制剂QSAR建模  |  独立项目", "2026.03 - 2026.06", [
                    "从ChEMBL 33获取6,896个EGFR激酶抑制剂的pIC50活性数据",
                    "使用RDKit提取12个2D分子描述符（MW, logP, TPSA, HBD, HBA等）构建特征矩阵",
                    "对比Random Forest与Gradient Boosting，最优RF达测试集R²=0.45、5-fold CV R²=0.49",
                    "特征重要性分析发现TPSA、logP和芳香环数为Top 3，与EGFR激酶口袋理化性质一致",
                    "完整pipeline自动化：数据清洗→特征提取→模型训练→交叉验证→可视化"
                ]),
                ...projectBlock("MD轨迹自动化分析工具  |  独立项目", "2026.05 - 2026.06", [
                    "基于Python + MDTraj开发AMBER/GROMACS轨迹自动化分析脚本",
                    "实现RMSD/RMSF/Rg/PCA本质动力学/自由能景观/配体-蛋白氢键分析",
                    "六合一分析图自动生成，替代GROMACS命令行+VMD手动操作流程",
                    "应用于500ns蛋白-配体复合物轨迹，识别与配体形成氢键占有率>50%的关键残基"
                ]),
            ]),

            // 核心技能
            ...section("核心技能", [
                new Table({
                    width: { size: 9000, type: WidthType.DXA },
                    columnWidths: [1800, 7200],
                    rows: [
                        skillRow("药学与实验技能", "熟练掌握四大药学（药理/药分/药化/药剂）、天然药化；熟悉有机合成、药物分析、药剂学、药理学等实验操作；能独立完成药理活性筛选、药物含量测定、制剂制备等常规实验"),
                        skillRow("化学基础", "系统掌握有机化学、无机化学、分析化学、物理化学、生物化学核心理论，具备扎实的化学键理论、热力学、动力学、光谱分析等基础知识"),
                        skillRow("计算化学", "分子对接（AutoDock Vina/Glide）、MD模拟（GROMACS/AMBER）、MD轨迹分析、量子力学计算（Gaussian）、自由能计算"),
                        skillRow("AI / ML", "QSAR建模、随机森林、梯度提升、交叉验证、过拟合诊断、特征工程"),
                        skillRow("深度学习", "PyTorch、MLP、训练循环"),
                        skillRow("编程与工具", "Python（pandas, matplotlib, RDKit, MDTraj, scikit-learn）、Git/GitHub、Jupyter"),
                    ]
                })
            ]),

            // 自我评价
            ...section("自我评价", [
                new Paragraph({ spacing: { before: 60, after: 60 }, children: [new TextRun({ text: "药物化学硕士，具备完整的计算化学研究能力，熟练掌握分子对接、分子动力学模拟、自由能计算与量子化学计算等多尺度模拟方法。自学Python与机器学习，能独立完成从数据获取、特征工程到模型训练评估的完整QSAR建模管线。擅长将传统CADD的物理化学理解与机器学习方法结合，开展有机制深度的计算研究。GitHub开源多个计算化学分析工具。", size: 21, font: "Microsoft YaHei" })] })
            ]),

            // 实习意向
            ...section("实习意向", [
                new Paragraph({ spacing: { before: 60, after: 60 }, children: [new TextRun({ text: "暑期实习  |  AI制药 / 计算化学 / CADD方向  |  2026.07 - 2026.09", size: 21, font: "Microsoft YaHei" })] })
            ]),
        ]
    }]
});

Packer.toBuffer(doc).then(buffer => {
    fs.writeFileSync("c:/Users/USER/Desktop/vscode claude code/interview-project/docs/黄君宝_简历.docx", buffer);
    console.log("Done: 黄君宝_简历.docx");
});
