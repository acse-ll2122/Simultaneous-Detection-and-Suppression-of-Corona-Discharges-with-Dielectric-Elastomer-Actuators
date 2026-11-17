# Simultaneous-Detection-and-Suppression-of-Corona-Discharges-with-Dielectric-Elastomer-Actuators

```text
Simultaneous-Detection-and-Suppression-of-Corona-Discharges-with-Dielectric-Elastomer-Actuators/
├── README.md
├── DEA_only/                    DEA 单独响应，论文 3.1 节
│   ├── code/                    模型、求解器与绘图代码
│   ├── results/
│   │   ├── step_damping/        电压阶跃与阻尼比较，图 2
│   │   └── pulse_response/      电压脉冲与参考半径比较，图 3
│   └── figure/                  图 2、3 的 .fig、.png
├── DEA_full/                    DEA 与电晕耦合模型
│   ├── baseline/                基准响应，论文 3.2 节
│   │   ├── code/                耦合模型、求解器与绘图代码
│   │   ├── results/             Spring 0–5 s、Gas 0–1 s 基准数据
│   │   └── figure/              图 4 的 .fig、.png
│   └── sensitive/               参数影响，论文 3.3 节
│       ├── code/                参数扫描与绘图代码，复用 baseline/code/
│       ├── results/
│       │   └── spring_radius/   弹簧自由半径 2、2.5、3、3.5 mm
│       └── figure/              图 5 的 .fig、.png
└── corona_validation/           电晕参数核验
    ├── code/                    C、n 的拟合、换算与核验代码
    ├── results/                 提取点、拟合表与核验记录
    └── figure/                  核验图 .fig、.png
```

代码为 `.m`，数据为 `.csv`、`.mat`、`.xlsx`，图件为 `.fig`、`.png`；保留模型所需的参数快照 `.json`。图 4 和图 5 共用基准 Spring 数据。文件包不放置 Word 文档和参考论文 PDF。

用 MATLAB R2025b 打开以下各 `code/` 中的 `run_all.m` 重新计算，或打开 `plot_all.m` 重绘已有结果，无需手动添加路径。运行入口会暂时切换到相应代码目录，完成后恢复原目录和路径。

| 内容 | 入口所在目录 |
|---|---|
| 图 2、3：DEA 单独响应 | `DEA_only/code/` |
| 图 4：完整耦合基准响应 | `DEA_full/baseline/code/` |
| 图 5：弹簧自由半径比较 | `DEA_full/sensitive/code/` |
| C、n 拟合与核验 | `corona_validation/code/` |

参数、物理方程和求解器分别在各模型的 `params.m`、`physics.m`、`solver.m`；参数扫描复用 `DEA_full/baseline/code/` 中的模型与参数。完整重算顺序为 DEA_only → baseline → sensitive → corona_validation，重算会覆盖对应输出。单独运行 sensitive 时，如缺少基准数据会先计算 baseline；若基准参数已改变，需先重算 baseline。重绘仅更新图件，不求解或覆盖数据。

图 5 通过 `DEA_full/sensitive/results/spring_radius/shared_baseline.csv` 引用图 4 的基准 Spring 数据；其 MAT 文件同样使用包内相对路径读取基准轨迹。C、n 核验使用代码内置的文献图上提取点，不需要 PDF，也不依赖 Python；提取点不是原始实验数据。

目录依据论文 `111626_E4` 版本整理。论文表 5 的电导率比较及正文中的耦合阻尼扫描，原包未找到独立数据与运行入口，未列为已备齐结果。

已在 MATLAB R2025b 完整重算 21 个工况并核对保存数据；四张结果图的 59 条曲线数据与线型均与原 FIG 一致。图 2–5 的重绘 PNG 与原图逐像素一致；图 2、3、4 的标题字重及图 3 的子图位置已在绘图代码中固定。
