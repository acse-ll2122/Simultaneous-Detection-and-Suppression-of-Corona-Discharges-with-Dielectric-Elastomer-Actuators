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
