# 实验对照用的 AMSET（重掺杂、300–600 K）

2026-10-31。目的：用 S 反推的载流子浓度，加上 AMSET 迁移率，预测 σ(T) 并与实验比较（论文报告第 5 节“如何进一步与实验对照”）。

现有的密网格结果（`../密网格AMSET/`）只覆盖 ±1e18–1e20，且三种材料只有 300 K。实验样品的反推浓度不在这个范围内：

| 材料 | 反推浓度（HSE / PBE 能带） | 本次 doping（cm⁻³，正为空穴） |
|---|---|---|
| Sn2Sb2Te5 | p 5.6e20 / 2.3e21 | 3.0e20, 6.0e20, 1.2e21, 2.4e21 |
| Sn2Bi2Te5 | p 7.2e20 / 1.2e21 | 3.0e20, 6.0e20, 1.2e21, 2.4e21 |
| Pb2Bi2Te5 | n 8.2e18 / 3.5e19 | −5.0e18, −1.0e19, −2.0e19, −4.0e19 |
| Pb2Sb2Te5 | 近补偿（净 1e16–1e17） | −1.0e17, 1.0e17 |

**温度**：每个材料拆成两个作业，A 为 [300, 400]，B 为 [500, 600]。同一材料的两个作业 doping 列表必须完全相同，合并时按温度拼接。

## 设置：除 doping、temperatures 外与 `step8_amset_full_f1` 完全相同

- **改动范围**：在集群 `~/work/ke_soc_dense_20261019/<材料>/ke-dft-cpu/` 下新建 `step8_amset_exp_A`、`step8_amset_exp_B`。从 `step8_amset_full_f1` 复制（或软链）全部输入：`settings.yaml`、`wavefunction.h5`、`deformation.h5`、能带输入文件等。`settings.yaml` 中只改 `doping` 和 `temperatures` 两项。以下各项都不动：
  - `interpolation_factor: 1`；
  - `scattering_type: [ADP, IMP, POP]`；
  - `bandgap`、介电常数、弹性常数、`pop_frequency`、`soc`。
- **资源**：hanhai25，128 进程，`--time=1-00:00:00`。内存按以前的实测：Sb 系 `--mem=300G`，Bi 系 `--mem=700G`。
- **预计耗时**：每个（浓度, 温度）组合 0.45–0.8 h，每个作业 2–7 h，8 个作业可以并行。
- **提交规矩**（运行记录第 1 节）：submit.sh 先给用户过目再提交。

## 要推送的文件

每个作业的结果放在 `14_散射机制分项迁移率/实验对照AMSET/<材料>/` 下：
- `transport_A.json` 和 `transport_B.json`（原文件名可以保留，前缀加 A/B 区分）；
- `settings_A.yaml`、`settings_B.yaml`；
- `meta.txt`：作业号、墙钟、峰值内存；
- `amset.log` 的最后 50 行。

**不要放进 `../密网格AMSET/`**：两处的 doping 列表不同，合并时会报错。

## 推上来以后

由 Claude 完成：
- 按 log(浓度) 插值出反推浓度处的 σ(T)；
- 同时用 AMSET 自己算的 S 再反推一次浓度，作为自洽对照；
- 与实验 σ(T) 比较，并计算实验的加权迁移率。

前提是用户提供四个样品的 σ(T)（带单位）。
