# HSE06+SOC 直接 BoltzTraP2（密网格 28×28×7/6）结果汇总

2026-11-03。四材料统一用 **A2B2Te5_SOC_dense 的密网格 step3_uniform 的 HSE06+SOC 能带**直接跑 BoltzTraP2（CRTA），不再用 PBE 比值法。

## 输入

- 能带：`A2B2Te5_SOC_dense/materials/<材料>/ke-dft-cpu/result/step3_uniform/vasprun.xml`（HSE06+SOC，PbSb/SnSb 为 28×28×7，PbBi/SnBi 为 28×28×6）
- 环境：hanhai25，Python 3.11.16 + numpy 1.26.4 + BoltzTraP2 24.1.1（`$HOME/venvs/bt2py311n1`）
- 脚本：`10_脚本/run_bt2_hse_grid.py`，调用 AutoZT `gen_step11_boltztrap.run_boltztrap_crta`
- 设置：温度 300–650 K 每 25 K；化学势 ±1.2 eV、600 点；DOS 6000 点；HSE **不加剪刀**

## 倍数收敛性（面内 S/µV·K⁻¹，300 K）

| 材料 | 浓度 cm⁻³ | mult=5 | mult=20 | mult=60 | |20−60|% |
|---|---|---:|---:|---:|---:|
| Pb₂Sb₂Te₅ | 3×10¹⁹ | 218.9 | 219.2 | 219.1 | 0.04 |
| Pb₂Sb₂Te₅ | 1×10²⁰ | 118.5 | 118.5 | 118.4 | 0.04 |
| Pb₂Sb₂Te₅ | 1×10²¹ | 27.1 | 27.1 | 27.1 | 0.02 |
| Pb₂Bi₂Te₅ | 3×10¹⁹ | 210.9 | 211.3 | 211.3 | 0.04 |
| Pb₂Bi₂Te₅ | 1×10²⁰ | 122.1 | 122.3 | 122.4 | 0.05 |
| Pb₂Bi₂Te₅ | 1×10²¹ | 31.5 | 31.6 | 31.6 | 0.06 |
| Sn₂Sb₂Te₅ | 3×10¹⁹ | 198.3 | 198.4 | 198.4 | 0.01 |
| Sn₂Sb₂Te₅ | 1×10²⁰ | 116.4 | 116.4 | 116.5 | 0.00 |
| Sn₂Sb₂Te₅ | 1×10²¹ | 28.3 | 28.2 | 28.2 | 0.03 |
| Sn₂Bi₂Te₅ | 3×10¹⁹ | 196.0 | 196.1 | 196.2 | 0.01 |
| Sn₂Bi₂Te₅ | 1×10²⁰ | 119.8 | 120.0 | 119.9 | 0.03 |
| Sn₂Bi₂Te₅ | 1×10²¹ | 35.8 | 35.7 | 35.7 | 0.05 |

**结论：|mult20 − mult60| 全部 < 0.06%，远低于 5% 阈值，插值已收敛。**

## 耗时（hanhai25）

| 材料 | mult=5 | mult=20 | mult=60（计算节点，~115 GB 峰值内存） |
|---|---:|---:|---:|
| Pb₂Sb₂Te₅ | 267 s | 650 s | ~35 min |
| Pb₂Bi₂Te₅ | 241 s | 532 s | ~35 min |
| Sn₂Sb₂Te₅ | 248 s | 657 s | ~32 min |
| Sn₂Bi₂Te₅ | 236 s | 543 s | ~28 min |

## 提醒（给 Claude 接图时）

- 这是**密网格 HSE** 结果，不是交接里说的 12×12×2。因为 PbSb 的 12×12×2 HSE 在已到期的 clus_xa 上、本地未拉回，而四材料密网格 vasprun 均在 `_归档_20261028...` 中，故统一用密网格（更密、更可靠）。
- `[WARN] 未做剪刀修正` 是预期行为：HSE 本就不加剪刀，脚本只是提示没取到 PBE 目标带隙。
- 每个 `boltztrap_crta.json` 的 `grid` 字段写了 12×12×2 是脚本默认文案，实际为密网格（见 `input_dir`/材料）。