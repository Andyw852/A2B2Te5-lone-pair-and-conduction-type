# 带边 s 占比与输运参数的因果检验（纯计算）

只读复用仓库已有数据，不含新的 DFT 计算。报告正文在 Claude Docs 文档《A₂B₂Te₅：带边 s 占比与热电输运参数的因果分析（纯计算）》。

| 脚本 | 检验内容 | 输入 |
| --- | --- | --- |
| `脚本/mechanism_tests.py` | 价带内 A/B/Te-s 权重随能量的分布；价带顶带 E(k) 与 s 权重 w(k) 的相关及两能级拟合 | `07_带边能窗/源数据/*/解析投影.npz`（HSE06+SOC 路径） |
| `脚本/chain_tests.py` | N_V、N_C（300 K）、带边态数、价/导带谷数、PBE VBM 位置；BoltzTraP2 与 AMSET 固定浓度输运；形变势；s 占比窗口 | dense `step3_uniform/EIGENVAL`、`step8.1_boltztrap/boltztrap_crta.json`、`05_SOC_AMSET补算` |
| `脚本/rank_table.py` | 三种 s 定义与能带量/输运量的 Spearman 秩相关（n = 4）；非简并公式预测的 S | `结果/chain.json`、`结果/mech.json` |

运行（仓库根目录为第一个参数）：

```
python 脚本/chain_tests.py <仓库根目录> 结果/chain.json
python 脚本/mechanism_tests.py <仓库根目录> 结果/mech.json
python 脚本/rank_table.py 结果
```

口径：谷数为周期网格上价带顶两组 Kramers 带的局部极大值（与 26 个近邻比较），受 28×28×7（Bi 系 28×28×6）网格分辨率限制。n = 4 时两个无关量 |ρ| ≥ 0.8 的概率为 1/3，秩相关只用于说明“不能据此下结论”。
