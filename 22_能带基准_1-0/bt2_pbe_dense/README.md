# PBE+SOC 密网格 BoltzTraP2（不加剪刀）

四种材料各一个 `boltztrap_crta.json`。

| 项 | 内容 |
|---|---|
| 能带 | PBE+SOC，`A2B2Te5_SOC_dense` 的 step3_uniform 密网格：Sb 系 28×28×7，Bi 系 28×28×6，全网格（ISYM = −1） |
| 剪刀 | 不加，用 PBE+SOC 本身的带隙 |
| 输运 | BoltzTraP2 常数弛豫时间（AutoZT `run_boltztrap_crta`），插值倍数 60 |
| 温度 | 300–650 K，每 25 K |
| 面内 S | 取 xx、yy |
| 载流子浓度 | + 为空穴、− 为电子 |

**来源**：分支 `bt2-hse-grid`（提交 9204897）。那批结果当时误标为 HSE06+SOC。核对后确认：输入 vasprun.xml 的 SHA-256 与 PBE 密网格逐字节相同，所以它实际是 PBE+SOC、不加剪刀。这里改正了 `functional` 字段，并在 `provenance` 中写明来历。

**插值倍数**：mult 20 与 60 相比，在 ±3×10¹⁹–2×10²¹ cm⁻³、300–600 K 范围内最多差 0.22 µV/K。

**与加剪刀的密网格结果（AutoZT 原结果，11-04 以前报告里的 PBE）比较**：在四个样品的霍尔浓度下几乎相同。只有 Pb₂Bi₂Te₅ 在 500 K 以上，|S| 比加剪刀的高 2–6 µV/K（带隙 0.143 eV 对 0.114 eV，双极弱一些）。

2026-11-04 起报告的输运全部用这批数据（`10_脚本/paper_figures.py` 的 `pbe_dense()`）：第 4 节的图 8、表 6，第 5 节的图 9、表 7、表 8。
