# 交接：直接用 HSE 能带跑 BoltzTraP2（2026-11-02）

给 AutoZT 的 BoltzTraP2 步骤（step8.1，`run_boltztrap_crta`）。只做后处理，不需要新的 VASP 计算。

## 为什么要跑

- 报告现在的“HSE 估计”是 PBE+SOC BoltzTraP2 的 S 乘以多谷模型的 HSE/PBE 比值。文献里没有这种做法，用户要求改为直接用 HSE。
- 多谷模型单独也不能当定量结果。在同一套 PBE 能带上和 BoltzTraP2 全区比较，它在高浓度下系统性偏高。下表是 300 K、p 型的 S（µV/K），左为 BoltzTraP2，右为多谷模型。n 型也大多偏高（只有 Sn₂Sb₂Te₅ 在 1×10²¹ 时接近）。

  | 浓度（cm⁻³） | Pb₂Sb₂Te₅ | Pb₂Bi₂Te₅ | Sn₂Sb₂Te₅ | Sn₂Bi₂Te₅ |
  |---|---|---|---|---|
  | 1×10²⁰ | 120 / 193 | 121 / 210 | 113 / 157 | 122 / 193 |
  | 1×10²¹ | 27 / 79 | 32 / 102 | 28 / 59 | 35 / 81 |

  原因：模型是椭球谷近似，只取一条价带、一条导带。约 10²¹ cm⁻³ 时费米能级深入带内 400–500 meV，这些近似都不成立。样品浓度是 8×10¹⁹–1.7×10²¹ cm⁻³，正好落在这个范围里。
- 已有的全区 HSE 数据，就是 1-0 基准的 HSE06+SOC 自洽计算：Γ 中心 12×12×2 网格，38 个不可约点，四种材料都有。文献的常规做法是用 BoltzTraP2 直接插值这些 HSE 本征值，不加剪刀修正。

## 要跑的作业（4 种材料 × 2）

| 材料 | HSE06+SOC（主结果） | PBE+SOC 同网格（检验网格够不够） |
|---|---|---|
| Pb₂Sb₂Te₅ | clus_xa `~/autozt_work/1_0/Pb2Sb2Te5/c1_scf/`（vasprun.xml 14 MB） | jzzn `~/ke_1_0_bandgap/Pb2Sb2Te5/b_pbe_parity2/` |
| Pb₂Bi₂Te₅ | hanhai 1-0 的 `c_hse_A`（作业 306976） | 1-0 的 `a_pbe`（作业 305553–305555 之一） |
| Sn₂Sb₂Te₅ | hanhai `c_hse_A`（作业 306978） | `a_pbe` |
| Sn₂Bi₂Te₅ | hanhai `c_hse_A`（作业 307002） | `a_pbe` |

八个目录都是同一 12×12×2 网格。`run_boltztrap_crta` 读目录里的 vasprun.xml。如果某个目录没有 vasprun.xml、只有 EIGENVAL，告诉我，另想办法。

## 设置

- 与现有密网格结果（`A2B2Te5_SOC_dense/materials/<材料>/ke-dft-cpu/result/step8.1_boltztrap/boltztrap_crta.json`）用同一函数、同一输出格式：CRTA，面内取 (xx+yy)/2，带符号载流子浓度（+ 为空穴）。
- **HSE**：`scissor_gap_eV=None`，用 HSE 自己的带隙。
- **PBE 同网格**：剪刀与密网格一致（`run_bt2_soc.py` 的 GAPS：PbBi 0.1141、PbSb 0.0654、SnBi 0.1366、SnSb 0.1352 eV）。这样它与密网格结果只差 k 网格。
- **温度**：300–650 K，每 25 K（实验点在 300–620 K；密网格只算了 300、500、700 K）。
- **插值倍数 `kmesh_mult`**：38 个点很少，密网格用的 2 不够。HSE、PBE 都跑 20 和 60 两档。300 K 下两档的 S 相差 < 5% 才算收敛，取 60 档。
- **化学势**：`mu_window_ev=2.0`，`mu_npts` ≥ 1000，`dos_npts` 3000。要能覆盖 ±10¹⁸ 到 ±3×10²¹ cm⁻³。

## 结果放哪

```
22_能带基准_1-0/<材料>/c1_HSE_SOC/bt2_mult20/boltztrap_crta.json
22_能带基准_1-0/<材料>/c1_HSE_SOC/bt2_mult60/boltztrap_crta.json
22_能带基准_1-0/<材料>/b_PBE_SOC/bt2_mult20/boltztrap_crta.json
22_能带基准_1-0/<材料>/b_PBE_SOC/bt2_mult60/boltztrap_crta.json
```

Pb₂Sb₂Te₅ 没有材料子目录，直接放在 `22_能带基准_1-0/c1_HSE_SOC/`、`22_能带基准_1-0/b_PBE_SOC/` 下面。

## 验收（我来做）

1. **网格够不够**：比较 PBE 12×12×2 与 PBE 28×28×7（密网格），看 300 K、p 型与 n 型、3×10¹⁹–2×10²¹ cm⁻³ 的 S。
   - 相差 ≤ 10%：12×12×2 够用，HSE 结果直接进报告。
   - 相差更大：需要更密的 HSE 网格（非自洽 HSE，费时）或 Wannier 插值，另议。
2. **价带顶是否插值出来**：价带顶在镜面内的一般位置（c4/c5 真实谷），不在网格点上。用插值后的 HSE 价带顶能量，对照直接计算的真实谷能量。如果 AutoZT 能输出插值能带（或 .bt2 文件），一并带回。

## 之后报告怎么改

- 所有 HSE 数值改为 HSE BoltzTraP2 的直接结果，删去比值法：图 8、表 6、图 9、表 7、图 S3、表 A1。
- 多谷模型只保留定性用途：说明哪些谷在导电（孤对电子谷的电导占比）。
- S(T) 对照建议改为在霍尔浓度下直接计算，没有可调参数；“需要多少浓度才能对上”只作诊断。待用户确认。
