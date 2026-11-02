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
- 与实验的对照，用户已确定（2026-11-02）：在室温霍尔浓度下直接算 S(T)，不调任何参数。

## 输入（集群上的 8 个作业目录，都是同一 12×12×2 网格）

| 材料 | HSE06+SOC（主结果） | PBE+SOC 同网格（检验网格够不够） | NELECT |
|---|---|---|---|
| Pb₂Sb₂Te₅ | clus_xa `~/autozt_work/1_0/Pb2Sb2Te5/c1_scf/`（vasprun.xml 14 MB） | jzzn `~/ke_1_0_bandgap/Pb2Sb2Te5/b_pbe_parity2/` | 68 |
| Pb₂Bi₂Te₅ | hanhai 1-0 的 `c_hse_A`（作业 306976） | 1-0 的 `a_pbe`（作业 305553–305555 之一） | 88 |
| Sn₂Sb₂Te₅ | hanhai `c_hse_A`（作业 306978） | `a_pbe` | 68 |
| Sn₂Bi₂Te₅ | hanhai `c_hse_A`（作业 307002） | `a_pbe` | 88 |

每个目录先核对四项：
- vasprun.xml 存在；
- KPOINTS 是 Gamma 12 12 2；
- INCAR 有 `LSORBIT = .TRUE.`；HSE 目录有 `LHFCALC = .TRUE.`，PBE 目录没有；
- OUTCAR 已收敛。

## 运行

脚本 `10_脚本/run_bt2_hse_grid.py`。它不改 AutoZT，只调用 AutoZT 的 `gen_step11_boltztrap.run_boltztrap_crta`（与密网格结果同一函数），设置如下：
- 温度 300–650 K，每 25 K；
- 化学势 ±1.2 eV、600 点，DOS 6000 点；
- HSE 不加剪刀；PBE 剪刀与密网格相同（PbBi 0.1141、PbSb 0.0654、SnBi 0.1366、SnSb 0.1352 eV），这样它与密网格只差 k 网格；
- NELECT 不符时报错退出。

每个目录跑三档插值倍数 `--mult 5 20 60`（38 个点很少，要看结果是否随倍数收敛），4 × 2 × 3 = 24 次。例：

```
python 10_脚本/run_bt2_hse_grid.py --autozt <AutoZT 的 step8.1 目录> --material Sn2Sb2Te5 --functional hse \
    --run-dir <c_hse_A 目录> --mult 20 --out 22_能带基准_1-0/Sn2Sb2Te5/c1_HSE_SOC/bt2_mult20 \
    --kcheck 22_能带基准_1-0/Sn2Sb2Te5/bt2_check_kpoints.txt
```

`--kcheck` 让脚本另外输出插值能带在检验 k 点上的能量（`bands_check.txt`）：
- 0–203 号是 c2 的 HSE 路径；
- Sn 系另有 c5 镜面内的两个真实价带谷。

## 结果放哪

```
22_能带基准_1-0/<材料>/c1_HSE_SOC/bt2_mult{5,20,60}/{boltztrap_crta.json, bands_check.txt}
22_能带基准_1-0/<材料>/b_PBE_SOC/bt2_mult{5,20,60}/{boltztrap_crta.json, bands_check.txt}
```

Pb₂Sb₂Te₅ 没有材料子目录，直接放在 `22_能带基准_1-0/c1_HSE_SOC/`、`22_能带基准_1-0/b_PBE_SOC/` 下。每个 JSON 约 0.5 MB。

## 验收（Claude 做）

1. **倍数收敛**：mult 20 与 60 的 S（300 K，±10²⁰、±10²¹ cm⁻³）相差 < 5%。
2. **网格够不够**：比较 PBE 12×12×2 与 PBE 28×28×7 密网格（`A2B2Te5_SOC_dense/.../step8.1_boltztrap/boltztrap_crta.json`），看 300 K、p 型与 n 型、3×10¹⁹–2×10²¹ cm⁻³ 的 S。
   - 相差 ≤ 10%：12×12×2 够用，HSE 结果直接进报告。
   - 相差更大：需要更密的 HSE 网格（非自洽 HSE，费时）或 Wannier 插值，另议。
3. **插值能带**：`bands_check.txt` 与直接计算对照。路径对 `c2_hse/path_eigen.txt`；Sn 系镜面内的真实价带谷对 c5 的结果，看插值有没有把价带顶放对。

## 之后报告怎么改

- 所有 HSE 数值改为 HSE BoltzTraP2 的直接结果，删去比值法：图 8、表 6、图 9、表 7、图 S3、表 A1。
- 图 S3 改为在室温霍尔浓度下直接计算的 S(T)，没有可调参数。表 A1 改为霍尔浓度下的计算值与偏差；“需要多少浓度才能对上”只在正文提一句作诊断。
- 多谷模型只保留定性用途：说明哪些谷在导电（孤对电子谷的电导占比）。
