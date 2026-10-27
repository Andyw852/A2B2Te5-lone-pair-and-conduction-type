# Pb₂Bi₂Te₅ · c2_hse（HSE06+SOC 能带路径）

2026-10-27 完成。

## 1. 计算设置

| 项 | 值 |
|---|---|
| 体系 | Pb₂Bi₂Te₅，9 原子，P-3m1 |
| 泛函 | HSE06（`LHFCALC=.TRUE.`） |
| SOC | 是（`LSORBIT=.TRUE.`，`vasp_ncl`） |
| **NELECT** | **88**（VBM = 第 88 带，CBM = 第 89 带；Bi 5d¹⁰ 半芯贡献 20 个电子） |
| ENCUT | 360 eV |
| k 网格（SCF） | 12×12×2，38 个不可约点 |
| **路径** | `KPOINTS_OPT`，**204 个点**（与另三种材料同一条路径） |
| 起点 | `c_hse_A` 的 HSE WAVECAR（`ISTART=1`） |

## 2. 带边结果

| 量 | 值 |
|---|---|
| VBM | 5.7296 eV @ k=97（坐标 (0, 0, −0.5) = **A 点**） |
| CBM | 5.8715 eV @ k=0（**Γ**） |
| **路径带隙** | **0.1419 eV** |
| `E-fermi` | 5.8006 eV |

**VBM 在 A 点、CBM 在 Γ —— 间接带隙沿 Γ–A。**与 SnSb 不同（SnSb 的 CBM 在 Γ–A 线上、VBM 在 Γ）。

## 3. 与 SCF 网格（`c_hse_A`，12×12×2）的交叉验证

| 量 | `c_hse_A` | `c2_hse` 路径 | 差 |
|---|---:|---:|---:|
| VBM | 5.729552 | 5.729600 | <1×10⁻⁵ eV |
| CBM | 5.871537 | 5.871500 | <4×10⁻⁵ eV |
| 带隙 | 0.1420 | **0.1419** | <1×10⁻⁴ eV |

**Γ 点 band 88/89 与 SCF 网格吻合到 10⁻⁵ eV 量级。**
（Γ 点全部 128 带的逐带最大差较大，来自 Bi 5d¹⁰ 深半芯带，不影响带边。）

## 4. 文件说明

| 文件 | 内容 |
|---|---|
| `path_eigen.txt` | 204 个 k 点 × 128 条带，三列 `k_index band_index energy_eV` |
| `KPOINTS_OPT` | 204 点路径定义 |
| `INCAR` / `KPOINTS` | SCF 设置 |
| `PROCAR_OPT.gz` | 路径上的投影（fat-band / 成分分析用） |

## 5. 提取方法

```
python ext_opt2.py vasprun.xml 204 path_eigen.txt
```

本征值在 `vasprun.xml` 的 `<eigenvalues_kpoints_opt>` 段，**不在** `EIGENVAL` 里。