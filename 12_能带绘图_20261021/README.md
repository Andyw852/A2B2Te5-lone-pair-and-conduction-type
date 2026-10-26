# 能带绘图源数据（2026-10-21）

> **2026-10-26 更正：本目录的 HSE 数据不能作为基准。** 源数据来自 `band_and_elf/huang/<材料>/step4_*`，其 HSE 自洽用的是 `ISYM = 0` 加 12×12×2 的 38 个不可约 k 点列表，电荷密度和交换算符没有按点群对称化，也就是 [`12_s占比与Seebeck有效质量/README.md`](../12_s占比与Seebeck有效质量/README.md) 第 3.2 节和 [`22_能带基准_1-0/README.md`](../22_能带基准_1-0/README.md) 指出的“旧 HSE”错误（各材料 `源数据/<材料>/INCAR_band`、`INCAR_fatband` 中可见 `ISYM = 0`）。正确的 HSE+SOC 结果见 `22_能带基准_1-0/`（`ISYM = 3`，四种材料同一网格）。
> 因此本目录的带隙（0.065 / 0.115 / 0.135 / 约 0.196 eV）、能带形状、轨道 fat-band 图（包括 PbSb 价带顶约 17% 的阳离子 s；新 HSE 为 0%）都只反映旧 HSE，**不要引用**。Bi 系路径上的单点尖峰与此是否有关不再深究。脚本与画法仍可复用于正确的 HSE 能带路径。

本目录是“能带图精细化”这一轮绘图（能带简并/带收敛 + 元素轨道参与度）所用的**源数据、解析结果与脚本索引**。
图本身在 `02_图表/`（或本目录 `图表_预览/`），脚本在 `10_脚本/`。

## 1. 数据来源

| 内容 | 原始位置（计算归档） |
|---|---|
| HSE+SOC 能带路径数据 `band.dat / band_klabels.txt / band_summary.json` | `band_and_elf/huang/<材料>/step4_band_plot/` |
| 路径 lm 投影 `PROCAR_OPT`（Sb 系）/`PROCAR`（Bi 系） | `band_and_elf/huang/<材料>/step4_fatband_plot/` |
| 计算输入 `INCAR_band / INCAR_fatband / KPOINTS / KPOINTS_OPT / POSCAR` | 同上两个目录 |

- 泛函：**HSE06 + SOC**；路径：Γ–M–K–Γ–A–L–H–A|L–M|H–K；每材料 **200–204 个路径点**，96（Sb 系）/120（Bi 系）条能带。
- 能量零点：VBM。
- 原始 `PROCAR`/`PROCAR_OPT` 体积 60–90 MB/材料，未直接入库；其大小与 SHA256 见
  `源数据/<材料>/原始来源.json` 与 `校验/原始大文件_SHA256.json`，可用于核对。
  同一批原始投影在仓库 `07_带边能窗/源数据/<材料>/` 里已有 gzip 版本（`PROCAR_OPT.gz` / `PROCAR.gz`），可直接复算。

## 2. 文件说明

```
源数据/<材料>/
  band.dat             # 第 1 列 k 累积距离 (1/Å)，其后 96/120 列 E-E_VBM (eV)
  band_klabels.txt     # 高对称点标签 + k 距离
  band_summary.json    # 带隙类型/数值、VBM/CBM 的 k 分数坐标
  POSCAR               # 结构（SOC 结构）
  INCAR_band           # 能带计算 INCAR（HSE+SOC）
  INCAR_fatband        # 投影/fatband 计算 INCAR
  KPOINTS / KPOINTS_OPT# 均匀网格 + 路径 k 点（Bi 系无 KPOINTS_OPT，路径在 KPOINTS 里）
  投影_解析.npz         # 绘图直接用的解析数组（见下）
  原始来源.json         # 原始 PROCAR 文件大小/SHA256 + 解析口径
摘要/
  band_summary_refined.{csv,md,json}   # 带隙、VBM/CBM 位置、kT 内谷数、高对称点能级表
校验/
  SHA256.json               # 本目录全部文件校验和
  原始大文件_SHA256.json     # 原始 PROCAR 文件校验和
```

`投影_解析.npz` 字段：

| 字段 | 形状 | 含义 |
|---|---|---|
| `kdist` | (nk,) | 路径累积 k 距离 (1/Å) |
| `E` | (nk, nb) | E − E_VBM (eV)，与 `band.dat` 完全一致 |
| `frac` | (nk, nb, ncat) | 各“元素-轨道”类别的**相对份额**（按该带投影总和归一） |
| `tot` | (nk, nb) | 该带投影总权重（原子球内份额，窗口内约 0.44–0.58） |
| `cats` | (ncat,) | 类别名，如 `Sb-s Sb-p Te-s Te-p Pb-s Pb-p`（= POSCAR 元素顺序 × s,p） |
| `nv` | (1,) | 占据带数（SOC 下 = NELECT） |

## 3. 解析口径（重要）

- SOC 下 PROCAR 每条带有 **4 个块**：**第 1 块是电荷投影**，后 3 块是自旋分量（可正可负，不用于轨道成分）。
- 每块 9 条 ion 行 = 该离子 `s py pz px dxy dyz dz2 dxz x2-y2` 权重 → 按元素归并为 **s / p** 两类（d 留作 other，未画）。
- 因此图里的“参与度”是**原子球内投影的相对份额**（fat-band 常规口径），不是绝对电子数。
- SOC + 反演对称 ⇒ 每条能带 **Kramers 二重简并**（绘图时每对只画一条，图例注明 ×2）。

## 4. 复现

在仓库根目录直接运行即可（脚本自动读 `12_能带绘图_20261021/源数据/`，输出到 `12_能带绘图_20261021/复现输出/`，该目录已在 `.gitignore` 中）：

```bash
pip install -r 10_脚本/requirements.txt              # numpy, matplotlib
python 10_脚本/band_refine.py        # bands_<材料>.png/pdf、4 材料对比、简并/收敛诊断 + band_summary_refined.*
python 10_脚本/orbital_fatband.py    # 每个材料一张 orbital_fatband_<材料>.png/pdf（气泡颜色 = 元素-轨道，面积 = 份额）
```

常用选项：

| 命令 | 作用 |
|---|---|
| `orbital_fatband.py --materials Sn2Bi2Te5` | 只画指定材料 |
| `orbital_fatband.py --emin -2 --emax 2` | 改能窗（默认 VBM −3…+3 eV） |
| `--keep-dup-kpts`（两个脚本都有） | 不剔除 Bi 系落在 SCF 网格上的路径点，复现剔除前的数值（见下文“Bi 系重复 k 点”） |
| `orbital_fatband.py --procar-root 07_带边能窗/源数据` | 不用 npz，从原始 `PROCAR(.gz)` 重新解析投影，并与 `投影_解析.npz` 逐元素比对 |
| `--data-root <dir> --out <dir>`（两个脚本都有） | 指定输入/输出目录；`<dir>/<材料>/step4_band_plot` 结构（原工作区）也能识别 |

复现核对（2026-10-21）：
- `band_refine.py` 生成的 `band_summary_refined.{csv,md,json}` 与 `摘要/` 中的文件**逐字节相同**；
  加 `--keep-dup-kpts` 时 csv/md 与剔除前的旧版逐字节相同（json 仅多出 `dropped_kpts`、`levels[].dk` 两个记录字段）；
- `--procar-root 07_带边能窗/源数据` 重新解析的份额与 `投影_解析.npz` 最大差 3×10⁻⁸（float32 精度）；
  能窗（VBM±3.2 eV）内 PROCAR 能量与 `band.dat` 逐点偏差 ≤ 9 meV（能窗外最高几条空带偏差大，不参与作图）。

`export_plot_source_data.py` 从原工作区 `band_and_elf/huang/<材料>/` 重新生成本目录的打包内容（需要原始 PROCAR，仓库外运行）。

### 轨道 fat-band 的读法

每个材料一张图（文献常用的气泡式 fat-band）：能带为细灰线，每个 (k, 能带) 上按 5 类画气泡，**气泡面积 ∝ 该元素-轨道的份额**：

| 类别 | 颜色 | 标记 |
|---|---|---|
| A-s（Pb-s / Sn-s） | 红 `#e34948` | 圆点 |
| A-p（Pb-p / Sn-p） | 琥珀 `#eda100` | 圆点 |
| B-s（Sb-s / Bi-s） | 绿 `#008300` | 菱形 |
| B-p（Sb-p / Bi-p） | 青 `#1baf7a` | 菱形 |
| Te（s+p） | 蓝 `#2a78d6` | 圆点 |

- Te-s 在 VBM±3 eV 内平均只占 2–9%，并入 Te；份额 < 3% 的气泡不画；气泡每隔 1 个 k 点画一个。
- 五色取自经色觉校验的分类色板（全两两对比正常视觉 ΔE ≥ 15.6）；红–绿、红–青两对对红绿色弱者偏近，故 A/B 再用圆点/菱形区分。
- 全部气泡按面积从大到小排序后绘制，小份额（如 Sn₂Bi₂Te₅ 价带顶约 25% 的 Sn-s）压在大的 Te 气泡上面，不会被遮住。
- Kramers 简并对只画一组气泡（能量、份额取两条平均；能窗内两条份额差 ≤ 0.03），细线仍画全部能带；路径断点 A|L、M|H 处不连线。
- Bi 系默认剔除落在 SCF 网格上的路径点（见下节），图例下方注明剔除点数。

### Bi 系重复 k 点：已剔除（2026-10-21）

Bi 系（Pb₂Bi₂Te₅、Sn₂Bi₂Te₅）用“带权重的 SCF 网格（12×12×2）+ 零权重路径点”的 HSE 方案。路径点里有 14 / 18 个**正好落在 SCF 网格上**，
这些点上 Fock 项 q→0 的处理与其它路径点不同：价带被整体抬高 0.03–0.47 eV（中位值），导带基本不动（非重复点上同一指标只有 ±0.01 eV），
在能带上表现为单点竖直尖峰。这是该方案的已知陷阱，不是物理特征。Γ、M、K、A、L、H 都在网格上，所以高对称点本身也在其中。

两个脚本**默认剔除这些点**（判据：路径点分数坐标 × 网格数 = 整数），画线时跨过被剔除的点、但仍不跨路径断点连线；
能量零点改为剩余点的 VBM，带隙、VBM/CBM 位置按剩余点重算；高对称点能级表取最近的剩余路径点（表中注明 Δk）。
Sb 系用 `KPOINTS_OPT`，不受影响，数值不变。

| 材料 | 剔除点数 | 剔除前（band_summary.json） | 剔除后 |
|---|---:|---|---|
| Pb₂Bi₂Te₅ | 14 | 0.1141 eV 间接；CBM 在 Γ | 0.1148 eV 间接；VBM (0, 0, 0.14)，CBM (0, 0, 0.07)（Γ 本身已剔除） |
| Sn₂Bi₂Te₅ | 18 | 0.1366 eV 间接；VBM 在 (1/3, 0, 1/2)（本点 0.000 eV，两侧邻点 −0.336 / −0.382 eV，是尖峰） | **0.1958 eV** 间接；VBM (0.045, 0, 0.5)（A 附近，比原零点低 0.059 eV），CBM (0, 0, 0.07) |

注意：剔除点只会让 VBM 变低、CBM 变高，且路径上缺了高对称点本身，所以上表带隙是该路径带隙的**上界估计**（带边若恰在高对称点上会被错过）；
`源数据/<材料>/band_summary.json` 保持原样未改。定量结论建议用 `KPOINTS_OPT`（同 Sb 系）重算 Bi 系路径。

## 5. 与本仓库其它目录的关系

- `07_带边能窗/`：带边能窗统计（1/25/50/100 meV 窗口内的轨道占比、入选态清单）与原始 PROCAR gzip；
  本轮绘图复用了同一批投影，但按“元素 × s/p”6 类重新解析并保留全部路径点，用于画 fat-band。
- `02_图表/`：本轮的 `bands_*.png/pdf`、`orbital_fatband_<材料>.png/pdf` 若定稿也放这里（目前在 `图表_预览/`）。
