# 能带绘图源数据（2026-10-21）

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
python 10_脚本/orbital_fatband.py    # orbital_fatband_4mats.png/pdf（颜色 = 元素-轨道，层厚 = 参与份额）
```

常用选项：

| 命令 | 作用 |
|---|---|
| `orbital_fatband.py --materials Sn2Bi2Te5` | 只画单材料 → `orbital_fatband_Sn2Bi2Te5.png/pdf` |
| `orbital_fatband.py --procar-root 07_带边能窗/源数据` | 不用 npz，从原始 `PROCAR(.gz)` 重新解析投影，并与 `投影_解析.npz` 逐元素比对 |
| `--data-root <dir> --out <dir>`（两个脚本都有） | 指定输入/输出目录；`<dir>/<材料>/step4_band_plot` 结构（原工作区）也能识别 |

复现核对（2026-10-21）：
- `band_refine.py` 生成的 `band_summary_refined.{csv,md,json}` 与 `摘要/` 中的文件**逐字节相同**；
- `--procar-root 07_带边能窗/源数据` 重新解析的份额与 `投影_解析.npz` 最大差 3×10⁻⁸（float32 精度）；
  能窗（VBM±3.2 eV）内 PROCAR 能量与 `band.dat` 逐点偏差 ≤ 9 meV（能窗外最高几条空带偏差大，不参与作图）。

`export_plot_source_data.py` 从原工作区 `band_and_elf/huang/<材料>/` 重新生成本目录的打包内容（需要原始 PROCAR，仓库外运行）。

### 轨道 fat-band 的读法

每条能带画成总线宽固定的一根“线缆”，由内到外依次是 **A-s、B-s、Te-s、A-p、B-p、Te-p** 六层，
每层厚度 ∝ 该元素-轨道在该 (k, 带) 的份额（份额 < 1% 的层不画）。颜色：Te 蓝、A(Pb/Sn) 红、B(Sb/Bi) 绿；s 浅色、p 深色。
因为逐层累积线宽、内层后画，小份额层（如 Sn₂Bi₂Te₅ 价带顶约 25% 的 Sn-s）不会被大份额的 Te-p 层盖住。
路径断点 A|L、M|H 处不连线。

## 5. 与本仓库其它目录的关系

- `07_带边能窗/`：带边能窗统计（1/25/50/100 meV 窗口内的轨道占比、入选态清单）与原始 PROCAR gzip；
  本轮绘图复用了同一批投影，但按“元素 × s/p”6 类重新解析并保留全部路径点，用于画 fat-band。
- `02_图表/`：本轮的 `bands_*.png/pdf`、`orbital_fatband_4mats.png/pdf` 若定稿也放这里（目前在 `图表_预览/`）。
