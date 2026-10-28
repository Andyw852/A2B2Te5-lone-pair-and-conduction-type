# Sn₂Bi₂Te₅ · c4_mirror + c5_star：镜面内离轴能谷的真实极值与质量张量（HSE06+SOC）

2026-10-28，用户同意补算。**请在 `claude/keen-sagan-15fv7v` 分支上工作并推送**（报告的最新版本在这个分支上，脚本也只在这里）：

```bash
git fetch origin claude/keen-sagan-15fv7v
git checkout -B claude/keen-sagan-15fv7v origin/claude/keen-sagan-15fv7v
```

## 为什么要算

c3_valley 表明，SnBi 价带在 Γ–M、A–L 线上的三组谷，在路径上的点都不是三维极值：#11 是鞍点，#101、#119 的极值偏离路径约 0.05 Å⁻¹（见 `../c3_valley/README.md`）。这两条线在镜面内、不在二重轴上，对称性不要求沿 c 的梯度为零，所以真实极值在镜面内的一般 kz 处。

后果是：SnBi 的真实价带顶比路径上的高，路径带隙 0.182 eV 只是上限；p 型 S 是暂定值。n 型中承担约一半电导的 A–L 导带口袋同样在这个镜面里。

本计算分两步：
1. **c4_mirror**：在镜面内均匀取点，找出真实的价带极大和导带极小。
2. **c5_star**：在找到的极值处做与 c3 相同的 11 点星形，拟合完整的质量张量。

## 设置（两步相同，都和 c3_valley 一样，只换 KPOINTS_OPT）

| 项 | 内容 |
|---|---|
| INCAR | 用本目录的 `INCAR`（与 `../c3_valley/INCAR` 逐字相同：HSE06+SOC，`ISTART = 1`，`ISYM = 3`，`LORBIT = 11`，`PRECFOCK = Fast`，`NKRED 2 2 1`，`TIME = 0.2`） |
| POSCAR、POTCAR、KPOINTS | 与 c_hse_A、c3_valley 相同（KPOINTS 为 Γ 12×12×2，本目录已附） |
| WAVECAR | 从已收敛的 `c_hse_A/WAVECAR` 起步，与 c3 相同。**不要用 PBE 的** |
| KPOINTS_OPT | 第 1 步用本目录的（208 点）；第 2 步用脚本生成的 `../c5_star/KPOINTS_OPT`（约 40 点） |
| 墙钟 | 第 1 步与 c2_hse（204 点）相当：`--time=3-00:00:00`，后台 `sleep 253800` 后写 `LABORT = .TRUE.`；第 2 步 1 天足够 |
| 内存 | 不加 `KPOINTS_OPT_NKBATCH`；内存不够再加 |
| 结尾检查 | `General timing` 与 `aborting loop because EDIFF is reached` 都要有，否则 `exit 1` |

## 第 1 步：c4_mirror（本目录）

`KPOINTS_OPT`：镜面上的点 (u, 0, w)，倒格矢分数坐标。u 从 0 到 −0.40，步长 0.025（17 个，u = 0 是 Γ–A 线）；w 从 −1/2 到 5/12，步长 1/12（12 个）。再加 4 个 c2 路径校验点（#0 Γ、#11、#101、#119）。共 208 点。

由于 (u, w) 与 (−u, −w) 等价，这已覆盖整个镜面中 |u| ≤ 0.40 的部分。点的几何见 `c4_grid.json`，由 `python 10_脚本/c4_mirror_scan.py grid Sn2Bi2Te5 --out 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror` 生成。

**验收**：
1. HSE 自洽从 c_hse_A 的波函数起步，应在几步内收敛，energy(sigma->0) 与 c_hse_A 相差 < 1 meV。超过 15 个电子步或能量对不上，说明 WAVECAR 没读进，停下汇报。
2. 算完后在集群上直接运行（只需 numpy）：
   ```bash
   python 10_脚本/c4_mirror_scan.py analyze 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror | tee 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror/analyze.log
   ```
   第一行是校验点与 c2 路径的最大差，应 < 0.5 meV；否则停下汇报。脚本会写出 `c4_extrema.json`，并生成第 2 步的 `../c5_star/KPOINTS_OPT` 与 `c5_points.json`。
3. 如果价带列表第一行注明“在 umin 边界上”，说明最高的极大可能在 u < −0.40 处。这时停下汇报，需要把网格扩到 −0.5 再算。

## 第 2 步：c5_star（`../c5_star/`）

第 1 步的 analyze 会为以下极值各写一组 11 点星形：离价带顶 60 meV 内的价带极大，以及离导带底 100 meV 内的导带极小。星形与 c3 相同：沿线、面内垂直、沿 c 各 ±δ，加 (线, c) 平面对角点，δ = 0.02、δz = 0.03 Å⁻¹。另外再加 4 个校验点。

**第 1 步验收通过后直接提交第 2 步，不用等 Claude。** 设置同上，只把 KPOINTS_OPT 换成 `../c5_star/KPOINTS_OPT`。算完后运行：

```bash
python 10_脚本/c4_mirror_scan.py fit 22_能带基准_1-0/Sn2Bi2Te5/c5_star | tee 22_能带基准_1-0/Sn2Bi2Te5/c5_star/fit.log
```

输出包括：
- 每个谷的梯度、主质量，以及驻点是否落在采样半径内。
- 真实价带顶、导带底和带隙。
- 结果写入 `c5_masses.json`。

如果离价带顶 30 meV 内有谷注明“未确认为真极值”，照常推送，并在汇报里写明。

## 要推送的文件（推到 `claude/keen-sagan-15fv7v`）

| 目录 | 文件 |
|---|---|
| `c4_mirror/` | `PROCAR_OPT.gz`（**必须**：分析用它的 8 位小数能量）、`INCAR`、`KPOINTS`、`KPOINTS_OPT`、`能量.txt`（energy(sigma->0)、E-fermi、EDIFF reached 行）、`analyze.log`、`c4_extrema.json` |
| `c5_star/` | `PROCAR_OPT.gz`（必须）、`INCAR`、`KPOINTS`、`KPOINTS_OPT`、`c5_points.json`、`能量.txt`、`fit.log`、`c5_masses.json` |

提交信息用中文，不写模型名。汇报“已入库”前用 `git ls-files` 确认。

## 之后（Claude）

- 把 c5 的谷（真实位置、能量与质量张量）接入 `10_脚本/multivalley_seebeck.py`，重算 SnBi 的 p 型和 n 型 S，去掉报告里的“暂定”。
- 用真实价带顶更新 SnBi 的带隙（表 2、图 2）和图 4 中孤对电子谷的深度。
