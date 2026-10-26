# c2：四种材料的 HSE+SOC 与 PBE+SOC 能带路径

用户 2026-10-26 选定（不做 Wannier + BoltzTraP2）。目的：

- 得到正确的 HSE 能带图和轨道 fat-band 图，替换 `12_能带绘图_20261021/` 中作废的旧 HSE 图；
- 补上 12×12×2 网格看不到的地方：网格点之间的极值位置（如 Γ–A 线中间、SnBi 两组相差 17 meV 的价带谷），以及带边有效质量；
- 用多谷抛物带模型估同浓度 S（p、n 分开），并与同一路径上的 PBE 对比。

## 路径

`KPOINTS_OPT`：Γ–M–K–Γ–A–L–H–A|L–M|H–K，204 点，倒格子分数坐标，与 12_ 中 Sb 系用过的 seekpath 路径相同。四种材料的 POSCAR 都是同一 120° 六方设置，所以这一份文件四种材料通用。

| 高对称点 | Γ | M | K | Γ | A | L | H | A\|L | M\|H | K |
|---|---|---|---|---|---|---|---|---|---|---|
| 路径点序号（从 0 起） | 0 | 33 | 52 | 90 | 97 | 130 | 149 | 187/188 | 195/196 | 203 |

其中 18 个点落在 12×12×2 自洽网格上（序号 0、11、22、33、52、71、90、97、108、119、130、149、168、187、188、195、196、203），用来与自洽部分的本征值核对。

## 作业设置（每种材料两个）

| 项 | HSE 路径（`c2_hse`） | PBE 路径（`c2_pbe`） |
|---|---|---|
| 起点 | 复制 c_hse_A 的 INCAR、POSCAR、POTCAR、KPOINTS（PbSb 用 c1 的，从 clus_xa 取） | 复制 a_pbe 的（PbSb 用 jzzn 上 b 的） |
| 新增 | 本目录的 `KPOINTS_OPT` | 同左 |
| WAVECAR | 已收敛的 `c_hse_A/WAVECAR`（PbSb 用 c1 的）；脚本只在目录里没有 WAVECAR 时复制。**不要用 a_pbe 的** | 不需要：a_pbe 的 INCAR 本来就是 `ISTART = 0`，PBE 自洽约 1 分钟 |
| INCAR 改动 | `TIME = 0.2`、`NELM = 250`（PbSb 的 c1 原为 0.4、100，也照改）；`LWAVE = .FALSE.`、`LCHARG = .FALSE.`；保留 `ISYM = 3`、`LORBIT = 11`；不加 `KPOINTS_OPT_NKBATCH`（内存不够再加） | `LWAVE = .FALSE.`、`LCHARG = .FALSE.`；保留 `LORBIT = 11` |
| 墙钟 | `--time=3-00:00:00`，后台 `sleep 253800` 后写 `LABORT = .TRUE.` | 默认 |
| 结尾检查 | General timing 和 `aborting loop because EDIFF is reached`，任一缺失就 `exit 1` | 同左 |

## 验收

1. HSE 自洽部分从已收敛的波函数起步，应在几步内收敛，energy(sigma->0) 与 c_hse_A 相差 < 1 meV（PbSb 与 c1 相差 < 2 meV）。超过 15 个电子步或能量对不上，说明 WAVECAR 没读进，停下汇报。
2. 18 个落在网格上的路径点：第 NELECT、NELECT+1 条带的能量与 c_hse_A（PbSb：c1）的 EIGENVAL 相差 ≤ 2 meV（Claude 复核）。
3. 路径上没有单点尖峰（Claude 复核）。

## 推送位置

`22_能带基准_1-0/<材料>/c2_hse/`、`c2_pbe/`（PbSb 在 `22_能带基准_1-0/c2_hse/`、`c2_pbe/`）：INCAR、KPOINTS、KPOINTS_OPT、`PROCAR_OPT.gz`、`能量.txt`（energy(sigma->0)、E-fermi、EDIFF reached 行）、EIGENVAL（若另有 EIGENVAL_OPT 也推）。`vasprun.xml` 压缩后不超过 50 MB 就一并推 `vasprun.xml.gz`。

## 之后（Claude）

核对上面第 2、3 条；画四种材料的 HSE 能带图和轨道 fat-band 图（脚本 `10_脚本/band_refine.py`、`orbital_fatband.py` 改读 PROCAR_OPT）；读出价带、导带各谷的位置、能量和有效质量；用多谷抛物带模型估同浓度 S，p、n 分开，HSE 与 PBE 对比。
