# Sn₂Sb₂Te₅ · c6_local：三个网格极大的局部细网格（HSE06+SOC）

2026-10-31。做法、设置与 [`../../Sn2Bi2Te5/c6_local/README.md`](../../Sn2Bi2Te5/c6_local/README.md) 完全相同（那边已跑通），只换了材料。

## 为什么要算

SnSb 的 c4 中有三个价带极大只停在网格点上：

| 标签 | c4 中心 (u, w) | c4 网格能量（相对价带顶） | 说明 |
|---|---|---:|---|
| g0 | (−0.1111, −1/4) | −123.5 meV | R₁′；已按网格粗估计入 S |
| g1 | (−0.1389, −1/12) | −129.8 meV | 粗网格上曲率不全为负，未计入 S |
| g2 | (−0.0556, 5/12) | −130.3 meV | 位置与 SnBi 的 R₁″ 相同；未计入 S |

- **为什么要细化**：3×10¹⁹ cm⁻³ 时 SnSb 的费米能级在价带顶下 113 meV，这三个谷正好在费米能级附近。
- **对 S 的影响**：p 型 S 现为 164 µV/K。完全不计 g0 时为 134 µV/K。细化后数值才能定下来。

## 网格

由下面的命令生成：

```
python 10_脚本/c4_mirror_scan.py local 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3 --out 22_能带基准_1-0/Sn2Sb2Te5/c6_local --win-v 0.20 --split
```

- **每个谷**：u 方向 7 个点、w 方向 9 个点，步长 1/60 × 1/36；另加 ±δ 两个镜面外的点。
- **覆盖**：沿线 ±0.083 Å⁻¹、沿 c ±0.041 Å⁻¹。
- **点数与分批**：加锚点与校验点共 202 点，`KPOINTS_OPT_NKBATCH = 36`，分 6 批（36 × 5 + 22）。
- **FOCK_FCC 检查**：所有点都通过奇数/24 检查。
- **不要手改 `KPOINTS_OPT`。**

## 计算

| 项 | 内容 |
|---|---|
| INCAR、KPOINTS | 本目录的，从 SnSb `c4_mirror_v3` 原样复制（`NBANDS = 84`，`NKRED 2 2 1`，`KPOINTS_OPT_NKBATCH = 36`） |
| POSCAR、POTCAR、WAVECAR | 与 SnSb `c_hse_A` 相同，`ISTART = 1` |
| KPOINTS_OPT | `b01`…`b06` 里各放对应的 `KPOINTS_OPT_bNN`，改名为 `KPOINTS_OPT` |
| 资源 | 每批 64 核，`--time=0-12:00:00` |

**每个作业的检查**：
- 自洽能量与 c4_mirror_v3 一致（−41.37891726 eV）；
- 没有 `internal error`；
- 有 `General timing`。

**要推送的输出**：每批的 `PROCAR_OPT.gz`、`OSZICAR`、`out.log` 尾部。

## 分析

```bash
python 10_脚本/c4_mirror_scan.py local-fit 22_能带基准_1-0/Sn2Sb2Te5/c6_local \
  --parts 22_能带基准_1-0/Sn2Sb2Te5/c6_local/b0{1,2,3,4,5,6} \
  | tee 22_能带基准_1-0/Sn2Sb2Te5/c6_local/local_fit.log
```

- `local-fit` 已改用整块补丁的双三次样条，SnBi 的 R₁″ 就是这样确认的。
- **验收**：校验点与路径的差 < 0.5 meV。
- **以下情况停下汇报**：
  - 最高点落在补丁的最外一圈（脚本会注明，并且不写入这个谷）；
  - 某个谷报“未确认”。
- **推上来以后**：`multivalley_seebeck.py` 会自动读入 `c6_masses.json`。g1、g2 若确认是极大，也会替换掉原来不计入的网格点。之后由 Claude 重算 S、重出图表。
