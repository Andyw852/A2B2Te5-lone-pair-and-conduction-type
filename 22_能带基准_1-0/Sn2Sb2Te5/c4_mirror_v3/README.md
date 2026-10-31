# Sn₂Sb₂Te₅ · c4_mirror_v3：镜面网格找真实谷（HSE06+SOC）

2026-10-30。方法、网格、批次与 [`../../Sn2Bi2Te5/c4_mirror_v3/README.md`](../../Sn2Bi2Te5/c4_mirror_v3/README.md) 完全相同（FOCK_FCC 根因也写在那里），只换了材料。

## 为什么要算

SnSb 的 HSE 价带顶在 Γ。孤对电子谷都在 Γ 下面 130–150 meV：
- Σ（#11，−131.5 meV）
- R₁（#101，−141.0 meV）
- R₂（#120，−144.6 meV）

这三个谷现在只有**路径上的值**。在 SnBi 中，同类谷离开路径进入镜面后升高了 22–38 meV，质量也和路径估计差很多。按 SnBi 的幅度估计，SnSb p 型 S（3e19，300 K）会从现在的 162 µV/K 变到 162–171 µV/K 左右；Fig 8 中 137 µV/K 那组也要相应改动。

导带同样要扫：导带底在 Γ–A 线（kz = 0.143）上，R₁ 位置的导带谷在 +65 meV 处。

## 网格

```
python 10_脚本/c4_mirror_scan.py grid Sn2Sb2Te5 --out 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3 --du 1/36 --umin=-15/36 --split
```

- 16 × 12 = 192 个网格点，加 6 个锚点、1 个 c2 校验点（#101）、4 个对称校验点，共 203 点，分 6 批（36 × 5 + 23）。
- 步长 1/36，12u = −k/3，没有奇数/24 点。
- 文件：
  - `KPOINTS_OPT`（6 批合一）、`KPOINTS_OPT_b01`…`_b06`；
  - `c4_grid.json`：每个点的全局序号，**不要手改 KPOINTS_OPT**。

## 计算：每批一个作业，6 个并行

| 项 | 内容 |
|---|---|
| INCAR | 本目录的：SnSb `c2_hse/INCAR`（`NBANDS = 84`、`TIME = 0.2`、`NELM = 250`、`NKRED 2 2 1`）+ `KPOINTS_OPT_NKBATCH = 36`。**NKRED 不要改** |
| KPOINTS | Gamma 12 12 2（与 c_hse_A 相同） |
| POSCAR、POTCAR | 与 SnSb `c_hse_A`、`c2_hse` 相同 |
| WAVECAR | 每个目录各复制一份已收敛的 SnSb `c_hse_A/WAVECAR`（`ISTART = 1`） |
| KPOINTS_OPT | `bNN` 里放 `KPOINTS_OPT_bNN`，改名为 `KPOINTS_OPT` |
| 资源 | 64 核快节点，`--time=0-12:00:00`，后台 `sleep 41400` 后写 `LABORT = .TRUE.`；SnSb 的 NBANDS 比 SnBi 少，应略快 |

检查与要推送的输出同 SnBi：自洽能量与 c_hse_A 的差 < 1 meV；没有 `internal error`；有 `General timing`；推送 `PROCAR_OPT.gz`、`OSZICAR`、`out.log` 尾部。

## 分析

```bash
python 10_脚本/c4_mirror_scan.py analyze 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3 \
  --parts 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3/b0{1,2,3,4,5,6} --split --win-v 0.20 \
  | tee 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3/analyze.log
```

- **`--win-v 0.20` 必须加**。默认 0.060 eV 会把离 Γ 130–150 meV 的孤对电子谷全部排除，c5 就没有点可算。
- **验收**：校验点与 c2 路径的差 < 0.5 meV，对称校验 < 0.1 meV。价带第一行注明“在 umin 边界上”时停下汇报。
- **通过后**：脚本写出 `c4_extrema.json`，并在 `../c5_star/` 写出 c5 的 KPOINTS_OPT（含分批文件）。c5 的设置同上，可以直接提交。
- 若有极大只停在网格点上、没能二次拟合（`网格点（未细化）`），用 `local` 子命令出 c6，做法同 SnBi 的 `c6_local`。
- **c5 拟合**：

  ```bash
  python 10_脚本/c4_mirror_scan.py fit 22_能带基准_1-0/Sn2Sb2Te5/c5_star --parts <c5 各批目录>
  ```

**注意**：`multivalley_seebeck.py` 的 `true_load` 目前只为 SnBi 读真实谷。SnSb 的 c4/c5 推上来后，由 Claude 接入并检查 Γ 谷（在 Γ–A 三重轴上，不进 c5 星形）的处理。

## 结果（2026-10-31，`analyze.log`、`../c5_star/fit.log`）

计算由 DeepSeek 完成（`9e55d1a`）。
- c4 校验 0.047 meV、对称校验 0.0016 meV；c5 校验 0.042 meV。
- 六批与 c5 的自洽能量都是 −41.37891726 eV。

**带边**：价带顶仍在 Γ，镜面内没有更高的谷；带隙 160.3 meV，与路径相同。导带底仍在 Γ–A 线上。c4/c5 的“导带谷” (±0.006, 0, ∓0.163) 高 0.7–0.9 meV，星形拟合的驻点退回到 u ≈ 0，就是路径上的 Γ–A 驼峰，不是新谷。

**孤对电子价带谷**（名称按最近的路径谷，与 SnBi 相同的规则）：

| 谷 | 中心 (u, 0, w) | 相对价带顶 | 阳离子 s | 主质量（mₑ） | 来源 |
|---|---|---:|---:|---|---|
| R₂′ | (−0.3459, 0, 0.4083) | −84.0 meV | 33.4% | 0.05、0.069、0.318 | c5，真极值 |
| R₁′ | (−0.1111, 0, −1/4) | ≥ −123.5 meV | 26.7% | 网格粗估 | c4 网格点，未细化 |
| Σ′ | (−0.1774, 0, 0.0625) | −127.1 meV | 26.1% | 0.039、0.08、2.23 | c5，真极值 |

- R₂′ 比路径上的 R₂（−146 meV）高 62 meV；Σ′ 只比路径上的 Σ（−131.5）高 4 meV。
- `analyze.log` 把 R₂′ 写成“Σ 谷”，是因为命名规则不同。按位置，它对应 SnBi 的 R₂′ (−0.335, 0, 0.396)。
- 另有两个网格点：(−0.1389, −0.0833) 在 −129.8 meV，(−0.0556, 5/12) 在 −130.3 meV。它们在粗网格上不是二维极大（曲率不全为负），`multivalley_seebeck.py` 不计入。

**对 S 的影响**：SnSb p 型 S（3e19，300 K）由 137（路径）改为 164 µV/K。不计 R₁′ 时为 134，所以要做 c6，见 [`../c6_local/README.md`](../c6_local/README.md)。
