# Sn₂Bi₂Te₅ · c4_mirror_w2：按 w 分组、逐批作业的镜面网格（HSE06+SOC）

2026-10-29。目的与 [`../c4_mirror/README.md`](../c4_mirror/README.md) 相同：在镜面内找价带真实极大、导带真实极小。网格不变（u 17 个 × w 12 个），**只改了点在批里的排法**，并且每批单独一个作业。INCAR、KPOINTS 与 c4_mirror 逐字相同（`NKRED 2 2 1` 不变，`KPOINTS_OPT_NKBATCH = 36`）。

## 为什么改

c4_mirror 的两种排布都在第 2 批崩溃（`internal error in FOCK_FCC: number of k-points incorrect 68 288 2 2 1`，作业 313107、315591）。

集群侧统计了已知的 12 个批：
- 崩溃的批都是 kz 取 12 个值、Γ–A 线上点很少的批；
- kz 取 12 个值而通过的两批，都含整条 Γ–A 线；
- kz 不超过 8 个值的批（c2 六批、c3）全部通过。

这只是经验规律，机制不明。同一组数据也能被“批内是否有 k 与 −k 成对”完全分开，因为含整条 Γ–A 线的批正好有这样的点对。数据分不清是哪一个特征起作用，所以这里不依赖任何一种解释，只做两件事：

1. **每批只含 2 个 w**（相差 1/2），加上锚点，批内 kz 只有 2–3 个值。这与从未出错的 c2、c3 同属“kz 少”的一类。
2. **每批单独一个作业**。崩溃发生在一批开始时（DAV 0 步），所以每个作业开跑约 1 小时就知道这一批能不能算。某批崩溃只损失这一批，其余照常；6 批并行，墙钟约 6–7 小时，原方案约 1.5 天。

## 点的排布

由 `python 10_脚本/c4_mirror_scan.py grid Sn2Bi2Te5 --out 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_w2 --wgroup 2 --split` 生成。

| 批 | 文件 | 点数 | 网格的 w | 批内 kz 值 | 锚点 |
|---|---|---|---|---|---|
| 1 | `KPOINTS_OPT_b01` | 36 | −1/2、0 | 2 个 | #97 A；另含校验点 #101 |
| 2 | `KPOINTS_OPT_b02` | 36 | −5/12、1/12 | 3 个 | #11 |
| 3 | `KPOINTS_OPT_b03` | 36 | −1/3、1/6 | 3 个 | #22 |
| 4 | `KPOINTS_OPT_b04` | 36 | −1/4、1/4 | 3 个 | #33 M |
| 5 | `KPOINTS_OPT_b05` | 36 | −1/6、1/3 | 3 个 | #97 A |
| 6 | `KPOINTS_OPT_b06` | 35 | −1/12、5/12 | 3 个 | #108 |

- 每批 = 1 个锚点 + 2 个 w × 17 个 u + 补位。
- 批 2–5 各有 1 个补位点：批内某个网格点的另一个旋转表示（kz 不变），兼作对称校验。
- 网格点仍轮流写成三个等价表示 (u,0,w)、(−u,u,w)、(0,−u,w)，与 c4_mirror 第二版相同。
- `KPOINTS_OPT` 是 6 批按顺序合在一起的文件：只在单个作业顺序跑全部 6 批时用，此时 `NKBATCH = 36` 正好切出这 6 批。**逐批作业用 `_b01`…`_b06`。**
- `c4_grid.json` 记录每个点的全局序号，分析脚本靠它取点。**不要手改任何 KPOINTS_OPT**，要改就改脚本重新生成。

## 计算（每批一个作业，b01–b06）

| 项 | 内容 |
|---|---|
| 目录 | 集群上建 `c4_mirror_w2/b01` … `b06` |
| INCAR、KPOINTS | 本目录的（与 c4_mirror、c3_valley 相同，只多 `KPOINTS_OPT_NKBATCH = 36`）。**NKRED 不要改** |
| POSCAR、POTCAR | 与 c_hse_A、c3_valley 相同 |
| WAVECAR | 从已收敛的 `c_hse_A/WAVECAR` 复制一份到每个目录（不要 6 个作业共用一个文件），不要用 PBE 的 |
| KPOINTS_OPT | `bNN` 目录里放 `KPOINTS_OPT_bNN`，改名为 `KPOINTS_OPT` |
| 资源 | 与 c4_mirror 相同（64 核快节点）。墙钟 `--time=0-12:00:00`，后台 `sleep 41400` 后写 `LABORT = .TRUE.` |

**每个作业的检查**
1. **自洽**：从 c_hse_A 的波函数起步，几步内收敛，energy(sigma->0) 与 c_hse_A 相差 < 1 meV。
2. **批开始后**：`out.log` 里 `k-point batch [1-36\36]`（b06 为 `[1-35\35]`）之后出现 DAV 迭代，说明过了 FOCK_FCC。出现 `internal error` 就记下报错原文，这个作业停下，其余继续。
3. **结尾**：`General timing` 与 `aborting loop because EDIFF is reached` 都要有。

**要推送的输出**（每个 `bNN/`）：`PROCAR_OPT.gz`（必须）、`OSZICAR`、`out.log` 中从 `k-point batch` 开始的部分。

## 分析（6 批都完成后，在集群上或本地）

```bash
python 10_脚本/c4_mirror_scan.py analyze 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_w2 \
  --parts 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_w2/b0{1,2,3,4,5,6} --c5-group --split \
  | tee 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_w2/analyze.log
```

- `--parts` 按批号顺序读 6 个 PROCAR_OPT，并逐批核对点数。
- 开头两行是验收数据：
  - 7 个校验点与 c2 路径的最大差 < 0.5 meV（6 个作业各自自洽，锚点同时检验作业之间是否一致）；
  - 对称校验 < 0.1 meV。

  任一不满足就停下汇报。
- 价带列表第一行若注明“在 umin 边界上”，停下汇报。

满足的话，脚本写出 `c4_extrema.json`，并在 `../c5_star/` 写出 c5 的 `KPOINTS_OPT_b01` …：每个谷一批（锚点 + 11 点星形 + 补位），批内 kz 只有 3–5 个值。

- **c5 的跑法**：每批同样单独一个作业，设置同上，只换 KPOINTS_OPT。
- **c5 的分析**：`python 10_脚本/c4_mirror_scan.py fit 22_能带基准_1-0/Sn2Bi2Te5/c5_star --parts <c5 各批目录，按批号>`。
- **c5 不用等 Claude**：第 1 步验收通过后直接提交。

脚本的 `--parts`、`--c5-group`、`--split` 已用合成数据测过：
- 构造满足 D3d 对称的合成能带；
- 逐批读入与单个 PROCAR_OPT 的结果逐字节相同；
- 少一批时报错；
- c5 星形能还原合成能带的极值与质量。
