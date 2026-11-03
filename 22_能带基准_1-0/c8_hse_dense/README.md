# Pb2Sb2Te5 · c8_hse_dense：HSE06+SOC 加密网格（28x28x7，515 个不可约点）

2026-11-03。用途：直接用 HSE 能带算输运（BoltzTraP2），替换报告里的比值法。网格为什么取 28x28x7，见
[`21_阶段报告与交接/交接_HSE直接BoltzTraP2_20261102.md`](../../21_阶段报告与交接/交接_HSE直接BoltzTraP2_20261102.md)。
本目录文件由 `python 10_脚本/hse_dense_kpoints.py` 生成，**不要手改**。

做法与 `Sn2Bi2Te5/c4_mirror_v3` 相同：在 12×12×2 的 HSE 自洽之上，用 `KPOINTS_OPT` 算更多 k 点（非自洽）。

| 文件 | 内容 |
|---|---|
| `INCAR` | 与 c2_hse 相同（`NKRED 2 2 1`、`PRECFOCK = Fast`、`ISTART = 1`），只改了 SYSTEM、`NBANDS = 96`，加 `KPOINTS_OPT_NKBATCH = 36` |
| `KPOINTS` | 12×12×2，与自洽相同 |
| `KPOINTS_OPT_b01` … `_b15` | 15 批，每批第一个点是 Γ（锚点，用来核对各批一致），其余为新点 |
| `KPOINTS_OPT` | 全部 515 个点合在一起（单个作业顺序跑时用） |
| `c8_grid.json` | 每个点的坐标与批号 |

所有点都不在 FOCK_FCC 判据边界上（12·k 不为半整数）。

## 0 只对 Pb₂Sb₂Te₅：先跑 b00（HSE 自洽，产生 WAVECAR）

- 目录 `b00/`：本目录的 INCAR，但改两处：`LWAVE = .TRUE.`；删掉 `KPOINTS_OPT_NKBATCH` 一行。**不放 KPOINTS_OPT。**
- KPOINTS、POSCAR、POTCAR 与 c2_hse 相同；WAVECAR 用 c2_hse 当时起步的那个 PBE WAVECAR（`ISTART = 1`）。
- c2_hse 的自洽从 PBE 起步，连同 204 个路径点共 17 小时 43 分；b00 只做自洽，预计 10 小时左右，`--time=1-00:00:00`。
- 检查：`aborting loop because EDIFF is reached`；带隙应与 c2_hse 一致（Γ 处约 51 meV）。
- 之后各批从 `b00/WAVECAR` 起步。

## 1 各批（建议每批一个作业，并行）

- 目录 `b01` … `b15`：本目录的 INCAR、KPOINTS；`KPOINTS_OPT_bNN` 改名为 `KPOINTS_OPT`。
- POSCAR、POTCAR 与 c2_hse 相同。
- WAVECAR：`b00/WAVECAR`，给每个目录各复制一份。**不要用 PBE 的。**
- 资源：64 核快节点，`--time=0-12:00:00`（c4 的经验：36 个点约 6–7 小时），后台 `sleep 41400` 后写 `LABORT = .TRUE.`。

**每个作业的检查**
1. 自洽：从 WAVECAR 起步，几步内收敛。
2. 出现 `k-point batch` 后有 DAV 迭代。出现 `internal error` 就贴报错原文和批号，这个作业停下，其余继续。
3. 结尾有 `General timing` 和 `aborting loop because EDIFF is reached`。

## 2 合并成 BoltzTraP2 输入，再跑 BoltzTraP2

```
python 10_脚本/vasprun_to_gene.py --opt --mesh 28 28 7 --material Pb2Sb2Te5 --functional hse \
    --vasprun b01/vasprun.xml … b15/vasprun.xml --out 22_能带基准_1-0/c8_hse_dense/gene
python 10_脚本/run_bt2_hse_grid.py --autozt <AutoZT step8.1 目录> --material Pb2Sb2Te5 --functional hse \
    --run-dir 22_能带基准_1-0/c8_hse_dense/gene --mult 60 --out 22_能带基准_1-0/c8_hse_dense/gene/bt2_mult60 --kcheck 22_能带基准_1-0/bt2_check_kpoints.txt
```

- 第二条再用 `--mult 120` 跑一次（输出 `bt2_mult120`），用来核对插值倍数已收敛。
- 合并时会打印锚点一致性（各批在 Γ 上的本征值最多差，应 < 1 meV）。还会核对恰好是一套 515 个不可约点。不符就报错，停下汇报。

**要推送的**：
- `gene/`（GENE 文件、`source.json`、两个 `bt2_mult*/`）；
- 每批的 `OSZICAR` 与 `out.log` 中从 `k-point batch` 开始的部分。

**vasprun.xml 不推送**，留在集群上。
