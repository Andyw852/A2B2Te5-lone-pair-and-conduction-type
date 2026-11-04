# SOC 机制检验的 HSE 版：HSE06 不加 SOC（2026-11-04）

**目的**：HSE 能带专题（`21_阶段报告与交接/HSE能带专题_20261104.html`）要求全部为 HSE 级别。图 H8（Γ 点能级与 SOC）要对比加、不加 SOC，加 SOC 的 HSE 已有（`<材料>/c2_hse`），**缺的是 HSE06 不加 SOC**。本目录就是这四个计算。PBE 版见 `../SOC机制检验/`。

## 四个作业

`e1_hse_nosoc_Pb2Sb2Te5`、`e1_hse_nosoc_Pb2Bi2Te5`、`e1_hse_nosoc_Sn2Sb2Te5`、`e1_hse_nosoc_Sn2Bi2Te5`。

- **优先** Sn₂Sb₂Te₅、Sn₂Bi₂Te₅：图 H8 只画这两种。
- Pb 系用于表里 Pb 系的拆分（标量 / SOC）。

每个作业目录内由 `submit.sh` 依次跑三步，中断后重交会跳过已完成的步：

| 步 | INCAR | 内容 | 产物 |
|---|---|---|---|
| 1 | `INCAR.1_pbe` | PBE 不加 SOC 自洽（与 `../SOC机制检验/d1_nosoc_<材料>` 相同，只改 `LWAVE = .TRUE.`、ENCUT 与 c2_hse 一致） | WAVECAR |
| 2 | `INCAR.2_hse_scf` | HSE06 不加 SOC 自洽，从第 1 步的 WAVECAR 续算（`ISTART = 1`） | WAVECAR |
| 3 | `INCAR.3_hse_path` | 从第 2 步续算，加 `KPOINTS_OPT`（c2 路径 204 点，`KPOINTS_OPT_NKBATCH = 36`） | **PROCAR_OPT**、vasprun.xml |

**与 `<材料>/c2_hse` 的差别**：INCAR 逐项相同，只做两处改动。
- 去掉 `LSORBIT`、`GGA_COMPAT`、`MAGMOM`、`AMIX_MAG`、`BMIX_MAG`，用 **vasp_std**。
- `NBANDS` 取 PBE 无 SOC 时用过的值：Pb₂Sb₂Te₅ 48、Pb₂Bi₂Te₅ 54、Sn₂Sb₂Te₅ 42、Sn₂Bi₂Te₅ 54。

其余照旧，包括 HSE06（`AEXX = 0.25`、`HFSCREEN = 0.2`）、`NKRED 2 2 1`、`PRECFOCK = Fast`、`ISYM = 3`、`ALGO = Damped`、`TIME = 0.2`、ENCUT、12×12×2。

**文件**
- `KPOINTS`：c2_hse 的 12×12×2。
- `KPOINTS_OPT.path`：`../c2_path/KPOINTS_OPT`，第 3 步前由脚本复制为 `KPOINTS_OPT`。这条路径以前在 c2_hse 上跑通过，FOCK_FCC 的“奇数/24”检查没有问题。**不要手改。**
- `POSCAR`：`A2B2Te5_SOC_dense/materials/<材料>/POSCAR`，与 PBE 的 d1 相同。
- `POTCAR`：从 c2_hse 的作业目录复制，元素顺序同 POSCAR 第 6 行（B 位、Te、A 位）。

`submit.sh` 里的分区、模块、vasp_std 路径照抄了 d1（jzzn）。若 c2_hse 当时在别的集群上跑，按那边的环境改。

## 检查

- 三步都有 `General timing` 与 `aborting loop because EDIFF is reached`（脚本已检查）。
- PROCAR_OPT 含 204 个 k 点。
- 第 2 步自洽与第 3 步路径在 Γ 点的最高价带、最低导带相差 < 1 meV。
- **能带反转应当消失**：不加 SOC 时 Γ 处的价带顶应是孤对电子 s 态，B 位 p 态在导带。这一点与 PBE 的 d1 相同。若不是，先报告，不要继续。

## 交回

每个目录交回：`PROCAR_OPT.gz`、`vasprun.xml.gz`（第 3 步）、`OUTCAR.1/2/3`（或其末尾 200 行）、`OSZICAR.1/2/3`、`out.log.3`、`能量.txt`（Γ 点、Σ 谷附近几条带的能量）。

## 收到后（Claude 做）

```
python 10_脚本/soc_mechanism.py --functional hse --fig <输出目录>     # 打印 HSE 的能级与拆分
python 10_脚本/hse_band_figures.py --only H8 T                        # 有了 e1 结果，图 H8 自动改用 HSE
python 10_脚本/build_hse_bands_report.py
```

Σ 谷取 Γ–M 段内最高价带的局部极大。2026-11-04 起这样取；原来取的是段上最大值，不加 SOC 时会落在紧挨 Γ 的点，见 `../SOC机制检验/README.md` 开头的改正。
