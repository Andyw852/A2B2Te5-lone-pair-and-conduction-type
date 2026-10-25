# 退火前后主导缺陷：D/V/O 缺陷计算的 VASP 流程与分析

2026-10-18。本目录把 [`补算实施方案_20261007`](../补算实施方案_20261007/开始阅读.txt) 中已有的 Pb₂Sb₂Te₅ D/V/O 几何（9 个母胞、234 个缺陷模板，全部未弛豫）接成一条能直接提交的计算流程，并补上续算审计 6.3 节列出的缺口：INCAR、KPOINTS、POTCAR、提交脚本、参考相、能量收集、形成能与判型分析。**本次没有运行 VASP，没有任何新能量**；脚本只用合成数据做过端到端测试，见第 6 节。

**2026-10-23 进展**：参考相 6 个已算完；母胞 O_ordered（hanhai 301880）在弛豫，D_seed1729（301881）排队。1b 的决定和成本估算见第 2.1 节。

## 1 要回答的问题与文献做法

师兄的要求是“和文献一样补上计算：退火前晶格未平整化时的主要缺陷类型，和退火后晶格平整化时的主要缺陷类型”。文献为 Evang 等，Adv. Mater. 34, 2106868 (2022)，第 5–10 页的做法已在 10 月 7 日的方案里逐页核对过。本流程逐项对应如下：

| 文献 | 本流程 |
|---|---|
| 无序岩盐：Te 占阴离子位，Pb/Sb/空位在阳离子位随机分布（Pb₂Sb₂Te₅ 含 20% 阳离子空位） | **D**：4 个随机种子母胞，144 原子（另含 16 个结构空位） |
| 有序：保持立方堆垛，沿 [111] 分成 Sb–Pb–Pb–Sb–空位层 | **O**：同一母格、同一组成 |
| —— | **V**：空位已成层但 Pb/Sb 仍混排。用来分开“空位成层”和“阳离子分层”两个效应（文献第 7 页认为两者都起作用） |
| 母胞先优化到理论密度；缺陷在母胞晶格固定下只弛豫原子 | 母胞 ISIF=3，缺陷 ISIF=2 |
| PBE+D3，无序 2×2×2、有序 2×2×1 k 网格 | PBE-D3(BJ)，`IVDW=12`；三族统一用 2×2×2（`--shift` 可改为 MP 偏移网格） |
| 12 类点缺陷，全部中性计算，用电子计数和 PDOS 判断施主/受主 | 同样 12 类；`analyze` 先做中性形成能和电子计数；带电态与电中性作为第二层（可选） |
| 无序 194 原子 × 16 母胞 × 8 位点 | 144 原子 × 4 母胞 × 4 位点；先做预筛，排序不稳定再扩大采样 |

文献结论（待本体系检验）：无序时 Pb 填空位、Sb 填空位和 Sb_Pb 容易形成，它们都是施主，所以是 n 型；有序后这三者代价上升，受主 Pb_Sb 占优，变为 p 型。

## 2 分期（`任务分期.csv`）

| 阶段 | 内容 | 作业数 | 回答什么 |
|---|---|---:|---|
| 0 | 参考相：Pb、Sb、Te、PbTe、Sb₂Te₃、9 原子层状 Pb₂Sb₂Te₅ | 6 | 化学势，以及层状相相对 PbTe+Sb₂Te₃ 的能量 |
| 1a | 9 个母胞弛豫 + 静态 | 9 | D/V/O 每个化学式单元的能量差，即有序化驱动力；弛豫后是否仍保持目标占位 |
| **1b** | **D×4 和 O×1 上的 5 类关键缺陷**：Sb_Pb、Pb_Sb、Pb_on_Vac、Sb_on_Vac、V_Te（D 每类 4 位点，O 取全部不等价位点） | **87** | **与文献相同的核心对比：退火前（D）和退火后（O）谁最低，是施主还是受主** |
| 2 | V×4 上的同 5 类 | 80 | 区分空位成层和 Pb/Sb 分层各自的作用 |
| 3 | 其余 7 类：V_Pb、V_Sb、Te_on_Vac、Pb_Te、Sb_Te、Te_Pb、Te_Sb | 67 | 排除被遗漏的低能者（旧数据中 Sb_Te 曾在贫 Te 条件下最低）。**做完这一步才能说“整体主导缺陷”** |
| 4（可选） | 竞争缺陷最低位点的带电态 + 母胞 DOSCAR | 依 1b–3 的结果而定 | 冻结缺陷电中性下的 E_F、n−p；先确认有序母胞是否开带隙 |

建议先用 1 个母胞和 2–3 个缺陷在集群上试跑，确定核数和耗时，再批量提交。

### 2.1 1b 的决定（2026-10-23）

- **缺陷做弛豫（ISIF=2），不做单点。** 缺陷是从未弛豫模板放入的，反位和填空位的弛豫能一般在 0.3–1 eV，与要分辨的 D–O 差别（文献中 Sb_Pb 随有序化上升 0.7–0.9 eV）同量级；“填入空位的原子是否移出原位 1 Å”也只能从弛豫后的结构判断。
- **分两批。** 87 个作业 = 4 个 D 母胞 × 5 类 × 4 位点（80）+ O 的 7 个。
  - **第一批 27 个**：O_ordered 的 7 个 + D_seed1729 的 20 个。`defects --families DO` 只为 `STATUS=DONE` 的母胞生成作业，所以这两个母胞算完后直接运行即可，不需要额外筛选。
  - **第二批**：其余 3 个 D 母胞（60 个），只有在第一批中 D 与 O 的差别没有明显超出 D 内部位点之间的离散时才补。
- **弛豫阶段放宽收敛**：`INCAR.defect_relax` 的 EDIFF 由 1E-6 改为 1E-5，EDIFFG 由 −0.02 改为 −0.03。最终能量仍取 `static/` 单点（EDIFF 1E-6、`LREAL=.FALSE.`），不受影响。母胞设置不变。
- **成本**：母胞 O_ordered 的前两个离子步共 64 min（每步约 17 个电子步）。这是最贵的情况：从未弛豫结构起步、ISIF=3、含启动开销。缺陷从已弛豫母胞出发，后续离子步一般只需 5–8 个电子步、共 15–40 步。单个缺陷的耗时等母胞第 5–10 步的实测值再定，不要用前两步外推。
- **集群**：13_ 全部改在 hanhai 上跑（jzzn 排不上）。所有生成命令加 `--cluster hanhai --ntasks 64 --vasp <hanhai 的 vasp_std>`，用 `inputs/submit.sh.hanhai.template`（CPU-256C768GB 分区，墙钟上限 15 天，默认 `--time 7-00:00:00`、`--mem 300G`）。它和 jzzn 模板的计算逻辑相同：relax（不收敛则续算）→ static → 写 `STATUS`。参考相此前在 jzzn 用 vasp 6.4.3 算过，要在 hanhai 用同一个 vasp 重算，保证形成能里的能量来自同一版本。
- **O_ordered 的补救**（hanhai 301880 用的是只有一次 mpirun 的脚本，没有 relax/static/STATUS）：作业结束后（收敛或被墙钟终止都一样），把 VASP 输出连同 INCAR 移进 `relax/`，再用 hanhai 模板生成的 submit.sh 提交。脚本检查 `relax/OUTCAR`：已收敛就只做 static，未收敛就从 CONTCAR 续算。
- **母胞弛豫步长**（2026-10-23）：O_ordered 用默认 POTIM（IBRION=2 时为 0.5）弛豫，前 9 步正常（能量 −583.66 → −589.79 eV，体积 4767 → 4995 Å³），第 10–12 步体积基本不变，能量却升到 −90.6 eV。`INCAR.parent_relax` 已加 `POTIM = 0.1`。
  - **原因已由单点对照确认**（hanhai 302292/302293，同一 INCAR、NSW=0）：第 9、10 步构型的单点 F 分别为 −589.960、−589.993 eV，而弛豫过程中第 10 步给出 −508.741 eV。两个构型体积只差 0.7%、最短键长正常（Pb–Te 3.19、Sb–Te 3.00 Å），所以结构没坏，是 ISIF=3 下晶胞变形后基组和 `LREAL=Auto` 的实空间投影没有重新优化造成的计算误差。`INCAR.parent_relax` 因此改为 `LREAL = .FALSE.`（80aaa99）。（早先从 OUTCAR 抓的体积表错位两步：OUTCAR 开头会先打印 `volume of cell`。以 ase 读 XDATCAR 为准。）
  - **收敛后要看压强**：母胞弛豫收敛后，检查 `static/OUTCAR` 的 `external pressure`。绝对值超过约 5 kbar 时，从 `relax/CONTCAR` 再弛豫一轮（体积变化大时，弛豫末尾的应力仍带基组误差）。
  - O_ordered 从第 9 步结构重新开始，D_seed1729 从模板重新开始；两个作业的 `relax/` 都要清空，否则模板会沿用旧的 `relax/INCAR`（不含 POTIM），或者从发散后的 CONTCAR 续算。
- **母胞分段弛豫**（2026-10-25）：D_seed1729 在 POTIM=0.1（第 34 步）和 0.05（重启后第 26 步）下仍出现能量突然上升，POTIM=0.02 才稳定。根源和 O_ordered 一样：ISIF=3 下体积变化后基组不再匹配晶胞，误差随步数累积。`INCAR.parent_relax` 改为 `POTIM = 0.02`、`NSW = 20`，生成母胞时加 `--restarts 15`：每 20 步从 CONTCAR 重启一次，基组随之重建，最多 300 步。模板在每段开始前删除 STOPCAR，所以运行中若发现能量上升（dE > +0.5 eV），可以写 `LSTOP = .TRUE.` 到 `relax/STOPCAR`，让 VASP 在当前离子步后停下，下一段自动从 CONTCAR 续算。
- **换行符**：模型包 `补算实施方案_20261007/rocksalt_*` 中的 POSCAR 是 CRLF 换行，VASP 拒读。`write_job` 现在写出时统一转成 LF（此前母胞作业因此失败过两次）。

## 3 操作

所有命令在本目录下执行，`WORK` 为任意新工作目录。脚本依赖 numpy、scipy、matplotlib；生成参考相需要 ase。

**集群上的目录布局**：脚本按“`scripts/` 的上两级 = 仓库根目录”找依赖，所以集群上要保持仓库的相对路径（本目录可以改名）：

```
<根>/<本目录>/inputs、scripts
<根>/补算实施方案_20261007/rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED/   # parents、defects 都要（约 5 MB）
<根>/补算实施方案_20261007/scripts/regenerate_defects_from_relaxed_parent.py   # defects 要（只依赖 numpy）
<根>/A2B2Te5_SOC_dense/materials/*/POSCAR、POTCAR             # POTCAR 拼接与 refs 要
```

只有 `refs` 需要 ase；集群上没有 ase 时，可以在有 ase 的机器上运行 `refs --cluster hanhai`，再把生成的作业目录拷过去。

```bash
# 0+1a：生成母胞与参考相作业（POTCAR 由仓库内同版本 PAW 拼接：Pb_d、Sb、Te）
python scripts/dvo_workflow.py parents --out WORK
python scripts/dvo_workflow.py refs    --out WORK
#   hanhai：以上及下面每条生成命令都加 --cluster hanhai --ntasks 64 --vasp <hanhai 的 vasp_std>
#   上传后在每个作业目录执行 sbatch submit.sh；脚本自动 relax（不收敛则续算，至多 3 轮）→ static，写 STATUS=DONE

# 1b：母胞 STATUS=DONE 后，用原工具从 relax/CONTCAR 重建缺陷，再写作业（默认 5 类关键缺陷、≤4 位点、q=0）
python scripts/dvo_workflow.py defects --parents WORK/parents --out WORK --families DO
# 2：再加 V 族；3：加其余 7 类
python scripts/dvo_workflow.py defects --parents WORK/parents --out WORK --families V
python scripts/dvo_workflow.py defects --parents WORK/parents --out WORK --types V_Pb V_Sb Te_on_Vac Pb_Te Sb_Te Te_Pb Te_Sb

# 收集与分析（任何阶段都可以运行，只统计 STATUS=DONE 的作业）
python scripts/dvo_workflow.py collect --root WORK --out WORK/energies.csv
python scripts/dvo_workflow.py analyze --energies WORK/energies.csv --out WORK/analysis --T-anneal 473 573

# 4（可选）带电态：charges.csv 每行写 parent_id,defect_id,"q1 q2"，例如 D_seed1729,Sb_on_Vac__s01,"1 2 3"
python scripts/dvo_workflow.py defects --parents WORK/parents --out WORK --charges charges.csv
python scripts/dvo_workflow.py analyze --energies WORK/energies.csv --out WORK/analysis --dos
```

常用选项：`--no-potcar`（不写 POTCAR，到集群上再拷）；`--ntasks/--nodes/--partition/--vasp`（默认沿用 jzzn 的 `cpu192` 与 vasp.6.4.3）；`--spin auto`（奇数电子缺陷用 ISPIN=2，用于抽查自旋影响）；`--mu mu.json`（手工指定化学势）。

**Pb₂Bi₂Te₅ 对照**（方案 B3）：先用原生成器生成 Bi 模型，再把 `--models` 指向它：

```bash
python ../补算实施方案_20261007/scripts/build_rocksalt_models.py --material Pb2Bi2Te5 --a 6.3 --repeat 4 4 1 --out BI_MODELS
python scripts/dvo_workflow.py parents --models BI_MODELS --out WORK_BI
python scripts/dvo_workflow.py refs --material Pb2Bi2Te5 --out WORK_BI   # Pb、Te、Bi、PbTe、Bi2Te3、层状 9 原子
```

## 4 计算设置（`inputs/`）

| 项 | 取值 | 理由 |
|---|---|---|
| 泛函 / vdW | PBE + D3(BJ)（`IVDW=12`） | 与文献一致，也与本项目 SOC 密网格一致 |
| POTCAR | Pb_d、Sb、Te（Bi 用 Bi_d），从仓库 `A2B2Te5_SOC_dense/materials/*/POTCAR` 拼接 | 与已有计算使用同一组数据集 |
| ENCUT | 360 eV | 与 Pb₂Sb₂Te₅ 密网格相同；约为最大 ENMAX 的 1.5 倍，ISIF=3 时可抑制 Pulay 应力 |
| 母胞 | ISIF=3，EDIFFG=−0.02，ISYM=0，POTIM=0.02，LREAL=.FALSE.，NSW=20 分段（`--restarts 15`） | 文献要求优化到理论密度；关闭对称，避免把 O 母胞锁在理想对称下 |
| 缺陷 | ISIF=2，晶格取母胞弛豫后的值；弛豫阶段 EDIFF=1E-5、EDIFFG=−0.03 | 与文献一致；形成能只做固定晶格比较；能量取自 static/ 单点 |
| 能量 | 统一由 `static/` 的 `energy(sigma->0)` 给出，`LREAL=.FALSE.`，`LORBIT=11` | 母胞与缺陷用同一套设置；DOSCAR 供第 4 阶段使用 |
| 展宽 | ISMEAR=0，SIGMA=0.05 | 窄隙或半金属宿主 |
| 自旋 / SOC | 默认 ISPIN=1，无 SOC | 文献级预筛；对决定排序的缺陷，再用 `--spin auto` 和 SOC 单点抽查 |
| 化学势 | 默认按 PbTe+Sb₂Te₃ 两相平衡取 Te-rich 和 Te-poor 两端（Pb₂Sb₂Te₅ 处在这条伪二元线上） | D/V/O 三族用同一组 μ。同种缺陷在不同族之间的能量差与 μ 无关 |

## 5 `analyze` 输出与读法

- `形成能_逐位点.csv`：每个位点、每种化学势条件下的 Ef；带电态还附带势对齐 ΔV 和像电荷修正。
- `形成能_家族统计.csv`：族 × 缺陷类型的位点数、最小值、中位数、最大值，以及构型平均 Ef_eff(T) = −kT ln⟨e^(−Ef/kT)⟩。另外记录“填入空位的原子是否移出原位 1 Å 以上”。
- `DVO_结论摘要.json`：
  - `dominant_by_condition_family_T`：各族 Ef_eff 最低的三类缺陷，以及按电子计数（Pb²⁺、Sb³⁺、Te²⁻，完全电离）得到的施主−受主净浓度，用来判断 n/p 倾向。这就是文献层级的“主要缺陷类型”判据。
  - `same_type_shift_between_families`：同一类缺陷 O−D、V−D 的 Ef 差。它与化学势无关，最适合直接对照研究记录中摘录的文献结论（Sb_Pb 形成能随有序化上升约 0.7–0.9 eV；Pb/Sb 填空位随空位成层变贵）。
  - `host_energy_per_formula_unit`：D、V 相对 O 每个化学式单元的能量差，即有序化驱动力。
  - `charged_neutrality`（第 4 阶段）：退火温度下平衡，冻结缺陷总量后在 300 K 重新求解，给出 E_F 和 n−p。只统计带电态的母胞；|n−p| 低于 1e17 cm⁻³ 时记为 undetermined。
- `DVO_中性形成能分布.png/.svg`：每个化学势条件一行，横轴为缺陷类型（括号内为电子计数电荷），D、V、O 三族的每个位点都画出。

**判读顺序**：先看 1a 的母胞，确认弛豫后仍是 D/V/O，并读出有序化能；再看 1b 中 D 与 O 各自最低的缺陷类型及其施主/受主属性；然后看 `same_type_shift` 是否重现“填空位施主变贵、Pb_Sb 相对变便宜”；最后用阶段 3 排除其他类型。如果排序差异小于位点之间的离散，结论应写成“在所采样模型中未分辨”。

## 6 已验证与未验证

- 已验证（仅软件层面）：在临时目录中，用玩具能量伪造 OUTCAR/DOSCAR，完整走通 `parents → refs → defects（234 个中性 + 11 个带电）→ collect（260）→ analyze（含冻结电中性）`。缺陷几何由原 `regenerate_defects_from_relaxed_parent.py` 生成，NELECT、POTCAR 顺序、Bi 对照分支都已检查。测试数据没有提交。
- 未验证：任何真实能量、弛豫是否保持目标相、位点和尺寸收敛、SOC 影响。中性形成能不等于载流子浓度；无序母胞可能没有带隙，带电态的参考取自母胞自身费米能。静态计算不给出退火温度或动力学。144 原子小于文献的 194/432 原子。
