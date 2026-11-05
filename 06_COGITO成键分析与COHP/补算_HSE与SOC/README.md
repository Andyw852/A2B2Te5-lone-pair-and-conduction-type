# COHP 补算：HSE06 与 SOC（2026-11-05）

**目的**：报告图 3-2 的 COHP 目前只有 PBE、不加 SOC（COGITO）。用户要求用更高的层次检验“孤对电子是填满的阳离子 s–Te p 反键态”这一结论。

- **f**：COGITO + HSE06（不加 SOC）。COGITO 的文档要求 `LSORBIT = .FALSE.`，所以 COGITO 只能做到这一步。
- **g**：LOBSTER，PBE 不加 SOC → PBE+SOC → HSE06+SOC。先确认 LOBSTER 接不接受 SOC 的波函数，再算 HSE06+SOC。

两组互相独立，可以同时交。四种材料都要算；如果排队紧张，先算 **Sn₂Bi₂Te₅、Sn₂Sb₂Te₅**。

| 目录 | 程序 | 步骤 | 大致耗时（每种材料） |
|---|---|---|---|
| `f_cogito_hse/<材料>/` | VASP（vasp_std）+ COGITO | ① PBE 不加 SOC 自洽 → ② HSE06 不加 SOC 自洽，各跑一遍 COGITO | ① 约 1 小时；② 几小时到一天；COGITO 每次约 25 分钟 |
| `g_lobster_soc/<材料>/` | VASP（vasp_std / vasp_ncl）+ LOBSTER | ① PBE 不加 SOC → ② PBE+SOC →（可选）③ HSE06+SOC，每步之后跑 LOBSTER | ①② 各约 1–3 小时；③ 很贵（全 k 网格、HSE+SOC），可能要几天 |

**通用约定**
- **结构**：`POSCAR` 取自 `A2B2Te5_SOC_dense/materials/<材料>/POSCAR`，与 c2 能带、BoltzTraP2 以及现有 COGITO 结果是同一个结构。不要换。
- **POTCAR**：从 `<材料>/c2_hse` 的作业目录复制。元素顺序同 POSCAR 第 6 行（B 位、Te、A 位）；Pb、Sn、Bi 用带 d 的势（`Pb_d`、`Sn_d`、`Bi_d`），Sb、Te 不带 d。NELECT 应为 68（Sb 系）、88（Bi 系）。
- **k 网格**：Γ 12×12×2，与 c2 自洽相同。
- **集群环境**：`submit*.sh` 里的分区、模块、vasp 路径、conda 环境照抄 `SOC机制检验_HSE/e1` 与 `06_COGITO成键分析与COHP/运行日志/<材料>/submit.sh`（jzzn）。换机器就改这几行。
- **AutoZT**：可以用 AutoZT 的 `skill/cohp-cogito` 跑 COGITO 三件套，但上游 VASP 必须用本目录的 INCAR。不要改 AutoZT 的代码。
- **检查不过时先报告**，不要自行改参数重算。

---

## f：COGITO + HSE06（不加 SOC）

每个材料目录内先交 `submit_vasp.sh`，跑完再交 `submit_cogito.sh`。

| 步 | INCAR | 内容 | 与已有计算的差别 |
|---|---|---|---|
| 1 | `INCAR.1_pbe` | PBE 不加 SOC 自洽，`vasp_std` | 由 `<材料>/c2_pbe/INCAR` 去掉 SOC（`LSORBIT`、`LNONCOLLINEAR`、`SAXIS`、`MAGMOM`）；`NBANDS = 144`、`ISYM = 2`、`LWAVE = .TRUE.`，即原来 COGITO 的上游设置 |
| 2 | `INCAR.2_hse` | HSE06 不加 SOC 自洽，从第 1 步的 WAVECAR 续算 | 与 `22_能带基准_1-0/SOC机制检验_HSE/e1_hse_nosoc_<材料>/INCAR.2_hse_scf` 相同，只改 `NBANDS = 144`、`SIGMA = 0.01`（与第 1 步一致）、`LWAVE = .TRUE.`；`ISYM = 3`、`NKRED 2 2 1`、`PRECFOCK = Fast` 照旧 |

`submit_vasp.sh` 每步之后把 COGITO 要读的文件（POSCAR、POTCAR、KPOINTS、INCAR、OUTCAR、vasprun.xml、WAVECAR、IBZKPT）复制到 `cogito_pbe/`、`cogito_hse/`。`submit_cogito.sh` 在这两个目录里依次跑 `COGITO`、`COGITOanalyze`、`COGITOpost`。

第 1 步顺便是对照：设置与原来的 PBE COGITO 相同，结果应当复现 `06_COGITO成键分析与COHP/结果/<材料>/bond_info.txt`。

用 48 核：VASP 会把 NBANDS 补成并行组数的整数倍，144 能被 48 整除，HSE 一步（不设 NCORE）就不会被补到 160。

**检查**
- 两步都有 `General timing` 和 `aborting loop because EDIFF is reached`（脚本已检查）。
- 两个 OUTCAR 中 `NBANDS = 144`；`IBZKPT.1` 与 `IBZKPT.2` 相同（c2 自洽为 38 个不可约 k 点）。
- 不加 SOC 时 Γ 处不应有带反转：Γ 处的价带顶是阳离子 s 态，B 位 p 态在导带。与 PBE 的 `SOC机制检验/d1` 相同。
- COGITO：`percent charge spilling` 应在 1% 以内（原来的 PBE 运行为 0.3%）。所有 `Warning`、`ERROR` 原样记下。原来的 PBE 运行也报过 “NRMSE … much larger in later step”，不是致命错误，但要报告。
- 对照：`cogito_pbe/bond_info.txt` 中最短的 A–Te、B–Te 键的 ICOHP，与 `结果/<材料>/bond_info.txt` 相差应小于 0.03 eV/键。超过就先报告。
- 如果 COGITO 读不了 HSE 的输出（报错、或 charge spilling 明显变大），原样报告错误，不要改 COGITO 或 AutoZT。

**交回**（每个材料）
- `cogito_pbe/`、`cogito_hse/` 中 COGITO 三件套写出的全部文件（不要 WAVECAR）。必须包括：
  - `bond_info.txt`、`all_atoms.json`、`all_unique_bonds.json`；
  - COGITOpost 生成的 COHP 交互图（`*.html`；我用 `10_脚本/cohp_extract.py` 从里面取轨道分辨的 COHP）；
  - `cogito_run.log`、`analyze.log`、`post.log`、`orb_converg_info.json`。
- `OUTCAR.1`、`OUTCAR.2` 末尾 200 行，`OSZICAR.1/2`，`IBZKPT.1/2`，`out.log.1/2`。

---

## g：LOBSTER，PBE → PBE+SOC → HSE06+SOC

### 第 0 步：先查 LOBSTER（不用算）

1. 集群上有没有 LOBSTER，哪个版本：`command -v lobster`、`ls /public/home/wangchao/software | grep -i lobster`、`module avail 2>&1 | grep -i lobster`。运行一次不带输入的 `lobster`，开头会打印版本号。
2. 手册（LOBSTER 安装包里的 PDF 或 README）中有没有讲非共线 / 自旋–轨道耦合（noncollinear、spin–orbit、LSORBIT）。原样摘录相关的句子。
3. 没有 LOBSTER，或手册明确说不支持非共线：只跑第 1 步（作同一程序内的对照，与 COGITO 互相印证），第 2、3 步不跑，报告即可。
4. 手册没写清楚：照常交 `submit.sh`。脚本在第 2 步之后先试 LOBSTER，它不接受 SOC 的 WAVECAR 就会停下，不会白跑第 3 步。

### 三步

| 步 | INCAR | 内容 | 与已有计算的差别 |
|---|---|---|---|
| 1 | `INCAR.1_pbe_nosoc` | PBE 不加 SOC，`vasp_std` | 由 `c2_pbe/INCAR` 去掉 SOC；`ISYM = 0`、`NBANDS = 64`、`LWAVE = .TRUE.` |
| 2 | `INCAR.2_pbe_soc` | PBE+SOC，`vasp_ncl` | 与 `c2_pbe/INCAR` 相同，只改 `ISYM = 0`、`NBANDS = 128`、`LWAVE = .TRUE.` |
| 3（可选） | `INCAR.3_hse_soc` | HSE06+SOC 自洽，`vasp_ncl`，从第 2 步的 WAVECAR 续算 | 与 `c2_hse/INCAR` 相同，只改 `ISYM = 0`、`NBANDS = 128`、`SIGMA = 0.01`、`LWAVE = .TRUE.` |

- **ISYM = 0 与 NBANDS**：LOBSTER 要用全布里渊区的波函数（只用时间反演），所以 `ISYM = 0`。若 LOBSTER 报错说 k 点不全或对称性不对，改成 `ISYM = -1` 重跑该步，并在报告里写明。NBANDS 不少于基组函数数：Sb 系 46（不加 SOC）、92（SOC），Bi 系 56、112。
- **`lobsterin`**：
  - 基组 `pbeVaspFit2015`：Pb 6s 6p 5d、Sn 5s 5p 4d、Bi 6s 6p 5d、Sb 5s 5p、Te 5s 5p。
  - 取 2.95–3.30 Å 内的全部键，即阳离子–Te 键（3.01–3.26 Å）；Te–Te 最短 3.58 Å，不会进来。
  - `orbitalWise` 给出轨道分辨的 COHP（s–p、p–p 等）。若安装的版本不认 `cohpGenerator … orbitalWise` 这种写法，按手册改成等价写法，并报告改了什么。
- **运行方式**：`submit.sh` 每步 VASP 之后在原地跑 LOBSTER，结果放在 `lobster_1_pbe_nosoc/`、`lobster_2_pbe_soc/`、`lobster_3_hse_soc/`。
  - 第 2 步的 LOBSTER 失败时脚本停下，不跑第 3 步。
  - 第 3 步只在 `sbatch --export=ALL,RUN_HSE=1 submit.sh` 时才跑。不跑第 3 步时保留第 2 步的 WAVECAR，以后可以直接续跑。
  - 找不到 LOBSTER 时先 `export LOBSTER=<路径>` 再交。
- **第 3 步**很贵（HSE+SOC，`ISYM = 0`，12×12×2 只用时间反演约 148 个 k 点）。只在第 2 步的 LOBSTER 成功之后再跑；先跑 Sn₂Bi₂Te₅、Sn₂Sb₂Te₅。

**检查**
- 每步 VASP 收敛（脚本已检查），OUTCAR 中的 NBANDS 与 INCAR 一致（脚本会提示）。
- `lobsterout` 中 `abs. total spilling` 与 `abs. charge spilling` 都在 2% 以内（一般 < 1.5%）。原样记下数值。
- `ICOHPLIST.lobster` 中有阳离子–Te 键（3.01–3.26 Å）。第 1 步最短的 B–Te 键 ICOHP 应为负，且与 COGITO 同号、同一量级。两种程序的基组不同，数值不必相等。
- 第 2 步 LOBSTER 若报错，把 `lobster_2_pbe_soc/lobster.out` 与 `lobsterout` 原样交回，这就是“LOBSTER 不支持 SOC”的证据。

**交回**（每个材料）
- `lobster_*/` 整个目录（含 `ICOHPLIST.lobster`、`COHPCAR.lobster.gz`、`lobsterout`、`lobster.out`、`lobsterin`、INCAR、OUTCAR）。
- `OSZICAR.*`、`out.log.*`，以及第 0 步查到的 LOBSTER 版本和手册摘录。

---

## 收到后（Claude 做）

- f：用 `10_脚本/cohp_extract.py` 从 `cogito_hse/` 的交互图取轨道分辨 COHP。图 3-2 加 HSE 与 PBE 的对比；表中给出两种泛函下 s–p、p–p、总 ICOHP 与价带顶以下 1 eV 内的 s–p ICOHP；先用 `cogito_pbe/` 核对能否复现原结果。
- g：同一程序（LOBSTER）内比较 PBE、PBE+SOC、HSE06+SOC，看 SOC、HSE 改变 s–p 反键多少，并与 COGITO 的 PBE 结果交叉检验。
- 报告第 3 节“局限”与第 7 节“待补”相应更新。
