# SOC 机制检验：Bi 为什么把 Γ 处的价带态压下去（2026-10-31）

> **结果（2026-10-31 晚）**：解释成立，Γ 谷被压低完全来自 Bi 6p 的 SOC。
> - **无 SOC 时**：Γ 价带顶是 s 态，B 位 p 态都在导带。
> - **有 SOC 时**：p₁/₂ 态被压到 s 态之下（带反转），压低的量 SnSb 262 meV、SnBi 714 meV。
> - **Sb → Bi 的拆分**：Sn 系 −194 = 标量 +257 + SOC −452；Pb 系 −177 = +291 − 467（meV）。
> - **晶格交换**：只差约 25 meV。
>
> 见 `analysis.log`、`结果说明_20261031.md` 末节，以及论文报告第 3 节与图 S2。旧版脚本按“B 位 p ≥ 25%”选带，在无 SOC 时选错了态，已改为按轨道成分追踪。

## 要检验的解释

Γ 处的价带态约 40% 是 B 位 p，p_z : p_xy ≈ 1 : 2.5，接近 j = 1/2 的 p₁/₂ 态（1 : 2）。孤对电子谷（Σ、R）的 B 位 p 只有 1–3%。Bi 6p 的自旋–轨道耦合比 Sb 5p 强得多，会把这个 p₁/₂ 型态压下去，孤对电子谷不受影响，于是 Γ 相对孤对电子谷下降。

现有 PBE+SOC 中，这一下降为 Sn 系 194 meV、Pb 系 177 meV；HSE 中为 226、197 meV。

## 两组计算（都是 PBE，6 个作业，每个约 1 小时）

| 目录 | 内容 | 与现有 `<材料>/c2_pbe` 的差别 |
|---|---|---|
| `d1_nosoc_Pb2Sb2Te5`、`d1_nosoc_Pb2Bi2Te5`、`d1_nosoc_Sn2Sb2Te5`、`d1_nosoc_Sn2Bi2Te5` | ① 关掉 SOC | INCAR 去掉 `LSORBIT`、`LNONCOLLINEAR`、`SAXIS`、`MAGMOM`，`NBANDS` 减半；用 `vasp_std`。其余逐字相同 |
| `d2_Sn2Bi2Te5_in_Sn2Sb2Te5_lattice` | ② SnBi 的成分放进 SnSb 的晶格和原子坐标 | INCAR 与 SnBi 的 c2_pbe 相同（含 SOC）；POSCAR 用 SnSb 的晶格和坐标，只把元素行改成 `Bi Te Sn` |
| `d2_Sn2Sb2Te5_in_Sn2Bi2Te5_lattice` | ② 反过来 | INCAR 与 SnSb 的 c2_pbe 相同；POSCAR 用 SnBi 的晶格和坐标，元素行 `Sb Te Sn` |

**共同设置**
- `KPOINTS`：Gamma 12×12×2；
- `KPOINTS_OPT`：c2 路径，204 点（`c2_path/KPOINTS_OPT`）；
- 从头自洽，`ISTART = 0`，不需要 WAVECAR。

**POSCAR 与 POTCAR**
- **POSCAR**：①用 `A2B2Te5_SOC_dense/materials/<材料>/POSCAR`。已和 c2_pbe 的 vasprun.xml 核对过，晶格与坐标一致。
- **POTCAR**：用成分对应材料的 a_pbe 的 POTCAR，元素顺序与 POSCAR 第 6 行一致（B 位、Te、A 位）。② 的 POTCAR 跟成分走，不跟晶格走。

## 检查与输出

- **检查**：
  - 自洽收敛，有 `aborting loop because EDIFF is reached`；
  - 有 `General timing`；
  - PROCAR_OPT 含 204 个 k 点。
- **NCORE**：c2_pbe 的 INCAR 写了 `NCORE = 4`，原样保留。若 VASP 因 KPOINTS_OPT 拒绝 NCORE，就删掉这一行（运行记录第 6 节），并在汇报中注明。
- **推送**（每个目录）：
  - `PROCAR_OPT.gz`（必须）；
  - `OSZICAR`；
  - `能量.txt`：energy(sigma->0)、E-fermi、EDIFF 行；
  - `vasprun.xml.gz`：不超过 50 MB 时才推；
  - `out.log` 尾部。

## 分析

```bash
python 10_脚本/soc_mechanism.py
```

脚本对每个计算读出三个能量（都相对 Σ 谷）：
- Γ 价带顶；
- Γ 处“B 位 p 态”：B 位 p ≥ 25% 的最高价带，同时给出它的 p_z 比例；
- Γ 导带。

然后算两个量：
- **双重差分**：Sb → Bi 的变化，有 SOC 减无 SOC。A 位和 Te 的 SOC 在两种材料中相同，相减后基本抵消，剩下的是 Bi 比 Sb 多出来的 SOC 的贡献。
- **晶格交换**：同一晶格中只换 B 位元素的变化。

**判读**

| 结果 | 结论 |
|---|---|
| 无 SOC 时 Sb → Bi 的下降基本消失，且晶格交换后下降仍在 | 主要是 Bi 6p 的 SOC |
| 只消失一部分 | SOC 是主因，其余来自 Bi 6s、6p 能级或体积，由 ② 拆分 |
| 基本不变 | SOC 解释不成立 |

**注意**：关掉 SOC 后，Γ 处的能级可能重新排序，反转甚至可能消失。所以要按轨道成分（B 位 p、p_z 比例）追踪同一个态，不能按带序号比较。脚本已经这样做。
