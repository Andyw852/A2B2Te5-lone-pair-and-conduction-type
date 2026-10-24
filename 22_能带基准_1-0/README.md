# 1-0 能带基准：PBE+SOC 与 HSE+SOC（Pb₂Sb₂Te₅）

2026-10-24 完成。这是问题一（孤对电子与电输运）的**能带基准**：此前所有依赖带隙的结论
（低浓度 S、双极热导、n 型变号阈值、剪刀修正）都以旧 HSE 的 0.065 eV 为依据，
而旧 HSE 的 INCAR 用了 `ISYM=0` + 不可约 k 点列表（密度未对称化，计算错误）。本页给出正确结果。

## 计算设置

| 项 | 作业 b（PBE+SOC） | 作业 c1（HSE+SOC） |
|---|---|---|
| 作业号 | jzzn `ke_1_0_bandgap/.../b_pbe_parity2` | clus_xa `72279386` → `c1_scf` |
| 泛函 | **PBE**（基线 INCAR） | **HSE06**（`LHFCALC=.TRUE.`、`AEXX=0.25`） |
| SOC | `LSORBIT=.TRUE.`，`vasp_ncl` | 同左 |
| KPOINTS | 自动 Γ **12×12×2** → **38 个不可约点** | 同左（同一网格，保证可比） |
| NBANDS | 96（NELECT=68） | 96 |
| ENCUT | 360 | 360 |
| 对称性 | 默认（`ISYM=2`） | **`ISYM=3`**（与 b 的不可约集相同） |
| HSE 加速 | — | `PRECFOCK=Fast`、`NKREDX/Y=2`、`NKREDZ=1` |
| 墙钟 | — | 21.2 h（33 个电子步） |

## 结果

| 方法 | VBM (eV) | CBM (eV) | **带隙 (eV)** | 带边 k 点 | 类型 |
|---|---:|---:|---:|---|---|
| **PBE+SOC**（b） | 6.2713 | 6.4980 | **0.2268** | VBM k#3、CBM k#1 | **间接** |
| **HSE+SOC**（c1） | 6.2370 | 6.2883 | **0.0513** | 均为 k#1 | **直接（Γ 点）** |

VBM/CBM 按能带序号读取：SOC 下每条带装 1 个电子，NELECT=68 → 第 68 条为价带顶、第 69 条为导带底。

## 三条结论

### 1. HSE 把 PbSb 的带隙压到 51 meV，只有 PBE 的 1/4.4

```
PBE+SOC  0.2268 eV
HSE+SOC  0.0513 eV     比值 0.226
```

### 2. PBE 给间接带隙，HSE 给 Γ 点直接带隙

PBE 的 VBM 在 k#3、CBM 在 k#1；HSE 的 VBM 与 CBM 都在 Γ（k#1）。
这解释了 12_ 里“PbSb 接近 s–p 反转”的描述——HSE 下带边确实汇聚到 Γ 点。

### 3. 旧 HSE 的 0.065 eV 数值上接近，但计算本身是错的

旧 HSE 用 `ISYM=0` 配不可约 k 点列表，电荷密度与交换算符未按对称性补全。
新 HSE（`ISYM=3`、与 PBE 同一 38 点网格）给出 **0.0513 eV**。
数值相近纯属巧合，**不能再用旧值做依据**。

## 对已有结论的影响

| 结论 | 原依据带隙 | 新 HSE 带隙 | 影响 |
|---|---:|---:|---|
| 低浓度 S 由带隙决定、p 型 S 峰值 | 0.065 | **0.0513** | 双极效应更强，S 峰值更低 |
| n 型 S 变号阈值（原 3.0e18） | 0.065 | 0.0513 | 阈值会**更小** |
| 双极热导分解（Lorenz 5.3 → 1.6） | 0.065 | 0.0513 | 双极贡献会**更大** |
| **AMSET 的剪刀目标（scissor）** | 0.065 | **0.0513** | **需要更新输入后重跑或至少记录偏差** |

**方向性结论不变**：PbSb 仍是四种材料里带隙最小的，双极效应最显著。
但**绝对值需要按 0.0513 eV 复核**。

## 待补

1. **Bi 系的 HSE**（PbBi、SnBi）：按 19_ 1-0 的方案，用 12×12×2 + `PRECFOCK=Normal` + KPOINTS_OPT；
   注意 PbBi/SnBi 的旧 HSE 路径用了零权重 k 点，带隙不可靠，必须重算。
2. **mBJ** 暂搁置（VASP 报 `XC_FOCK_READER: Only one tag`，mBJ 与 `GGA` 冲突，且 `STOP 1` 原因未定）。
3. 带边 k 点的 HSE 带隙（1-8 里的“密网格带边点”）——本轮只给了 38 点网格上的值。

## 数据位置

- PBE+SOC：jzzn `~/ke_1_0_bandgap/Pb2Sb2Te5/b_pbe_parity2/`（EIGENVAL 137 KB、WAVECAR 269 MB）
- HSE+SOC：clus_xa `~/autozt_work/1_0/Pb2Sb2Te5/c1_scf/`（EIGENVAL 137 KB、vasprun.xml 14 MB）