# 执行进展：密网格 AMSET 的 wave 步（2026-10-19）

配合 `README.md` 的 P1.1。本文件记录实际执行路径、踩到的坑和已产生的结果。

## 1 找到了密网格那批数据的 AutoZT 配置

之前 `az progress` 里看不到 A2B2Te5，是因为项目配置不在全局 `setting/`，而在项目自己的目录里：

`/mnt/d/tf_data/work_AutoZT/A2B2Te5_SOC_dense/project_setting/tf_A2B2Te5_SOC_dense.yaml`

```yaml
ke-dft-cpu:
  local_root: /mnt/d/tf_data/work_AutoZT/A2B2Te5_SOC_dense/materials
  work_dir:   /public/home/wangchao/ke_soc_dense_20261008
  hpc: jzzn
  max_jobs: 8
```

原 steps 只有 5 步：`step0_seed / step3_uniform / step8.1_boltztrap / step8.2_mass / step8.3_mass_fit`。
**没有 `step4_wave`、也没有 `step8_amset`**——这就是 `amset wave` 从未跑过的原因（不是失败，是没定义）。
集群上 `ke-dft-cpu/` 里的 `gen_dense_soc.py`、`run_bt2_soc.py` 等就是这个 YAML 的 `gen` 命令。

## 2 已补 step4_wave

原文件已备份为 `tf_A2B2Te5_SOC_dense.yaml.bak-before-s4`。新增：

```yaml
      - seq: 4
        name: step4_wave
        label: S4_wave_dense_soc
        src: step4_wave
        check: marker
        marker: "wavefunction.h5:"
        gen: "gen_wave_dense.py --material {mat}"
        gen_need: []
        needs: [step3_uniform]
```

已用 YAML 解析校验：6 步、依赖链 `step3_uniform → step4_wave` 正确。

**生成器** `project_setting/gen_wave_dense.py`：与技能版 `gen_step6_wave.py` 的唯一区别是
生成期**不检查 WAVECAR**，改由作业内对 `../step3_uniform/WAVECAR` 建相对软链
（密网格 WAVECAR 46–51 GB，本项目的 `fetch_files` 从不回传它；技能版会因
`src.is_file()` 直接退出）。

## 3 一个必须记住的坑：QOS 墙钟

首次提交一直 PENDING，原因是 `QOSMaxWallDurationPerJobLimit`：

| QOS | MaxWall |
|---|---|
| **normal**（账号默认） | **00:01:00** ← 1 分钟 |
| debug | 01:00:00 |
| **regular** | **1-00:00:00** |
| premium | 2-00:00:00 |

所以提交脚本必须显式写 `#SBATCH --qos=regular`，否则 8 h 的作业永远排不上。
这很可能也是项目历史上那次 `00:01:00` TIMEOUT 的真因（当时归因给了 `#SBATCH` 占位符问题）。
`submit_amset.tpl` 一类模板建议也补上这一行。

## 4 已完成的结果

| 项 | 值 |
|---|---|
| 作业 | 4226390（Pb₂Sb₂Te₅） |
| 状态 | COMPLETED，ExitCode 0:0 |
| 墙钟 | 1677 s ≈ **28 min** |
| 峰值内存 | **115 847 476 KB ≈ 110 GB** |
| 纳入带 | `Including bands 59—78`（20 条） |
| `wavefunction.h5` | 6 337 997 405 B ≈ **5.9 GB** |

Sn₂Sb₂Te₅ 已提交：作业 4227066（同一套脚本，带 `--qos=regular`）。

## 5 下一步

1. 等 SnSb 完成，核 `wavefunction.h5` 尺寸与带数。
2. `step8_amset` 还不能直接加：技能版它的 `needs` 要求
   `step5_dielect / step6_elastic / step7b_deform_read / step2 能带`，本项目一套都没有。
   试点阶段按既定方案**沿用 12×12×2 的 ε∞/ε₀/弹性/形变势**，需要为该 YAML 放宽 needs
   或另写一个 amset 步骤定义。
3. memory 110 GB 说明密网格 wave 很吃内存：后续材料与 AMSET 主步需要按此估资源。

## 6 energy_cutoff 与形变势带号（2026-10-19 实测更正）

**三个校验，前两个是硬报错，第三个要靠自己核：**

| 对象 | AMSET 校验什么 | 后果 |
|---|---|---|
| `wavefunction.h5` | 带数（`set_overlap_calculator`） | 不匹配直接报错 |
| `deformation.h5` | 带数（`Deformation potential file does not contain the correct number of bands`） | 不匹配直接报错 |
| 两者 | **带号** | **不校验**：带数恰好相同而带号错开时会静默算错 |

> 更正：此前写的「AMSET 根本不让你跑」只对一半。它校验**带数**，不校验**带号**。

**选带规则**（由 `10_脚本/amset_band_window.py` 验证）：某条带只要在
`[VBM - e, CBM + e]` 内有能量就保留（未被剪刀修正的 DFT 带边）。

**关键结论：用默认 `energy_cutoff = 1.5 eV` 时，密网格选带与旧 deformation.h5 逐个带号重合：**

| 材料 | e = 1.0（原做法） | e = 1.5（默认，旧计算口径） | 旧 deformation.h5 |
|---|---|---|---|
| Pb2Sb2Te5 | 59-78（20）不匹配 | **53-82（30）** 匹配 | 53-82（30） |
| Pb2Bi2Te5 | 77-94（18）不匹配 | **67-100（34）** 匹配 | 67-100（34） |
| Sn2Sb2Te5 | 61-80（20）不匹配 | **57-82（26）** 匹配 | 57-82（26） |
| Sn2Bi2Te5 | 81-96（16）不匹配 | **73-102（30）** 匹配 | 73-102（30） |

所以**不需要重做形变势（28 h），只需用 `-e 1.5` 重跑 wave（约 40 min/材料）**。
同时 `energy_cutoff` 回到与旧计算相同的 1.5 eV，新旧对比中唯一的实质变化就是 DFT 网格，
700-900 K 下 1.0 eV 偏紧的问题也随之消失。

**已执行**：取消三个用 `energy_cutoff=1.0` + 旧形变势的作业（4245803 / 4246175 / 4227148，
它们一启动就会报形变势带数错）；用 `-e 1.5` 重提 wave：Pb2Sb2Te5 = 4247041、Sn2Sb2Te5 = 4247042。
提交后须核对日志里的 `Including bands` 是否为 53-82 / 57-82。

**残余局限**：形变势仍来自 12x12x2（ADP 的 k 分辨率没变）。ADP 弱于 POP 3-17 倍，
先按局限注明；若密网格结果显示 ADP 变重要，再考虑在密网格上重做。

## 7 各步骤泛函一致性审计（2026-10-19）

| AMSET 输入 | 来源 | 泛函 |
|---|---|---|
| 能带 / 波函数 | 密网格 step3_uniform | PBE-D3(BJ) + SOC |
| 形变势 ADP | step7_deform（18 个应变单点） | PBE-D3(BJ) + SOC |
| 弹性常数 Cij | step6_elastic（IBRION=6, ISYM=2） | PBE-D3(BJ)，无 SOC |
| 介电 eps_inf / eps_0 / Born Z* / omega_LO | step5_dielect SOC DFPT（9x9x2 / v2） | PBE + SOC（DFPT 不含 D3） |
| POP 频率 | 同上 DFPT 的 Gamma 声子导出 | 同上 |
| 带隙 | 历史 HSE+SOC 剪刀值 | HSE + SOC |

**口径说明（更正此前"介电缺 D3、影响未量化"的误判）：**

> 介电与 Born 有效电荷由 DFPT（PBE+SOC）在 PBE-D3(BJ) 弛豫结构上计算；D3 为经典对势、
> 不参与电子响应，VASP 的 DFPT 亦不支持 vdW 修正（生成脚本 KEEP_D3_IN_DFPT=False 即按此设计），
> 故 eps_inf、Z* 不受影响。Gamma 点力常数（eps_0 离子部分与 omega_LO）未含 D3 的贡献，
> 对层内光学模影响预期为百分之几量级；**S5 无虚频**。

**Gamma 虚频核查（只读，四材料 DFPT OUTCAR）**：

| 材料 | 模态数 | 最小频率 / THz | 3 个近零模（声学） | 负频 | 最低光学模 / THz |
|---|---:|---:|---|---:|---:|
| Pb2Sb2Te5 | 27 | +0.0009 | 0.0009 / 0.0019 / 0.0068 | 0 | 0.888 |
| Sn2Sb2Te5 | 27 | +0.0015 | 0.0015 / 0.0017 / 0.0121 | 0 | 0.940 |
| Pb2Bi2Te5 | 27 | +0.0014 | 0.0014 / 0.0063 / 0.0173 | 0 | 0.747 |
| Sn2Bi2Te5 | 27 | +0.0036 | 0.0036 / 0.0057 / 0.0235 | 0 | 0.783 |

四材料均无虚频（无 f/i 标记），最小三个模即 Gamma 声学模（数值应严格为零）。
最低光学模与 18_ 的 lowest_optical_THz 一致（0.888 / 0.747 / 0.940 / 0.783），
说明"无 D3 力常数 + 带 D3 结构"没有在层间模上引入软化，eps_0 与 omega_LO 可用。

**代码注释建议（未改，属共享技能代码，待批准）**：
`AutoZT/skill/ke-dft-cpu/step5_dielect/gen_step8_dielect.py` 第 121 行附近注释
"只影响总能，不影响 eps" 应改为"不影响 eps_inf / Z*；对 Gamma 点力常数（进而 eps_0 离子部分、
omega_LO）的贡献在 DFPT 中被略去"。
