# 执行进展：密网格 AMSET 的 wave 步（2026-12）

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
