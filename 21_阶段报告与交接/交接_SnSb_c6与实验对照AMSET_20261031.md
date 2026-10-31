# 交接（2026-10-31）：③ SnSb c6_local　④ 实验对照 AMSET　⑤ 收尾 Sb 系 400–700 K AMSET

三项互相独立。③ 在 clus_xa（VASP），④⑤ 在 hanhai25（AMSET）。

## ③ Sn₂Sb₂Te₅ c6_local（必须做）

- **目录**：`22_能带基准_1-0/Sn2Sb2Te5/c6_local/`，共 202 点、6 批。设置与刚跑通的 SnBi c6 相同，详见目录内 `README.md`。
- **每批**：在 `bNN/` 里放本目录的 INCAR、KPOINTS 和 `KPOINTS_OPT_bNN`（改名为 `KPOINTS_OPT`）；POSCAR、POTCAR、WAVECAR 与 SnSb `c_hse_A` 相同。
- **资源**：64 核，`--time=0-12:00:00`。
- **检查**：
  - 自洽能量为 −41.37891726 eV；
  - 没有 `internal error`；
  - 有 `General timing`。
- **推送**：每批的 `PROCAR_OPT.gz`、`OSZICAR`、`out.log` 尾部。
- **分析**：

  ```bash
  python 10_脚本/c4_mirror_scan.py local-fit 22_能带基准_1-0/Sn2Sb2Te5/c6_local \
    --parts 22_能带基准_1-0/Sn2Sb2Te5/c6_local/b0{1,2,3,4,5,6} | tee 22_能带基准_1-0/Sn2Sb2Te5/c6_local/local_fit.log
  ```

- **验收**：校验点与路径的差 < 0.5 meV。某个谷报“最高点在边界”或“未确认”时，停下汇报，不要自行扩网格。

## ④ 实验对照 AMSET（重掺杂、300–600 K）

**前提**：用户能拿到四个样品的电导率 σ(T)。拿不到就先不跑。

做法、doping 列表和推送约定见 `14_散射机制分项迁移率/实验对照AMSET/README.md`。要点：

| 材料 | doping（cm⁻³，正为 p） | 作业 A / B 的温度 | 内存 |
|---|---|---|---|
| Sn2Sb2Te5 | 3.0e20, 6.0e20, 1.2e21, 2.4e21 | [300, 400] / [500, 600] | 300G |
| Sn2Bi2Te5 | 3.0e20, 6.0e20, 1.2e21, 2.4e21 | 同上 | 700G |
| Pb2Bi2Te5 | −5.0e18, −1.0e19, −2.0e19, −4.0e19 | 同上 | 700G |
| Pb2Sb2Te5 | −1.0e17, 1.0e17 | 同上 | 300G |

- **settings.yaml**：从各材料的 `step8_amset_full_f1` 复制，只改 `doping`、`temperatures`，其余（`interpolation_factor: 1`、散射机制、`bandgap` 等）一律不动。`wavefunction.h5`、`deformation.h5` 直接沿用。
- **资源**：128 进程，24 h 墙钟，预计每个作业 2–7 h。
- **提交前**：把 submit.sh 给用户过目。
- **推送位置**：`14_散射机制分项迁移率/实验对照AMSET/<材料>/`。不要放进 `密网格AMSET/`。

## ⑤ Sb 系 400–700 K AMSET 收尾（查一下即可）

运行记录中 PbSb、SnSb 的 400–700 K 作业（300783–300786）重提过，但 `14_/密网格AMSET/` 里这两种材料仍然只有 300 K。

- **作业已完成**：把 `transport_*.json` 按 `密网格AMSET/README.md` 的约定入库，doping 列表与 300 K 的相同，可以合并；
- **作业未完成**：报告状态即可，不用重跑。

## 推上来以后 Claude 做什么

- **③**：重算 Sn₂Sb₂Te₅ 的 S（目前范围 134–164 µV/K），更新图 6、图 8、图 9，表 6、表 7 和正文。
- **④**：拿到用户的 σ(T) 后：
  - 用同一个反推浓度预测 σ(T)，与实验比较；
  - 比较两个 Sn 样品的 σ 之比（约 1.3 说明差别来自能带，约 2.5 说明来自浓度）；
  - 计算实验的加权迁移率。
