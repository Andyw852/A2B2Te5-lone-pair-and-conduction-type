# Sn2Bi2Te5 · c7_true_path_pbe：真实谷路径的 PBE+SOC（与 c7_true_path 同一组 k 点）

说明、资源与验收见 [`21_阶段报告与交接/交接_PBE真实谷路径_20261101.md`](../../../21_阶段报告与交接/交接_PBE真实谷路径_20261101.md)。

- INCAR = c2_pbe 去掉 NCORE；KPOINTS = c2_pbe；KPOINTS_OPT 与 `../c7_true_path/KPOINTS_OPT` 逐点相同（一次算完，不分批）；POSCAR 与 a_pbe 相同。
- 分析：`python 10_脚本/true_valley_path.py check 22_能带基准_1-0/Sn2Bi2Te5/c7_true_path_pbe --side c2_pbe`

## 结果（2026-11-01，`check.log`）

DeepSeek 完成（`67b3de9`），自洽能量与 c2_pbe 一致到 1e-7 eV。锚点与 c2_pbe 差 0.048 meV；本路径上价带最高在第 21 点（Σ′ 附近，(−0.188, 0, 0.118)），比 PBE 路径价带顶高 +41.9 meV；Σ′ 处 +39.5、R₂′ 处 +21.3 meV。
比 AutoZT 第 3 步 PBE 密网格上的最高点还高（同一 k 点附近），说明这段路径抓住了 PBE 的真实极大。已接进论文报告图 2(c)(d)、表 2、图 7。
