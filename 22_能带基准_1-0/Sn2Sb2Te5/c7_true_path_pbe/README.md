# Sn2Sb2Te5 · c7_true_path_pbe：真实谷路径的 PBE+SOC（与 c7_true_path 同一组 k 点）

说明、资源与验收见 [`21_阶段报告与交接/交接_PBE真实谷路径_20261101.md`](../../../21_阶段报告与交接/交接_PBE真实谷路径_20261101.md)。

- INCAR = c2_pbe 去掉 NCORE；KPOINTS = c2_pbe；KPOINTS_OPT 与 `../c7_true_path/KPOINTS_OPT` 逐点相同（一次算完，不分批）；POSCAR 与 a_pbe 相同。
- 分析：`python 10_脚本/true_valley_path.py check 22_能带基准_1-0/Sn2Sb2Te5/c7_true_path_pbe --side c2_pbe`
