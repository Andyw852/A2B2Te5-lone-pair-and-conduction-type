# Sn2Sb2Te5 · c7_true_path：穿过真实价带极大的附加能带路径（HSE06+SOC）

由 `10_脚本/true_valley_path.py make Sn2Sb2Te5 --out 22_能带基准_1-0/Sn2Sb2Te5/c7_true_path --split` 生成。说明、资源与检查见 [`21_阶段报告与交接/交接_真实谷能带路径_20261031.md`](../../../21_阶段报告与交接/交接_真实谷能带路径_20261031.md)。

**路径**：Γ → Σ′ → R₂′ → L，全部在含 c 轴的镜面 (u, 0, w) 内。折点 Σ′、R₂′ 是 c5 拟合出的真实极大，坐标见 `true_path.json`。

**k 点**：步长 ≤ 0.015 Å⁻¹，每批开头放一个 c2 路径上的点作锚点；共 2 批（`KPOINTS_OPT_b01`、`_b02`）。

**INCAR、KPOINTS**：与 `c4_mirror_v3` 相同（含 `KPOINTS_OPT_NKBATCH = 36`、`NKRED 2 2 1`）。所有点都已通过 FOCK_FCC 奇数/24 检查。**不要手改 KPOINTS_OPT。**
