# 11_续算审计_DeepSeek — 只读诊断与后处理

本目录是对 `A2B2Te5_SOC_dense` 四材料 SOC 密网格、BoltzTraP2 CRTA、有效质量拟合、
旧 HSE s 投影和 `补算实施方案_20261007` 中 D/V/O 几何包的**只读审计**。没有改写
`A2B2Te5_SOC_dense` 下任何原始文件，没有提交/恢复任何远端作业，没有运行 VASP、
BoltzTraP2 或 AMSET。

## 入口

- 主报告：`审计报告_20261017.md`
- 复现：`bash scripts/run_all.sh`
- 来源与哈希：`json/source_manifest.json`、`json/hash_check_summary.json`、`csv/hash_check.csv`
- 审计产物清单：`json/audit_output_manifest.json`

## 2026-10-18 仓库发布说明

本次保留原审计报告和 CSV/JSON 数值，未重新运行全套审计。原报告中的绝对路径及
`json/source_manifest.json`、`json/audit_output_manifest.json` 记录的是原运行环境和当时哈希，
不是当前克隆目录的路径。发布后的文件哈希见仓库根目录 `发布清单_20261018.json`。

配套方案已放在仓库根目录 `补算实施方案_20261007/`；`scripts/audit_dvo.py` 现在优先从
该目录读取，兼容原仓库外相邻目录。只适配输入路径，未改变能量、拟合或输运算法。

完整复现仍需要自行补齐 `.gitignore` 列出的 8 个超大 `vasprun.xml`，以及原审计使用的
外部旧 81 原子项目 `/mnt/d/tf_data`。未上传这些原始数据，不能将单次 Git 克隆描述为
完整计算环境。已有审计 CSV/JSON 可以直接阅读；`run_all.sh` 会写本目录的后处理产物，
在另外的副本中复现可保留本次发布的历史结果。

原报告的诊断分级不等于最终物理结论。尤其是带隙中点不必等于电荷中性点，
中隙处净浓度非零本身不能证明浓度零点错误；应结合 DOS 和电子/空穴平衡复核该解释。

## 目录

```
11_续算审计_DeepSeek/
├── 审计报告_20261017.md
├── README_审计说明.md
├── scripts/
│   ├── audit_common.py            独立 EIGENVAL/POSCAR/KPOINTS 解析与二次拟合
│   ├── audit_inventory.py         四材料设置、输出、provenance、哈希清单
│   ├── audit_mass.py              125 点局部网格的独立质量复算与窗口/边界诊断
│   ├── audit_bt2.py               BoltzTraP2 JSON 审计 + PCHIP 固定浓度提取
│   ├── audit_s_projection.py      旧 HSE s 投影/50 meV 窗口与 PBE 密网格对照
│   ├── audit_dvo.py               D/V/O 几何包与旧 81/72 原子缺陷数据核查
│   ├── make_manifest.py           来源哈希与审计产物哈希
│   └── run_all.sh                 依次运行以上脚本
├── csv/
│   ├── inventory.csv
│   ├── mass_fit_recheck.csv
│   ├── mass_fit_variants.csv
│   ├── mass_variant_summary.csv
│   ├── mass_point_diagnostics.csv
│   ├── band_edge_check.csv
│   ├── band_tracking.csv
│   ├── bt2_checks.csv
│   ├── bt2_mu_grid_summary.csv
│   ├── bt2_interpolated_S_sigma.csv
│   ├── bt2_interpolation_sensitivity.csv
│   ├── bt2_reliable_S_sigma.csv
│   ├── bt2_marginal_S_sigma.csv
│   ├── bt2_unreliable_S_sigma.csv
│   ├── s_projection_audit.csv
│   ├── s_projection_source_check.csv
│   ├── dvo_status.csv
│   ├── dvo_old_defects.csv
│   ├── dvo_old_81atom.csv
│   ├── dvo_old_72atom_zips.csv
│   └── hash_check.csv
├── json/
│   ├── inventory_summary.json
│   ├── mass_audit_summary.json
│   ├── bt2_audit_summary.json
│   ├── s_projection_summary.json
│   ├── dvo_audit_summary.json
│   ├── hash_check_summary.json
│   ├── source_manifest.json
│   └── audit_output_manifest.json
└── logs/
    ├── audit_inventory.log
    ├── audit_mass.log
    ├── audit_s_projection.log
    ├── audit_bt2.log
    ├── audit_dvo.log
    ├── make_manifest.log
    ├── step3_vasprun_sha256.txt
    ├── test_effective_mass.log
    └── run_all_final.log
```

## 口径

- 质量单位：`m_e`，`ħ²/m_e = 7.619964 eV·Å²`；`q` 为 Cartesian `Å⁻¹`。
- 质量顺序：电子能量 Hessian 本征值升序；空穴质量按同一顺序取 `ħ²/(−λ)`，因此空穴
  质量列表不是升序，和旧 `mass_tensor.json` 一致。
- 浓度：`carrier_conc_cm-3` 负 = n 型电子，正 = p 型空穴；由 `(N + nelect)/V` 得到。
- `mu_rel_efermi_eV = 0` 是 scissor 后重新对中的带隙中央，不是原始 PBE 费米能。
- `sigma/tau`、`kappa_e/tau` 是 per-τ 量；本审计不报告绝对电导率、热导率或迁移率。
- 标量张量值是 3D 对角平均 `(xx+yy+zz)/3`；层状材料的面内值看 `*_xx/yy`。
- 固定浓度提取只用 PCHIP/线性插值，所有目标都在 `c(μ)` 范围内；不外推。
- 可靠性分层：A = 单载流子近似下可用；B = 边缘/接近交叉/插值敏感；C = 双极或符号不可靠。

## 未执行

- 没有运行 VASP/DFT、BoltzTraP2、AMSET 或任何 HPC/远端任务。
- 没有弛豫 D/V/O 母胞，没有计算任何 D/V/O 或 81 原子缺陷能量。
- 没有把旧 81 原子或 72 原子缺陷数据当作 D/V/O；没有恢复旧缺陷作业。
- 没有读取或输出任何凭据。
