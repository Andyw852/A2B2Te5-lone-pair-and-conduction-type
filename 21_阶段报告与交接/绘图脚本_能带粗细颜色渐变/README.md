# 能带“粗细渐变 + 颜色渐变”绘图包（2026-11-06）

报告图 10-1（以及图 3-1）的画法：能带的线宽和颜色同时随阳离子 s 成分（孤对电子 s 态的权重）变化，权重越大线越粗、颜色越深；水平虚线为 300 K 的化学势。整个文件夹可以单独拷走使用。

## 文件

| 文件 | 内容 |
|---|---|
| `plot_fatband_gradient.py` | 绘图脚本，只依赖 numpy 与 matplotlib |
| `data/HSE06+SOC/` | 报告图 10-1 的数据（HSE06+SOC） |
| `data/PBE+SOC/` | 同一画法的 PBE+SOC 数据（截图里那张旧版图 10 就是它：Sn₂Sb₂Te₅ E_g = 146 meV、Sn₂Bi₂Te₅ 75 meV） |
| `fatband_HSE06+SOC.png/.pdf` | 脚本默认输出，与报告图 10-1 逐像素一致 |
| `fatband_PBE+SOC.png/.pdf`、`fatband_PBE+SOC_Sn_only.png/.pdf` | PBE+SOC 四种材料、只画两种 Sn 材料（同截图） |

每个 `data/<泛函>/` 中：

| 文件 | 列 | 用途 |
|---|---|---|
| `<材料>_main_path.csv` | `k_index, k_distance_1/A, label, band, E-EVBM_eV, cation_s_percent` | 高对称路径 Γ–M–K–Γ–A–L–H–A 上价带顶以下 14 条、以上 6 条能带（SOC 下 Kramers 对各一行）。`label` 非空的行是高对称点；`band` 中 `_VB_top`、`_CB_bottom` 标出价带顶、导带底 |
| `<材料>_off_path_G-S-R2-L.csv` | 同上（`i` 代替 `k_index`） | 只有 Sn 系：穿过镜面内真实价带极大 Σ′、R₂′ 的附加路径，不在高对称路径上，图中单独画在右侧浅紫底小图，横向放大 1.8 倍（k 距离存的是原始值） |
| `chemical_potential_300K.csv` | `material, case, carrier_conc_signed_cm-3, mu-EVBM_eV, Hall_flagged_by_instrument` | 300 K 的化学势：`Hall` 为样品霍尔浓度（橙色虚线），`3e19` 为同类型 3×10¹⁹ cm⁻³（灰色虚线）。浓度正为空穴、负为电子 |
| `band_gap.csv` | `material, Eg_meV (incl. off path), …` | 子图标题里的带隙 |
| `作图数据_<泛函>.xlsx` | 同上，每个 CSV 一页 | 用 Origin 等软件作图时用 |

能量零点都是本泛函的真实价带顶（Sn₂Bi₂Te₅ 的价带顶在附加路径的 Σ′ 上）。数据来源：HSE06+SOC 与 PBE+SOC 路径能带及其 PROCAR 投影（`22_能带基准_1-0/<材料>/c2_hse`、`c2_pbe`、`c7_true_path*`），化学势取自同一泛函的 BoltzTraP2 结果；由仓库中 `10_脚本/export_fatband_kit.py` 导出，与报告图 10-1 用的是同一个取数函数。

## 用法

```
python plot_fatband_gradient.py                                 # HSE06+SOC 四种材料（报告图 10-1）
python plot_fatband_gradient.py --functional PBE+SOC            # PBE+SOC
python plot_fatband_gradient.py --mats Sn2Sb2Te5 Sn2Bi2Te5      # 只画其中几种
python plot_fatband_gradient.py --no-mu                         # 不画化学势（即图 3-1 的样式）
python plot_fatband_gradient.py --refine 4                      # 每段细分 4 份，渐变更平滑
python plot_fatband_gradient.py --weight-sat 25 --lw-min 0.4 --lw-add 3.0
python plot_fatband_gradient.py --colors "#FDE0DD,#FA9FB5,#C51B8A,#7A0177"   # 换渐变色（浅 → 深）
python plot_fatband_gradient.py --out my_figure                  # 输出 my_figure.png 与 my_figure.pdf
```

## 渐变是怎么做的

1. 每条能带拆成相邻 k 点之间的小线段，交给 `matplotlib.collections.LineCollection`。
2. 每段的权重 w 取两端阳离子 s 成分的平均。同一个 w 同时决定：
   - **颜色**：`LineCollection(..., cmap=渐变色, norm=Normalize(0, 满刻度))`，再 `set_array(w)`；
   - **线宽**：`linewidths = lw_min + lw_add × min(w / 满刻度, 1)`。
   默认满刻度 30%、`lw_min = 0.5`、`lw_add = 2.4` pt，与报告相同；超过满刻度的按满刻度画。
3. `capstyle='round'` 让相邻两段线宽不同时接口处不留缺口。
4. `--refine N` 把每段再细分 N 份，能量与权重线性内插，粗细与颜色的过渡更平滑（数据点约 0.025 Å⁻¹ 一个，默认不细分已足够）。

核心只有 `gradient_line()` 一个函数。换成自己的数据时，准备同样列名的长表 CSV（每行一个 k 点、一条带），`cation_s_percent` 换成任何 0–100 的权重即可。
