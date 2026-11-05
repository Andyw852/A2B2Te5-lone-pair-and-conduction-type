"""作图数据的索引与 Excel 汇总（2026-11-04）：由 build_paper_report.py 在打包时调用。

对压缩包里的 data/ 目录：
1. 写 data/作图数据索引.csv：每幅图、每个子图 → 数据文件 → 横轴、纵轴（或颜色、分组）取哪一列 → 说明；
2. 每幅图单独一个 Excel（用户 2026-11-04 要求“每个图数据单独放”），放在 作图数据/<图目录名>.xlsx：
   第一页“说明”是这幅图的索引行（另加工作表名），其后每个 CSV 一页，数值存为数字，可直接在 Excel、Origin 里作图；
   正文表格单独一个 作图数据/tables.xlsx；
3. 核对：索引里写到的文件都存在，data/ 里的每个文件都在索引里，缺一样就报错，保证不漏、不多。

单独运行（检查用）：python 10_脚本/plot_data_index.py <解压后的 data 目录> <输出目录>
"""
import csv, shutil, fnmatch, os, re, sys

# (图/表, 子图, 内容, 文件（相对 data/，可用通配符）, 横轴, 纵轴 / 颜色 / 分组, 说明)
INDEX = [
    ('图 1', '(a)(b)', '144 原子超胞的格点占位（O 有序、D 无序）', 'Fig1_structures/*_sites.csv',
     'position_in_plane_row', 'plane_from_bottom_along_111', '颜色按 species；sublattice 为 vacancy 的画空心圆'),
    ('图 1', '(c)(d)', '侧视：各原子面的高度与堆垛位置', 'Fig1_structures/side_view_planes.csv',
     'stacking_position_0_1_2', 'z_A_relative_to_central_Te', '按 structure 分两列；颜色按 species'),
    ('图 1、表 1', '', '三种结构的晶格、体积与能量', 'Fig1_structures/structure_comparison.csv', '', '', '即表 1'),
    ('图 1', '', '结构文件（POSCAR，未弛豫的 O、D 母胞与所给的层状结构）', 'Fig1_structures/*.vasp', '', '', 'VASP 格式，不进 Excel'),
    ('图 2-1', '(a)–(d)', 'HSE06+SOC 路径能带（零点为真实价带顶）', 'Fig2-1_HSE_bands/*_HSE06+SOC_path.csv', 'k_distance_1/A', 'E-EVBM_eV；按 band 分线', '价带蓝、导带红'),
    ('图 2-1', '(c)(d) 右侧小图', '附加路径 Γ–Σ′–R₂′–L（不在高对称路径上）', 'Fig2-1_HSE_bands/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A（图中 ×1.8）', 'E-EVBM_eV', '只有 Sn 系'),
    ('图 2-1', '标题', '带隙（真实值与路径值）', 'Fig2-1_HSE_bands/band_gaps_HSE06+SOC.csv', '', '', ''),
    ('图 2-2', '(a)–(d)', 'PBE+SOC（虚线）路径能带', 'Fig2-2_bands_PBE_vs_HSE/*2Te5_PBE+SOC.csv',
     'k_distance_1/A（label 列为高对称点）', 'band_N_E-EVBM_eV（每条带一列）', 'SOC 下每条带二重简并，相邻两列相同；能量相对本泛函的价带顶'),
    ('图 2-2', '(a)–(d)', 'HSE06+SOC（实线）路径能带', 'Fig2-2_bands_PBE_vs_HSE/*2Te5_HSE06+SOC.csv', 'k_distance_1/A', 'band_N_E-EVBM_eV', ''),
    ('图 2-2', '(c)(d) 右侧小图', '附加路径 Γ–Σ′–R₂′–L', 'Fig2-2_bands_PBE_vs_HSE/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A', 'band_N_E-EVBM_eV', '图中横向放大 1.8 倍'),
    ('图 2-2、表 2', '', '路径带隙、Γ 处能级与反转带隙', 'Fig2-2_bands_PBE_vs_HSE/band_gaps_on_path.csv', '', '', '单位 meV'),
    ('图 3-1', '(a)–(d)', 'HSE06+SOC 能带与阳离子 s 成分', 'Fig3-1_HSE_bands_cation_s/*_HSE06+SOC_cation_s.csv',
     'k_distance_1/A', 'E-EVBM_eV；颜色、线宽 = cation_s_percent', '长表：按 band 分组连线'),
    ('图 3-1', '(c)(d) 右侧小图', '附加路径 Γ–Σ′–R₂′–L', 'Fig3-1_HSE_bands_cation_s/*_true_valley_path_G-S-R2-L.csv',
     'k_distance_1/A', 'E-EVBM_true_eV；颜色、线宽 = cation_s_percent', '能量相对真实价带顶'),
    ('图 3-2', '(a)–(h)', '阳离子–Te 键的轨道分辨 COHP（−COHP，正为成键）', 'Fig3-2_COHP/*_COHP.csv',
     '-COHP_…（每个轨道对一列）', 'E-EF_eV', '上排 A 位–Te、下排 B 位–Te；着色为阳离子 s–Te p'),
    ('图 3-2', '框内数值', '积分到 E_F 的 ICOHP（eV/键）', 'Fig3-2_COHP/ICOHP_summary.csv', '', '', '正值为反键'),
    ('图 4-1', '(a)–(d)', '阳离子 s 来自哪个原子（A 位、B 位）', 'Fig4-1_cation_s_by_atom/*_s_by_atom.csv',
     'k_distance_1/A', 'E-EVBM_eV；线宽 = A_s + B_s，颜色 = A_s/(A_s + B_s)', 'A_s + B_s < 15% 时画灰'),
    ('图 4-1、表 4', '', '关键态的位置与成分', 'Fig4-1_cation_s_by_atom/key_states.csv', '', '', '即表 4'),
    ('图 4-2、图 S3a–d', '全部子图', '各轨道成分：A 位 s、B 位 s、A 位 p、B 位 p、Te p（另有 Te s 与阳离子 s 合计）', 'Fig4-2_HSE_orbitals/*_HSE06+SOC_path.csv',
     'k_distance_1/A', 'E-EVBM_eV；圆点面积（图 4-2）或颜色、线宽（图 S3）按各成分（%）', ''),
    ('图 4-2、图 S3c–d', 'Sn 系右侧小图', '同上，附加路径', 'Fig4-2_HSE_orbitals/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A（图中 ×1.8）', '同上', ''),
    ('图 5', '(a)–(d)', 'Sn₂Sb₂Te₅、Sn₂Bi₂Te₅ 的轨道投影', 'Fig5_SnSb_vs_SnBi/*_projections.csv',
     'k_distance_1/A', 'E-EVBM_eV；线宽 = Sn_s_percent（a、b）或 B_p_percent（c、d）', '长表：按 band 分组连线'),
    ('图 5', '(a)–(d) ▼', '路径上的价带谷', 'Fig5_SnSb_vs_SnBi/*_valence_valleys.csv', 'path_k_index', 'E-EVBM_meV', ''),
    ('图 5', '(e)', '价带能级（真实值，相对真实价带顶）', 'Fig5_SnSb_vs_SnBi/levels_true.csv', 'material（两列）', 'E_*_meV', ''),
    ('图 5、表 5', '', '路径上的能级与 Γ 处的态', 'Fig5_SnSb_vs_SnBi/levels.csv', '', '', '相对路径价带顶'),
    ('图 5、表 5', '', 'Γ 处最高价带、最低导带的能量与成分', 'Fig5_SnSb_vs_SnBi/Gamma_states.csv', '', '', ''),
    ('图 6-1', '', '四种材料价带上的谷：能量、等价口袋数、k、阳离子 s、主轴有效质量', 'Fig6-1_valence_levels/valence_valleys_HSE06+SOC.csv',
     'material', 'E-EVBM_meV；颜色按 kind', ''),
    ('图 6-2', '(a)–(d)', 'Sn 系价带口袋在布里渊区中的位置', 'Fig6-2_pockets_BZ/pocket_positions.csv', 'kx_1/A（俯视）；三维另用 kz_1/A', 'ky_1/A；形状按 valley', ''),
    ('图 7', '', '最高的孤对电子价带谷相对价带顶的能量', 'Fig7_lone_pair_valley/highest_lone_pair_valley.csv',
     'material', 'E-EVBM_meV；按 functional 分组', ''),
    ('图 8-1', '(a)(b)', '300 K 面内 S 随载流子浓度（HSE06+SOC + BoltzTraP2）', 'Fig8-1_Seebeck_HSE/S_vs_concentration_300K_HSE.csv',
     'net_concentration_cm-3（对数）', 'S…（取绝对值作图；为负的一段即 S 已变号，画点线）', '按 material、carrier_type（p → a，n → b）分组'),
    ('图 8-2', '(a)–(d)', '同上，HSE 与 PBE', 'Fig8-2_Seebeck_HSE_vs_PBE/S_vs_concentration_300K_HSE_and_PBE.csv',
     'net_concentration_cm-3（对数）', '|S|；按 functional（颜色）、carrier_type（线型）分组', '每种材料一个子图'),
    ('表 6', '', '300 K、3×10¹⁹ 与 1×10²⁰ cm⁻³ 的 S（HSE 与 PBE）', 'Table6_S_300K/S_300K_HSE_and_PBE.csv', '', '', '即表 6'),
    ('图 9-1、9-2', '点', '实验 T、σ、S', 'Fig9-*/experiment_T_sigma_S.csv', 'T_K', 'S_exp_uV/K（σ 另列，S/m）', '按 material 分四个子图'),
    ('图 9-1、9-2', '灰字', '室温霍尔', 'Fig9-*/hall_room_temperature.csv', '', '', ''),
    ('图 9-1', '线', '霍尔浓度下 HSE06+SOC 计算的 S(T)（不调参数）', 'Fig9-1_S_vs_T_HSE/calculated_S_T_at_Hall_conc.csv', 'T_K', 'S_BoltzTraP2_uV/K', 'T 即计算温度（每 25 K）'),
    ('图 9-2', '线', '同上，HSE 与 PBE', 'Fig9-2_S_vs_T_HSE_vs_PBE/calculated_S_T_at_Hall_conc.csv', 'T_K', 'S_BoltzTraP2_uV/K；按 functional 分线', ''),
    ('表 7', '', '霍尔浓度下的 S：实验、HSE、PBE，最低与最高测量温度、均方根偏差', 'Table7_S_T_at_Hall/S_T_at_Hall_HSE_and_PBE.csv', '', '', '即表 7'),
    ('表 8', '', 'Pisarenko 有效质量：实验、HSE、PBE 在霍尔浓度下的 S 与读出的 m*（不作图）', 'Table8_Pisarenko_mass/samples_Pisarenko_mass_HSE_and_PBE.csv',
     '', '', '即表 8；Pb₂Sb₂Te₅ 的实验 m* 留空（补偿）'),
    ('图 10-1', '虚线', 'HSE06+SOC 能带上 300 K 的化学势 µ（霍尔浓度、3×10¹⁹ cm⁻³）', 'Fig10-1_Fermi_window_HSE/chemical_potential_300K_HSE.csv',
     '', 'mu-EVBM_meV', '能带与阳离子 s 成分同图 3-1（Fig3-1_HSE_bands_cation_s/）'),
    ('图 10-2', '虚线', 'HSE 与 PBE 各自在霍尔浓度下的 µ', 'Fig10-2_Fermi_window_HSE_vs_PBE/chemical_potential_300K_Hall_HSE_and_PBE.csv',
     '', 'mu-EVBM_meV（各自的价带顶为零点）', '能带同图 2-2（Fig2-2_bands_PBE_vs_HSE/）'),
    ('图 S1', '(a)(b)', '镜面内最高价带的能量网格', 'FigS1_mirror_plane/*_mirror_plane_grid_meV.csv', '行：x=-u（面内）', '列：w（kz）；值 = E − VBM（meV）', '矩阵，画等高图'),
    ('图 S1', '标记', '路径上的谷与镜面内的真实极大', 'FigS1_mirror_plane/valleys_in_mirror_plane.csv', 'x=-u', 'w', ''),
    ('图 S1、第 2 节', '', '路径上的谷与真实极大的能量差', 'FigS1_mirror_plane/path_vs_true_maxima.csv', '', 'E_path_meV、E_true_meV、rise_meV', '不作图'),
    ('图 S2', '', 'Γ 点能级（PBE，有、无 SOC）', 'FigS2_Gamma_levels_SOC/gamma_levels_PBE.csv', 'material、SOC', 'E_minus_Sigma_meV；按 state 分组', ''),
    ('图 S2', '箭头', 'Sb → Bi 变化的标量 / SOC 拆分', 'FigS2_Gamma_levels_SOC/SbBi_decomposition_meV.csv', '', '', ''),
    ('图 11', '(a)(b)', '12 类缺陷逐位点的中性形成能', 'Fig11_formation_energies/formation_energy_per_site.csv',
     'Ef_Te-rich_eV（a）或 Ef_Te-poor_eV（b）', 'defect_type（括号内 nominal_charge）；颜色按 phase', ''),
    ('图 11', '', '两端化学势', 'Fig11_formation_energies/chemical_potentials.csv', '', '', ''),
    ('图 12', '点', '有序化引起的形成能变化（逐位点）', 'Fig12_ordering/delta_Ef_per_site.csv', 'defect_type', 'dEf_O-D_eV', ''),
    ('图 12', '条', '同上，中位数', 'Fig12_ordering/delta_Ef_median.csv', 'defect_type', 'median_dEf_O-D_eV', ''),
    ('表 9、表 10', 'L 列', '层状相缺陷计算（3×3×1，81 原子）全部顶点的形成能排名', 'L_layered_defects/*.csv', '', '', ''),
    ('表 9、表 10', 'L 列', '同上，原始结果', 'L_layered_defects/*.json', '', '', 'JSON，不进 Excel'),
    ('正文表格', '', '表 1–11（从报告 HTML 抽出）', 'tables/*.csv', '', '', ''),
]

# 每幅图的数据来源与绘制方法：键为 data/ 下的目录名。build() 写进每个 Excel 的“说明”页，methods_html() 生成报告中的一节。
METHODS = {
    'Fig1_structures': ('图 1',
        '(a)(b) O、D 两个 144 原子母胞的格点表与未弛豫 POSCAR（13_DVO缺陷_VASP流程与分析 的结构生成结果）；(c) 所给的 9 原子层状 POSCAR；'
        '(d) O 的理想岩盐格点按弛豫后的体积缩放；表 1 的体积与能量来自 PBE-D3 全弛豫。',
        '(a)(b) 把 160 个岩盐格点按沿 [111] 的 10 个原子面分行，每行 16 个格点按超胞平移排序；(c)(d) 取各原子面相对中心 Te 面的高度 z 与堆垛位置（0、1、2）。',
        '(a)(b) 横轴为面内序号、纵轴为原子面；Pb、Sb 为实心圆，Te 为小圆，结构空位为空心圆。(c)(d) 横轴为堆垛位置、纵轴为 z（Å），标出范德华间隙与空位面的 Te–Te 距离。',
        'paper_figures.fig1'),
    'Fig2-1_HSE_bands': ('图 2-1',
        '22_能带基准_1-0/<材料>/c2_hse（PROCAR_OPT、path_eigen）；Sn 系附加路径 <材料>/c7_true_path。', '零点为真实价带顶（Sn₂Bi₂Te₅ 取镜面搜索的 Σ′，其余为路径最高点）；SOC 下每条带是二重简并的一对，源数据逐条列出。',
        '2×2，每个子图有自己的纵坐标；价带蓝、导带红；两条虚线为价带顶与导带底；▼ 标出价带谷；Sn 系的附加路径单独画在右侧小图。', 'hse_band_figures.h1'),
    'Fig2-2_bands_PBE_vs_HSE': ('图 2-2、表 2',
        'c2_pbe、c2_hse：PBE+SOC、HSE06+SOC 在各自 12×12×2 自洽之上用 KPOINTS_OPT 算的 204 点高对称路径；'
        'Sn 系另有 c7_true_path（HSE）、c7_true_path_pbe（PBE）：穿过镜面内真实价带极大的附加路径 Γ–Σ′–R₂′–L。',
        '每种泛函以两段路径上的最高价带为能量零点；横轴为累计 k 距离（Å⁻¹），路径断开处不累加；附加路径横向放大 1.8 倍，单独画在右侧小图（不在高对称路径上）。',
        '实线 HSE06+SOC、虚线 PBE+SOC，同色；纵轴 −0.66 至 0.66 eV；▼ 标出各泛函价带顶的位置；子图标题为带隙。',
        'paper_figures.fig2'),
    'Fig3-1_HSE_bands_cation_s': ('图 3-1、表 3',
        'c2_hse 的 PROCAR_OPT：路径上每个 k 点、每条带的原子轨道投影（HSE06+SOC）；Sn 系附加路径 c7_true_path 的 PROCAR_OPT。',
        '阳离子 s 成分 = （A 位 s + B 位 s）/ 全部投影；能量零点为真实价带顶（Sn 系取镜面内三维搜索得到的价带顶）。',
        '每条带画成分段线，颜色深浅与线宽都随阳离子 s 成分变化；两条虚线为价带顶与导带底；▼ 标出路径上的价带谷。',
        'paper_figures.fig3'),
    'Fig3-2_COHP': ('图 3-2',
        '06_COGITO成键分析与COHP：COGITO（PBE、不加 SOC、所给结构）的轨道分辨 COHP；由交互图 <材料>_COHP交互图.html 取出（10_脚本/cohp_extract.py）。',
        '交互图里每根键一组堆叠曲线，相邻相减得各轨道对；总 COHP 积分到 E_F 与 bond_info.txt 的 ICOHP 恰差 1.50 倍（含两根等长键的组为 3.0 倍），换算成每根键；取 −COHP（正为成键）。ICOHP 正值为反键。',
        '2×4：上排 A 位–Te、下排 B 位–Te，每种材料取最短的键；横轴 −COHP、纵轴 E − E_F；阳离子 s–Te p 着色（反键一侧填色），p–p 灰，总 COHP 黑细线；框内为 ICOHP。', 'paper_figures.fig3_2、cohp_extract.py'),
    'Fig4-1_cation_s_by_atom': ('图 4-1、表 4',
        '同图 3-1（c2_hse、c7_true_path 的 PROCAR_OPT）。',
        '分别求 A 位（Pb、Sn）、B 位（Sb、Bi）、Te 的 s 成分；颜色取 A 位 s /（A 位 s + B 位 s）。',
        '线宽同图 3-1；颜色从紫（A 位）到橙（B 位），阳离子 s 低于 15% 时渐变为灰。能量零点同图 3-1。',
        'paper_figures.fig4'),
    'Fig4-2_HSE_orbitals': ('图 4-2、图 S3a–d',
        '同图 3-1（c2_hse、c7_true_path 的 PROCAR_OPT）。', '零点为真实价带顶（Sn₂Bi₂Te₅ 取镜面搜索的 Σ′，其余为路径最高点）；SOC 下每条带是二重简并的一对，源数据逐条列出。 p 成分为 p_x、p_y、p_z 投影之和。',
        '图 4-2：每种材料一行，左为 A 位 s 与 B 位 s 叠在一张能带图上，右为 A 位 p、B 位 p、Te p 叠在一起；圆点面积 ∝ 轨道比例，颜色按原子所在的位（A 紫、B 橙、Te 绿），同一 k 点大的点先画；▼ 标出价带谷。'
        '图 S3a–d：同一组数据，每种材料一张 3×2，五个成分各一个子图，按位着色的渐变，满刻度 50%。', 'hse_band_figures.h4、h4_each'),
    'Fig5_SnSb_vs_SnBi': ('图 5、表 5',
        '同图 3，只取 Sn₂Sb₂Te₅、Sn₂Bi₂Te₅；(e) 另用 c5_star 镜面内真实极大的能量。',
        '(a)(b) 取 Sn 5s 投影，(c)(d) 取 B 位 p 投影；(e) 各价带谷相对真实价带顶的能量。',
        '(a)–(d) 线宽表示投影权重，▼ 标出谷；(e) 能级图：橙色为 Γ 单谷，紫色细线为孤对电子谷的真实极大，虚线连接两种材料的同一个态。',
        'paper_figures.fig5'),
    'Fig6-1_valence_levels': ('图 6-1、表 3',
        '22_能带基准_1-0/<材料>/c2_hse；Sn 系 c4_mirror_v3、c5_star、c6_local（镜面内的真实极大）。',
        'Γ–A 轴上的谷取路径上 Γ–A 段的最高价带；孤对电子谷：Sn 系取镜面搜索的真实极大，Pb 系取 Σ、R 线上阳离子 s ≥ 15% 的局部极大（路径值）。等价口袋数由 D3d 点群与时间反演展开。',
        '每种材料一列，线为能级，右侧标名称、能量与 ×n；浅蓝底为 2k_BT（300 K）。', 'hse_band_figures.h5'),
    'Fig6-2_pockets_BZ': ('图 6-2',
        '同图 6-1（Sn 系）。', '谷展开成第一布里渊区内的全部等价点，换到笛卡尔坐标。',
        '左右两列对齐：(a)(b) 沿 c 俯视，(c)(d) 三维透视（k_z 放大 4 倍），下方为各组谷的能量与等价口袋数；形状区分谷，离价带顶超过 2k_BT 的画浅色。', 'hse_band_figures.h6'),
    'Fig7_lone_pair_valley': ('图 7',
        '每种材料、每种泛函中最高的孤对电子价带谷（阳离子 s ≥ 15%）：Sn 系 HSE 取 c5_star 镜面搜索，PBE 取 c7_true_path_pbe；Pb 系取路径上的值。',
        '能量相对真实价带顶（meV）。',
        '横轴为材料，两种泛函并排；纵轴为能量；阴影为 2kBT（300 K，52 meV）。',
        'paper_figures.fig7'),
    'Fig8-1_Seebeck_HSE': ('图 8-1',
        '22_能带基准_1-0/<材料>/c8_hse_dense/gene/bt2_mult60：HSE06+SOC 28×28×7（Sb 系，515 个不可约点）、28×28×6（Bi 系，450 个）KPOINTS_OPT 非自洽能带（分 13–15 批，锚点一致到 0.1 meV），BoltzTraP2 常数弛豫时间、不加剪刀、插值倍数 60（与 120 比 S 最多差 0.3 µV/K）。',
        '取 300 K；面内 S =（S_xx + S_yy）/2；1×10¹⁸–3×10²¹ cm⁻³ 取 73 个对数等距的浓度，在化学势网格上按 log|N| 线性取值。不用任何模型或比值估计。',
        '与图 8-2 同一组数据；双对数，横轴浓度、纵轴 S（p 型）或 |S|（n 型），每种材料一个颜色和符号（同实验第一版 Origin 图）；点线为 3×10¹⁹ cm⁻³；S 变号的一段画点线。', 'transport_figures.fig8_1'),
    'Fig8-2_Seebeck_HSE_vs_PBE': ('图 8-2',
        'HSE：同图 8-1；PBE：22_能带基准_1-0/bt2_pbe_dense（PBE+SOC 28×28×7/6 全网格，BoltzTraP2，不加剪刀，插值倍数 60）。',
        '同图 8-1。', '每种材料一个子图；蓝 HSE、橙 PBE；实线 p 型、虚线 n 型（|S|）；变号的一段为细点线；灰竖线为样品霍尔浓度。', 'transport_figures.fig8_2'),
    'Table6_S_300K': ('表 6', '同图 8-2。', '300 K，±3×10¹⁹、±1×10²⁰ cm⁻³ 的 S。', '（不作图）', 'transport_figures.tables_hse_pbe'),
    'Fig9-1_S_vs_T_HSE': ('图 9-1',
        '实验：四个样品的原始 T–σ–S 表与室温霍尔（21_阶段报告与交接/实验数据/；Pb₂Sb₂Te₅ 用原始表）。计算：22_能带基准_1-0/<材料>/c8_hse_dense/gene/bt2_mult60：HSE06+SOC 28×28×7（Sb 系，515 个不可约点）、28×28×6（Bi 系，450 个）KPOINTS_OPT 非自洽能带（分 13–15 批，锚点一致到 0.1 meV），BoltzTraP2 常数弛豫时间、不加剪刀、插值倍数 60（与 120 比 S 最多差 0.3 µV/K）。',
        '浓度固定为霍尔值，不调任何参数；取 BoltzTraP2 的计算温度（300–650 K，每 25 K）；面内 S。均方根偏差在 8 个测量温度上比较，计算值按 25 K 网格线性内插。',
        '四个子图（每个样品一个），空心圆为实验，蓝线为 HSE；灰字为霍尔浓度（标红的注明 flagged），蓝字为均方根偏差。', 'transport_figures.fig9_1'),
    'Fig9-2_S_vs_T_HSE_vs_PBE': ('图 9-2',
        '同图 9-1，另加 PBE（bt2_pbe_dense）。', '同图 9-1。', '空心圆为实验，蓝实线 HSE、橙虚线 PBE，两种泛函的均方根偏差标在左上（n 型在左下）。', 'transport_figures.fig9_2'),
    'Table7_S_T_at_Hall': ('表 7', '同图 9-2。', '最低、最高测量温度的实验与计算值；300–620 K 的均方根偏差。', '（不作图）', 'transport_figures.tables_hse_pbe'),
    'Table8_Pisarenko_mass': ('表 8',
        '实验：最低测量温度（约 301 K）的 S 与室温霍尔浓度。计算：HSE（c8_hse_dense）与 PBE（bt2_pbe_dense）在霍尔浓度、300 K 的 S。参照：HSE06+SOC 能带上单个谷的态密度质量（能带曲率拟合）。',
        '单抛物带、声学声子散射：S =（kB/e）(2F₁/F₀ − η)，n = 4π(2m*kBT/h²)^{3/2}·F_{1/2}(η)。由 |S| 解出 η，再由霍尔浓度解出 m*，实验与两种计算各自读出。Pb₂Sb₂Te₅（电子、空穴补偿）的实验 m* 不适用，留空。',
        '（不作图）', 'transport_figures.tables_hse_pbe、paper_figures.spb_mass'),
    'Fig10-1_Fermi_window_HSE': ('图 10-1',
        '能带：c2_hse 路径与 PROCAR 的阳离子 s 成分（同图 3-1），Sn 系附加路径 c7_true_path；化学势：c8_hse_dense 的 BoltzTraP2（300 K）；对齐用 GENE 网格（c8_hse_dense/gene 的 .energy 与 source.json）。',
        'µ 取化学势网格上带符号浓度等于霍尔值（或 ±3×10¹⁹ cm⁻³）处的值，按 log|N| 线性取值。BoltzTraP2 的 µ 相对 GENE 的费米能（网格带隙中点）；GENE 网格与路径都含 Γ，用 Γ 处最高价带与最低导带把网格移到路径的零点（两者给出的偏移相差 < 0.3 meV）。零点为真实价带顶。',
        '四个子图（每种材料一个，各有纵坐标），能带颜色与线宽表示阳离子 s 成分；橙色虚线为霍尔浓度下的 µ，灰色虚线为同类型 3×10¹⁹ cm⁻³ 下的 µ；Sn 系附加路径单独画在右侧小图。', 'transport_figures.fig10_1'),
    'Fig10-2_Fermi_window_HSE_vs_PBE': ('图 10-2',
        'HSE 同图 10-1；PBE：c2_pbe 路径（c7_true_path_pbe）与 bt2_pbe_dense，对齐用 A2B2Te5_SOC_dense 的 step3_uniform（EIGENVAL、OUTCAR 的 E-fermi）。',
        '两种泛函各以本泛函的真实价带顶为零点；µ 的取法与对齐同图 10-1。', '蓝实线 HSE 能带、橙虚线 PBE 能带；同色长虚线为各自在霍尔浓度下 300 K 的 µ。', 'transport_figures.fig10_2'),
    'FigS1_mirror_plane': ('图 S1',
        'c4_mirror_v3：Sn 系含 c 轴的镜面内 HSE06+SOC 网格（面内步长 1/36，kz 取 12 个值，KPOINTS_OPT 分 6 批）；真实极大取 c5_star、c6_local。',
        '最高价带（二重简并的一对取平均），相对真实价带顶；用 (u, w) ≡ (−u, −w) 补全另一半，三次样条插值后画图。path_vs_true_maxima.csv 为路径上的谷与真实极大的能量差，不作图。',
        '色图加等高线（每 25 meV），低于 −250 meV 的按 −250 meV 画；白色虚线为能带路径经过的位置；空心为路径上的谷，红色实心为真实极大。', 'hse_band_figures.h7'),
    'FigS2_Gamma_levels_SOC': ('图 S2',
        '22_能带基准_1-0/SOC机制检验：PBE 下有、无 SOC 的 Γ 点能级（d1：不加 SOC；d2：互换晶格）。',
        '能量相对 Σ 谷（Γ–M 段内最高价带的局部极大；2026-11-04 改正，原取段上最大值，不加 SOC 时取错）；B 位 p 态分成 pz、pxy 与加 SOC 后的 p1/2；Sb → Bi 的变化拆成标量与 SOC 两部分。',
        '能级图：两种材料各两列（无 SOC、有 SOC），紫色为孤对电子 s 态，橙色为 B 位 p 态，箭头上为 SOC 的压低量。',
        'soc_mechanism.py'),
    'Fig11_formation_energies': ('图 11、表 9',
        '13_DVO缺陷_VASP流程与分析/结果_第一批/analysis：D、O 两个 144 原子母胞中 12 类缺陷逐位点的中性形成能（PBE-D3），以及两端化学势。',
        'E_f = E(缺陷) − E(母胞) − Σ nᵢμᵢ；分富 Te、贫 Te 两端；每个位点一个值。',
        '横轴为形成能，纵轴为缺陷类型（上半施主、灰底，下半受主），绿色为 O、紫红为 D；(a) 富 Te，(b) 贫 Te。',
        'paper_figures.fig11'),
    'Fig12_ordering': ('图 12',
        '同图 11。',
        'ΔE_f = E_f(O，该类缺陷的最低位点) − E_f(D，每个位点)，与化学势无关；每类取中位数。',
        '点为 D 的各位点，条为中位数；正值表示退火后更难形成。',
        'paper_figures.fig12'),
    'L_layered_defects': ('表 9、表 10 的 L 列',
        '早先一批在所给层状结构上（3×3×1，81 原子）的缺陷计算，原始结果取自 03_缺陷形成能与PN型数据/Pb2Sb2Te5.zip。',
        '各化学势顶点下的形成能排名与原始结果，原样保存。', '（不作图）', 'build_paper_report.package'),
    'tables': ('正文表格',
        '报告 HTML 中带 data-csv 属性的表。', '原样抽出（下标写成 _x，上标写成 ^x）。', '（不作图）', 'build_paper_report.Tables'),
}


def methods_html():
    """报告“附：每幅图的数据与绘制方法”一节的表格（只列作图的目录）。"""
    import html
    out = ['<table>', '<thead><tr><th>图</th><th>数据来源</th><th>处理</th><th>画法</th><th>程序</th></tr></thead>', '<tbody>']
    for k, (fig, src, proc, draw, prog) in METHODS.items():
        if k in ('L_layered_defects', 'tables') or k.startswith('Table'):
            continue
        out.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><code>%s</code></td></tr>'
                   % tuple(html.escape(x) for x in (fig, src, proc, draw, prog)))
    return '\n'.join(out + ['</tbody></table>'])


HEAD = ['图/表', '子图', '内容', '数据文件（data/ 下）', '横轴', '纵轴 / 颜色 / 分组', '说明']
SHORT = {'Pb2Sb2Te5': 'PbSb', 'Pb2Bi2Te5': 'PbBi', 'Sn2Sb2Te5': 'SnSb', 'Sn2Bi2Te5': 'SnBi'}


def files_of(data, pat):
    out = []
    for dp, _, fs in os.walk(data):
        for f in fs:
            rel = os.path.relpath(os.path.join(dp, f), data).replace(os.sep, '/')
            if fnmatch.fnmatch(rel, pat):
                out.append(rel)
    return sorted(out, key=lambda t: [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', t)])   # Table2 排在 Table10 前


def sheet_name(rel, used):
    stem = os.path.splitext(rel.split('/', 1)[1])[0]
    for k, v in SHORT.items():
        stem = stem.replace(k, v)
    stem = stem.replace('_HSE06+SOC', '_HSE').replace('_PBE+SOC', '_PBE').replace('true_valley_path_G-S-R2-L', 'truepath')
    name = re.sub(r'[\[\]:*?/\\]', '_', stem)[:31]
    base, i = name, 2
    while name.lower() in {u.lower() for u in used}:     # Excel 的工作表名不分大小写
        name = base[:28] + '_%d' % i; i += 1
    used.add(name)
    return name


def number(x):
    try:
        v = float(x)
        return v if x.strip().lower() not in ('', 'nan', 'inf', '-inf', 'infinity', '-infinity') else x
    except ValueError:
        return x


CATS = {   # 主报告压缩包的数据分类（2026-11-04 用户要求分门别类）：data/<类别>/<图目录>/，作图数据/<类别>/<图目录>.xlsx
    'Fig1_structures': '01_结构',
    'Fig2-1_HSE_bands': '02_能带与轨道成分', 'Fig2-2_bands_PBE_vs_HSE': '02_能带与轨道成分', 'Fig3-1_HSE_bands_cation_s': '02_能带与轨道成分',
    'Fig3-2_COHP': '02_能带与轨道成分', 'Fig4-1_cation_s_by_atom': '02_能带与轨道成分', 'Fig4-2_HSE_orbitals': '02_能带与轨道成分',
    'Fig5_SnSb_vs_SnBi': '02_能带与轨道成分', 'FigS2_Gamma_levels_SOC': '02_能带与轨道成分',
    'Fig6-1_valence_levels': '03_价带谷与口袋', 'Fig6-2_pockets_BZ': '03_价带谷与口袋', 'Fig7_lone_pair_valley': '03_价带谷与口袋', 'FigS1_mirror_plane': '03_价带谷与口袋',
    'Fig8-1_Seebeck_HSE': '04_输运与实验对照', 'Fig8-2_Seebeck_HSE_vs_PBE': '04_输运与实验对照', 'Table6_S_300K': '04_输运与实验对照',
    'Fig9-1_S_vs_T_HSE': '04_输运与实验对照', 'Fig9-2_S_vs_T_HSE_vs_PBE': '04_输运与实验对照', 'Table7_S_T_at_Hall': '04_输运与实验对照',
    'Table8_Pisarenko_mass': '04_输运与实验对照', 'Fig10-1_Fermi_window_HSE': '04_输运与实验对照', 'Fig10-2_Fermi_window_HSE_vs_PBE': '04_输运与实验对照',
    'Fig11_formation_energies': '05_缺陷与导电类型', 'Fig12_ordering': '05_缺陷与导电类型', 'L_layered_defects': '05_缺陷与导电类型',
    'tables': '06_正文表格',
}


def build(data, out_dir, cats=None):
    """data：解压后的 data 目录；out_dir：放每幅图 Excel 的目录。写 data/作图数据索引.csv 与 out_dir/<图目录名>.xlsx。
    返回 (文件数, Excel 个数, 工作表数)。"""
    listed, rows = set(), []
    for fig, panel, what, pat, x, y, note in INDEX:
        fs = files_of(data, pat)
        if not fs:
            raise SystemExit('[ERROR] 索引里的 %s 没有对应文件' % pat)
        listed.update(fs); rows.extend((fig, panel, what, rel, x, y, note) for rel in fs)
    every = set(files_of(data, '*')) - {'作图数据索引.csv'}
    if every - listed:
        raise SystemExit('[ERROR] 这些文件不在索引里：%s' % sorted(every - listed))
    cats = cats or {}
    cat = lambda folder: cats.get(folder, '')
    if cats:                                                   # 每个图目录移到 data/<类别>/ 下，索引与 Excel 里的路径跟着改
        miss = sorted({r[3].split('/', 1)[0] for r in rows} - set(cats))
        if miss: raise SystemExit('[ERROR] 这些目录没有分类：%s' % miss)
        for folder in sorted({r[3].split('/', 1)[0] for r in rows}):
            os.makedirs(os.path.join(data, cat(folder)), exist_ok=True)
            shutil.move(os.path.join(data, folder), os.path.join(data, cat(folder), folder))
    full = lambda rel: (cat(rel.split('/', 1)[0]) + '/' if cat(rel.split('/', 1)[0]) else '') + rel
    with open(os.path.join(data, '作图数据索引.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(HEAD); w.writerows([r[:3] + (full(r[3]),) + r[4:] for r in rows])

    os.makedirs(out_dir, exist_ok=True); nbook = nsheet = 0
    for folder in sorted({r[3].split('/', 1)[0] for r in rows}):
        mine = [r for r in rows if r[3].split('/', 1)[0] == folder]
        if folder not in METHODS:
            raise SystemExit('[ERROR] %s 没有写绘制方法（METHODS）' % folder)
        os.makedirs(os.path.join(out_dir, cat(folder)), exist_ok=True)
        path = os.path.join(out_dir, cat(folder), folder + '.xlsx')
        nsheet += write_book(path, folder, mine, data, full); nbook += 1
        check_book(path)
    return len(listed), nbook, nsheet


def write_book(path, folder, mine, data, full):
    """用 xlsxwriter 写一个 Excel（openpyxl 写的文件在用户的 Excel 里报损坏，2026-11-05 改）。
    不合并单元格；数值用 write_number（NaN、inf 写成空格）、文字用 write_string（不会被当成公式）；工作表名大小写不敏感地去重。"""
    import xlsxwriter
    wb = xlsxwriter.Workbook(path, {'strings_to_formulas': False, 'strings_to_urls': False, 'nan_inf_to_errors': False})
    bold = wb.add_format({'bold': True}); title = wb.add_format({'bold': True, 'font_size': 12})
    wrap = wb.add_format({'text_wrap': True, 'valign': 'top'}); wrapb = wb.add_format({'text_wrap': True, 'valign': 'top', 'bold': True})
    fig_, src_, proc_, draw_, prog_ = METHODS[folder]
    ix = wb.add_worksheet('说明'); used = {'说明'}; sheets = {}
    ix.write_string(0, 0, '%s：数据来源与绘制方法' % fig_, title)
    r = 1
    for k_, v_ in (('数据来源', src_), ('处理', proc_), ('画法', draw_), ('程序（10_脚本/）', prog_)):
        ix.write_string(r, 0, k_, wrapb); ix.write_string(r, 1, v_, wrap); r += 1
    r += 1; ix.write_string(r, 0, '各数据表', title); r += 1
    for j, h in enumerate(HEAD + ['工作表']): ix.write_string(r, j, h, bold)
    r += 1
    for fig, panel, what, rel, x, y, note in mine:
        if rel.endswith('.csv') and rel not in sheets:
            sheets[rel] = sheet_name(rel, used)
        for j, v in enumerate([fig, panel, what, 'data/' + full(rel), x, y, note, sheets.get(rel, '（不是表格，见 data/ 中的原文件）')]):
            ix.write_string(r, j, v, wrap)
        r += 1
    ix.set_column(0, 0, 14); ix.set_column(1, 1, 90); ix.set_column(2, 2, 34); ix.set_column(3, 3, 52)
    ix.set_column(4, 6, 34); ix.set_column(7, 7, 24)
    for rel, name in sheets.items():
        ws = wb.add_worksheet(name); ws.write_string(0, 0, '来源：data/' + full(rel))
        with open(os.path.join(data, full(rel)), encoding='utf-8-sig') as f:
            for i, row in enumerate(csv.reader(f)):
                for j, v in enumerate(row):
                    if i == 0: ws.write_string(1, j, v, bold); continue
                    x_ = number(v)
                    if isinstance(x_, float):
                        if x_ != x_ or x_ in (float('inf'), float('-inf')): ws.write_string(i + 1, j, v)
                        else: ws.write_number(i + 1, j, x_)
                    else: ws.write_string(i + 1, j, x_)
        ws.freeze_panes(2, 0)
    wb.close()
    return len(sheets)


def check_book(path):
    """回读检查：openpyxl 能打开、每页都有内容、没有公式。"""
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True)
    for ws in wb.worksheets:
        n = 0
        for row in ws.iter_rows(values_only=True):
            n += 1
            for v in row:
                if isinstance(v, str) and v.startswith('='): raise SystemExit('[ERROR] %s 的 %s 有公式样的内容' % (path, ws.title))
        if n == 0: raise SystemExit('[ERROR] %s 的 %s 是空的' % (path, ws.title))
    wb.close()


if __name__ == '__main__':
    print(build(sys.argv[1], sys.argv[2]))
