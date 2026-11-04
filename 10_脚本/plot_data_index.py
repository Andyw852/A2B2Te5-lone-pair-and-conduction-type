"""作图数据的索引与 Excel 汇总（2026-11-04）：由 build_paper_report.py 在打包时调用。

对压缩包里的 data/ 目录：
1. 写 data/作图数据索引.csv：每幅图、每个子图 → 数据文件 → 横轴、纵轴（或颜色、分组）取哪一列 → 说明；
2. 每幅图单独一个 Excel（用户 2026-11-04 要求“每个图数据单独放”），放在 作图数据/<图目录名>.xlsx：
   第一页“说明”是这幅图的索引行（另加工作表名），其后每个 CSV 一页，数值存为数字，可直接在 Excel、Origin 里作图；
   正文表格单独一个 作图数据/tables.xlsx；
3. 核对：索引里写到的文件都存在，data/ 里的每个文件都在索引里，缺一样就报错，保证不漏、不多。

单独运行（检查用）：python 10_脚本/plot_data_index.py <解压后的 data 目录> <输出目录>
"""
import csv, fnmatch, os, re, sys

# (图/表, 子图, 内容, 文件（相对 data/，可用通配符）, 横轴, 纵轴 / 颜色 / 分组, 说明)
INDEX = [
    ('图 1', '(a)(b)', '144 原子超胞的格点占位（O 有序、D 无序）', 'Fig1_structures/*_sites.csv',
     'position_in_plane_row', 'plane_from_bottom_along_111', '颜色按 species；sublattice 为 vacancy 的画空心圆'),
    ('图 1', '(c)(d)', '侧视：各原子面的高度与堆垛位置', 'Fig1_structures/side_view_planes.csv',
     'stacking_position_0_1_2', 'z_A_relative_to_central_Te', '按 structure 分两列；颜色按 species'),
    ('图 1、表 1', '', '三种结构的晶格、体积与能量', 'Fig1_structures/structure_comparison.csv', '', '', '即表 1'),
    ('图 1', '', '结构文件（POSCAR，未弛豫的 O、D 母胞与所给的层状结构）', 'Fig1_structures/*.vasp', '', '', 'VASP 格式，不进 Excel'),
    ('图 2', '(a)–(d)', 'PBE+SOC（虚线）与 HSE06+SOC（实线）路径能带', 'Fig2_bands_PBE_vs_HSE/*2Te5_PBE+SOC.csv',
     'k_distance_1/A（label 列为高对称点）', 'band_N_E-EVBM_eV（每条带一列）', 'SOC 下每条带二重简并，相邻两列相同；能量相对本泛函的价带顶'),
    ('图 2', '(a)–(d)', '同上，HSE06+SOC', 'Fig2_bands_PBE_vs_HSE/*2Te5_HSE06+SOC.csv', 'k_distance_1/A', 'band_N_E-EVBM_eV', ''),
    ('图 2', '(c)(d) 右侧', '附加路径 Γ–Σ′–R₂′–L（穿过镜面内真实价带极大）', 'Fig2_bands_PBE_vs_HSE/*_true_valley_path_G-S-R2-L.csv',
     'k_distance_1/A', 'band_N_E-EVBM_eV', '图中横向放大 1.8 倍'),
    ('图 2、表 2', '', '路径带隙、Γ 处能级与反转带隙', 'Fig2_bands_PBE_vs_HSE/band_gaps_on_path.csv', '', '', '单位 meV'),
    ('图 3', '(a)–(d)', 'HSE06+SOC 能带与阳离子 s 成分', 'Fig3_HSE_bands_cation_s/*_HSE06+SOC_cation_s.csv',
     'k_distance_1/A', 'E-EVBM_eV；颜色、线宽 = cation_s_percent', '长表：按 band 分组连线'),
    ('图 3', '(c)(d) 右侧', '附加路径 Γ–Σ′–R₂′–L', 'Fig3_HSE_bands_cation_s/*_true_valley_path_G-S-R2-L.csv',
     'k_distance_1/A', 'E-EVBM_true_eV；颜色、线宽 = cation_s_percent', '能量相对真实价带顶'),
    ('图 4', '(a)–(d)', '阳离子 s 来自哪个原子（A 位、B 位）', 'Fig4_cation_s_by_atom/*_s_by_atom.csv',
     'k_distance_1/A', 'E-EVBM_eV；线宽 = A_s + B_s，颜色 = A_s/(A_s + B_s)', 'A_s + B_s < 15% 时画灰'),
    ('图 4、表 4', '', '关键态的位置与成分', 'Fig4_cation_s_by_atom/key_states.csv', '', '', '即表 4'),
    ('图 5', '(a)–(d)', 'Sn₂Sb₂Te₅、Sn₂Bi₂Te₅ 的轨道投影', 'Fig5_SnSb_vs_SnBi/*_projections.csv',
     'k_distance_1/A', 'E-EVBM_eV；线宽 = Sn_s_percent（a、b）或 B_p_percent（c、d）', '长表：按 band 分组连线'),
    ('图 5', '(a)–(d) ▼', '路径上的价带谷', 'Fig5_SnSb_vs_SnBi/*_valence_valleys.csv', 'path_k_index', 'E-EVBM_meV', ''),
    ('图 5', '(e)', '价带能级（真实值，相对真实价带顶）', 'Fig5_SnSb_vs_SnBi/levels_true.csv', 'material（两列）', 'E_*_meV', ''),
    ('图 5、表 5', '', '路径上的能级与 Γ 处的态', 'Fig5_SnSb_vs_SnBi/levels.csv', '', '', '相对路径价带顶'),
    ('图 5、表 5', '', 'Γ 处最高价带、最低导带的能量与成分', 'Fig5_SnSb_vs_SnBi/Gamma_states.csv', '', '', ''),
    ('图 6', '(a)–(d)', '价带口袋在布里渊区中的位置', 'Fig6_valence_pockets/pocket_positions.csv',
     'kx_1/A（俯视）；斜视另用 kz_1/A', 'ky_1/A；形状按 valley，颜色按谷的类型', 'plane 列：区内实心、kz = ½ 区界面空心'),
    ('图 6', '(e)', '四种材料的价带能级', 'Fig6_valence_pockets/levels_HSE.csv', 'material', 'E-EVBM_meV；颜色按 kind', ''),
    ('图 6', '(e) 上方', '谷间能量差 ΔE', 'Fig6_valence_pockets/valley_offset_HSE.csv', 'material', 'dE_valley_offset_meV', ''),
    ('图 7', '', '最高的孤对电子价带谷相对价带顶的能量', 'Fig7_lone_pair_valley/highest_lone_pair_valley.csv',
     'material', 'E-EVBM_meV；按 functional 分组', ''),
    ('图 8', '(a)(b)', '300 K 面内 S 随载流子浓度（PBE+SOC + BoltzTraP2）', 'Fig8_Seebeck/S_vs_concentration_300K.csv',
     'concentration_cm-3（对数）', '|S|_PBE+SOC_BoltzTraP2_uV/K', '按 material、carrier_type（p → a，n → b）分组'),
    ('图 8、表 6', '', '表 6：3×10¹⁹、1×10²⁰ cm⁻³ 的 S', 'Fig8_Seebeck/S_300K_table6.csv', '', '', ''),
    ('图 8', '圆点', '实验值：约 301 K 的 S、室温霍尔浓度', 'Fig8_Seebeck/experiment_points.csv',
     '|Hall_conc_signed_cm-3|', 'S_exp_uV/K（n 型取绝对值）', 'Hall_flagged = yes 的画空心'),
    ('图 9', '灰线', '单抛物带 Pisarenko 参考线', 'Fig9_Pisarenko/single_parabolic_band_lines.csv',
     'concentration_cm-3（对数）', '|S|_uV/K（对数）；按 m_DOS_me 分组', ''),
    ('图 9、表 7', '点', '样品：霍尔浓度、实验与计算 S、有效质量', 'Fig9_Pisarenko/samples_Pisarenko_mass.csv',
     '|Hall_conc_signed_cm-3|', '|S_exp_300K|（实心圆）、|S_PBE+SOC_BoltzTraP2_at_Hall_conc|（空心方块）', '即表 7；Pb₂Sb₂Te₅ 的实验 m* 留空（补偿，见 note）'),
    ('图 S1', '(a)(b)', '镜面内最高价带的能量网格', 'FigS1_mirror_plane/*_mirror_plane_grid_meV.csv',
     '行：x=-u（面内）', '列：w（kz）；值 = E − VBM（meV）', '矩阵，画等高图'),
    ('图 S1', '(c)(d)', '真实极大与沿 kz 的切线', 'FigS1_mirror_plane/valleys_and_cuts.csv', 'w', 'E_minus_VBM_meV；按 valley_or_cut 分组', ''),
    ('图 S2', '', 'Γ 点能级（PBE，有、无 SOC）', 'FigS2_Gamma_levels_SOC/gamma_levels_PBE.csv',
     'material、SOC', 'E_minus_Sigma_meV；按 state 分组', ''),
    ('图 S2', '箭头', 'Sb → Bi 变化的标量 / SOC 拆分', 'FigS2_Gamma_levels_SOC/SbBi_decomposition_meV.csv', '', '', ''),
    ('图 S3', '点', '实验 T、σ、S', 'FigS3_S_vs_T/experiment_T_sigma_S.csv', 'T_K', 'S_exp_uV/K（σ 另列，S/m）', '按 material 分四个子图'),
    ('图 S3', '线', '霍尔浓度下计算的 S(T)（不调参数）', 'FigS3_S_vs_T/calculated_S_T_at_Hall_conc.csv',
     'T_K', 'S_PBE+SOC_BoltzTraP2_uV/K', 'PBE+SOC 密网格（28×28×7/6）BoltzTraP2，不加剪刀；T 即计算温度（每 25 K）'),
    ('图 S3', '灰字', '室温霍尔', 'FigS3_S_vs_T/hall_room_temperature.csv', '', '', ''),
    ('图 S3、表 A1', '', '最低、最高温的实验与计算值和偏差', 'FigS3_S_vs_T/summary_at_Hall_conc.csv', '', '', '即表 A1；计算值为 PBE+SOC'),
    ('图 10', '(a)(b)', '12 类缺陷逐位点的中性形成能', 'Fig10_formation_energies/formation_energy_per_site.csv',
     'Ef_Te-rich_eV（a）或 Ef_Te-poor_eV（b）', 'defect_type（括号内 nominal_charge）；颜色按 phase', ''),
    ('图 10', '', '两端化学势', 'Fig10_formation_energies/chemical_potentials.csv', '', '', ''),
    ('图 11', '点', '有序化引起的形成能变化（逐位点）', 'Fig11_ordering/delta_Ef_per_site.csv', 'defect_type', 'dEf_O-D_eV', ''),
    ('图 11', '条', '同上，中位数', 'Fig11_ordering/delta_Ef_median.csv', 'defect_type', 'median_dEf_O-D_eV', ''),
    ('表 9、表 10', 'L 列', '层状相缺陷计算（3×3×1，81 原子）全部顶点的形成能排名', 'L_layered_defects/*.csv', '', '', ''),
    ('表 9、表 10', 'L 列', '同上，原始结果', 'L_layered_defects/*.json', '', '', 'JSON，不进 Excel'),
    ('正文表格', '', '表 1–10、表 A1（从报告 HTML 抽出）', 'tables/*.csv', '', '', ''),
]

# 每幅图的数据来源与绘制方法：键为 data/ 下的目录名。build() 写进每个 Excel 的“说明”页，methods_html() 生成报告中的一节。
METHODS = {
    'Fig1_structures': ('图 1',
        '(a)(b) O、D 两个 144 原子母胞的格点表与未弛豫 POSCAR（13_DVO缺陷_VASP流程与分析 的结构生成结果）；(c) 所给的 9 原子层状 POSCAR；'
        '(d) O 的理想岩盐格点按弛豫后的体积缩放；表 1 的体积与能量来自 PBE-D3 全弛豫。',
        '(a)(b) 把 160 个岩盐格点按沿 [111] 的 10 个原子面分行，每行 16 个格点按超胞平移排序；(c)(d) 取各原子面相对中心 Te 面的高度 z 与堆垛位置（0、1、2）。',
        '(a)(b) 横轴为面内序号、纵轴为原子面；Pb、Sb 为实心圆，Te 为小圆，结构空位为空心圆。(c)(d) 横轴为堆垛位置、纵轴为 z（Å），标出范德华间隙与空位面的 Te–Te 距离。',
        'paper_figures.fig1'),
    'Fig2_bands_PBE_vs_HSE': ('图 2、表 2',
        'c2_pbe、c2_hse：PBE+SOC、HSE06+SOC 在各自 12×12×2 自洽之上用 KPOINTS_OPT 算的 204 点高对称路径；'
        'Sn 系另有 c7_true_path（HSE）、c7_true_path_pbe（PBE）：穿过镜面内真实价带极大的附加路径 Γ–Σ′–R₂′–L。',
        '每种泛函以两段路径上的最高价带为能量零点；横轴为累计 k 距离（Å⁻¹），路径断开处不累加；附加路径横向放大 1.8 倍，接在右侧。',
        '实线 HSE06+SOC、虚线 PBE+SOC，同色；纵轴 −0.66 至 0.66 eV；▼ 标出各泛函价带顶的位置；子图标题为带隙。',
        'paper_figures.fig2'),
    'Fig3_HSE_bands_cation_s': ('图 3、表 3',
        'c2_hse 的 PROCAR_OPT：路径上每个 k 点、每条带的原子轨道投影（HSE06+SOC）；Sn 系附加路径 c7_true_path 的 PROCAR_OPT。',
        '阳离子 s 成分 = （A 位 s + B 位 s）/ 全部投影；能量零点为真实价带顶（Sn 系取镜面内三维搜索得到的价带顶）。',
        '每条带画成分段线，颜色深浅与线宽都随阳离子 s 成分变化；两条虚线为价带顶与导带底；▼ 标出路径上的价带谷。',
        'paper_figures.fig3'),
    'Fig4_cation_s_by_atom': ('图 4、表 4',
        '同图 3（c2_hse、c7_true_path 的 PROCAR_OPT）。',
        '分别求 A 位（Pb、Sn）、B 位（Sb、Bi）、Te 的 s 成分；颜色取 A 位 s /（A 位 s + B 位 s）。',
        '线宽同图 3；颜色从紫（A 位）到橙（B 位），阳离子 s 低于 15% 时渐变为灰。能量零点同图 3。',
        'paper_figures.fig4'),
    'Fig5_SnSb_vs_SnBi': ('图 5、表 5',
        '同图 3，只取 Sn₂Sb₂Te₅、Sn₂Bi₂Te₅；(e) 另用 c5_star 镜面内真实极大的能量。',
        '(a)(b) 取 Sn 5s 投影，(c)(d) 取 B 位 p 投影；(e) 各价带谷相对真实价带顶的能量。',
        '(a)–(d) 线宽表示投影权重，▼ 标出谷；(e) 能级图：橙色为 Γ 单谷，紫色细线为孤对电子谷的真实极大，虚线连接两种材料的同一个态。',
        'paper_figures.fig5'),
    'Fig6_valence_pockets': ('图 6',
        '路径上的价带谷（c2_hse）；Sn 系另用镜面内的真实谷（c4_mirror_v3 网格 + c5_star 星形拟合）。',
        '每个谷用点群 D3d 的对称操作展开成全部等价口袋，换成笛卡尔坐标（Å⁻¹）；能量相对价带顶；ΔE = 最高的孤对电子谷 − Γ–A 轴上的谷。',
        '(a)(b) 沿 c 俯视，(c)(d) 斜视（kz 放大 4 倍），画出第一布里渊区边框；形状区分谷，颜色区分 Γ–A 轴上的谷与孤对电子谷；'
        '区内实心、区界面上空心，离价带顶超过 2kBT 的淡色。(e) 四种材料的价带能级，阴影为 2kBT，ΔE 标在列上方。',
        'paper_figures.fig6'),
    'Fig7_lone_pair_valley': ('图 7',
        '每种材料、每种泛函中最高的孤对电子价带谷（阳离子 s ≥ 15%）：Sn 系 HSE 取 c5_star 镜面搜索，PBE 取 c7_true_path_pbe；Pb 系取路径上的值。',
        '能量相对真实价带顶（meV）。',
        '横轴为材料，两种泛函并排；纵轴为能量；阴影为 2kBT（300 K，52 meV）。',
        'paper_figures.fig7'),
    'Fig8_Seebeck': ('图 8、表 6',
        '22_能带基准_1-0/bt2_pbe_dense：PBE+SOC 28×28×7（Sb 系）、28×28×6（Bi 系）全网格，BoltzTraP2 常数弛豫时间，不加剪刀，插值倍数 60。',
        '取 300 K；面内 S =（S_xx + S_yy）/2；在 1×10¹⁸–3×10²¹ cm⁻³ 取 73 个对数等距的浓度，'
        '每个浓度在 BoltzTraP2 的化学势网格上按 log|N| 线性取值。表 6 取 3×10¹⁹、1×10²⁰ cm⁻³。不用任何模型或比值估计。',
        '双对数：横轴为浓度，纵轴为 S（p 型）或 |S|（n 型）；每种材料一个颜色（分类色板前四格，已做色盲校验）；'
        '圆点为实验（约 301 K 的 S 画在霍尔浓度处，霍尔标红的空心）；点线为 3×10¹⁹ cm⁻³；n 型低浓度 S 变号时曲线向下离开画框。',
        'paper_figures.fig8、pbe_dense'),
    'Fig9_Pisarenko': ('图 9、表 7',
        '实验：最低测量温度（约 301 K）的 S 与室温霍尔浓度（21_阶段报告与交接/实验数据/）。计算：bt2_pbe_dense 在霍尔浓度、300 K 的 S。'
        '参照：HSE06+SOC 能带上单个谷的态密度质量（能带曲率拟合）。',
        '单抛物带、声学声子散射：S =（kB/e）(2F₁/F₀ − η)，n = 4π(2m*kBT/h²)^{3/2}·F_{1/2}(η)。给定 m* 扫 η 得到参考线；'
        '每个样品由 |S| 解出 η，再由霍尔浓度解出 m*，实验点与计算点各自读出。',
        '双对数；灰线为 m* = 0.1、0.3、1、3 me；实心圆为实验，空心方块为计算，同一样品以虚线相连，颜色同图 8；点旁标实验的 m*。'
        'Pb₂Sb₂Te₅（补偿）只画点、不读 m*。',
        'paper_figures.fig9、spb_mass'),
    'FigS1_mirror_plane': ('图 S1',
        'c4_mirror_v3：Sn 系含 c 轴的镜面内 HSE06+SOC 网格（面内步长 1/36，kz 取 12 个值，KPOINTS_OPT 分 6 批）；c5_star：真实极大附近的星形采样与拟合。',
        '最高价带的能量相对真实价带顶（meV）；沿 kz 穿过真实极大取切线。',
        '(a)(b) 色图加等高线，横轴为沿 Γ→M 的面内 k、纵轴为 kz；白色虚线为能带路径的位置；空心标记为路径谷，红色实心为真实极大。'
        '(c)(d) 切线，空心为路径上的值，箭头标路径低估的量。',
        'paper_figures.figS1'),
    'FigS2_Gamma_levels_SOC': ('图 S2',
        '22_能带基准_1-0/SOC机制检验：PBE 下有、无 SOC 的 Γ 点能级（d1：不加 SOC；d2：互换晶格）。',
        '能量相对 Σ 谷；B 位 p 态分成 pz、pxy 与加 SOC 后的 p1/2；Sb → Bi 的变化拆成标量与 SOC 两部分。',
        '能级图：两种材料各两列（无 SOC、有 SOC），紫色为孤对电子 s 态，橙色为 B 位 p 态，箭头上为 SOC 的压低量。',
        'soc_mechanism.py'),
    'FigS3_S_vs_T': ('图 S3、表 A1',
        '实验：四个样品的原始 T–σ–S 表与室温霍尔（21_阶段报告与交接/实验数据/；Pb₂Sb₂Te₅ 用原始表）。计算：bt2_pbe_dense。',
        '浓度固定为霍尔值，不调任何参数；取 BoltzTraP2 的计算温度（300–650 K，每 25 K）；面内 S。'
        '均方根偏差在 8 个测量温度上比较，计算值按 25 K 网格线性内插。',
        '四个子图，空心圆为实验，橙线为计算；灰字为霍尔浓度（标红的注明 flagged），橙字为均方根偏差。',
        'paper_figures.figS3'),
    'Fig10_formation_energies': ('图 10、表 8',
        '13_DVO缺陷_VASP流程与分析/结果_第一批/analysis：D、O 两个 144 原子母胞中 12 类缺陷逐位点的中性形成能（PBE-D3），以及两端化学势。',
        'E_f = E(缺陷) − E(母胞) − Σ nᵢμᵢ；分富 Te、贫 Te 两端；每个位点一个值。',
        '横轴为形成能，纵轴为缺陷类型（上半施主、灰底，下半受主），绿色为 O、紫红为 D；(a) 富 Te，(b) 贫 Te。',
        'paper_figures.fig10'),
    'Fig11_ordering': ('图 11',
        '同图 10。',
        'ΔE_f = E_f(O，该类缺陷的最低位点) − E_f(D，每个位点)，与化学势无关；每类取中位数。',
        '点为 D 的各位点，条为中位数；正值表示退火后更难形成。',
        'paper_figures.fig11'),
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
        if k in ('L_layered_defects', 'tables'):
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
    while name in used:
        name = base[:28] + '_%d' % i; i += 1
    used.add(name)
    return name


def number(x):
    try:
        v = float(x)
        return v if x.strip() not in ('', 'nan') else x
    except ValueError:
        return x


def build(data, out_dir):
    """data：解压后的 data 目录；out_dir：放每幅图 Excel 的目录。写 data/作图数据索引.csv 与 out_dir/<图目录名>.xlsx。
    返回 (文件数, Excel 个数, 工作表数)。"""
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font
    listed, rows = set(), []
    for fig, panel, what, pat, x, y, note in INDEX:
        fs = files_of(data, pat)
        if not fs:
            raise SystemExit('[ERROR] 索引里的 %s 没有对应文件' % pat)
        listed.update(fs); rows.extend((fig, panel, what, rel, x, y, note) for rel in fs)
    every = set(files_of(data, '*')) - {'作图数据索引.csv'}
    if every - listed:
        raise SystemExit('[ERROR] 这些文件不在索引里：%s' % sorted(every - listed))
    with open(os.path.join(data, '作图数据索引.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(HEAD); w.writerows(rows)

    os.makedirs(out_dir, exist_ok=True); nbook = nsheet = 0
    for folder in sorted({r[3].split('/', 1)[0] for r in rows}):
        mine = [r for r in rows if r[3].split('/', 1)[0] == folder]
        wb = Workbook(); ix = wb.active; ix.title = '说明'; used = {'说明'}; sheets = {}
        if folder not in METHODS:
            raise SystemExit('[ERROR] %s 没有写绘制方法（METHODS）' % folder)
        fig_, src_, proc_, draw_, prog_ = METHODS[folder]
        ix.append(['%s：数据来源与绘制方法' % fig_]); ix['A1'].font = Font(bold=True, size=12)
        for k_, v_ in (('数据来源', src_), ('处理', proc_), ('画法', draw_), ('程序（10_脚本/）', prog_)):
            ix.append([k_, v_]); ix.merge_cells(start_row=ix.max_row, start_column=2, end_row=ix.max_row, end_column=8)
            ix.cell(ix.max_row, 1).font = Font(bold=True); ix.row_dimensions[ix.max_row].height = 48
        ix.append([]); ix.append(['各数据表']); ix.cell(ix.max_row, 1).font = Font(bold=True, size=12)
        ix.append(HEAD + ['工作表']); head_row = ix.max_row
        for c in ix[head_row]: c.font = Font(bold=True)
        for fig, panel, what, rel, x, y, note in mine:
            if rel.endswith('.csv') and rel not in sheets:
                sheets[rel] = sheet_name(rel, used)
            ix.append([fig, panel, what, 'data/' + rel, x, y, note, sheets.get(rel, '（不是表格，见 data/ 中的原文件）')])
        for col, wd in zip('ABCDEFGH', (11, 10, 34, 52, 30, 40, 40, 24)):
            ix.column_dimensions[col].width = wd
        for row in ix.iter_rows(min_row=2):
            for c in row: c.alignment = Alignment(wrap_text=True, vertical='top')
        for rel, name in sheets.items():
            ws = wb.create_sheet(name); ws.append(['来源：data/' + rel])
            with open(os.path.join(data, rel), encoding='utf-8-sig') as f:
                for i, r in enumerate(csv.reader(f)):
                    ws.append(r if i == 0 else [number(v) for v in r])
            for c in ws[2]: c.font = Font(bold=True)
            ws.freeze_panes = 'A3'
        wb.save(os.path.join(out_dir, folder + '.xlsx')); nbook += 1; nsheet += len(sheets)
    return len(listed), nbook, nsheet


if __name__ == '__main__':
    print(build(sys.argv[1], sys.argv[2]))
