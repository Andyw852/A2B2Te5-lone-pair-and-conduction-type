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
    ('图 8', '(a)(b)', '300 K 面内 S 随载流子浓度', 'Fig8_Seebeck/S_vs_concentration_300K.csv',
     'concentration_cm-3（对数）', '|S|_uV/K', '按 material、carrier_type（p → a，n → b）、bands 分组；PBE+SOC 四种材料取最小、最大作灰带。HSE 为暂用估计'),
    ('图 8、表 6', '', '3×10¹⁹ cm⁻³ 的 S', 'Fig8_Seebeck/S_300K_3e19.csv', '', '', 'HSE 为暂用估计'),
    ('表 6 附', '', '逐谷电导份额', 'Fig8_Seebeck/valley_conductivity_shares.csv', '', '', ''),
    ('表 6 附', '', '孤对电子谷电导份额汇总', 'Fig8_Seebeck/lone_pair_share_summary.csv', '', '', ''),
    ('表 6 附', '', '孤对电子谷份额随空穴浓度', 'Fig8_Seebeck/lone_pair_share_vs_concentration.csv',
     'hole_conc_cm-3', 'lone_pair_valley_share_percent', ''),
    ('1.4 节', '', 'Sn₂Bi₂Te₅ 价带谷 c3 拟合检查', 'Fig8_Seebeck/Sn2Bi2Te5_c3_valley_check.csv', '', '', ''),
    ('图 9', '灰线', '单抛物带 Pisarenko 参考线', 'Fig9_Pisarenko/single_parabolic_band_lines.csv',
     'concentration_cm-3（对数）', '|S|_uV/K（对数）；按 m_DOS_me 分组', ''),
    ('图 9、表 7', '点', '样品：霍尔浓度、实验与计算 S、有效质量', 'Fig9_Pisarenko/samples_Pisarenko_mass.csv',
     '|Hall_conc_signed_cm-3|', '|S_exp_300K|、|S_HSE_at_Hall_conc|、|S_PBE_at_Hall_conc|', '即表 7；HSE 为暂用估计'),
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
        wb = Workbook(); ix = wb.active; ix.title = '说明'; used = {'说明'}
        ix.append(HEAD + ['工作表']); sheets = {}
        for c in ix[1]: c.font = Font(bold=True)
        for fig, panel, what, rel, x, y, note in mine:
            if rel.endswith('.csv') and rel not in sheets:
                sheets[rel] = sheet_name(rel, used)
            ix.append([fig, panel, what, 'data/' + rel, x, y, note, sheets.get(rel, '（不是表格，见 data/ 中的原文件）')])
        for col, wd in zip('ABCDEFGH', (11, 10, 34, 52, 30, 40, 40, 24)):
            ix.column_dimensions[col].width = wd
        for row in ix.iter_rows(min_row=2):
            for c in row: c.alignment = Alignment(wrap_text=True, vertical='top')
        ix.freeze_panes = 'A2'
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
