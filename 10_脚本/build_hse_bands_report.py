#!/usr/bin/env python3
"""HSE06+SOC 能带专题的单文件 HTML 与压缩包（2026-11-04）。

用法：先 python 10_脚本/hse_band_figures.py（图与源数据），再 python 10_脚本/build_hse_bands_report.py
输入：21_阶段报告与交接/HSE能带_20261104/HSE能带_20261104.src.html、…/图/（含 source_data/）
输出：21_阶段报告与交接/HSE能带专题_20261104.html（图片嵌入）、21_阶段报告与交接/A2B2Te5_HSE_bands_20261104.zip
占位符：{{STYLE}} 主报告的样式；{{FIG:名称}} 图；{{TABLE:名称}} 由 source_data/tables/名称.csv 生成的表；{{PLOTMETHODS}} 绘制方法表。
"""
import base64, csv, glob, html, io, os, re, shutil, sys, tempfile, zipfile
from html.parser import HTMLParser
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import plot_data_index as pdi
D = os.path.join(ROOT, '21_阶段报告与交接'); HD = os.path.join(D, 'HSE能带_20261104')
FIGS = os.path.join(HD, '图'); DATA = os.path.join(FIGS, 'source_data')
NAME = 'A2B2Te5_HSE_bands_20261104'

INDEX = [
    ('图 H1', '(a)–(d)', 'HSE06+SOC 路径能带（零点为真实价带顶）', 'H1_bands/*_HSE06+SOC_path.csv', 'k_distance_1/A', 'E-EVBM_eV；按 band 分线', '价带蓝、导带红'),
    ('图 H1', '(c)(d) 右段', '附加路径 Γ–Σ′–R₂′–L', 'H1_bands/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A（图中 ×1.8）', 'E-EVBM_eV', '只有 Sn 系'),
    ('图 H1、表 1', '标题', '带隙（真实值与路径值）', 'H1_bands/band_gaps_HSE06+SOC.csv', '', '', ''),
    ('图 H2', '(a)–(d)', '阳离子 s 成分', 'H2_cation_s/*_HSE06+SOC_path.csv', 'k_distance_1/A', 'E-EVBM_eV；颜色与线宽按 cation_s_percent', ''),
    ('图 H2', '(c)(d) 右段', '同上，附加路径', 'H2_cation_s/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A（图中 ×1.8）', '同上', ''),
    ('图 H3', '(a)–(d)', '阳离子 s 来自 A 位还是 B 位', 'H3_s_by_site/*_HSE06+SOC_path.csv', 'k_distance_1/A',
     'E-EVBM_eV；颜色按 A_s/(A_s+B_s)，线宽按 cation_s_percent', ''),
    ('图 H3', '(c)(d) 右段', '同上，附加路径', 'H3_s_by_site/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A（图中 ×1.8）', '同上', ''),
    ('图 H4、H4a–d', 'H4 左列 s、右列 p；H4a–d 每个成分一个子图', '每种材料的 A 位 s、B 位 s（左），A 位 p、B 位 p、Te p（右）', 'H4_orbitals/*_HSE06+SOC_path.csv', 'k_distance_1/A',
     'E-EVBM_eV；圆点面积按 A_s、B_s 或 A_p、B_p、Te_p（%）', '另列 cation_s、Te_s'),
    ('图 H4', 'Sn 系右侧小图', '同上，附加路径', 'H4_orbitals/*_true_valley_path_G-S-R2-L.csv', 'k_distance_1/A（图中 ×1.8）', '同上', ''),
    ('图 H5、表 3', '', '价带上的谷：能量、等价口袋数、k、阳离子 s、主轴有效质量', 'H5_valence_levels/valence_valleys_HSE06+SOC.csv',
     'material', 'E-EVBM_meV；颜色按 kind', ''),
    ('图 H6', '(a)(b)', '口袋中心在第一布里渊区中的位置', 'H6_pockets/pocket_positions.csv', 'kx_1/A', 'ky_1/A；形状按 valley', ''),
    ('图 H7', '(a)(b)', '镜面内最高价带的能量网格', 'H7_mirror_plane/*_mirror_plane_grid_meV.csv', '行：x = −u', '列：w；值 = E − VBM（meV）', '矩阵，画等高图'),
    ('图 H7', '标记', '路径上的谷与镜面内的真实极大', 'H7_mirror_plane/valleys_in_mirror_plane.csv', 'x=-u', 'w', ''),
    ('图 H7、第 5 节', '', '路径上的谷与真实极大的能量差', 'H7_mirror_plane/path_vs_true_maxima.csv', '', 'E_path_meV、E_true_meV、rise_meV', '不作图'),
    ('图 H8', '', 'Γ 点能级（有、无 SOC；文件名注明泛函）', 'H8_Gamma_levels_SOC/gamma_levels_*.csv', 'material、SOC', 'E_minus_Sigma_meV；按 state', ''),
    ('图 H8', '箭头', 'Sb → Bi 变化的标量 / SOC 拆分', 'H8_Gamma_levels_SOC/SbBi_decomposition_meV.csv', '', '', ''),
    ('表 4', '', '加 SOC 时 Γ 点能级，PBE 与 HSE', 'H8_Gamma_levels_SOC/Gamma_SOC_levels_PBE_vs_HSE.csv', '', '', ''),
    ('正文表格', '', '表 1–4', 'tables/*.csv', '', '', ''),
]
P = '22_能带基准_1-0/<材料>/c2_hse（PROCAR_OPT、path_eigen）'
BAND = '零点为真实价带顶（Sn₂Bi₂Te₅ 取镜面搜索的 Σ′，其余为路径最高点）；SOC 下每条带是二重简并的一对，源数据逐条列出。'
METHODS = {
    'H1_bands': ('图 H1、表 1', P + '；Sn 系附加路径 <材料>/c7_true_path。', BAND, '2×2，每个子图有自己的纵坐标；价带蓝、导带红；两条虚线为价带顶与导带底；Sn 系的附加路径单独画在右侧小图。', 'hse_band_figures.h1'),
    'H2_cation_s': ('图 H2', P + '；c7_true_path。', BAND + ' 阳离子 s = 离子 1–2（B 位）与 8–9（A 位）的 s 投影之和 ÷ 球内总投影。',
                    '颜色深浅与线宽都按阳离子 s；四个子图共用一个普通色柱，满刻度 30%，超过的按 30% 画。', 'hse_band_figures.h2'),
    'H3_s_by_site': ('图 H3', P + '；c7_true_path。', BAND + ' A 位 s、B 位 s 分别取离子 8–9、1–2 的 s 投影。',
                     '线宽按阳离子 s 总量；颜色按 A 位所占比例（紫 A、橙 B），阳离子 s < 15% 时变灰。', 'hse_band_figures.h3'),
    'H4_orbitals': ('图 H4、H4a–d', P + '；c7_true_path。', BAND + ' p 成分为 p_x、p_y、p_z 投影之和。',
                    '每种材料一行：左为 s（A 位 s、B 位 s 叠在同一张能带图上），右为 p（A 位 p、B 位 p、Te p 叠在一起）；圆点面积 ∝ 轨道比例，颜色按原子所在的位（A 紫、B 橙、Te 绿），同一 k 点上大的点先画；▼ 标出价带上的谷。H4a–d 为同一组数据：每种材料一张 3×2，A 位 s、B 位 s、A 位 p、B 位 p、Te p 各一个子图，按位着色的渐变（A 紫、B 橙、Te 绿），满刻度 50%。', 'hse_band_figures.h4、h4_each'),
    'H5_valence_levels': ('图 H5、表 3', P + '；Sn 系 c4_mirror_v3、c5_star、c6_local。',
                          'Γ–A 轴上的谷取路径上 Γ–A 段的最高价带；孤对电子谷：Sn 系取镜面搜索的真实极大，Pb 系取 Σ、R 线上阳离子 s ≥ 15% 的局部极大（路径值）。等价口袋数由 D3d 点群与时间反演展开。',
                          '每种材料一列，线为能级，右侧标名称、能量与 ×n；浅蓝底为 2k_BT（300 K）。', 'hse_band_figures.h5、valley_list'),
    'H6_pockets': ('图 H6', '同图 H5（Sn 系）。', '谷展开成第一布里渊区内的全部等价点，换到笛卡尔坐标，沿 c 投影。',
                   '(a)(b) 六边形为第一布里渊区沿 c 的投影，右侧列表为能量与 ×n；(c)(d) 三维透视，k_z 放大 4 倍；形状区分谷，离价带顶超过 2k_BT 的画浅色。', 'hse_band_figures.h6'),
    'H7_mirror_plane': ('图 H7', '22_能带基准_1-0/<材料>/c4_mirror_v3（6 批）；真实极大取 c5_star、c6_local。',
                        '最高价带（二重简并的一对取平均），相对真实价带顶；用 (u, w) ≡ (−u, −w) 补全另一半，三次样条插值后画图。',
                        '色图加等高线（每 25 meV），低于 −250 meV 的按 −250 meV 画；白色虚线为能带路径经过的位置；空心为路径上的谷，红色实心为真实极大。path_vs_true_maxima.csv 为两者的能量差（路径值取 c2_hse，真实值取 c5/c6），不作图。', 'hse_band_figures.h7、paper_figures.mirror_map'),
    'H8_Gamma_levels_SOC': ('图 H8、表 4', 'PBE：22_能带基准_1-0/SOC机制检验/d1_nosoc_<材料>（不加 SOC）与 <材料>/c2_pbe（加 SOC）；表 4 另取 c2_hse。',
                            'Γ 点的态按轨道成分追踪；能量相对 Σ 谷 = Γ–M 段内最高价带的局部极大（2026-11-04 改正，原取段上最大值）。',
                            '能级图：两种材料各两列（不加、加 SOC）；箭头为 SOC 把 p 态压低的量。', 'soc_mechanism.figure、hse_band_figures.h8'),
    'tables': ('正文表格', '本页 HTML 的表 1，以及 hse_band_figures.tables 生成的表 2–4。', '原样。', '（不作图）', 'hse_band_figures.tables'),
}

CATS = {'H1_bands': '01_能带与带隙', 'H2_cation_s': '02_轨道成分', 'H3_s_by_site': '02_轨道成分', 'H4_orbitals': '02_轨道成分',
        'H5_valence_levels': '03_价带谷与口袋', 'H6_pockets': '03_价带谷与口袋', 'H7_mirror_plane': '04_镜面搜索',
        'H8_Gamma_levels_SOC': '05_Γ点能级与SOC', 'tables': '06_正文表格'}


class Tables(HTMLParser):
    """抽出带 data-csv 属性的手写表格（表 1）；下标写成 _x，上标写成 ^x。"""
    def __init__(self):
        super().__init__(); self.tables = []; self.cur = None; self.cell = None

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'table' and 'data-csv' in a: self.cur = {'name': a['data-csv'], 'rows': []}
        elif self.cur is not None and tag == 'tr': self.cur['rows'].append([])
        elif self.cur is not None and tag in ('td', 'th'): self.cell = ''
        elif tag in ('sub', 'sup') and self.cell is not None: self.cell += '_' if tag == 'sub' else '^'

    def handle_endtag(self, tag):
        if tag == 'table' and self.cur is not None: self.tables.append(self.cur); self.cur = None
        elif self.cur is not None and tag in ('td', 'th') and self.cell is not None:
            self.cur['rows'][-1].append(' '.join(self.cell.split())); self.cell = None

    def handle_data(self, data):
        if self.cell is not None: self.cell += data


def csv_table(m):
    name = m.group(1)
    rows = list(csv.reader(open(os.path.join(DATA, 'tables', name + '.csv'), encoding='utf-8-sig')))
    def num(v):
        try: float(v); return True
        except ValueError: return False
    out = [f'<table data-csv-generated="{name}">', '<thead><tr>' + ''.join(f'<th>{html.escape(h)}</th>' for h in rows[0]) + '</tr></thead>', '<tbody>']
    for r in rows[1:]:
        out.append('<tr>' + ''.join(f'<td class="num">{v.replace("-", "−")}</td>' if num(v) else f'<td>{html.escape(v)}</td>' for v in r) + '</tr>')
    return '\n'.join(out + ['</tbody></table>'])


def embed(m):
    im = Image.open(os.path.join(FIGS, m.group(1) + '.png')).convert('RGB')
    if im.width > 1800:
        im = im.resize((1800, round(im.height * 1800 / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, 'PNG', optimize=True)
    return 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode()


README = """A₂B₂Te₅ 的 HSE06+SOC 能带（专题），2026-11-04
=====================================================
HSE能带专题.html          报告（图片用相对路径，需与 figures/ 同目录）
HSE能带专题_单文件.html    同上，图片已嵌入，可单独转发
figures/png/  figures/pdf/  图 H1–H8、H4a–d（600 dpi PNG 与矢量 PDF），文件名即图号
data/<类别>/<图目录>/      每幅图的源数据（CSV，UTF-8），目录与图的对应见 data/作图数据索引.csv。类别：
                           01_能带与带隙（H1）；02_轨道成分（H2、H3、H4、H4a–d）；03_价带谷与口袋（H5、H6）；
                           04_镜面搜索（H7）；05_Γ点能级与SOC（H8、表 4）；06_正文表格（表 1–4）
作图数据/<类别>/<图目录>.xlsx  同样的数据，每幅图一个 Excel；第一页“说明”为数据来源、处理与画法

全部为 HSE06+SOC，只有图 H8（SOC 机制）目前是 PBE：HSE 不加 SOC 的计算见仓库 22_能带基准_1-0/SOC机制检验_HSE/。
能量零点为真实价带顶；Sn 系的孤对电子谷取镜面内三维搜索的真实极大。
程序：10_脚本/hse_band_figures.py、10_脚本/build_hse_bands_report.py、10_脚本/soc_mechanism.py
"""


def main():
    main_src = open(os.path.join(D, '论文报告_20261028.src.html'), encoding='utf-8').read()
    style = main_src[main_src.index('<style>'):main_src.index('</style>') + len('</style>')]
    pdi.INDEX, pdi.METHODS = INDEX, METHODS
    src = open(os.path.join(HD, 'HSE能带_20261104.src.html'), encoding='utf-8').read()
    src = src.replace('{{STYLE}}', style).replace('{{PLOTMETHODS}}', pdi.methods_html())
    src = re.sub(r'\{\{TABLE:([^}]+)\}\}', csv_table, src)
    p = Tables(); p.feed(src)
    for t in p.tables:                                          # 手写的表（表 1）也存成 CSV
        with open(os.path.join(DATA, 'tables', t['name'] + '.csv'), 'w', newline='', encoding='utf-8-sig') as f:
            csv.writer(f).writerows(t['rows'])
    standalone = re.sub(r'\{\{FIG:([^}]+)\}\}', embed, src)
    out = os.path.join(D, 'HSE能带专题_20261104.html'); open(out, 'w', encoding='utf-8').write(standalone)
    print(out, f'{len(standalone) / 1e6:.1f} MB')
    tmp = tempfile.mkdtemp(); top = os.path.join(tmp, NAME)
    for ext in ('png', 'pdf'):
        os.makedirs(os.path.join(top, 'figures', ext))
        for f in sorted(glob.glob(os.path.join(FIGS, f'H*.{ext}'))): shutil.copy(f, os.path.join(top, 'figures', ext))
    shutil.copytree(DATA, os.path.join(top, 'data'))
    nfile, nbook, nsheet = pdi.build(os.path.join(top, 'data'), os.path.join(top, '作图数据'), CATS)
    print(f'作图数据：{nfile} 个文件全部在索引里；{nbook} 个 Excel、{nsheet} 个工作表')
    open(os.path.join(top, 'HSE能带专题.html'), 'w', encoding='utf-8').write(re.sub(r'\{\{FIG:([^}]+)\}\}', lambda m: f'figures/png/{m.group(1)}.png', src))
    open(os.path.join(top, 'HSE能带专题_单文件.html'), 'w', encoding='utf-8').write(standalone)
    open(os.path.join(top, 'README.txt'), 'w', encoding='utf-8-sig', newline='').write(README.replace('\n', '\r\n'))
    z = os.path.join(D, NAME + '.zip')
    with zipfile.ZipFile(z, 'w', zipfile.ZIP_DEFLATED) as zf:
        for dp, _, fs in sorted(os.walk(top)):
            for f in sorted(fs):
                full = os.path.join(dp, f); zf.write(full, os.path.relpath(full, tmp))
    n = sum(len(fs) for _, _, fs in os.walk(top)); shutil.rmtree(tmp)
    print(z, f'{os.path.getsize(z) / 1e6:.1f} MB，{n} 个文件')


if __name__ == '__main__':
    main()
