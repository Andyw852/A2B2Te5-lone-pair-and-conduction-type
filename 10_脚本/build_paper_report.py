"""由论文报告的 HTML 模板生成两样东西：

1. 单文件 HTML：{{FIG:名称}} 换成嵌入的图片（缩到 1800 px 宽的 PNG，base64），可直接打开或转发；
2. 压缩包：HTML 两份（图片用相对路径的一份；图片已嵌入的单文件一份）+ 全部图（PNG、PDF）+ 每幅图的源数据 + 正文表格（CSV），分目录存放：
       A2B2Te5_report_20261028/
         A2B2Te5_report.html     A2B2Te5_report_单文件.html     README.txt
         figures/png/  figures/pdf/
         data/Fig1_structures/ … data/Fig11_ordering/   （paper_figures.py 写出的源数据）
         data/tables/                                 （正文表 1–10，从 HTML 中抽出）
         data/L_layered_defects/                      （早先一批层状相缺陷计算的原始结果，取自 03_缺陷形成能与PN型数据/Pb2Sb2Te5.zip）

用法：先运行 10_脚本/paper_figures.py 生成图与源数据，再 python 10_脚本/build_paper_report.py
输入：21_阶段报告与交接/论文报告_20261028.src.html、21_阶段报告与交接/论文图_20261028/（含 source_data/）
输出：21_阶段报告与交接/论文报告_20261028.html、21_阶段报告与交接/A2B2Te5_report_20261028.zip
"""
import base64, csv, glob, io, os, re, shutil, tempfile, zipfile
from html.parser import HTMLParser
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, '21_阶段报告与交接')
FIGS = os.path.join(D, '论文图_20261028')
NAME = 'A2B2Te5_report_20261028'
src = open(os.path.join(D, '论文报告_20261028.src.html'), encoding='utf-8').read()


def embed(m):
    im = Image.open(os.path.join(FIGS, m.group(1) + '.png')).convert('RGB')
    if im.width > 1800:
        im = im.resize((1800, round(im.height * 1800 / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, 'PNG', optimize=True)
    return 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode()


class Tables(HTMLParser):
    """抽出带 data-csv 属性的表格；下标写成 _x，上标写成 ^x；表题取紧挨在前面的 <p class="tcap">。"""
    def __init__(self):
        super().__init__(); self.tables = []; self.cap = None; self.in_cap = False; self.cur = None; self.cell = None

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'p' and a.get('class') == 'tcap': self.cap, self.in_cap = '', True
        elif tag == 'table' and 'data-csv' in a: self.cur = {'name': a['data-csv'], 'caption': self.cap or '', 'rows': []}
        elif self.cur is not None and tag == 'tr': self.cur['rows'].append([])
        elif self.cur is not None and tag in ('td', 'th'): self.cell = ''
        elif tag in ('sub', 'sup') and self.cell is not None: self.cell += '_' if tag == 'sub' else '^'
        elif tag in ('sub', 'sup') and self.in_cap: self.cap += '_' if tag == 'sub' else '^'

    def handle_endtag(self, tag):
        if tag == 'p' and self.in_cap: self.in_cap = False
        elif tag == 'table' and self.cur is not None: self.tables.append(self.cur); self.cur = None; self.cap = None
        elif self.cur is not None and tag in ('td', 'th') and self.cell is not None:
            self.cur['rows'][-1].append(' '.join(self.cell.split())); self.cell = None

    def handle_data(self, data):
        if self.in_cap: self.cap += data
        elif self.cell is not None: self.cell += data


def readme(tables):
    t = '\n'.join(f'      {x["name"]}.csv    {" ".join(x["caption"].split())}' for x in tables)
    return f'''A2B2Te5（A = Pb, Sn；B = Sb, Bi）中的孤对电子与导电类型：论文写作用报告与数据
生成日期 2026-10-28。全部为理论计算；实验 S(T) 为用户提供的实验值。
结构：问题一（PBE、HSE 能带与 S）用所给的四种 9 原子层状结构。问题二只算 Pb2Sb2Te5，用了三套缺陷计算：
D（无序立方）、O（有序立方）为 144 原子岩盐超胞；L 为所给层状结构的 3×3×1 超胞（81 原子，早先一批）。见报告 1.1 节、图 1、表 1。

A2B2Te5_report.html    报告正文（用浏览器打开；图片从 figures/png 读取，请先整体解压再打开）
A2B2Te5_report_单文件.html  同一报告，图片已嵌入，单独打开或转发都能显示
figures/png/           图 1–11，600 dpi PNG（文件名的 FigN 即正文的图 N）
figures/pdf/           图 1–11，矢量 PDF（字体已嵌入，英文 Times New Roman，公式 STIX）
data/                  每幅图的源数据，即图中画出的数值，目录名与图文件名对应
  Fig1_structures/            O、D 两个 144 原子母胞的格点占位（每个格点的原子面、种类、坐标）与未弛豫 POSCAR；
                              所给的 9 原子层状结构 POSCAR（问题一的能带与 S 都用它）；图 1c、d 各原子面的高度与堆垛位置；
                              三种结构的晶格、体积与能量对比（structure_comparison.csv，即表 1）
  Fig2_bands_PBE_vs_HSE/      四种材料 PBE+SOC、HSE06+SOC 的路径能带（相对价带顶，eV）；路径带隙及 Γ 处最高价带、
                              最低导带的能量与反转带隙（band_gaps_on_path.csv）
  Fig3_HSE_bands_cation_s/    HSE06+SOC 路径能带及每个 k 点、每条带的阳离子 s 成分（%）
  Fig4_cation_s_by_atom/      HSE06+SOC 路径能带上每个 k 点、每条带的 A（Pb/Sn）、B（Sb/Bi）、Te 的 s 成分（%）；关键态（表 4）
  Fig5_SnSb_vs_SnBi/          Sn2Sb2Te5 与 Sn2Bi2Te5 路径能带上的 Sn-s、B-s、Sn-p、B-p、Te-p 投影；Γ 单谷与三组孤对电子谷的能量与成分；
                              Γ 处最高价带与最低导带的能量与成分（Gamma_states.csv）（表 5）
  Fig6_valence_pockets/       Sn2Sb2Te5、Sn2Bi2Te5 价带口袋在布里渊区中的位置（pocket_positions.csv：每个口袋一行，
                              分数坐标、笛卡尔坐标、所在面、能量、等价口袋数）；四种材料的价带能级（levels_HSE.csv）与谷间能量差（HSE06+SOC）
                              ΔE = E(Σ 孤对电子谷) − E(Γ–A 轴上最高价带态)（valley_offset_HSE.csv）
  Fig7_lone_pair_valley/      最高孤对电子价带谷相对价带顶的能量（meV）、所在 k 点与成分
  Fig8_Seebeck/               300 K、3×10^19 cm^-3 的面内 S：PBE+SOC（BoltzTraP2）与 HSE 估计；表 6 的逐谷电导份额
                              （valley_conductivity_shares.csv：每个谷的位置、离带边能量、阳离子 s 成分、份额；lone_pair_share_summary.csv：汇总）
                              lone_pair_share_vs_concentration.csv：孤对电子谷份额随空穴浓度的变化；
                              Sn2Bi2Te5_c3_valley_check.csv：SnBi 价带谷 c3 拟合的梯度、主质量、驻点位移与是否采用
  Fig9_S_vs_T/                实验 S(T)、计算 S(T)（HSE 与 PBE 能带）、拟合浓度与最大偏差
  Fig10_formation_energies/   12 类缺陷逐位点的中性形成能（富 Te、贫 Te），两端化学势
  Fig11_ordering/             有序化引起的形成能变化 ΔE_f = E_f(O) − E_f(D)，逐位点与中位数
  L_layered_defects/          早先一批在所给层状结构上（3×3×1，81 原子）的缺陷计算，即正文表 9、表 10 的 L 列：
                              L_defect_ranking_all_vertices.csv  5 个化学势顶点下全部缺陷在平衡费米能级处的形成能排名
                              L_formation_energy_results.json    原始结果（各电荷态的形成能、转变能级、化学势、载流子浓度）
  tables/                     正文表格：
{t}

约定：能量单位 eV（另注明者除外），S 单位 µV/K，载流子浓度单位 cm^-3（正为空穴、负为电子），
k 点坐标为倒格矢分数坐标；谷的名称为所在高对称点或高对称线的标准记号（Σ：Γ–M，R：A–L，Δ：Γ–A，Λ：Γ–K，
Q：A–H，T：M–K；同一条线上的几个谷按离 Γ–A 轴由近到远加下标，如 R₁、R₂），见报告 1.2 节；能带 CSV 中 band_N 为 VASP 能带序号（SOC，_VB_top 为价带顶，_CB_bottom 为导带底）。
CSV 为 UTF-8（带 BOM），可直接用 Excel、Origin 打开。

生成脚本（在仓库中）：10_脚本/paper_figures.py（图与源数据）、10_脚本/multivalley_seebeck.py（S 的模型）、
10_脚本/build_paper_report.py（本压缩包）。
'''.replace('\n', '\r\n')


def package():
    p = Tables(); p.feed(src)
    tmp = tempfile.mkdtemp(); top = os.path.join(tmp, NAME)
    for ext in ('png', 'pdf'):
        os.makedirs(os.path.join(top, 'figures', ext))
        for f in sorted(glob.glob(os.path.join(FIGS, f'Fig*.{ext}'))):
            shutil.copy(f, os.path.join(top, 'figures', ext))
    shutil.copytree(os.path.join(FIGS, 'source_data'), os.path.join(top, 'data'))
    L = os.path.join(top, 'data', 'L_layered_defects'); os.makedirs(L)
    with zipfile.ZipFile(os.path.join(ROOT, '03_缺陷形成能与PN型数据', 'Pb2Sb2Te5.zip')) as zold:
        for suffix, name in (('全部顶点缺陷排名.csv', 'L_defect_ranking_all_vertices.csv'),
                             ('formation_energy_results.json', 'L_formation_energy_results.json')):
            src_name = next(n for n in zold.namelist() if n.endswith(suffix))
            open(os.path.join(L, name), 'wb').write(zold.read(src_name))
    os.makedirs(os.path.join(top, 'data', 'tables'))
    for t in p.tables:
        with open(os.path.join(top, 'data', 'tables', t['name'] + '.csv'), 'w', newline='', encoding='utf-8-sig') as f:
            csv.writer(f).writerows(t['rows'])
    html = re.sub(r'\{\{FIG:([^}]+)\}\}', lambda m: f'figures/png/{m.group(1)}.png', src)
    open(os.path.join(top, 'A2B2Te5_report.html'), 'w', encoding='utf-8').write(html)
    open(os.path.join(top, 'A2B2Te5_report_单文件.html'), 'w', encoding='utf-8').write(standalone)
    open(os.path.join(top, 'README.txt'), 'w', encoding='utf-8-sig', newline='').write(readme(p.tables))
    out = os.path.join(D, NAME + '.zip')
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for dp, _, fs in sorted(os.walk(top)):
            for f in sorted(fs):
                full = os.path.join(dp, f); z.write(full, os.path.relpath(full, tmp))
    n = sum(len(fs) for _, _, fs in os.walk(top)); shutil.rmtree(tmp)
    print(out, f'{os.path.getsize(out) / 1e6:.1f} MB，{n} 个文件，表格 {len(p.tables)} 个')


standalone = re.sub(r'\{\{FIG:([^}]+)\}\}', embed, src)
out = os.path.join(D, '论文报告_20261028.html')
open(out, 'w', encoding='utf-8').write(standalone)
print(out, f'{len(standalone) / 1e6:.1f} MB')
package()
