"""由论文报告的 HTML 模板生成两样东西：

1. 单文件 HTML：{{FIG:名称}} 换成嵌入的图片（缩到 1800 px 宽的 PNG，base64），可直接打开或转发；
2. 压缩包：HTML 两份（图片用相对路径的一份；图片已嵌入的单文件一份）+ 全部图（PNG、PDF）+ 每幅图的源数据 + 正文表格（CSV），分目录存放：
       A2B2Te5_report_20261028/
         A2B2Te5_report.html     A2B2Te5_report_单文件.html     README.txt
         figures/png/  figures/pdf/
         data/<类别>/<图目录>/    （paper_figures.py 写出的源数据，图 1–12（含 -1/-2 两幅）、S1–S3、表 6–8，按 plot_data_index.CATS 分类）
         data/06_正文表格/tables/                       （正文表 1–11，从 HTML 中抽出）
         data/05_缺陷与导电类型/L_layered_defects/      （早先一批层状相缺陷计算的原始结果，取自 03_缺陷形成能与PN型数据/Pb2Sb2Te5.zip）
         data/作图数据索引.csv                          （每幅图、每个子图用哪个文件、横轴纵轴取哪一列；10_脚本/plot_data_index.py）
         作图数据/<类别>/<图目录>.xlsx                  （每幅图单独一个 Excel，xlsxwriter 写出并读回检查：第一页说明，其后每个数据表一页）

用法：先运行 10_脚本/paper_figures.py 生成图与源数据，再 python 10_脚本/build_paper_report.py
输入：21_阶段报告与交接/论文报告_20261028.src.html、21_阶段报告与交接/论文图_20261028/（含 source_data/）
输出：21_阶段报告与交接/论文报告_20261028.html、21_阶段报告与交接/A2B2Te5_report_20261028.zip
"""
import base64, csv, glob, io, os, re, shutil, sys, tempfile, zipfile
from html.parser import HTMLParser
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_data_index
D = os.path.join(ROOT, '21_阶段报告与交接')
FIGS = os.path.join(D, '论文图_20261028')
NAME = 'A2B2Te5_report_20261028'
src = open(os.path.join(D, '论文报告_20261028.src.html'), encoding='utf-8').read()
src = src.replace('{{PLOTMETHODS}}', plot_data_index.methods_html())      # 每幅图的绘制方法，来自 plot_data_index.METHODS


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
生成日期 2026-10-28，2026-11-05 整合版：HSE06+SOC 能带专题并入本报告；HSE06+SOC 与 PBE+SOC 的 S 都是同一 28×28×7/6 网格上
BoltzTraP2 的直接计算（不加剪刀）；新增 COHP（图 3-2 为 HSE06，图 3-3 为 HSE06 与 PBE 对比，另用 LOBSTER 交叉检验）。改动见报告文首的版本说明。计算全部为理论计算；实验 S(T)、σ(T) 与室温霍尔为用户提供的实验值。
结构：问题一（能带、COHP 与 S）用所给的四种 9 原子层状结构。问题二只算 Pb2Sb2Te5，用了三套缺陷计算：
D（无序立方）、O（有序立方）为 144 原子岩盐超胞；L 为所给层状结构的 3×3×1 超胞（81 原子，早先一批）。见报告 1.1 节、图 1、表 1。

A2B2Te5_report.html    报告正文（用浏览器打开；图片从 figures/png 读取，请先整体解压再打开）
A2B2Te5_report_单文件.html  同一报告，图片已嵌入，单独打开或转发都能显示
figures/png/           全部图，600 dpi PNG（文件名的 FigN、FigN-1、FigN-2 即正文的图 N、N-1、N-2）
figures/pdf/           同上，矢量 PDF（字体已嵌入，英文 Times New Roman，公式 STIX）
                       输运图都分两幅：“-1”只画 HSE06+SOC，“-2”为 HSE06+SOC 与 PBE+SOC 对比（图 8、9、10）
作图数据/<类别>/       每幅图单独一个 Excel（文件名同 data/ 下的图目录名；第一页“说明”写明每个子图用哪一页、横轴纵轴取哪一列，
                       其后每个数据表一页，数值为数字格式，可直接在 Excel、Origin 中作图）；正文表格在 06_正文表格/tables.xlsx
data/作图数据索引.csv  同一份索引的 CSV（全部图）
data/<类别>/<图目录>/   每幅图的源数据，即图中画出的数值，图目录名与图文件名对应。类别：
                       01_结构（图 1）；02_能带与轨道成分（图 2-1 至 5、S2、S3）；03_价带谷与口袋（图 6-1、6-2、7、S1）；
                       04_输运与实验对照（图 8-1 至 10-2、表 6–8）；05_缺陷与导电类型（图 11、12，L 层状相缺陷）；06_正文表格（表 1–11）
                       下面按图目录逐个说明（能量相对真实价带顶，另注明者除外）：
  01_结构/
  Fig1_structures/            O、D 两个 144 原子母胞的格点占位与未弛豫 POSCAR；所给的 9 原子层状结构 POSCAR；
                              图 1c、d 各原子面的高度与堆垛位置；三种结构的晶格、体积与能量对比（structure_comparison.csv，即表 1）
  02_能带与轨道成分/
  Fig2-1_HSE_bands/           HSE06+SOC 路径能带；Sn 系的附加路径 Γ–Σ′–R₂′–L（*_true_valley_path_G-S-R2-L.csv）；真实带隙与路径带隙
  Fig2-2_bands_PBE_vs_HSE/    PBE+SOC、HSE06+SOC 的路径能带（各相对本泛函的价带顶）；路径带隙及 Γ 处能级与反转带隙（band_gaps_on_path.csv）
  Fig3-1_HSE_bands_cation_s/  HSE06+SOC 路径能带及每个 k 点、每条带的阳离子 s 成分（%）；Sn 系的附加路径
  Fig3-2_COHP_HSE/            HSE06（不加 SOC，COGITO）最短 A–Te、B–Te 键的轨道分辨 −COHP（阳离子 s–Te p、阳离子 p–Te p、总和；正为成键，
                              能量相对价带顶）；ICOHP（全部占据态与价带顶以下 1 eV 内，ICOHP_summary.csv，正为反键）
  Fig3-3_COHP_HSE_vs_PBE/     同一批键在 HSE06 与 PBE 下的 −COHP 与 ICOHP 的变化；LOBSTER（PBE）的交叉检验（按能量窗口的 s–p ICOHP）
  Fig4-1_cation_s_by_atom/    HSE06+SOC 路径能带上每个 k 点、每条带的 A（Pb/Sn）、B（Sb/Bi）、Te 的 s 成分（%）；关键态（表 4）
  Fig4-2_HSE_orbitals/        HSE06+SOC 路径能带上的 A 位 s、B 位 s、A 位 p、B 位 p、Te p 成分（图 4-2 与图 S3a–d 共用）
  Fig5_SnSb_vs_SnBi/          Sn2Sb2Te5 与 Sn2Bi2Te5 的轨道投影、价带谷、价带能级（levels_true.csv 为真实值；其余相对路径价带顶）、
                              Γ 处最高价带与最低导带的能量与成分（表 5）
  FigS2_Gamma_levels_SOC/     Γ 点能级（PBE，有、无 SOC）与 Sb → Bi 变化的标量 / SOC 拆分（10_脚本/soc_mechanism.py）
  03_价带谷与口袋/
  Fig6-1_valence_levels/      四种材料价带上的谷：能量、等价口袋数、k、阳离子 s、主轴有效质量（HSE06+SOC）
  Fig6-2_pockets_BZ/          Sn 系全部价带口袋在布里渊区中的位置（pocket_positions.csv：每个口袋一行）
  Fig7_lone_pair_valley/      最高孤对电子价带谷相对价带顶的能量（meV）、所在 k 点与成分（两种泛函）
  FigS1_mirror_plane/         Sn 系镜面内最高价带的能量网格（meV）；路径上的谷与真实极大；两者的能量差（path_vs_true_maxima.csv）
  04_输运与实验对照/          （BoltzTraP2：HSE06+SOC 与 PBE+SOC，同一 28×28×7/6 网格，不加剪刀，常数弛豫时间）
  Fig8-1_Seebeck_HSE/         300 K 面内 S 随载流子浓度（HSE06+SOC）
  Fig8-2_Seebeck_HSE_vs_PBE/  同上，两种泛函（按 functional 列区分）
  Table6_S_300K/              表 6：300 K、±3×10¹⁹ 与 ±1×10²⁰ cm⁻³ 的 S
  Fig9-1_S_vs_T_HSE/          四个样品的实验 T、σ（S/m）、S 与室温霍尔；霍尔浓度下直接计算的 S(T)（HSE06+SOC，不调参数）
  Fig9-2_S_vs_T_HSE_vs_PBE/   同上，计算含两种泛函
  Table7_S_T_at_Hall/         表 7：最低、最高测量温度的实验与计算值，300–620 K 的均方根偏差
  Table8_Pisarenko_mass/      表 8（不作图）：样品的霍尔浓度、实验 S，同一浓度下 HSE、PBE 的 S，以及由此读出的态密度有效质量
  Fig10-1_Fermi_window_HSE/   300 K 的化学势 µ（霍尔浓度与 3×10¹⁹ cm⁻³，HSE06+SOC；µ 取自 BoltzTraP2，用 Γ 处本征值对齐到能带路径）
  Fig10-2_Fermi_window_HSE_vs_PBE/  两种泛函各自在霍尔浓度下的 µ
  05_缺陷与导电类型/
  Fig11_formation_energies/   12 类缺陷逐位点的中性形成能（富 Te、贫 Te），两端化学势
  Fig12_ordering/             有序化引起的形成能变化 ΔE_f = E_f(O) − E_f(D)，逐位点与中位数
  L_layered_defects/          早先一批在所给层状结构上（3×3×1，81 原子）的缺陷计算，即正文表 9、表 10 的 L 列：
                              L_defect_ranking_all_vertices.csv  5 个化学势顶点下全部缺陷在平衡费米能级处的形成能排名
                              L_formation_energy_results.json    原始结果（各电荷态的形成能、转变能级、化学势、载流子浓度）
  06_正文表格/
  tables/                     正文表格：
{t}

约定：能量单位 eV（另注明者除外），S 单位 µV/K，载流子浓度单位 cm^-3（正为空穴、负为电子），
k 点坐标为倒格矢分数坐标；谷的名称为所在高对称点或高对称线的标准记号（Σ：Γ–M，R：A–L，Δ：Γ–A，Λ：Γ–K，
Q：A–H，T：M–K；同一条线上的几个谷按离 Γ–A 轴由近到远加下标，如 R₁、R₂），见报告 1.2 节；能带 CSV 中 band_N 为 VASP 能带序号（SOC，_VB_top 为价带顶，_CB_bottom 为导带底）。
CSV 为 UTF-8（带 BOM），可直接用 Excel、Origin 打开。

生成脚本（在仓库中）：10_脚本/paper_figures.py（图与源数据；输运图 transport_figures.py，HSE 能带与轨道图 hse_band_figures.py，
COHP 提取 cohp_extract.py --run pbe/hse/lobster）；BoltzTraP2 结果 22_能带基准_1-0/<材料>/c8_hse_dense/gene/bt2_mult60/（HSE06+SOC）、
22_能带基准_1-0/bt2_pbe_dense/（PBE+SOC）；10_脚本/build_paper_report.py（本压缩包）。
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
    nfile, nbook, nsheet = plot_data_index.build(os.path.join(top, 'data'), os.path.join(top, '作图数据'), plot_data_index.CATS)
    print(f'作图数据：{nfile} 个文件全部在索引里；{nbook} 个 Excel、{nsheet} 个工作表')
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
