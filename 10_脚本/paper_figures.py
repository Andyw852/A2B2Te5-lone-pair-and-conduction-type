"""Journal-style figures for the manuscript (English labels; 600-dpi PNG + vector PDF), numbered as in 论文报告_20261028.

Fig. 1  Structures: O/D 144-atom supercells (a, b) and side views of the user's layered phase vs model O (c, d)
                                                                               补算实施方案_20261007/.../parents/<id>/site_map.json;
                                                                               A2B2Te5_SOC_dense/materials/Pb2Sb2Te5/POSCAR; 13_…/energies.csv
Fig. 2  PBE+SOC vs HSE06+SOC band structures (c2 path)                       22_能带基准_1-0/<mat>/c2_{pbe,hse}/path_eigen.txt
Fig. 3  HSE06+SOC bands weighted by cation-s (lone-pair) character             22_能带基准_1-0/<mat>/c2_hse/PROCAR_OPT.gz
Fig. 4  HSE06+SOC bands: which cation (A = Pb/Sn or B = Sb/Bi) carries the s weight   same
Fig. 5  Sn2Sb2Te5 vs Sn2Bi2Te5: Sn-5s and B-site p (Sb 5p / Bi 6p) projections near the VBM, and level diagram   same
Fig. 6  Valence pockets in the Brillouin zone (Sn2Sb2Te5 vs Sn2Bi2Te5) and the valence-valley offset ΔE = E(Σ) − E(Γ–A), HSE only
                                                                               same + A2B2Te5_SOC_dense/materials/<mat>/POSCAR (spglib)
Fig. 7  Energy of the highest lone-pair valence valley, PBE vs HSE            same
Fig. 8  Pisarenko plots at 300 K (HSE estimate | PBE BoltzTraP2; p | n) with the samples at their Hall concentrations
                                                                               10_脚本/multivalley_seebeck.py; boltztrap_crta.json; 21_…/实验数据/
        (+ Table 6 at 3e19 cm-3, per-valley conductivity shares and cation-s character)
Fig. 9  S(T) of all four samples: experiment vs calculation at the room-temperature Hall concentration (no free parameter)
                                                                               21_阶段报告与交接/实验数据/（raw T–σ–S, Hall）;
                                                                               A2B2Te5_SOC_dense/materials/<mat>/.../boltztrap_crta.json
Fig. 10 Neutral formation energies of 12 defect types in O and D               13_DVO缺陷_VASP流程与分析/结果_第一批/analysis/形成能_逐位点.csv
Fig. 11 Ordering-induced change of the formation energy, O − D                  same
Fig. S1 Top valence band in the mirror plane (c4 grid) for Sn2Sb2Te5 / Sn2Bi2Te5, with k_z cuts   22_能带基准_1-0/<mat>/c4_mirror_v3/b0*/PROCAR_OPT.gz

Valleys off the high-symmetry points are named after the line they lie on (standard hexagonal notation, Bradley & Cracknell):
Γ–M Σ, M–K T, K–Γ Λ, Γ–A Δ, A–L R, L–H S, H–A Q; the two A–L valleys are R1 (near A) and R2 (near L). See kname().

Every figure also writes its plotted numbers as CSV (UTF-8 with BOM) to --data/<figure name>/.

Style: Times New Roman (Liberation Serif, same metrics, where TNR is not installed), STIX math, 7–8 pt text, 0.8-pt axes,
inward ticks on all four sides, panel labels (a), (b), …, panel titles above the axes;
double column 7.0 in, single column 3.35 in; Okabe–Ito colors (HSE #0072B2, PBE #D55E00, ordered #009E73, disordered #CC79A7).

Usage: python 10_脚本/paper_figures.py [--out 21_阶段报告与交接/论文图_20261028] [--data <out>/source_data] [--only 1 2 ...]
"""
import argparse, csv, json, os, shutil, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import spglib
from scipy import integrate, optimize

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import c2_valleys as cv
import multivalley_seebeck as mv
import report_figures_20261028 as rf

MATS = rf.MATS
TEX = {m: m.replace('2', '$_2$').replace('5', '$_5$') for m in MATS}
NEL = rf.NEL
C_HSE, C_PBE, C_O, C_D = '#0072B2', '#D55E00', '#009E73', '#CC79A7'
C_PB, C_SB, C_TE = '#333333', '#E69F00', '#56B4E9'
GRID = '#D9D9D9'
LAB = rf.LAB
W2, W1 = 7.0, 3.35
KEND = 187          # 路径截到第二个 A（标准六方路径 Γ–M–K–Γ–A–L–H–A），区界上沿 c 的 L–M、H–K 两小段不画
DASH = (0, (3, 1.6))
DATA = None         # 源数据目录；None 时不写


def style():
    # 英文 Times New Roman；本机没有时用度量完全相同的 Liberation Serif 代替
    for f in font_manager.findSystemFonts():
        n = os.path.basename(f).lower().replace('-', '').replace(' ', '').replace('_', '')
        if n.startswith(('times', 'liberationserif')):
            font_manager.fontManager.addfont(f)
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman', 'Liberation Serif', 'DejaVu Serif'],
        'mathtext.fontset': 'stix',
        'font.size': 7.5, 'axes.labelsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'legend.fontsize': 7, 'axes.linewidth': 0.8, 'lines.linewidth': 1.0,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
        'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5, 'legend.frameon': False,
        'savefig.dpi': 600, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'axes.unicode_minus': True, 'figure.facecolor': 'white'})


def panel(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')


def save(fig, out, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(out, f'{name}.{ext}'))
    plt.close(fig)


def dump(sub, name, header, rows):
    if DATA is None: return
    d = os.path.join(DATA, sub); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def r(v, n):
    return round(float(v), n)


def sci(N, unit=True):
    a, e = f'{abs(N):.1e}'.split('e')
    return rf'{a}$\times10^{{{int(e)}}}$' + (r' cm$^{-3}$' if unit else '')


def kticks(ax, x):
    t, l = [], []
    for i in sorted(LAB):
        if i > KEND: break
        t.append(x[i]); l.append(LAB[i])
    for v in t[1:-1]: ax.axvline(v, color=GRID, lw=0.5, zorder=0)
    ax.set_xticks(t, l); ax.set_xlim(x[0], x[KEND]); ax.tick_params(axis='x', top=False, length=0, pad=3)


def vmark(ax, xp, yp, name, lead=None):
    """在能带图上标出一个谷。默认在谷顶上方画 ▼ 和名称；lead=(tx, ty) 时把名称放在 (tx, ty)，用细线指向谷顶，
    用于 Γ–A 这类很窄、两个标记会挤在一起的地方。"""
    lab = MATH.get(name, name)
    if lead is None:
        ax.plot(xp, yp + 0.035, 'v', ms=3.2, color='k', zorder=4)
        ax.text(xp, yp + 0.062, lab, ha='center', va='bottom', fontsize=7, zorder=4, bbox=dict(fc='white', ec='none', pad=0.3))
    else:
        ax.annotate(lab, xy=(xp, yp + 0.01), xytext=lead, ha='center', va='bottom', fontsize=7, zorder=4,
                    bbox=dict(fc='white', ec='none', pad=0.3),
                    arrowprops=dict(arrowstyle='-|>', lw=0.5, color='k', mutation_scale=5, shrinkA=0, shrinkB=0))


def kpath():
    return cv.read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))


def bands(m, side):
    E = cv.read_path_eigen(os.path.join(rf.base(m), side, 'path_eigen.txt'))
    nel = NEL[m]; vbm = np.nanmax(E[:, nel - 1]); gap = np.nanmin(E[:, nel]) - vbm
    return E - vbm, gap


# 路径上离开对称点的位置按所在的高对称线命名（六方布里渊区标准记号，Bradley & Cracknell）
SEGS = [(0, 33, 'Σ'), (33, 52, 'T'), (52, 90, 'Λ'), (90, 97, 'Δ'), (97, 130, 'R'), (130, 149, 'S'), (149, 187, 'Q')]
MATH = {'R₁': r'R$_1$', 'R₂': r'R$_2$', 'R₁′': r'R$_1$′', 'R₁″': r'R$_1$″', 'R₂′': r'R$_2$′'}   # 图中用 mathtext 写下标


def kname(i):
    """c2 路径序号 → 位置名：对称点用点名；线上用线名，A–L（R）线按离 A 的远近分 R₁（u < 1/6）、R₂。"""
    if i in LAB and i <= KEND: return LAB[i]
    for a, b, n in SEGS:
        if a < i < b:
            return n if n != 'R' else ('R₁' if (i - a) / (2 * (b - a)) < 1 / 6 else 'R₂')
    return f'k{i}'


AXD = {'Σ': lambda i: i, 'Λ': lambda i: 90 - i, 'R': lambda i: i - 97, 'Q': lambda i: 187 - i}   # 离 Γ–A 轴的远近（路径步数）
SUB = '₁₂₃₄'


def knames(idx):
    """一组谷（同一材料、同一条带）的位置名。同一条线上有几个谷时按离 Γ–A 轴由近到远加下标；R 线始终用 R₁、R₂（与表 5 一致）。"""
    idx = [i for i in idx if i is not None]
    base = {i: kname(i) for i in idx}; out = dict(base)
    for n in set(base.values()):
        same = sorted([i for i in idx if base[i] == n], key=AXD.get(n, lambda i: i))
        if len(same) > 1 and n in AXD: out.update({i: n + SUB[k] for k, i in enumerate(same)})
    return out


def bandname(b, nel):
    return f'band_{b + 1}' + ('_VB_top' if b == nel - 1 else '_CB_bottom' if b == nel else '')


# ------------------------------------------------------------------ Fig. 1
SUPER = os.path.join(ROOT, '补算实施方案_20261007', 'rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED', 'parents')
LAYERED = os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', 'Pb2Sb2Te5', 'POSCAR')   # 用户结构，问题一的能带与 S 都用它
ENERGIES = os.path.join(ROOT, '13_DVO缺陷_VASP流程与分析', '结果_第一批', 'energies.csv')
ZMAX = 11.3


def read_poscar(path):
    L = [l.split() for l in open(path).read().splitlines()]
    lat = np.array(L[2:5], float) * float(L[1][0]); sym, cnt = L[5], [int(c) for c in L[6]]
    frac = np.array([l[:3] for l in L[8:8 + sum(cnt)]], float)
    return lat, frac, [e for e, c in zip(sym, cnt) for _ in range(c)]


def host_energies():
    """13_ 的母胞与参考相：每 f.u. 的能量与体积（PBE-D3、无 SOC，缺陷计算的同一套设置）。"""
    E = {}
    for row in csv.DictReader(open(ENERGIES, encoding='utf-8-sig')):
        if row['status'] == 'DONE' and row['kind'] in ('parent', 'ref'):
            E[row['job'].split('/')[1]] = (float(row['energy_eV']), float(row['volume_A3']), int(row['n_atoms']))
    fu = {'O_ordered': 16, 'D_seed1729': 16, 'Pb2Sb2Te5_layered9': 1}
    hull = E['PbTe'][0] * 2 + E['Sb2Te3'][0]          # 2 PbTe + Sb2Te3，每 f.u. Pb2Sb2Te5
    return {k: dict(E=E[k][0] / n, V=E[k][1] / n, hull=E[k][0] / n - hull) for k, n in fu.items()}


def planes_layered():
    """层状相各原子面：高度（以两 Pb 面之间的中心 Te 为 0）、面内堆垛位置 0/1/2（单位 p/3，p = a√3/2）。"""
    lat, frac, sp = read_poscar(LAYERED)
    a, c = np.linalg.norm(lat[0]), lat[2, 2]; p = a * np.sqrt(3) / 2
    out = []
    for m in (-1, 0, 1):
        for f, e in zip(frac, sp):
            z = (f[2] + m) * c
            out.append((z, e, round((-(f @ lat)[1] % p) / p * 3) % 3))
    out = sorted({(round(z, 4), e, w) for z, e, w in out if abs(z) <= ZMAX})
    return out, a, c, p


def planes_O(a_eff):
    """有序模型 O：理想岩盐格点（a = 6.2 Å 的起始坐标）按弛豫后的体积缩放到 a_eff；仓库里没有弛豫后的坐标。"""
    d = json.load(open(os.path.join(SUPER, 'O_ordered', 'site_map.json')))
    n = np.ones(3) / np.sqrt(3); u = np.array([1, 1, -2]) / np.sqrt(6)
    L3 = np.array(d['lattice_vectors_A_rows'])[2]; a0 = 6.2; p0 = a0 / np.sqrt(2) * np.sqrt(3) / 2
    pl = {}
    for s in d['sites']:
        r_ = np.array(s['cartesian_coordinates_A'])
        pl[s['q_times_2_mod_10']] = (r_ @ n, 'Vac' if s['is_structural_ghost_vacancy'] else s['occupancy'], (r_ @ u) % p0 / p0 * 3)
    zc = pl[3][0]                                    # 第 3 面：两 Pb 面之间的中心 Te
    assert pl[3][1] == 'Te' and pl[2][1] == pl[4][1] == 'Pb'
    out = []
    for m in (-1, 0, 1):
        for h, e, w in pl.values():
            z = (h + m * (L3 @ n) - zc) * a_eff / a0
            if abs(z) <= ZMAX: out.append((round(z, 4), e, round(w + m * (L3 @ u) / p0 * 3 - pl[3][2]) % 3))
    return sorted(out), a_eff / np.sqrt(2) * np.sqrt(3) / 2


def side_view(ax, planes, p, title, sub):
    ax.axhspan(min(z for z, e, _ in planes if abs(z) < 8 and e == 'Te'), max(z for z, e, _ in planes if abs(z) < 8 and e == 'Te'),
               color='#EFEFEF', lw=0, zorder=0)
    for z, e, w in planes:
        for k in range(3):
            x = (w / 3 + k) * p
            if e == 'Te':
                ax.plot(x, z, 'o', ms=3.6, color=C_TE, mec='none', zorder=3)
            elif e == 'Vac':
                ax.plot(x, z, 'o', ms=5.0, mfc='white', mec='k', mew=0.7, zorder=3)
            else:
                ax.plot(x, z, 'o', ms=5.0, color=C_PB if e == 'Pb' else C_SB, mec='white', mew=0.4, zorder=3)
        ax.text(-0.9, z, e, ha='right', va='center', fontsize=6)
    ax.set_xlim(-3.0, 3 * p + 5.6); ax.set_ylim(-ZMAX - 0.7, ZMAX + 0.7); ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([]); ax.tick_params(top=False, right=False)
    ax.set_title(title, fontsize=7.5, pad=3); ax.set_xlabel(sub, fontsize=6.5, labelpad=4)


def gap(ax, x, z0, z1, text):
    ax.annotate('', xy=(x, z0), xytext=(x, z1), arrowprops=dict(arrowstyle='<->', lw=0.6, shrinkA=0, shrinkB=0))
    ax.text(x + 0.35, (z0 + z1) / 2, text, ha='left', va='center', fontsize=6, linespacing=1.15)


def fig1(out):
    H = host_energies()
    a_eff = {k: (H[k]['V'] * 16 / 20) ** (1 / 3) for k in ('O_ordered', 'D_seed1729')}   # 160 个岩盐格点 = 20 个立方晶胞
    fig = plt.figure(figsize=(W2, 5.75))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, 3.35], hspace=0.2, bottom=0.075, top=0.97)
    axs = [fig.add_subplot(g) for g in gs[0].subgridspec(1, 2, wspace=0.08)]
    for ax, pid, title in zip(axs, ('O_ordered', 'D_seed1729'), ('Ordered model O (after annealing)', 'Disordered model D (before annealing)')):
        sites = json.load(open(os.path.join(SUPER, pid, 'site_map.json')))['sites']
        rows = []
        for q in range(10):
            row = sorted([s for s in sites if s['q_times_2_mod_10'] == q], key=lambda s: tuple(s['replica']))
            y = 9 - q
            for i, s in enumerate(row):
                sp = 'Vac' if s['is_structural_ghost_vacancy'] else s['occupancy']
                rows.append([s['site_id'], y, i, s['sublattice'], sp] + [r(v, 5) for v in s['fractional_coordinates']]
                            + [r(v, 4) for v in s['cartesian_coordinates_A']])
                if s['sublattice'] == 'anion':
                    ax.plot(i, y, 'o', ms=3.2, color=C_TE, mec='none')
                elif sp == 'Vac':
                    ax.plot(i, y, 'o', ms=5.2, mfc='white', mec='k', mew=0.7)
                else:
                    ax.plot(i, y, 'o', ms=5.2, color=C_PB if sp == 'Pb' else C_SB, mec='white', mew=0.4)
        dump('Fig1_structures', f'{pid}_sites.csv',
             ['site_id', 'plane_from_bottom_along_111', 'position_in_plane_row', 'sublattice', 'species',
              'frac_a', 'frac_b', 'frac_c', 'x_A', 'y_A', 'z_A'], sorted(rows, key=lambda t: (t[1], t[2])))
        if DATA is not None:
            shutil.copy(os.path.join(SUPER, pid, 'POSCAR'), os.path.join(DATA, 'Fig1_structures', f'{pid}_POSCAR_unrelaxed.vasp'))
        ax.set_xlim(-0.8, 15.8); ax.set_ylim(-0.8, 9.8); ax.set_xticks([]); ax.set_yticks([])
        ax.tick_params(top=False, right=False)
        ax.set_title(title, fontsize=8, pad=3)
    axs[0].annotate('', xy=(-2.9, 9.4), xytext=(-2.9, -0.4), xycoords='data', annotation_clip=False,
                    arrowprops=dict(arrowstyle='->', lw=0.8))
    axs[0].text(-3.3, 4.5, '[111]', rotation=90, ha='right', va='center', fontsize=7)
    for y, lab in zip(range(10), ['Sb', 'Te', 'Pb', 'Te', 'Pb', 'Te', 'Sb', 'Te', 'Vac', 'Te'][::-1]):
        axs[0].text(-1.1, y, lab, ha='right', va='center', fontsize=6.5)
    panel(axs[0], '(a)', x=-0.09); panel(axs[1], '(b)')

    # (c)(d) 侧视：用户的层状相与有序模型 O
    PL, aL, cL, pL = planes_layered()
    PO, pO = planes_O(a_eff['O_ordered'])
    bot = gs[1].subgridspec(1, 2, wspace=0.0)
    ac, ad = fig.add_subplot(bot[0]), fig.add_subplot(bot[1])
    HL, HO = H['Pb2Sb2Te5_layered9'], H['O_ordered']
    side_view(ac, PL, pL, 'Layered phase, $P\\bar{3}m1$, 9 atoms (user structure)\nused for bands, $S$ (Figs. 2–9) and defect set L',
              f'$a_\\mathrm{{hex}}$ = {aL:.3f} Å,  $V$ = {HL["V"]:.1f} Å$^3$/f.u.,  $\\Delta E$ = 0')
    side_view(ad, PO, pO, 'Ordered model O, 144 atoms (rocksalt supercell)\nhost for the defect calculations (Figs. 10, 11)',
              f'$a_\\mathrm{{hex}}^\\mathrm{{eq}}$ = {a_eff["O_ordered"] / np.sqrt(2):.3f} Å,  $V$ = {HO["V"]:.1f} Å$^3$/f.u.,  '
              f'$\\Delta E$ = +{HO["E"] - HL["E"]:.3f} eV/f.u.')
    ac.set_ylabel('$z$ (Å)')
    ad.tick_params(labelleft=True)
    tL = sorted(z for z, e, _ in PL if e == 'Te' and z > 0)
    gL = (tL[-2], tL[-1]); gap(ac, 3 * pL + 0.5, *gL, f'vdW gap\n{gL[1] - gL[0]:.2f} Å')
    vO = [z for z, e, _ in PO if e == 'Vac' and z > 0][0]; tO = sorted(z for z, e, _ in PO if e == 'Te' and z > 0)
    gO = (max(t for t in tO if t < vO), min(t for t in tO if t > vO))
    gap(ad, 3 * pO + 0.5, *gO, f'vacancy\nplane\n{gO[1] - gO[0]:.2f} Å')
    ac.text(3 * pL + 0.85, 0, '9-plane\nblock', ha='left', va='center', fontsize=6, color='#555555')
    ad.text(3 * pO + 0.85, 0, '9-plane\nblock', ha='left', va='center', fontsize=6, color='#555555')
    panel(ac, '(c)', x=0.02); panel(ad, '(d)', x=0.02)

    h = [plt.Line2D([], [], ls='', marker='o', ms=5, color=C_PB), plt.Line2D([], [], ls='', marker='o', ms=5, color=C_SB),
         plt.Line2D([], [], ls='', marker='o', ms=5, mfc='white', mec='k'), plt.Line2D([], [], ls='', marker='o', ms=3.5, color=C_TE)]
    fig.legend(h, ['Pb', 'Sb', 'Structural vacancy', 'Te'], loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.01))

    if DATA is not None:
        shutil.copy(LAYERED, os.path.join(DATA, 'Fig1_structures', 'layered_Pb2Sb2Te5_POSCAR_user.vasp'))
    dump('Fig1_structures', 'side_view_planes.csv',
         ['structure', 'z_A_relative_to_central_Te', 'species', 'stacking_position_0_1_2', 'geometry'],
         [['layered_P-3m1_9atom', r(z, 3), e, w, 'user POSCAR (as used for bands)'] for z, e, w in PL]
         + [['O_ordered', r(z, 3), e, w, f'ideal rocksalt sites scaled to relaxed volume (a_eff = {a_eff["O_ordered"]:.4f} A)'] for z, e, w in PO])
    VL = abs(np.linalg.det(read_poscar(LAYERED)[0]))
    dump('Fig1_structures', 'structure_comparison.csv',
         ['structure', 'used_for', 'how_obtained', 'atoms_in_cell', 'structural_vacancies',
          'a_hex_A', 'V_per_fu_A3', 'E_per_fu_eV', 'E_minus_layered_eV_per_fu', 'E_minus_2PbTe+Sb2Te3_eV_per_fu',
          'vertical_Te-Te_spacing_across_gap_or_vacancy_plane_A', 'note'],
         [['layered_P-3m1_9atom', 'bands (PBE+SOC, HSE06+SOC), orbital character, topology, Seebeck; reference phase in defect study',
           "user structure (PBE-D3); for the energy re-relaxed with the defect settings (ISIF=3)", 9, 0,
           r(aL, 4), r(HL['V'], 2), r(HL['E'], 5), 0.0, r(HL['hull'], 4), r(gL[1] - gL[0], 3),
           f'a_hex and spacing from the user POSCAR (V = {VL:.2f} A3/f.u.); V and E from the re-relaxation'],
          ['O_ordered', 'host for defects, ordered (after annealing)',
           'built from ideal rocksalt (a = 6.2 A start) following Evang et al. 2022; PBE-D3 full relaxation (ISIF=3)', 144, 16,
           r(a_eff['O_ordered'] / np.sqrt(2), 4), r(HO['V'], 2), r(HO['E'], 5), r(HO['E'] - HL['E'], 4), r(HO['hull'], 4),
           r(gO[1] - gO[0], 3), 'a_hex equivalent from volume, a_eff/sqrt(2) with a_eff = (V/20)^(1/3); spacing for ideal sites at a_eff'],
          ['D_seed1729', 'host for defects, disordered (before annealing)',
           'built from ideal rocksalt (a = 6.2 A start) following Evang et al. 2022; PBE-D3 full relaxation (ISIF=3)', 144, 16,
           r(a_eff['D_seed1729'] / np.sqrt(2), 4), r(H['D_seed1729']['V'], 2), r(H['D_seed1729']['E'], 5),
           r(H['D_seed1729']['E'] - HL['E'], 4), r(H['D_seed1729']['hull'], 4), '', 'vacancies random on the cation sublattice; a_hex equivalent from volume']])
    save(fig, out, 'Fig1_structures')
    return {'a_eff_O': r(a_eff['O_ordered'], 4), 'a_eff_D': r(a_eff['D_seed1729'], 4), 'gap_layered': r(gL[1] - gL[0], 3),
            'vac_O': r(gO[1] - gO[0], 3), 'E': {k: {kk: r(vv, 4) for kk, vv in v.items()} for k, v in H.items()}}


# ------------------------------------------------------------------ Fig. 2
def fig2(out):
    K = kpath(); gaps = []
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.0), sharey=True, gridspec_kw=dict(hspace=0.34, wspace=0.06))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        x, _ = rf.path_x(m); nel = NEL[m]; g = {}; gt = {}; sh = {}; gam = []; X0 = x[KEND]; tpx = None; xv = {}
        for side, nm, col, ls, lw, z in (('c2_pbe', 'PBE', C_PBE, DASH, 0.9, 2), ('c2_hse', 'HSE', C_HSE, '-', 1.1, 3)):
            E, g[nm] = bands(m, side); bs = range(nel - 10, nel + 8)
            gam += [r(1000 * E[0, nel - 1], 1), r(1000 * E[0, nel], 1), r(1000 * (E[0, nel] - E[0, nel - 1]), 1)]   # Γ 处最高价带、最低导带、两者之差（相对路径价带顶）
            # Sn 系接上附加路径 Γ–Σ′–R₂′–L（c7，两种泛函同一组 k 点）；零点改为两段路径上的最高价带：
            # HSE 取 c5/c6 拟合的真实价带顶（与附加路径上的最高点相差 ≤ 0.1 meV），PBE 取附加路径上的最高点
            tp = true_path_bands(m, side)
            if tp: tp['x'] = TP_XSCALE * tp['x']; tpx = tp
            if side == 'c2_hse':
                T5 = mv.true_load(m); dv = (T5['vbm'] - T5['path_vbm']) if T5 else 0.0
            else:
                dv = max(0.0, float(np.nanmax(tp['E'][:, nel - 1]))) if tp else 0.0
                if tp: tp['E'] = tp['E'] - dv
            E = E - dv; sh[nm] = dv
            gt[nm] = min(np.nanmin(E[:, nel]), np.nanmin(tp['E'][:, nel]) if tp else np.inf)
            for b in bs:
                ax.plot(x[:KEND + 1], E[:KEND + 1, b], color=col, ls=ls, lw=lw, zorder=z)
                if tp: ax.plot(X0 + tp['x'], tp['E'][:, b], color=col, ls=ls, lw=lw, zorder=z)
            # ▼（颜色同线）：该泛函的价带顶在哪里。与路径上的最高点相差 0.5 meV 以内算同一高度；
            # 价带顶在 Γ 或 A 时（路径首尾、附加路径起点都有 Γ），标在中间的 Γ–A 段上
            vb = E[:KEND + 1, nel - 1]; top = np.nanmax(vb)
            if tp and np.nanmax(tp['E'][:, nel - 1]) > top + 5e-4:
                xv[nm] = X0 + tp['x'][int(np.nanargmax(tp['E'][:, nel - 1]))]
            else:
                cand = [i for i in range(KEND + 1) if vb[i] > top - 5e-4]
                xv[nm] = x[next((i for i in cand if i in AXIS), cand[0])]
            dump('Fig2_bands_PBE_vs_HSE', f'{m}_{nm.replace("HSE", "HSE06")}+SOC.csv',
                 ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label'] + [bandname(b, nel) + '_E-EVBM_eV' for b in bs],
                 [[i] + [r(v, 6) for v in K[i]] + [r(x[i], 5), LAB.get(i, '')] + [r(E[i, b], 5) for b in bs] for i in range(KEND + 1)])
            if tp:
                dump('Fig2_bands_PBE_vs_HSE', f'{m}_{nm.replace("HSE", "HSE06")}+SOC_true_valley_path_G-S-R2-L.csv',
                     ['i', 'k_distance_1/A', 'label'] + [bandname(b, nel) + '_E-EVBM_eV' for b in bs],
                     [[i, r(tp['x'][i] / TP_XSCALE, 5), tp['labels'].get(i, '')] + [r(tp['E'][i, b], 5) for b in bs] for i in range(len(tp['x']))])
        gaps.append([m, r(1000 * g['PBE'], 1), r(1000 * g['HSE'], 1)] + gam
                    + [r(1000 * gt['PBE'], 1), r(1000 * gt['HSE'], 1), r(1000 * sh['PBE'], 1), r(1000 * sh['HSE'], 1)])
        ax.axhline(0, color='k', lw=0.5, ls=':', zorder=1)
        kticks(ax, x); ax.set_ylim(-0.66, 0.66); ax.set_yticks(np.arange(-0.6, 0.61, 0.2))
        if tpx: tp_decorate(ax, tpx, X0, nel, mark=())
        stack = abs(xv['PBE'] - xv['HSE']) < 0.012 * (ax.get_xlim()[1] - ax.get_xlim()[0])   # 两个 ▼ 重叠时上下叠放，PBE 在上
        for nm, col in (('HSE', C_HSE), ('PBE', C_PBE)):
            ax.plot(xv[nm], 0.085 if (stack and nm == 'PBE') else 0.04, 'v', ms=3.6, color=col, mec='white', mew=0.3, zorder=5)
        ax.set_title(f'{TEX[m]}     $E_g$: PBE {1000 * gt["PBE"]:.0f} meV, HSE {1000 * gt["HSE"]:.0f} meV', fontsize=8, pad=4)
        panel(ax, f'({"abcd"[n]})')
    dump('Fig2_bands_PBE_vs_HSE', 'band_gaps_on_path.csv',
         ['material', 'Eg_PBE+SOC_meV', 'Eg_HSE06+SOC_meV'] +
         [f'{s}_Gamma_{q}_meV' for s in ('PBE+SOC', 'HSE06+SOC') for q in ('VB', 'CB', 'gap')] +
         ['Eg_PBE+SOC_incl_true_path_meV', 'Eg_HSE06+SOC_incl_true_path_meV', 'VBM_above_path_PBE+SOC_meV', 'VBM_above_path_HSE06+SOC_meV'], gaps)
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    h = [plt.Line2D([], [], color=C_HSE, lw=1.1), plt.Line2D([], [], color=C_PBE, lw=0.9, ls=DASH)]
    fig.legend(h, ['HSE06+SOC', 'PBE+SOC'], loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.905))
    save(fig, out, 'Fig2_bands_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 3
TP_XSCALE = 1.8      # 图 3 中附加路径的横向放大倍数（它只有 0.85 Å⁻¹，原尺度下刻度挤在一起）


def true_path_bands(m, side='c2_hse'):
    """c7_true_path：穿过镜面内真实价带极大的附加路径 Γ–Σ′–R₂′–L（HSE06+SOC，单独计算，分批）。
    side='c2_pbe'：同一组 k 点的 PBE+SOC（c7_true_path_pbe，一次算完），能量相对 PBE 高对称路径上的价带顶。
    返回 dict(x[Å⁻¹，从 0 起], E[k, band]（HSE 相对真实价带顶，eV）, fs[k, band]（阳离子 s）, labels{序号: 名称})；没有结果时返回 None。"""
    import c4_mirror_scan as c4
    pbe = side == 'c2_pbe'
    d = os.path.join(ROOT, '22_能带基准_1-0', m, 'c7_true_path' + ('_pbe' if pbe else ''))
    f = os.path.join(d, 'true_path.json')
    if not os.path.exists(f): return None
    meta = json.load(open(f, encoding='utf-8'))
    if pbe:
        if not any(os.path.exists(os.path.join(d, q)) for q in ('PROCAR_OPT.gz', 'PROCAR_OPT')): return None
        parts, ref = None, np.nanmax(cv.read_path_eigen(os.path.join(rf.base(m), 'c2_pbe', 'path_eigen.txt'))[:, NEL[m] - 1])
    else:
        if not os.path.isdir(os.path.join(d, 'b01')): return None
        parts, ref = [os.path.join(d, f'b{i + 1:02d}') for i in range(len(meta['batches']))], meta['vbm']
    E, fs = c4.read_parts(d, parts, meta)
    W = np.concatenate([cv.read_procar(c4.find_procar(q))[1] for q in (parts or [d])], axis=0)   # [k, band, ion, orb]，图 4、图 5 用
    k = np.array(meta['path_k_index']); Kq = np.array(c4.read_kpoints_opt(os.path.join(d, 'KPOINTS_OPT')))
    return dict(x=np.array(meta['x_invA']), E=E[k] - ref, fs=fs[k], W=W[k], k=Kq[k], labels={int(i): v for i, v in meta['labels'].items()})


def tp_decorate(ax, tp, X0, nel, mark=('Σ′', 'R₂′')):
    """把附加路径 Γ–Σ′–R₂′–L 接在主路径右边之后的装饰：浅色底、粗竖线断开、刻度 A|Γ、Σ′、R₂′、L，真实极大标 ▼ 加名称。
    先画好主路径并调用 kticks、set_ylim，再调用本函数。"""
    xe = X0 + tp['x'][-1]
    ax.axvspan(X0, xe, color='#F1EEF7', zorder=0, lw=0); ax.axvline(X0, color='k', lw=0.8, zorder=1)
    tk, tl = list(ax.get_xticks()), [t.get_text() for t in ax.get_xticklabels()]
    tl[-1] = tl[-1] + '|Γ'
    for i, nm in sorted(tp['labels'].items()):
        if i == 0: continue
        tk.append(X0 + tp['x'][i]); tl.append(MATH.get(nm, nm))
        if nm != 'L': ax.axvline(X0 + tp['x'][i], color=GRID, lw=0.5, zorder=0)
        if nm in mark: vmark(ax, X0 + tp['x'][i], tp['E'][i, nel - 1], nm)
    ax.set_xticks(tk, tl); ax.set_xlim(ax.get_xlim()[0], xe)
    ax.text((X0 + xe) / 2, ax.get_ylim()[0] + 0.015, f'off path\n(k ×{TP_XSCALE:g})', ha='center', va='bottom', fontsize=5.6,
            linespacing=1.0, color='#444444', bbox=dict(fc='#F1EEF7', ec='none', pad=0.5))


def fig3(out):
    cmap = LinearSegmentedColormap.from_list('s', ['#DCE8F4', '#9CC3E4', '#4A90C8', '#1F5FA0', '#08306B'])
    norm = Normalize(0, 0.30); K = kpath()
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.0), sharey=True, gridspec_kw=dict(hspace=0.34, wspace=0.06))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        x, _ = rf.path_x(m); nel = NEL[m]
        E, gap = bands(m, 'c2_hse')
        T5 = mv.true_load(m); dv = (T5['vbm'] - T5['path_vbm']) if T5 else 0.0
        E, gap = E - dv, gap - dv                                # 能量零点取真实价带顶（Sn₂Bi₂Te₅ 比路径上的最高点高 38 meV）
        tp = true_path_bands(m); X0 = x[KEND]                    # 附加路径 Γ–Σ′–R₂′–L 接在主路径右边，横向放大 TP_XSCALE 倍
        if tp: tp['x'] = TP_XSCALE * tp['x']
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fs = rf.cation_s(Wp); rows = []
        for b in range(nel - 8, nel + 6):
            segs, cols = [], []
            for i in range(KEND):
                segs.append([(x[i], E[i, b]), (x[i + 1], E[i + 1, b])]); cols.append((fs[i, b] + fs[i + 1, b]) / 2)
            cols = np.array(cols)
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=0.5 + 2.4 * np.clip(cols / 0.3, 0, 1), capstyle='round')
            lc.set_array(cols); ax.add_collection(lc)
            rows += [[i] + [r(v, 6) for v in K[i]] + [r(x[i], 5), LAB.get(i, ''), bandname(b, nel), r(E[i, b], 5), r(100 * fs[i, b], 2)]
                     for i in range(KEND + 1)]
            if tp:
                xs, es, fc = X0 + tp['x'], tp['E'][:, b], tp['fs'][:, b]
                cc = (fc[:-1] + fc[1:]) / 2
                lc = LineCollection([[(xs[i], es[i]), (xs[i + 1], es[i + 1])] for i in range(len(xs) - 1)], cmap=cmap, norm=norm,
                                    linewidths=0.5 + 2.4 * np.clip(cc / 0.3, 0, 1), capstyle='round')
                lc.set_array(cc); ax.add_collection(lc)
        dump('Fig3_HSE_bands_cation_s', f'{m}_HSE06+SOC_cation_s.csv',
             ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV', 'cation_s_percent'], rows)
        ax.axhline(0, color='k', lw=0.5, ls=':'); ax.axhline(gap, color='k', lw=0.5, ls=':')
        for nm, _E, kind, i in levels(m, true=False):        # 标出路径上的各谷（名称与图 5、图 6 相同）
            if kind == 'ref': continue
            # Γ–A 段很窄：Δ，以及带隙很小时的 Γ（▼ 会压到导带上）改用引线，名称放在左上方
            lead = (x[90] - 0.28, 0.10) if nm == 'Δ' or (nm == 'Γ' and gap < 0.1) else None
            vmark(ax, x[i], E[i, nel - 1], nm, lead)
        kticks(ax, x); ax.set_ylim(-0.5, 0.5)
        if tp:
            tp_decorate(ax, tp, X0, nel)
            dump('Fig3_HSE_bands_cation_s', f'{m}_true_valley_path_G-S-R2-L.csv',
                 ['i', 'k_distance_1/A', 'label', 'band', 'E-EVBM_true_eV', 'cation_s_percent'],
                 [[i, r(tp['x'][i] / TP_XSCALE, 5), tp['labels'].get(i, ''), bandname(b, nel), r(tp['E'][i, b], 5), r(100 * tp['fs'][i, b], 2)]
                  for b in range(nel - 8, nel + 6) for i in range(len(tp['x']))])
        ax.set_title(f'{TEX[m]}     $E_g$ = {1000 * gap:.0f} meV' + (f' (path only: {1000 * (gap + dv):.0f} meV)' if dv > 5e-4 else ''),
                     fontsize=8, pad=4)
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.7, pad=0.015, aspect=30, extend='max')
    cb.set_label('Cation-s weight (%)'); cb.set_ticks([0, .1, .2, .3], labels=['0', '10', '20', '30'])
    cb.ax.tick_params(direction='in', length=2, width=0.6); cb.outline.set_linewidth(0.6)
    save(fig, out, 'Fig3_HSE_bands_cation_s')


# ------------------------------------------------------------------ Fig. 4
SITE_CMAP = LinearSegmentedColormap.from_list('ab', ['#E66100', '#F2B27A', '#E3E3E3', '#B7A3D6', '#5D3A9B'])   # B 位（橙）→ A 位（紫）
GREY = np.array([0.80, 0.80, 0.80, 1.0])


def s_by_atom(W):
    """W[k, band, ion, orb]，离子顺序 B×2、Te×5、A×2（四种材料的 POSCAR 都是 B Te A）。返回 A、B、Te 的 s 占球内总投影的比例。"""
    tot = np.maximum(W.sum((-1, -2)), 1e-9)
    return W[:, :, [7, 8], 0].sum(-1) / tot, W[:, :, [0, 1], 0].sum(-1) / tot, W[:, :, 2:7, 0].sum(-1) / tot


def fig4(out):
    norm = Normalize(0, 1); K = kpath(); key = []
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.0), sharey=True, gridspec_kw=dict(hspace=0.34, wspace=0.06))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        x, _ = rf.path_x(m); nel = NEL[m]; A_, B_ = m[:2], m[3:5]
        E, gap = bands(m, 'c2_hse')
        T5 = mv.true_load(m); dv = (T5['vbm'] - T5['path_vbm']) if T5 else 0.0
        E, gap = E - dv, gap - dv                                # 能量零点取真实价带顶（与图 3 相同）
        tp = true_path_bands(m); X0 = x[KEND]
        if tp: tp['x'] = TP_XSCALE * tp['x']; ta, tb, _tt = s_by_atom(tp['W'])
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fa, fb, ft = s_by_atom(Wp); rows = []
        for b in range(nel - 8, nel + 6):
            segs, cols, wid = [], [], []
            for i in range(KEND):
                a2, b2 = (fa[i, b] + fa[i + 1, b]) / 2, (fb[i, b] + fb[i + 1, b]) / 2
                tot = a2 + b2; w = np.clip(tot / 0.15, 0, 1)      # 阳离子 s 少于 15% 时逐渐变灰
                segs.append([(x[i], E[i, b]), (x[i + 1], E[i + 1, b])])
                cols.append(w * np.array(SITE_CMAP(norm(a2 / max(tot, 1e-9)))) + (1 - w) * GREY)
                wid.append(0.5 + 2.4 * np.clip(tot / 0.3, 0, 1))
            if tp:                                               # 附加路径 Γ–Σ′–R₂′–L，同样着色
                xs, es = X0 + tp['x'], tp['E'][:, b]
                for i in range(len(xs) - 1):
                    a2, b2 = (ta[i, b] + ta[i + 1, b]) / 2, (tb[i, b] + tb[i + 1, b]) / 2
                    tot = a2 + b2; w = np.clip(tot / 0.15, 0, 1)
                    segs.append([(xs[i], es[i]), (xs[i + 1], es[i + 1])])
                    cols.append(w * np.array(SITE_CMAP(norm(a2 / max(tot, 1e-9)))) + (1 - w) * GREY)
                    wid.append(0.5 + 2.4 * np.clip(tot / 0.3, 0, 1))
            ax.add_collection(LineCollection(segs, colors=cols, linewidths=wid, capstyle='round'))
            rows += [[i] + [r(v, 6) for v in K[i]] + [r(x[i], 5), LAB.get(i, ''), bandname(b, nel), r(E[i, b], 5),
                     r(100 * fa[i, b], 2), r(100 * fb[i, b], 2), r(100 * ft[i, b], 2)] for i in range(KEND + 1)]
        dump('Fig4_cation_s_by_atom', f'{m}_HSE06+SOC_s_by_atom.csv',
             ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV',
              f'A_{A_}_s_percent', f'B_{B_}_s_percent', 'Te_s_percent'], rows)
        # 关键态：价带顶、导带底、最高孤对电子价带谷（Kramers 对平均）
        iv = int(np.nanargmax(E[:KEND + 1, nel - 1])); ic = int(np.nanargmin(E[:KEND + 1, nel])); lp = rf.lone_pair_offset(m, 'c2_hse')
        vb, cb = [nel - 2, nel - 1], [nel, nel + 1]
        P = lambda i, bb: [i] + [r(v, 4) for v in K[i]] + [r(1000 * E[i, bb].mean(), 1), A_, r(100 * fa[i, bb].mean(), 1),
                                                          B_, r(100 * fb[i, bb].mean(), 1), r(100 * ft[i, bb].mean(), 1)]
        Q = lambda i, bb: [f'c7:{i} {tp["labels"].get(i, "")}'.strip()] + [r(v, 4) for v in tp['k'][i]] + [
            r(1000 * tp['E'][i, bb].mean(), 1), A_, r(100 * ta[i, bb].mean(), 1), B_, r(100 * tb[i, bb].mean(), 1), r(100 * _tt[i, bb].mean(), 1)]
        kv, kl = P(iv, vb), P(lp[1], vb)
        if tp:                                                   # Sn 系：价带顶、孤对电子谷取附加路径上的真实极大（与图 7 相同）
            jv = int(np.nanargmax(tp['E'][:, nel - 1]))
            if tp['E'][jv, nel - 1] > E[iv, nel - 1] + 5e-4: kv = Q(jv, vb); kv[4] = 0.0           # 它就是真实价带顶
            tv = true_valleys(m); t0 = tv['valleys'][0]; jl = next((i for i, nm in tp['labels'].items() if nm == t0['name']), None)
            if jl is not None: kl = Q(jl, vb); kl[4] = r(t0['E'], 1)       # 成分取附加路径上的这一点，能量取三维搜索值
        for st, row in (('VBM', kv), ('lone-pair valence valley', kl), ('CBM', P(ic, cb))):
            key.append([m, st] + row)
        ax.axhline(0, color='k', lw=0.5, ls=':'); ax.axhline(gap, color='k', lw=0.5, ls=':')
        kticks(ax, x); ax.set_ylim(-0.5, 0.5)
        if tp: tp_decorate(ax, tp, X0, nel, mark=())            # 本图主路径也不标谷
        ax.set_title(f'{TEX[m]}     A = {A_}, B = {B_}', fontsize=8, pad=4)
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SITE_CMAP), ax=axs, shrink=0.7, pad=0.015, aspect=30)
    cb.set_label('Share of cation-s on the A site'); cb.set_ticks([0, 0.5, 1], labels=['B only\n(Sb, Bi)', 'equal', 'A only\n(Pb, Sn)'])
    cb.ax.tick_params(direction='in', length=2, width=0.6, labelsize=6.5); cb.outline.set_linewidth(0.6)
    dump('Fig4_cation_s_by_atom', 'key_states.csv', ['material', 'state', 'k_index (c7: true-valley path)', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'E-EVBM_meV',
                                                     'A', 'A_s_percent', 'B', 'B_s_percent', 'Te_s_percent'], key)
    save(fig, out, 'Fig4_cation_s_by_atom')
    return {f'{k[0]}|{k[1]}': (k[8], k[10]) for k in key}


# ------------------------------------------------------------------ Fig. 5
C_A, C_B = '#5D3A9B', '#E66100'          # A 位（Sn）紫、B 位（Sb/Bi）橙，与图 4 的两端一致
SNVAL = {'Γ': 0, 'Σ': 11, 'R₁': 101, 'R₂': 119}   # 价带：Γ 单谷与三组孤对电子谷（路径序号）


def fig5(out):
    """Sn₂Sb₂Te₅ 与 Sn₂Bi₂Te₅ 对比：上排 Sn 5s 投影（孤对电子谷），下排 B 位 p 投影（Γ 单谷），右侧为两者的能级。"""
    mats = ['Sn2Sb2Te5', 'Sn2Bi2Te5']
    fig = plt.figure(figsize=(W2, 4.7))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.55], hspace=0.36, wspace=0.12)
    lv = {}; gam = []
    for c, m in enumerate(mats):
        x, _ = rf.path_x(m); nel = NEL[m]; B_ = m[3:5]
        E, gap = bands(m, 'c2_hse')                              # E：相对路径价带顶，表 5 与 valence_valleys.csv 用它
        T5 = mv.true_load(m); dv = (T5['vbm'] - T5['path_vbm']) if T5 else 0.0
        Ep = E - dv                                              # 画图：能量零点取真实价带顶（与图 3 相同）
        tp = true_path_bands(m); X0 = x[KEND]
        if tp:
            tp['x'] = TP_XSCALE * tp['x']; tW = tp['W']; ttot = np.maximum(tW.sum((-1, -2)), 1e-9)
            tproj = (tW[:, :, [7, 8], 0].sum(-1) / ttot, tW[:, :, [0, 1], 1:4].sum((-1, -2)) / ttot)   # Sn 5s、B 位 p
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        tot = np.maximum(Wp.sum((-1, -2)), 1e-9)
        sn_s = Wp[:, :, [7, 8], 0].sum(-1) / tot; b_s = Wp[:, :, [0, 1], 0].sum(-1) / tot
        b_p = Wp[:, :, [0, 1], 1:4].sum((-1, -2)) / tot; sn_p = Wp[:, :, [7, 8], 1:4].sum((-1, -2)) / tot
        te_p = Wp[:, :, 2:7, 1:4].sum((-1, -2)) / tot
        rows = []
        for row, (w, wmax, col, lab) in enumerate(((sn_s, 0.25, C_A, 'Sn 5s'), (b_p, 0.45, C_B, f'{B_} {"5p" if B_ == "Sb" else "6p"}'))):
            ax = fig.add_subplot(gs[row, c])
            for b in range(nel - 8, nel + 6):
                ax.plot(x[:KEND + 1], Ep[:KEND + 1, b], color='#C8C8C8', lw=0.5, zorder=1)
                segs = [[(x[i], Ep[i, b]), (x[i + 1], Ep[i + 1, b])] for i in range(KEND)]
                wid = [3.2 * np.clip((w[i, b] + w[i + 1, b]) / 2 / wmax, 0, 1) for i in range(KEND)]
                if tp:                                           # 附加路径 Γ–Σ′–R₂′–L，同一投影
                    xs, es, tw = X0 + tp['x'], tp['E'][:, b], tproj[row]
                    ax.plot(xs, es, color='#C8C8C8', lw=0.5, zorder=1)
                    segs += [[(xs[i], es[i]), (xs[i + 1], es[i + 1])] for i in range(len(xs) - 1)]
                    wid += [3.2 * np.clip((tw[i, b] + tw[i + 1, b]) / 2 / wmax, 0, 1) for i in range(len(xs) - 1)]
                ax.add_collection(LineCollection(segs, colors=col, linewidths=wid, capstyle='round', zorder=2))
            ax.axhline(0, color='k', lw=0.5, ls=':')
            kticks(ax, x); ax.set_ylim(-0.42, 0.30)
            if tp: tp_decorate(ax, tp, X0, nel, mark=('Σ′', 'R₂′') if row == 0 else ())
            if c == 0: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
            else: ax.tick_params(labelleft=False)
            ax.set_title(f'{TEX[m]}: {lab}', fontsize=8, pad=4, color=col)
            panel(ax, f'({"abcd"[2 * row + c]})')
            for nm, i in SNVAL.items():                      # 标出四个谷（名称与图 3、图 6 相同）：上排标孤对电子谷，下排标 Γ 单谷
                if (row == 0) == (nm != 'Γ'):
                    vmark(ax, x[90] if nm == 'Γ' else x[i], Ep[i, nel - 1], nm)
        for nm, i in SNVAL.items():
            bb = [nel - 2, nel - 1]; g = lambda X: r(100 * X[i, bb].mean(), 1)
            lv[(m, nm)] = 1000 * E[i, bb].mean()
            rows.append([nm, i] + [r(v, 4) for v in kpath()[i]] + [r(1000 * E[i, bb].mean(), 1), g(sn_s), g(b_s), g(sn_p), g(b_p), g(te_p)])
        for st, bb in (('Gamma_VB', [nel - 2, nel - 1]), ('Gamma_CB', [nel, nel + 1])):   # Γ 处最高价带、最低导带（表 5 前两行）
            gam.append([m, st, r(1000 * E[0, bb].mean(), 1)] + [r(100 * X[0, bb].mean(), 1) for X in (sn_s, b_s, sn_p, b_p, te_p)])
        dump('Fig5_SnSb_vs_SnBi', f'{m}_valence_valleys.csv',
             ['valley', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'E-EVBM_meV', 'Sn_s_percent', f'{B_}_s_percent',
              'Sn_p_percent', f'{B_}_p_percent', 'Te_p_percent'], rows)
        dump('Fig5_SnSb_vs_SnBi', f'{m}_HSE06+SOC_projections.csv',
             ['k_index', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV', 'Sn_s_percent', f'{B_}_s_percent', 'Sn_p_percent', f'{B_}_p_percent', 'Te_p_percent'],
             [[i, r(x[i], 5), LAB.get(i, ''), bandname(b, nel), r(E[i, b], 5), r(100 * sn_s[i, b], 2), r(100 * b_s[i, b], 2),
               r(100 * sn_p[i, b], 2), r(100 * b_p[i, b], 2), r(100 * te_p[i, b], 2)] for b in range(nel - 8, nel + 6) for i in range(KEND + 1)])
    # (e) 能级：各自以真实价带顶为零点。Γ 单谷取 Γ 点本身（路径值减价带顶的移动），孤对电子谷取三维搜索的真实极大
    ax = fig.add_subplot(gs[:, 2]); X = (0.0, 1.7); hw = 0.3
    lt = {}
    for m in mats:
        tv = true_valleys(m); sh = tv['shift'] if tv else 0.0
        lt[(m, 'Γ')] = lv[(m, 'Γ')] - sh
        for t in (tv['valleys'] if tv else []): lt[(m, t['name'])] = t['E']
    names = [n for n in ('Σ′', 'R₁′', 'R₂′', 'R₁″') if all((m, n) in lt for m in mats)]
    lpt = [max(lt[(m, n)] for n in names) for m in mats]; lpb = [min(lt[(m, n)] for n in names) for m in mats]
    for k in range(2):                               # Sn 5s 孤对电子谷：色块标出范围，细线为各个谷
        ax.add_patch(plt.Rectangle((X[k] - hw, lpb[k]), 2 * hw, lpt[k] - lpb[k], color=C_A, alpha=0.13, lw=0, zorder=1))
        for n in names: ax.plot([X[k] - hw, X[k] + hw], [lt[(mats[k], n)]] * 2, color=C_A, lw=1.0, solid_capstyle='butt', zorder=2)
        ax.plot([X[k] - hw, X[k] + hw], [lt[(mats[k], 'Γ')]] * 2, color=C_B, lw=2.2, solid_capstyle='butt', zorder=3)
    for y, col in ((lpt, C_A), ([lt[(m, 'Γ')] for m in mats], C_B)):   # 只连两组的最高点与 Γ：一升一降
        ax.plot([X[0] + hw, X[1] - hw], y, color=col, lw=0.8, ls=(0, (2, 1.5)), zorder=2)
    kw = dict(fontsize=6.5, ha='center')
    ax.text(X[0], lt[(mats[0], 'Γ')] + 3, 'Γ (Sb p)', va='bottom', color=C_B, **kw)
    ax.text(X[0], lpb[0] - 3, 'Sn 5s\n(4 valleys)', va='top', color=C_A, linespacing=1.1, **kw)
    ax.text(X[1], lt[(mats[1], 'Γ')] - 3, 'Γ (Bi p)', va='top', color=C_B, **kw)
    ax.text(X[1], lpb[1] - 3, 'Sn 5s\n(4 valleys)', va='top', color=C_A, linespacing=1.1, **kw)
    dump('Fig5_SnSb_vs_SnBi', 'levels_true.csv', ['state', 'E_Sn2Sb2Te5_meV', 'E_Sn2Bi2Te5_meV'],
         [[nm, r(lt[('Sn2Sb2Te5', nm)], 1), r(lt[('Sn2Bi2Te5', nm)], 1)] for nm in ['Γ'] + names])
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xlim(-0.6, 2.3); ax.set_ylim(-160, 14)
    ax.set_xticks(X, [TEX[m] for m in mats], fontsize=6.5); ax.tick_params(axis='x', top=False, length=0)
    ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (meV)'); ax.yaxis.set_label_position('right'); ax.yaxis.tick_right()
    ax.tick_params(axis='y', left=True, labelleft=False, right=True, labelright=True)
    ax.set_title('Valence-band levels', fontsize=8, pad=4)
    panel(ax, '(e)', x=-0.04)
    dump('Fig5_SnSb_vs_SnBi', 'Gamma_states.csv',
         ['material', 'state', 'E-EVBM_meV', 'Sn_s_percent', 'B_s_percent', 'Sn_p_percent', 'B_p_percent', 'Te_p_percent'], gam)
    dump('Fig5_SnSb_vs_SnBi', 'levels.csv', ['valley', 'E_Sn2Sb2Te5_meV', 'E_Sn2Bi2Te5_meV'],
         [[nm, r(lv[('Sn2Sb2Te5', nm)], 1), r(lv[('Sn2Bi2Te5', nm)], 1)] for nm in SNVAL])
    save(fig, out, 'Fig5_SnSb_vs_SnBi_orbitals')
    return {f'{a}|{b}': r(v, 1) for (a, b), v in lv.items()}


# ------------------------------------------------------------------ Fig. 6
KT300 = 1000 * mv.KB * 300                      # 300 K 的 k_BT（meV）
AXIS = range(90, 98)                            # 路径上的 Γ–A 段（Δ 线）


def bz_lattice(m):
    """所给层状结构的倒格矢（行为 b1, b2, b3，Å⁻¹，含 2π）与 12 个点群操作（D3d，spglib）。"""
    lat, frac, sp = read_poscar(os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'POSCAR'))
    num = [sorted(set(sp)).index(e) for e in sp]
    rots = spglib.get_symmetry((lat, frac, num), symprec=1e-3)['rotations']
    rots = np.unique(rots.reshape(len(rots), -1), axis=0).reshape(-1, 3, 3)
    assert len(rots) == 12, len(rots)
    return 2 * np.pi * np.linalg.inv(lat).T, rots


def star(k, rots):
    """k 的全部等价点（分数坐标，模倒格矢去重；含时间反演，有反演时不增加新点）。"""
    out = []
    for R in rots:
        for sg in (1, -1):
            q = sg * (np.asarray(k, float) @ R); q -= np.round(q); q[np.isclose(q, -0.5)] = 0.5
            if not any(np.allclose(q, p, atol=1e-6) for p in out): out.append(q)
    return out


def ws(q, B):
    """分数坐标 q 换到第一布里渊区（面内取最短的等价像），返回笛卡尔坐标（Å⁻¹）。"""
    c = [(q + np.array([i, j, 0])) @ B for i in (-1, 0, 1) for j in (-1, 0, 1)]
    return min(c, key=lambda v: np.hypot(v[0], v[1]))


def top_vb(m, side):
    """路径上最高价带（Kramers 对取平均，与图 5、表 5 一致），相对价带顶，meV。"""
    E, _ = bands(m, side); nel = NEL[m]
    return 1000 * E[:, nel - 2:nel].mean(1)


PATH_LP = {'Σ': (-1 / 6, 0.0), 'R₁': (-2 / 33, 0.5), 'R₂': (-1 / 3, 0.5)}     # 路径上的三组孤对电子谷 (u, w)


def true_valleys(m):
    """镜面内的真实价带谷（c4_mirror_v3 + c5_star；目前只有 Sn₂Bi₂Te₅）。按离哪个路径谷最近命名并加撇：Σ′、R₂′；
    同一个路径谷演化出两个极大时，能量高的为 ′、低的为 ″（R₁′、R₁″）。能量相对真实价带顶（meV）。没有结果时返回 None。"""
    T = mv.true_load(m)
    if T is None: return None
    B, _ = bz_lattice(m); nb1, cz = np.linalg.norm(B[0]), abs(B[2, 2])
    def dist(t, uw):                                            # 笛卡尔距离（Å⁻¹）；分数坐标下沿 c 的距离被放大约 4.6 倍
        dw = abs(t['w'] - uw[1]) % 1; dw = min(dw, 1 - dw)
        return np.hypot((t['u'] - uw[0]) * nb1, dw * cz)
    groups = {}
    for t in T['valleys']:
        groups.setdefault(min(PATH_LP, key=lambda n: dist(t, PATH_LP[n])), []).append(t)
    out = []
    for n, ts in groups.items():
        for k, t in enumerate(sorted(ts, key=lambda t: -t['E'])):
            out.append(dict(name=n + '′″'[k], kfrac=(t['u'], 0.0, t['w']), E=1000 * (t['E'] - T['vbm']), s=t['s'],
                            grid='c4 网格' in t['src'], m=t['m'], principal=t['principal']))
    return dict(shift=1000 * (T['vbm'] - T['path_vbm']), valleys=sorted(out, key=lambda d: -d['E']))


def crossover():
    """价带谷间能量差 ΔE（HSE06+SOC）：最高孤对电子谷（与图 7 同一判据）减 Γ–A 轴（Γ、Δ 线、A）上最高的价带态，四种材料。"""
    K = kpath(); rows = []
    for m in MATS:
        _, rots = bz_lattice(m)
        for side, fn in (('c2_hse', 'HSE06+SOC'),):
            v = top_vb(m, side); ia = max(AXIS, key=lambda i: v[i])
            lp = rf.lone_pair_offset(m, side); il = lp[1]
            tv = true_valleys(m)
            if tv:                                                  # 有镜面内的真实谷：价带顶与孤对电子谷取真实值
                v = v - tv['shift']; t0 = tv['valleys'][0]
                rows.append(dict(m=m, fn=fn, EG=v[0], Ea=v[ia], ia=ia, na=len(star(K[ia], rots)), El=t0['E'], il=None, lname=t0['name'],
                                 kl=t0['kfrac'], nl=len(star(t0['kfrac'], rots)), s=t0['s'], d=t0['E'] - v[ia]))
                continue
            rows.append(dict(m=m, fn=fn, EG=v[0], Ea=v[ia], ia=ia, na=len(star(K[ia], rots)),
                             El=v[il], il=il, lname=kname(il), kl=tuple(K[il]), nl=len(star(K[il], rots)), s=lp[2], d=v[il] - v[ia]))
    return rows


def pockets(m):
    """HSE06+SOC：Γ–A 轴上最高的价带态与三组孤对电子谷（Σ、R₁、R₂），展开成第一布里渊区内的全部等价口袋。"""
    B, rots = bz_lattice(m); K = kpath(); v = top_vb(m, 'c2_hse')
    ia = max(AXIS, key=lambda i: v[i]); tv = true_valleys(m)
    if tv:                                                      # 镜面内的真实谷（价带顶取真实值；轴上的态仍取路径，它在三重轴上）
        pts = star(K[ia], rots); out = [dict(name=kname(ia), i=ia, E=v[ia] - tv['shift'], col=C_B, frac=pts, cart=[ws(q, B) for q in pts], grid=False)]
        for t in tv['valleys']:
            pts = star(t['kfrac'], rots)
            out.append(dict(name=t['name'], i=None, E=t['E'], col=C_A, frac=pts, cart=[ws(q, B) for q in pts], grid=t['grid']))
        return B, out
    out = []
    for i, col in [(ia, C_B)] + [(i, C_A) for nm, i in SNVAL.items() if nm != 'Γ']:
        pts = star(K[i], rots)
        out.append(dict(name=kname(i), i=i, E=v[i], col=col, frac=pts, cart=[ws(q, B) for q in pts], grid=False))
    return B, out


def hexagon(B):
    K0 = np.array([-2 / 3, 1 / 3, 0]) @ B; a = np.arctan2(K0[1], K0[0]) + np.arange(7) * np.pi / 3
    return np.hypot(*K0[:2]) * np.cos(a), np.hypot(*K0[:2]) * np.sin(a)


MK = {'Γ': 'o', 'Δ': 'o', 'Σ': 'o', 'R₁': 's', 'R₂': 'D', 'Σ′': 'o', 'Σ″': 'h', 'R₁′': 's', 'R₁″': 'p', 'R₂′': 'D'}   # 形状区分谷：Σ 圆、R₁ 方、R₂ 菱形、R₁″ 五边形；实心/空心区分所在面
MS = {'s': 0.55, 'D': 0.8, 'p': 0.6, 'h': 0.65}                                    # 方块、菱形看起来比圆大，略缩小


def near(E):
    return E >= -2 * KT300                      # 价带顶 2k_BT 以内：画大、画实


def on_face(f):
    return np.isclose(abs(f[2]), 0.5)


def bz_top(ax, m, B, P, letter):
    hx, hy = hexagon(B)
    ax.plot(hx, hy, color='k', lw=0.8, zorder=1)
    r0 = 0.5 * np.linalg.norm(B[0][:2]); a0 = np.arctan2(B[0][1], B[0][0])
    M = [(r0 * np.cos(a0 + t * np.pi / 3), r0 * np.sin(a0 + t * np.pi / 3)) for t in range(6)]   # 6 个 M，在边中点
    for q in M: ax.plot([0, q[0]], [0, q[1]], color='#BDBDBD', lw=0.5, ls=DASH, zorder=0)
    for g in P:
        big = near(g['E'])
        for c, f in zip(g['cart'], g['frac']):
            ax.plot(c[0], c[1], MK[g['name']], ms=(5.6 if big else 3.4) * MS.get(MK[g['name']], 1), mew=0.9, color=g['col'], mfc='white' if on_face(f) else g['col'],
                    alpha=1.0 if big else 0.45, zorder=3)
    Mq = max(M, key=lambda c: c[0] + 0.1 * c[1]); kv = max(zip(hx, hy), key=lambda t: (round(t[1], 3), t[0]))
    ax.text(Mq[0] + 0.04, Mq[1], 'M (L)', ha='left', va='center', fontsize=6.5)
    ax.text(kv[0] + 0.03, kv[1] + 0.02, 'K (H)', ha='left', va='bottom', fontsize=6.5)
    ax.annotate('Γ (A)', (0, 0), (-0.62, -0.50), fontsize=6.5, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.5, color='#555555', shrinkA=1, shrinkB=3))
    y = -0.66                                     # 右下角小表：位置名、相对价带顶的能量（meV）、等价口袋数
    for g in P:
        al = 1.0 if near(g['E']) else 0.45
        ax.plot(0.74, y, MK[g['name']], ms=3.6 * (0.9 if MK[g['name']] == 's' else MS.get(MK[g['name']], 1)), mew=0.8, color=g['col'], mfc='white' if all(map(on_face, g['frac'])) else g['col'], alpha=al)
        ax.text(0.83, y, MATH.get(g['name'], g['name']), ha='left', va='center', fontsize=6.3)
        ax.text(1.66, y, f"{g['E']:.0f}".replace('-', '−') + ('*' if g.get('grid') else '') + f"  ×{len(g['frac'])}", ha='right', va='center', fontsize=6.3)
        y -= 0.115
    n = sum(len(g['frac']) for g in P if near(g['E']))
    ax.set_title(f"{TEX[m]}: {n} pocket{'s' if n > 1 else ''} within $2k_\\mathrm{{B}}T$", fontsize=8, pad=2)
    ax.set_aspect('equal'); ax.set_xlim(-1.0, 1.68); ax.set_ylim(-1.19, 0.98); ax.axis('off')
    ax.text(0.0, 1.06, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')


def bz_3d(ax, m, B, P, zs):
    hx, hy = hexagon(B); h = 0.5 * np.linalg.norm(B[2]) * zs
    for z in (-h, h): ax.plot(hx, hy, [z] * 7, color='k', lw=0.6)
    for x, y in zip(hx[:-1], hy[:-1]): ax.plot([x, x], [y, y], [-h, h], color='k', lw=0.6)
    ax.plot([0, 0], [0, 0], [-h, h], color='#BDBDBD', lw=0.5, ls=DASH)
    for g in P:
        big = near(g['E'])
        for c, f in zip(g['cart'], g['frac']):
            zz = (h, -h) if on_face(f) else (c[2] * zs,)          # 区界面上的口袋：顶面、底面是同一个，两处都画
            for z in zz:
                ax.plot([c[0]], [c[1]], [z], MK[g['name']], ms=(3.8 if big else 2.6) * MS.get(MK[g['name']], 1), color=g['col'], mfc='white' if on_face(f) else g['col'],
                        mew=0.8, alpha=1.0 if big else 0.45)
    # 每组谷标一次名称，与 (a)(b) 角上的小表对照：R₂ 标在右前方的口袋旁，R₁、Σ 标在左前方，Δ 标在上方的那个点旁
    for g in P:
        if g['name'] == 'Γ': continue
        cs = [np.array([c[0], c[1], h if on_face(f) else c[2] * zs]) for c, f in zip(g['cart'], g['frac'])]
        if g['name'] == 'R₂': c, dx, ha = max(cs, key=lambda c: (c[0], -c[1])), 0.07, 'left'
        elif g['name'] == 'Δ': c, dx, ha = max(cs, key=lambda c: c[2]), 0.06, 'left'
        else: c, dx, ha = min(cs, key=lambda c: (c[0], c[1])), -0.07, 'right'
        ax.text(c[0] + dx, c[1], c[2], MATH.get(g['name'], g['name']), fontsize=6.5, ha=ha, va='center', color=g['col'],
                alpha=1.0 if near(g['E']) else 0.6)
    # 对称点标签，方向与 (a)(b) 的俯视图一致：Γ、A 在轴上，M/L 在 30° 方向的边中点，K/H 在 60° 方向的顶点
    r0 = 0.5 * np.linalg.norm(B[0][:2]); a0 = np.arctan2(B[0][1], B[0][0]); kv = max(zip(hx, hy), key=lambda t: (round(t[1], 3), t[0]))
    for z in (0, h, -h): ax.plot([0], [0], [z], marker='+', ms=5, mew=0.9, color='k', zorder=5)     # Γ（区中心）与 A（顶、底面中心）
    ax.text(0.03, 0, h + 0.05, 'A', fontsize=6.5, ha='left', va='bottom'); ax.text(-0.07, 0, 0.02, 'Γ', fontsize=7, ha='right', va='bottom')
    for z, (lm, lk) in ((h, ('L', 'H')), (-h, ('M', 'K'))):
        ax.text(1.08 * r0 * np.cos(a0), 1.08 * r0 * np.sin(a0), z, lm, fontsize=6.5, ha='left', va='center')
        ax.text(1.05 * kv[0], 1.05 * kv[1], z, lk, fontsize=6.5, ha='left', va='bottom')
    # 从 −y 方向斜上方看：kx 水平向右，与俯视图 (a)(b) 的方向相同
    ax.view_init(elev=30, azim=-90); ax.set_box_aspect((1, 1, 2 * h / 1.9), zoom=1.3); ax.set_axis_off()


def levels(m, true=True):
    """HSE06+SOC 价带上的能级（相对价带顶，meV）：Γ–A 轴上的谷（Γ、Δ 或 A）、Γ 点本身（不是谷时作参照）、
    孤对电子谷（Sn 系用表 5 的 Σ、R₁、R₂；Pb 系取 Σ、R 线上阳离子 s ≥ 15% 的局部极大，260 meV 以内）。"""
    v = top_vb(m, 'c2_hse'); ia = max(AXIS, key=lambda i: v[i])
    tv = true_valleys(m) if true else None
    if tv:                                                      # 真实谷（相对真实价带顶）；轴上的态与 Γ 点仍取路径值，减去价带顶的移动
        v = v - tv['shift']
        out = [(kname(ia), v[ia], 'axis', ia)] + ([('Γ', v[0], 'ref', 90)] if ia != 90 else [])
        return out + [(t['name'], t['E'], 'lp', None) for t in tv['valleys']]
    out = [(kname(ia), v[ia], 'axis', ia)]
    if ia != 90: out.append(('Γ', v[0], 'ref', 90))
    if m.startswith('Sn'):
        out += [(nm, v[i], 'lp', i) for nm, i in SNVAL.items() if nm != 'Γ']
    else:
        x, _ = rf.path_x(m); nel = NEL[m]
        E = cv.read_path_eigen(os.path.join(rf.base(m), 'c2_hse', 'path_eigen.txt'))
        _, W = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz')); fs = rf.cation_s(W)[:, nel - 2:nel].mean(1)
        for i, _m in cv.extrema(x, E[:, nel - 1], 'max', 600):
            if (0 < i < 33 or 97 < i < 130) and fs[i] >= LP_MIN and v[i] > -260: out.append((kname(i), v[i], 'lp', i))
    return out


def fig6(out):
    """左列 Sn₂Sb₂Te₅、右列 Sn₂Bi₂Te₅：(a)(b) 价带口袋沿 c 俯视，(c)(d) 三维透视；(e) 四种材料的价带谷间能量差 ΔE。"""
    data = {m: pockets(m) for m in ('Sn2Sb2Te5', 'Sn2Bi2Te5')}
    fig = plt.figure(figsize=(W2, 7.2))
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 0.82], height_ratios=[1, 1.0, 1.35], wspace=0.04, hspace=0.10)
    zs = 4
    for k, m in enumerate(data):
        at = fig.add_subplot(gs[0, k]); bz_top(at, m, *data[m], f'({"ab"[k]})')
        # 三维面板与俯视面板对齐：中心放在俯视图六边形中心的正下方，宽度相同
        at.apply_aspect(); pt = at.get_position(); x0, x1 = at.get_xlim()
        cx = pt.x0 + (0 - x0) / (x1 - x0) * pt.width; c1 = gs[1, k].get_position(fig)
        wd = 2.1 / (x1 - x0) * pt.width * 1.25
        ax = fig.add_axes([cx - wd / 2, c1.y0, wd, c1.height], projection='3d')
        bz_3d(ax, m, *data[m], zs=zs)
        ax.text2D(0.5, 1.04, f'{TEX[m]}: 3D view ($k_z$ ×{zs})', transform=ax.transAxes, fontsize=8, ha='center', va='bottom')
        ax.text2D(-0.08, 1.04, f'({"cd"[k]})', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')
    # (e) 四种材料的价带能级（HSE06+SOC）：Sb → Bi 时 Γ 相对孤对电子谷下降；Pb 系价带顶移到 A（仍是单谷），Sn 系移到孤对电子多谷
    rows = crossover(); by = {d['m']: d for d in rows}
    ax = fig.add_subplot(gs[2, :2]); lv_rows = []
    X = {'Pb2Sb2Te5': 0.0, 'Pb2Bi2Te5': 1.0, 'Sn2Sb2Te5': 2.9, 'Sn2Bi2Te5': 3.9}; hw = 0.26
    ax.axhspan(-2 * KT300, 0, color='#EEF3F8', zorder=0, lw=0)
    L = {m: levels(m) for m in MATS}
    for m in MATS:
        x0 = X[m]
        for nm, E, kind, _i in L[m]:
            col = C_A if kind == 'lp' else C_B; ref = kind == 'ref'
            ax.plot([x0 - hw, x0 + hw], [E, E], color=col, lw=0.9 if ref else 1.6, ls=(0, (2, 1.5)) if ref else '-', zorder=2,
                    alpha=0.7 if ref else 1)
            if kind != 'ref':
                ax.plot(x0, E, MK.get(nm, 'o'), ms={'o': 4.2, 's': 3.6, 'D': 3.4, 'p': 4.0, 'h': 4.0}[MK.get(nm, 'o')],
                        color=col, mfc='white' if nm in ('R₁', 'R₂') else col, mew=0.8, zorder=3)
            lv_rows.append([m, nm, kind, r(E, 1)])
        for nm, E, kind, _i in L[m]:                               # 轴上的谷标在左侧；不是谷的 Γ 点标在右侧
            if kind == 'axis': ax.text(x0 - hw - 0.04, E, nm, ha='right', va='center', fontsize=6.5, color=C_B)
            if kind == 'ref': ax.text(x0 + hw + 0.04, E, 'Γ (not a max)', ha='left', va='center', fontsize=6.0, color=C_B)
        # 孤对电子谷（紫）标名称：Sb 化合物标在左侧、Bi 化合物标在右侧（避开连线）；相距太近的几条线把名称上下错开，用短线连到各自的能级
        lp = sorted([(E, nm) for nm, E, kind, _i in L[m] if kind == 'lp'], reverse=True)
        ys = [E for E, _ in lp]; gap_min = 13.0
        for _ in range(50):
            for k in range(1, len(ys)):
                if ys[k - 1] - ys[k] < gap_min:
                    mid = (ys[k - 1] + ys[k]) / 2; ys[k - 1], ys[k] = mid + gap_min / 2, mid - gap_min / 2
        side = -1 if 'Sb' in m else 1
        xe = x0 + side * hw; xt = xe + side * 0.16
        for (E, nm), yt in zip(lp, ys):
            ax.plot([xe, xe + side * 0.10], [E, yt], color=C_A, lw=0.5, zorder=2)
            ax.text(xt, yt, MATH.get(nm, nm), ha='left' if side > 0 else 'right', va='center', fontsize=6.5, color=C_A)
        d = by[m]['d']
        ax.text(x0, 30, f"ΔE = {d:+.0f}".replace('-', '−'), ha='center', va='center', fontsize=6.5, color=C_A if d > 0 else C_B)
    for (a1, a2) in (('Pb2Sb2Te5', 'Pb2Bi2Te5'), ('Sn2Sb2Te5', 'Sn2Bi2Te5')):   # 同一个态连虚线：Γ 点、Σ 谷
        for nm, col in (('Γ', C_B), ('Σ', C_A)):
            e1 = next(E for n, E, k, _ in L[a1] if n.rstrip('′″') == nm); e2 = next(E for n, E, k, _ in L[a2] if n.rstrip('′″') == nm)
            ax.plot([X[a1] + hw, X[a2] - hw], [e1, e2], color=col, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xticks(list(X.values()), [TEX[m] for m in X]); ax.tick_params(axis='x', top=False, length=0)
    ax.set_xlim(-0.6, 5.0); ax.set_ylim(-262, 52); ax.set_yticks(range(-250, 1, 50))
    ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (meV)')
    for xm, t in ((0.5, 'Pb series: Sb → Bi'), (3.4, 'Sn series: Sb → Bi')):
        ax.text(xm, -0.085, t, ha='center', va='top', fontsize=7, transform=ax.get_xaxis_transform())   # 框外，材料名下方
    ax.axvline(2.25, color=GRID, lw=0.6)
    ax.text(4.9, -2 * KT300 / 2, r'$2k_\mathrm{B}T$', ha='right', va='center', fontsize=6.3)
    ax.set_title('Valence-band levels of the four compounds (HSE06+SOC)', fontsize=8, pad=4)
    dump('Fig6_valence_pockets', 'levels_HSE.csv', ['material', 'state', 'kind (axis: valley on the Γ–A axis; lp: lone-pair valley; '
         'ref: Γ point, not a valley)', 'E-EVBM_meV'], lv_rows)
    panel(ax, '(e)')
    # 右侧：(a)–(d) 的图例
    lg = fig.add_subplot(gs[0:2, 2]); lg.axis('off')
    h = [Line2D([], [], ls='', marker='o', ms=5, color=C_B, label='Γ / Δ: Γ–A axis state (p orbitals)'),
         Line2D([], [], ls='', marker='o', ms=5, color=C_A, label='Σ: lone-pair valley on Γ–M'),
         Line2D([], [], ls='', marker='s', ms=4.5, color=C_A, mfc='white', label=r'R$_1$: lone-pair valley on A–L, near A'),
         Line2D([], [], ls='', marker='D', ms=4.0, color=C_A, mfc='white', label=r'R$_2$: lone-pair valley on A–L, near L'),
         Line2D([], [], ls='', marker='p', ms=4.5, color=C_A, label=r'′, ″: true maximum off the line (mirror-plane search)'),
         Line2D([], [], ls='', marker='o', ms=5, color='#555555', label=r'filled: $k_z$ = 0 plane or inside'),
         Line2D([], [], ls='', marker='o', ms=5, color='#555555', mfc='white', label=r'open: $k_z$ = ½ zone face'),
         Line2D([], [], ls='', marker='o', ms=3.4, color=C_A, alpha=0.45, label=r'faded: more than $2k_\mathrm{B}T$ below VBM')]
    lg.legend(handles=h, loc='upper left', bbox_to_anchor=(0.04, 0.88), title='(a)–(d): reciprocal (k) space;\neach symbol is a valley, not an atom', title_fontsize=7, handletextpad=0.3, labelspacing=0.6, borderaxespad=0)
    lg.text(0.06, 0.12, 'Labels: line names of the hexagonal zone\n(Σ: Γ–M, R: A–L, Δ: Γ–A); energies in meV\nrelative to the VBM; ×n = equivalent pockets;\n* c4 grid point, not refined (lower bound).',
            transform=lg.transAxes, fontsize=6.3, va='bottom', ha='left', color='#333333')
    save(fig, out, 'Fig6_valence_pockets_crossover')
    prow = []
    for m, (B, P) in data.items():
        for g in P:
            for f, c in zip(g['frac'], g['cart']):
                prow.append([m, g['name'], g['i'], *[r(x, 4) for x in f], *[r(x, 4) for x in c],
                             'mirror plane, general kz' if (g['i'] is None) else ('kz=1/2 zone face' if on_face(f) else ('kz=0' if np.isclose(f[2], 0) else 'Γ–A axis (Δ line)')),
                             r(g['E'], 1), len(g['frac']), 'c4 grid point (lower bound)' if g.get('grid') else ('c5 true extremum' if g['i'] is None else 'band path')])
    dump('Fig6_valence_pockets', 'pocket_positions.csv', ['material', 'valley', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3',
         'kx_1/A', 'ky_1/A', 'kz_1/A', 'plane', 'E-EVBM_meV (HSE06+SOC)', 'n_equivalent_pockets', 'source'], prow)
    dump('Fig6_valence_pockets', 'valley_offset_HSE.csv',
         ['material', 'functional', 'E_Gamma_meV', 'E_GammaA_axis_max_meV', 'axis_position', 'axis_path_k_index', 'axis_n_pockets',
          'E_lone_pair_max_meV', 'lone_pair_valley', 'lone_pair_k_frac_1', 'lone_pair_k_frac_3', 'lone_pair_n_pockets', 'lone_pair_cation_s_percent',
          'dE_valley_offset_meV (E_lone_pair - E_axis)', 'E_lone_pair_minus_E_Gamma_meV', 'source'],
         [[d['m'], d['fn'], r(d['EG'], 1), r(d['Ea'], 1), kname(d['ia']), d['ia'], d['na'], r(d['El'], 1), d['lname'], r(d['kl'][0], 4), r(d['kl'][2], 4),
           d['nl'], r(100 * d['s'], 1), r(d['d'], 1), r(d['El'] - d['EG'], 1),
           'true VBM from c4_mirror_v3 + c5_star (mirror plane)' if d['il'] is None else 'band path c2 (lone-pair valley is a lower bound)'] for d in rows])
    return {m: [(g['name'], len(g['frac']), r(g['E'], 1)) for g in P] for m, (B, P) in data.items()} | \
           {f"{d['m']}|{d['fn']}": r(d['d'], 1) for d in rows}


# ------------------------------------------------------------------ Fig. 7
def fig7(out):
    K = kpath(); rows = []
    fig, ax = plt.subplots(figsize=(W1, 2.1))
    kT2 = 2 * 1000 * mv.KB * 300
    ax.axvspan(-kT2, 0, color='#EEF3F8', zorder=0, lw=0)
    ax.text(-kT2 / 2, 3.55, r'$2k_\mathrm{B}T$', ha='center', va='center', fontsize=6.5)
    for j, m in enumerate(MATS[::-1]):
        nel = NEL[m]; v = {}
        for side, nm in (('c2_pbe', 'PBE+SOC'), ('c2_hse', 'HSE06+SOC')):
            lp = rf.lone_pair_offset(m, side); c = (lp[0], f'path {lp[1]}', tuple(K[lp[1]]), lp[2])
            tq = true_path_bands(m, side)
            if tq and side == 'c2_hse':                          # Sn 系 HSE：镜面内三维搜索的真实谷，相对真实价带顶
                t0 = true_valleys(m)['valleys'][0]; c = (t0['E'], t0['name'] + ' (mirror-plane search)', t0['kfrac'], t0['s'])
            elif tq:                                             # Sn 系 PBE：附加路径上的孤对电子极大也参与比较，零点取两段路径上的最高价带
                e = 1000 * tq['E'][:, nel - 1]; sh = max(0.0, float(np.nanmax(e))); c = (c[0] - sh,) + c[1:]
                fq = tq['fs'][:, nel - 2:nel].mean(1)
                for i in range(1, len(e) - 1):
                    if e[i] >= e[i - 1] and e[i] >= e[i + 1] and fq[i] >= 0.15 and e[i] - sh > c[0]:
                        c = (e[i] - sh, f'c7:{i} {tq["labels"].get(i, "")}'.strip() + ' (true-valley path)', tuple(tq['k'][i]), fq[i])
            v[nm] = c; rows.append([m, nm, r(c[0], 1), c[1]] + [r(q, 4) for q in c[2]] + [r(100 * c[3], 1)])
        p, h = v['PBE+SOC'], v['HSE06+SOC']
        ax.plot([p[0], h[0]], [j, j], color='#BDBDBD', lw=1.2, zorder=1)
        ax.plot(p[0], j + 0.08, 'o', ms=4.5, color=C_PBE, mec='white', mew=0.6, zorder=3, label='PBE+SOC' if j == 0 else None)
        ax.plot(h[0], j - 0.08, 'o', ms=4.5, color=C_HSE, mec='white', mew=0.6, zorder=3, label='HSE06+SOC' if j == 0 else None)
        ax.text(h[0] - 7, j - 0.08, f'{h[0]:.0f}'.replace('-', '−'), ha='right', va='center', fontsize=6.3)
    dump('Fig7_lone_pair_valley', 'highest_lone_pair_valley.csv',
         ['material', 'functional', 'E-EVBM_meV', 'where', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'cation_s_percent'], rows[::-1])
    ax.set_yticks(range(4), [TEX[m] for m in MATS[::-1]]); ax.set_ylim(-0.6, 3.85)
    ax.set_xlim(-275, 15); ax.set_xlabel('Highest lone-pair valence valley, $E - E_{\\mathrm{VBM}}$ (meV)')
    ax.tick_params(axis='y', length=0, right=False); ax.legend(loc='lower left', handletextpad=0.2, borderaxespad=0.2)
    save(fig, out, 'Fig7_lone_pair_valley_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 8
# 图 8(c)(d) 的配色：色相按 A 位（Sn 紫，与图 5 的 Sn 5s 孤对电子同色；Pb 灰），深浅按 B 位（Bi 深、Sb 浅）
# 图 8、图 9：每种材料一个颜色，取自用户用 Origin 画的第一版 S(T) 图（实验数据/四种材料S_T_第一版.webp，2026-11-04 用户要求）。
# 红、绿两色对红绿色弱区分不足（validate_palette.js：deutan ΔE 5.7），所以同时用 Origin 式的“线 + 符号”，符号与第一版图相同。
MCOL = {'Pb2Sb2Te5': '#F24040', 'Pb2Bi2Te5': '#525252', 'Sn2Sb2Te5': '#39AE6B', 'Sn2Bi2Te5': '#1B72DF'}
MMARK = {'Pb2Sb2Te5': 'o', 'Pb2Bi2Te5': 's', 'Sn2Sb2Te5': 'v', 'Sn2Bi2Te5': '^'}


def fig8(out):
    """300 K 面内 S 随载流子浓度（理论）：四种材料 PBE+SOC 密网格（28×28×7/6）BoltzTraP2、不加剪刀的直接结果
    （22_能带基准_1-0/bt2_pbe_dense，pbe_dense()），不用任何估计法（2026-11-04 用户要求）。(a) p 型，(b) n 型。
    浓度在 BoltzTraP2 的化学势网格上按 log|N| 线性取值。HSE06+SOC 的直接结果（c8_hse_dense）完成后加入。
    只画计算，不画实验值（与实验的对照在图 9、图 S3；2026-11-04 用户要求）。"""
    T0 = 300; rows, tab = [], []
    fig, A = plt.subplots(1, 2, figsize=(W2, 2.75), gridspec_kw=dict(wspace=0.26))
    Ns = np.logspace(18, np.log10(3e21), 73)
    for i, (ax, sign) in enumerate(zip(A, (1, -1))):
        for m in MATS:
            T, f = pbe_dense(m); it = list(T).index(T0)
            v = np.array([sign * f(sign * N, it) for N in Ns])
            ax.plot(Ns, v, color=MCOL[m], lw=1.1, marker=MMARK[m], ms=4.2, markevery=(3 * MATS.index(m) % 6, 6), label=TEX[m], zorder=3)
            rows += [[m, 'pn'[i], f'{N:.3e}', r(x, 2)] for N, x in zip(Ns, v)]
            if i == 0:
                tab.append([m] + [r(f(N, it), 1) for N in (3e19, 1e20, -3e19, -1e20)])
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(Ns[0], Ns[-1]); ax.set_ylim(3, 600)   # 对数纵轴：样品的 S 只有 5–90 µV/K
        ax.set_yticks([5, 10, 20, 50, 100, 200, 500], ['5', '10', '20', '50', '100', '200', '500']); ax.minorticks_off()
        ax.set_xlabel(('Hole' if sign > 0 else 'Electron') + r' concentration (cm$^{-3}$)')
        ax.set_ylabel(r'$S$ (μV K$^{-1}$)' if sign > 0 else r'$|S|$ (μV K$^{-1}$)')
        ax.set_title(f'{"p" if sign > 0 else "n"}-type, 300 K, PBE+SOC + BoltzTraP2', fontsize=8, pad=4)
        ax.axvline(3e19, color='#999999', lw=0.6, ls=(0, (1, 1.5)), zorder=0)
        panel(ax, f'({"ab"[i]})')
    h, l = A[0].get_legend_handles_labels()
    A[0].legend(h, l, loc='lower left', fontsize=6.5, handlelength=2.2, labelspacing=0.3)
    d = os.path.join(DATA, 'Fig8_Seebeck') if DATA else None
    if d and os.path.isdir(d): shutil.rmtree(d)          # 旧版（HSE 估计、多谷模型的谷份额）的文件一并清掉
    dump('Fig8_Seebeck', 'S_vs_concentration_300K.csv',
         ['material', 'carrier_type', 'concentration_cm-3', '|S|_PBE+SOC_BoltzTraP2_uV/K (28x28x7/6 dense grid, no scissor)'], rows)
    dump('Fig8_Seebeck', 'S_300K_table6.csv', ['material', 'S_p_3e19_uV/K', 'S_p_1e20_uV/K', 'S_n_3e19_uV/K', 'S_n_1e20_uV/K'], tab)
    save(fig, out, 'Fig8_Seebeck_vs_concentration')
    return tab


LP_MIN = 0.15        # 孤对电子谷的判据：谷极值点处阳离子 s 成分 ≥ 15%（与图 4 相同）


def valley_index(name):
    """multivalley_seebeck 的谷名 → c2 路径上的 k 序号（取成分用）。"""
    if name.startswith('mirror'): return None                       # 镜面内的真实谷（c4/c5），不在路径上
    if name.startswith('#'): return int(name[1:].split()[0])
    if name.startswith('Γ–A'): return 90 + round(float(name.split('±')[1]) * 14)
    return mv.POINTS[name][0][0]


def valley_shares():
    """表 4 的电导份额：多谷模型（HSE、CRTA、300 K、面内、3e19）中每个谷占总电导（含少数载流子）的比例，
    以及谷极值点处的阳离子 s 成分；p 型孤对电子谷 = 价带谷中阳离子 s ≥ 15% 者，n 型“Γ 导带”= Γ 或 Γ–A 驼峰处的导带谷。"""
    K = kpath(); rows, summary, scan = [], [], []
    for m in MATS:
        c3 = mv.c3_usable(m)[0]
        V, gap = mv.build(m, 'v', mv.WIN, c3); C, _ = mv.build(m, 'c', mv.WIN)
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fs = rf.cation_s(Wp); nel = NEL[m]
        nm = {'v': knames([valley_index(v['name']) for v in V]), 'c': knames([valley_index(v['name']) for v in C])}
        tv = {tuple(np.round(t['kfrac'], 4)): t['name'] for t in (true_valleys(m) or {'valleys': []})['valleys']}
        for N, lab in ((3e19, 'p'), (-3e19, 'n')):
            _, EF, parts = mv.solve(V, C, gap, N, 300)
            lp = 0.0
            for (side, name, w), v in zip(parts, V + C):
                i = valley_index(name); b = nel - 1 if side == 'v' else nel
                cs = v['s'] if i is None else fs[i, b]; kf = v['kfrac'] if i is None else K[i]
                lab_v = tv.get(tuple(np.round(kf, 4)), name) if i is None else nm[side][i]
                flag = (side == 'v' and lab == 'p' and cs >= LP_MIN) or (side == 'c' and lab == 'n' and (name in ('Γ',) or name.startswith('Γ–A')))
                lp += w if flag else 0.0
                rows.append([m, lab, r(1000 * EF, 1), 'valence' if side == 'v' else 'conduction', lab_v, '' if i is None else i]
                            + [r(x, 4) for x in kf] + [r(1000 * v['off'], 1), v['g'], r(100 * cs, 1), r(100 * w, 2), 'yes' if flag else 'no'])
            summary.append([m, lab, r(100 * lp, 1)])
        # p 型孤对电子谷份额随空穴浓度的变化（判据同上）
        lpv = [v['name'] for v in V if (v['s'] if valley_index(v['name']) is None else fs[valley_index(v['name']), nel - 1]) >= LP_MIN]
        for N in (1e19, 3e19, 1e20, 3e20):
            S, EF, pp = mv.solve(V, C, gap, N, 300)
            scan.append([m, f'{N:.0e}', r(1000 * EF, 1), r(100 * sum(w for side, name, w in pp if side == 'v' and name in lpv), 1), r(S, 1)])
    dump('Fig8_Seebeck', 'valley_conductivity_shares.csv',
         ['material', 'carrier_type', 'EF_minus_VBM_meV', 'band', 'valley (line/point name, see kname)', 'k_index_on_path', 'k_frac_1', 'k_frac_2', 'k_frac_3',
          'offset_from_band_edge_meV', 'degeneracy_g', 'cation_s_percent', 'conductivity_share_percent',
          'counted (p: lone-pair valence valley, cation s >= 15%; n: Gamma conduction valley)'], rows)
    dump('Fig8_Seebeck', 'lone_pair_share_summary.csv',
         ['material', 'carrier_type', 'share_percent (p: lone-pair valence valleys; n: Gamma conduction valley)'], summary)
    dump('Fig8_Seebeck', 'lone_pair_share_vs_concentration.csv',
         ['material', 'hole_conc_cm-3', 'EF_minus_VBM_meV', 'lone_pair_valley_share_percent', 'S_model_raw_uV/K (HSE path model, not calibrated)'], scan)
    print('shares', summary)
    c3 = mv.c3_load('Sn2Bi2Te5')
    if c3:
        rows = []
        for k, v in c3.items():
            H = np.array(v['H_eVA2_line_perp_c']); g = np.array(v['gradient'])
            ok, why = mv.c3_verdict(v, 'v')
            dk = -np.linalg.solve(H, g); dE = -0.5 * g @ np.linalg.solve(H, g)
            rows.append([kname(int(k)), k, r(v['E0'], 5)] + [r(x, 4) for x in g] + [r(x, 4) for x in v['principal']]
                        + [r(x, 4) for x in dk] + [r(1000 * dE, 1), 'yes' if ok else 'no', why])
        dump('Fig8_Seebeck', 'Sn2Bi2Te5_c3_valley_check.csv',
             ['valley', 'valley_k_index', 'E_center_eV', 'grad_line_eVA', 'grad_perp_eVA', 'grad_c_eVA',
              'principal_mass_1_me', 'principal_mass_2_me', 'principal_mass_3_me',
              'stationary_shift_line_1/A', 'stationary_shift_perp_1/A', 'stationary_shift_c_1/A',
              'quadratic_extrapolated_dE_meV', 'used_in_model', 'reason_not_used'], rows)


# ------------------------------------------------------------------ Fig. 9
# 实验数据（用户提供，21_阶段报告与交接/实验数据/）：四个样品的原始 T、σ（10⁴ S/m）、S（2026-11-02），与室温霍尔
EXPDIR = os.path.join(ROOT, '21_阶段报告与交接', '实验数据')


def read_exp(m):
    """T（K）、σ（S/m）、S（µV/K）。Pb₂Sb₂Te₅ 的文件表头写“单位未给出”，σ 同为 10⁴ S/m（实验数据/README.md）。"""
    rows = list(csv.reader(open(os.path.join(EXPDIR, f'{m}_T_sigma_S.csv'), encoding='utf-8-sig')))[1:]
    return [(float(t), float(s) * 1e4, float(S)) for t, s, S in rows]


def read_hall():
    """室温霍尔：N（cm⁻³，正为空穴、负为电子）、μ_H（cm²/V·s）、仪器是否把结果标红。"""
    out = {}
    for d in csv.DictReader(open(os.path.join(EXPDIR, '霍尔_20261102.csv'), encoding='utf-8-sig')):
        out[d['sample']] = dict(N=(1 if d['carrier_type'] == 'P' else -1) * float(d['n_cm-3']), mu=float(d['mu_H_cm2/Vs']),
                                flagged=d['instrument_flag'] != 'none')
    return out


EXP = {m: [(t, S) for t, _, S in read_exp(m)] for m in MATS}
HALL = read_hall()


def model(m, side, tau='crta'):
    V, gap = mv.build(m, 'v', mv.WIN, None, side); C, _ = mv.build(m, 'c', mv.WIN, side=side)
    return lambda N, T: mv.solve(V, C, gap, N, T, tau)[0]


C_ADP = '#56B4E9'          # 图 9：HSE 能带 + 声学声子散射（浅蓝细线）
SCOL = {'H': (C_HSE, '-', 1.3), 'A': (C_ADP, '-', 0.9), 'P': (C_PBE, DASH, 1.1)}


def solve_N(f, target, T, lo, hi):
    for _ in range(70):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if (f(mid, T) - target) * (f(lo, T) - target) > 0 else (lo, mid)
    return (lo + hi) / 2


def solve_log(f, data, lo, hi, sign=1):
    """在 log10 浓度 [lo, hi] 内找 f(sign·10^x, T0) = S0 的全部根（T0、S0 取 data 的第一个实验点）。
    双极区 S 随浓度不单调，可能有两个根；取与整条实验曲线最大偏差最小的那个。"""
    T0, S0 = data[0]
    g = lambda x: f(sign * 10.0 ** x, T0) - S0
    xs = np.linspace(lo, hi, 31); gs = [g(x) for x in xs]; roots = []
    for a, b, ga, gb in zip(xs[:-1], xs[1:], gs[:-1], gs[1:]):
        if ga * gb > 0: continue
        for _ in range(50):
            mid = (a + b) / 2; gm = g(mid)
            a, b, ga = (mid, b, gm) if gm * ga > 0 else (a, mid, ga)
        roots.append(sign * 10.0 ** ((a + b) / 2))
    return min(roots, key=lambda N: max(abs(f(N, t) - s) for t, s in data))


def s_methods(m, sign):
    """与图 8 同口径的三种算法，返回 S(N, T)（N 带符号，正为空穴）：'P' PBE+SOC BoltzTraP2 全区（常数弛豫时间）；
    'H' HSE 估计 = BoltzTraP2 × 多谷模型 HSE/PBE 比值（常数弛豫时间）；'A' 同上，比值的分子换成声学声子散射。"""
    b = bt2_signed(m, sign); gh, gp, ga = model(m, 'c2_hse'), model(m, 'c2_pbe'), model(m, 'c2_hse', 'adp')
    return {'H': lambda N, t: b(N, t) * gh(N, t) / gp(N, t), 'A': lambda N, t: b(N, t) * ga(N, t) / gp(N, t), 'P': lambda N, t: b(N, t)}


SNAME = {'H': 'HSE06+SOC', 'A': 'HSE06+SOC, acoustic phonon', 'P': 'PBE+SOC'}


# 单抛物带（声学声子散射，λ = 0）的 Pisarenko 关系：文献中常用来由 S 与霍尔浓度读出态密度有效质量（May & Snyder；Pei 等）
KBe = 86.17333262     # k_B/e，µV/K


def fd(j, eta):
    return integrate.quad(lambda x: x ** j / (1 + np.exp(x - eta)), 0, max(60.0, eta + 60.0))[0]


def spb_S(eta):
    return KBe * (2 * fd(1, eta) / fd(0, eta) - eta)


def spb_n(eta, ms, T=300):
    """cm⁻³：n = 4π(2 m* k_B T / h²)^{3/2} F_{1/2}(η)。"""
    return 4 * np.pi * (2 * ms * 9.1093837e-31 * 1.380649e-23 * T / 6.62607015e-34 ** 2) ** 1.5 * fd(0.5, eta) * 1e-6


def spb_mass(S, n, T=300):
    """由 |S|（µV/K）与浓度（cm⁻³）反推单抛物带态密度有效质量（m_e）。"""
    eta = optimize.brentq(lambda x: spb_S(x) - abs(S), -15, 60)
    return (abs(n) / spb_n(eta, 1.0, T)) ** (2 / 3)


SPB_ETA = np.linspace(-3, 40, 120)
SPB_CURVE = [(spb_S(e), e) for e in SPB_ETA]


# 单个谷的态密度质量（HSE，路径或镜面搜索拟合的三个主质量的几何平均），作 Pisarenko 图中的“单谷”参考线
SINGLE = {'Pb2Sb2Te5': (0.10, r'$\Gamma$ valley'), 'Pb2Bi2Te5': (0.23, r'$\Gamma$ valley'),
          'Sn2Sb2Te5': (0.14, r'$\Gamma$ valley'), 'Sn2Bi2Te5': (0.13, 'one lone-pair pocket')}


def fig9(out):
    """与实验对照（文献常用的 Pisarenko 图，单个面板）：300 K 的 |S| 对室温霍尔浓度（双对数）。灰线为单抛物带（声学声子散射）
    在 m* = 0.1、0.3、1、3 m_e 下的关系，作读数用的“刻度”；实心圆为样品；空心方块为同一霍尔浓度下 PBE+SOC 密网格 BoltzTraP2
    （不加剪刀）的直接计算，用同样方法读出 m*。不画任何估计值；HSE06+SOC 直接结果完成后加入。
    Pb₂Sb₂Te₅ 为补偿体系，单一载流子的关系不适用，不画（附录 A）。"""
    T = 300; lines, samples = [], []
    fig, ax = plt.subplots(figsize=(4.2, 3.3))
    ss = np.array([S for S, _ in SPB_CURVE]); lo, hi, ylo, yhi = 1e19, 1e22, 4, 300
    for ms in (0.1, 0.3, 1.0, 3.0):
        nn = np.array([spb_n(e, ms, T) for e in SPB_ETA]); ok = (nn >= lo) & (nn <= hi) & (ss >= ylo) & (ss <= yhi)
        ax.plot(nn[ok], ss[ok], color='#A6A6A6', lw=0.8, zorder=1)
        lines += [[r(ms, 3), f'{a_:.3e}', r(b_, 2)] for a_, b_ in zip(nn, ss)]
        k = np.where(ok)[0][0]                                   # 标在线段进入画框的一端
        ax.text(nn[k] * 1.08, ss[k] * (0.93 if ss[k] > 0.9 * yhi else 1.04), f'{ms:g} $m_e$', fontsize=6.5, color='#7F7F7F',
                ha='left', va='top' if ss[k] > 0.9 * yhi else 'bottom')
    lab = {'Pb2Bi2Te5': ('n', (8, 2)), 'Sn2Bi2Te5': ('p', (-8, 6)), 'Sn2Sb2Te5': ('p', (8, -2)), 'Pb2Sb2Te5': ('p', (-8, -2))}
    for m, (typ, off) in lab.items():
        h = HALL[m]; sign = 1 if h['N'] > 0 else -1; c = MCOL[m]; Tg, f = pbe_dense(m)
        n0 = abs(h['N']); S0 = abs(EXP[m][0][1]); sp = abs(f(h['N'], list(Tg).index(T)))
        ms_exp, ms_p = spb_mass(S0, n0, T), spb_mass(sp, n0, T)
        ax.plot([n0, n0], [S0, sp], color=c, lw=0.8, ls=(0, (1.5, 1.5)), zorder=2)   # 同一样品的实验点与计算点相连
        ax.plot(n0, sp, 's', ms=8.5, mfc='white', mec=c, mew=1.2, zorder=3)          # 计算画得大、实验画得小，重合时两者都看得见
        ax.plot(n0, S0, 'o', ms=5, mfc=c, mec=c, zorder=4)
        comp = m == 'Pb2Sb2Te5'                                   # 电子、空穴补偿：单抛物带的 m* 不适用，只画点
        ax.annotate(f'{TEX[m]} ({typ})' + ('*' if h['flagged'] else '') + ('\ncompensated, no $m^*$' if comp else f'\n$m^*$ = {ms_exp:.2f} $m_e$'),
                    (n0, S0), xytext=off, textcoords='offset points', ha='left' if off[0] > 0 else 'right', va='center', fontsize=6.5, color='#333333',
                    bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.85), zorder=6)
        samples.append([m, typ, f'{h["N"]:.3e}', 'yes' if h['flagged'] else 'no', r(EXP[m][0][1], 2), '' if comp else r(ms_exp, 3),
                        r(sign * sp, 1), r(ms_p, 3), SINGLE[m][0],
                        'electrons and holes compensate; single-band m* not applicable to the experiment' if comp else ''])
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, hi); ax.set_ylim(ylo, yhi)
    ax.set_yticks([5, 10, 20, 50, 100, 200], ['5', '10', '20', '50', '100', '200']); ax.minorticks_off()
    ax.set_xlabel(r'Hall carrier concentration (cm$^{-3}$)'); ax.set_ylabel(r'$|S|$ at 300 K (μV K$^{-1}$)')
    hh = [plt.Line2D([], [], ls='', marker='o', ms=5, mfc='#555555', mec='#555555'),
          plt.Line2D([], [], ls='', marker='s', ms=7.5, mfc='white', mec='#555555', mew=1.2),
          plt.Line2D([], [], color='#A6A6A6', lw=0.8)]
    ax.legend(hh, ['Experiment', 'PBE+SOC + BoltzTraP2 (calc.)', 'Single parabolic band'], loc='upper right', fontsize=6.5, handletextpad=0.3)
    d = os.path.join(DATA, 'Fig9_Pisarenko') if DATA else None
    if d and os.path.isdir(d): shutil.rmtree(d)
    dump('Fig9_Pisarenko', 'single_parabolic_band_lines.csv', ['m_DOS_me', 'concentration_cm-3', '|S|_uV/K'], lines)
    dump('Fig9_Pisarenko', 'samples_Pisarenko_mass.csv', ['material', 'carrier_type', 'Hall_conc_signed_cm-3', 'Hall_flagged', 'S_exp_300K_uV/K', 'm_DOS_exp_SPB_me',
                                                        'S_PBE+SOC_BoltzTraP2_at_Hall_conc_uV/K', 'm_DOS_PBE_SPB_me', 'single_valley_m_DOS_HSE_band_fit_me', 'note'], samples)
    save(fig, out, 'Fig9_Pisarenko_experiment')
    return samples


def st_methods(m):
    """附录 A 拟合 S(T) 用的算法，返回 {'H': S_HSE(N, T), 'P': S_PBE(N, T)} 与浓度的搜索方式。N 带符号（正为空穴）。
    Sn 系：重掺杂 p 型，PBE = BoltzTraP2，HSE = BoltzTraP2 × 多谷模型 HSE/PBE 比值（与图 8 同口径）。
    Pb₂Bi₂Te₅：PBE = BoltzTraP2；HSE 拟合出的浓度低于 5×10¹⁸ cm⁻³（BoltzTraP2 插值不可靠），直接用多谷模型。
    Pb₂Sb₂Te₅：近补偿，净掺杂 N 很小、可正可负，两种能带都直接用多谷模型（含电子与空穴）。"""
    if m.startswith('Sn'):
        f = s_methods(m, 1); return {'H': f['H'], 'P': f['P']}, ('log', 1, 19.0, 21.8)
    if m == 'Pb2Bi2Te5':
        return {'H': model(m, 'c2_hse'), 'P': s_methods(m, -1)['P']}, ('log', -1, 17.0, 20.5)
    return {'H': model(m, 'c2_hse'), 'P': model(m, 'c2_pbe')}, ('lin', 0, -3e19, 3e19)


def fit_N(f, data, how):
    """只调一个常数——净载流子浓度 N，使 300–620 K 的计算 S 与实验的均方根偏差最小。先在网格上找，再局部细化。"""
    kind, sign, lo, hi = how
    to_N = (lambda x: sign * 10.0 ** x) if kind == 'log' else (lambda x: x)
    rms = lambda x: np.sqrt(np.mean([(f(to_N(x), t) - S) ** 2 for t, S in data]))
    xs = np.linspace(lo, hi, 61 if kind == 'log' else 241); vs = [rms(x) for x in xs]; k = int(np.argmin(vs))
    res = optimize.minimize_scalar(rms, bounds=(xs[max(k - 1, 0)], xs[min(k + 1, len(xs) - 1)]), method='bounded', options=dict(xatol=1e-4 if kind == 'log' else 1e13))
    return to_N(res.x), float(res.fun)


PBEDENSE = os.path.join(ROOT, '22_能带基准_1-0', 'bt2_pbe_dense')


def pbe_dense(m):
    """PBE+SOC 密网格（28×28×7/6）BoltzTraP2、不加剪刀（22_能带基准_1-0/bt2_pbe_dense）：返回计算温度 T 与 f(N, i) =
    第 i 个计算温度、带符号浓度 N（+ 空穴）处的面内 S（µV/K）。浓度方向在 BoltzTraP2 的化学势网格上按 log|N| 线性取值，温度不插值。"""
    d = json.load(open(os.path.join(PBEDENSE, m, 'boltztrap_crta.json'), encoding='utf-8'))
    T = np.array(d['temperatures_K']); C = np.array(d['carrier_conc_cm-3'], float)
    S = (np.array(d['seebeck_xx_V_per_K'], float) + np.array(d['seebeck_yy_V_per_K'], float)) / 2 * 1e6

    def f(N, i):
        sg = 1 if N > 0 else -1; ok = sg * C[i] > 0; x = np.log10(sg * C[i][ok]); o = np.argsort(x)
        assert x.min() <= np.log10(abs(N)) <= x.max()
        return float(np.interp(np.log10(abs(N)), x[o], S[i][ok][o]))
    return T, f


def figS3(out):
    """附录 A：实验 S(T) 与计算，计算直接取室温霍尔浓度，不调任何参数（2026-11-03 用户确定）。
    计算用 PBE+SOC 密网格 BoltzTraP2、不加剪刀，曲线就是 BoltzTraP2 的计算温度（每 25 K），不用任何估计法（2026-11-04 用户要求）；
    均方根偏差在测量温度上比较，计算值按 25 K 网格线性内插。HSE06+SOC 的直接结果（c8_hse_dense）完成后再加。"""
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.9), gridspec_kw=dict(hspace=0.42, wspace=0.24))
    ylim = {'Pb2Sb2Te5': (0, 80), 'Pb2Bi2Te5': (-175, 0), 'Sn2Sb2Te5': (0, 60), 'Sn2Bi2Te5': (0, 90)}
    exp_rows, calc_rows, rows = [], [], []
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        h = HALL[m]; N = h['N']; data = EXP[m]; T, f = pbe_dense(m)
        v = np.array([f(N, i) for i in range(len(T))]); at = lambda t: float(np.interp(t, T, v))
        Te, Se = zip(*data); exp_rows += [[m, r(t, 2), r(sg, 1), r(x, 3)] for t, sg, x in read_exp(m)]
        ax.plot(Te, Se, 'o', ms=4.5, mfc='white', mec='k', mew=0.9, ls='', zorder=4)
        ax.plot(T, v, color=C_PBE, lw=1.4, zorder=3)
        calc_rows += [[m, f'{N:.3e}', int(t), r(x, 2)] for t, x in zip(T, v)]
        dev = [at(t) - S for t, S in data]; rmsv = float(np.sqrt(np.mean(np.square(dev))))
        rows.append([m, 'p' if N > 0 else 'n', f'{N:.3e}', 'yes' if h['flagged'] else 'no', r(data[0][0], 1), r(data[0][1], 1), r(at(data[0][0]), 1),
                     r(data[-1][0], 1), r(data[-1][1], 1), r(at(data[-1][0]), 1), r(rmsv, 1), r(max(dev, key=abs), 1)])
        ax.set_xlim(280, 660); ax.set_ylim(*ylim[m])
        ax.set_title(f'{TEX[m]} ({"n" if N < 0 else "p"}-type)', fontsize=8, pad=4)
        ax.text(0.04, 0.95, f'Hall: {"$p$" if N > 0 else "$n$"} = {sci(abs(N))}' + (' (flagged)' if h['flagged'] else ''),
                transform=ax.transAxes, ha='left', va='top', fontsize=6.5, color='#555555')
        ax.text(0.04, 0.865, f'rms deviation {rmsv:.0f} μV K$^{{-1}}$', transform=ax.transAxes, ha='left', va='top', fontsize=6.5, color=C_PBE)
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$S$ (μV K$^{-1}$)')
    for ax in axs[1]: ax.set_xlabel('$T$ (K)')
    hh = [plt.Line2D([], [], ls='', marker='o', ms=4.5, mfc='white', mec='k', mew=0.9), plt.Line2D([], [], color=C_PBE, lw=1.4)]
    fig.legend(hh, ['Experiment', 'PBE+SOC + BoltzTraP2 at the Hall concentration (no adjustable parameter)'],
               loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.935), columnspacing=1.4, handlelength=1.8)
    d = os.path.join(DATA, 'FigS3_S_vs_T') if DATA else None
    if d and os.path.isdir(d): shutil.rmtree(d)          # 旧版文件一并清掉
    dump('FigS3_S_vs_T', 'experiment_T_sigma_S.csv', ['material', 'T_K', 'sigma_exp_S/m', 'S_exp_uV/K'], exp_rows)
    dump('FigS3_S_vs_T', 'hall_room_temperature.csv', ['material', 'carrier_conc_signed_cm-3 (+ holes, - electrons)', 'Hall_mobility_cm2/Vs', 'flagged_by_instrument'],
         [[k, f'{x["N"]:.3e}', x['mu'], 'yes' if x['flagged'] else 'no'] for k, x in HALL.items()])
    dump('FigS3_S_vs_T', 'calculated_S_T_at_Hall_conc.csv',
         ['material', 'Hall_conc_signed_cm-3', 'T_K', 'S_PBE+SOC_BoltzTraP2_uV/K (28x28x7/6 dense grid, no scissor)'], calc_rows)
    dump('FigS3_S_vs_T', 'summary_at_Hall_conc.csv', ['material', 'type', 'Hall_conc_signed_cm-3', 'Hall_flagged', 'T_low_K', 'S_exp_low_uV/K', 'S_PBE_low_uV/K',
                                                      'T_high_K', 'S_exp_high_uV/K', 'S_PBE_high_uV/K', 'rms_deviation_uV/K', 'max_deviation_calc-exp_uV/K'], rows)
    save(fig, out, 'FigS3_S_vs_T_Hall')
    return rows


# ------------------------------------------------------------------ Fig. 10, 11
NAMES = {'Pb_on_Vac': r'Pb$_{\mathrm{Vac}}$', 'Sb_on_Vac': r'Sb$_{\mathrm{Vac}}$', 'Te_on_Vac': r'Te$_{\mathrm{Vac}}$',
         'Sb_Pb': r'Sb$_{\mathrm{Pb}}$', 'Pb_Sb': r'Pb$_{\mathrm{Sb}}$', 'Te_Pb': r'Te$_{\mathrm{Pb}}$', 'Te_Sb': r'Te$_{\mathrm{Sb}}$',
         'Sb_Te': r'Sb$_{\mathrm{Te}}$', 'Pb_Te': r'Pb$_{\mathrm{Te}}$', 'V_Te': r'$V_{\mathrm{Te}}$', 'V_Pb': r'$V_{\mathrm{Pb}}$', 'V_Sb': r'$V_{\mathrm{Sb}}$'}
ROWS = [(t, q) for t, q in rf.DONORS] + [(t, q.replace('-', '−')) for t, q in rf.ACCEPT]
DEFDIR = os.path.join(ROOT, '13_DVO缺陷_VASP流程与分析', '结果_第一批', 'analysis')


def ef_sites():
    """逐位点形成能：[(相, 缺陷类型, 名义电荷, 位点, E_f 富 Te, E_f 贫 Te)]。"""
    q = dict(rf.DONORS + rf.ACCEPT)
    out = []
    for row in csv.DictReader(open(os.path.join(DEFDIR, '形成能_逐位点.csv'), encoding='utf-8-sig')):
        out.append(('O' if row['parent_id'] == 'O_ordered' else 'D', row['defect_type'], q[row['defect_type']], row['defect_id'],
                    float(row['Ef0_Te-rich']), float(row['Ef0_Te-poor'])))
    return out


def fig10(out):
    d = rf.read_ef()
    fig, axs = plt.subplots(1, 2, figsize=(W2, 3.3), sharey=True, gridspec_kw=dict(wspace=0.06))
    for n, (ax, c) in enumerate(zip(axs, ('Te-rich', 'Te-poor'))):
        ax.axvline(0, color='k', lw=0.5)
        ax.axhspan(len(rf.DONORS) - 0.5, len(ROWS) - 0.5, color='#F2F2F2', zorder=0, lw=0)
        for j, (t, q) in enumerate(ROWS):
            y = len(ROWS) - 1 - j
            for fam, col, dy, lab in (('O', C_O, 0.16, 'Ordered (O)'), ('D', C_D, -0.16, 'Disordered (D)')):
                v = d.get((c, fam, t), [])
                ax.plot(v, [y + dy] * len(v), 'o', ms=3.6, color=col, mec='white', mew=0.4, ls='', label=lab if j == 0 else None, zorder=3)
        ax.set_xlim(-0.45, 3.05); ax.set_xlabel(r'Neutral formation energy (eV)')
        ax.text(0.98, 0.985, {'Te-rich': r'Te-rich ($\Delta\mu_{\mathrm{Te}}$ = 0)', 'Te-poor': r'Te-poor ($\Delta\mu_{\mathrm{Sb}}$ = 0)'}[c],
                transform=ax.transAxes, ha='right', va='top', fontsize=7)
        ax.tick_params(axis='y', length=0, right=False); panel(ax, f'({"ab"[n]})')
    axs[0].set_yticks(range(len(ROWS)), [f'{NAMES[t]} ({q})' for t, q in ROWS[::-1]])
    axs[0].text(2.95, len(ROWS) - 1.05, 'donors', ha='right', fontsize=6.5, style='italic')
    axs[0].text(2.95, len(rf.ACCEPT) - 1.05, 'acceptors', ha='right', fontsize=6.5, style='italic')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=2, bbox_to_anchor=(0.55, 1.04))
    order = {t: i for i, (t, _) in enumerate(ROWS)}
    dump('Fig10_formation_energies', 'formation_energy_per_site.csv',
         ['phase', 'defect_type', 'nominal_charge', 'site', 'Ef_Te-rich_eV', 'Ef_Te-poor_eV'],
         [[a, b, c, s, r(x, 4), r(y, 4)] for a, b, c, s, x, y in sorted(ef_sites(), key=lambda t: (order[t[1]], t[0] != 'O', t[3]))])
    mu = json.load(open(os.path.join(DEFDIR, 'DVO_结论摘要.json')))['chemical_potentials']['conditions']
    dump('Fig10_formation_energies', 'chemical_potentials.csv', ['condition', 'mu_Pb_eV', 'mu_Sb_eV', 'mu_Te_eV', 'dmu_Pb_eV', 'dmu_Sb_eV', 'dmu_Te_eV'],
         [[c] + [r(mu[c][e], 4) for e in ('Pb', 'Sb', 'Te')] + [r(mu[c]['delta'][e], 4) for e in ('Pb', 'Sb', 'Te')] for c in ('Te-rich', 'Te-poor')])
    save(fig, out, 'Fig10_defect_formation_energies')


def fig11(out):
    sites = ef_sites(); rows, med = [], []
    fig, ax = plt.subplots(figsize=(W1, 3.0))
    ax.axvline(0, color='k', lw=0.5)
    ax.axhspan(len(rf.DONORS) - 0.5, len(ROWS) - 0.5, color='#F2F2F2', zorder=0, lw=0)
    for j, (t, q) in enumerate(ROWS):
        y = len(ROWS) - 1 - j
        o = min(s[4] for s in sites if s[0] == 'O' and s[1] == t)
        ds = sorted((s for s in sites if s[0] == 'D' and s[1] == t), key=lambda s: s[3])
        sh = o - np.array([s[4] for s in ds])
        rows += [[t, s[2], s[3], r(o, 4), r(s[4], 4), r(v, 4)] for s, v in zip(ds, sh)]
        med.append([t, ds[0][2], len(ds), r(np.median(sh), 3)])
        ax.barh(y, np.median(sh), height=0.55, color=C_D, alpha=0.35, lw=0)
        ax.plot(sh, [y] * len(sh), 'o', ms=3.2, color=C_D, mec='white', mew=0.4, zorder=3)
    ax.set_yticks(range(len(ROWS)), [f'{NAMES[t]} ({q})' for t, q in ROWS[::-1]]); ax.tick_params(axis='y', length=0, right=False)
    ax.set_xlim(-0.8, 2.0); ax.set_xlabel(r'$\Delta E_f = E_f^{\mathrm{O}} - E_f^{\mathrm{D}}$ (eV)')
    ax.text(-0.76, len(ROWS) - 0.95, 'donors', ha='left', va='center', fontsize=6.5, style='italic')
    ax.text(-0.76, len(rf.ACCEPT) - 0.95, 'acceptors', ha='left', va='center', fontsize=6.5, style='italic')
    dump('Fig11_ordering', 'delta_Ef_per_site.csv', ['defect_type', 'nominal_charge', 'D_site', 'Ef_O_lowest_site_eV', 'Ef_D_eV', 'dEf_O-D_eV'], rows)
    dump('Fig11_ordering', 'delta_Ef_median.csv', ['defect_type', 'nominal_charge', 'n_D_sites', 'median_dEf_O-D_eV'], med)
    save(fig, out, 'Fig11_ordering_induced_change')


# ------------------------------------------------------------------ Fig. S1
def mirror_map(m):
    """c4_mirror_v3 网格上最高价带（Kramers 对平均）的能量，相对真实价带顶（meV）。返回 x（= −u，指向 M）、w、E[x, w]，
    用 (u, w) ≡ (−u, −w) 补全 x < 0 的一半。"""
    import c4_mirror_scan as c4
    d = os.path.join(ROOT, '22_能带基准_1-0', m, 'c4_mirror_v3'); meta = json.load(open(os.path.join(d, 'c4_grid.json')))
    E, _ = c4.read_parts(d, [os.path.join(d, f'b0{i}') for i in range(1, 7)], meta); nel = meta['nelect']
    us, ws = np.array(meta['us']), np.array(meta['ws']); gk = np.array(meta['grid_k_index'])
    G = 1000 * (E[gk, nel - 2:nel].mean(-1) - mv.true_load(m)['vbm'])          # [u, w]
    x = -us; jneg = [int(np.argmin(np.abs(((-w - ws) + 0.5) % 1 - 0.5))) for w in ws]   # −w 在网格中的序号
    X = np.concatenate([-x[:0:-1], x]); GX = np.concatenate([G[:0:-1][:, jneg], G])
    return X, ws, GX


def figS1(out):
    from scipy.interpolate import RectBivariateSpline
    fig = plt.figure(figsize=(W2, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.42, wspace=0.22)
    cmap = plt.get_cmap('viridis'); lo = -250
    cuts = {}; rows = []
    for k, m in enumerate(('Sn2Sb2Te5', 'Sn2Bi2Te5')):
        X, ws, GX = mirror_map(m)
        we = np.concatenate([ws[-2:] - 1, ws, ws[:3] + 1]); Ge = np.concatenate([GX[:, -2:], GX, GX[:, :3]], axis=1)   # w 方向周期延拓
        S = RectBivariateSpline(X, we, Ge, kx=3, ky=3, s=0)
        xf, wf = np.linspace(X.min(), X.max(), 241), np.linspace(-0.5, 0.5, 201); F = S(xf, wf)
        ax = fig.add_subplot(gs[0, k])
        im = ax.pcolormesh(xf, wf, np.clip(F.T, lo, 0), cmap=cmap, vmin=lo, vmax=0, shading='auto', rasterized=True)
        ax.contour(xf, wf, F.T, levels=np.arange(-300, 1, 25), colors='white', linewidths=0.35, alpha=0.7)
        for w0 in (-0.5, 0.0, 0.5): ax.axhline(w0, color='w', lw=1.6, ls=(0, (1, 1)))
        ax.axvline(0, color='w', lw=1.0, ls=(0, (1, 1)))
        ax.text(0.39, 0.02, 'Γ–M (path)', color='w', fontsize=6.5, ha='right', va='bottom')
        ax.text(0.39, 0.48, 'A–L (path)', color='w', fontsize=6.5, ha='right', va='top')
        ax.text(0.01, -0.33, 'Γ–A', color='w', fontsize=6.5, ha='left', rotation=90)
        for nm, (u0, w0) in PATH_LP.items():                         # 路径上的谷（空心）
            for xx, ww in ((-u0, w0), (u0, -w0)):
                for wv in ((ww, ww - 1) if np.isclose(abs(ww), 0.5) else (ww,)):
                    if abs(wv) <= 0.5 + 1e-9: ax.plot(xx, wv, MK[nm], ms=5 * MS.get(MK[nm], 1), mfc='none', mec='k', mew=0.9, zorder=4, clip_on=False)
        for t in true_valleys(m)['valleys']:                             # 三维搜索的真实极大（实心）
            u0, w0 = t['kfrac'][0], t['kfrac'][2]
            for xx, ww in ((-u0, w0), (u0, -w0)):
                ax.plot(xx, ww, MK[t['name']], ms=5.2 * MS.get(MK[t['name']], 1), color='#FF5A5F', mec='k', mew=0.5, zorder=5)
            right = -u0 > 0.28
            ax.annotate(f"{MATH.get(t['name'], t['name'])} {t['E']:.0f}".replace('-', '−') + ('*' if t['grid'] else ''), (-u0, w0),
                        (-6 if right else 6, -11 if (right or w0 > 0.35) else 5), textcoords='offset points', fontsize=6.3, color='w', fontweight='bold',
                        ha='right' if right else 'left')
            rows.append([m, t['name'], round(-u0, 4), round(w0, 4), round(t['E'], 1)])
            if t['name'] in ('Σ′', 'R₂′'): cuts[(m, t['name'])] = (-u0, w0, t['E'], S(np.array([-u0]), wf)[0])
        ax.set_xlim(X.min(), X.max()); ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel(r'in-plane $k$ along Γ→M ($-u$, units of $b_1$)'); ax.set_ylabel(r'$k_z$ ($w$, units of $b_3$)' if k == 0 else '')
        ax.set_title(f'{TEX[m]}: top valence band in the mirror plane', fontsize=8, pad=3); panel(ax, f'({"ab"[k]})')
        for nm, (u0, w0) in PATH_LP.items():
            cuts[(m, nm)] = S(np.array([-u0]), np.array([w0]))[0, 0]          # 路径上该谷的能量（同一样条，作对照）
        dump('FigS1_mirror_plane', f'{m}_mirror_plane_grid_meV.csv', ['x=-u'] + [f'w={w:.4f}' for w in ws], [[round(xx, 4)] + [round(v, 1) for v in row] for xx, row in zip(X, GX)])
    cb = fig.colorbar(im, ax=fig.axes[:2], fraction=0.025, pad=0.02, extend='min'); cb.set_label(r'$E - E_\mathrm{VBM}$ (meV)')
    # (c)(d) 沿 c 的切线，经过三维搜索的真实极大；路径只在 w = 0（Γ–M）和 w = ±½（A–L）两处经过这条线附近
    for k, (tn, pn) in enumerate((('Σ′', 'Σ'), ('R₂′', 'R₂'))):
        ax = fig.add_subplot(gs[1, k])
        for m, col in (('Sn2Sb2Te5', C_B), ('Sn2Bi2Te5', C_A)):
            x0, w0, E0, cut = cuts[(m, tn)]
            ax.plot(wf, cut, color=col, lw=1.3, label=f'{TEX[m]} ($-u$ = {x0:.3f})')
            ax.plot(w0, E0, MK[tn], ms=5.5 * MS.get(MK[tn], 1), color=col, mec='k', mew=0.5, zorder=5)
            wp = PATH_LP[pn][1]; Ep = cuts[(m, pn)]
            for ww in ((wp,) if wp == 0 else (-0.5, 0.5)):
                ax.plot(ww, Ep, MK[pn], ms=5.5 * MS.get(MK[pn], 1), mfc='white', mec=col, mew=1.0, zorder=5, clip_on=False)
            ax.annotate('', (w0, E0), (wp if wp == 0 else 0.5, Ep), arrowprops=dict(arrowstyle='->', lw=0.7, color=col, shrinkA=4, shrinkB=4))
            ax.text(w0 + 0.04 if wp == 0 else 0.03, (E0 + Ep) / 2 if wp == 0 else E0 - 12, f'+{E0 - Ep:.0f} meV', color=col, fontsize=6.5,
                    ha='left', va='center')
            rows += [[m, f'cut through {tn}', round(x0, 4), round(w, 4), round(e, 1)] for w, e in zip(wf[::5], cut[::5])]
        for ww in ((0.0,) if pn == 'Σ' else (-0.5, 0.5)): ax.axvline(ww, color='#999999', lw=0.8, ls=DASH)
        ax.text(0.0 if pn == 'Σ' else 0.5, 1.01, 'Γ–M plane (path)' if pn == 'Σ' else 'A–L plane (path)', transform=ax.get_xaxis_transform(),
                ha='center' if pn == 'Σ' else 'right', va='bottom', fontsize=6.5, color='#555555')
        ax.set_xlim(-0.5, 0.5); ax.set_ylim(-300, 15); ax.set_xlabel(r'$k_z$ ($w$, units of $b_3$)')
        ax.set_ylabel(r'$E - E_\mathrm{VBM}$ (meV)' if k == 0 else '')
        ax.set_title(f'Cut along $k_z$ through {MATH.get(tn, tn)} (filled) vs path {MATH.get(pn, pn)} (open)', fontsize=8, pad=11)
        ax.legend(loc='upper left' if pn == 'Σ' else 'center', fontsize=6.3, ncol=1); panel(ax, f'({"cd"[k]})', y=1.06)
    dump('FigS1_mirror_plane', 'valleys_and_cuts.csv', ['material', 'valley_or_cut', 'x=-u', 'w', 'E_minus_VBM_meV'], rows)
    save(fig, out, 'FigS1_mirror_plane_maps')


def bt2_signed(m, sign):
    """BoltzTraP2（PBE+SOC，全区，CRTA，面内）的 S(N, T)：sign = +1 空穴、−1 电子，N 取绝对值的对数插值，T 方向二次插值。"""
    d = json.load(open(os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'ke-dft-cpu', 'result', 'step8.1_boltztrap', 'boltztrap_crta.json')))
    T = np.array(d['temperatures_K']); c = np.array(d['carrier_conc_cm-3'])
    S = (np.array(d['seebeck_xx_V_per_K']) + np.array(d['seebeck_yy_V_per_K'])) / 2 * 1e6

    def f(N, t):
        v = []
        for it in range(len(T)):
            ok = sign * c[it] > 0; x = np.log10(sign * c[it][ok]); o = np.argsort(x)
            v.append(np.interp(np.log10(abs(N)), x[o], S[it][ok][o]))
        return np.polyval(np.polyfit(T, v, 2), t)
    return f


def main():
    global DATA
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.join(ROOT, '21_阶段报告与交接', '论文图_20261028'))
    ap.add_argument('--data'); ap.add_argument('--only', nargs='+'); a = ap.parse_args()
    DATA = a.data or os.path.join(a.out, 'source_data')
    os.makedirs(a.out, exist_ok=True); style()
    todo = a.only or [str(i) for i in range(1, 12)] + ['S1', 'S3']
    info = {}
    for i in todo:
        res = {'1': fig1, '2': fig2, '3': fig3, '4': fig4, '5': fig5, '6': fig6, '7': fig7, '8': fig8, '9': fig9, '10': fig10, '11': fig11, 'S1': figS1, 'S3': figS3}[i](a.out)
        if res is not None: info[i] = res
    print(info); print('→', a.out, '；源数据 →', DATA)


if __name__ == '__main__':
    main()
