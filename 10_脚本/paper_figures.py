"""Journal-style figures for the manuscript (English labels; 600-dpi PNG + vector PDF), numbered as in 论文报告_20261028.

Fig. 1  Structures: O/D 144-atom supercells (a, b) and side views of the user's layered phase vs model O (c, d)
                                                                               补算实施方案_20261007/.../parents/<id>/site_map.json;
                                                                               A2B2Te5_SOC_dense/materials/Pb2Sb2Te5/POSCAR; 13_…/energies.csv
Fig. 2  PBE+SOC vs HSE06+SOC band structures (c2 path)                       22_能带基准_1-0/<mat>/c2_{pbe,hse}/path_eigen.txt
Fig. 3  HSE06+SOC bands weighted by cation-s (lone-pair) character             22_能带基准_1-0/<mat>/c2_hse/PROCAR_OPT.gz
Fig. 4  HSE06+SOC bands: which cation (A = Pb/Sn or B = Sb/Bi) carries the s weight   same
Fig. 5  Sn2Sb2Te5 vs Sn2Bi2Te5: Sn-5s and B-site p (Sb 5p / Bi 6p) projections near the VBM, and level diagram   same
Fig. 6  Valence pockets in the Brillouin zone (Sn2Sb2Te5 vs Sn2Bi2Te5) and the lone-pair / Γ–A crossover energy, HSE only
                                                                               same + A2B2Te5_SOC_dense/materials/<mat>/POSCAR (spglib)
Fig. 7  Energy of the highest lone-pair valence valley, PBE vs HSE            same
Fig. 8  Seebeck coefficient at 3e19 cm-3: PBE (BoltzTraP2) vs HSE estimate     10_脚本/multivalley_seebeck.py
        (+ per-valley conductivity shares and cation-s character for Table 6)
Fig. 9  S(T) of all four compounds: experiment vs calculation (HSE and PBE)    experimental S(T) from the user;
                                                                               A2B2Te5_SOC_dense/materials/<mat>/.../boltztrap_crta.json
Fig. 10 Neutral formation energies of 12 defect types in O and D               13_DVO缺陷_VASP流程与分析/结果_第一批/analysis/形成能_逐位点.csv
Fig. 11 Ordering-induced change of the formation energy, O − D                  same

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


def kpath():
    return cv.read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))


def bands(m, side):
    E = cv.read_path_eigen(os.path.join(rf.base(m), side, 'path_eigen.txt'))
    nel = NEL[m]; vbm = np.nanmax(E[:, nel - 1]); gap = np.nanmin(E[:, nel]) - vbm
    return E - vbm, gap


# 路径上离开对称点的位置按所在的高对称线命名（六方布里渊区标准记号，Bradley & Cracknell）
SEGS = [(0, 33, 'Σ'), (33, 52, 'T'), (52, 90, 'Λ'), (90, 97, 'Δ'), (97, 130, 'R'), (130, 149, 'S'), (149, 187, 'Q')]
MATH = {'R₁': r'R$_1$', 'R₂': r'R$_2$'}          # 图中用 mathtext 写下标


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
        x, _ = rf.path_x(m); nel = NEL[m]; g = {}; gam = []
        for side, nm, col, ls, lw, z in (('c2_pbe', 'PBE', C_PBE, DASH, 0.9, 2), ('c2_hse', 'HSE', C_HSE, '-', 1.1, 3)):
            E, g[nm] = bands(m, side); bs = range(nel - 10, nel + 8)
            gam += [r(1000 * E[0, nel - 1], 1), r(1000 * E[0, nel], 1), r(1000 * (E[0, nel] - E[0, nel - 1]), 1)]   # Γ 处最高价带、最低导带、两者之差
            for b in bs:
                ax.plot(x[:KEND + 1], E[:KEND + 1, b], color=col, ls=ls, lw=lw, zorder=z)
            dump('Fig2_bands_PBE_vs_HSE', f'{m}_{nm.replace("HSE", "HSE06")}+SOC.csv',
                 ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label'] + [bandname(b, nel) + '_E-EVBM_eV' for b in bs],
                 [[i] + [r(v, 6) for v in K[i]] + [r(x[i], 5), LAB.get(i, '')] + [r(E[i, b], 5) for b in bs] for i in range(KEND + 1)])
        gaps.append([m, r(1000 * g['PBE'], 1), r(1000 * g['HSE'], 1)] + gam)
        ax.axhline(0, color='k', lw=0.5, ls=':', zorder=1)
        kticks(ax, x); ax.set_ylim(-0.66, 0.66); ax.set_yticks(np.arange(-0.6, 0.61, 0.2))
        ax.set_title(f'{TEX[m]}     $E_g$: PBE {1000 * g["PBE"]:.0f} meV, HSE {1000 * g["HSE"]:.0f} meV', fontsize=8, pad=4)
        panel(ax, f'({"abcd"[n]})')
    dump('Fig2_bands_PBE_vs_HSE', 'band_gaps_on_path.csv',
         ['material', 'Eg_PBE+SOC_meV', 'Eg_HSE06+SOC_meV'] +
         [f'{s}_Gamma_{q}_meV' for s in ('PBE+SOC', 'HSE06+SOC') for q in ('VB', 'CB', 'gap')], gaps)
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    h = [plt.Line2D([], [], color=C_HSE, lw=1.1), plt.Line2D([], [], color=C_PBE, lw=0.9, ls=DASH)]
    fig.legend(h, ['HSE06+SOC', 'PBE+SOC'], loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0.905))
    save(fig, out, 'Fig2_bands_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 3
def fig3(out):
    cmap = LinearSegmentedColormap.from_list('s', ['#DCE8F4', '#9CC3E4', '#4A90C8', '#1F5FA0', '#08306B'])
    norm = Normalize(0, 0.30); K = kpath()
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.0), sharey=True, gridspec_kw=dict(hspace=0.34, wspace=0.06))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        x, _ = rf.path_x(m); nel = NEL[m]
        E, gap = bands(m, 'c2_hse')
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
        dump('Fig3_HSE_bands_cation_s', f'{m}_HSE06+SOC_cation_s.csv',
             ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV', 'cation_s_percent'], rows)
        ax.axhline(0, color='k', lw=0.5, ls=':'); ax.axhline(gap, color='k', lw=0.5, ls=':')
        kticks(ax, x); ax.set_ylim(-0.5, 0.5)
        ax.set_title(f'{TEX[m]}     $E_g$ = {1000 * gap:.0f} meV', fontsize=8, pad=4)
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
            ax.add_collection(LineCollection(segs, colors=cols, linewidths=wid, capstyle='round'))
            rows += [[i] + [r(v, 6) for v in K[i]] + [r(x[i], 5), LAB.get(i, ''), bandname(b, nel), r(E[i, b], 5),
                     r(100 * fa[i, b], 2), r(100 * fb[i, b], 2), r(100 * ft[i, b], 2)] for i in range(KEND + 1)]
        dump('Fig4_cation_s_by_atom', f'{m}_HSE06+SOC_s_by_atom.csv',
             ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV',
              f'A_{A_}_s_percent', f'B_{B_}_s_percent', 'Te_s_percent'], rows)
        # 关键态：价带顶、导带底、最高孤对电子价带谷（Kramers 对平均）
        iv = int(np.nanargmax(E[:KEND + 1, nel - 1])); ic = int(np.nanargmin(E[:KEND + 1, nel])); lp = rf.lone_pair_offset(m, 'c2_hse')
        for st, i, bb in (('VBM', iv, [nel - 2, nel - 1]), ('lone-pair valence valley', lp[1], [nel - 2, nel - 1]), ('CBM', ic, [nel, nel + 1])):
            key.append([m, st, i] + [r(v, 4) for v in K[i]] + [r(1000 * E[i, bb].mean(), 1), A_, r(100 * fa[i, bb].mean(), 1),
                       B_, r(100 * fb[i, bb].mean(), 1), r(100 * ft[i, bb].mean(), 1)])
        ax.axhline(0, color='k', lw=0.5, ls=':'); ax.axhline(gap, color='k', lw=0.5, ls=':')
        kticks(ax, x); ax.set_ylim(-0.5, 0.5)
        ax.set_title(f'{TEX[m]}     A = {A_}, B = {B_}', fontsize=8, pad=4)
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=SITE_CMAP), ax=axs, shrink=0.7, pad=0.015, aspect=30)
    cb.set_label('Share of cation-s on the A site'); cb.set_ticks([0, 0.5, 1], labels=['B only\n(Sb, Bi)', 'equal', 'A only\n(Pb, Sn)'])
    cb.ax.tick_params(direction='in', length=2, width=0.6, labelsize=6.5); cb.outline.set_linewidth(0.6)
    dump('Fig4_cation_s_by_atom', 'key_states.csv', ['material', 'state', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'E-EVBM_meV',
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
        E, gap = bands(m, 'c2_hse')
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        tot = np.maximum(Wp.sum((-1, -2)), 1e-9)
        sn_s = Wp[:, :, [7, 8], 0].sum(-1) / tot; b_s = Wp[:, :, [0, 1], 0].sum(-1) / tot
        b_p = Wp[:, :, [0, 1], 1:4].sum((-1, -2)) / tot; sn_p = Wp[:, :, [7, 8], 1:4].sum((-1, -2)) / tot
        te_p = Wp[:, :, 2:7, 1:4].sum((-1, -2)) / tot
        rows = []
        for row, (w, wmax, col, lab) in enumerate(((sn_s, 0.25, C_A, 'Sn 5s'), (b_p, 0.45, C_B, f'{B_} {"5p" if B_ == "Sb" else "6p"}'))):
            ax = fig.add_subplot(gs[row, c])
            for b in range(nel - 8, nel + 6):
                ax.plot(x[:KEND + 1], E[:KEND + 1, b], color='#C8C8C8', lw=0.5, zorder=1)
                segs = [[(x[i], E[i, b]), (x[i + 1], E[i + 1, b])] for i in range(KEND)]
                wid = [3.2 * np.clip((w[i, b] + w[i + 1, b]) / 2 / wmax, 0, 1) for i in range(KEND)]
                ax.add_collection(LineCollection(segs, colors=col, linewidths=wid, capstyle='round', zorder=2))
            ax.axhline(0, color='k', lw=0.5, ls=':')
            kticks(ax, x); ax.set_ylim(-0.42, 0.30)
            if c == 0: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
            else: ax.tick_params(labelleft=False)
            ax.set_title(f'{TEX[m]}: {lab}', fontsize=8, pad=4, color=col)
            panel(ax, f'({"abcd"[2 * row + c]})')
            for nm, i in SNVAL.items():                      # 标出四个谷：上排标孤对电子谷，下排标 Γ 单谷
                if (row == 0) == (nm != 'Γ'):
                    xi = x[90] if nm == 'Γ' else x[i]
                    ax.plot(xi, E[i, nel - 1] + 0.035, 'v', ms=3.2, color='k', zorder=3)
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
    # (e) 能级：各自以价带顶为零点
    ax = fig.add_subplot(gs[:, 2]); X = (0.0, 1.7); hw = 0.3
    for nm, col, lw in (('Γ', C_B, 2.0), ('Σ', C_A, 1.4), ('R₁', C_A, 1.4), ('R₂', C_A, 1.4)):
        y = [lv[(m, nm)] for m in mats]
        for k in range(2): ax.plot([X[k] - hw, X[k] + hw], [y[k]] * 2, color=col, lw=lw, solid_capstyle='butt')
        ax.plot([X[0] + hw, X[1] - hw], y, color=col, lw=0.6, ls=(0, (2, 1.5)))
    for k, m in enumerate(mats):                     # 标签放在各自能级上方
        B_ = m[3:5]; top = max(lv[(m, v)] for v in ('Σ', 'R₁', 'R₂'))
        ax.text(X[k], lv[(m, 'Γ')] + 2.5, f'Γ ({B_} p)', ha='center', va='bottom', fontsize=6.5, color=C_B)
        ax.text(X[k], top + 2.5 + (9 if k == 1 else 0), 'Sn 5s', ha='center', va='bottom', fontsize=6.5, color=C_A)
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xlim(-0.75, 2.45); ax.set_ylim(-165, 26)
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


def crossover():
    """转变能量（HSE06+SOC）：最高孤对电子谷（与图 7 同一判据）减 Γ–A 轴（Γ、Δ 线、A）上最高的价带态，四种材料。"""
    K = kpath(); rows = []
    for m in MATS:
        _, rots = bz_lattice(m)
        for side, fn in (('c2_hse', 'HSE06+SOC'),):
            v = top_vb(m, side); ia = max(AXIS, key=lambda i: v[i])
            lp = rf.lone_pair_offset(m, side); il = lp[1]
            rows.append(dict(m=m, fn=fn, EG=v[0], Ea=v[ia], ia=ia, na=len(star(K[ia], rots)),
                             El=v[il], il=il, nl=len(star(K[il], rots)), s=lp[2], d=v[il] - v[ia]))
    return rows


def pockets(m):
    """HSE06+SOC：Γ–A 轴上最高的价带态与三组孤对电子谷（Σ、R₁、R₂），展开成第一布里渊区内的全部等价口袋。"""
    B, rots = bz_lattice(m); K = kpath(); v = top_vb(m, 'c2_hse')
    ia = max(AXIS, key=lambda i: v[i])
    out = []
    for i, col in [(ia, C_B)] + [(i, C_A) for nm, i in SNVAL.items() if nm != 'Γ']:
        pts = star(K[i], rots)
        out.append(dict(name=kname(i), i=i, E=v[i], col=col, frac=pts, cart=[ws(q, B) for q in pts]))
    return B, out


def hexagon(B):
    K0 = np.array([-2 / 3, 1 / 3, 0]) @ B; a = np.arctan2(K0[1], K0[0]) + np.arange(7) * np.pi / 3
    return np.hypot(*K0[:2]) * np.cos(a), np.hypot(*K0[:2]) * np.sin(a)


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
            ax.plot(c[0], c[1], 'o', ms=5.6 if big else 3.4, mew=0.9, color=g['col'], mfc='white' if on_face(f) else g['col'],
                    alpha=1.0 if big else 0.45, zorder=3)
    Mq = max(M, key=lambda c: c[0] + 0.1 * c[1]); kv = max(zip(hx, hy), key=lambda t: (round(t[1], 3), t[0]))
    ax.text(Mq[0] + 0.04, Mq[1], 'M (L)', ha='left', va='center', fontsize=6.5)
    ax.text(kv[0] + 0.03, kv[1] + 0.02, 'K (H)', ha='left', va='bottom', fontsize=6.5)
    ax.annotate('Γ (A)', (0, 0), (-0.62, -0.50), fontsize=6.5, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.5, color='#555555', shrinkA=1, shrinkB=3))
    y = -0.66                                     # 右下角小表：位置名、相对价带顶的能量（meV）、等价口袋数
    for g in P:
        al = 1.0 if near(g['E']) else 0.45
        ax.plot(0.74, y, 'o', ms=3.4, mew=0.8, color=g['col'], mfc='white' if all(map(on_face, g['frac'])) else g['col'], alpha=al)
        ax.text(0.83, y, MATH.get(g['name'], g['name']), ha='left', va='center', fontsize=6.3)
        ax.text(1.66, y, f"{g['E']:.0f}".replace('-', '−') + f"  ×{len(g['frac'])}", ha='right', va='center', fontsize=6.3)
        y -= 0.115
    n = sum(len(g['frac']) for g in P if near(g['E']))
    ax.set_title(f"{TEX[m]}: {n} pocket{'s' if n > 1 else ''} within $2k_\\mathrm{{B}}T$", fontsize=8, pad=2)
    ax.set_aspect('equal'); ax.set_xlim(-1.0, 1.68); ax.set_ylim(-1.12, 0.98); ax.axis('off')
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
                ax.plot([c[0]], [c[1]], [z], 'o', ms=4.6 if big else 3.0, color=g['col'], mfc='white' if on_face(f) else g['col'],
                        mew=0.9, alpha=1.0 if big else 0.45)
    ax.text(-0.10, 0, h + 0.06, 'A', fontsize=6.5, ha='center', va='bottom'); ax.text(0.12, -0.10, -0.08, 'Γ', fontsize=6.5, ha='left', va='top')
    ax.view_init(elev=22, azim=-62); ax.set_box_aspect((1, 1, 2 * h / 1.9)); ax.set_axis_off()


def fig6(out):
    """(a)(b) Sn₂Sb₂Te₅、Sn₂Bi₂Te₅ 价带口袋沿 c 俯视，(c) Sn₂Bi₂Te₅ 三维透视；(d) 四种材料的转变能量。"""
    data = {m: pockets(m) for m in ('Sn2Sb2Te5', 'Sn2Bi2Te5')}
    fig = plt.figure(figsize=(W2, 4.25))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.05], height_ratios=[1, 0.92], wspace=0.04, hspace=0.16)
    for k, m in enumerate(data):
        bz_top(fig.add_subplot(gs[0, k]), m, *data[m], f'({"ab"[k]})')
    ax = fig.add_subplot(gs[0, 2], projection='3d'); zs = 4
    bz_3d(ax, 'Sn2Bi2Te5', *data['Sn2Bi2Te5'], zs=zs)
    ax.text2D(0.5, 0.97, f'{TEX["Sn2Bi2Te5"]}: 3D view ($k_z$ ×{zs})', transform=ax.transAxes, fontsize=8, ha='center', va='bottom')
    ax.text2D(0.04, 0.97, '(c)', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')
    # (d) 转变能量（HSE06+SOC）：正为孤对电子多谷在上（紫），负为 Γ–A 轴上的态在上（橙），与 (a)–(c) 同色
    rows = crossover(); by = {d['m']: d for d in rows}
    ax = fig.add_subplot(gs[1, :2])
    ax.axvspan(-KT300, KT300, color='#EEF3F8', zorder=0, lw=0); ax.axvline(0, color='k', lw=0.6, zorder=2)
    for j, m in enumerate(MATS[::-1]):
        d = by[m]; v = d['d']; col = C_A if v > 0 else C_B
        ax.barh(j, v, height=0.5, color=col, lw=0, zorder=1)
        lab = f'{v:+.0f}'.replace('-', '−') + f"  ({kname(d['il'])} − {kname(d['ia'])})"
        ax.text(v + (6 if v > 0 else -6), j, lab, ha='left' if v > 0 else 'right', va='center', fontsize=6.5)
    ax.set_yticks(range(4), [TEX[m] for m in MATS[::-1]]); ax.set_ylim(-0.6, 4.1); ax.tick_params(axis='y', length=0, right=False)
    ax.set_xlim(-330, 170)
    ax.set_xlabel(r'$E_{\mathrm{lone\ pair}} - E_{\Gamma\mathrm{-A\ axis}}$ (meV)')
    ax.text(-165, 3.62, 'Γ–A axis state on top (1–2 pockets)', ha='center', va='center', fontsize=6.5, color=C_B)
    ax.text(85, 3.62, 'lone-pair valleys on top' '\n' '(6 pockets each)', ha='center', va='center', fontsize=6.5, color=C_A)
    ax.text(0, -0.45, r'$\pm k_\mathrm{B}T$', ha='center', va='center', fontsize=6.3, bbox=dict(fc='#EEF3F8', ec='none', pad=0.6), zorder=3)
    ax.set_title('Crossover energy (HSE06+SOC)', fontsize=8, pad=4)
    panel(ax, '(d)')
    # 右下：(a)–(c) 的图例
    lg = fig.add_subplot(gs[1, 2]); lg.axis('off')
    h = [Line2D([], [], ls='', marker='o', ms=5, color=C_B, label='B-site p state on the Γ–A axis'),
         Line2D([], [], ls='', marker='o', ms=5, color=C_A, label='Sn 5s lone-pair valley'),
         Line2D([], [], ls='', marker='o', ms=5, color='#555555', label=r'filled: $k_z$ = 0 plane or inside'),
         Line2D([], [], ls='', marker='o', ms=5, color='#555555', mfc='white', label=r'open: $k_z$ = ½ zone face'),
         Line2D([], [], ls='', marker='o', ms=3.4, color=C_A, alpha=0.45, label=r'faded: more than $2k_\mathrm{B}T$ below VBM')]
    lg.legend(handles=h, loc='upper left', bbox_to_anchor=(0, 0.97), title='(a)–(c)', title_fontsize=7, handletextpad=0.3, labelspacing=0.6, borderaxespad=0)
    lg.text(0.02, 0.02, 'Labels: line names of the hexagonal zone\n(Σ: Γ–M, R: A–L, Δ: Γ–A); energies in meV\nrelative to the VBM; ×n = equivalent pockets.',
            transform=lg.transAxes, fontsize=6.3, va='bottom', ha='left', color='#333333')
    save(fig, out, 'Fig6_valence_pockets_crossover')
    prow = []
    for m, (B, P) in data.items():
        for g in P:
            for f, c in zip(g['frac'], g['cart']):
                prow.append([m, g['name'], g['i'], *[r(x, 4) for x in f], *[r(x, 4) for x in c],
                             'kz=1/2 zone face' if on_face(f) else ('kz=0' if np.isclose(f[2], 0) else 'Γ–A axis (Δ line)'),
                             r(g['E'], 1), len(g['frac'])])
    dump('Fig6_valence_pockets', 'pocket_positions.csv', ['material', 'valley', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3',
         'kx_1/A', 'ky_1/A', 'kz_1/A', 'plane', 'E-EVBM_meV (HSE06+SOC)', 'n_equivalent_pockets'], prow)
    dump('Fig6_valence_pockets', 'crossover_energy_HSE.csv',
         ['material', 'functional', 'E_Gamma_meV', 'E_GammaA_axis_max_meV', 'axis_position', 'axis_path_k_index', 'axis_n_pockets',
          'E_lone_pair_max_meV', 'lone_pair_valley', 'lone_pair_path_k_index', 'lone_pair_n_pockets', 'lone_pair_cation_s_percent',
          'E_lone_pair_minus_E_axis_meV', 'E_lone_pair_minus_E_Gamma_meV'],
         [[d['m'], d['fn'], r(d['EG'], 1), r(d['Ea'], 1), kname(d['ia']), d['ia'], d['na'], r(d['El'], 1), kname(d['il']), d['il'], d['nl'],
           r(100 * d['s'], 1), r(d['d'], 1), r(d['El'] - d['EG'], 1)] for d in rows])
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
        p, h = rf.lone_pair_offset(m, 'c2_pbe'), rf.lone_pair_offset(m, 'c2_hse')
        for nm, v in (('PBE+SOC', p), ('HSE06+SOC', h)):
            rows.append([m, nm, r(v[0], 1), v[1]] + [r(c, 4) for c in K[v[1]]] + [r(100 * v[2], 1)])
        ax.plot([p[0], h[0]], [j, j], color='#BDBDBD', lw=1.2, zorder=1)
        ax.plot(p[0], j + 0.08, 'o', ms=4.5, color=C_PBE, mec='white', mew=0.6, zorder=3, label='PBE+SOC' if j == 0 else None)
        ax.plot(h[0], j - 0.08, 'o', ms=4.5, color=C_HSE, mec='white', mew=0.6, zorder=3, label='HSE06+SOC' if j == 0 else None)
        ax.text(h[0] - 7, j - 0.08, f'{h[0]:.0f}'.replace('-', '−'), ha='right', va='center', fontsize=6.3)
    dump('Fig7_lone_pair_valley', 'highest_lone_pair_valley.csv',
         ['material', 'functional', 'E-EVBM_meV', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'cation_s_percent'], rows[::-1])
    ax.set_yticks(range(4), [TEX[m] for m in MATS[::-1]]); ax.set_ylim(-0.6, 3.85)
    ax.set_xlim(-275, 15); ax.set_xlabel('Highest lone-pair valence valley, $E - E_{\\mathrm{VBM}}$ (meV)')
    ax.tick_params(axis='y', length=0, right=False); ax.legend(loc='lower left', handletextpad=0.2, borderaxespad=0.2)
    save(fig, out, 'Fig7_lone_pair_valley_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 8
def fig8(out):
    est = rf.hse_estimates(); prov = mv.c3_provisional('Sn2Bi2Te5')
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.5), gridspec_kw=dict(wspace=0.18))
    xs = np.arange(4); w = 0.36
    for k, ax in enumerate(axs):
        pb = [abs(rf.BTP[m][k]) for m in MATS]; hs = [abs(est[m][k]) for m in MATS]
        ax.bar(xs - w / 2, pb, w * 0.95, color=C_PBE, label='PBE+SOC (BoltzTraP2)', lw=0)
        bars = ax.bar(xs + w / 2, hs, w * 0.95, color=C_HSE, label='HSE06+SOC (this work)', lw=0)
        for i, v in enumerate(hs):
            ax.text(xs[i] + w / 2, v + 4, f'{v:.0f}' + ('*' if (k == 0 and i == 3 and prov) else ''), ha='center', fontsize=6.5)
        if k == 0 and prov: bars[3].set_hatch('/////'); bars[3].set_edgecolor('white')
        ax.set_xticks(xs, [TEX[m] for m in MATS]); ax.tick_params(axis='x', top=False)
        ax.set_ylim(0, 270); ax.set_ylabel(r'$|S|$ (μV K$^{-1}$)' if k else r'$S$ (μV K$^{-1}$)')
        ax.text(0.04, 0.97, ('p-type' if k == 0 else 'n-type') + r', $3\times10^{19}$ cm$^{-3}$, 300 K', transform=ax.transAxes, ha='left', va='top', fontsize=7)
        panel(ax, f'({"ab"[k]})')
    dump('Fig8_Seebeck', 'S_300K_3e19.csv', ['material', 'carrier_type', 'carrier_conc_cm-3', 'S_PBE+SOC_BoltzTraP2_uV/K',
                                             'S_HSE06+SOC_estimate_uV/K', 'provisional'],
         [[m, 'pn'[k], 3e19, rf.BTP[m][k], r(est[m][k], 1), 'yes' if (k == 0 and m == 'Sn2Bi2Te5' and prov) else 'no']
          for k in range(2) for m in MATS])
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.07), handlelength=1.2)
    save(fig, out, 'Fig8_Seebeck_PBE_vs_HSE')
    valley_shares()
    return est


LP_MIN = 0.15        # 孤对电子谷的判据：谷极值点处阳离子 s 成分 ≥ 15%（与图 4 相同）


def valley_index(name):
    """multivalley_seebeck 的谷名 → c2 路径上的 k 序号（取成分用）。"""
    if name.startswith('#'): return int(name[1:].split()[0])
    if name.startswith('Γ–A'): return 90 + round(float(name.split('±')[1]) * 14)
    return mv.POINTS[name][0][0]


def valley_shares():
    """表 4 的电导份额：多谷模型（HSE、CRTA、300 K、面内、3e19）中每个谷占总电导（含少数载流子）的比例，
    以及谷极值点处的阳离子 s 成分；p 型孤对电子谷 = 价带谷中阳离子 s ≥ 15% 者，n 型“Γ 导带”= Γ 或 Γ–A 驼峰处的导带谷。"""
    K = kpath(); rows, summary, scan = [], [], []
    for m in MATS:
        c3 = mv.c3_usable(m)[0]
        V, gap = mv.build(m, 'v', 300, c3); C, _ = mv.build(m, 'c', 300)
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fs = rf.cation_s(Wp); nel = NEL[m]
        nm = {'v': knames([valley_index(v['name']) for v in V]), 'c': knames([valley_index(v['name']) for v in C])}
        for N, lab in ((3e19, 'p'), (-3e19, 'n')):
            _, EF, parts = mv.solve(V, C, gap, N, 300)
            lp = 0.0
            for (side, name, w), v in zip(parts, V + C):
                i = valley_index(name); b = nel - 1 if side == 'v' else nel; cs = fs[i, b]
                flag = (side == 'v' and lab == 'p' and cs >= LP_MIN) or (side == 'c' and lab == 'n' and (name in ('Γ',) or name.startswith('Γ–A')))
                lp += w if flag else 0.0
                rows.append([m, lab, r(1000 * EF, 1), 'valence' if side == 'v' else 'conduction', nm[side][i], i]
                            + [r(x, 4) for x in K[i]] + [r(1000 * v['off'], 1), v['g'], r(100 * cs, 1), r(100 * w, 2), 'yes' if flag else 'no'])
            summary.append([m, lab, r(100 * lp, 1)])
        # p 型孤对电子谷份额随空穴浓度的变化（判据同上）
        lpv = [name for name, v in zip([x[1] for x in parts[:len(V)]], V) if fs[valley_index(name), nel - 1] >= LP_MIN]
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
# 实验 S(T)，从所给曲线图读取（µV/K）
# 实验 S(T)（µV/K），用户提供。Pb₂Sb₂Te₅ 为原始数据（2026-10-29，见 21_阶段报告与交接/实验数据/），其余三种从所给曲线图读取
EXP = {'Pb2Sb2Te5': [(301.58, 5.52), (319.32, 5.39), (370.16, 5.60), (418.62, 4.94), (468.56, 5.09), (519.25, 6.93), (568.17, 10.67), (619.36, 17.10)],
       'Pb2Bi2Te5': [(302, -89), (318, -89), (368, -101), (415, -116), (463, -127), (513, -133), (562, -132), (611, -127)],
       'Sn2Sb2Te5': [(300, 16), (320, 16), (370, 19), (418, 23), (467, 28), (518, 33), (567, 37), (618, 42)],
       'Sn2Bi2Te5': [(302, 30), (318, 31), (368, 39), (417, 46), (466, 55), (518, 62), (567, 66), (619, 67)]}


def model(m, side):
    V, gap = mv.build(m, 'v', 300, None, side); C, _ = mv.build(m, 'c', 300, side=side)
    return lambda N, T: mv.solve(V, C, gap, N, T, 'crta')[0]


def bt2(m):
    """PBE+SOC BoltzTraP2（全布里渊区，CRTA，面内 xx、yy 平均）的 S(p, T)，p > 0 为空穴。只算了 300、500、700 K，T 方向二次插值。"""
    d = json.load(open(os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'ke-dft-cpu', 'result', 'step8.1_boltztrap', 'boltztrap_crta.json')))
    T = np.array(d['temperatures_K']); c = np.array(d['carrier_conc_cm-3'])
    S = (np.array(d['seebeck_xx_V_per_K']) + np.array(d['seebeck_yy_V_per_K'])) / 2 * 1e6

    def f(p, t):
        v = []
        for it in range(len(T)):
            ok = c[it] > 0; o = np.argsort(c[it][ok])
            v.append(np.interp(np.log10(p), np.log10(c[it][ok][o]), S[it][ok][o]))
        return np.polyval(np.polyfit(T, v, 2), t)
    return f


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


def fig9(out):
    Ts = np.arange(300, 641, 10)
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.9), gridspec_kw=dict(hspace=0.40, wspace=0.24))
    titles = {'Pb2Sb2Te5': 'near-compensated', 'Pb2Bi2Te5': 'n-type', 'Sn2Sb2Te5': 'p-type', 'Sn2Bi2Te5': 'p-type'}
    ylim = {'Pb2Sb2Te5': (-140, 30), 'Pb2Bi2Te5': (-245, -58), 'Sn2Sb2Te5': (0, 57), 'Sn2Bi2Te5': (0, 90)}
    # 拟合出的载流子浓度标在各图空白处：(x, y, 对齐)，轴坐标
    note = {'Pb2Sb2Te5': (0.97, 0.52, 'right'), 'Pb2Bi2Te5': (0.97, 0.93, 'right'), 'Sn2Sb2Te5': (0.04, 0.93, 'left'), 'Sn2Bi2Te5': (0.04, 0.93, 'left')}
    exp_rows, calc_rows, fit_rows = [], [], []
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        Te, Se = zip(*EXP[m]); T0, S0 = Te[0], Se[0]
        exp_rows += [[m, t, s] for t, s in EXP[m]]
        curves = []
        if m == 'Pb2Sb2Te5':        # 近补偿：用 301.6 K 的 S 定净掺杂（N 在 ±10¹⁹ 内单调，二分）；比较两种带隙
            for side, nm in (('c2_hse', 'HSE06+SOC'), ('c2_pbe', 'PBE+SOC')):
                g = model(m, side); N = solve_N(g, S0, T0, -1e19, 1e19)
                gap = mv.build(m, 'v', 300, None, side)[1]
                curves.append((nm, 'path multivalley model', g, N, f'{nm[:3]} ($E_g$ = {1000 * gap:.0f} meV): $p-n$ = {"+" if N > 0 else "\u2212"}{sci(N)}'))
        elif m == 'Pb2Bi2Te5':      # n 型：302 K 的实验值定电子浓度
            for side, nm in (('c2_hse', 'HSE06+SOC'), ('c2_pbe', 'PBE+SOC')):
                g = model(m, side); N = solve_log(g, EXP[m], 17, 20, -1)
                curves.append((nm, 'path multivalley model', g, N, f'{nm[:3]}: $n$ = {sci(N)}'))
        else:                       # Sn 系重掺杂 p 型：BoltzTraP2 全区；HSE = BoltzTraP2 × 模型 HSE/PBE 比值（1.7 节第 5 步）
            b = bt2(m); gh, gp = model(m, 'c2_hse'), model(m, 'c2_pbe')
            hse = lambda p, t, b=b, gh=gh, gp=gp: b(p, t) * gh(p, t) / gp(p, t)
            for nm, meth, g in (('HSE06+SOC', 'BoltzTraP2 (PBE) x model HSE/PBE ratio', hse), ('PBE+SOC', 'BoltzTraP2 (PBE, full BZ)', b)):
                N = solve_log(g, EXP[m], 20, 21.8)
                curves.append((nm, meth, g, N, f'{nm[:3]}: $p$ = {sci(N)}'))
        ax.plot(Te, Se, 'o', ms=4, mfc='white', mec='k', mew=0.8, ls='', zorder=3)
        for nm, meth, g, N, lab in curves:
            v = [g(N, t) for t in Ts]
            ax.plot(Ts, v, color=C_HSE if nm[0] == 'H' else C_PBE, ls='-' if nm[0] == 'H' else DASH)
            calc_rows += [[m, nm, meth, f'{N:.3e}', int(t), r(s, 2)] for t, s in zip(Ts, v)]
            dev = [g(N, t) - s for t, s in EXP[m]]
            fit_rows.append([m, nm, meth, f'{N:.3e}', T0, S0, r(max(dev, key=abs), 1)])
        if m == 'Pb2Sb2Te5': ax.axhline(0, color='k', lw=0.5, ls=':')
        ax.set_xlim(280, 640); ax.set_ylim(*ylim[m])
        ax.set_title(f'{TEX[m]} ({titles[m]})', fontsize=8, pad=4)
        x0, y0, ha = note[m]
        for k, (nm, _, _, _, lab) in enumerate(curves):
            ax.text(x0, y0 - 0.075 * k, lab, transform=ax.transAxes, ha=ha, va='top', fontsize=6.5, color=C_HSE if nm[0] == 'H' else C_PBE)
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$S$ (μV K$^{-1}$)')
    for ax in axs[1]: ax.set_xlabel('$T$ (K)')
    h = [plt.Line2D([], [], ls='', marker='o', ms=4, mfc='white', mec='k', mew=0.8), plt.Line2D([], [], color=C_HSE),
         plt.Line2D([], [], color=C_PBE, ls=DASH)]
    fig.legend(h, ['Experiment', 'HSE06+SOC bands', 'PBE+SOC bands'], loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.935))
    dump('Fig9_S_vs_T', 'experiment_S_T.csv', ['material', 'T_K', 'S_exp_uV/K'], exp_rows)
    dump('Fig9_S_vs_T', 'calculated_S_T.csv', ['material', 'bands', 'method', 'net_doping_p-n_cm-3', 'T_K', 'S_calc_uV/K'], calc_rows)
    dump('Fig9_S_vs_T', 'fit_summary.csv', ['material', 'bands', 'method', 'net_doping_p-n_cm-3', 'fit_T_K', 'fit_S_exp_uV/K',
                                             'max_deviation_calc-exp_uV/K'], fit_rows)
    save(fig, out, 'Fig9_S_vs_T_experiment')
    return {f'{a}|{b}': (c, g) for a, b, _, c, _, _, g in fit_rows}


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


def main():
    global DATA
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.join(ROOT, '21_阶段报告与交接', '论文图_20261028'))
    ap.add_argument('--data'); ap.add_argument('--only', nargs='+'); a = ap.parse_args()
    DATA = a.data or os.path.join(a.out, 'source_data')
    os.makedirs(a.out, exist_ok=True); style()
    todo = a.only or [str(i) for i in range(1, 12)]
    info = {}
    for i in todo:
        res = {'1': fig1, '2': fig2, '3': fig3, '4': fig4, '5': fig5, '6': fig6, '7': fig7, '8': fig8, '9': fig9, '10': fig10, '11': fig11}[i](a.out)
        if res is not None: info[i] = res
    print(info); print('→', a.out, '；源数据 →', DATA)


if __name__ == '__main__':
    main()
