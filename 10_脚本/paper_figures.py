"""Journal-style figures for the manuscript (English labels; 600-dpi PNG + vector PDF), numbered as in 论文报告_20261028.

Fig. 1  Structures: O/D 144-atom supercells (a, b) and side views of the user's layered phase vs model O (c, d)
                                                                               补算实施方案_20261007/.../parents/<id>/site_map.json;
                                                                               A2B2Te5_SOC_dense/materials/Pb2Sb2Te5/POSCAR; 13_…/energies.csv
Fig. 2  PBE+SOC vs HSE06+SOC band structures (c2 path)                       22_能带基准_1-0/<mat>/c2_{pbe,hse}/path_eigen.txt
Fig. 3  HSE06+SOC bands weighted by cation-s (lone-pair) character             22_能带基准_1-0/<mat>/c2_hse/PROCAR_OPT.gz
Fig. 4  Energy of the highest lone-pair valence valley, PBE vs HSE            same
Fig. 5  Seebeck coefficient at 3e19 cm-3: PBE (BoltzTraP2) vs HSE estimate     10_脚本/multivalley_seebeck.py
        (+ per-valley conductivity shares and cation-s character for Table 4)
Fig. 6  S(T) of all four compounds: experiment vs calculation (HSE and PBE)    experimental points digitized from the user's plot;
                                                                               A2B2Te5_SOC_dense/materials/<mat>/.../boltztrap_crta.json
Fig. 7  Neutral formation energies of 12 defect types in O and D               13_DVO缺陷_VASP流程与分析/结果_第一批/analysis/形成能_逐位点.csv
Fig. 8  Ordering-induced change of the formation energy, O − D                  same

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
    side_view(ac, PL, pL, 'Layered phase, $P\\bar{3}m1$, 9 atoms (user structure)\nused for bands, $S$ (Figs. 2–6) and defect set L',
              f'$a_\\mathrm{{hex}}$ = {aL:.3f} Å,  $V$ = {HL["V"]:.1f} Å$^3$/f.u.,  $\\Delta E$ = 0')
    side_view(ad, PO, pO, 'Ordered model O, 144 atoms (rocksalt supercell)\nhost for the defect calculations (Figs. 7, 8)',
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
        x, _ = rf.path_x(m); nel = NEL[m]; g = {}
        for side, nm, col, ls, lw, z in (('c2_pbe', 'PBE', C_PBE, DASH, 0.9, 2), ('c2_hse', 'HSE', C_HSE, '-', 1.1, 3)):
            E, g[nm] = bands(m, side); bs = range(nel - 10, nel + 8)
            for b in bs:
                ax.plot(x[:KEND + 1], E[:KEND + 1, b], color=col, ls=ls, lw=lw, zorder=z)
            dump('Fig2_bands_PBE_vs_HSE', f'{m}_{nm.replace("HSE", "HSE06")}+SOC.csv',
                 ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label'] + [bandname(b, nel) + '_E-EVBM_eV' for b in bs],
                 [[i] + [r(v, 6) for v in K[i]] + [r(x[i], 5), LAB.get(i, '')] + [r(E[i, b], 5) for b in bs] for i in range(KEND + 1)])
        gaps.append([m, r(1000 * g['PBE'], 1), r(1000 * g['HSE'], 1)])
        ax.axhline(0, color='k', lw=0.5, ls=':', zorder=1)
        kticks(ax, x); ax.set_ylim(-0.6, 0.6)
        ax.set_title(f'{TEX[m]}     $E_g$: PBE {1000 * g["PBE"]:.0f} meV, HSE {1000 * g["HSE"]:.0f} meV', fontsize=8, pad=4)
        panel(ax, f'({"abcd"[n]})')
    dump('Fig2_bands_PBE_vs_HSE', 'band_gaps_on_path.csv', ['material', 'Eg_PBE+SOC_meV', 'Eg_HSE06+SOC_meV'], gaps)
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
def fig4(out):
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
    dump('Fig4_lone_pair_valley', 'highest_lone_pair_valley.csv',
         ['material', 'functional', 'E-EVBM_meV', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'cation_s_percent'], rows[::-1])
    ax.set_yticks(range(4), [TEX[m] for m in MATS[::-1]]); ax.set_ylim(-0.6, 3.85)
    ax.set_xlim(-275, 15); ax.set_xlabel('Highest lone-pair valence valley, $E - E_{\\mathrm{VBM}}$ (meV)')
    ax.tick_params(axis='y', length=0, right=False); ax.legend(loc='lower left', handletextpad=0.2, borderaxespad=0.2)
    save(fig, out, 'Fig4_lone_pair_valley_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 5
def fig5(out):
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
    dump('Fig5_Seebeck', 'S_300K_3e19.csv', ['material', 'carrier_type', 'carrier_conc_cm-3', 'S_PBE+SOC_BoltzTraP2_uV/K',
                                             'S_HSE06+SOC_estimate_uV/K', 'provisional'],
         [[m, 'pn'[k], 3e19, rf.BTP[m][k], r(est[m][k], 1), 'yes' if (k == 0 and m == 'Sn2Bi2Te5' and prov) else 'no']
          for k in range(2) for m in MATS])
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.07), handlelength=1.2)
    save(fig, out, 'Fig5_Seebeck_PBE_vs_HSE')
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
    K = kpath(); rows, summary = [], []
    for m in MATS:
        c3 = mv.c3_usable(m)[0]
        V, gap = mv.build(m, 'v', 300, c3); C, _ = mv.build(m, 'c', 300)
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fs = rf.cation_s(Wp); nel = NEL[m]
        for N, lab in ((3e19, 'p'), (-3e19, 'n')):
            _, EF, parts = mv.solve(V, C, gap, N, 300)
            lp = 0.0
            for (side, name, w), v in zip(parts, V + C):
                i = valley_index(name); b = nel - 1 if side == 'v' else nel; cs = fs[i, b]
                flag = (side == 'v' and lab == 'p' and cs >= LP_MIN) or (side == 'c' and lab == 'n' and (name in ('Γ',) or name.startswith('Γ–A')))
                lp += w if flag else 0.0
                rows.append([m, lab, r(1000 * EF, 1), 'valence' if side == 'v' else 'conduction', name.split(' [')[0], i]
                            + [r(x, 4) for x in K[i]] + [r(1000 * v['off'], 1), v['g'], r(100 * cs, 1), r(100 * w, 2), 'yes' if flag else 'no'])
            summary.append([m, lab, r(100 * lp, 1)])
    dump('Fig5_Seebeck', 'valley_conductivity_shares.csv',
         ['material', 'carrier_type', 'EF_minus_VBM_meV', 'band', 'valley', 'k_index_on_path', 'k_frac_1', 'k_frac_2', 'k_frac_3',
          'offset_from_band_edge_meV', 'degeneracy_g', 'cation_s_percent', 'conductivity_share_percent',
          'counted (p: lone-pair valence valley, cation s >= 15%; n: Gamma conduction valley)'], rows)
    dump('Fig5_Seebeck', 'lone_pair_share_summary.csv',
         ['material', 'carrier_type', 'share_percent (p: lone-pair valence valleys; n: Gamma conduction valley)'], summary)
    print('shares', summary)
    c3 = mv.c3_load('Sn2Bi2Te5')
    if c3:
        rows = []
        for k, v in c3.items():
            H = np.array(v['H_eVA2_line_perp_c']); g = np.array(v['gradient'])
            ok, why = mv.c3_verdict(v, 'v')
            dk = -np.linalg.solve(H, g); dE = -0.5 * g @ np.linalg.solve(H, g)
            rows.append([k, r(v['E0'], 5)] + [r(x, 4) for x in g] + [r(x, 4) for x in v['principal']]
                        + [r(x, 4) for x in dk] + [r(1000 * dE, 1), 'yes' if ok else 'no', why])
        dump('Fig5_Seebeck', 'Sn2Bi2Te5_c3_valley_check.csv',
             ['valley_k_index', 'E_center_eV', 'grad_line_eVA', 'grad_perp_eVA', 'grad_c_eVA',
              'principal_mass_1_me', 'principal_mass_2_me', 'principal_mass_3_me',
              'stationary_shift_line_1/A', 'stationary_shift_perp_1/A', 'stationary_shift_c_1/A',
              'quadratic_extrapolated_dE_meV', 'used_in_model', 'reason_not_used'], rows)


# ------------------------------------------------------------------ Fig. 6
# 实验 S(T)，从所给曲线图读取（µV/K）
# 实验 S(T)（µV/K）。Pb₂Sb₂Te₅ 为用户给的原始数据（2026-10-29），其余三种从所给曲线图读取
EXP = {'Pb2Sb2Te5': [(301.58, 5.52), (319.32, 5.39), (370.16, 5.60), (418.62, 4.94), (468.56, 5.09), (519.25, 6.93), (568.17, 10.67), (619.36, 17.10)],
       'Pb2Bi2Te5': [(302, -89), (318, -89), (368, -101), (415, -116), (463, -127), (513, -133), (562, -132), (611, -127)],
       'Sn2Sb2Te5': [(300, 16), (320, 16), (370, 19), (418, 23), (467, 28), (518, 33), (567, 37), (618, 42)],
       'Sn2Bi2Te5': [(302, 30), (318, 31), (368, 39), (417, 46), (466, 55), (518, 62), (567, 66), (619, 67)]}
EXP_SRC = {'Pb2Sb2Te5': 'raw data from user', 'Pb2Bi2Te5': 'digitized from plot', 'Sn2Sb2Te5': 'digitized from plot', 'Sn2Bi2Te5': 'digitized from plot'}
# Pb₂Sb₂Te₅ 的实验电导率（与 S 同一组温度；单位未给出，只用相对 302 K 的比值）
SIG_PBSB = [21.82167, 20.56344, 17.39372, 14.98943, 13.02654, 11.30616, 10.02781, 9.0085]
MARK = {'Pb2Sb2Te5': 'o', 'Pb2Bi2Te5': 's', 'Sn2Sb2Te5': 'v', 'Sn2Bi2Te5': '^'}


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


def fig6(out):
    """上排：Pb₂Sb₂Te₅（实测原始数据）的 S(T) 与 σ(T)；下排：其余三种材料的 S(T)。"""
    Ts = np.arange(300, 641, 10)
    fig = plt.figure(figsize=(W2, 5.0))
    gs = fig.add_gridspec(2, 6, hspace=0.45, wspace=1.3)
    axs = [fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])] + [fig.add_subplot(gs[1, 2 * k:2 * k + 2]) for k in range(3)]
    order = ['Pb2Sb2Te5', None, 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
    titles = {'Pb2Sb2Te5': 'near-compensated', 'Pb2Bi2Te5': 'n-type', 'Sn2Sb2Te5': 'p-type', 'Sn2Bi2Te5': 'p-type'}
    ylim = {'Pb2Sb2Te5': (-140, 30), 'Pb2Bi2Te5': (-305, -62), 'Sn2Sb2Te5': (0, 57), 'Sn2Bi2Te5': (0, 90)}
    legloc = {'Pb2Sb2Te5': dict(loc='center right', bbox_to_anchor=(1.0, 0.42)), 'Pb2Bi2Te5': dict(loc='lower left', fontsize=6),
              'Sn2Sb2Te5': dict(loc='upper left', fontsize=6), 'Sn2Bi2Te5': dict(loc='upper left', fontsize=6)}
    exp_rows, calc_rows, fit_rows, sig_rows = [], [], [], []
    for ax, m in zip(axs, order):
        if m is None: continue
        Te, Se = zip(*EXP[m]); T0, S0 = Te[0], Se[0]
        exp_rows += [[m, t, s, EXP_SRC[m]] + ([SIG_PBSB[k]] if m == 'Pb2Sb2Te5' else ['']) for k, (t, s) in enumerate(EXP[m])]
        curves = []
        if m == 'Pb2Sb2Te5':        # 近补偿：用 301.6 K 的 S 定净掺杂（N 在 ±10¹⁹ 内单调，二分）；比较两种带隙
            for side, nm in (('c2_hse', 'HSE06+SOC'), ('c2_pbe', 'PBE+SOC')):
                g = model(m, side); N = solve_N(g, S0, T0, -1e19, 1e19)
                V, gp = mv.build(m, 'v', 300, None, side); C, _ = mv.build(m, 'c', 300, side=side)
                st = [mv.sigma_tau(V, C, gp, N, t) for t in Te]
                sig_rows += [[nm, f'{N:.3e}', t, r(SIG_PBSB[k] / SIG_PBSB[0], 4), r(st[k][0] / st[0][0], 4), f'{st[k][1]:.3e}', f'{st[k][2]:.3e}']
                             for k, t in enumerate(Te)]
                curves.append((nm, 'path multivalley model', g, N, f'{nm[:3]} bands ($E_g$ = {1000 * gp:.0f} meV)', st))
        elif m == 'Pb2Bi2Te5':      # n 型：302 K 的实验值定电子浓度
            for side, nm in (('c2_hse', 'HSE06+SOC'), ('c2_pbe', 'PBE+SOC')):
                g = model(m, side); N = solve_log(g, EXP[m], 17, 20, -1)
                curves.append((nm, 'path multivalley model', g, N, f'{nm[:3]}, $n$ = {sci(N)}', None))
        else:                       # Sn 系重掺杂 p 型：BoltzTraP2 全区；HSE = BoltzTraP2 × 模型 HSE/PBE 比值（1.7 节第 5 步）
            b = bt2(m); gh, gp_ = model(m, 'c2_hse'), model(m, 'c2_pbe')
            hse = lambda p, t, b=b, gh=gh, gp_=gp_: b(p, t) * gh(p, t) / gp_(p, t)
            for nm, meth, g in (('HSE06+SOC', 'BoltzTraP2 (PBE) x model HSE/PBE ratio', hse), ('PBE+SOC', 'BoltzTraP2 (PBE, full BZ)', b)):
                N = solve_log(g, EXP[m], 20, 21.8)
                curves.append((nm, meth, g, N, f'{nm[:3]}, $p$ = {sci(N)}', None))
        ax.plot(Te, Se, MARK[m], ms=4, mfc='white', mec='k', mew=0.8, ls='', label='Exp.', zorder=3)
        for nm, meth, g, N, lab, _ in curves:
            v = [g(N, t) for t in Ts]
            ax.plot(Ts, v, color=C_HSE if nm[0] == 'H' else C_PBE, ls='-' if nm[0] == 'H' else DASH, label=lab)
            calc_rows += [[m, nm, meth, f'{N:.3e}', int(t), r(s, 2)] for t, s in zip(Ts, v)]
            dev = [g(N, t) - s for t, s in EXP[m]]
            fit_rows.append([m, nm, meth, f'{N:.3e}', T0, S0, r(max(dev, key=abs), 1)])
        if m == 'Pb2Sb2Te5': ax.axhline(0, color='k', lw=0.5, ls=':')
        ax.set_xlim(280, 640); ax.set_ylim(*ylim[m]); ax.set_xlabel('$T$ (K)')
        ax.set_title(f'{TEX[m]} ({titles[m]})', fontsize=8, pad=4)
        ax.legend(**legloc[m])
        if m == 'Pb2Sb2Te5': pbsb = curves
    for ax in (axs[0], axs[2]): ax.set_ylabel(r'$S$ (μV K$^{-1}$)')
    # (b) Pb₂Sb₂Te₅ 的 σ(T)：实测值与模型的 σ/τ，都归一到 301.6 K
    ax = axs[1]; Te = [t for t, _ in EXP['Pb2Sb2Te5']]
    ax.plot(Te, np.array(SIG_PBSB) / SIG_PBSB[0], 'o', ms=4, mfc='white', mec='k', mew=0.8, ls='', label=r'Exp. $\sigma$', zorder=3)
    for nm, _, _, N, _, st in pbsb:
        y = [x[0] / st[0][0] for x in st]
        ax.plot(Te, y, color=C_HSE if nm[0] == 'H' else C_PBE, ls='-' if nm[0] == 'H' else DASH, marker='', label=f'{nm[:3]} bands, $\\sigma/\\tau$')
    k = np.polyfit(np.log(Te), np.log(SIG_PBSB), 1)[0]
    ax.text(300, 0.40, f'Exp.: $\\sigma \\propto T^{{{k:.2f}}}$'.replace('-', '\u2212'), ha='left', va='center', fontsize=7)
    ax.set_yscale('log'); ax.set_ylim(0.3, 60); ax.set_xlim(280, 640); ax.set_xlabel('$T$ (K)')
    ax.set_ylabel(r'$\sigma(T)/\sigma(302\ \mathrm{K})$')
    ax.set_title(f'{TEX["Pb2Sb2Te5"]}: electrical conductivity', fontsize=8, pad=4)
    ax.legend(loc='upper left')
    for k, ax in enumerate(axs): panel(ax, f'({"abcde"[k]})')
    dump('Fig6_S_vs_T', 'experiment_S_T.csv', ['material', 'T_K', 'S_exp_uV/K', 'source', 'sigma_exp_(unit_not_given)'], exp_rows)
    dump('Fig6_S_vs_T', 'calculated_S_T.csv', ['material', 'bands', 'method', 'net_doping_p-n_cm-3', 'T_K', 'S_calc_uV/K'], calc_rows)
    dump('Fig6_S_vs_T', 'fit_summary.csv', ['material', 'bands', 'method', 'net_doping_p-n_cm-3', 'fit_T_K', 'fit_S_exp_uV/K',
                                             'max_deviation_calc-exp_uV/K'], fit_rows)
    dump('Fig6_S_vs_T', 'Pb2Sb2Te5_conductivity.csv', ['bands', 'net_doping_p-n_cm-3', 'T_K', 'sigma_exp_ratio_to_302K',
                                                       'sigma_over_tau_model_ratio_to_302K', 'p_model_cm-3', 'n_model_cm-3'], sig_rows)
    save(fig, out, 'Fig6_S_vs_T_experiment')
    return {f'{a}|{b}': (c, g) for a, b, _, c, _, _, g in fit_rows}


# ------------------------------------------------------------------ Fig. 7, 8
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


def fig7(out):
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
    dump('Fig7_formation_energies', 'formation_energy_per_site.csv',
         ['phase', 'defect_type', 'nominal_charge', 'site', 'Ef_Te-rich_eV', 'Ef_Te-poor_eV'],
         [[a, b, c, s, r(x, 4), r(y, 4)] for a, b, c, s, x, y in sorted(ef_sites(), key=lambda t: (order[t[1]], t[0] != 'O', t[3]))])
    mu = json.load(open(os.path.join(DEFDIR, 'DVO_结论摘要.json')))['chemical_potentials']['conditions']
    dump('Fig7_formation_energies', 'chemical_potentials.csv', ['condition', 'mu_Pb_eV', 'mu_Sb_eV', 'mu_Te_eV', 'dmu_Pb_eV', 'dmu_Sb_eV', 'dmu_Te_eV'],
         [[c] + [r(mu[c][e], 4) for e in ('Pb', 'Sb', 'Te')] + [r(mu[c]['delta'][e], 4) for e in ('Pb', 'Sb', 'Te')] for c in ('Te-rich', 'Te-poor')])
    save(fig, out, 'Fig7_defect_formation_energies')


def fig8(out):
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
    dump('Fig8_ordering', 'delta_Ef_per_site.csv', ['defect_type', 'nominal_charge', 'D_site', 'Ef_O_lowest_site_eV', 'Ef_D_eV', 'dEf_O-D_eV'], rows)
    dump('Fig8_ordering', 'delta_Ef_median.csv', ['defect_type', 'nominal_charge', 'n_D_sites', 'median_dEf_O-D_eV'], med)
    save(fig, out, 'Fig8_ordering_induced_change')


def main():
    global DATA
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.join(ROOT, '21_阶段报告与交接', '论文图_20261028'))
    ap.add_argument('--data'); ap.add_argument('--only', nargs='+'); a = ap.parse_args()
    DATA = a.data or os.path.join(a.out, 'source_data')
    os.makedirs(a.out, exist_ok=True); style()
    todo = a.only or [str(i) for i in range(1, 9)]
    info = {}
    for i in todo:
        res = {'1': fig1, '2': fig2, '3': fig3, '4': fig4, '5': fig5, '6': fig6, '7': fig7, '8': fig8}[i](a.out)
        if res is not None: info[i] = res
    print(info); print('→', a.out, '；源数据 →', DATA)


if __name__ == '__main__':
    main()
