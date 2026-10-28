"""Journal-style figures for the manuscript (English labels; 600-dpi PNG + vector PDF).

Fig. 1  PBE+SOC vs HSE06+SOC band structures (c2 path)                       22_能带基准_1-0/<mat>/c2_{pbe,hse}/path_eigen.txt
Fig. 2  HSE06+SOC bands weighted by cation-s (lone-pair) character             22_能带基准_1-0/<mat>/c2_hse/PROCAR_OPT.gz
Fig. 3  Energy of the highest lone-pair valence valley, PBE vs HSE            same
Fig. 4  Seebeck coefficient at 3e19 cm-3: PBE (BoltzTraP2) vs HSE estimate     10_脚本/multivalley_seebeck.py
Fig. 5  S(T) of Pb2Bi2Te5 and Pb2Sb2Te5: experiment vs HSE-band model          experimental points digitized from the user's plot
Fig. 6  144-atom supercells: ordered (O) and disordered (D) cation sublattice   补算实施方案_20261007/.../parents/<id>/site_map.json
Fig. 7  Neutral formation energies of 12 defect types in O and D               13_DVO缺陷_VASP流程与分析/结果_第一批/analysis/形成能_逐位点.csv
Fig. 8  Ordering-induced change of the formation energy, O − D                  same

Style: Times New Roman (Liberation Serif, same metrics, where TNR is not installed), STIX math, 7–8 pt text, 0.8-pt axes, inward ticks on all four sides, panel labels (a), (b), …;
double column 7.0 in, single column 3.35 in; Okabe–Ito colors (HSE #0072B2, PBE #D55E00, ordered #009E73, disordered #CC79A7).

Usage: python 10_脚本/paper_figures.py [--out 21_阶段报告与交接/论文图_20261028] [--only 1 2 ...]
"""
import argparse, json, os, sys
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


def kticks(ax, x):
    t, l = [], []
    for i in sorted(LAB):
        if i > KEND: break
        t.append(x[i]); l.append(LAB[i])
    for v in t[1:-1]: ax.axvline(v, color=GRID, lw=0.5, zorder=0)
    ax.set_xticks(t, l); ax.set_xlim(x[0], x[KEND]); ax.tick_params(axis='x', top=False, length=0, pad=3)


def bands(m, side):
    E = cv.read_path_eigen(os.path.join(rf.base(m), side, 'path_eigen.txt'))
    nel = NEL[m]; vbm = np.nanmax(E[:, nel - 1]); gap = np.nanmin(E[:, nel]) - vbm
    return E - vbm, gap


def segments(x, y):
    return [(x[:KEND + 1], y[:KEND + 1])]


# ------------------------------------------------------------------ Fig. 1
def fig1(out):
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.9), sharey=True, gridspec_kw=dict(hspace=0.28, wspace=0.06))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        x, _ = rf.path_x(m); nel = NEL[m]
        for side, col, ls, lw, z in (('c2_pbe', C_PBE, (0, (3, 1.6)), 0.9, 2), ('c2_hse', C_HSE, '-', 1.1, 3)):
            E, gap = bands(m, side)
            for b in range(nel - 10, nel + 8):
                for xs, ys in segments(x, E[:, b]):
                    ax.plot(xs, ys, color=col, ls=ls, lw=lw, zorder=z)
            if side == 'c2_pbe': gp = gap
            else: gh = gap
        ax.axhline(0, color='k', lw=0.5, ls=':', zorder=1)
        kticks(ax, x); ax.set_ylim(-0.6, 0.6)
        ax.text(0.015, 0.97, TEX[m], transform=ax.transAxes, ha='left', va='top', fontsize=8,
                bbox=dict(fc='white', ec='none', pad=1.2))
        ax.text(0.985, 0.03, f'$E_g$: PBE {1000 * gp:.0f} meV, HSE {1000 * gh:.0f} meV', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=6.8, bbox=dict(fc='white', ec='none', pad=1.2))
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    h = [plt.Line2D([], [], color=C_HSE, lw=1.1), plt.Line2D([], [], color=C_PBE, lw=0.9, ls=(0, (3, 1.6)))]
    fig.legend(h, ['HSE06+SOC', 'PBE+SOC'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.0))
    save(fig, out, 'Fig1_bands_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 2
def fig2(out):
    cmap = LinearSegmentedColormap.from_list('s', ['#DCE8F4', '#9CC3E4', '#4A90C8', '#1F5FA0', '#08306B'])
    norm = Normalize(0, 0.30)
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.9), sharey=True, gridspec_kw=dict(hspace=0.28, wspace=0.06))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        x, _ = rf.path_x(m); nel = NEL[m]
        E, gap = bands(m, 'c2_hse')
        _, Wp = cv.read_procar(os.path.join(rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fs = rf.cation_s(Wp)
        for b in range(nel - 8, nel + 6):
            segs, cols = [], []
            for i in range(KEND):
                segs.append([(x[i], E[i, b]), (x[i + 1], E[i + 1, b])]); cols.append((fs[i, b] + fs[i + 1, b]) / 2)
            cols = np.array(cols)
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=0.5 + 2.4 * np.clip(cols / 0.3, 0, 1), capstyle='round')
            lc.set_array(cols); ax.add_collection(lc)
        ax.axhline(0, color='k', lw=0.5, ls=':'); ax.axhline(gap, color='k', lw=0.5, ls=':')
        kticks(ax, x); ax.set_ylim(-0.5, 0.5)
        ax.text(0.015, 0.97, f'{TEX[m]}  ($E_g$ = {1000 * gap:.0f} meV)', transform=ax.transAxes, ha='left', va='top',
                fontsize=8, bbox=dict(fc='white', ec='none', pad=1.2))
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.7, pad=0.015, aspect=30, extend='max')
    cb.set_label('Cation-s weight (%)'); cb.set_ticks([0, .1, .2, .3], labels=['0', '10', '20', '30'])
    cb.ax.tick_params(direction='in', length=2, width=0.6); cb.outline.set_linewidth(0.6)
    save(fig, out, 'Fig2_HSE_bands_cation_s')


# ------------------------------------------------------------------ Fig. 3
def fig3(out):
    fig, ax = plt.subplots(figsize=(W1, 2.1))
    kT2 = 2 * 1000 * mv.KB * 300
    ax.axvspan(-kT2, 0, color='#EEF3F8', zorder=0, lw=0)
    ax.text(-kT2 / 2, 3.55, r'$2k_\mathrm{B}T$', ha='center', va='center', fontsize=6.5)
    for j, m in enumerate(MATS[::-1]):
        p, h = rf.lone_pair_offset(m, 'c2_pbe'), rf.lone_pair_offset(m, 'c2_hse')
        ax.plot([p[0], h[0]], [j, j], color='#BDBDBD', lw=1.2, zorder=1)
        ax.plot(p[0], j + 0.08, 'o', ms=4.5, color=C_PBE, mec='white', mew=0.6, zorder=3, label='PBE+SOC' if j == 0 else None)
        ax.plot(h[0], j - 0.08, 'o', ms=4.5, color=C_HSE, mec='white', mew=0.6, zorder=3, label='HSE06+SOC' if j == 0 else None)
        ax.text(h[0] - 7, j - 0.08, f'{h[0]:.0f}'.replace('-', '\u2212'), ha='right', va='center', fontsize=6.3)
    ax.set_yticks(range(4), [TEX[m] for m in MATS[::-1]]); ax.set_ylim(-0.6, 3.85)
    ax.set_xlim(-275, 15); ax.set_xlabel('Highest lone-pair valence valley, $E - E_{\\mathrm{VBM}}$ (meV)')
    ax.tick_params(axis='y', length=0, right=False); ax.legend(loc='lower left', handletextpad=0.2, borderaxespad=0.2)
    save(fig, out, 'Fig3_lone_pair_valley_PBE_vs_HSE')


# ------------------------------------------------------------------ Fig. 4
def fig4(out):
    est = rf.hse_estimates(); prov = not os.path.exists(os.path.join(rf.base('Sn2Bi2Te5'), 'c3_valley', 'c3_masses.json'))
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
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.07), handlelength=1.2)
    save(fig, out, 'Fig4_Seebeck_PBE_vs_HSE')
    return est


# ------------------------------------------------------------------ Fig. 5
EXP = {'Pb2Bi2Te5': [(302, -89), (318, -89), (368, -101), (415, -116), (463, -127), (513, -133), (562, -132), (611, -127)],
       'Pb2Sb2Te5': [(302, 0), (318, 0), (370, -1), (418, -1), (468, 0.5), (520, 4.5), (568, 10), (620, 18)]}


def model(m, side):
    V, gap = mv.build(m, 'v', 300, None, side); C, _ = mv.build(m, 'c', 300, side=side)
    return lambda N, T: mv.solve(V, C, gap, N, T, 'crta')[0]


def solve_N(f, target, T, lo, hi):
    for _ in range(70):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if (f(mid, T) - target) * (f(lo, T) - target) > 0 else (lo, mid)
    return (lo + hi) / 2


def fig5(out):
    Ts = np.arange(300, 641, 20)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.5), gridspec_kw=dict(wspace=0.25))
    # (a) Pb2Bi2Te5: electron density fitted to S(302 K)
    f = model('Pb2Bi2Te5', 'c2_hse')
    lgN = np.log10(np.logspace(17, 20, 31)); vals = [f(-10 ** g, 302) for g in lgN]
    j = next(i for i in range(len(vals) - 1) if (vals[i] + 89) * (vals[i + 1] + 89) <= 0)
    lo, hi = lgN[j], lgN[j + 1]
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if (f(-10 ** mid, 302) + 89) * (f(-10 ** lo, 302) + 89) > 0 else (lo, mid)
    n_fit = 10 ** ((lo + hi) / 2)
    ax = axs[0]
    T, S = zip(*EXP['Pb2Bi2Te5']); ax.plot(T, S, 's', ms=4, mfc='white', mec='k', mew=0.8, label='Exp.', zorder=3)
    ax.plot(Ts, [f(-n_fit, t) for t in Ts], '-', color=C_HSE, label=f'HSE bands, $n$ = {n_fit / 1e18:.1f}' + r'$\times10^{18}$ cm$^{-3}$')
    ax.set_ylabel(r'$S$ (μV K$^{-1}$)'); ax.set_xlabel('$T$ (K)'); ax.set_ylim(-150, -60); ax.set_xlim(280, 640)
    ax.text(0.03, 0.05, TEX['Pb2Bi2Te5'], transform=ax.transAxes, fontsize=8); ax.legend(loc='upper right'); panel(ax, '(a)')
    # (b) Pb2Sb2Te5: compensated sample (S(302 K) = 0) with the HSE gap vs the PBE gap
    ax = axs[1]
    T, S = zip(*EXP['Pb2Sb2Te5']); ax.plot(T, S, 'o', ms=4, mfc='white', mec='k', mew=0.8, label='Exp.', zorder=3)
    for side, col, ls, lab in (('c2_hse', C_HSE, '-', 'HSE bands ($E_g$ = 51 meV)'), ('c2_pbe', C_PBE, (0, (3, 1.6)), 'PBE bands ($E_g$ = 216 meV)')):
        g = model('Pb2Sb2Te5', side)
        N0 = solve_N(g, 0.0, 302, -1e19, 1e19)
        ax.plot(Ts, [g(N0, t) for t in Ts], color=col, ls=ls, label=lab)
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xlabel('$T$ (K)'); ax.set_ylabel(r'$S$ (μV K$^{-1}$)'); ax.set_xlim(280, 640); ax.set_ylim(-140, 40)
    ax.text(0.03, 0.05, TEX['Pb2Sb2Te5'] + ', compensated ($S$ = 0 at 302 K)', transform=ax.transAxes, fontsize=7)
    ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.45)); panel(ax, '(b)')
    save(fig, out, 'Fig5_S_vs_T_experiment')
    return n_fit


# ------------------------------------------------------------------ Fig. 6
def fig6(out):
    base = os.path.join(ROOT, '补算实施方案_20261007', 'rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED', 'parents')
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.55), gridspec_kw=dict(wspace=0.08))
    for ax, pid, title in zip(axs, ('O_ordered', 'D_seed1729'), ('Ordered (O)', 'Disordered (D)')):
        sites = json.load(open(os.path.join(base, pid, 'site_map.json')))['sites']
        for q in range(10):
            row = sorted([s for s in sites if s['q_times_2_mod_10'] == q], key=lambda s: tuple(s['replica']))
            y = 9 - q
            for i, s in enumerate(row):
                if s['sublattice'] == 'anion':
                    ax.plot(i, y, 'o', ms=3.2, color=C_TE, mec='none')
                elif s['is_structural_ghost_vacancy']:
                    ax.plot(i, y, 'o', ms=5.2, mfc='white', mec='k', mew=0.7)
                else:
                    ax.plot(i, y, 'o', ms=5.2, color=C_PB if s['occupancy'] == 'Pb' else C_SB, mec='white', mew=0.4)
        ax.set_xlim(-0.8, 15.8); ax.set_ylim(-0.8, 9.8); ax.set_xticks([]); ax.set_yticks([])
        ax.tick_params(top=False, right=False)
        ax.set_title(title, fontsize=8, pad=3)
    axs[0].annotate('', xy=(-2.9, 9.4), xytext=(-2.9, -0.4), xycoords='data', annotation_clip=False,
                    arrowprops=dict(arrowstyle='->', lw=0.8))
    axs[0].text(-3.3, 4.5, '[111]', rotation=90, ha='right', va='center', fontsize=7)
    for y, lab in zip(range(10), ['Sb', 'Te', 'Pb', 'Te', 'Pb', 'Te', 'Sb', 'Te', 'Vac', 'Te'][::-1]):
        axs[0].text(-1.1, y, lab, ha='right', va='center', fontsize=6.5)
    panel(axs[0], '(a)', x=-0.09); panel(axs[1], '(b)')
    h = [plt.Line2D([], [], ls='', marker='o', ms=5, color=C_PB), plt.Line2D([], [], ls='', marker='o', ms=5, color=C_SB),
         plt.Line2D([], [], ls='', marker='o', ms=5, mfc='white', mec='k'), plt.Line2D([], [], ls='', marker='o', ms=3.5, color=C_TE)]
    fig.legend(h, ['Pb', 'Sb', 'Structural vacancy', 'Te'], loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.06))
    save(fig, out, 'Fig6_supercells_O_D')


# ------------------------------------------------------------------ Fig. 7, 8
NAMES = {'Pb_on_Vac': r'Pb$_{\mathrm{Vac}}$', 'Sb_on_Vac': r'Sb$_{\mathrm{Vac}}$', 'Te_on_Vac': r'Te$_{\mathrm{Vac}}$',
         'Sb_Pb': r'Sb$_{\mathrm{Pb}}$', 'Pb_Sb': r'Pb$_{\mathrm{Sb}}$', 'Te_Pb': r'Te$_{\mathrm{Pb}}$', 'Te_Sb': r'Te$_{\mathrm{Sb}}$',
         'Sb_Te': r'Sb$_{\mathrm{Te}}$', 'Pb_Te': r'Pb$_{\mathrm{Te}}$', 'V_Te': r'$V_{\mathrm{Te}}$', 'V_Pb': r'$V_{\mathrm{Pb}}$', 'V_Sb': r'$V_{\mathrm{Sb}}$'}
ROWS = [(t, q) for t, q in rf.DONORS] + [(t, q.replace('-', '−')) for t, q in rf.ACCEPT]


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
    save(fig, out, 'Fig7_defect_formation_energies')


def fig8(out):
    d = rf.read_ef()
    fig, ax = plt.subplots(figsize=(W1, 3.0))
    ax.axvline(0, color='k', lw=0.5)
    ax.axhspan(len(rf.DONORS) - 0.5, len(ROWS) - 0.5, color='#F2F2F2', zorder=0, lw=0)
    for j, (t, q) in enumerate(ROWS):
        y = len(ROWS) - 1 - j
        o = min(d[('Te-rich', 'O', t)]); sh = o - np.array(d[('Te-rich', 'D', t)])
        ax.barh(y, np.median(sh), height=0.55, color=C_D, alpha=0.35, lw=0)
        ax.plot(sh, [y] * len(sh), 'o', ms=3.2, color=C_D, mec='white', mew=0.4, zorder=3)
    ax.set_yticks(range(len(ROWS)), [f'{NAMES[t]} ({q})' for t, q in ROWS[::-1]]); ax.tick_params(axis='y', length=0, right=False)
    ax.set_xlim(-0.8, 2.0); ax.set_xlabel(r'$\Delta E_f = E_f^{\mathrm{O}} - E_f^{\mathrm{D}}$ (eV)')
    ax.text(-0.76, len(ROWS) - 0.95, 'donors', ha='left', va='center', fontsize=6.5, style='italic')
    ax.text(-0.76, len(rf.ACCEPT) - 0.95, 'acceptors', ha='left', va='center', fontsize=6.5, style='italic')
    save(fig, out, 'Fig8_ordering_induced_change')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.join(ROOT, '21_阶段报告与交接', '论文图_20261028'))
    ap.add_argument('--only', nargs='+'); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); style()
    todo = a.only or [str(i) for i in range(1, 9)]
    info = {}
    for i in todo:
        r = {'1': fig1, '2': fig2, '3': fig3, '4': fig4, '5': fig5, '6': fig6, '7': fig7, '8': fig8}[i](a.out)
        if r is not None: info[i] = r
    print(info); print('→', a.out)


if __name__ == '__main__':
    main()
