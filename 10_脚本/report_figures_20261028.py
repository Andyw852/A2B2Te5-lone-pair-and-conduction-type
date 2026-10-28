"""结论报告（2026-10-28）的插图，全部从仓库数据生成。

图 1  HSE 能带按阳离子 s 成分着色（孤对电子态在哪里）           22_能带基准_1-0/<材料>/c2_hse/{path_eigen.txt, PROCAR_OPT.gz}
图 2  最高的孤对电子价带谷相对价带顶的能量：PBE → HSE           同上 + c2_pbe
图 3  同浓度 Seebeck 系数：PBE（BoltzTraP2）与 HSE 估计          10_脚本/multivalley_seebeck.py
图 4  12 类缺陷的中性形成能：有序相（退火后）与无序相（退火前）    13_DVO缺陷_VASP流程与分析/结果_第一批/analysis/形成能_逐位点.csv
图 5  有序化带来的形成能变化（O − D，与化学势无关）               同上

用法：python 10_脚本/report_figures_20261028.py [--out 21_阶段报告与交接/图_20261028]
"""
import argparse, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c2_valleys as cv
import multivalley_seebeck as mv

MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
TEX = {m: m.replace('2', '$_2$').replace('5', '$_5$') for m in MATS}
NEL = {m: 88 if m[3:5] == 'Bi' else 68 for m in MATS}
C_PBE, C_HSE = '#eb6834', '#2a78d6'           # 参考调色板 1、2 号（已验校）
C_O, C_D = '#1baf7a', '#6250d6'               # 3、7 号：有序相、无序相（绿色对比度 < 3:1，图中配直接标注）
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d9d8d3'
RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
LAB = {0: 'Γ', 33: 'M', 52: 'K', 90: 'Γ', 97: 'A', 130: 'L', 149: 'H', 187: 'A', 188: 'L', 195: 'M', 196: 'H', 203: 'K'}


def setup_fonts():
    for f in font_manager.findSystemFonts():
        if 'wqy' in f.lower() or 'notosanscjk' in f.lower().replace('-', ''):
            font_manager.fontManager.addfont(f)
            plt.rcParams['font.family'] = font_manager.FontProperties(fname=f).get_name(); break
    plt.rcParams.update({'mathtext.fontset': 'dejavusans', 'axes.unicode_minus': False, 'axes.edgecolor': INK2,
                         'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.8,
                         'font.size': 9, 'savefig.dpi': 200, 'savefig.bbox': 'tight', 'figure.facecolor': 'white'})


def base(m):
    return os.path.join(ROOT, '22_能带基准_1-0', '' if m == 'Pb2Sb2Te5' else m)


def path_x(m):
    K = cv.read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    A = cv.read_lat(os.path.join(ROOT, '08_介电复核', '源数据', m, 'POSCAR'))
    return cv.kdist(K, 2 * np.pi * np.linalg.inv(A).T), A


def cation_s(W):
    """W[k, band, ion, orb]：阳离子（离子 1–2 为 B，8–9 为 A）s 成分占球内总投影的比例。"""
    return W[:, :, [0, 1, 7, 8], 0].sum(-1) / np.maximum(W.sum((-1, -2)), 1e-9)


def ticks(ax, x):
    t, l = [], []
    for i in sorted(LAB):
        if (i - 1, i) in cv.BREAKS: l[-1] += '|' + LAB[i]; continue
        t.append(x[i]); l.append(LAB[i])
    for v in t[1:-1]: ax.axvline(v, color=GRID, lw=0.6, zorder=0)
    ax.set_xticks(t, l, fontsize=7.5); ax.set_xlim(x[0], x[-1])


# ---------------------------------------------------------------- 图 1
def fig1(out):
    cmap = LinearSegmentedColormap.from_list('s', RAMP); norm = Normalize(0, 0.30)
    fig, axs = plt.subplots(2, 2, figsize=(11.0, 6.8), sharey=True)
    for ax, m in zip(axs.flat, MATS):
        x, _ = path_x(m); nel = NEL[m]
        E = cv.read_path_eigen(os.path.join(base(m), 'c2_hse', 'path_eigen.txt'))
        _, W = cv.read_procar(os.path.join(base(m), 'c2_hse', 'PROCAR_OPT.gz'))
        fs = cation_s(W); vbm = np.nanmax(E[:, nel - 1]); gap = np.nanmin(E[:, nel]) - vbm
        for b in range(nel - 6, nel + 4):
            y = E[:, b] - vbm; f = fs[:, b]
            segs = [[(x[i], y[i]), (x[i + 1], y[i + 1])] for i in range(len(x) - 1) if (i, i + 1) not in cv.BREAKS]
            cols = [(f[i] + f[i + 1]) / 2 for i in range(len(x) - 1) if (i, i + 1) not in cv.BREAKS]
            cols = np.array(cols)
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=0.9 + 2.6 * np.clip(cols / 0.30, 0, 1), capstyle='round'); lc.set_array(cols)
            ax.add_collection(lc)
        ax.axhline(0, color=INK2, lw=0.6, ls=':'); ax.axhline(gap, color=INK2, lw=0.6, ls=':')
        ticks(ax, x); ax.set_ylim(-0.5, 0.5)
        ax.set_title(f'{TEX[m]}（HSE06+SOC，$E_g$ = {1000 * gap:.0f} meV）', fontsize=10, color=INK)
    for ax in axs[:, 0]: ax.set_ylabel(r'$E - E_\mathrm{VBM}$ (eV)')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.75, pad=0.02, extend='max')
    cb.set_label('阳离子 s 成分（%，孤对电子 s 态；线越粗越深越多）'); cb.set_ticks([0, .1, .2, .3], labels=['0', '10', '20', '≥30'])
    fig.savefig(os.path.join(out, '图1_HSE能带_阳离子s成分.png')); plt.close(fig)


# ---------------------------------------------------------------- 图 2
def lone_pair_offset(m, side):
    x, _ = path_x(m); nel = NEL[m]
    E = cv.read_path_eigen(os.path.join(base(m), side, 'path_eigen.txt'))
    _, W = cv.read_procar(os.path.join(base(m), side, 'PROCAR_OPT.gz'))
    v = E[:, nel - 1]; fs = cation_s(W)[:, nel - 2:nel].mean(1)
    best = None
    for i, (dE, _m) in cv.extrema(x, v, 'max', 600):
        if i in (0, 90) and v[91] > v[90]: continue                 # Γ 沿 c 不是极大：鞍点
        if i in (97, 187) and v[96] > v[97]: continue               # A 同理
        if i not in (0, 90, 97, 187) and i <= 187:
            j = i + 97 if i <= 90 else i - 97                       # 同一面内位置的另一个 kz
            if 0 <= j < len(v) and v[j] > v[i]: continue            # 沿 c 不是极大：鞍点
        if fs[i] >= 0.15 and (best is None or dE > best[0]): best = (dE, i, fs[i])
    return best


def fig2(out):
    fig, ax = plt.subplots(figsize=(6.4, 2.9))
    kT2 = 2 * 1000 * mv.KB * 300
    ax.axvspan(-kT2, 0, color='#eef3fb', zorder=0); ax.text(-kT2 / 2, 3.62, '2$k_BT$（300 K）', ha='center', va='bottom', fontsize=7.5, color=INK2)
    for j, m in enumerate(MATS[::-1]):
        p, h = lone_pair_offset(m, 'c2_pbe'), lone_pair_offset(m, 'c2_hse')
        ax.plot([p[0], h[0]], [j, j], color=GRID, lw=2, zorder=1)
        ax.scatter([p[0]], [j + 0.09], s=46, color=C_PBE, zorder=3, edgecolor='white', linewidth=1.2, label='PBE+SOC' if j == 0 else None)
        ax.scatter([h[0]], [j - 0.09], s=46, color=C_HSE, zorder=3, edgecolor='white', linewidth=1.2, label='HSE06+SOC' if j == 0 else None)
        ax.text(h[0] - 6, j - 0.36, f'{h[0]:.0f}', ha='right', fontsize=7.5, color=INK2)
        ax.text(p[0] + 4, j + 0.24, f'{p[0]:.0f}', ha='left', fontsize=7.5, color=INK2)
    ax.set_yticks(range(4), [TEX[m] for m in MATS[::-1]]); ax.set_ylim(-0.6, 3.9)
    ax.set_xlim(-262, 22); ax.set_xlabel('最高的孤对电子价带谷相对价带顶的能量（meV）')
    ax.grid(axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.legend(loc='lower left', frameon=False, fontsize=8)
    fig.savefig(os.path.join(out, '图2_孤对电子谷位置_PBE与HSE.png')); plt.close(fig)


# ---------------------------------------------------------------- 图 3
BTP = {'Pb2Sb2Te5': (204, -125), 'Pb2Bi2Te5': (208, -128), 'Sn2Sb2Te5': (202, -159), 'Sn2Bi2Te5': (204, -186)}   # 2.4 节，300 K，3e19


def hse_estimates():
    est = {}
    for m in MATS:
        r = {}
        for side in ('c2_pbe', 'c2_hse'):
            c3 = mv.c3_usable(m)[0] if side == 'c2_hse' else None
            V, gap = mv.build(m, 'v', 300, c3, side); C, _ = mv.build(m, 'c', 300, side=side)
            r[side] = [mv.solve(V, C, gap, N, 300, 'crta')[0] for N in (3e19, -3e19)]
        est[m] = [BTP[m][i] * r['c2_hse'][i] / r['c2_pbe'][i] for i in range(2)]
    return est


def fig3(out, est, snbi_provisional):
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.2))
    xs = np.arange(4); w = 0.36
    for k, (ax, title) in enumerate(zip(axs, ('p 型，3×10$^{19}$ cm$^{-3}$', 'n 型，3×10$^{19}$ cm$^{-3}$（取 |S|）'))):
        pb = [abs(BTP[m][k]) for m in MATS]; hs = [abs(est[m][k]) for m in MATS]
        ax.bar(xs - w / 2 - 0.01, pb, w, color=C_PBE, label='PBE+SOC（BoltzTraP2，2.4 节）')
        bars = ax.bar(xs + w / 2 + 0.01, hs, w, color=C_HSE, label='HSE06+SOC（多谷模型，比值修正）')
        if k == 0 and snbi_provisional:
            bars[3].set_hatch('////'); bars[3].set_edgecolor('white')
            ax.text(xs[3] + w / 2, hs[3] + 5, f'{hs[3]:.0f}（暂定）', ha='center', fontsize=7.5, color=INK2)
        for n, (xi, v) in enumerate(zip(xs + w / 2 + 0.01, hs)):
            if k == 0 and snbi_provisional and n == 3: continue
            ax.text(xi, v / 2, f'{v:.0f}', ha='center', va='center', fontsize=7.5, color='white')
        ax.set_xticks(xs, [TEX[m] for m in MATS]); ax.set_title(title, fontsize=9.5, color=INK)
        ax.set_ylim(0, 245); ax.grid(axis='y', color=GRID, lw=0.5); ax.set_axisbelow(True)
        for s in ('top', 'right'): ax.spines[s].set_visible(False)
    axs[0].set_ylabel('Seebeck 系数（µV/K，300 K，面内）')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=2, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, '图3_同浓度Seebeck_PBE与HSE.png')); plt.close(fig)


# ---------------------------------------------------------------- 图 4、5
DONORS = [('Pb_on_Vac', '+2'), ('Sb_on_Vac', '+3'), ('Sb_Pb', '+1'), ('Te_Pb', '+2'), ('Te_Sb', '+1'), ('V_Te', '+2')]
ACCEPT = [('Pb_Sb', '-1'), ('Sb_Te', '-1'), ('V_Pb', '-2'), ('Te_on_Vac', '-2'), ('V_Sb', '-3'), ('Pb_Te', '-2')]


def tex_def(t):
    a, b = t.split('_', 1)
    b = {'on_Vac': '空位'}.get(b, b)
    return (f'{a} 填空位' if b == '空位' else (f'$V_{{\\rm {b}}}$' if a == 'V' else f'{a}$_{{\\rm {b}}}$'))


def read_ef():
    f = os.path.join(ROOT, '13_DVO缺陷_VASP流程与分析', '结果_第一批', 'analysis', '形成能_逐位点.csv')
    d = {}
    for r in csv.DictReader(open(f, encoding='utf-8-sig')):
        fam = 'O' if r['parent_id'] == 'O_ordered' else 'D'
        for c in ('Te-rich', 'Te-poor'):
            d.setdefault((c, fam, r['defect_type']), []).append(float(r[f'Ef0_{c}']))
    return d


def fig4(out, d):
    rows = DONORS + ACCEPT
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 5.0), sharey=True)
    for ax, c in zip(axs, ('Te-rich', 'Te-poor')):
        ax.axvline(0, color=INK2, lw=0.7)
        ax.axhspan(len(DONORS) - 0.5, len(rows) - 0.5, color='#f6f5f2', zorder=0)
        for j, (t, q) in enumerate(rows):
            y = len(rows) - 1 - j
            for fam, col, dy in (('O', C_O, 0.17), ('D', C_D, -0.17)):
                v = d.get((c, fam, t), [])
                ax.scatter(v, [y + dy] * len(v), s=34, color=col, edgecolor='white', linewidth=1, zorder=3,
                           label=({'O': '有序相 O（退火后）', 'D': '无序相 D（退火前）'}[fam] if j == 0 else None))
        ax.set_xlim(-0.45, 3.05); ax.set_xlabel('中性形成能（eV）')
        ax.set_title({'Te-rich': '富 Te（Δμ$_{\\rm Te}$ = 0）', 'Te-poor': '贫 Te（Δμ$_{\\rm Sb}$ = 0）'}[c], fontsize=9.5, color=INK)
        ax.grid(axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True)
        for s in ('top', 'right'): ax.spines[s].set_visible(False)
    axs[0].set_yticks(range(len(rows)), [f'{tex_def(t)}（{q}）' for t, q in rows[::-1]])
    axs[0].text(-0.43, len(rows) - 0.75, '施主', fontsize=8, color=INK2, va='top')
    axs[0].text(-0.43, len(ACCEPT) - 0.75, '受主', fontsize=8, color=INK2, va='top')
    axs[1].legend(loc='lower right', frameon=True, framealpha=0.95, fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, '图4_缺陷形成能_有序与无序.png')); plt.close(fig)


def fig5(out, d):
    rows = DONORS + ACCEPT
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.axvline(0, color=INK2, lw=0.7)
    ax.axhspan(len(DONORS) - 0.5, len(rows) - 0.5, color='#f6f5f2', zorder=0)
    for j, (t, q) in enumerate(rows):
        y = len(rows) - 1 - j
        o = min(d[('Te-rich', 'O', t)]); dv = np.array(d[('Te-rich', 'D', t)])      # 对有序相的最低位点（与 D_O对比.md 同口径）
        shift = o - dv                                   # 有序化带来的升高（正 = 退火后变贵）
        ax.barh(y, np.median(shift), height=0.55, color=C_D, alpha=0.35, zorder=1)
        ax.scatter(shift, [y] * len(shift), s=26, color=C_D, edgecolor='white', linewidth=1, zorder=3)
        ax.text(max(shift.max(), 0) + 0.05, y, f'{np.median(shift):+.2f}' + ('' if len(shift) > 1 else '（1 个位点）'), va='center', fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(rows)), [f'{tex_def(t)}（{q}）' for t, q in rows[::-1]])
    ax.set_xlim(-0.9, 2.3); ax.set_xlabel('有序化（退火）带来的形成能变化 O - D（eV；O 取最低位点；正 = 退火后更难形成）')
    ax.text(-0.88, len(rows) - 0.75, '施主', fontsize=8, color=INK2, va='top'); ax.text(-0.88, len(ACCEPT) - 0.75, '受主', fontsize=8, color=INK2, va='top')
    ax.grid(axis='x', color=GRID, lw=0.5); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(out, '图5_有序化带来的形成能变化.png')); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.join(ROOT, '21_阶段报告与交接', '图_20261028'))
    ap.add_argument('--only', nargs='+')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); setup_fonts()
    todo = a.only or ['1', '2', '3', '4', '5']
    if '1' in todo: fig1(a.out)
    if '2' in todo: fig2(a.out)
    if '3' in todo:
        est = hse_estimates()
        snbi_prov = mv.c3_provisional('Sn2Bi2Te5')
        print('HSE 估计（p 3e19, n 3e19）：', {m: [round(v) for v in est[m]] for m in MATS}, '；SnBi p 暂定' if snbi_prov else '')
        fig3(a.out, est, snbi_prov)
    d = read_ef()
    if '4' in todo: fig4(a.out, d)
    if '5' in todo: fig5(a.out, d)
    print('→', a.out)


if __name__ == '__main__':
    main()
