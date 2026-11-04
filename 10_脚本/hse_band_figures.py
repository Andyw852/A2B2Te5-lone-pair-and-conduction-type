#!/usr/bin/env python3
"""HSE06+SOC 能带专题（2026-11-04 用户要求单独成册）：所有 HSE 能带、轨道成分、价带谷与口袋、镜面搜索、Γ 点能级。

与论文报告的关系：数据与零点约定完全相同（复用 paper_figures 的读数与判据），只是重新组织、每张图只做一件事，
每个子图有自己的纵坐标（用户要求）。图 6 拆成 H5（四种材料的价带能级）与 H6（Sn 系口袋在布里渊区中的位置，去掉三维图）；
图 S1(c)(d) 删去（用户意见），路径值与真实极大之差并入 H7 的数据；Sn 系的附加路径 Γ–Σ′–R₂′–L 单独画在能带图右侧。

  H1  HSE06+SOC 能带（四种材料）
  H2  阳离子 s 成分（孤对电子 s 态的权重）
  H3  阳离子 s 来自 A 位还是 B 位
  H4a–d  每种材料的 s、p 成分：阳离子 s、A 位 p、B 位 p、Te p
  H5  四种材料的价带能级（谷的名称、能量、等价口袋数）
  H6  Sn 系价带口袋在布里渊区中的位置（沿 c 俯视）
  H7  Sn 系含 c 轴的镜面内最高价带的能量
  H8  Γ 点能级与 SOC（目前只有 PBE；HSE 不加 SOC 的计算见 22_能带基准_1-0/SOC机制检验_HSE/）

能量零点：真实价带顶（Sn 系取镜面内三维搜索的结果，Pb 系取能带路径上的最高点）。
用法：python 10_脚本/hse_band_figures.py [--out 21_阶段报告与交接/HSE能带_20261104/图] [--only H1 H5 ...]
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

MATS, TEX, NEL, KEND, r = pf.MATS, pf.TEX, pf.NEL, pf.KEND, pf.r
SN = ('Sn2Sb2Te5', 'Sn2Bi2Te5')
W2 = pf.W2
CM = {   # 每种成分一个单色渐变（浅 → 深）
    'cation_s': LinearSegmentedColormap.from_list('s', ['#DCE8F4', '#9CC3E4', '#4A90C8', '#1F5FA0', '#08306B']),
    'A_p': LinearSegmentedColormap.from_list('ap', ['#EDE7F6', '#C3B1E1', '#8E6CC4', '#5D3A9B', '#331A66']),
    'B_p': LinearSegmentedColormap.from_list('bp', ['#FDEBDD', '#F8C39A', '#F08C4A', '#D45A0F', '#8C3A05']),
    'Te_p': LinearSegmentedColormap.from_list('tp', ['#E5F5E0', '#A8DCA0', '#5DB461', '#2C8A3F', '#0B5A24']),
}
WMAX = {'cation_s': 0.30, 'A_p': 0.60, 'B_p': 0.60, 'Te_p': 0.80}
LONG = {'cation_s': 'cation s (A-s + B-s)', 'A_p': 'A-site p', 'B_p': 'B-site p', 'Te_p': 'Te p'}


def comp(W):
    """W[k, band, ion, orb]（离子 B×2、Te×5、A×2；轨道 s、p_y、p_z、p_x、d…）→ 各成分占球内总投影的比例。"""
    t = np.maximum(W.sum((-1, -2)), 1e-9)
    c = {'A_s': W[:, :, 7:9, 0].sum(-1) / t, 'B_s': W[:, :, 0:2, 0].sum(-1) / t, 'Te_s': W[:, :, 2:7, 0].sum(-1) / t,
         'A_p': W[:, :, 7:9, 1:4].sum((-1, -2)) / t, 'B_p': W[:, :, 0:2, 1:4].sum((-1, -2)) / t, 'Te_p': W[:, :, 2:7, 1:4].sum((-1, -2)) / t}
    c['cation_s'] = c['A_s'] + c['B_s']
    return c


def hse(m):
    """HSE06+SOC 路径能带（零点为真实价带顶，同论文图 3）、附加路径（Sn 系）、两段的轨道成分。"""
    x, _ = pf.rf.path_x(m); nel = NEL[m]
    E, _ = pf.bands(m, 'c2_hse')
    T5 = pf.mv.true_load(m); dv = (T5['vbm'] - T5['path_vbm']) if T5 else 0.0
    E = E - dv
    _, W = pf.cv.read_procar(os.path.join(pf.rf.base(m), 'c2_hse', 'PROCAR_OPT.gz'))
    tp = pf.true_path_bands(m)
    if tp:
        tp['x'] = pf.TP_XSCALE * tp['x']; tp['c'] = comp(tp['W'])
    cbm = min(np.nanmin(E[:KEND + 1, nel]), np.nanmin(tp['E'][:, nel]) if tp else np.inf)
    return dict(m=m, x=x, X0=x[KEND], nel=nel, E=E, c=comp(W), tp=tp, gap=cbm, gap_path=np.nanmin(E[:KEND + 1, nel]) - np.nanmax(E[:KEND + 1, nel - 1]),
                shift=dv)


def segs(d, b, w=None, part='main'):
    """第 b 条带在主路径（part='main'）或附加路径（part='tp'）上的线段，以及每段两端成分 w（dict 键）的平均。"""
    if part == 'main':
        x, E, c, n = d['x'][:KEND + 1], d['E'][:KEND + 1, b], d['c'], KEND + 1
    else:
        tp = d['tp']
        if tp is None or b >= tp['E'].shape[1]: return [], np.array([])
        x, E, c, n = tp['x'], tp['E'][:, b], tp['c'], len(tp['x'])
    S = [[(x[i], E[i]), (x[i + 1], E[i + 1])] for i in range(n - 1)]
    V = (c[w][:n - 1, b] + c[w][1:n, b]) / 2 if w else np.array([])
    return S, V


def axes_for(fig, cell, d):
    """一个子图：主路径一块；Sn 系右边另加一小块画附加路径 Γ–Σ′–R₂′–L（不在高对称路径上，单独画开，2026-11-04 用户意见）。"""
    if not d['tp']: return [(fig.add_subplot(cell), 'main')]
    sub = cell.subgridspec(1, 2, width_ratios=[d['x'][KEND] - d['x'][0], d['tp']['x'][-1]], wspace=0.12)
    a = fig.add_subplot(sub[0]); t = fig.add_subplot(sub[1], sharey=a)
    return [(a, 'main'), (t, 'tp')]


def frame(axs, d, ylim, title, n=None):
    """主路径：高对称点刻度、纵轴、标题；附加路径：浅紫底、Γ Σ′ R₂′ L 刻度、不重复纵轴。两块都画带边虚线。"""
    for ax, part in axs:
        ax.axhline(0, color='k', lw=0.5, ls=':', zorder=1); ax.axhline(d['gap'], color='k', lw=0.5, ls=':', zorder=1)
        if part == 'main':
            pf.kticks(ax, d['x']); ax.set_ylim(*ylim)
            ax.set_title(title, fontsize=8, pad=4); ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)', labelpad=2)
            if n is not None: pf.panel(ax, f'({"abcdefgh"[n]})')
        else:
            tp = d['tp']; lab = sorted(tp['labels'].items())
            ax.set_facecolor('#F1EEF7'); ax.set_xlim(0, tp['x'][-1])
            ax.set_xticks([tp['x'][i] for i, _ in lab], [pf.MATH.get(nm, nm) for _, nm in lab])
            for i, _ in lab[1:-1]: ax.axvline(tp['x'][i], color=pf.GRID, lw=0.5, zorder=0)
            ax.tick_params(axis='x', top=False, length=0, pad=3); ax.tick_params(axis='y', labelleft=False)
            ax.set_title('off path\n(mirror plane)', fontsize=6.2, pad=3, color='#444444', linespacing=1.0)


def figure_cells(nrow=2, ncol=2, h=5.4, ws=0.30):
    fig = plt.figure(figsize=(W2, h)); return fig, fig.add_gridspec(nrow, ncol, hspace=0.42, wspace=ws)


def band_rows(d, keys):
    """源数据：主路径与附加路径上每个 k 点、每条带的能量与成分（%）。"""
    K = pf.kpath(); nel = d['nel']; bs = range(nel - 24, min(nel + 16, d['E'].shape[1]))
    main = [[i] + [r(v, 6) for v in K[i]] + [r(d['x'][i], 5), pf.LAB.get(i, ''), pf.bandname(b, nel), r(d['E'][i, b], 5)]
            + [r(100 * d['c'][k][i, b], 2) for k in keys] for b in bs for i in range(KEND + 1)]
    tp = d['tp']; extra = []
    if tp:
        extra = [[i, r(tp['x'][i] / pf.TP_XSCALE, 5), tp['labels'].get(i, ''), pf.bandname(b, nel), r(tp['E'][i, b], 5)]
                 + [r(100 * tp['c'][k][i, b], 2) for k in keys] for b in bs if b < tp['E'].shape[1] for i in range(len(tp['x']))]
    return main, extra


def dump_bands(sub, d, keys):
    main, extra = band_rows(d, keys); m = d['m']
    pf.dump(sub, f'{m}_HSE06+SOC_path.csv', ['k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV']
            + [f'{k}_percent' for k in keys], main)
    if extra:
        pf.dump(sub, f'{m}_HSE06+SOC_true_valley_path_G-S-R2-L.csv', ['i', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV']
                + [f'{k}_percent' for k in keys], extra)


# ------------------------------------------------------------------ H1
def h1(out):
    fig, gs = figure_cells(); rows = []; every = []
    for n, m in enumerate(MATS):
        d = hse(m); nel = d['nel']; axs = axes_for(fig, gs[n // 2, n % 2], d); every += [a for a, _ in axs]
        for ax, part in axs:
            for b in range(nel - 24, min(nel + 16, d['E'].shape[1])):
                S, _ = segs(d, b, part=part)
                if S: ax.add_collection(LineCollection(S, colors='#B2182B' if b >= nel else '#1F3A93', linewidths=0.9, zorder=3))
        frame(axs, d, (-1.0, 1.0), f'{TEX[m]}     $E_g$ = {1000 * d["gap"]:.0f} meV', n); mark_valleys(axs, d)
        axs[0][0].set_yticks(np.arange(-1.0, 1.01, 0.5))
        rows.append([m, r(1000 * d['gap'], 1), r(1000 * d['gap_path'], 1), r(1000 * d['shift'], 1)])
        dump_bands('H1_bands', d, [])
    h = [Line2D([], [], color='#1F3A93', lw=1.0), Line2D([], [], color='#B2182B', lw=1.0), Line2D([], [], color='k', lw=0.6, ls=':')]
    fig.legend(h, ['valence bands', 'conduction bands', 'band edges (VBM = 0, CBM)'], loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.93))
    pf.dump('H1_bands', 'band_gaps_HSE06+SOC.csv', ['material', 'Eg_true_meV', 'Eg_on_path_meV', 'true_VBM_above_path_VBM_meV'], rows)
    pf.save(fig, out, 'H1_HSE_bands')
    return rows


# ------------------------------------------------------------------ H2、H4：单一成分的 fat band
CMAP = CM['cation_s']          # 所有成分用同一种颜色（2026-11-04 用户要求统一）


def fat(axs, d, key, ylim, title, wmax, n=None, cmap=None):
    norm = Normalize(0, wmax); nel = d['nel']; cmap = cmap or CMAP
    for ax, part in axs:
        for b in range(nel - 24, min(nel + 16, d['E'].shape[1])):
            S, V = segs(d, b, key, part)
            if not S: continue
            lc = LineCollection(S, cmap=cmap, norm=norm, linewidths=0.4 + 2.6 * np.clip(V / wmax, 0, 1), capstyle='round', zorder=3)
            lc.set_array(np.clip(V, 0, wmax)); ax.add_collection(lc)
    frame(axs, d, ylim, title, n)


def colorbar(fig, axs, wmax, label):
    """整张图一个普通色柱（两端不加箭头）；超过满刻度的按满刻度画。"""
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, wmax), cmap=CMAP), ax=axs, shrink=0.7, pad=0.02, aspect=30)
    t = np.arange(0, wmax + 1e-9, 0.1); cb.set_ticks(t, labels=[f'{100 * v:.0f}' for v in t]); cb.set_label(label)
    cb.ax.tick_params(direction='in', length=2, width=0.6, labelsize=6.5); cb.outline.set_linewidth(0.6)


def h2(out):
    fig, gs = figure_cells(); every = []
    for n, m in enumerate(MATS):
        d = hse(m); axs = axes_for(fig, gs[n // 2, n % 2], d); every += [a for a, _ in axs]
        fat(axs, d, 'cation_s', (-0.5, 0.5), f'{TEX[m]}     $E_g$ = {1000 * d["gap"]:.0f} meV', 0.30, n); mark_valleys(axs, d)
        dump_bands('H2_cation_s', d, ['cation_s'])
    colorbar(fig, every, 0.30, 'Cation-s weight (%)')
    pf.save(fig, out, 'H2_cation_s')


SITE_COL = {'A': pf.C_A, 'B': pf.C_B, 'Te': '#1A9850'}     # 按原子所在的位着色：A 位紫、B 位橙、Te 绿（与图 H3 一致）
AREA = 70                                                     # 圆点面积（pt²）= AREA × 该轨道的比例


def dots(ax, part, d, items):
    """圆点 fat band：先画细灰线（能带本身），再叠圆点，面积 ∝ 轨道比例。所有成分的点按面积从大到小画，小点压在大点上面，
    加白色细边，同一 k 点上几种轨道的贡献都看得见。"""
    nel = d['nel']; bs = range(nel - 24, min(nel + 16, d['E'].shape[1]))
    if part == 'main': x, E, c, sl = d['x'][:KEND + 1], d['E'][:KEND + 1], d['c'], slice(0, KEND + 1)
    else: tp = d['tp']; x, E, c, sl = tp['x'], tp['E'], tp['c'], slice(None)
    bs = [b for b in bs if b < E.shape[1]]
    for b in bs: ax.plot(x, E[:, b], color='#B5B5B5', lw=0.5, zorder=1)
    X, Y, S, C = [], [], [], []
    for key, col in items:
        W = c[key][sl]
        for b in bs:
            X.append(x); Y.append(E[:, b]); S.append(W[:, b]); C += [col] * len(x)
    X, Y, S = np.concatenate(X), np.concatenate(Y), np.concatenate(S); C = np.array(C)
    o = np.argsort(-S); o = o[S[o] > 0.02]
    ax.scatter(X[o], Y[o], s=AREA * np.clip(S[o], 0, 1), c=C[o], alpha=0.75, edgecolors='white', linewidths=0.25, zorder=2, rasterized=True)


def mark_valleys(axs, d):
    """在最高价带上标出各谷（▼ 加名称，名称与图 H5、H6 相同；同论文图 3、图 5）。附加路径上标 Σ′、R₂′。"""
    m, nel, x = d['m'], d['nel'], d['x']
    for ax, part in axs:
        if part == 'main':
            for nm, _E, kind, i in pf.levels(m, true=False):
                if kind == 'ref': continue
                lead = (x[90] - 0.28, 0.10) if nm == 'Δ' or (nm == 'Γ' and d['gap'] < 0.1) else None
                pf.vmark(ax, x[i], d['E'][i, nel - 1], nm, lead)
        else:
            for i, nm in d['tp']['labels'].items():
                if nm in ('Σ′', 'R₂′'): pf.vmark(ax, d['tp']['x'][i], d['tp']['E'][i, nel - 1], nm)


def h4(out):
    """一张图：每种材料一行，左 s（A 位 s、B 位 s 叠在同一张能带图上），右 p（A 位 p、B 位 p、Te p 叠在一起）。
    2026-11-04 用户要求：不同轨道的贡献画在同一张能带图里，s 单独一张。"""
    fig = plt.figure(figsize=(W2, 10.6)); gs = fig.add_gridspec(4, 2, hspace=0.42, wspace=0.30, top=0.93)
    for j, m in enumerate(MATS):
        d = hse(m); A_, B_ = m[:2], m[3:5]
        ns = lambda e: {'Pb': '6s', 'Sn': '5s', 'Sb': '5s', 'Bi': '6s'}[e]
        for col, (what, items) in enumerate((('s', [('A_s', SITE_COL['A']), ('B_s', SITE_COL['B'])]),
                                             ('p', [('Te_p', SITE_COL['Te']), ('B_p', SITE_COL['B']), ('A_p', SITE_COL['A'])]))):
            axs = axes_for(fig, gs[j, col], d)
            for ax, part in axs: dots(ax, part, d, items)
            title = (f'{TEX[m]}: s ({A_} {ns(A_)}, {B_} {ns(B_)})' if what == 's' else f'{TEX[m]}: p ({A_} p, {B_} p, Te p)')
            frame(axs, d, (-0.8, 0.8), title, 2 * j + col)
            if what == 's': mark_valleys(axs, d)
        dump_bands('H4_orbitals', d, ['cation_s', 'A_s', 'B_s', 'Te_s', 'A_p', 'B_p', 'Te_p'])
    hc = [Line2D([], [], ls='', marker='o', ms=6, color=SITE_COL['A'], alpha=0.6), Line2D([], [], ls='', marker='o', ms=6, color=SITE_COL['B'], alpha=0.6),
          Line2D([], [], ls='', marker='o', ms=6, color=SITE_COL['Te'], alpha=0.6)]
    hs = [Line2D([], [], ls='', marker='o', ms=np.sqrt(AREA * w), color='#777777', alpha=0.5) for w in (0.1, 0.3, 0.6)]
    fig.legend(hc + hs, ['A site (Pb, Sn)', 'B site (Sb, Bi)', 'Te', '10%', '30%', '60%'], loc='lower center', ncol=6, bbox_to_anchor=(0.5, 0.945),
               title='colour = atomic site;  dot area = orbital weight', title_fontsize=7, handletextpad=0.2, columnspacing=1.2)
    pf.save(fig, out, 'H4_orbitals')


SITE_CM = {'A': CM['A_p'], 'B': CM['B_p'], 'Te': CM['Te_p']}   # H4a–d：按原子所在的位着色的渐变（A 位紫、B 位橙、Te 绿）
W4 = 0.50                                                     # H4a–d 五个子图同一刻度，超过的按满刻度画


def h4_each(out):
    """H4a–d：每种材料一张 3×2，(a) A 位 s，(b) B 位 s，(c) A 位 p，(d) B 位 p，(e) Te p，(f) 色柱；s 子图标出价带谷。
    与总览图 H4 同一组数据，用户要求两种画法都保留（2026-11-04）。"""
    for j, m in enumerate(MATS):
        d = hse(m); A_, B_ = m[:2], m[3:5]
        ns = lambda e: {'Pb': '6s', 'Sn': '5s', 'Sb': '5s', 'Bi': '6s'}[e]
        fig = plt.figure(figsize=(W2, 7.9)); gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.30)
        panels = [('A_s', f'A-site s ({A_} {ns(A_)})', 'A'), ('B_s', f'B-site s ({B_} {ns(B_)})', 'B'),
                  ('A_p', f'A-site p ({A_})', 'A'), ('B_p', f'B-site p ({B_})', 'B'), ('Te_p', 'Te p', 'Te')]
        for n, (key, lab, site) in enumerate(panels):
            axs = axes_for(fig, gs[n // 2, n % 2], d)
            fat(axs, d, key, (-0.8, 0.8), lab, W4, n, SITE_CM[site])
            if key.endswith('_s'): mark_valleys(axs, d)
        al = fig.add_subplot(gs[2, 1]); al.axis('off')
        for i, (site, lab) in enumerate((('A', f'A site ({A_}): s and p'), ('B', f'B site ({B_}): s and p'), ('Te', 'Te: p'))):
            cax = al.inset_axes([0.12, 0.78 - 0.32 * i, 0.76, 0.08])
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, W4), cmap=SITE_CM[site]), cax=cax, orientation='horizontal')
            t = np.arange(0, W4 + 1e-9, 0.1); cb.set_ticks(t, labels=[f'{100 * v:.0f}' for v in t])
            cb.ax.tick_params(direction='in', length=2, width=0.6, labelsize=6.5); cb.outline.set_linewidth(0.6)
            cax.set_title(lab + ' — orbital weight (%)', fontsize=7, pad=2)
        fig.suptitle(f'{TEX[m]} (HSE06+SOC): orbital weights', fontsize=9, y=0.94)
        pf.save(fig, out, f'H4{"abcd"[j]}_orbitals_{m}')


# ------------------------------------------------------------------ H3：s 来自 A 位还是 B 位
def h3(out):
    norm = Normalize(0, 1); fig, gs = figure_cells(); every = []
    for n, m in enumerate(MATS):
        d = hse(m); nel = d['nel']; axs = axes_for(fig, gs[n // 2, n % 2], d); every += [a for a, _ in axs]
        for ax, part in axs:
            for b in range(nel - 24, min(nel + 16, d['E'].shape[1])):
                S, A = segs(d, b, 'A_s', part); _, Bs = segs(d, b, 'B_s', part)
                if not S: continue
                tot = A + Bs; w = np.clip(tot / 0.15, 0, 1)[:, None]
                cols = w * np.array([pf.SITE_CMAP(norm(a / max(t, 1e-9))) for a, t in zip(A, tot)]) + (1 - w) * pf.GREY
                ax.add_collection(LineCollection(S, colors=cols, linewidths=0.5 + 2.4 * np.clip(tot / 0.3, 0, 1), capstyle='round', zorder=3))
        frame(axs, d, (-0.5, 0.5), f'{TEX[m]}     A = {m[:2]}, B = {m[3:5]}', n); mark_valleys(axs, d)
        dump_bands('H3_s_by_site', d, ['A_s', 'B_s', 'cation_s'])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=pf.SITE_CMAP), ax=every, shrink=0.7, pad=0.02, aspect=30)
    cb.set_label('Share of cation-s on the A site'); cb.set_ticks([0, 0.5, 1], labels=['B only\n(Sb, Bi)', 'equal', 'A only\n(Pb, Sn)'])
    cb.ax.tick_params(direction='in', length=2, width=0.6, labelsize=6.5); cb.outline.set_linewidth(0.6)
    pf.save(fig, out, 'H3_s_by_site')


# ------------------------------------------------------------------ H5：四种材料的价带能级
def valley_list(m):
    """价带上的谷（HSE06+SOC，相对真实价带顶）：名称、种类、能量、k（分数）、等价口袋数、阳离子 s、来源、主轴有效质量。"""
    K = pf.kpath(); _, rots = pf.bz_lattice(m); nel = NEL[m]
    _, W = pf.cv.read_procar(os.path.join(pf.rf.base(m), 'c2_hse', 'PROCAR_OPT.gz')); cs = comp(W)['cation_s'][:, nel - 2:nel].mean(1)
    tv = {t['name']: t for t in (pf.true_valleys(m) or {'valleys': []})['valleys']}
    out = []
    for nm, E, kind, i in pf.levels(m):
        if nm in tv:
            t = tv[nm]
            out.append(dict(name=nm, kind=kind, E=E, k=t['kfrac'], n=len(pf.star(t['kfrac'], rots)), s=t['s'],
                            src='c4 mirror grid (not refined)' if t['grid'] else 'c5 star fit (true maximum)', mass=sorted(t['m'])))
        else:
            k = tuple(K[i if i is not None else 0]); ii = i if i is not None else 0
            out.append(dict(name=nm, kind=kind, E=E, k=k, n=len(pf.star(k, rots)), s=cs[ii], src='band path c2' + (' (Γ point, not a maximum)' if kind == 'ref' else ''),
                            mass=None))
    return out


def h5(out):
    fig, ax = plt.subplots(figsize=(W2, 3.6))
    X = {'Pb2Sb2Te5': 0.0, 'Pb2Bi2Te5': 1.25, 'Sn2Sb2Te5': 2.9, 'Sn2Bi2Te5': 4.15}; hw = 0.22; rows = []
    ax.axhspan(-2 * pf.KT300, 0, color='#EEF3F8', zorder=0, lw=0)
    for m in MATS:
        x0 = X[m]; V = valley_list(m); lab = []
        for v in V:
            col = pf.C_B if v['kind'] in ('axis', 'ref') else pf.C_A; ref = v['kind'] == 'ref'
            ax.plot([x0 - hw, x0 + hw], [v['E']] * 2, color=col, lw=0.9 if ref else 2.0, ls=(0, (2, 1.5)) if ref else '-', solid_capstyle='butt', zorder=3)
            path_lp = v['kind'] == 'lp' and v['mass'] is None
            txt = ('Γ' if ref else pf.MATH.get(v['name'], v['name']) + f'  {v["E"]:.0f}'.replace('-', '−') + f'  ×{v["n"]}' + (' †' if path_lp else ''))
            lab.append((v['E'], txt, col))
            rows.append([m, v['name'], v['kind'], r(v['E'], 1), v['n'], r(v['k'][0], 4), r(v['k'][1], 4), r(v['k'][2], 4), r(100 * v['s'], 1), v['src']]
                        + ([r(q, 3) for q in v['mass']] if v['mass'] else ['', '', '']))
        lab.sort(key=lambda t: -t[0]); ys = [t[0] for t in lab]
        for _ in range(60):                                   # 名称上下错开，至少相隔 11 meV
            for k in range(1, len(ys)):
                if ys[k - 1] - ys[k] < 11: mid = (ys[k - 1] + ys[k]) / 2; ys[k - 1], ys[k] = mid + 5.5, mid - 5.5
        for (e, txt, col), y in zip(lab, ys):
            ax.plot([x0 + hw, x0 + hw + 0.08], [e, y], color=col, lw=0.5, zorder=2)
            ax.text(x0 + hw + 0.11, y, txt, ha='left', va='center', fontsize=6.3, color='#222222')
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.set_xticks(list(X.values()), [TEX[m] for m in X]); ax.tick_params(axis='x', top=False, length=0)
    ax.set_xlim(-0.45, 5.35); ax.set_ylim(-265, 18); ax.set_yticks(range(-250, 1, 50)); ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (meV)')
    ax.axvline(2.25, color=pf.GRID, lw=0.6)
    for xm, t in ((0.62, 'Pb series'), (3.52, 'Sn series')):
        ax.text(xm, -0.09, t, ha='center', va='top', fontsize=7, transform=ax.get_xaxis_transform())
    ax.text(-0.42, -pf.KT300, r'$2k_\mathrm{B}T$ (300 K)', ha='left', va='center', fontsize=6.3, color='#555555')
    h = [Line2D([], [], color=pf.C_B, lw=2), Line2D([], [], color=pf.C_A, lw=2), Line2D([], [], color=pf.C_B, lw=0.9, ls=(0, (2, 1.5)))]
    ax.legend(h, ['valley on the Γ–A axis (cation/Te p)', 'lone-pair valley (cation s–Te p)', 'Γ point (reference, not a valley)'],
              loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6.5, handlelength=1.8, columnspacing=1.2)
    pf.dump('H5_valence_levels', 'valence_valleys_HSE06+SOC.csv',
            ['material', 'valley', 'kind (axis: on the Γ–A axis; lp: lone-pair valley; ref: Γ point, not a maximum)', 'E-EVBM_meV',
             'n_equivalent_pockets', 'k_frac_1', 'k_frac_2', 'k_frac_3', 'cation_s_percent', 'source', 'm1_me', 'm2_me', 'm3_me'], rows)
    pf.save(fig, out, 'H5_valence_levels')
    return rows


# ------------------------------------------------------------------ H6：Sn 系口袋在布里渊区中的位置
def h6(out):
    """两列（Sn₂Sb₂Te₅、Sn₂Bi₂Te₅）上下对齐：(a)(b) 沿 c 俯视，(c)(d) 三维透视（k_z 放大 4 倍），下方为每组谷的能量与等价口袋数。"""
    fig = plt.figure(figsize=(W2, 6.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.15, 0.22], wspace=0.08, hspace=0.06, bottom=0.07)
    rows = []
    for k, m in enumerate(SN):
        B, P = pf.pockets(m); hx, hy = pf.hexagon(B)
        ax = fig.add_subplot(gs[0, k])
        ax.plot(hx, hy, color='k', lw=0.8, zorder=1)
        r0 = 0.5 * np.linalg.norm(B[0][:2]); a0_ = np.arctan2(B[0][1], B[0][0])
        Ms = [(r0 * np.cos(a0_ + t * np.pi / 3), r0 * np.sin(a0_ + t * np.pi / 3)) for t in range(6)]
        for q in Ms: ax.plot([0, q[0]], [0, q[1]], color='#BDBDBD', lw=0.5, ls=pf.DASH, zorder=0)
        for g in P:
            big = pf.near(g['E'])
            for c, f in zip(g['cart'], g['frac']):
                ax.plot(c[0], c[1], pf.MK[g['name']], ms=(6.6 if big else 4.2) * pf.MS.get(pf.MK[g['name']], 1), mew=0.9, color=g['col'],
                        mfc='white' if pf.on_face(f) else g['col'], alpha=1.0 if big else 0.4, zorder=3)
                rows.append([m, g['name'], r(g['E'], 1), len(g['frac']), *[r(v, 4) for v in f], *[r(v, 4) for v in c], 'yes' if big else 'no'])
        Mq = max(Ms, key=lambda c: c[0] + 0.1 * c[1]); kv = max(zip(hx, hy), key=lambda t: (round(t[1], 3), t[0]))
        ax.text(Mq[0] + 0.04, Mq[1], 'M', ha='left', va='center', fontsize=7); ax.text(kv[0] + 0.03, kv[1] + 0.02, 'K', ha='left', va='bottom', fontsize=7)
        ax.text(0.03, -0.06, 'Γ', fontsize=7, ha='left', va='top')
        n = sum(len(g['frac']) for g in P if pf.near(g['E']))
        ax.set_title(f"{TEX[m]}: {n} pocket{'s' if n > 1 else ''} within $2k_\\mathrm{{B}}T$", fontsize=8, pad=4)
        ax.set_aspect('equal'); ax.set_xlim(-1.0, 1.0); ax.set_ylim(-0.92, 0.92); ax.axis('off')
        ax.text(0.0, 1.02, f'({"ab"[k]})', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='left', va='bottom')
        a3 = fig.add_subplot(gs[1, k], projection='3d'); pf.bz_3d(a3, m, B, P, zs=4)
        a3.text2D(0.5, 1.0, f'{TEX[m]}: 3D view ($k_z$ ×4)', transform=a3.transAxes, fontsize=8, ha='center', va='bottom')
        a3.text2D(0.0, 1.0, f'({"cd"[k]})', transform=a3.transAxes, fontsize=9, fontweight='bold', ha='left', va='bottom')
        al = fig.add_subplot(gs[2, k]); al.axis('off'); al.set_xlim(0, 1); al.set_ylim(0, 1)
        for j, g in enumerate(sorted(P, key=lambda g: -g['E'])):        # 每列一个小表：形状、名称、能量、等价口袋数（两栏）
            x0, y0 = (0.06 if j < 3 else 0.54), 0.90 - 0.36 * (j % 3); big = pf.near(g['E'])
            al.plot(x0, y0, pf.MK[g['name']], ms=5 * pf.MS.get(pf.MK[g['name']], 1), mew=0.8, color=g['col'], alpha=1 if big else 0.4,
                    mfc='white' if all(map(pf.on_face, g['frac'])) else g['col'], clip_on=False)
            al.text(x0 + 0.05, y0, pf.MATH.get(g['name'], g['name']), ha='left', va='center', fontsize=7)
            al.text(x0 + 0.42, y0, f"{g['E']:.0f}".replace('-', '−') + f" meV  ×{len(g['frac'])}", ha='right', va='center', fontsize=7)
    h = [Line2D([], [], ls='', marker='o', ms=5, color=pf.C_B, label='Γ / Δ: on the Γ–A axis'),
         Line2D([], [], ls='', marker='o', ms=5, color=pf.C_A, label='Σ′: lone-pair valley near Γ–M'),
         Line2D([], [], ls='', marker='s', ms=4.2, color=pf.C_A, label=r'R$_1$′: near A–L, close to A'),
         Line2D([], [], ls='', marker='p', ms=4.6, color=pf.C_A, label=r'R$_1$″: near A–L, closest to A'),
         Line2D([], [], ls='', marker='D', ms=3.8, color=pf.C_A, label=r'R$_2$′: near A–L, close to L'),
         Line2D([], [], ls='', marker='o', ms=3.6, color=pf.C_A, alpha=0.4, label=r'faded: more than $2k_\mathrm{B}T$ below the VBM')]
    fig.legend(handles=h, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.065), fontsize=6.5, handletextpad=0.3, columnspacing=1.2)
    pf.dump('H6_pockets', 'pocket_positions.csv', ['material', 'valley', 'E-EVBM_meV', 'n_equivalent_pockets', 'k_frac_1', 'k_frac_2', 'k_frac_3',
                                                  'kx_1/A', 'ky_1/A', 'kz_1/A', 'within_2kBT'], rows)
    pf.save(fig, out, 'H6_pockets_BZ')


# ------------------------------------------------------------------ H7：镜面内最高价带
def spline(m):
    from scipy.interpolate import RectBivariateSpline
    X, ws, GX = pf.mirror_map(m)
    we = np.concatenate([ws[-2:] - 1, ws, ws[:3] + 1]); Ge = np.concatenate([GX[:, -2:], GX, GX[:, :3]], axis=1)   # w 方向周期延拓
    return X, ws, GX, RectBivariateSpline(X, we, Ge, kx=3, ky=3, s=0)


def h7(out):
    fig, axs = plt.subplots(1, 2, figsize=(W2, 3.2), gridspec_kw=dict(wspace=0.30))
    lo = -250; rows = []
    for k, (ax, m) in enumerate(zip(axs, SN)):
        X, ws, GX, S = spline(m)
        xf, wf = np.linspace(X.min(), X.max(), 241), np.linspace(-0.5, 0.5, 201); F = S(xf, wf)
        im = ax.pcolormesh(xf, wf, np.clip(F.T, lo, 0), cmap='viridis', vmin=lo, vmax=0, shading='auto', rasterized=True)
        ax.contour(xf, wf, F.T, levels=np.arange(-300, 1, 25), colors='white', linewidths=0.35, alpha=0.7)
        for w0 in (-0.5, 0.0, 0.5): ax.axhline(w0, color='w', lw=1.6, ls=(0, (1, 1)))
        ax.axvline(0, color='w', lw=1.0, ls=(0, (1, 1)))
        ax.text(X.max() - 0.01, 0.015, 'Γ–M (path)', color='w', fontsize=6.3, ha='right', va='bottom')
        ax.text(0.20, 0.485, 'A–L (path)', color='w', fontsize=6.3, ha='center', va='top')
        ax.text(0.012, -0.30, 'Γ–A (path)', color='w', fontsize=6.3, ha='left', rotation=90)
        for nm, (u0, w0) in pf.PATH_LP.items():                       # 路径上的谷：空心
            for xx, ww in ((-u0, w0), (u0, -w0)):
                for wv in ((ww, ww - 1) if np.isclose(abs(ww), 0.5) else (ww,)):
                    if abs(wv) <= 0.5 + 1e-9:
                        ax.plot(xx, wv, pf.MK[nm], ms=5 * pf.MS.get(pf.MK[nm], 1), mfc='none', mec='k', mew=0.9, zorder=4, clip_on=False)
        for t in pf.true_valleys(m)['valleys']:                       # 镜面内的真实极大：红色实心
            u0, w0 = t['kfrac'][0], t['kfrac'][2]
            for xx, ww in ((-u0, w0), (u0, -w0)):
                ax.plot(xx, ww, pf.MK[t['name']], ms=5.2 * pf.MS.get(pf.MK[t['name']], 1), color='#FF5A5F', mec='k', mew=0.5, zorder=5)
            dx, dy, ha, va = {'R₂′': (0, -8, 'center', 'top'), 'R₁″': (6, -3, 'left', 'top')}.get(t['name'], (6, 4, 'left', 'bottom'))
            ax.annotate(f"{pf.MATH.get(t['name'], t['name'])} {t['E']:.0f}".replace('-', '−'), (-u0, w0), (dx, dy), textcoords='offset points',
                        fontsize=6.3, color='w', fontweight='bold', ha=ha, va=va)
        ax.set_xlim(X.min(), X.max()); ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel(r'in-plane $k$ along Γ→M ($-u$, units of $b_1$)'); ax.set_ylabel(r'$k_z$ ($w$, units of $b_3$)')
        ax.set_title(f'{TEX[m]}: top valence band in the mirror plane', fontsize=8, pad=3); pf.panel(ax, f'({"ab"[k]})')
        pf.dump('H7_mirror_plane', f'{m}_mirror_plane_grid_meV.csv', ['x=-u'] + [f'w={w:.4f}' for w in ws],
                [[round(xx, 4)] + [round(v, 1) for v in row] for xx, row in zip(X, GX)])
        rows += [[m, t['name'], r(-t['kfrac'][0], 4), r(t['kfrac'][2], 4), r(t['E'], 1)] for t in pf.true_valleys(m)['valleys']]
        rows += [[m, nm + ' (path)', r(-u0, 4), r(w0, 4), r(float(S(np.array([-u0]), np.array([w0]))[0, 0]), 1)] for nm, (u0, w0) in pf.PATH_LP.items()]
    cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02); cb.set_label(r'$E - E_\mathrm{VBM}$ (meV)')
    pf.dump('H7_mirror_plane', 'valleys_in_mirror_plane.csv', ['material', 'valley', 'x=-u', 'w', 'E-EVBM_meV'], rows)
    pf.dump('H7_mirror_plane', 'path_vs_true_maxima.csv', ['material', 'path_valley', 'true_maximum', 'E_path_meV', 'E_true_meV', 'rise_meV',
                                                           'k_distance_1/A', 'path_x=-u', 'path_w', 'true_x=-u', 'true_w'], sum((path_vs_true(m) for m in SN), []))
    pf.save(fig, out, 'H7_mirror_plane_maps')


ORDER = ('Σ′', 'R₁′', 'R₁″', 'R₂′')


def path_vs_true(m):
    """每个真实谷与它所来自的路径谷：路径值（c2_hse）、真实极大（c5/c6）、差值、k 空间距离。图 H8 已按用户意见删去，数据并入图 H7。"""
    tv = pf.true_valleys(m); v = pf.top_vb(m, 'c2_hse') - tv['shift']
    B, _ = pf.bz_lattice(m); nb1, cz = np.linalg.norm(B[0]), abs(B[2, 2]); out = []
    for t in sorted(tv['valleys'], key=lambda t: ORDER.index(t['name'])):
        pn = t['name'].rstrip('′″'); Ep, Et = v[pf.SNVAL[pn]], t['E']; up, wp = pf.PATH_LP[pn]; ut, wt = t['kfrac'][0], t['kfrac'][2]
        cand = [(up, wp + n) for n in (-1, 0, 1)] + [(-up, -wp + n) for n in (-1, 0, 1)]
        dk = min(float(np.hypot((c[0] - ut) * nb1, (c[1] - wt) * cz)) for c in cand)
        out.append([m, pn, t['name'], r(Ep, 1), r(Et, 1), r(Et - Ep, 1), r(dk, 4), r(-up, 4), r(wp, 4), r(-ut, 4), r(wt, 4)])
    return out


# ------------------------------------------------------------------ H8：Γ 点能级与 SOC
def h8(out):
    """有 HSE 不加 SOC 的结果（SOC机制检验_HSE/e1_hse_nosoc_<材料>）就画 HSE，否则画 PBE（图名后缀 _HSE / _PBE）。"""
    import soc_mechanism as sm
    fn = 'hse' if all(os.path.exists(os.path.join(sm.R, 'SOC机制检验_HSE', f'e1_hse_nosoc_{m}', 'PROCAR_OPT.gz')) for m in SN) else 'pbe'
    res = {}
    for key, (d, comp_) in sm.runs(fn).items():
        q = sm.gamma_states(d, comp_) if os.path.isdir(d) else None
        if q: res[key] = q
    dec = {}
    for A_ in ('Pb', 'Sn'):
        sb, bi = f'{A_}2Sb2Te5', f'{A_}2Bi2Te5'
        if not all((mm, mm, s_) in res for mm in (sb, bi) for s_ in (True, False)): continue
        lo = {mm: res[(mm, mm, True)]['p12'] - res[(mm, mm, False)]['pc'] for mm in (sb, bi)}
        dec[A_] = dict(total=res[(bi, bi, True)]['p12'] - res[(sb, sb, True)]['p12'], scalar=res[(bi, bi, False)]['pc'] - res[(sb, sb, False)]['pc'],
                       soc=lo[bi] - lo[sb], low_sb=lo[sb], low_bi=lo[bi])
    tag = 'HSE' if fn == 'hse' else 'PBE'
    d0, s0 = pf.dump, pf.save
    pf.dump = lambda sub, name, header, rows: d0('H8_Gamma_levels_SOC' if sub == 'FigS2_Gamma_levels_SOC' else sub, name, header, rows)
    pf.save = lambda fig, o, name: s0(fig, o, f'H8_Gamma_levels_SOC_{tag}' if name.startswith('FigS2_Gamma_levels_SOC') else name)
    try:
        data = pf.DATA; sm.figure(res, dec, out, fn); pf.DATA = data
    finally:
        pf.dump, pf.save = d0, s0
    rows = []                                                   # 加 SOC 时 PBE 与 HSE 的 Γ 点（表 4 的同一判据）
    for mm in sm.NELECT:
        for f_, side in (('PBE+SOC', 'c2_pbe'), ('HSE06+SOC', 'c2_hse')):
            d = os.path.join(sm.R, side) if mm == 'Pb2Sb2Te5' else os.path.join(sm.R, mm, side)
            q = sm.gamma_states(d, mm)
            rows.append([mm, f_, r(q['p12'], 1), r(q['s'], 1), q['s_in']])
    pf.dump('H8_Gamma_levels_SOC', 'Gamma_SOC_levels_PBE_vs_HSE.csv',
            ['material', 'functional', 'Gamma_VB_p1/2_minus_Sigma_meV', 'Gamma_s_state_minus_Sigma_meV', 's_state_in'], rows)
    return fn


# ------------------------------------------------------------------ 表格（专题报告里的表直接由这里的 CSV 生成，首行为表头）
def tables(out):
    K = pf.kpath(); key = []
    keys = ('A_s', 'B_s', 'Te_s', 'A_p', 'B_p', 'Te_p')
    for m in MATS:
        d = hse(m); nel = d['nel']; tp = d['tp']; vb, cb = [nel - 2, nel - 1], [nel, nel + 1]
        def at(i, bb, on_tp=False, name=''):
            src = tp if on_tp else d; E = src['E']; c = src['c']
            kk = tuple(tp['k'][i]) if on_tp else tuple(K[i])
            where = name or pf.kname(i)
            return [where, '(' + ', '.join(f'{v:.3f}'.replace('-', '−') for v in kk) + ')', r(1000 * E[i, bb].mean(), 1)] + [r(100 * c[q][i, bb].mean(), 1) for q in keys]
        # 价带顶
        iv = int(np.nanargmax(d['E'][:KEND + 1, nel - 1]))
        if tp and np.nanmax(tp['E'][:, nel - 1]) > d['E'][iv, nel - 1] + 5e-4:
            j = int(np.nanargmax(tp['E'][:, nel - 1])); row_v = at(j, vb, True, tp['labels'].get(j, 'off path') + '（附加路径）')
        else:
            ia = max(range(90, 98), key=lambda i: d['E'][i, nel - 1]); row_v = at(ia if d['E'][ia, nel - 1] > d['E'][iv, nel - 1] - 5e-4 else iv, vb)
        # 最高的孤对电子价带谷
        lp = max((v for v in valley_list(m) if v['kind'] == 'lp'), key=lambda v: v['E'])
        jl = next((i for i, nm in (tp['labels'].items() if tp else []) if nm == lp['name']), None)
        if jl is not None: row_l = at(jl, vb, True, lp['name'] + '（附加路径）')
        else:
            il = next(i for nm, E, kind, i in pf.levels(m) if nm == lp['name']); row_l = at(il, vb, name=lp['name'])
        # 导带底
        ic = int(np.nanargmin(d['E'][:KEND + 1, nel]))
        ia = min(range(90, 98), key=lambda i: d['E'][i, nel]); row_c = at(ia if d['E'][ia, nel] < d['E'][ic, nel] + 5e-4 else ic, cb)
        for st, row in (('价带顶', row_v), ('最高的孤对电子价带谷', row_l), ('导带底', row_c)):
            key.append([TEX[m].replace('$_2$', '₂').replace('$_5$', '₅'), st] + row)
    pf.dump('tables', 'T2_key_states.csv', ['材料', '态', '位置', 'k（分数坐标）', 'E − E_VBM (meV)', 'A 位 s (%)', 'B 位 s (%)', 'Te s (%)',
                                            'A 位 p (%)', 'B 位 p (%)', 'Te p (%)'], key)
    nm = lambda m: m.replace('2', '₂').replace('5', '₅')
    kind = {'axis': 'Γ–A 轴上的谷', 'lp': '孤对电子谷', 'ref': 'Γ 点（不是极大）'}
    rows = []
    for m in MATS:
        for v in valley_list(m):
            rows.append([nm(m), v['name'], kind[v['kind']], r(v['E'], 1), v['n'], '(' + ', '.join(f'{q:.3f}'.replace('-', '−') for q in v['k']) + ')',
                         r(100 * v['s'], 1), ' / '.join(f'{q:.3f}' for q in v['mass']) if v['mass'] else '—',
                         {'c5 star fit (true maximum)': '镜面搜索 + 星形拟合（真实极大）', 'c4 mirror grid (not refined)': '镜面网格点（未细化）'}.get(v['src'], '能带路径')])
    pf.dump('tables', 'T3_valence_valleys.csv', ['材料', '谷', '类型', 'E − E_VBM (meV)', '等价口袋数', 'k（分数坐标）', '阳离子 s (%)',
                                                 '主轴有效质量 (mₑ)', '来源'], rows)
    import soc_mechanism as sm
    g = []
    for m in sm.NELECT:
        for fn, side in (('PBE+SOC', 'c2_pbe'), ('HSE06+SOC', 'c2_hse')):
            q = sm.gamma_states(os.path.join(sm.R, side) if m == 'Pb2Sb2Te5' else os.path.join(sm.R, m, side), m)
            g.append([nm(m), fn, int(round(q['p12'])), int(round(q['s'])), q['s_in']])
    pf.dump('tables', 'T4_Gamma_levels.csv', ['材料', '泛函', 'Γ 处最高价带（p₁/₂ 型）− Σ 谷 (meV)', 'Γ 处 s 态 − Σ 谷 (meV)', 's 态位于'], g)
    return key


FIGS = {'H1': h1, 'H2': h2, 'H3': h3, 'H4': h4, 'H4x': h4_each, 'H5': h5, 'H6': h6, 'H7': h7, 'H8': h8, 'T': tables}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(pf.ROOT, '21_阶段报告与交接', 'HSE能带_20261104', '图'))
    ap.add_argument('--only', nargs='+'); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); pf.style(); pf.DATA = os.path.join(a.out, 'source_data')
    for k in a.only or FIGS:
        print(k, FIGS[k](a.out) or '')
    print('→', a.out)


if __name__ == '__main__':
    main()
