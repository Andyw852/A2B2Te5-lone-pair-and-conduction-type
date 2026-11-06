#!/usr/bin/env python3
"""能带“粗细渐变 + 颜色渐变”（fat band）独立绘图脚本 —— 报告图 10-1（以及图 3-1）的画法。

只依赖 numpy 与 matplotlib；数据在同目录的 data/<泛函>/ 下（由仓库 10_脚本/export_fatband_kit.py 从计算结果导出）。

画法要点
  1. 每条能带拆成相邻 k 点之间的小线段，用 matplotlib.collections.LineCollection 一次画出。
  2. 每段的权重 = 两端阳离子 s 成分的平均；同一个权重同时决定
       颜色：映射到单色渐变色图（浅蓝 → 深蓝），Normalize(0, 满刻度)；
       线宽：lw_min + lw_add × min(权重 / 满刻度, 1)。
     所以权重越大，线越粗、颜色越深；超过满刻度（默认 30%）的按满刻度画。
  3. capstyle='round' 让相邻两段线宽不同时接口处不留缺口；--refine N 把每段再细分 N 份（线性内插能量与权重），渐变更平滑。
  4. Sn 系另有不在高对称路径上的附加路径 Γ–Σ′–R₂′–L，单独画在右侧浅紫底小图，横向放大 1.8 倍，与主图共用纵轴。
  5. 水平虚线为 300 K 的化学势 µ：橙色为样品霍尔浓度，灰色为同类型 3×10¹⁹ cm⁻³。

用法
  python plot_fatband_gradient.py                               # HSE06+SOC 四种材料（即报告图 10-1）
  python plot_fatband_gradient.py --functional PBE+SOC          # PBE+SOC（截图里的旧版图 10）
  python plot_fatband_gradient.py --mats Sn2Sb2Te5 Sn2Bi2Te5    # 只画两种材料（一行两个子图）
  python plot_fatband_gradient.py --refine 4 --weight-sat 25 --lw-min 0.4 --lw-add 3.0
  python plot_fatband_gradient.py --no-mu                       # 不画化学势，即图 3-1 的样式
  python plot_fatband_gradient.py --colors "#FDE0DD,#FA9FB5,#C51B8A,#7A0177"   # 换一个渐变色
输出：<--out>.png 与 <--out>.pdf（默认 fatband_<泛函>）
"""
import argparse, csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
TEX = {m: m.replace('2', '$_2$').replace('5', '$_5$') for m in MATS}
COLORS = '#DCE8F4,#9CC3E4,#4A90C8,#1F5FA0,#08306B'     # 报告用的单色渐变（浅 → 深）
C_HALL, C_REF = '#E69F00', '#555555'                   # 霍尔浓度、3×10¹⁹ 的 µ 虚线
DASH = (0, (3, 1.6))
OFF_BG, GRID = '#F1EEF7', '#D9D9D9'
OFF_XSCALE = 1.8                                       # 附加路径横向放大倍数（与报告相同）


def style():
    """Times New Roman（没有时用度量相同的 Liberation Serif），公式 STIX，期刊尺寸的字号与刻度。"""
    for f in font_manager.findSystemFonts():
        n = os.path.basename(f).lower().replace('-', '').replace(' ', '').replace('_', '')
        if n.startswith(('times', 'liberationserif')): font_manager.fontManager.addfont(f)
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman', 'Liberation Serif', 'DejaVu Serif'], 'mathtext.fontset': 'stix',
        'font.size': 7.5, 'axes.labelsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'axes.linewidth': 0.8, 'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'legend.frameon': False,
        'savefig.dpi': 600, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02, 'pdf.fonttype': 42, 'axes.unicode_minus': True})


def read_bands(path):
    """长表 CSV（每行一个 k 点、一条带）→ x（k 距离，Å⁻¹）、{序号: 高对称点名}、[(带名, E, 权重)]（权重为 0–1）。"""
    rows = list(csv.DictReader(open(path, encoding='utf-8-sig')))
    kcol = 'k_index' if 'k_index' in rows[0] else 'i'
    order, data = [], {}
    for r in rows:
        b = r['band']
        if b not in data: order.append(b); data[b] = []
        data[b].append((int(r[kcol]), float(r['k_distance_1/A']), r['label'], float(r['E-EVBM_eV']), float(r['cation_s_percent']) / 100))
    first = sorted(data[order[0]])
    x = np.array([p[1] for p in first]); labels = {i: p[2] for i, p in enumerate(first) if p[2]}
    bands = [(b, np.array([p[3] for p in sorted(data[b])]), np.array([p[4] for p in sorted(data[b])])) for b in order]
    return x, labels, bands


def gradient_line(ax, x, y, w, cmap, norm, wsat, lw_min, lw_add, refine=1, zorder=3):
    """一条能带的双渐变线：颜色与线宽都随权重 w 变化。"""
    if refine > 1:                                     # 每段细分 refine 份，能量与权重线性内插
        t = np.linspace(0, 1, refine + 1)[:-1]
        f = lambda a: np.concatenate([a[i] + (a[i + 1] - a[i]) * t for i in range(len(a) - 1)] + [a[-1:]])
        x, y, w = f(x), f(y), f(w)
    pts = np.column_stack([x, y]); segs = np.stack([pts[:-1], pts[1:]], axis=1)
    c = (w[:-1] + w[1:]) / 2                           # 每段取两端的平均
    lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=lw_min + lw_add * np.clip(c / wsat, 0, 1),
                        capstyle='round', joinstyle='round', zorder=zorder)
    lc.set_array(np.clip(c, 0, wsat)); ax.add_collection(lc)


def sci(N):
    a, e = f'{abs(N):.1e}'.split('e'); return rf'{a}$\times10^{{{int(e)}}}$'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--functional', default='HSE06+SOC', choices=['HSE06+SOC', 'PBE+SOC'])
    ap.add_argument('--data', help='数据目录（默认 data/<泛函>）')
    ap.add_argument('--mats', nargs='+', default=MATS, choices=MATS)
    ap.add_argument('--out', help='输出文件名（不含扩展名），默认 fatband_<泛函>')
    ap.add_argument('--weight-sat', type=float, default=30, help='满刻度（%%）：权重达到它时颜色最深、线最粗，默认 30')
    ap.add_argument('--lw-min', type=float, default=0.5, help='权重为 0 时的线宽（pt），默认 0.5')
    ap.add_argument('--lw-add', type=float, default=2.4, help='满刻度时在 lw-min 之上增加的线宽（pt），默认 2.4')
    ap.add_argument('--colors', default=COLORS, help='渐变色（逗号分隔的十六进制，浅 → 深）')
    ap.add_argument('--refine', type=int, default=1, help='每段细分份数，>1 时渐变更平滑，默认 1（与报告相同）')
    ap.add_argument('--ylim', type=float, nargs=2, default=(-0.62, 0.52))
    ap.add_argument('--no-mu', action='store_true', help='不画化学势虚线')
    ap.add_argument('--no-offpath', action='store_true', help='不画右侧附加路径小图')
    a = ap.parse_args()
    style()
    d = a.data or os.path.join(HERE, 'data', a.functional); wsat = a.weight_sat / 100
    cmap = LinearSegmentedColormap.from_list('grad', a.colors.split(',')); norm = Normalize(0, wsat)
    gap = {r['material']: float(r['Eg_meV (incl. off path)']) for r in csv.DictReader(open(os.path.join(d, 'band_gap.csv'), encoding='utf-8-sig'))}
    mu = {}
    for r in csv.DictReader(open(os.path.join(d, 'chemical_potential_300K.csv'), encoding='utf-8-sig')):
        mu.setdefault(r['material'], []).append((r['case'], float(r['carrier_conc_signed_cm-3']), float(r['mu-EVBM_eV'])))

    n = len(a.mats); nc = 1 if n == 1 else 2; nr = (n + 1) // 2
    fig = plt.figure(figsize=(7.0 if nc == 2 else 3.6, 2.6 * nr)); gs = fig.add_gridspec(nr, nc, hspace=0.38, wspace=0.30); every = []
    for k, m in enumerate(a.mats):
        x, lab, bands = read_bands(os.path.join(d, f'{m}_main_path.csv'))
        fo = os.path.join(d, f'{m}_off_path_G-S-R2-L.csv'); off = None if (a.no_offpath or not os.path.exists(fo)) else read_bands(fo)
        cell = gs[k // nc, k % nc]
        if off:                                        # 主图 + 右侧小图，宽度按各自的 k 长度
            sub = cell.subgridspec(1, 2, width_ratios=[x[-1] - x[0], OFF_XSCALE * off[0][-1]], wspace=0.12)
            ax = fig.add_subplot(sub[0]); at = fig.add_subplot(sub[1], sharey=ax)
        else:
            ax, at = fig.add_subplot(cell), None
        every += [q for q in (ax, at) if q is not None]
        for _, E, w in bands: gradient_line(ax, x, E, w, cmap, norm, wsat, a.lw_min, a.lw_add, a.refine)
        if off:
            xo = OFF_XSCALE * off[0]
            for _, E, w in off[2]: gradient_line(at, xo, E, w, cmap, norm, wsat, a.lw_min, a.lw_add, a.refine)
        if not a.no_mu:
            for case, N, u in mu.get(m, []):
                hall = case == 'Hall'
                for q in (ax, at):
                    if q is not None: q.axhline(u, color=C_HALL if hall else C_REF, lw=1.4 if hall else 1.0, ls=DASH, zorder=4)
                up = (N < 0) == hall                    # 文字放在虚线上方还是下方（与报告相同）
                ax.text(0.985, u + (0.012 if up else -0.012), f'{"$p$" if N > 0 else "$n$"} = {sci(N)}' + (' (Hall)' if hall else ''),
                        transform=ax.get_yaxis_transform(), ha='right', va='bottom' if up else 'top', fontsize=6.3, zorder=6,
                        bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.85))
        # 主图：高对称点刻度、竖线、纵轴
        ticks = sorted(lab); ax.set_xticks([x[i] for i in ticks], [lab[i] for i in ticks])
        for i in ticks[1:-1]: ax.axvline(x[i], color=GRID, lw=0.5, zorder=0)
        ax.set_xlim(x[0], x[-1]); ax.set_ylim(*a.ylim); ax.set_yticks(np.arange(-0.5, a.ylim[1] + 1e-9, 0.25))
        ax.tick_params(axis='x', top=False, length=0, pad=3); ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)', labelpad=2)
        ax.set_title(f'{TEX[m]}     $E_g$ = {gap[m]:.0f} meV', fontsize=8, pad=4)
        ax.text(-0.02, 1.02, f'({"abcd"[k]})', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='right', va='bottom')
        if off:                                        # 右侧小图：浅紫底、Γ Σ′ R₂′ L 刻度、不重复纵轴
            lo = sorted(off[1]); at.set_facecolor(OFF_BG); at.set_xlim(0, xo[-1])
            at.set_xticks([xo[i] for i in lo], [off[1][i].replace('R₂′', r'R$_2$′') for i in lo])
            for i in lo[1:-1]: at.axvline(xo[i], color=GRID, lw=0.5, zorder=0)
            at.tick_params(axis='x', top=False, length=0, pad=3); at.tick_params(axis='y', labelleft=False)
            at.set_title('off path\n(mirror plane)', fontsize=6.2, pad=3, color='#444444', linespacing=1.0)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=every, shrink=0.7, pad=0.015, aspect=30)
    t = np.arange(0, wsat + 1e-9, 0.1); cb.set_ticks(t, labels=[f'{100 * v:.0f}' for v in t]); cb.set_label('Cation-s weight (%)')
    cb.ax.tick_params(direction='in', length=2, width=0.6); cb.outline.set_linewidth(0.6)
    if not a.no_mu:
        fig.legend([Line2D([], [], color=C_HALL, lw=1.4, ls=DASH), Line2D([], [], color=C_REF, lw=1.0, ls=DASH)],
                   [r'$\mu$ at the sample Hall concentration (300 K)', r'$\mu$ at $3\times10^{19}$ cm$^{-3}$, same carrier type'],
                   loc='lower center', ncol=2, bbox_to_anchor=(0.45, 0.925 if nr == 2 else 0.97), columnspacing=1.6, handlelength=2.6)
    out = a.out or os.path.join(HERE, f'fatband_{a.functional}')
    for ext in ('png', 'pdf'): fig.savefig(f'{out}.{ext}')
    print('→', out + '.png / .pdf')


if __name__ == '__main__':
    main()
