#!/usr/bin/env python3
"""讲解图：为什么 Pb₂Bi₂Te₅ 样品的 |S| 最高（2026-11-11，用户问 n 型 Pb₂Bi₂Te₅ 为何性能最好）。

(a) 同一浓度下四种材料的 S 差不多：PBE+SOC BoltzTraP2 在三个相同浓度（Pb₂Bi₂Te₅ 样品的 n 8.1×10¹⁹、同浓度的 p 型、
    另外三个样品附近的 p 9.1×10²⁰）下的 300 K |S|。组内四根柱子差不多，组间相差几倍：S 的差别来自浓度。
(b) 浓度低为什么 S 大：Pb₂Bi₂Te₅ 的 PBE+SOC 能带，样品浓度（8.1×10¹⁹）与假想的 9.1×10²⁰ 下 300 K 的化学势 µ，
    阴影为 µ ± 3k_BT（导电的能量窗口）。µ 离导带底越远，窗口两侧的贡献越对称、越抵消，S 越小。
(c) 为什么只有它浓度低：四个样品的室温霍尔浓度排成 2×2（A 位 × B 位）。含 Sn 或含 Sb 的三个都是 ~10²¹ 的 p 型，
    只有两者都不含的 Pb₂Bi₂Te₅ 是低浓度 n 型。原因是推测（母体二元化合物的已知规律），本工作尚未对 Pb₂Bi₂Te₅ 做缺陷计算。
数值全部来自已有计算（bt2_pbe_dense、c2_pbe）与实验数据，不做新的估计。
输出：21_阶段报告与交接/讲解图_20261111/讲解_Pb2Bi2Te5的S为何最高.png、.pdf 与 data/*.csv
用法：python 10_脚本/explain_pbbi_high_S.py
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import transport_figures as tf
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm, LinearSegmentedColormap

OUT = os.path.join(pf.ROOT, '21_阶段报告与交接', '讲解图_20261111')
NAME = {'Pb2Sb2Te5': 'Pb₂Sb₂Te₅', 'Pb2Bi2Te5': 'Pb₂Bi₂Te₅', 'Sn2Sb2Te5': 'Sn₂Sb₂Te₅', 'Sn2Bi2Te5': 'Sn₂Bi₂Te₅'}
KT = 1000 * pf.mv.KB * 300


def style():
    pf.style()
    for f in font_manager.findSystemFonts():
        if 'wqy' in f.lower() or 'wenquanyi' in f.lower(): font_manager.fontManager.addfont(f)
    have = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams['font.family'] = [n for n in ('Times New Roman', 'Liberation Serif') if n in have][:1] + ['WenQuanYi Zen Hei', 'DejaVu Serif']


def dump(name, header, rows):
    os.makedirs(os.path.join(OUT, 'data'), exist_ok=True)
    with open(os.path.join(OUT, 'data', name), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def main():
    style(); os.makedirs(OUT, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 2.75)); gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.95, 0.95], wspace=0.42)

    # (a) 同一浓度下四种材料的 |S|（PBE+SOC，300 K）：每组内四根柱子差不多，组与组之间（浓度）差得多
    ax = fig.add_subplot(gs[0]); groups = ((-8.11e19, 'n 型\n8.1×10¹⁹'), (8.11e19, 'p 型\n8.1×10¹⁹'), (9.07e20, 'p 型\n9.1×10²⁰'))
    rows = []; wd = 0.19
    for g, (N, lab) in enumerate(groups):
        for k, m in enumerate(pf.MATS):
            T, f = tf.LOAD['PBE'](m); v = abs(f(N, list(T).index(300)))
            ax.bar(g + (k - 1.5) * (wd + 0.01), v, wd, color=pf.MCOL[m], zorder=3, label=NAME[m] if g == 0 else None)
            rows.append([m, 'n' if N < 0 else 'p', f'{abs(N):.3e}', round(v, 1)])
    ax.set_xticks(range(3), [g[1] for g in groups], fontsize=6.6); ax.tick_params(axis='x', length=0, top=False)
    ax.set_ylim(0, 175); ax.set_ylabel('|S|，300 K（μV/K）'); ax.grid(axis='y', color=pf.GRID, lw=0.4, zorder=0)
    ax.text(0, 112, 'Pb₂Bi₂Te₅ 样品\n就在这个浓度', ha='center', fontsize=5.8, color='#333333')
    ax.text(2, 47, '另外三个样品\n在 ~10²¹', ha='center', fontsize=5.8, color='#333333')
    ax.set_xlabel('相同的载流子浓度（cm⁻³）')
    ax.set_title('(a) 同一浓度下，四种材料的 S 差不多', fontsize=7.6, pad=4, loc='left')
    ax.legend(fontsize=5.6, loc='upper right', handlelength=1.0, labelspacing=0.2, ncol=1)
    ax.text(0.98, 0.60, 'PBE+SOC 计算', transform=ax.transAxes, fontsize=5.6, color='#555555', ha='right')
    dump('a_S_same_concentration_PBE_300K.csv', ['material', 'carrier_type', 'concentration_cm-3', 'abs_S_PBE+SOC_uV/K'], rows)
    dump('a_samples.csv', ['material', 'type', 'Hall_conc_cm-3', 'abs_S_exp_300K_uV/K'],
         [[m, 'p' if pf.HALL[m]['N'] > 0 else 'n', f"{abs(pf.HALL[m]['N']):.3e}", round(abs(pf.EXP[m][0][1]), 2)] for m in pf.MATS])

    # (b) Pb₂Bi₂Te₅ 能带上两种浓度的 µ
    ax = fig.add_subplot(gs[1]); B = tf.bands_fn('Pb2Bi2Te5', 'PBE'); x, E, nel = B['x'], B['E'], B['nel']; K = pf.KEND
    for b in range(nel - 6, nel + 6):
        ax.plot(x[:K + 1], E[:K + 1, b], color='#B2182B' if b >= nel else '#1F3A93', lw=0.8, zorder=2)
    T, f = tf.LOAD['PBE']('Pb2Bi2Te5'); it = list(T).index(300); brow = []
    for N, col, tag in ((-8.11e19, '#E69F00', '样品 8.1×10¹⁹'), (-9.07e20, '#7A7A7A', '假如 9.1×10²⁰')):
        u = B['mu'](N); dep = 1000 * (u - B['gap']); S = f(N, it)
        ax.axhspan(u - 3 * KT / 1000, u + 3 * KT / 1000, color=col, alpha=0.18, lw=0, zorder=1)
        ax.axhline(u, color=col, lw=1.3, ls=pf.DASH, zorder=3)
        ax.text(0.03, u + 3 * KT / 1000 + 0.01, f'{tag}：S = {S:.0f} μV/K\nµ 在导带底上方 {dep:.0f} meV（{dep / KT:.0f} kBT）'.replace('-', '−'),
                transform=ax.get_yaxis_transform(), ha='left', va='bottom', fontsize=5.6, color='#222222', linespacing=1.1,
                bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.9), zorder=5)
        brow.append([f'{abs(N):.2e}', round(u, 4), round(dep, 1), round(dep / KT, 2), round(S, 1)])
    ax.axhline(B['gap'], color='k', lw=0.5, ls=':'); ax.axhline(0, color='k', lw=0.5, ls=':')
    pf.kticks(ax, x); ax.set_ylim(-0.25, 0.80); ax.set_ylabel(r'$E - E_\mathrm{VBM}$ (eV)')
    ax.set_title('(b) 浓度低 → µ 离带边近 → S 大', fontsize=7.6, pad=4, loc='left')
    ax.text(0.02, 0.02, 'Pb₂Bi₂Te₅，PBE+SOC；阴影为 µ ± 3kBT', transform=ax.transAxes, fontsize=5.6, color='#555555', va='bottom',
            bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.85))
    dump('b_Pb2Bi2Te5_mu_300K_PBE.csv', ['electron_conc_cm-3', 'mu-EVBM_eV', 'mu_above_CBM_meV', 'mu_above_CBM_kBT', 'S_PBE+SOC_uV/K'], brow)

    # (c) 2×2：A 位 × B 位的霍尔浓度
    ax = fig.add_subplot(gs[2]); grid = [['Sn2Sb2Te5', 'Sn2Bi2Te5'], ['Pb2Sb2Te5', 'Pb2Bi2Te5']]
    cmap = LinearSegmentedColormap.from_list('n', ['#F2F0F7', '#BCBDDC', '#756BB1', '#3F007D']); norm = LogNorm(5e19, 2e21); crow = []
    for i, r in enumerate(grid):
        for j, m in enumerate(r):
            N = pf.HALL[m]['N']; dark = abs(N) > 4e20
            ax.add_patch(plt.Rectangle((j, 1 - i), 0.96, 0.96, color=cmap(norm(abs(N))), lw=0))
            ax.text(j + 0.48, 1 - i + 0.66, NAME[m], ha='center', va='center', fontsize=7, color='white' if dark else 'k', fontweight='bold')
            ax.text(j + 0.48, 1 - i + 0.36, f'{"p" if N > 0 else "n"} 型\n{abs(N):.1e}'.replace('e+', '×10^').replace('×10^20', '×10²⁰').replace('×10^21', '×10²¹').replace('×10^19', '×10¹⁹'),
                    ha='center', va='center', fontsize=6.3, color='white' if dark else 'k', linespacing=1.1)
            crow.append([m, 'Sn' if m[:2] == 'Sn' else 'Pb', m[3:5], 'p' if N > 0 else 'n', f'{abs(N):.3e}', 'yes' if pf.HALL[m]['flagged'] else 'no'])
    ax.set_xlim(-0.75, 2.0); ax.set_ylim(-0.05, 2.0); ax.set_aspect('equal'); ax.axis('off')
    ax.text(0.48, 2.03, 'B = Sb\n易生成受主\n（如 Sb₂Te₃）', ha='center', va='bottom', fontsize=5.6, linespacing=1.05)
    ax.text(1.48, 2.03, 'B = Bi\n接近本征\n（如 Bi₂Te₃）', ha='center', va='bottom', fontsize=5.6, linespacing=1.05)
    ax.text(-0.08, 1.48, 'A = Sn\n易生成受主\n（如 SnTe）', ha='right', va='center', fontsize=5.6, linespacing=1.05)
    ax.text(-0.08, 0.48, 'A = Pb\n接近本征\n（如 PbTe）', ha='right', va='center', fontsize=5.6, linespacing=1.05)
    ax.set_title('(c) 为什么只有它浓度低（推测）', fontsize=7.6, pad=34, loc='center')
    ax.text(0.62, -0.12, '格中为室温霍尔浓度（cm⁻³），颜色越深越高。\nSn：Sn 空位；Sb：Sb 占 Te 位。均为母体化合物的\n已知规律，尚未对 Pb₂Bi₂Te₅ 做缺陷计算', ha='center', va='top', fontsize=5.4, color='#555555')
    dump('c_Hall_2x2.csv', ['material', 'A_site', 'B_site', 'type', 'Hall_conc_cm-3', 'Hall_flagged_by_instrument'], crow)
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(OUT, f'讲解_Pb2Bi2Te5的S为何最高.{ext}'))
    print('→', OUT)


if __name__ == '__main__':
    main()
