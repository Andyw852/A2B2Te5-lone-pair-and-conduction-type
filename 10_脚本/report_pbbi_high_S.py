#!/usr/bin/env python3
"""专题报告：Pb₂Bi₂Te₅ 样品的 Seebeck 系数为什么最高（2026-11-11）。

生成四幅图（每幅图的数据写成 CSV 与一个 Excel）、报告 HTML（正文中的数值由同一批数据填入），并打包：
  图 1  实验：四个样品的 S(T) 与功率因子 S²σ(T)
  图 2  计算：相同载流子浓度下四种材料的 |S|（PBE+SOC、HSE06+SOC，300 K）
  图 3  计算：Pb₂Bi₂Te₅ 的 |S| 随电子浓度；能带上两种浓度的化学势（PBE+SOC）
  图 4  实验：四个样品的室温霍尔浓度按 A 位 × B 位排成 2×2（浓度为何低：推测）
数值全部来自已有计算（bt2_pbe_dense、c8_hse_dense、c2_pbe）与实验原始数据，不做新的估计。
输出：21_阶段报告与交接/专题_Pb2Bi2Te5的S为何最高_20261111/（图、数据、Excel），
      21_阶段报告与交接/专题_Pb2Bi2Te5的S为何最高_20261111.html（单文件）与 .zip（报告 + 图 + 作图数据）
用法：python 10_脚本/report_pbbi_high_S.py
"""
import base64, csv, io, os, re, shutil, sys, tempfile, zipfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import transport_figures as tf
import plot_data_index as pdi
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from PIL import Image

TAG = '专题_Pb2Bi2Te5的S为何最高_20261111'
D = os.path.join(pf.ROOT, '21_阶段报告与交接'); OUT = os.path.join(D, TAG)
NAME = {'Pb2Sb2Te5': 'Pb₂Sb₂Te₅', 'Pb2Bi2Te5': 'Pb₂Bi₂Te₅', 'Sn2Sb2Te5': 'Sn₂Sb₂Te₅', 'Sn2Bi2Te5': 'Sn₂Bi₂Te₅'}
HN = {m: NAME[m].replace('₂', '<sub>2</sub>').replace('₅', '<sub>5</sub>') for m in NAME}
KT = 1000 * pf.mv.KB * 300
DATA = {}          # 图目录 → [(文件名, 表头, 行, 说明)]
V = {}             # 正文里的数值


def style():
    pf.style()
    for f in font_manager.findSystemFonts():
        if 'wqy' in f.lower() or 'wenquanyi' in f.lower(): font_manager.fontManager.addfont(f)
    have = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams['font.family'] = [n for n in ('Times New Roman', 'Liberation Serif') if n in have][:1] + ['WenQuanYi Zen Hei', 'DejaVu Serif']


def dump(folder, name, header, rows, note):
    DATA.setdefault(folder, []).append((name, header, rows, note))
    d = os.path.join(OUT, 'data', folder); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def save(fig, name):
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(OUT, 'figures', f'{name}.{ext}'))
    plt.close(fig)


def sup(N):
    """2.0e+20 → 2.0×10²⁰（HTML 用 <sup>）"""
    a, e = f'{abs(N):.1e}'.split('e'); return f'{a}×10<sup>{int(e)}</sup>'


# ------------------------------------------------------------------ 图 1：实验
def fig1():
    fig, A = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw=dict(wspace=0.30)); rows = []
    for m in pf.MATS:
        e = np.array(pf.read_exp(m)); T, sig, S = e[:, 0], e[:, 1], e[:, 2]; P = (S * 1e-6) ** 2 * sig * 1e4     # µW cm⁻¹ K⁻²
        h = pf.HALL[m]; c = pf.MCOL[m]; lab = f'{NAME[m]}（{"p" if h["N"] > 0 else "n"} {abs(h["N"]):.1e}）'.replace('e+', '×10^')
        lab = re.sub(r'\^(\d+)', lambda q: ''.join('⁰¹²³⁴⁵⁶⁷⁸⁹'[int(ch)] for ch in q.group(1)), lab)
        A[0].plot(T, S, '-o', color=c, ms=3.5, lw=1.1, mfc='white' if m == 'Pb2Sb2Te5' else c, label=lab)
        A[1].plot(T, P, '-o', color=c, ms=3.5, lw=1.1, mfc='white' if m == 'Pb2Sb2Te5' else c)
        rows += [[m, round(t, 2), round(s / 1e4, 4), round(x, 3), round(p, 3)] for t, s, x, p in zip(T, sig, S, P)]
        i = int(np.argmax(P)); V[f'{m}_PF300'] = P[0]; V[f'{m}_PFmax'] = P[i]; V[f'{m}_TPFmax'] = T[i]
        j = int(np.argmax(abs(S))); V[f'{m}_Smax'] = S[j]; V[f'{m}_TSmax'] = T[j]; V[f'{m}_S300'] = S[0]; V[f'{m}_sig300'] = sig[0] / 1e4
    A[0].axhline(0, color='k', lw=0.6); A[0].set_xlabel('T（K）'); A[0].set_ylabel('S（μV/K）'); A[0].set_ylim(-150, 80)
    A[0].set_title('(a) 实验 S(T)：只有 Pb₂Bi₂Te₅ 是 n 型，|S| 最大', fontsize=7.6, loc='left', pad=4)
    A[0].legend(fontsize=5.8, loc='upper left', bbox_to_anchor=(0.01, 0.63), title='样品（霍尔浓度，cm⁻³）', title_fontsize=5.8, labelspacing=0.2)
    A[1].set_xlabel('T（K）'); A[1].set_ylabel('功率因子 S²σ（μW cm⁻¹ K⁻²）'); A[1].set_ylim(0, 9.5)
    A[1].set_title('(b) 实验功率因子：Pb₂Bi₂Te₅ 最高', fontsize=7.6, loc='left', pad=4)
    for a in A: a.grid(color=pf.GRID, lw=0.4, zorder=0)
    A[1].text(0.03, 0.97, '空心：Pb₂Sb₂Te₅ 疑为两相样品', transform=A[1].transAxes, ha='left', va='top', fontsize=5.6, color='#555555')
    dump('Fig1_measured_S_PF', 'measured_T_sigma_S_PF.csv', ['material', 'T_K', 'sigma_1e4_S/m', 'S_uV/K', 'PF_uW/cm/K2 (= S^2 sigma)'], rows,
         '四个样品的实验 T、σ、S（实验原始表），功率因子 S²σ 由同一行的 S 与 σ 相乘得到。')
    dump('Fig1_measured_S_PF', 'Hall_room_temperature.csv', ['material', 'carrier_type', 'Hall_conc_cm-3', 'Hall_mobility_cm2/Vs', 'flagged_by_instrument'],
         [[m, 'p' if pf.HALL[m]['N'] > 0 else 'n', f"{abs(pf.HALL[m]['N']):.3e}", pf.HALL[m]['mu'], 'yes' if pf.HALL[m]['flagged'] else 'no'] for m in pf.MATS],
         '室温霍尔（实验原始表）。Pb₂Sb₂Te₅、Sn₂Sb₂Te₅ 的霍尔被仪器标红。')
    save(fig, 'Fig1_measured_S_PF')


# ------------------------------------------------------------------ 图 2：相同浓度
GROUPS = ((-8.11e19, 'n 型\n8.1×10¹⁹'), (8.11e19, 'p 型\n8.1×10¹⁹'), (9.07e20, 'p 型\n9.1×10²⁰'))


def fig2():
    fig, A = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw=dict(wspace=0.28)); rows = []; wd = 0.19
    for ax, fn in zip(A, ('PBE', 'HSE')):
        for g, (N, lab) in enumerate(GROUPS):
            for k, m in enumerate(pf.MATS):
                T, f = tf.LOAD[fn](m); v = abs(f(N, list(T).index(300))); V[f'{fn}_{m}_{g}'] = v
                ax.bar(g + (k - 1.5) * (wd + 0.01), v, wd, color=pf.MCOL[m], zorder=3, label=NAME[m] if g == 0 else None)
                rows.append([tf.FN[fn], m, 'n' if N < 0 else 'p', f'{abs(N):.3e}', round(v, 1)])
        ax.set_xticks(range(3), [g[1] for g in GROUPS], fontsize=6.6); ax.tick_params(axis='x', length=0, top=False)
        ax.set_ylim(0, 175); ax.set_ylabel('|S|，300 K（μV/K）'); ax.grid(axis='y', color=pf.GRID, lw=0.4, zorder=0)
        ax.set_xlabel('相同的载流子浓度（cm⁻³）')
        ax.set_title(f'({"ab"[fn == "HSE"]}) {tf.FN[fn]}：S 主要随浓度变，材料间差别小', fontsize=7.6, loc='left', pad=4)
        ax.text(0, 118 if fn == 'PBE' else 92, 'Pb₂Bi₂Te₅ 样品\n在这个浓度', ha='center', fontsize=5.8, color='#333333')
        ax.text(2, 47, '另外三个样品\n在 ~10²¹', ha='center', fontsize=5.8, color='#333333')
    A[0].legend(fontsize=5.6, loc='upper right', handlelength=1.0, labelspacing=0.2)
    dump('Fig2_same_concentration', 'S_at_same_concentration_300K.csv', ['functional', 'material', 'carrier_type', 'concentration_cm-3', 'abs_S_uV/K'], rows,
         '300 K、相同载流子浓度下的 |S|：BoltzTraP2（PBE+SOC：bt2_pbe_dense；HSE06+SOC：c8_hse_dense，mult 60），常数弛豫时间、不加剪刀，'
         '浓度方向在化学势网格上按 log|N| 线性取值。')
    save(fig, 'Fig2_same_concentration')


# ------------------------------------------------------------------ 图 3：浓度与化学势（Pb₂Bi₂Te₅）
def fig3():
    fig = plt.figure(figsize=(7.0, 2.8)); gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=0.30)
    ax = fig.add_subplot(gs[0]); NS = np.logspace(19, np.log10(3e21), 61); rows = []
    for fn, ls in (('PBE', '-'), ('HSE', (0, (3, 1.4)))):
        T, f = tf.LOAD[fn]('Pb2Bi2Te5'); it = list(T).index(300); v = np.array([-f(-N, it) for N in NS])
        ax.plot(NS, np.where(v > 0, v, np.nan), color=pf.C_PBE if fn == 'PBE' else pf.C_HSE, lw=1.3, ls=ls, label=f'计算 {tf.FN[fn]}')
        rows += [[tf.FN[fn], f'{N:.4e}', round(x, 2)] for N, x in zip(NS, v)]
        for N in (8.11e19, 3e20, 9.07e20, 1.72e21): V[f'{fn}_PbBi_S_{N:.2e}'] = f(-N, it)
    Nh = abs(pf.HALL['Pb2Bi2Te5']['N']); Se = abs(pf.EXP['Pb2Bi2Te5'][0][1])
    ax.plot(Nh, Se, 'o', ms=6, color='#333333', zorder=5, label='样品（实验，约 300 K）')
    T, f = tf.LOAD['PBE']('Pb2Bi2Te5'); it = list(T).index(300); Sh = abs(f(-9.07e20, it))
    ax.plot(9.07e20, Sh, 's', ms=5, mfc='white', mec=pf.C_PBE, mew=1.2, zorder=5, label='假如电子浓度为 9.1×10²⁰')
    ax.annotate('', (9.07e20, Sh + 4), (Nh * 1.15, Se - 4), arrowprops=dict(arrowstyle='->', lw=0.9, color='#555555'))
    Sl = abs(f(-8.11e19, it)); ax.text(2.6e20, 72, f'浓度 ×{9.07e20 / 8.11e19:.0f}\n|S| 降到 {Sh / Sl:.0%}', fontsize=6.3, color='#333333')
    ax.set_xscale('log'); ax.set_xlim(1e19, 3e21); ax.set_ylim(0, 200); ax.grid(color=pf.GRID, lw=0.4)
    ax.set_xlabel('电子浓度（cm⁻³）'); ax.set_ylabel('|S|，300 K（μV/K）')
    ax.set_title('(a) Pb₂Bi₂Te₅：电子越多，|S| 越小', fontsize=7.6, loc='left', pad=4); ax.legend(fontsize=5.8, loc='upper right', labelspacing=0.2)
    dump('Fig3_concentration_and_mu', 'Pb2Bi2Te5_S_vs_electron_conc_300K.csv', ['functional', 'electron_conc_cm-3', 'abs_S_uV/K'], rows,
         'Pb₂Bi₂Te₅ 的 n 型 |S|（300 K）随电子浓度，两种泛函的 BoltzTraP2。实验点：霍尔浓度 8.11×10¹⁹、约 302 K 的 S。')
    # (b) 能带与化学势
    ax = fig.add_subplot(gs[1]); B = tf.bands_fn('Pb2Bi2Te5', 'PBE'); x, E, nel = B['x'], B['E'], B['nel']; K = pf.KEND; brow = []
    for b in range(nel - 6, nel + 6): ax.plot(x[:K + 1], E[:K + 1, b], color='#B2182B' if b >= nel else '#1F3A93', lw=0.8, zorder=2)
    for N, col, tag in ((-8.11e19, '#E69F00', '样品 8.1×10¹⁹'), (-9.07e20, '#7A7A7A', '假如 9.1×10²⁰')):
        u = B['mu'](N); dep = 1000 * (u - B['gap']); S = f(N, it)
        ax.axhspan(u - 3 * KT / 1000, u + 3 * KT / 1000, color=col, alpha=0.18, lw=0, zorder=1); ax.axhline(u, color=col, lw=1.3, ls=pf.DASH, zorder=3)
        ax.text(0.03, u + 3 * KT / 1000 + 0.01, f'{tag}：S = {S:.0f} μV/K\nµ 在导带底上方 {dep:.0f} meV（{dep / KT:.0f} kBT）'.replace('-', '−'),
                transform=ax.get_yaxis_transform(), ha='left', va='bottom', fontsize=5.6, color='#222222', linespacing=1.1,
                bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.9), zorder=5)
        brow.append([f'{abs(N):.3e}', round(u, 4), round(dep, 1), round(dep / KT, 2), round(S, 1)]); V[f'mu_{abs(N):.2e}'] = (dep, dep / KT, S)
    ax.axhline(B['gap'], color='k', lw=0.5, ls=':'); ax.axhline(0, color='k', lw=0.5, ls=':')
    pf.kticks(ax, x); ax.set_ylim(-0.25, 0.80); ax.set_ylabel(r'$E - E_\mathrm{VBM}$ (eV)')
    ax.set_title('(b) 浓度低 → µ 离带边近 → S 大', fontsize=7.6, loc='left', pad=4)
    ax.text(0.02, 0.02, 'Pb₂Bi₂Te₅，PBE+SOC；阴影为 µ ± 3kBT', transform=ax.transAxes, fontsize=5.6, color='#555555', va='bottom',
            bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.85))
    dump('Fig3_concentration_and_mu', 'Pb2Bi2Te5_mu_300K_PBE.csv', ['electron_conc_cm-3', 'mu-EVBM_eV', 'mu_above_CBM_meV', 'mu_above_CBM_kBT', 'S_PBE+SOC_uV/K'], brow,
         '300 K 的化学势 µ（BoltzTraP2，对齐到能带路径的真实价带顶），以及 µ 高出导带底多少。')
    dump('Fig3_concentration_and_mu', 'Pb2Bi2Te5_PBE+SOC_bands.csv', ['k_index', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV'],
         [[i, pf.r(x[i], 5), pf.LAB.get(i, ''), pf.bandname(b, nel), pf.r(E[i, b], 5)] for b in range(nel - 6, nel + 6) for i in range(K + 1)],
         'Pb₂Bi₂Te₅ 的 PBE+SOC 路径能带（价带顶以下 6 条、以上 6 条），能量相对价带顶。')
    save(fig, 'Fig3_concentration_and_mu')


# ------------------------------------------------------------------ 图 4：2×2 霍尔浓度
def fig4():
    fig, ax = plt.subplots(figsize=(3.6, 3.0)); grid = [['Sn2Sb2Te5', 'Sn2Bi2Te5'], ['Pb2Sb2Te5', 'Pb2Bi2Te5']]
    cmap = LinearSegmentedColormap.from_list('n', ['#F2F0F7', '#BCBDDC', '#756BB1', '#3F007D']); norm = LogNorm(5e19, 2e21); rows = []
    for i, r in enumerate(grid):
        for j, m in enumerate(r):
            N = pf.HALL[m]['N']; dark = abs(N) > 4e20; tc = 'white' if dark else 'k'
            ax.add_patch(plt.Rectangle((j, 1 - i), 0.96, 0.96, color=cmap(norm(abs(N))), lw=0))
            ax.text(j + 0.48, 1 - i + 0.70, NAME[m], ha='center', va='center', fontsize=7.5, color=tc, fontweight='bold')
            n_ = f'{abs(N):.1e}'.split('e'); nn = n_[0] + '×10' + ''.join('⁰¹²³⁴⁵⁶⁷⁸⁹'[int(c)] for c in str(int(n_[1])))
            ax.text(j + 0.48, 1 - i + 0.36, f'{"p" if N > 0 else "n"} 型 {nn}\nS = {pf.EXP[m][0][1]:+.0f} μV/K'.replace('-', '−'),
                    ha='center', va='center', fontsize=6.4, color=tc, linespacing=1.2)
            rows.append([m, m[:2], m[3:5], 'p' if N > 0 else 'n', f'{abs(N):.3e}', round(pf.EXP[m][0][1], 2), 'yes' if pf.HALL[m]['flagged'] else 'no'])
    ax.set_xlim(-0.8, 2.0); ax.set_ylim(-0.05, 2.0); ax.set_aspect('equal'); ax.axis('off')
    for j, t in enumerate(('B = Sb\n易生成受主\n（如 Sb₂Te₃）', 'B = Bi\n接近本征\n（如 Bi₂Te₃）')): ax.text(j + 0.48, 2.03, t, ha='center', va='bottom', fontsize=6, linespacing=1.05)
    for i, t in enumerate(('A = Sn\n易生成受主\n（如 SnTe）', 'A = Pb\n接近本征\n（如 PbTe）')): ax.text(-0.06, 1.48 - i, t, ha='right', va='center', fontsize=6, linespacing=1.05)
    ax.set_title('为什么只有 Pb₂Bi₂Te₅ 浓度低（推测）', fontsize=8, pad=36)
    ax.text(0.6, -0.1, '格中为室温霍尔浓度（cm⁻³）与约 300 K 的 S，颜色越深浓度越高。\nSn：Sn 空位；Sb：Sb 占 Te 位。均为母体化合物的已知规律，\n尚未对 Pb₂Bi₂Te₅ 做缺陷计算。',
            ha='center', va='top', fontsize=5.6, color='#555555')
    dump('Fig4_Hall_2x2', 'Hall_and_S_by_site.csv', ['material', 'A_site', 'B_site', 'type', 'Hall_conc_cm-3', 'S_exp_300K_uV/K', 'Hall_flagged_by_instrument'], rows,
         '四个样品的室温霍尔浓度与约 300 K 的 S（实验），按 A 位（Sn、Pb）× B 位（Sb、Bi）排列。')
    save(fig, 'Fig4_Hall_2x2')


# ------------------------------------------------------------------ 报告
def table(head, rows, num_from=1):
    h = ''.join(f'<th{" class=num" if j >= num_from else ""}>{x}</th>' for j, x in enumerate(head))
    b = ''.join('<tr>' + ''.join(f'<td{" class=num" if j >= num_from else ""}>{x}</td>' for j, x in enumerate(r)) + '</tr>' for r in rows)
    return f'<table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>'


def html():
    css = open(os.path.join(D, '论文报告_20261028.src.html'), encoding='utf-8').read()
    css = css[css.index('<style>'):css.index('</style>') + 8]
    t8 = {r['material']: r for r in csv.DictReader(open(os.path.join(D, '论文图_20261028', 'source_data', 'Table8_Pisarenko_mass',
                                                                     'samples_Pisarenko_mass_HSE_and_PBE.csv'), encoding='utf-8-sig'))}
    m = 'Pb2Bi2Te5'; f1 = lambda v: f'{v:.1f}'.replace('-', '−'); f0 = lambda v: f'{v:.0f}'.replace('-', '−')
    T1 = table(['样品', '霍尔浓度（cm⁻³）', '霍尔迁移率（cm²/V·s）', 'S，约 300 K（µV/K）', 'σ，约 300 K（10⁴ S/m）', 'S²σ，约 300 K（µW cm⁻¹ K⁻²）', 'S²σ 最大（温度）'],
               [[HN[k], ('p ' if pf.HALL[k]['N'] > 0 else 'n ') + sup(pf.HALL[k]['N']) + ('（标红）' if pf.HALL[k]['flagged'] else ''), f"{pf.HALL[k]['mu']:.1f}",
                 f1(V[f'{k}_S300']), f1(V[f'{k}_sig300']), f'{V[f"{k}_PF300"]:.2f}', f'{V[f"{k}_PFmax"]:.1f}（{V[f"{k}_TPFmax"]:.0f} K）'] for k in pf.MATS])
    T2 = table(['相同浓度（cm⁻³）', '泛函'] + [HN[k] for k in pf.MATS],
               [[g[1].replace('\n', ' ').replace('10¹⁹', '10<sup>19</sup>').replace('10²⁰', '10<sup>20</sup>'), tf.FN[fn]] + [f0(V[f'{fn}_{k}_{gi}']) for k in pf.MATS]
                for gi, g in enumerate(GROUPS) for fn in ('PBE', 'HSE')], num_from=2)
    mu1, mu2 = V['mu_8.11e+19'], V['mu_9.07e+20']
    rr = sorted(V[f'{fn}_{k}_1'] / V[f'{fn}_{k}_2'] for fn in ('PBE', 'HSE') for k in pf.MATS)
    S = lambda N, fn='PBE': f1(V[f'{fn}_PbBi_S_{N:.2e}'])
    rng = lambda fn, g: (f0(min(V[f'{fn}_{k}_{g}'] for k in pf.MATS)), f0(max(V[f'{fn}_{k}_{g}'] for k in pf.MATS)))
    fig = lambda n, cap, w='': f'<figure><img src="{{{{FIG:{n}}}}}" alt="{n}"{w}><figcaption>{cap}</figcaption></figure>'
    body = f'''<h1>Pb<sub>2</sub>Bi<sub>2</sub>Te<sub>5</sub> 样品的 Seebeck 系数为什么最高</h1>
<div class="meta">A₂B₂Te₅（A = Pb, Sn；B = Sb, Bi）专题说明 · 2026-11-11 · 实验数据与理论计算（PBE+SOC、HSE06+SOC + BoltzTraP2）对照 · 每幅图的数据见压缩包的 <code>data/</code> 与 <code>作图数据/</code></div>

<div class="box"><b>结论</b>
<ol>
<li><b>Pb₂Bi₂Te₅ 样品的 |S| 和功率因子在四个样品中最高</b>（约 300 K：S = {f1(V[f"{m}_S300"])} µV/K，S²σ = {V[f"{m}_PF300"]:.1f} µW cm⁻¹ K⁻²；图 1）。<b>但在相同的载流子浓度下，它的 S 与另外三种材料相近，并不更高</b>（图 2）。</li>
<li><b>它的 S 高，是因为样品的载流子浓度低一个数量级</b>：电子 8.1×10<sup>19</sup> cm⁻³，另外三个样品是约 10<sup>21</sup> cm⁻³ 的空穴。浓度低，化学势只比导带底高 {mu1[0]:.0f} meV（{mu1[1]:.0f} k<sub>B</sub>T），S 就大；同一材料若有 9.1×10<sup>20</sup> 个电子，S 只有 {f0(mu2[2])} µV/K（图 3）。</li>
<li><b>浓度为什么低（推测，尚未计算）</b>：另外三个样品都含有容易产生受主的 Sn（Sn 空位）或 Sb（Sb 占 Te 位），都是重掺杂的 p 型；只有 Pb₂Bi₂Te₅ 两者都不含，受主少，较弱的施主占上风，成为浓度不高的 n 型（图 4）。</li>
<li><b>启示</b>：n 型并不本质上更好，Pb₂Bi₂Te₅ 只是退火后浓度正好落在 10<sup>20</sup> 附近，离最佳区最近。计算表明，同一浓度下 p 型的 S 甚至更大（8.1×10<sup>19</sup> 时 {rng("PBE", 1)[0]}–{rng("PBE", 1)[1]} µV/K）；若能把 Sn 系的空穴降到 ~10<sup>20</sup>，它们的 S 也会大幅升高。</li>
</ol></div>

<h2>1 实验事实：Pb₂Bi₂Te₅ 的 |S| 与功率因子最高</h2>
{fig("Fig1_measured_S_PF", "<b>图 1</b> 四个退火样品的实验 S(T) 与功率因子 S²σ(T)。图例中为室温霍尔浓度。Pb₂Sb₂Te₅ 的 S 异常小，疑为 n、p 两相共存（主报告 5.3 节），用空心点。")}
<p class="tcap"><b>表 1</b> 四个样品的室温霍尔与约 300 K 的输运（实验）。“标红”：测试仪器把这次霍尔结果标红（可信度较低）。</p>
{T1}
<ul>
<li>只有 Pb₂Bi₂Te₅ 是 n 型，浓度 8.1×10<sup>19</sup> cm⁻³，比另外三个样品低一个数量级；它的霍尔迁移率也最高（{pf.HALL[m]["mu"]:.0f} cm²/V·s，另外三个 {min(pf.HALL[k]["mu"] for k in pf.MATS if k != m):.1f}–{max(pf.HALL[k]["mu"] for k in pf.MATS if k != m):.1f}）。</li>
<li>它的电导率只有另外三个的 {min(V[f"{m}_sig300"] / V[f"{k}_sig300"] for k in pf.MATS if k != m):.0%}–{max(V[f"{m}_sig300"] / V[f"{k}_sig300"] for k in pf.MATS if k != m):.0%}，但 |S| 是它们的 {min(abs(V[f"{m}_S300"] / V[f"{k}_S300"]) for k in pf.MATS if k != m):.0f}–{max(abs(V[f"{m}_S300"] / V[f"{k}_S300"]) for k in pf.MATS if k != m):.0f} 倍，S²σ 最高；|S| 在 {V[f"{m}_TSmax"]:.0f} K 达到最大 {f1(abs(V[f"{m}_Smax"]))} µV/K，S²σ 最大 {V[f"{m}_PFmax"]:.1f} µW cm⁻¹ K⁻²。</li>
<li>ZT 还需要热导率，本工作没有热导率数据。σ 低时电子热导也低（Wiedemann–Franz 关系），对 ZT 有利，但这里不作定量判断。</li>
</ul>

<h2>2 同一浓度下，Pb₂Bi₂Te₅ 的 S 并不更高</h2>
<p>用同一套 BoltzTraP2 计算（常数弛豫时间，不加剪刀），在三个相同的浓度下比较四种材料的 300 K |S|：Pb₂Bi₂Te₅ 样品的电子浓度（n 8.1×10<sup>19</sup>）、同样的空穴浓度（p 8.1×10<sup>19</sup>），以及另外三个样品附近的空穴浓度（p 9.1×10<sup>20</sup>）。</p>
{fig("Fig2_same_concentration", "<b>图 2</b> 相同载流子浓度下四种材料的 |S|（300 K）。(a) PBE+SOC，(b) HSE06+SOC。p 型浓度从 8.1×10<sup>19</sup> 升到 9.1×10<sup>20</sup>，每种材料的 |S| 都降为原来的 1/" + f'{rr[-1]:.1f}' + "–1/" + f'{rr[0]:.1f}' + "；同一浓度下材料之间的差别小得多：S 的差别主要来自浓度，而不是材料。")}
<p class="tcap"><b>表 2</b> 相同浓度下的 |S|（µV/K，300 K）。</p>
{T2}
<ul>
<li>电子都是 8.1×10<sup>19</sup> 时，PBE 下四种材料为 {rng("PBE", 0)[0]}–{rng("PBE", 0)[1]} µV/K，Pb₂Bi₂Te₅（{f0(V["PBE_Pb2Bi2Te5_0"])}）是其中最小的；HSE 下它同样最小（{f0(V["HSE_Pb2Bi2Te5_0"])}，HSE 的 n 型 Pb₂Bi₂Te₅ 偏低，见下一条）。</li>
<li>空穴都是 8.1×10<sup>19</sup> 时为 {rng("PBE", 1)[0]}–{rng("PBE", 1)[1]}（PBE）、{rng("HSE", 1)[0]}–{rng("HSE", 1)[1]} µV/K（HSE）；都是 9.1×10<sup>20</sup> 时只有 {rng("PBE", 2)[0]}–{rng("PBE", 2)[1]}、{rng("HSE", 2)[0]}–{rng("HSE", 2)[1]} µV/K。</li>
<li>n 型 Pb₂Bi₂Te₅ 以 PBE 为准：在霍尔浓度下 PBE 给出 {f1(-V["PBE_Pb2Bi2Te5_0"])} µV/K，与实验 {f1(V[f"{m}_S300"])} 相符；HSE 的导带只剩 Γ 单谷，只有 {f1(-V["HSE_Pb2Bi2Te5_0"])}（主报告 5.1 节）。</li>
</ul>

<h2>3 浓度低，化学势离带边近，S 就大</h2>
{fig("Fig3_concentration_and_mu", "<b>图 3</b> (a) Pb₂Bi₂Te₅ 的 n 型 |S| 随电子浓度（300 K），圆点为样品；方块为同一材料在 9.1×10<sup>20</sup>（Sn₂Bi₂Te₅ 样品的空穴浓度）下的计算值。(b) Pb₂Bi₂Te₅ 的 PBE+SOC 能带（价带蓝、导带红）及两种电子浓度下 300 K 的化学势 µ（虚线），阴影为 µ ± 3k<sub>B</sub>T，即参与导电的能量窗口。")}
<ul>
<li><b>S 衡量每个载流子比化学势多带了多少能量。</b>参与导电的只是 µ 附近几个 k<sub>B</sub>T 内的态（−∂f/∂E 的权重）。
  <ul>
  <li>µ 离导带底近（非简并、弱简并）：窗口下方是带隙，没有态，导电的电子几乎都在 µ 之上，不能抵消，S 大。</li>
  <li>µ 深入能带（强简并）：µ 上下都有很多态，高于和低于 µ 的电子对热流的贡献符号相反、大部分抵消，S 小。</li>
  </ul></li>
<li><b>电子越多，就要把导带填得越深。</b>Pb₂Bi₂Te₅ 样品（8.1×10<sup>19</sup>）的 µ 只比导带底高 {mu1[0]:.0f} meV（{mu1[1]:.0f} k<sub>B</sub>T），S = {f0(mu1[2])} µV/K；电子增加到 9.1×10<sup>20</sup> 时，µ 高出导带底 {mu2[0]:.0f} meV（{mu2[1]:.0f} k<sub>B</sub>T），S 只有 {f0(mu2[2])}；到 1.7×10<sup>21</sup> 只有 {S(1.72e21)} µV/K（PBE）。</li>
<li><b>能带的作用：导带是多谷的。</b>在 8.1×10<sup>19</sup> 下 S 能达到 {f0(mu1[2])} µV/K，靠的是导带态密度大：PBE 下导带底 23 meV 内有三个谷（含以 Bi 6p 为主的 A 谷），实验读出的态密度有效质量为 {float(t8[m]["m_DOS_exp_SPB_me"]):.2f} m<sub>e</sub>，HSE 的 Γ 单谷只有 {float(t8[m]["single_valley_m_DOS_HSE_me (geometric mean)"]):.2f} m<sub>e</sub>。但这一点另外三种材料的 n 型也差不多（图 2），所以它解释了 S 的大小，不解释 Pb₂Bi₂Te₅ 为什么最高。</li>
</ul>

<h2>4 浓度为什么低（推测）</h2>
{fig("Fig4_Hall_2x2", "<b>图 4</b> 四个样品的室温霍尔浓度与约 300 K 的 S，按 A 位 × B 位排列。行、列标签为母体二元化合物的已知规律。", ' style="max-width:min(100%,430px);display:block;margin:0 auto"')}
<p>样品是退火后直接测量的，载流子浓度由退火时形成的本征缺陷与成分偏差决定。<b>本工作只对 Pb₂Sb₂Te₅ 做过缺陷计算</b>（退火后的有序层状相中，受主 Pb<sub>Sb</sub> 的中性形成能在各类缺陷中最低：富 Te 0.14 eV、贫 Te 0.25 eV，呈 p 型；主报告第 6 节表 9），没有算 Pb₂Bi₂Te₅。下面的解释是推测：</p>
<ul>
<li><b>母体化合物的已知规律</b>：SnTe 中 Sn 空位极易形成，本征就是重 p 型；Sb₂Te₃ 中 Sb 占 Te 位的反位缺陷极易形成，也是强 p 型；PbTe、Bi₂Te₃ 接近本征（Bi₂Te₃ 偏富 Te 时为 n 型）。</li>
<li><b>四个样品正好对上</b>：含 Sn 或含 Sb 的三个都是约 10<sup>21</sup> 的 p 型，同时含两者的 Sn₂Sb₂Te₅ 浓度最高；只有两者都不含的 Pb₂Bi₂Te₅ 是低浓度 n 型。</li>
<li><b>可能的施主</b>：Bi 占 Pb 位（Bi<sub>Pb</sub>）与退火时 Te 挥发留下的 Te 空位。Bi³⁺ 与 Pb²⁺ 大小接近，施主 Bi<sub>Pb</sub> 与受主 Pb<sub>Bi</sub> 可能都不少、互相补偿，净载流子少，符号由两者的微小差别决定。</li>
</ul>

<h2>5 建议</h2>
<ol>
<li><b>计算</b>：按 Pb₂Sb₂Te₅ 的同一流程，算 Pb₂Bi₂Te₅ 层状相的本征缺陷形成能与平衡载流子浓度；同时算 Sn₂Bi₂Te₅ 作对照，检验 Sn 相关的受主是不是使它成为 p 型的原因。</li>
<li><b>实验</b>：EPMA/EDS 测成分，看是否缺 Te；比较不同退火时间、温度下的霍尔浓度；测热导率以得到 ZT。</li>
<li><b>掺杂</b>：Sn 系若能把空穴降到 8.1×10<sup>19</sup> cm⁻³，计算的 S 为 {f0(min(V[f"{fn}_{k}_1"] for fn in ("PBE", "HSE") for k in ("Sn2Sb2Te5", "Sn2Bi2Te5")))}–{f0(max(V[f"{fn}_{k}_1"] for fn in ("PBE", "HSE") for k in ("Sn2Sb2Te5", "Sn2Bi2Te5")))} µV/K（图 2、表 2，两种泛函），可检验“性能主要由浓度决定”的结论。</li>
</ol>

<h2>附：数据来源与方法</h2>
<ul>
<li><b>实验</b>：四个退火样品的 T、σ、S 原始表与室温霍尔（<code>21_阶段报告与交接/实验数据/</code>）。功率因子 S²σ 由同一温度下的 S 与 σ 相乘。Pb₂Sb₂Te₅、Sn₂Sb₂Te₅ 的霍尔被仪器标红。</li>
<li><b>计算</b>：所给 9 原子层状结构；BoltzTraP2，常数弛豫时间，不加剪刀，插值倍数 60，面内 S 取 (xx + yy)/2。两种泛函用同一网格（Sb 系 28×28×7、Bi 系 28×28×6）：PBE+SOC 为自洽后在该网格上的本征值（结果在 <code>22_能带基准_1-0/bt2_pbe_dense</code>），HSE06+SOC 在 HSE 自洽之上用 KPOINTS_OPT 算该网格的本征值（<code>22_能带基准_1-0/&lt;材料&gt;/c8_hse_dense</code>）。给定浓度时在化学势网格上按 log|N| 线性取值，温度不插值。图 3(b) 的化学势用 Γ 处的本征值对齐到能带路径的真实价带顶。</li>
<li><b>不做模型或估计</b>：文中所有计算值都直接取自 BoltzTraP2 结果；有效质量取自主报告表 8。</li>
<li>生成程序：<code>10_脚本/report_pbbi_high_S.py</code>（图、数据、本报告与压缩包）。</li>
</ul>
'''
    return f'<!doctype html>\n<html lang="zh-CN">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n<title>Pb₂Bi₂Te₅ 的 S 为何最高</title>\n{css}\n</head>\n<body>\n<main>\n{body}\n</main>\n</body>\n</html>\n'


def embed(mo):
    im = Image.open(os.path.join(OUT, 'figures', mo.group(1) + '.png')).convert('RGB')
    if im.width > 1800: im = im.resize((1800, round(im.height * 1800 / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, 'PNG', optimize=True); return 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode()


def excel(folder, path):
    import xlsxwriter
    wb = xlsxwriter.Workbook(path, {'strings_to_formulas': False, 'strings_to_urls': False, 'nan_inf_to_errors': False})
    b = wb.add_format({'bold': True}); wr = wb.add_format({'text_wrap': True, 'valign': 'top'})
    ix = wb.add_worksheet('说明'); ix.write_string(0, 0, f'{folder}：每个数据表一页', b)
    for r, (name, _, _, note) in enumerate(DATA[folder]): ix.write_string(r + 2, 0, name[:-4][:31]); ix.write_string(r + 2, 1, note, wr)
    ix.set_column(0, 0, 36); ix.set_column(1, 1, 110)
    for name, header, rows, note in DATA[folder]:
        ws = wb.add_worksheet(name[:-4][:31])
        for j, h in enumerate(header): ws.write_string(0, j, h, b)
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                x = pdi.number(str(v)); (ws.write_number if isinstance(x, float) else ws.write_string)(i + 1, j, x)
        ws.freeze_panes(1, 0); ws.set_column(0, len(header), 18)
    wb.close(); pdi.check_book(path)


def main():
    style(); shutil.rmtree(OUT, ignore_errors=True); os.makedirs(os.path.join(OUT, 'figures'))
    fig1(); fig2(); fig3(); fig4()
    for folder in DATA: excel(folder, os.path.join(OUT, '作图数据_' + folder + '.xlsx'))
    src = html(); open(os.path.join(OUT, 'report.src.html'), 'w', encoding='utf-8').write(src)
    single = re.sub(r'\{\{FIG:([^}]+)\}\}', embed, src)
    open(os.path.join(D, TAG + '.html'), 'w', encoding='utf-8').write(single)
    # 压缩包
    tmp = tempfile.mkdtemp(); top = os.path.join(tmp, TAG)
    os.makedirs(os.path.join(top, 'figures', 'png')); os.makedirs(os.path.join(top, 'figures', 'pdf')); os.makedirs(os.path.join(top, '作图数据'))
    for f in os.listdir(os.path.join(OUT, 'figures')): shutil.copy(os.path.join(OUT, 'figures', f), os.path.join(top, 'figures', f.rsplit('.', 1)[1]))
    shutil.copytree(os.path.join(OUT, 'data'), os.path.join(top, 'data'))
    for folder in DATA: shutil.copy(os.path.join(OUT, '作图数据_' + folder + '.xlsx'), os.path.join(top, '作图数据', folder + '.xlsx'))
    open(os.path.join(top, '报告.html'), 'w', encoding='utf-8').write(re.sub(r'\{\{FIG:([^}]+)\}\}', lambda mo: f'figures/png/{mo.group(1)}.png', src))
    open(os.path.join(top, '报告_单文件.html'), 'w', encoding='utf-8').write(single)
    readme = ['Pb₂Bi₂Te₅ 样品的 Seebeck 系数为什么最高：专题说明与作图数据（2026-11-11）', '',
              '报告.html            报告（图片从 figures/png 读取，请整体解压后打开）',
              '报告_单文件.html     同一报告，图片已嵌入，可单独打开或转发',
              'figures/png、pdf     图 1–4（600 dpi PNG 与矢量 PDF）',
              'data/<图>/           每幅图的源数据（CSV，UTF-8 带 BOM）',
              '作图数据/<图>.xlsx   同样的数据，每幅图一个 Excel，第一页“说明”写明每个表的内容', '']
    for folder, items in DATA.items():
        readme.append(folder + '/')
        readme += [f'  {name}  {note}' for name, _, _, note in items]
    readme += ['', '数值全部来自已有计算（BoltzTraP2：PBE+SOC bt2_pbe_dense、HSE06+SOC c8_hse_dense）与实验原始数据，不做模型或估计。',
               '生成程序（仓库中）：10_脚本/report_pbbi_high_S.py']
    open(os.path.join(top, 'README.txt'), 'w', encoding='utf-8-sig', newline='').write('\r\n'.join(readme) + '\r\n')
    out = os.path.join(D, TAG + '.zip')
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for dp, _, fs in sorted(os.walk(top)):
            for f in sorted(fs): z.write(os.path.join(dp, f), os.path.relpath(os.path.join(dp, f), tmp))
    n = sum(len(fs) for _, _, fs in os.walk(top)); shutil.rmtree(tmp)
    print(out, f'{os.path.getsize(out) / 1e6:.1f} MB，{n} 个文件')


if __name__ == '__main__':
    main()
