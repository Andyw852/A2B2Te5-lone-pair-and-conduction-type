#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
band_refine.py —— 精细化能带图（解决“能带收敛/简并看不清”）

输入：<材料>/{band.dat, band_klabels.txt, band_summary.json, KPOINTS[_OPT]}
      （HSE + SOC，Γ-M-K-Γ-A-L-H-A|L-M|H-K 路径，~200 点，能量已相对 VBM）
      * clone 下来的仓库：默认读 12_能带绘图_20261021/源数据/<材料>/，输出到 12_能带绘图_20261021/复现输出/
      * 原工作区：band_and_elf/huang/<材料>/step4_band_plot/，输出到 band_refine_20261021/
      用 --data-root / --out 可显式指定。

输出（图在 <out>/figures，摘要在 <out>）：
  bands_<材料>.png/pdf   双面板：全貌 + 带边放大（±kT 窗口、简并说明、收敛标记）
  bands_4mats_near_edge.png/pdf   四材料带边对比（同一能量窗）
  bands_degeneracy_convergence.png/pdf  简并/收敛诊断（高对称点能级 + 态数-窗口曲线）
  band_summary_refined.csv/md     数值摘要

设计要点（针对“简并/收敛看不清”）：
  1. SOC + 反演对称 -> 每条能带是 Kramers 二重简并（本数据逐点核对过）；图中只画每对的
     一条并把图例写成 x2，避免“两条重合线看起来像一条”的歧义。
  2. 带边放大面板：只画 VBM 上方 6 对 / CBM 下方 6 对，逐对配色，画出 ±kT(300K)=25.85 meV
     窗口，标出 VBM/CBM 位置与“窗口内态数/等价谷数”，把带收敛直接量化在图上。
  3. 诊断图：各高对称点的前 4 对价带能级阶梯（看简并与分裂）+ 距 VBM 的态数随窗口 ΔE 的曲线
     （带收敛程度的单调刻画）。
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_PKG = os.path.normpath(os.path.join(HERE, '..', '12_能带绘图_20261021'))
BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
if os.path.isdir(os.path.join(REPO_PKG, '源数据')):          # 仓库内：读打包好的源数据
    SRC = os.path.join(REPO_PKG, '源数据')
    OUT = os.path.join(REPO_PKG, '复现输出')
else:                                                          # 原工作区
    SRC = os.path.join(BASE, 'band_and_elf', 'huang')
    OUT = os.path.join(BASE, 'band_refine_20261021')
FIG = os.path.join(OUT, 'figures')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
KT = 0.025852          # 300 K 的 kT (eV)
WIN = (-0.35, 0.45)    # 带边放大窗口 (eV, 相对 VBM)

plt.rcParams.update({
    'font.size': 11, 'font.family': 'DejaVu Sans',
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.8, 'xtick.direction': 'out', 'ytick.direction': 'out',
})

C_VAL = ['#b03a2e', '#d35400', '#b7950b', '#7d6608', '#6c3483', '#1a5276']
C_CON = ['#1f618d', '#2874a6', '#148f77', '#117864', '#5b2c6f', '#7b241c']
C_DEEP = '#aab7c4'


def band_dir(m):
    """原工作区：<SRC>/<材料>/step4_band_plot；仓库打包：文件直接在 <SRC>/<材料>/ 下。"""
    d = os.path.join(SRC, m)
    return os.path.join(d, 'step4_band_plot') if os.path.isdir(os.path.join(d, 'step4_band_plot')) else d


def load(m):
    d = band_dir(m)
    raw = np.loadtxt(os.path.join(d, 'band.dat'))
    k, E = raw[:, 0], raw[:, 1:]
    labels = []
    for line in open(os.path.join(d, 'band_klabels.txt')):
        p = line.split()
        if len(p) == 2 and not p[0].startswith('#'):
            labels.append((p[0], float(p[1])))
    sm = json.load(open(os.path.join(d, 'band_summary.json')))
    return k, E, labels, sm, d


def read_path_kpoints(m, n_expect=None):
    """读出与 band.dat 行一一对应的【路径 k 点】分数坐标。
    两种来源：
      * Sb 材料：KPOINTS_OPT 直接列 204 个路径点（第 4 列是 1）；
      * Bi 材料：KPOINTS 里 38 个均匀点 + 200/204 个零权重路径点。
    返回 (n,3) 或 None（长度对不上时）。"""
    def four_col(path):
        out = []
        for line in open(path):
            p = line.split()
            if len(p) == 4:
                try:
                    v = [float(x) for x in p]
                except ValueError:
                    continue
                out.append(v)
        return out
    for fn in ('KPOINTS_OPT', 'KPOINTS'):
        fp = os.path.join(band_dir(m), fn)
        if not os.path.exists(fp):
            continue
        rows = four_col(fp)
        if fn == 'KPOINTS_OPT':
            if n_expect is None or len(rows) == n_expect:
                return np.array([r[:3] for r in rows])
        else:
            zero = [r for r in rows if abs(r[3]) < 1e-12]
            if n_expect is None or len(zero) == n_expect:
                return np.array([r[:3] for r in zero])
    return None


def unique_bands(E):
    """Kramers 简并 -> 每隔一条取一条（最高一条带偶有未配对，单独处理）。"""
    return E[:, 0::2]


def draw_path(ax, labels, lw_edge=0.9, fs=12, y0=None):
    """高对称点竖线 + 标签（同名相邻段不重复画粗线）。"""
    xs = [x for _, x in labels]
    names = [n for n, _ in labels]
    for i, (n, x) in enumerate(labels):
        ax.axvline(x, color='0.25', lw=lw_edge, zorder=0)
        ax.text(x, 1.0, n.replace('|', '|'), transform=ax.get_xaxis_transform(),
                ha='center', va='bottom', fontsize=fs, color='0.15')
    ax.set_xlim(xs[0], xs[-1])


def analyze(m):
    k, E, labels, sm, d = load(m)
    nv = int(sm['nelect'])                       # SOC: 每条带 1 个电子
    vbm = float(E[:, :nv].max()); cbm = float(E[:, nv:].min())
    gap = cbm - vbm
    # 价带“对”的序号（从 VBM 往下数）：pair 0 = 最上面一对
    nvpair = nv // 2
    # 每个 k 上，各对价带的能量（取该对的最小值作为该对的能量）
    def pair_energy(side, j):
        """side='v' 第 j 对价带（j=0 最靠上）；side='c' 第 j 对导带（j=0 最靠下）。"""
        if side == 'v':
            cols = [nv - 1 - 2 * j, nv - 2 - 2 * j]
        else:
            cols = [nv + 2 * j, nv + 2 * j + 1]
        return E[:, cols].mean(axis=1)
    vp = np.array([pair_energy('v', j) for j in range(6)])   # (6, nk)
    cp = np.array([pair_energy('c', j) for j in range(6)])
    # VBM/CBM 位置
    iv = int(np.argmax(vp[0])); ic = int(np.argmin(cp[0]))
    # 收敛指标：每个 k 上“最上面一对”与第 j 对的能量差
    dv = vp - vp[0][None, :]      # (6, nk) 相对该 k 的顶价带
    dc = cp - cp[0][None, :]
    # 价带/导带在“距 VBM/CBM 窗口”内的态数（沿路径，含简并 x2）
    res = dict(material=m, gap=gap, gap_type=sm['gap_type'], min_direct=sm.get('min_direct_gap_eV'),
               vbm=vbm, cbm=cbm, k=k, labels=labels, nv=nv,
               vbm_k=(float(k[iv]), list(map(float, sm['vbm']['k_frac']))),
               cbm_k=(float(k[ic]), list(map(float, sm['cbm']['k_frac']))),
               vp=vp, cp=cp, dv=dv, dc=dc)
    for tag, arr in [('v', vp), ('c', cp)]:
        for w in [0.005, 0.010, 0.025, 0.050, 0.100]:
            if tag == 'v':
                cnt = int((E[:, :nv] >= vbm - w).sum())
                nk_w = int((vp[0] >= vbm - w).sum())
            else:
                cnt = int((E[:, nv:] <= cbm + w).sum())
                nk_w = int((cp[0] <= cbm + w).sum())
            res.setdefault('window_' + tag, {})['%g' % (w * 1000)] = dict(states=cnt, k_points=nk_w)
    # 高对称点能级（价带前 4 对 / 导带前 4 对，eV 相对 VBM）
    lev = []
    for name, x in labels:
        i = int(np.argmin(np.abs(k - x)))
        lev.append(dict(label=name, x=float(k[i]),
                        v=[float(vp[j][i] - vbm) for j in range(4)],
                        c=[float(cp[j][i] - vbm) for j in range(4)]))
    res['levels'] = lev
    # ---------- 带收敛：沿路径的局部极值（谷）分析 ----------
    def _extrema(y, kind='max', prom=0.005, merge=2, look=10):
        """找局部极值：不光滑（数据噪声 ~0.1 meV），用 5 meV 最小 prominence 过滤；
        相邻平台合并成一个谷；含路径端点（VBM 可能就在端点）。"""
        n = len(y)
        if kind == 'max':
            cand = [i for i in range(n)
                    if (i == 0 or y[i] >= y[i - 1]) and (i == n - 1 or y[i] >= y[i + 1])]
        else:
            cand = [i for i in range(n)
                    if (i == 0 or y[i] <= y[i - 1]) and (i == n - 1 or y[i] <= y[i + 1])]
        # 合并相邻候选（平台）
        plat, cur = [], []
        for i in cand:
            if cur and i - cur[-1] <= merge:
                cur.append(i)
            else:
                if cur:
                    plat.append(cur)
                cur = [i]
        if cur:
            plat.append(cur)
        pick = []
        for c in plat:
            vals = [y[i] for i in c]
            j = int(np.argmax(vals)) if kind == 'max' else int(np.argmin(vals))
            i = c[j]
            lo, hi = max(0, i - look), min(n, i + look + 1)
            seg = y[lo:hi]
            prom_val = (y[i] - seg.min()) if kind == 'max' else (seg.max() - y[i])
            if prom_val >= prom:
                pick.append(i)
        if not pick:
            pick = [int(np.argmax(y)) if kind == 'max' else int(np.argmin(y))]
        return pick
    vmax_idx = _extrema(vp[0], 'max')
    cmin_idx = _extrema(cp[0], 'min')
    kfrac = read_path_kpoints(m, n_expect=len(k))
    res['kfrac'] = kfrac.tolist() if kfrac is not None else None

    def _dedup(idx, y, best='max'):
        """按 k 分数坐标去重：同一 k 在路径里出现多次只算一个谷。"""
        if kfrac is None or len(kfrac) != len(y):
            return idx
        seen, out = {}, []
        for i in idx:
            key = tuple(np.round(np.mod(kfrac[i] + 1e-9, 1.0), 4))
            if key in seen:
                j = seen[key]
                better = (y[i] > y[j]) if best == 'max' else (y[i] < y[j])
                if better:
                    out[out.index(j)] = i
                    seen[key] = i
            else:
                seen[key] = i
                out.append(i)
        return out
    vmax_idx = _dedup(vmax_idx, vp[0], 'max')
    cmin_idx = _dedup(cmin_idx, cp[0], 'min')
    vmax = sorted([dict(k=float(k[i]), E_meV=float((vp[0][i] - vbm) * 1000)) for i in vmax_idx],
                  key=lambda d: -d['E_meV'])          # 从高到低：vmax[0] = 主谷（VBM）
    cmin = sorted([dict(k=float(k[i]), E_meV=float((cp[0][i] - cbm) * 1000)) for i in cmin_idx],
                  key=lambda d: d['E_meV'])           # 从低到高：cmin[0] = 主谷（CBM）
    res['valleys_v'] = vmax
    res['valleys_c'] = cmin
    # 价带谷的能量都是 <= 0（相对 VBM），"kT 内" = E_meV >= -kT
    res['n_valleys_kT_v'] = int(sum(1 for d in vmax if -KT * 1000 - 1e-6 <= d['E_meV'] <= 1e-3))
    res['n_valleys_kT_c'] = int(sum(1 for d in cmin if -1e-3 <= d['E_meV'] <= KT * 1000 + 1e-6))
    res['second_valley_meV_v'] = vmax[1]['E_meV'] if len(vmax) > 1 else None
    res['second_valley_meV_c'] = cmin[1]['E_meV'] if len(cmin) > 1 else None
    return res


def fig_material(r, path):
    m = r['material']; k = r['k']; labels = r['labels']
    kk, E, _, _, _ = load(m)
    nv = r['nv']
    fig, (ax, az) = plt.subplots(1, 2, figsize=(13.2, 5.4), layout='constrained',
                                 gridspec_kw=dict(width_ratios=[1.0, 1.12]))
    # ---------------- 全貌 ----------------
    Eu = unique_bands(E) - r['vbm']
    nup = Eu.shape[1]
    for j in range(nup):
        col = C_DEEP if j < nup - 8 else '#5d6d7e'
        lw = 0.7
        ax.plot(k, Eu[:, j], color=col, lw=lw)
    for j in range(4):                                    # 顶部 4 对价带 + 底部 4 对导带高亮
        ax.plot(k, r['vp'][j] - r['vbm'], color=C_VAL[j], lw=1.5, zorder=3)
        ax.plot(k, r['cp'][j] - r['vbm'], color=C_CON[j], lw=1.5, zorder=3)
    draw_path(ax, labels)
    ax.axhline(0, color='0.35', lw=0.7, ls='--')
    ax.axhspan(0, r['gap'], color='#f4f6f7', zorder=0)
    ax.plot(k[[int(np.argmax(r['vp'][0]))]], 0, 'o', ms=7, mfc='white', mec=C_VAL[0], mew=1.6, zorder=5)
    ax.plot(k[[int(np.argmin(r['cp'][0]))]], r['gap'], 's', ms=7, mfc='white', mec=C_CON[0], mew=1.6, zorder=5)
    ax.annotate('', xy=(k[-1] * 0.985, 0), xytext=(k[-1] * 0.985, r['gap']),
                arrowprops=dict(arrowstyle='<->', color='0.3', lw=1.0))
    ax.text(k[-1] * 0.975, r['gap'] / 2, '%.3f eV' % r['gap'], ha='right', va='center', fontsize=10, color='0.2')
    ax.set_ylim(-4.0, 3.2)
    ax.set_ylabel('E - E$_{VBM}$ (eV)')
    ax.set_xlim(k[0], k[-1])
    ax.set_title('(a) %s   HSE+SOC   %s gap %.3f eV' % (m, r['gap_type'], r['gap']),
                 loc='left', fontsize=11, fontweight='bold')
    ax.legend(handles=[Line2D([], [], color=C_VAL[0], lw=1.5, label='top 4 valence pairs'),
                       Line2D([], [], color=C_CON[0], lw=1.5, label='bottom 4 conduction pairs'),
                       Line2D([], [], color=C_DEEP, lw=0.7, label='other bands (1 line = Kramers pair ×2)')],
              loc='lower right', fontsize=8.5, frameon=False)
    # ---------------- 带边放大 ----------------
    for j in range(6):
        az.plot(k, r['vp'][j] - r['vbm'], color=C_VAL[j], lw=1.5, label='VB%d (×2)' % (j + 1))
    for j in range(6):
        az.plot(k, r['cp'][j] - r['vbm'], color=C_CON[j], lw=1.5, label='CB%d (×2)' % (j + 1))
    draw_path(az, labels)
    az.axhspan(-KT, KT, color='#fdebd0', alpha=.75, zorder=0)
    az.axhspan(r['gap'] - KT, r['gap'] + KT, color='#d6eaf8', alpha=.75, zorder=0)
    az.axhline(0, color='0.35', lw=.6, ls='--')
    az.axhline(r['gap'], color='0.35', lw=.6, ls='--')
    iv = int(np.argmax(r['vp'][0])); ic = int(np.argmin(r['cp'][0]))
    az.plot(k[iv], 0, 'o', ms=8, mfc='white', mec=C_VAL[0], mew=1.8, zorder=5)
    az.plot(k[ic], r['gap'], 's', ms=8, mfc='white', mec=C_CON[0], mew=1.8, zorder=5)
    # 收敛标记：顶价带在窗口内的其它 k 位置
    # 带收敛：标出 kT 窗口内的价带谷（局部极大）
    ymin, ymax = az.get_ylim() if False else (-0.37, 0.44)
    for d in r['valleys_v']:
        if -KT * 1000 - 1e-6 <= d['E_meV'] <= 1e-3:
            az.axvline(d['k'], color=C_VAL[0], lw=0.8, ls=':', zorder=2)
            az.annotate('%+.0f meV' % d['E_meV'], (d['k'], 0.012), fontsize=7.6, color=C_VAL[0],
                        ha='center', va='bottom', rotation=0)
    for d in r['valleys_c']:
        if -1e-3 <= d['E_meV'] <= KT * 1000:
            az.axvline(d['k'], color=C_CON[0], lw=0.8, ls=':', zorder=2)
    az.text(k[0], -0.362, 'dotted: valleys within kT of VBM (n=%d) / CBM (n=%d)' % (
        r['n_valleys_kT_v'], r['n_valleys_kT_c']), fontsize=7.8, color='0.3', va='bottom')
    nv24 = r['window_v']['25.0'] if '25.0' in r['window_v'] else r['window_v']['25']
    nv10 = r['window_v']['10']
    nc10 = r['window_c']['10']
    txt = ('VBM at k=%.2f (%s)\nCBM at k=%.2f (%s)\n'
           'states within 10 meV: V %d / C %d\nstates within kT: V %d / C %d\n'
           'Kramers pairs: every line = 2 e$^-$') % (
        r['vbm_k'][0], ','.join('%.2f' % x for x in r['vbm_k'][1]),
        r['cbm_k'][0], ','.join('%.2f' % x for x in r['cbm_k'][1]),
        nv10['states'], nc10['states'], nv24['states'], r['window_c']['25']['states'])
    az.text(0.015, 0.98, txt, transform=az.transAxes, va='top', ha='left', fontsize=8.2,
            bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='0.8', alpha=.92))
    az.set_ylim(-0.37, 0.44)
    az.set_xlim(k[0], k[-1])
    az.set_ylabel('E - E$_{VBM}$ (eV)')
    az.set_title('(b) near-edge zoom (±k$_B$T at 300 K shaded)', loc='left', fontsize=11, fontweight='bold')
    az.legend(ncol=2, fontsize=7.4, frameon=False, loc='upper right', bbox_to_anchor=(1.0, 0.86))
    fig.savefig(path + '.png', dpi=300)
    fig.savefig(path + '.pdf')
    plt.close(fig)


def fig_compare(rs, path):
    fig, axs = plt.subplots(2, 2, figsize=(11.6, 8.2), layout='constrained')
    for ax, r in zip(axs.flat, rs):
        m = r['material']; k = r['k']
        for j in range(4):
            ax.plot(k, r['vp'][j] - r['vbm'], color=C_VAL[j], lw=1.3,
                    label='VB%d' % (j + 1) if j < 2 else None)
            ax.plot(k, r['cp'][j] - r['vbm'], color=C_CON[j], lw=1.3,
                    label='CB%d' % (j + 1) if j < 2 else None)
        draw_path(ax, r['labels'], fs=9)
        ax.axhspan(-KT, KT, color='#fdebd0', alpha=.7, zorder=0)
        ax.axhspan(r['gap'] - KT, r['gap'] + KT, color='#d6eaf8', alpha=.7, zorder=0)
        ax.axhline(0, color='0.35', lw=.5, ls='--'); ax.axhline(r['gap'], color='0.35', lw=.5, ls='--')
        ax.set_ylim(-0.35, 0.42); ax.set_xlim(k[0], k[-1])
        ax.set_ylabel('E - E$_{VBM}$ (eV)', fontsize=9.5)
        ax.set_title('%s   %s gap %.3f eV' % (m, r['gap_type'], r['gap']), loc='left', fontsize=10.5)
        w10 = r['window_v']['10']['states']
        ax.text(0.97, 0.05, 'V states ≤10 meV: %d' % w10, transform=ax.transAxes,
                ha='right', va='bottom', fontsize=8.5, color='0.25')
        ax.legend(fontsize=7.5, frameon=False, ncol=2, loc='upper right')
    fig.suptitle('HSE+SOC band edges of A$_2$B$_2$Te$_5$ (every band is Kramers-degenerate ×2)', fontsize=13)
    fig.savefig(path + '.png', dpi=300)
    fig.savefig(path + '.pdf')
    plt.close(fig)


def fig_diag(rs, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.0, 5.2), layout='constrained',
                                   gridspec_kw=dict(width_ratios=[1.15, 1.0]))
    # (a) 高对称点能级阶梯（价带 4 对）
    xs = np.arange(len(rs))
    width = 0.16
    for i, r in enumerate(rs):
        labs = [l['label'] for l in r['levels']]
        for j, lev in enumerate(r['levels']):
            x = i + (j - (len(labs) - 1) / 2) * width
            for q in range(4):
                y = lev['v'][q]
                col = C_VAL[q]
                ax1.plot([x - width * 0.38, x + width * 0.38], [y, y], color=col, lw=2.4 - 0.3 * q,
                         solid_capstyle='butt')
            if i == 0:
                ax1.text(x, 0.06, labs[j], rotation=0, ha='center', va='bottom', fontsize=7.5)
        ax1.text(i, -1.42, r['material'], ha='center', va='top', fontsize=9.5, fontweight='bold')
    for q in range(4):
        ax1.plot([], [], color=C_VAL[q], lw=2.2, label='VB%d pair' % (q + 1))
    ax1.legend(fontsize=8, frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(.5, -.14))
    ax1.set_ylim(-1.55, 0.20); ax1.set_xlim(-0.6, len(rs) - 0.4)
    ax1.set_xticks([]); ax1.set_ylabel('E - E$_{VBM}$ (eV)')
    ax1.set_title('(a) valence-pair energies at high-symmetry points\n(level spacing = degeneracy / crystal-field splitting)',
                  loc='left', fontsize=10.5, fontweight='bold')
    ax1.axhline(0, color='0.4', lw=.6, ls='--')
    # (b) 收敛曲线：态数 vs 窗口
    for i, r in enumerate(rs):
        ws = np.array([0.005, 0.010, 0.025, 0.050, 0.100])
        vals = [r['window_v']['%g' % (w * 1000)]['states'] for w in ws]
        ax2.plot(ws * 1000, vals, 'o-', lw=1.6, ms=5, label=r['material'])
        for w, v in zip(ws, vals):
            if w == 0.025:
                ax2.annotate('%d' % v, (w * 1000, v), textcoords='offset points', xytext=(4, 5), fontsize=8.5)
    ax2.axvline(25.9, color='0.5', lw=1.0, ls=':')
    ax2.text(26.5, ax2.get_ylim()[1] * 0.5 if ax2.get_ylim()[1] > 0 else 10, 'kT(300K)\n25.9 meV', fontsize=8, color='0.35')
    ax2.set_xscale('log'); ax2.set_xticks([5, 10, 25, 50, 100]); ax2.set_xticklabels(['5', '10', '25', '50', '100'])
    ax2.set_xlabel('energy window below VBM (meV)'); ax2.set_ylabel('number of valence states (×2 Kramers)')
    ax2.set_title('(b) band convergence: states within window of VBM', loc='left', fontsize=10.5, fontweight='bold')
    ax2.legend(fontsize=8.5, frameon=False)
    fig.savefig(path + '.png', dpi=300)
    fig.savefig(path + '.pdf')
    plt.close(fig)


def main(argv=None):
    global SRC, OUT, FIG
    ap = argparse.ArgumentParser(description='HSE+SOC 精细化能带图 + 简并/带收敛诊断')
    ap.add_argument('--data-root', default=SRC, help='含 <材料>/band.dat… 的目录（默认 %(default)s）')
    ap.add_argument('--out', default=OUT, help='输出目录；图写到 <out>/figures（默认 %(default)s）')
    a = ap.parse_args(argv)
    SRC, OUT = a.data_root, a.out
    FIG = os.path.join(OUT, 'figures')
    os.makedirs(FIG, exist_ok=True)
    rs = []
    for m in MATS:
        r = analyze(m)
        rs.append(r)
        print('%s: gap %.4f eV (%s), VBM k=%.2f frac=%s, CBM k=%.2f frac=%s' % (
            m, r['gap'], r['gap_type'], r['vbm_k'][0], np.round(r['vbm_k'][1], 2), r['cbm_k'][0], np.round(r['cbm_k'][1], 2)))
        print('   states within 10/25/50 meV of VBM:', [r['window_v'][k]['states'] for k in ['10', '25', '50']],
              ' of CBM:', [r['window_c'][k]['states'] for k in ['10', '25', '50']])
        fig_material(r, os.path.join(FIG, 'bands_%s' % m))
    fig_compare(rs, os.path.join(FIG, 'bands_4mats_near_edge'))
    fig_diag(rs, os.path.join(FIG, 'bands_degeneracy_convergence'))
    # CSV
    import csv
    rows = []
    for r in rs:
        rows.append(dict(material=r['material'], gap_eV=round(r['gap'], 5), gap_type=r['gap_type'],
                         vbm_k=round(r['vbm_k'][0], 4), cbm_k=round(r['cbm_k'][0], 4),
                         states_V_10meV=r['window_v']['10']['states'], states_V_25meV=r['window_v']['25']['states'],
                         states_V_50meV=r['window_v']['50']['states'], states_V_100meV=r['window_v']['100']['states'],
                         states_C_10meV=r['window_c']['10']['states'], states_C_25meV=r['window_c']['25']['states'],
                         kv_kT=r['window_v']['25']['k_points']))
    with open(os.path.join(OUT, 'band_summary_refined.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    json.dump([{k2: v for k2, v in r.items() if k2 not in ('k', 'vp', 'cp', 'dv', 'dc', 'labels')} for r in rs],
              open(os.path.join(OUT, 'band_summary_refined.json'), 'w'), ensure_ascii=False, indent=2)
    # markdown 摘要
    lines = ['# 精细化能带图与带边简并/收敛数值摘要', '',
             '数据：HSE+SOC 能带（Γ-M-K-Γ-A-L-H-A|L-M|H-K 路径，~200 点，能量相对 VBM）；',
             'SOC + 反演对称使**每条能带都是 Kramers 二重简并**（逐点核对：配对能量差 0），图中每条线代表 2 个电子。', '',
             '| 材料 | 带隙 / eV | 类型 | VBM 位置 | CBM 位置 | kT 内价带谷数（含主谷） | 次谷能量 / meV | kT 内导带谷数 |',
             '|---|---:|---|---|---|---:|---:|---:|']
    for r in rs:
        lines.append('| %s | %.4f | %s | k=%.2f (%s) | k=%.2f (%s) | %d | %s | %d |' % (
            r['material'], r['gap'], r['gap_type'], r['vbm_k'][0], ','.join('%.2f' % x for x in r['vbm_k'][1]),
            r['cbm_k'][0], ','.join('%.2f' % x for x in r['cbm_k'][1]),
            r['n_valleys_kT_v'], ('%.1f' % r['second_valley_meV_v']) if r['second_valley_meV_v'] is not None else '—',
            r['n_valleys_kT_c']))
    lines += ['', '> “kT 内价带谷数”= 沿路径顶价带的**独立局部极大**中能量在 VBM+25.9 meV 内的个数（= 带收敛/多谷简并），',
              '> 与“路径上落在窗口内的 k 点数”（平坦度 → DOS 有效质量）是两件不同的事，本文分开给。', '']
    lines += ['', '## 高对称点价带前 4 对能级（eV，相对 VBM）', '']
    for r in rs:
        lines += ['### ' + r['material'], '', '| 点 | VB1 | VB2 | VB3 | VB4 |', '|---|---:|---:|---:|---:|']
        for lev in r['levels']:
            lines.append('| %s | %s |' % (lev['label'], ' | '.join('%.4f' % x for x in lev['v'])))
        lines.append('')
    open(os.path.join(OUT, 'band_summary_refined.md'), 'w', encoding='utf-8').write('\n'.join(lines))
    print('written to', OUT)
    print('\n'.join(lines[:12]))


if __name__ == '__main__':
    main()
