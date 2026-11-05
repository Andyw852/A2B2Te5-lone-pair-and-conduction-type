#!/usr/bin/env python3
"""HSE06+SOC 密网格 BoltzTraP2（c8_hse_dense，28×28×7/6，mult 60）的输运图与表（2026-11-05）。
用户要求：原来只有 PBE 的图，每张都给一张 HSE 单独图和一张 HSE 与 PBE 的对比图。

  图 8-1  300 K 的 S 随载流子浓度（HSE06+SOC）               图 8-2  同上，HSE 与 PBE 对比（每种材料一个子图）
  图 9-1  霍尔浓度下的 S(T)：实验与 HSE06+SOC                图 9-2  实验、HSE、PBE 三者
  图 10-1 HSE06+SOC 能带上 300 K 的化学势                     图 10-2 HSE 与 PBE 能带及各自霍尔浓度下的化学势
  表 6、表 7、表 8 的 HSE 与 PBE 数值（tables_hse_pbe）

数据：22_能带基准_1-0/<材料>/c8_hse_dense/gene/bt2_mult60/boltztrap_crta.json（不加剪刀，常数弛豫时间），
PBE 对照仍为 bt2_pbe_dense（paper_figures.pbe_dense）。浓度方向在化学势网格上按 log|N| 线性取值，温度不插值。
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

MATS, TEX, NEL, KEND, r, dump, save, panel = pf.MATS, pf.TEX, pf.NEL, pf.KEND, pf.r, pf.dump, pf.save, pf.panel
C_HSE, C_PBE, DASH, W2 = pf.C_HSE, pf.C_PBE, pf.DASH, pf.W2
MULT = 60


def hse_json(m, mult=MULT):
    return os.path.join(pf.ROOT, '22_能带基准_1-0', m, 'c8_hse_dense', 'gene', f'bt2_mult{mult}', 'boltztrap_crta.json')


def load(path):
    """BoltzTraP2 结果：T、f(N, i)（第 i 个温度、带符号浓度 N 处的面内 S，µV/K）与原始的 µ、浓度网格。"""
    d = json.load(open(path, encoding='utf-8'))
    T = np.array(d['temperatures_K']); C = np.array(d['carrier_conc_cm-3'], float); mu = np.array(d['mu_rel_efermi_eV'], float)
    S = (np.array(d['seebeck_xx_V_per_K'], float) + np.array(d['seebeck_yy_V_per_K'], float)) / 2 * 1e6

    def f(N, i):
        sg = 1 if N > 0 else -1; ok = sg * C[i] > 0; x = np.log10(sg * C[i][ok]); o = np.argsort(x)
        assert x.min() <= np.log10(abs(N)) <= x.max()
        return float(np.interp(np.log10(abs(N)), x[o], S[i][ok][o]))

    def mu_at(N, i):
        sg = 1 if N > 0 else -1; ok = sg * C[i] > 0; x = np.log10(sg * C[i][ok]); o = np.argsort(x)
        return float(np.interp(np.log10(abs(N)), x[o], mu[ok][o]))
    return T, f, mu_at, d


def hse_dense(m):
    T, f, _, d = load(hse_json(m))
    assert d['functional'] == 'HSE06+SOC' and d['kmesh_mult'] == MULT and 'no scissor' in d['protocol'], (m, d['protocol'])
    return T, f


LOAD = {'HSE': hse_dense, 'PBE': pf.pbe_dense}
FN = {'HSE': 'HSE06+SOC', 'PBE': 'PBE+SOC'}
COL = {'HSE': C_HSE, 'PBE': C_PBE}
NS = np.logspace(18, np.log10(3e21), 73)


# ------------------------------------------------------------------ 图 8-1、8-2
def fig8_1(out):
    """300 K 面内 S 随载流子浓度，HSE06+SOC（与原图 8 的 PBE 版同一画法）。"""
    T0 = 300; rows = []
    fig, A = plt.subplots(1, 2, figsize=(W2, 2.75), gridspec_kw=dict(wspace=0.26))
    for i, (ax, sign) in enumerate(zip(A, (1, -1))):
        for m in MATS:
            T, f = hse_dense(m); it = list(T).index(T0)
            v = np.array([sign * f(sign * N, it) for N in NS])
            ax.plot(NS, np.where(v > 0, v, np.nan), color=pf.MCOL[m], lw=1.1, marker=pf.MMARK[m], ms=4.2,
                    markevery=(3 * MATS.index(m) % 6, 6), label=TEX[m], zorder=3)
            if (v < 0).any(): ax.plot(NS, np.where(v < 0, -v, np.nan), color=pf.MCOL[m], lw=1.1, ls=(0, (1, 1.4)), zorder=3)
            rows += [[m, 'pn'[i], f'{N:.3e}', r(x, 2)] for N, x in zip(NS, v)]
        style_sn(ax, sign, f'{"p" if sign > 0 else "n"}-type, 300 K, HSE06+SOC + BoltzTraP2'); panel(ax, f'({"ab"[i]})')
    A[0].legend(*A[0].get_legend_handles_labels(), loc='lower left', fontsize=6.5, handlelength=2.2, labelspacing=0.3)
    if any(float(x[3]) < 0 for x in rows):
        A[1].text(0.34, 0.04, 'dotted: $S>0$ despite net electrons\n(bipolar, holes dominate $S$)', transform=A[1].transAxes,
                  ha='left', va='bottom', fontsize=6.3, color='#222222', linespacing=1.15, bbox=dict(fc='white', ec='none', pad=0.5))
    dump('Fig8-1_Seebeck_HSE', 'S_vs_concentration_300K_HSE.csv', ['material', 'carrier_type', 'net_concentration_cm-3',
         'S_HSE06+SOC_BoltzTraP2_uV/K, sign of majority carrier positive (negative = S has changed sign, bipolar)'], rows)
    save(fig, out, 'Fig8-1_Seebeck_vs_concentration_HSE')


def style_sn(ax, sign, title):
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(NS[0], NS[-1]); ax.set_ylim(3, 600)
    ax.set_yticks([5, 10, 20, 50, 100, 200, 500], ['5', '10', '20', '50', '100', '200', '500']); ax.minorticks_off()
    ax.set_xlabel(('Net hole' if sign > 0 else 'Net electron') + r' concentration (cm$^{-3}$)' if sign else r'Net carrier concentration (cm$^{-3}$)')
    ax.set_ylabel(r'$S$ (μV K$^{-1}$)' if sign > 0 else r'$|S|$ (μV K$^{-1}$)')
    ax.set_title(title, fontsize=8, pad=4); ax.axvline(3e19, color='#999999', lw=0.6, ls=(0, (1, 1.5)), zorder=0)


def fig8_2(out):
    """HSE 与 PBE 对比：每种材料一个子图，蓝 HSE、橙 PBE；实线 p 型、虚线 n 型（|S|）；变号的一段为细点线；灰点为样品霍尔浓度。"""
    T0 = 300; rows = []
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.0), gridspec_kw=dict(hspace=0.42, wspace=0.26))
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        for fn in ('PBE', 'HSE'):
            T, f = LOAD[fn](m); it = list(T).index(T0)
            for sign, ls in ((1, '-'), (-1, (0, (4, 1.6)))):
                v = np.array([sign * f(sign * N, it) for N in NS])
                ax.plot(NS, np.where(v > 0, v, np.nan), color=COL[fn], lw=1.3 if fn == 'HSE' else 1.0, ls=ls, zorder=3 if fn == 'HSE' else 2)
                if (v < 0).any(): ax.plot(NS, np.where(v < 0, -v, np.nan), color=COL[fn], lw=0.8, ls=(0, (1, 1.4)), zorder=2)
                rows += [[m, fn, 'pn'[sign < 0], f'{N:.3e}', r(x, 2)] for N, x in zip(NS, v)]
        h = pf.HALL[m]; ax.axvline(abs(h['N']), color='#BBBBBB', lw=0.8, zorder=0)
        ax.text(abs(h['N']) * 0.9, 4, 'Hall', rotation=90, ha='right', va='bottom', fontsize=6.0, color='#777777')
        style_sn(ax, 0, f'{TEX[m]}, 300 K'); ax.set_ylabel(r'$|S|$ (μV K$^{-1}$)'); panel(ax, f'({"abcd"[n]})')
    hh = [Line2D([], [], color=C_HSE, lw=1.3), Line2D([], [], color=C_PBE, lw=1.0), Line2D([], [], color='#555555', lw=1.0),
          Line2D([], [], color='#555555', lw=1.0, ls=(0, (4, 1.6))), Line2D([], [], color='#555555', lw=0.8, ls=(0, (1, 1.4)))]
    fig.legend(hh, ['HSE06+SOC', 'PBE+SOC', 'p-type', 'n-type (|S|)', 'S has changed sign (bipolar)'], loc='lower center', ncol=5,
               bbox_to_anchor=(0.5, 0.935), columnspacing=1.2, handlelength=2.4)
    dump('Fig8-2_Seebeck_HSE_vs_PBE', 'S_vs_concentration_300K_HSE_and_PBE.csv', ['material', 'functional', 'carrier_type', 'net_concentration_cm-3',
         'S_uV/K, sign of majority carrier positive (negative = S has changed sign, bipolar)'], rows)
    save(fig, out, 'Fig8-2_Seebeck_vs_concentration_HSE_vs_PBE')


# ------------------------------------------------------------------ 图 9-1、9-2
YL = {'Pb2Sb2Te5': (0, 80), 'Pb2Bi2Te5': (-175, 0), 'Sn2Sb2Te5': (0, 60), 'Sn2Bi2Te5': (0, 90)}


def st_curve(m, fn):
    T, f = LOAD[fn](m); N = pf.HALL[m]['N']; v = np.array([f(N, i) for i in range(len(T))])
    dev = [float(np.interp(t, T, v)) - s for t, s in pf.EXP[m]]
    return T, v, dev, float(np.sqrt(np.mean(np.square(dev))))


def fig9(out, fns, name, folder):
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.9), gridspec_kw=dict(hspace=0.42, wspace=0.24))
    rows = []
    for n, (ax, m) in enumerate(zip(axs.flat, MATS)):
        h = pf.HALL[m]; N = h['N']; Te, Se = zip(*pf.EXP[m])
        ax.plot(Te, Se, 'o', ms=4.5, mfc='white', mec='k', mew=0.9, ls='', zorder=5)
        y = 0.865 if N > 0 else 0.05 + 0.085 * (len(fns) - 1)          # n 型（纵轴为负）文字放左下角，避开曲线
        for fn in fns:
            T, v, dev, rms = st_curve(m, fn)
            ax.plot(T, v, color=COL[fn], lw=1.4, ls='-' if fn == 'HSE' else DASH, zorder=4 if fn == 'HSE' else 3)
            ax.text(0.04, y, f'{FN[fn].split("+")[0]}: rms {rms:.0f} μV K$^{{-1}}$', transform=ax.transAxes, ha='left', va='top' if N > 0 else 'bottom',
                    fontsize=6.5, color=COL[fn])
            y -= 0.085
            rows += [[m, FN[fn], f'{N:.3e}', int(t), r(x, 2)] for t, x in zip(T, v)]
        ax.set_xlim(280, 660); ax.set_ylim(*YL[m]); ax.set_title(f'{TEX[m]} ({"n" if N < 0 else "p"}-type)', fontsize=8, pad=4)
        ax.text(0.04, 0.95, f'Hall: {"$p$" if N > 0 else "$n$"} = {pf.sci(abs(N))}' + (' (flagged)' if h['flagged'] else ''),
                transform=ax.transAxes, ha='left', va='top', fontsize=6.5, color='#555555')
        panel(ax, f'({"abcd"[n]})')
    for ax in axs[:, 0]: ax.set_ylabel(r'$S$ (μV K$^{-1}$)')
    for ax in axs[1]: ax.set_xlabel('$T$ (K)')
    hh = [Line2D([], [], ls='', marker='o', ms=4.5, mfc='white', mec='k', mew=0.9)] + \
         [Line2D([], [], color=COL[fn], lw=1.4, ls='-' if fn == 'HSE' else DASH) for fn in fns]
    fig.legend(hh, ['Experiment'] + [f'{FN[fn]} + BoltzTraP2 at the Hall concentration' for fn in fns], loc='lower center', ncol=1 + len(fns),
               bbox_to_anchor=(0.5, 0.935), columnspacing=1.2, handlelength=1.8, fontsize=6.6)
    dump(folder, 'experiment_T_sigma_S.csv', ['material', 'T_K', 'sigma_exp_S/m', 'S_exp_uV/K'],
         [[m, r(t, 2), r(sg, 1), r(x, 3)] for m in MATS for t, sg, x in pf.read_exp(m)])
    dump(folder, 'hall_room_temperature.csv', ['material', 'carrier_conc_signed_cm-3 (+ holes, - electrons)', 'Hall_mobility_cm2/Vs', 'flagged_by_instrument'],
         [[k, f'{x["N"]:.3e}', x['mu'], 'yes' if x['flagged'] else 'no'] for k, x in pf.HALL.items()])
    dump(folder, 'calculated_S_T_at_Hall_conc.csv', ['material', 'functional', 'Hall_conc_signed_cm-3', 'T_K', 'S_BoltzTraP2_uV/K (no scissor)'], rows)
    save(fig, out, name)


def fig9_1(out): fig9(out, ('HSE',), 'Fig9-1_S_vs_T_Hall_HSE', 'Fig9-1_S_vs_T_HSE')
def fig9_2(out): fig9(out, ('HSE', 'PBE'), 'Fig9-2_S_vs_T_Hall_HSE_vs_PBE', 'Fig9-2_S_vs_T_HSE_vs_PBE')


# ------------------------------------------------------------------ 表 6、7、8（HSE 与 PBE）
def tables_hse_pbe(out):
    T0 = 300; t6, t7, t8 = [], [], []
    for m in MATS:
        row = [m]
        for fn in ('HSE', 'PBE'):
            T, f = LOAD[fn](m); it = list(T).index(T0)
            row += [r(f(N, it), 1) for N in (3e19, 1e20, -3e19, -1e20)]
        t6.append(row)
        h = pf.HALL[m]; d = pf.EXP[m]; rw = [m, 'p' if h['N'] > 0 else 'n', f'{h["N"]:.3e}', 'yes' if h['flagged'] else 'no', r(d[0][0], 1), r(d[0][1], 1)]
        for fn in ('HSE', 'PBE'):
            T, v, dev, rms = st_curve(m, fn); rw += [r(np.interp(d[0][0], T, v), 1)]
        rw += [r(d[-1][0], 1), r(d[-1][1], 1)]
        for fn in ('HSE', 'PBE'):
            T, v, dev, rms = st_curve(m, fn); rw += [r(np.interp(d[-1][0], T, v), 1)]
        for fn in ('HSE', 'PBE'):
            rw += [r(st_curve(m, fn)[3], 1)]
        t7.append(rw)
        n0 = abs(h['N']); S0 = abs(d[0][1]); comp = m == 'Pb2Sb2Te5'; rm = [m, f'{h["N"]:.3e}', r(d[0][1], 2), '' if comp else r(pf.spb_mass(S0, n0, T0), 3)]
        for fn in ('HSE', 'PBE'):
            T, f = LOAD[fn](m); s = f(h['N'], list(T).index(T0)); rm += [r(s, 1), r(pf.spb_mass(abs(s), n0, T0), 3)]
        rm += [pf.SINGLE[m][0]]; t8.append(rm)
    dump('Table6_S_300K', 'S_300K_HSE_and_PBE.csv', ['material'] + [f'{fn}_{c}' for fn in ('HSE', 'PBE') for c in
         ('S_p_3e19_uV/K', 'S_p_1e20_uV/K', 'S_n_3e19_uV/K', 'S_n_1e20_uV/K')], t6)
    dump('Table7_S_T_at_Hall', 'S_T_at_Hall_HSE_and_PBE.csv', ['material', 'type', 'Hall_conc_signed_cm-3', 'Hall_flagged', 'T_low_K', 'S_exp_low_uV/K',
         'S_HSE_low_uV/K', 'S_PBE_low_uV/K', 'T_high_K', 'S_exp_high_uV/K', 'S_HSE_high_uV/K', 'S_PBE_high_uV/K', 'rms_HSE_uV/K', 'rms_PBE_uV/K'], t7)
    dump('Table8_Pisarenko_mass', 'samples_Pisarenko_mass_HSE_and_PBE.csv', ['material', 'Hall_conc_signed_cm-3', 'S_exp_300K_uV/K', 'm_DOS_exp_SPB_me',
         'S_HSE_at_Hall_uV/K', 'm_DOS_HSE_SPB_me', 'S_PBE_at_Hall_uV/K', 'm_DOS_PBE_SPB_me', 'single_valley_m_DOS_HSE_band_fit_me'], t8)
    return t6, t7, t8


# ------------------------------------------------------------------ 图 10-1、10-2：能带上的化学势
def gene_gamma(m):
    """HSE GENE 网格（c8_hse_dense/gene）：Γ 处最高价带、最低导带与 GENE 的费米能（BoltzTraP2 化学势的零点），绝对能量 eV。"""
    g = os.path.join(pf.ROOT, '22_能带基准_1-0', m, 'c8_hse_dense', 'gene')
    src = json.load(open(os.path.join(g, 'source.json'), encoding='utf-8'))
    L = open(os.path.join(g, f'{m}.energy')).read().split('\n'); nk, _, ef = L[1].split(); i = 2; RY = 13.605693122994
    iv = src['vbm_band_written_1based'] - 1
    for _ in range(int(nk)):
        p = L[i].split(); nb = int(p[3]); k = np.array([float(x) for x in p[:3]])
        e = [float(L[i + 1 + b]) * RY for b in range(nb)]
        if np.allclose(k, 0): return dict(vb_g=e[iv], cb_g=e[iv + 1], efermi=float(ef) * RY)
        i += 1 + nb
    raise SystemExit('[ERROR] %s 的 GENE 里没有 Γ' % m)


def bands_fn(m, fn):
    """路径能带（零点：本泛函的真实价带顶，与图 2 相同）、附加路径、阳离子 s 成分、Γ 对齐用的偏移、化学势函数。"""
    x, _ = pf.rf.path_x(m); nel = NEL[m]; side = 'c2_hse' if fn == 'HSE' else 'c2_pbe'
    E, _ = pf.bands(m, side); tp = pf.true_path_bands(m, side)
    if fn == 'HSE':
        T5 = pf.mv.true_load(m); dv = (T5['vbm'] - T5['path_vbm']) if T5 else 0.0
    else:
        dv = max(0.0, float(np.nanmax(tp['E'][:, nel - 1]))) if tp else 0.0
        if tp: tp['E'] = tp['E'] - dv
    E = E - dv
    if tp: tp['x'] = pf.TP_XSCALE * tp['x']
    _, Wp = pf.cv.read_procar(os.path.join(pf.rf.base(m), side, 'PROCAR_OPT.gz')); fs = pf.rf.cation_s(Wp)
    if fn == 'HSE': G = gene_gamma(m); T, _, mu_at, _ = load(hse_json(m))
    else: G = pf.dense_gamma(m); T, _, mu_at, _ = load(os.path.join(pf.PBEDENSE, m, 'boltztrap_crta.json'))
    o_v, o_c = E[0, nel - 1] - G['vb_g'], E[0, nel] - G['cb_g']; off = (o_v + o_c) / 2
    it = list(T).index(300)
    mu = lambda N: mu_at(N, it) + G['efermi'] + off
    gap = min(np.nanmin(E[:KEND + 1, nel]), np.nanmin(tp['E'][:, nel]) if tp else np.inf)
    return dict(x=x, E=E, tp=tp, fs=fs, nel=nel, mu=mu, gap=gap, resid=o_v - o_c)


CMAP = LinearSegmentedColormap.from_list('s', ['#DCE8F4', '#9CC3E4', '#4A90C8', '#1F5FA0', '#08306B'])


def fig10_1(out):
    """HSE06+SOC 能带（阳离子 s 成分）上 300 K 的化学势：样品霍尔浓度（橙）与 3×10¹⁹ cm⁻³（灰），与原图 10 的 PBE 版同一画法。"""
    norm = Normalize(0, 0.30); kT = pf.mv.KB * 300; rows = []
    fig = plt.figure(figsize=(W2, 5.2)); gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.30); every = []
    for n, m in enumerate(MATS):
        B = bands_fn(m, 'HSE'); x, E, tp, fs, nel = B['x'], B['E'], B['tp'], B['fs'], B['nel']
        ax, at = pf.band_axes(fig, gs[n // 2, n % 2], tp, x); every += [a for a in (ax, at) if a is not None]
        for b in range(nel - 14, nel + 6):
            cc = (fs[:KEND, b] + fs[1:KEND + 1, b]) / 2
            lc = LineCollection([[(x[i], E[i, b]), (x[i + 1], E[i + 1, b])] for i in range(KEND)], cmap=CMAP, norm=norm,
                                linewidths=0.5 + 2.4 * np.clip(cc / 0.3, 0, 1), capstyle='round', zorder=3); lc.set_array(cc); ax.add_collection(lc)
            if tp:
                xs, es, fc = tp['x'], tp['E'][:, b], tp['fs'][:, b]; c2 = (fc[:-1] + fc[1:]) / 2
                lc = LineCollection([[(xs[i], es[i]), (xs[i + 1], es[i + 1])] for i in range(len(xs) - 1)], cmap=CMAP, norm=norm,
                                    linewidths=0.5 + 2.4 * np.clip(c2 / 0.3, 0, 1), capstyle='round', zorder=3); lc.set_array(c2); at.add_collection(lc)
        NH = pf.HALL[m]['N']; NL = np.sign(NH) * pf.N_REF
        for N, col, lw, lab in ((NL, pf.C_WL, 1.0, 'reference'), (NH, pf.C_WH, 1.4, 'Hall')):
            u = B['mu'](N)
            for a in (ax, at):
                if a is not None: a.axhline(u, color=col, lw=lw, ls=DASH, zorder=4)
            depth = (u - B['gap']) if N < 0 else -u
            rows.append([m, lab, f'{N:.3e}', r(1000 * u, 1), r(1000 * depth, 1), r(depth / kT, 1), r(1000 * B['gap'], 1), r(1000 * B['resid'], 2)])
            ax.text(0.985, u + (0.012 if (N < 0) == (lab == 'Hall') else -0.012), f'{"$p$" if N > 0 else "$n$"} = {pf.sci(abs(N), unit=False)}' + (' (Hall)' if lab == 'Hall' else ''),
                    transform=ax.get_yaxis_transform(), ha='right', va='bottom' if (N < 0) == (lab == 'Hall') else 'top', fontsize=6.3, zorder=6,
                    bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.85))
        pf.kticks(ax, x); ax.set_ylim(-0.62, 0.52); ax.set_yticks(np.arange(-0.5, 0.51, 0.25)); ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)', labelpad=2)
        if tp: pf.tp_frame(at, tp, nel, mark=())
        ax.set_title(f'{TEX[m]}     $E_g$ = {1000 * B["gap"]:.0f} meV', fontsize=8, pad=4); panel(ax, f'({"abcd"[n]})')
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), ax=every, shrink=0.7, pad=0.015, aspect=30)
    cb.set_label('Cation-s weight (%)'); cb.set_ticks([0, .1, .2, .3], labels=['0', '10', '20', '30'])
    cb.ax.tick_params(direction='in', length=2, width=0.6); cb.outline.set_linewidth(0.6)
    fig.legend([Line2D([], [], color=pf.C_WH, lw=1.4, ls=DASH), Line2D([], [], color=pf.C_WL, lw=1.0, ls=DASH)],
               [r'$\mu$ at the sample Hall concentration (300 K)', r'$\mu$ at $3\times10^{19}$ cm$^{-3}$, same carrier type'],
               loc='lower center', ncol=2, bbox_to_anchor=(0.45, 0.925), columnspacing=1.6, handlelength=2.6)
    dump('Fig10-1_Fermi_window_HSE', 'chemical_potential_300K_HSE.csv', ['material', 'case', 'carrier_conc_signed_cm-3', 'mu-EVBM_meV',
         'depth_into_band_meV (p: below VBM, n: above CBM)', 'depth_in_kBT', 'Eg_HSE06+SOC_meV', 'alignment_residual_VB-CB_meV'], rows)
    save(fig, out, 'Fig10-1_Fermi_window_HSE')
    return rows


def fig10_2(out):
    """HSE（蓝实线）与 PBE（橙虚线）能带，各以本泛函的价带顶为零点；水平虚线为各自在样品霍尔浓度下 300 K 的化学势。"""
    rows = []
    fig = plt.figure(figsize=(W2, 5.2)); gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
    for n, m in enumerate(MATS):
        Bh, Bp = bands_fn(m, 'HSE'), bands_fn(m, 'PBE'); x, nel = Bh['x'], Bh['nel']; NH = pf.HALL[m]['N']
        ax, at = pf.band_axes(fig, gs[n // 2, n % 2], Bh['tp'], x)
        for fn, B, ls, lw, z in (('PBE', Bp, DASH, 0.9, 2), ('HSE', Bh, '-', 1.1, 3)):
            for b in range(nel - 10, nel + 8):
                ax.plot(x[:KEND + 1], B['E'][:KEND + 1, b], color=COL[fn], ls=ls, lw=lw, zorder=z)
                if B['tp'] is not None and at is not None: at.plot(B['tp']['x'], B['tp']['E'][:, b], color=COL[fn], ls=ls, lw=lw, zorder=z)
            u = B['mu'](NH)
            for a in (ax, at):
                if a is not None: a.axhline(u, color=COL[fn], lw=1.6, ls=(0, (5, 2)), zorder=5, alpha=0.9)
            rows.append([m, FN[fn], f'{NH:.3e}', r(1000 * u, 1), r(1000 * B['gap'], 1)])
        for a in (ax, at):
            if a is not None: a.axhline(0, color='k', lw=0.5, ls=':', zorder=1)
        pf.kticks(ax, x); ax.set_ylim(-0.62, 0.52); ax.set_yticks(np.arange(-0.5, 0.51, 0.25)); ax.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)', labelpad=2)
        if Bh['tp']: pf.tp_frame(at, Bh['tp'], nel, mark=())
        ax.set_title(f'{TEX[m]}  (Hall: {"$p$" if NH > 0 else "$n$"} = {pf.sci(abs(NH), unit=False)})', fontsize=8, pad=4); panel(ax, f'({"abcd"[n]})')
    fig.legend([Line2D([], [], color=C_HSE, lw=1.1), Line2D([], [], color=C_PBE, lw=0.9, ls=DASH),
                Line2D([], [], color='#555555', lw=1.6, ls=(0, (5, 2)))],
               ['HSE06+SOC bands', 'PBE+SOC bands', r'$\mu$ at the Hall concentration, 300 K (colour = functional)'],
               loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.925), columnspacing=1.4, handlelength=2.6)
    dump('Fig10-2_Fermi_window_HSE_vs_PBE', 'chemical_potential_300K_Hall_HSE_and_PBE.csv',
         ['material', 'functional', 'Hall_conc_signed_cm-3', 'mu-EVBM_meV (own VBM)', 'Eg_meV'], rows)
    save(fig, out, 'Fig10-2_Fermi_window_HSE_vs_PBE')
    return rows


FIGS = {'8-1': fig8_1, '8-2': fig8_2, '9-1': fig9_1, '9-2': fig9_2, '10-1': fig10_1, '10-2': fig10_2, 'T': tables_hse_pbe}

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=os.path.join(pf.ROOT, '21_阶段报告与交接', '论文图_20261028'))
    ap.add_argument('--data'); ap.add_argument('--only', nargs='+'); a = ap.parse_args()
    pf.DATA = a.data or os.path.join(a.out, 'source_data'); pf.style(); os.makedirs(a.out, exist_ok=True)
    for k in a.only or FIGS:
        res = FIGS[k](a.out)
        if k == 'T':
            for t in res:
                for row in t: print(row)
        elif res:
            for row in res: print(row)
