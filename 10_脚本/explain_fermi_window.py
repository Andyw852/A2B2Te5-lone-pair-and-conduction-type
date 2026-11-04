#!/usr/bin/env python3
"""讲解图（2026-11-04，回答用户“黄线远离价带顶，怎么和谷联系起来”）：以一个材料为例，把图 10 的一个子图拆成三栏，
纵轴都是能量，横着对齐读：
  (a) PBE+SOC 能带：µ（霍尔浓度，橙色虚线）以上的部分涂色 = 装着空穴的口袋（沿这条路径的截面）；
      µ 与能带的交点 = 费米面（输运发生的地方），交点处能带的斜率 = ħv；箭头 = µ 在谷顶以下的深度 D。
  (b) BoltzTraP2 的空穴浓度 p(µ)：µ 以上到价带顶的空态全部计入，所以由霍尔浓度定出 µ 的位置。
  (c) BoltzTraP2 的 σ/τ(µ)：S = −(π²/3)(k_B/e)k_BT·d ln σ/dE，切线斜率给出 S，与 BoltzTraP2 直接算的 S 对照。
数据全部取自 bt2_pbe_dense 与 c2_pbe 路径（同图 10，零点对齐方法也相同），不用模型。只作讲解，不进报告。

用法：python 10_脚本/explain_fermi_window.py [--material Sn2Bi2Te5] [--out 21_阶段报告与交接/讲解图_20261104]
"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import matplotlib.pyplot as plt
from matplotlib import font_manager

C_H, C_L, C_FILL, C_BAND = '#E69F00', '#555555', '#F6C66B', '#3A3A3A'


def cjk():
    for f in font_manager.findSystemFonts():
        if 'wqy' in f.lower() or 'wenquanyi' in f.lower():
            font_manager.fontManager.addfont(f)
    have = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams['font.family'] = [n for n in ('Times New Roman', 'Liberation Serif') if n in have][:1] + ['WenQuanYi Zen Hei', 'DejaVu Serif']


def usci(N):
    """9.1×10²⁰ 这样的写法，不用 mathtext（mathtext 里的中文会变成方框）。"""
    a, e = f'{abs(N):.1e}'.split('e')
    return f'{a}×10' + ''.join('⁰¹²³⁴⁵⁶⁷⁸⁹'[int(c)] for c in str(int(e)))


def load(m):
    """与 paper_figures.fig10 相同：路径能带（零点为两段路径上的最高价带）与对齐到同一零点的 BoltzTraP2 µ、p、σ/τ、S（300 K）。"""
    x, _ = pf.rf.path_x(m); nel = pf.NEL[m]; X0 = x[pf.KEND]
    E, _ = pf.bands(m, 'c2_pbe'); tp = pf.true_path_bands(m, 'c2_pbe')
    dv = max(0.0, float(np.nanmax(tp['E'][:, nel - 1]))) if tp else 0.0
    E = E - dv
    xs, Es = x[:pf.KEND + 1], E[:pf.KEND + 1]
    if tp:
        tp['E'] = tp['E'] - dv; tp['x'] = pf.TP_XSCALE * tp['x']
    D = pf.dense_gamma(m); off = ((E[0, nel - 1] - D['vb_g']) + (E[0, nel] - D['cb_g'])) / 2
    bt = json.load(open(os.path.join(pf.PBEDENSE, m, 'boltztrap_crta.json'), encoding='utf-8'))
    it = bt['temperatures_K'].index(300)
    mu = np.array(bt['mu_rel_efermi_eV']) + D['efermi'] + off
    C = np.array(bt['carrier_conc_cm-3'][it], float)
    sig = (np.array(bt['sigma_over_tau_xx'][it], float) + np.array(bt['sigma_over_tau_yy'][it], float)) / 2
    S = (np.array(bt['seebeck_xx_V_per_K'][it], float) + np.array(bt['seebeck_yy_V_per_K'][it], float)) / 2 * 1e6
    return dict(x=x, xs=xs, Es=Es, X0=X0, tp=tp, nel=nel, mu=mu, C=C, sig=sig, S=S)


def at(d, N):
    """带符号浓度 N 处的 µ、σ/τ、S、d ln σ/dE（按 log|N| 在化学势网格上线性取值）。"""
    sg = np.sign(N); ok = sg * d['C'] > 0; xx = np.log10(sg * d['C'][ok]); o = np.argsort(xx); q = np.log10(abs(N))
    dl = np.gradient(np.log(d['sig']), d['mu'])
    f = lambda y: float(np.interp(q, xx[o], y[ok][o]))
    return dict(mu=f(d['mu']), sig=f(d['sig']), S=f(d['S']), dlns=f(dl))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--material', default='Sn2Bi2Te5')
    ap.add_argument('--out', default=os.path.join(pf.ROOT, '21_阶段报告与交接', '讲解图_20261104')); a = ap.parse_args()
    pf.style(); cjk(); m = a.material; d = load(m); nel = d['nel']; kT = pf.mv.KB * 300
    NH = pf.HALL[m]['N']; NL = np.sign(NH) * 3e19; H, L = at(d, NH), at(d, NL); u = H['mu']
    assert NH > 0, '讲解图按 p 型写（空穴口袋）'

    fig, (A, B, Cx) = plt.subplots(1, 3, figsize=(7.4, 3.6), sharey=True, gridspec_kw=dict(width_ratios=[3.1, 1, 1], wspace=0.08))
    ylo, yhi = -0.45, 0.15

    # ---------- (a) 能带：µ 以上的口袋涂色，交点为费米面，箭头为深度 D
    segs = [(d['xs'], d['Es'])] + ([(d['X0'] + d['tp']['x'], d['tp']['E'])] if d['tp'] else [])
    tops, ncross, npocket = [], 0, 0
    for xx, EE in segs:
        for b in range(nel - 16, nel + 4):
            e = EE[:, b]
            A.plot(xx, e, color=C_BAND, lw=0.8, zorder=3)
            if b >= nel or (nel - 1 - b) % 2: continue              # SOC 下每条带是 Kramers 二重简并的一对，只数一条
            A.fill_between(xx, u, e, where=e > u, interpolate=True, color=C_FILL, alpha=0.45, lw=0, zorder=2)
            above = e > u; npocket += int(np.sum(above[1:] & ~above[:-1]) + above[0])
            for i in np.where(above[1:] != above[:-1])[0]:          # 与 µ 的交点（线性内插）
                xc = xx[i] + (u - e[i]) / (e[i + 1] - e[i]) * (xx[i + 1] - xx[i]); ncross += 1
                A.plot(xc, u, 'o', ms=3.2, mfc='white', mec='#B36B00', mew=0.9, zorder=5)
            for i in range(1, len(e) - 1):                         # µ 以上的谷顶（局部极大）
                if e[i] > u + 0.004 and e[i] >= e[i - 1] and e[i] >= e[i + 1]:
                    tops.append((e[i], xx[i]))
    A.axhline(u, color=C_H, lw=1.4, ls=pf.DASH, zorder=4); A.axhline(L['mu'], color=C_L, lw=1.0, ls=pf.DASH, zorder=4)
    pf.kticks(A, d['x']); A.set_ylim(ylo, yhi)
    if d['tp']: pf.tp_decorate(A, d['tp'], d['X0'], nel, mark=())
    # 深度箭头：最深（价带顶所在的口袋）、中等、很浅的各一个
    tops.sort()
    pick = [tops[-1], min(tops, key=lambda t: abs(t[0] - u - 0.15)), min((t for t in tops if t[1] < d['X0']), key=lambda t: abs(t[0] - u - 0.022))]
    for e0, x0 in pick:
        A.annotate('', xy=(x0, u), xytext=(x0, e0), arrowprops=dict(arrowstyle='<->', lw=0.8, color='#B36B00', shrinkA=0, shrinkB=0,
                   mutation_scale=6), zorder=6)
        A.text(x0 + 0.03, max((e0 + u) / 2, u + 0.03), f'D = {1000 * (e0 - u):.0f} meV', fontsize=6.5, va='bottom' if e0 - u < 0.05 else 'center',
               ha='left', zorder=7, bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.9))
    A.set_ylabel(r'$E - E_{\mathrm{VBM}}$ (eV)')
    A.set_title(f'(a) {m.replace("2", "₂").replace("5", "₅")}，PBE+SOC 能带', fontsize=8, loc='left')
    fig.text(0.075, -0.075, f'(a) 涂色：µ 以上装着空穴的口袋（沿这条路径的截面），此路径上 {npocket} 段。○：µ 与能带的交点（{ncross} 个）= 费米面，'
             f'输运只在这里发生；交点处能带的斜率 = ħv。\n     箭头 D：µ 在谷顶以下的深度。SOC 下每条带是二重简并的一对，只数一条。'
             f'(b)(c) 为 BoltzTraP2（PBE+SOC 密网格，300 K）的结果，纵轴与 (a) 对齐。',
             fontsize=6.6, ha='left', va='top', linespacing=1.4)

    # ---------- (b) 空穴浓度 p(µ)
    ok = d['C'] > 0
    B.plot(d['C'][ok], d['mu'][ok], color=C_BAND, lw=1.1)
    for N, Q, col in ((NH, H, C_H), (NL, L, C_L)):
        B.axhline(Q['mu'], color=col, lw=1.4 if col == C_H else 1.0, ls=pf.DASH)
        B.plot(abs(N), Q['mu'], 'o', ms=4, mfc=col, mec='white', mew=0.6, zorder=5)
    B.set_xscale('log'); B.set_xlim(3e18, 5e21); B.set_xticks([1e19, 1e20, 1e21])
    B.set_xlabel(r'$p$ (cm$^{-3}$)')
    B.set_title('(b) 空穴数 p(µ)', fontsize=8, loc='left')
    B.text(4e18, u + 0.008, f'µ：霍尔 p = {usci(NH)}', fontsize=6.3, ha='left', va='bottom', bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.9))
    B.text(4.5e21, L['mu'] - 0.008, f'µ：p = {usci(NL)}', fontsize=6.3, ha='right', va='top', bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.9))
    B.text(0.05, 0.97, 'µ 以上到价带顶的\n空态全部计入；\n由霍尔浓度定出 µ', transform=B.transAxes, fontsize=6.3, va='top', linespacing=1.3)

    # ---------- (c) σ/τ(µ) 与切线 → S
    Cx.plot(d['sig'][ok], d['mu'][ok], color=C_BAND, lw=1.1)
    for Q, col in ((H, C_H), (L, C_L)):
        Cx.axhline(Q['mu'], color=col, lw=1.4 if col == C_H else 1.0, ls=pf.DASH)
        ee = np.linspace(Q['mu'] - 0.05, Q['mu'] + 0.05, 20)
        Cx.plot(Q['sig'] * np.exp(Q['dlns'] * (ee - Q['mu'])), ee, color='#C0392B', lw=1.0, zorder=5)
        Smott = -(np.pi ** 2 / 3) * pf.mv.KB * kT * Q['dlns'] * 1e6
        Cx.text(Q['sig'] * 1.25, Q['mu'] + 0.018, f'$S$ = {Smott:.0f}\n(BoltzTraP2 {Q["S"]:.1f})',
                fontsize=6.2, va='bottom', ha='left', linespacing=1.15,
                bbox=dict(fc='white', ec='none', pad=0.3, alpha=0.9))
    Cx.set_xscale('log'); Cx.set_xlim(2e18, 1e21)
    Cx.set_xlabel(r'$\sigma/\tau$ (S m$^{-1}$ s$^{-1}$)')
    Cx.set_title('(c) σ(µ) 与切线', fontsize=8, loc='left')
    Cx.text(0.05, 0.97, '红线：切线。越接近水平\n（σ 随能量变得越快），\nS 越大：\n$S=-\\frac{\\pi^2}{3}\\frac{k_\\mathrm{B}}{e}k_\\mathrm{B}T\\,\\frac{d\\ln\\sigma}{dE}$',
            transform=Cx.transAxes, fontsize=6.3, va='top', linespacing=1.3)
    for ax in (B, Cx): ax.tick_params(axis='y', labelleft=False)
    os.makedirs(a.out, exist_ok=True)
    f = os.path.join(a.out, f'讲解_费米面口袋与S_{m}.png'); fig.savefig(f, dpi=300); plt.close(fig)
    print(f, f'µ_H = {1000 * u:.0f} meV, S_H: Mott {-(np.pi ** 2 / 3) * pf.mv.KB * kT * H["dlns"] * 1e6:.1f} / BT2 {H["S"]:.1f};',
          f'µ_L = {1000 * L["mu"]:.0f} meV, S_L: Mott {-(np.pi ** 2 / 3) * pf.mv.KB * kT * L["dlns"] * 1e6:.1f} / BT2 {L["S"]:.1f};',
          f'pockets {npocket}, crossings {ncross}; tops above µ (meV):', sorted(round(1000 * (t[0] - u)) for t in tops))


if __name__ == '__main__':
    main()
