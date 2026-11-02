"""霍尔数据与电导率（2026-11-02 用户提供）对照计算，结果写在 21_阶段报告与交接/实验数据/霍尔与电导率对照_20261102.md。

1. σ 的单位：表头为 10⁴ S/m；用霍尔测量的电阻率独立核对（Pb₂Bi₂Te₅、Pb₂Sb₂Te₅）。
2. Pb₂Bi₂Te₅（n 型，霍尔 n = 8.11×10¹⁹ cm⁻³）：在实测浓度下各种算法的 S(T)；单抛物带 Pisarenko 有效质量。
3. Pb₂Sb₂Te₅（霍尔 p = 1.04×10²¹ cm⁻³，仪器标红）：在该浓度下的 S；按报告的近补偿解释，载流子数与实测 σ 是否相容。
4. 两个 Sn 样品（无霍尔）：σ(SnSb)/σ(SnBi) 的实验值与三种反推浓度下的预测（同一弛豫时间）；Snyder 加权迁移率。

用法：python 10_脚本/hall_check.py      只读仓库内的数据，不写文件。
"""
import csv, os, sys
import numpy as np
from scipy import integrate, optimize

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import paper_figures as pf
import multivalley_seebeck as mv

DATA = os.path.join(ROOT, '21_阶段报告与交接', '实验数据')
KB, E, H, ME = 1.380649e-23, 1.602176634e-19, 6.62607015e-34, 9.1093837e-31


def read_tss(m):
    """T (K)、σ (S/m)、S (µV/K)。Pb₂Sb₂Te₅ 的旧文件表头写“单位未给出”，单位同为 10⁴ S/m。"""
    rows = list(csv.reader(open(os.path.join(DATA, f'{m}_T_sigma_S.csv'), encoding='utf-8-sig')))[1:]
    return [(float(t), float(s) * 1e4, float(S)) for t, s, S in rows]


def read_hall():
    return {r['sample']: r for r in csv.DictReader(open(os.path.join(DATA, '霍尔_20261102.csv'), encoding='utf-8-sig'))}


def F(j, eta):
    return integrate.quad(lambda x: x ** j / (1 + np.exp(x - eta)), 0, max(60, eta + 60))[0]


def spb_mass(S, T, n_m3, r):
    """单抛物带：由 |S| 定约化费米能 η，再由浓度定态密度有效质量。r = −1/2 声学声子，0 常数弛豫时间。"""
    lam = r + 0.5
    eta = optimize.brentq(lambda x: (KB / E) * ((2 + lam) * F(1 + lam, x) / ((1 + lam) * F(lam, x)) - x) * 1e6 - abs(S), -10, 40)
    return eta, (n_m3 / (4 * np.pi * F(0.5, eta))) ** (2 / 3) * H ** 2 / (2 * KB * T) / ME


def weighted_mobility(S, sigma, T):
    """Snyder 等，Adv. Mater. 32, 2001537 (2020)：由 S 与 σ 求加权迁移率（cm²/V·s）。"""
    s = abs(S) * 1e-6 / (KB / E)
    pre = 3 * H ** 3 * sigma / (8 * np.pi * E * (2 * ME * KB * T) ** 1.5) * 1e4
    return pre * (np.exp(s - 2) / (1 + np.exp(-5 * (s - 1))) + (3 / np.pi ** 2) * s / (1 + np.exp(5 * (s - 1))))


def sigma_over_tau(m, side, N, T, tau):
    """多谷模型的 σ/τ（任意单位，与 mv.solve 同口径），只用于两种材料之比。"""
    V, gap = mv.build(m, 'v', 300, None, side); C, _ = mv.build(m, 'c', 300, side=side)
    kT = mv.KB * T; occ = mv._occ(T); EF, _ = mv.fermi_level(V, C, gap, N, T); den = 0.0
    for sgn, L in ((1, V), (-1, C)):
        for v in L:
            eta = ((-v['off'] - EF) if sgn > 0 else (EF - gap - v['off'])) / kT
            den += v['g'] * np.trapezoid(v['vx2'] * (v['dos'] if tau == 'crta' else 1.0) * occ(eta)[1], mv.EPS)
    return den


def main():
    hall = read_hall(); exp = {m: read_tss(m) for m in pf.MATS}
    print('# 1 σ 的单位（表头 10⁴ S/m）：霍尔测得的 1/ρ 与表中最低温 σ 之比')
    for m in ('Pb2Bi2Te5', 'Pb2Sb2Te5'):
        s_h = 1 / float(hall[m]['rho_Ohm_m']); s_t = exp[m][0][1]
        print(f'  {m}: 霍尔 {s_h:.3g} S/m，表中 {s_t:.3g} S/m（{exp[m][0][0]:.0f} K），比值 {s_t / s_h:.2f}')

    m = 'Pb2Bi2Te5'; n = float(hall[m]['n_cm-3']); N = -n; data = [(t, S) for t, _, S in exp[m]]
    print(f'\n# 2 {m}：霍尔 n = {n:.3g} cm⁻³，μ_H = {hall[m]["mu_H_cm2/Vs"]} cm²/V·s')
    b = pf.bt2_signed(m, -1); gh, gp, ga = pf.model(m, 'c2_hse'), pf.model(m, 'c2_pbe'), pf.model(m, 'c2_hse', 'adp')
    meth = {'PBE+SOC BoltzTraP2（全区，常数τ）': lambda N, t: b(N, t),
            'HSE = BoltzTraP2 × HSE/PBE（常数τ）': lambda N, t: b(N, t) * gh(N, t) / gp(N, t),
            'HSE = BoltzTraP2 × HSE(声学声子)/PBE': lambda N, t: b(N, t) * ga(N, t) / gp(N, t)}
    print('  T/K    实验  ' + '  '.join(meth))
    for t, S in data:
        print(f'  {t:5.0f} {S:7.1f} ' + ' '.join(f'{f(N, t):9.1f}' for f in meth.values()))
    for k, f in meth.items():
        try: Nf = pf.solve_log(f, data, 17.5, 20.5, -1); print(f'  {k}：由 302 K 的 S 反推 n = {abs(Nf):.2e} cm⁻³')
        except ValueError: print(f'  {k}：10^17.5–10^20.5 内没有解')
    for r, nm in ((-0.5, '声学声子'), (0.0, '常数弛豫时间')):
        eta, ms = spb_mass(data[0][1], data[0][0], n * 1e6, r)
        print(f'  单抛物带（{nm}）：η = {eta:.2f}，m*_DOS = {ms:.2f} m_e')

    m = 'Pb2Sb2Te5'; p = float(hall[m]['n_cm-3'])
    print(f'\n# 3 {m}：霍尔 p = {p:.3g} cm⁻³（仪器标红，相位 {hall[m]["phase_deg"]}°）')
    b = pf.bt2_signed(m, 1); gh, gp = pf.model(m, 'c2_hse'), pf.model(m, 'c2_pbe')
    for t, _, S in exp[m][::7]:
        print(f'  {t:.0f} K：实验 {S:+.1f}；PBE BoltzTraP2 {b(p, t):+.1f}；HSE = BoltzTraP2 × 比值 {b(p, t) * gh(p, t) / gp(p, t):+.1f} µV/K')
    s302 = exp[m][0][1]
    for mu in (100, 1000):
        print(f'  σ(302 K) = {s302:.3g} S/m 若迁移率为 {mu} cm²/V·s，需载流子 {s302 / (E * mu * 1e-4) * 1e-6:.1e} cm⁻³'
              f'（报告的近补偿解释：净掺杂约 2×10¹⁶，本征约 1×10¹⁷ cm⁻³）')

    print('\n# 4 两个 Sn 样品：σ(SnSb)/σ(SnBi)')
    print('  实验：' + ' '.join(f'{a[1] / b_[1]:.2f}' for a, b_ in zip(exp['Sn2Sb2Te5'], exp['Sn2Bi2Te5'])) + '（300 → 620 K）')
    for nm, side, tau, pa, pb in (('HSE 常数τ', 'c2_hse', 'crta', 6.157e20, 7.221e20), ('HSE 声学声子', 'c2_hse', 'adp', 5.411e20, 3.534e20),
                                  ('PBE 常数τ', 'c2_pbe', 'crta', 2.308e21, 1.242e21)):
        r = sigma_over_tau('Sn2Sb2Te5', side, pa, 300, tau) / sigma_over_tau('Sn2Bi2Te5', side, pb, 300, tau)
        print(f'  {nm}：反推浓度 SnSb {pa:.1e}、SnBi {pb:.1e} cm⁻³ → 预测 {r:.2f}（同一 τ）')
    print('\n# 加权迁移率 μ_w（cm²/V·s，300 K 附近 → 620 K 附近）')
    for m in pf.MATS:
        v = [weighted_mobility(S, s, t) for t, s, S in exp[m]]
        print(f'  {m}: {v[0]:.0f} → {v[-1]:.0f}' + ('（近补偿或补偿时不适用）' if m == 'Pb2Sb2Te5' else ''))


if __name__ == '__main__':
    main()
