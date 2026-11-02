"""论文报告第 4、5 节引用的输运数值，一次算齐，供逐条核对正文（2026-11-02，取谷范围 mv.WIN = 700 meV）。

用法：python 10_脚本/report_numbers_check.py      只读，不写文件。
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import paper_figures as pf
import multivalley_seebeck as mv
import report_figures_20261028 as rf

T = 300


def main():
    meth = {(m, s): pf.s_methods(m, s) for m in pf.MATS for s in (1, -1)}

    print('# 表 6：300 K 面内 S（µV/K），PBE BoltzTraP2 → HSE 估计')
    for m in pf.MATS:
        cells = []
        for s, N in ((1, 3e19), (1, 1e20), (-1, 3e19), (-1, 1e20)):
            f = meth[(m, s)]; cells.append(f'{f["P"](s * N, T):7.1f} → {f["H"](s * N, T):6.1f}')
        print(f'  {m}: p3e19 {cells[0]} | p1e20 {cells[1]} | n3e19 {cells[2]} | n1e20 {cells[3]}')
    est = rf.hse_estimates()
    print('  rf.hse_estimates（表 6 第 1、3 列的原始来源，BTP 为 2.4 节的 BoltzTraP2 值）:', {m: [round(x, 1) for x in v] for m, v in est.items()})
    print('  BTP:', rf.BTP)

    print('\n# Sn₂Bi₂Te₅ / Sn₂Sb₂Te₅ 的 S 之比（300 K，同一空穴浓度）')
    for N in (1e19, 3e19, 1e20, 3e20, 6e20, 1e21, 1.72e21):
        h = [meth[(m, 1)]['H'](N, T) for m in ('Sn2Bi2Te5', 'Sn2Sb2Te5')]; p = [meth[(m, 1)]['P'](N, T) for m in ('Sn2Bi2Te5', 'Sn2Sb2Te5')]
        print(f'  {N:.2e}: HSE {h[0]:6.1f}/{h[1]:6.1f} = {h[0] / h[1]:.2f}   PBE {p[0]:6.1f}/{p[1]:6.1f} = {p[0] / p[1]:.2f}')

    print('\n# 同一空穴浓度下四种材料（300 K）')
    for N in (3e19, 1e21):
        for k in ('H', 'P'):
            print(f'  {N:.0e} {pf.SNAME[k]:10s}: ' + '  '.join(f'{m} {meth[(m, 1)][k](N, T):6.1f}' for m in pf.MATS))

    print('\n# n 型：Pb₂Bi₂Te₅ 霍尔浓度下四种材料的 |S|（300 K）')
    Nn = pf.HALL['Pb2Bi2Te5']['N']
    for k in ('H', 'P'):
        print(f'  {abs(Nn):.2e} {pf.SNAME[k]:10s}: ' + '  '.join(f'{m} {-meth[(m, -1)][k](Nn, T):6.1f}' for m in pf.MATS))

    print('\n# 样品：霍尔浓度下的 S、偏差、反推浓度（同 Fig9_S_vs_T/fit_summary.csv）')
    for m in pf.MATS:
        h = pf.HALL[m]; N = h['N']; s = 1 if N > 0 else -1; d = pf.EXP[m]
        lo, hi = (17.5, 20.5) if s < 0 else (19, 21.8)
        print(f'  {m}: 霍尔 {N:+.3e}{"（标红）" if h["flagged"] else ""}，实验 {d[0][0]:.0f} K {d[0][1]:+.1f}、{d[-1][0]:.0f} K {d[-1][1]:+.1f}')
        for k in ('H', 'A', 'P'):
            f = meth[(m, s)][k]; dev = [f(N, t) - S for t, S in d]
            try: Nf = abs(pf.solve_log(f, d, lo, hi, s)); inv = f'反推 {Nf:.2e}（{Nf / abs(N):.2f} 倍）'
            except ValueError: inv = '反推无解'
            print(f'     {pf.SNAME[k]:28s} {f(N, d[0][0]):+7.1f} … {f(N, d[-1][0]):+7.1f}  偏差 {min(dev):+6.1f}…{max(dev):+6.1f}  '
                  f'均方根 {np.sqrt(np.mean(np.square(dev))):5.1f}  {inv}')

    print('\n# Sn 两个样品之比（约 300 K）')
    ex = pf.EXP['Sn2Bi2Te5'][0][1] / pf.EXP['Sn2Sb2Te5'][0][1]
    nb, ns = pf.HALL['Sn2Bi2Te5']['N'], pf.HALL['Sn2Sb2Te5']['N']
    print(f'  实验 {ex:.2f}；实验 S 范围比 ' + ', '.join(f'{a[1] / b[1]:.2f}' for a, b in zip(pf.EXP['Sn2Bi2Te5'], pf.EXP['Sn2Sb2Te5'])))
    for k in ('H', 'P'):
        fb, fs = meth[('Sn2Bi2Te5', 1)][k], meth[('Sn2Sb2Te5', 1)][k]
        print(f'  {pf.SNAME[k]}: 同浓度（SnBi 的霍尔值）{fb(nb, 301) / fs(nb, 301):.2f}；各取霍尔值 {fb(nb, 301) / fs(ns, 301):.2f}')

    print('\n# 三维搜索（真实谷）对 3×10¹⁹ cm⁻³ p 型 S 的影响（HSE 估计）')
    for m in ('Sn2Sb2Te5', 'Sn2Bi2Te5'):
        vals = []
        for use in (False, True):
            mv.USE_TRUE = use
            V, gap = mv.build(m, 'v', mv.WIN, mv.c3_usable(m)[0], 'c2_hse'); C, _ = mv.build(m, 'c', mv.WIN)   # 与 rf.hse_estimates 相同
            Vp, gp_ = mv.build(m, 'v', mv.WIN, None, 'c2_pbe'); Cp, _ = mv.build(m, 'c', mv.WIN, side='c2_pbe')
            vals.append(rf.BTP[m][0] * mv.solve(V, C, gap, 3e19, T)[0] / mv.solve(Vp, Cp, gp_, 3e19, T)[0])
        mv.USE_TRUE = True
        print(f'  {m}: 只用路径 {vals[0]:.0f} → 真实谷 {vals[1]:.0f}')


if __name__ == '__main__':
    main()
