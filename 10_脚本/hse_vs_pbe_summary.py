#!/usr/bin/env python3
"""HSE06 与 PBE 对比汇总（2026-11-05，用户问：HSE 是否比 PBE 更明显、更支持孤对电子的结论）。

不做新的计算，也不估计：数值全部取自主报告的源数据（论文图_20261028/source_data：表 6、7、8，图 7，图 3-3），
只有 (d) 的 S(Sn₂Bi₂Te₅)/S(Sn₂Sb₂Te₅) 随浓度的曲线直接由两种泛函的 BoltzTraP2 结果算出（与表 6 同一取值方法）。

输出 21_阶段报告与交接/HSE与PBE对比_20261105/：
  HSE_vs_PBE_summary.png / .pdf      六个子图，每个子图回答一个问题
  data/*.csv                          每个子图一份源数据，另有 verdict.csv（逐条结论：PBE 值、HSE 值、哪个更支持、实验能否检验）
  HSE_vs_PBE_作图数据.xlsx            同样的数据，第一页为说明与结论表（xlsxwriter 写出并读回检查）
用法：python 10_脚本/hse_vs_pbe_summary.py
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import transport_figures as tf
import plot_data_index as pdi
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

SD = os.path.join(pf.ROOT, '21_阶段报告与交接', '论文图_20261028', 'source_data')
OUT = os.path.join(pf.ROOT, '21_阶段报告与交接', 'HSE与PBE对比_20261105')
MATS, TEX = pf.MATS, pf.TEX
SH = {'Pb2Sb2Te5': 'PbSb', 'Pb2Bi2Te5': 'PbBi', 'Sn2Sb2Te5': 'SnSb', 'Sn2Bi2Te5': 'SnBi'}
CH, CP = pf.C_HSE, pf.C_PBE
rd = lambda rel: list(csv.DictReader(open(os.path.join(SD, rel), encoding='utf-8-sig')))
f1 = lambda v: round(float(v), 1)
SHEETS = {}                                   # 名称 → (说明, 表头, 行)


def put(name, note, header, rows):
    SHEETS[name] = (note, header, rows)
    os.makedirs(os.path.join(OUT, 'data'), exist_ok=True)
    with open(os.path.join(OUT, 'data', name + '.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def bars(ax, labels, ph, pp, ylabel, fmt=None):
    x = np.arange(len(labels)); wd = 0.36
    ax.bar(x - wd / 2 - 0.01, pp, wd, color=CP, label='PBE+SOC', zorder=3)
    ax.bar(x + wd / 2 + 0.01, ph, wd, color=CH, label='HSE06+SOC', zorder=3)
    if fmt:
        for xi, v in zip(x - wd / 2 - 0.01, pp): ax.text(xi, v, fmt(v), ha='center', va='bottom' if v >= 0 else 'top', fontsize=5.6, color='#333333')
        for xi, v in zip(x + wd / 2 + 0.01, ph): ax.text(xi, v, fmt(v), ha='center', va='bottom' if v >= 0 else 'top', fontsize=5.6, color='#333333')
    ax.set_xticks(x, labels, fontsize=6.5); ax.tick_params(axis='x', length=0, top=False); ax.set_ylabel(ylabel)
    ax.axhline(0, color='k', lw=0.6, zorder=2); ax.grid(axis='y', color=pf.GRID, lw=0.4, zorder=0)


def main():
    pf.style(); os.makedirs(OUT, exist_ok=True)
    t6 = {r['material']: r for r in rd('Table6_S_300K/S_300K_HSE_and_PBE.csv')}
    t7 = {r['material']: r for r in rd('Table7_S_T_at_Hall/S_T_at_Hall_HSE_and_PBE.csv')}
    t8 = {r['material']: r for r in rd('Table8_Pisarenko_mass/samples_Pisarenko_mass_HSE_and_PBE.csv')}
    f7 = {(r['material'], r['functional'][:3]): r for r in rd('Fig7_lone_pair_valley/highest_lone_pair_valley.csv')}
    co = {(r['material'], r['site']): r for r in rd('Fig3-3_COHP_HSE_vs_PBE/ICOHP_HSE_vs_PBE.csv')}
    lab = [SH[m] for m in MATS]                    # 横轴用简称（PbSb = Pb₂Sb₂Te₅ …），全称写在图注与数据里
    fig, A = plt.subplots(2, 3, figsize=(pf.W2, 5.3), gridspec_kw=dict(hspace=0.62, wspace=0.42))

    # (a) 最高的孤对电子谷在不在价带顶
    pp = [float(f7[(m, 'PBE')]['E-EVBM_meV']) for m in MATS]; ph = [float(f7[(m, 'HSE')]['E-EVBM_meV']) for m in MATS]
    ax = A[0, 0]; ax.axhspan(-2 * pf.KT300, 0, color='#EEF3F8', zorder=0, lw=0)
    bars(ax, lab, ph, pp, r'$E - E_\mathrm{VBM}$ (meV)', lambda v: f'{v:.0f}'.replace('-', '−'))
    ax.text(3.45, -2 * pf.KT300 + 4, r'$2k_\mathrm{B}T$', ha='right', va='bottom', fontsize=6, color='#555555')
    ax.set_ylim(-265, 25); ax.set_title('(a) Is the highest lone-pair valley\nat the VBM?', fontsize=7.4, pad=3)
    put('a_lone_pair_valley', '最高的孤对电子价带谷相对价带顶的能量（meV），即图 7 的数值；阴影为 2k_BT（300 K）。',
        ['material', 'E_PBE+SOC_meV', 'E_HSE06+SOC_meV'], [[m, f1(a), f1(b)] for m, a, b in zip(MATS, pp, ph)])

    # (b) 低浓度 p 型 S
    pp = [float(t6[m]['PBE_S_p_3e19_uV/K']) for m in MATS]; ph = [float(t6[m]['HSE_S_p_3e19_uV/K']) for m in MATS]
    ax = A[0, 1]; bars(ax, lab, ph, pp, r'$S$ (μV K$^{-1}$)', lambda v: f'{v:.0f}')
    ax.set_ylim(0, 260); ax.set_title('(b) Do the four materials differ at\n' + r'$p$ = 3$\times$10$^{19}$ cm$^{-3}$ (300 K)?', fontsize=7.4, pad=3)
    put('b_S_p_3e19', '300 K、空穴 3×10¹⁹ cm⁻³ 的面内 S（µV/K），即表 6。',
        ['material', 'S_PBE+SOC_uV/K', 'S_HSE06+SOC_uV/K'], [[m, a, b] for m, a, b in zip(MATS, pp, ph)])

    # (c) A 位：同一 B 位下 Sn 比 Pb 高多少
    rows, cats, ph, pp = [], [], [], []
    for B in ('Sb', 'Bi'):
        for N, key in (('3e19', '3×10$^{19}$'), ('1e20', '10$^{20}$')):
            sn, pb = f'Sn2{B}2Te5', f'Pb2{B}2Te5'
            dh = float(t6[sn][f'HSE_S_p_{N}_uV/K']) - float(t6[pb][f'HSE_S_p_{N}_uV/K'])
            dp = float(t6[sn][f'PBE_S_p_{N}_uV/K']) - float(t6[pb][f'PBE_S_p_{N}_uV/K'])
            cats.append(f'{B} series\n{key}'); ph.append(dh); pp.append(dp); rows.append([f'{sn} − {pb}', N, round(dp, 1), round(dh, 1)])
    ax = A[0, 2]; bars(ax, cats, ph, pp, r'$S$(Sn) $-$ $S$(Pb) (μV K$^{-1}$)', lambda v: '0' if abs(v) < 0.5 else f'{v:+.0f}'.replace('-', '−'))
    ax.set_ylim(-25, 80); ax.set_title('(c) Does Sn (A site) raise S over Pb\nat the same hole concentration?', fontsize=7.4, pad=3)
    put('c_A_site_dS', '同一 B 位、同一空穴浓度、300 K 下 S(Sn 化合物) − S(Pb 化合物)（µV/K），由表 6 相减。',
        ['pair', 'hole_conc_cm-3', 'dS_PBE+SOC_uV/K', 'dS_HSE06+SOC_uV/K'], rows)

    # (d) B 位：同一浓度下 S(Sn₂Bi₂Te₅)/S(Sn₂Sb₂Te₅)，直接由 BoltzTraP2 算
    NS = np.logspace(19, np.log10(2e21), 41); cur = {}
    for fn in ('PBE', 'HSE'):
        Tb, fb = tf.LOAD[fn]('Sn2Bi2Te5'); Ts, fs = tf.LOAD[fn]('Sn2Sb2Te5'); ib, is_ = list(Tb).index(300), list(Ts).index(300)
        cur[fn] = np.array([fb(N, ib) / fs(N, is_) for N in NS])
        cur[fn + '_fixed'] = {N: fb(N, ib) / fs(N, is_) for N in (3e19, 1e20, 3e20, 1e21)}
    def s_low(m, fn):                                            # 霍尔浓度、最低测量温度下的计算 S（不取表 7 的四舍五入值）
        T, v, _, _ = tf.st_curve(m, fn); return float(np.interp(pf.EXP[m][0][0], T, v))
    hall = {fn: s_low('Sn2Bi2Te5', fn) / s_low('Sn2Sb2Te5', fn) for fn in ('HSE', 'PBE')}
    hexp = pf.EXP['Sn2Bi2Te5'][0][1] / pf.EXP['Sn2Sb2Te5'][0][1]
    ax = A[1, 0]
    ax.plot(NS, cur['PBE'], color=CP, lw=1.2, ls=pf.DASH); ax.plot(NS, cur['HSE'], color=CH, lw=1.4)
    ax.axhline(1, color='k', lw=0.6); ax.set_xscale('log'); ax.set_xlim(NS[0], NS[-1]); ax.set_ylim(0.80, 1.50)
    ax.set_xlabel(r'Hole concentration (cm$^{-3}$)'); ax.set_ylabel(r'$S$(Sn$_2$Bi$_2$Te$_5$) / $S$(Sn$_2$Sb$_2$Te$_5$)')
    ax.grid(color=pf.GRID, lw=0.4, zorder=0)
    ax.text(0.97, 0.04, f'at the two Hall conc. (~300 K):\nexp {hexp:.2f}, HSE {hall["HSE"]:.2f}, PBE {hall["PBE"]:.2f}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=5.8, color='#333333', bbox=dict(fc='white', ec='none', pad=0.5, alpha=0.9))
    ax.set_title('(d) Does Bi (B site) raise S over Sb\nat the same concentration (300 K)?', fontsize=7.4, pad=3)
    put('d_B_site_ratio_vs_conc', 'S(Sn₂Bi₂Te₅)/S(Sn₂Sb₂Te₅)，同一空穴浓度、300 K，由两种泛函的 BoltzTraP2 结果直接算出（浓度方向按 log|N| 线性取值，同表 6）。'
        '另：两个样品在各自的霍尔浓度、约 300 K 下的比值（实验、HSE、PBE，取自表 7）写在 d_B_site_ratio_at_Hall。',
        ['hole_conc_cm-3', 'ratio_PBE+SOC', 'ratio_HSE06+SOC'], [[f'{N:.6e}', round(a, 4), round(b, 4)] for N, a, b in zip(NS, cur['PBE'], cur['HSE'])])
    put('d_B_site_ratio_fixed_conc', '同上，取 3×10¹⁹、1×10²⁰、3×10²⁰、1×10²¹ cm⁻³ 四个浓度（即报告 4.2 节的表）。',
        ['hole_conc_cm-3', 'ratio_PBE+SOC', 'ratio_HSE06+SOC'], [[f'{N:.0e}', round(cur['PBE_fixed'][N], 3), round(cur['HSE_fixed'][N], 3)] for N in (3e19, 1e20, 3e20, 1e21)])
    put('d_B_site_ratio_at_Hall', '两个 Sn 样品各自在霍尔浓度下的 S 之比（约 300 K，表 7）：不同浓度的比较，含浓度的贡献。',
        ['quantity', 'experiment', 'PBE+SOC', 'HSE06+SOC'], [['S(Sn2Bi2Te5)/S(Sn2Sb2Te5) at Hall', round(hexp, 3), round(hall['PBE'], 3), round(hall['HSE'], 3)]])

    # (e) COHP：价带顶以下 1 eV 内的 s–p 反键
    keys = [(m, 'A') for m in MATS] + [(m, 'B') for m in MATS]
    cats = [f"{SH[m]}\n{co[(m, s)]['cation']}" for m, s in keys]
    pp = [float(co[k]['ICOHP_s-p_within_1eV_below_VBM_PBE_eV']) for k in keys]; ph = [float(co[k]['ICOHP_s-p_within_1eV_below_VBM_HSE_eV']) for k in keys]
    ax = A[1, 1]; bars(ax, cats, ph, pp, 's–p ICOHP, top 1 eV of VB\n(eV/bond; > 0 antibonding)')
    for t in ax.get_xticklabels(): t.set_fontsize(5.4)
    ax.axvline(3.5, color='#999999', lw=0.6); ax.text(1.5, 1.33, 'A site–Te', ha='center', fontsize=6.2); ax.text(5.5, 1.33, 'B site–Te', ha='center', fontsize=6.2)
    ax.set_ylim(0, 1.45); ax.set_title('(e) Is the filled s–p antibonding near the\nVBM stronger in HSE? (COHP, no SOC)', fontsize=7.4, pad=3)
    put('e_COHP_s-p_near_VBM', '每种材料最短的阳离子–Te 键，价带顶以下 1 eV 内的阳离子 s–Te p ICOHP（eV/键，正为反键；COGITO，不加 SOC，各以本泛函的价带顶为零点），即图 3-3 的数据；'
        '另列全部占据态的 s–p、p–p 与总 ICOHP。',
        ['material', 'site', 'cation', 'bond_lengths_A', 's-p_top1eV_PBE', 's-p_top1eV_HSE', 's-p_all_PBE', 's-p_all_HSE', 'p-p_PBE', 'p-p_HSE', 'total_PBE', 'total_HSE'],
        [[m, s, co[(m, s)]['cation'], co[(m, s)]['bond_lengths_averaged_A'], co[(m, s)]['ICOHP_s-p_within_1eV_below_VBM_PBE_eV'], co[(m, s)]['ICOHP_s-p_within_1eV_below_VBM_HSE_eV'],
          co[(m, s)]['ICOHP_s-p_PBE_eV'], co[(m, s)]['ICOHP_s-p_HSE_eV'], co[(m, s)]['ICOHP_p-p_PBE_eV'], co[(m, s)]['ICOHP_p-p_HSE_eV'],
          co[(m, s)]['ICOHP_total_PBE_eV'], co[(m, s)]['ICOHP_total_HSE_eV']] for m, s in keys])

    # (f) 与样品：300–620 K 的均方根偏差
    pp = [float(t7[m]['rms_PBE_uV/K']) for m in MATS]; ph = [float(t7[m]['rms_HSE_uV/K']) for m in MATS]
    cats = [f"{SH[m]}\n({t7[m]['type']}{', flagged' if t7[m]['Hall_flagged'] == 'yes' else ''})" for m in MATS]
    ax = A[1, 2]; bars(ax, cats, ph, pp, r'rms deviation, 300–620 K (μV K$^{-1}$)', lambda v: f'{v:.1f}')
    for t in ax.get_xticklabels(): t.set_fontsize(5.6)
    ax.set_ylim(0, 58); ax.set_title('(f) Which functional matches the samples?\n(at the Hall concentration, no fitting)', fontsize=7.4, pad=3)
    put('f_rms_vs_experiment', '霍尔浓度下计算的 S(T) 与实验的均方根偏差（300–620 K，8 个测量温度），即表 7；另列约 300 K 的实验与计算值，以及表 8 的态密度有效质量。'
        'Pb₂Sb₂Te₅ 为两相共存的推测样品（第 5.3 节），不作单相比较。',
        ['material', 'type', 'Hall_flagged', 'rms_PBE_uV/K', 'rms_HSE_uV/K', 'S_exp_300K', 'S_PBE_300K', 'S_HSE_300K', 'm_exp_me', 'm_PBE_me', 'm_HSE_me'],
        [[m, t7[m]['type'], t7[m]['Hall_flagged'], t7[m]['rms_PBE_uV/K'], t7[m]['rms_HSE_uV/K'], t7[m]['S_exp_low_uV/K'], t7[m]['S_PBE_low_uV/K'], t7[m]['S_HSE_low_uV/K'],
          t8[m]['m_DOS_exp_SPB_me'], t8[m]['m_DOS_PBE_SPB_me'], t8[m]['m_DOS_HSE_SPB_me']] for m in MATS])

    fig.legend([Patch(color=CP), Patch(color=CH), Line2D([], [], color=CP, ls=pf.DASH, lw=1.2), Line2D([], [], color=CH, lw=1.4)],
               ['PBE+SOC (bars)', 'HSE06+SOC (bars)', 'PBE+SOC (line, d)', 'HSE06+SOC (line, d)'], loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.955), fontsize=6.6)
    fig.text(0.5, 0.005, 'PbSb = Pb$_2$Sb$_2$Te$_5$, PbBi = Pb$_2$Bi$_2$Te$_5$, SnSb = Sn$_2$Sb$_2$Te$_5$, SnBi = Sn$_2$Bi$_2$Te$_5$.  '
             'Bars and lines are taken from the report (Tables 6–8, Figs. 3-3, 7); (d) is computed from the BoltzTraP2 results.', ha='center', va='bottom', fontsize=5.8, color='#444444')
    pf.save(fig, OUT, 'HSE_vs_PBE_summary')

    # 结论表
    s = lambda m, k: float(t6[m][k])
    spread = lambda fn: (min(s(m, f'{fn}_S_p_3e19_uV/K') for m in MATS), max(s(m, f'{fn}_S_p_3e19_uV/K') for m in MATS))
    V = [
        ['孤对电子谷是否在价带顶', '最高孤对电子谷 − 价带顶（meV）', 'PbSb −1.7、PbBi −16.8、SnSb 0、SnBi 0：四种都在 17 meV 内', 'PbSb −236、PbBi −91、SnSb −84、SnBi 0：只有 SnBi 在价带顶',
         'HSE 能区分材料（PBE 下四种都是多谷）', '无直接实验'],
        ['低浓度 S 随谷结构分开', '300 K、3×10¹⁹ cm⁻³ 的 p 型 S 范围（µV/K）', '%.0f–%.0f（差 %.0f）' % (*spread('PBE'), spread('PBE')[1] - spread('PBE')[0]),
         '%.0f–%.0f（差 %.0f）' % (*spread('HSE'), spread('HSE')[1] - spread('HSE')[0]), 'HSE 明显得多', '样品都在约 10²¹ cm⁻³，检验不了；需要 3×10¹⁹ 的 p 型样品'],
        ['A 位：Sn 比 Pb 提高 S', 'S(Sn) − S(Pb)，3×10¹⁹（µV/K）：Sb 系、Bi 系', '%+.0f、%+.0f（Pb 反而略高）' % (rows[0][2], rows[2][2]),
         '%+.0f、%+.0f' % (rows[0][3], rows[2][3]), '只有 HSE 支持', '同上，检验不了'],
        ['B 位：Bi 比 Sb 提高 S（同浓度）', 'S(Sn₂Bi₂Te₅)/S(Sn₂Sb₂Te₅)，3×10¹⁹', '%.2f' % cur['PBE_fixed'][3e19], '%.2f' % cur['HSE_fixed'][3e19],
         'HSE 更明显', '霍尔浓度下的比值：实验 %.2f、HSE %.2f、PBE %.2f，HSE 更接近' % (hexp, hall['HSE'], hall['PBE'])],
        ['孤对电子 = 填满的 s–p 反键（COHP）', '价带顶以下 1 eV 内 s–p ICOHP（eV/键）：A 位、B 位范围',
         '%s；%s' % (rng(co, 'A', 'PBE'), rng(co, 'B', 'PBE')), '%s；%s' % (rng(co, 'A', 'HSE'), rng(co, 'B', 'HSE')),
         '两者相同：都是填满的反键；A 位 HSE 略强 %.2f–%.2f' % dA(co), 'COHP 无直接实验；另一个程序（LOBSTER）给出同样的反键与趋势'],
        ['p 型样品（两个 Sn）', '300–620 K 均方根偏差（µV/K）：SnSb、SnBi', '%s、%s' % (t7['Sn2Sb2Te5']['rms_PBE_uV/K'], t7['Sn2Bi2Te5']['rms_PBE_uV/K']),
         '%s、%s' % (t7['Sn2Sb2Te5']['rms_HSE_uV/K'], t7['Sn2Bi2Te5']['rms_HSE_uV/K']), 'HSE 更接近', 'm*：实验 1.13、1.35；HSE 1.31、1.59；PBE 1.51、1.73'],
        ['n 型样品（Pb₂Bi₂Te₅）：导带多谷', '均方根偏差（µV/K）；约 300 K 的 S', '%s；%s（实验 %s）' % (t7['Pb2Bi2Te5']['rms_PBE_uV/K'], t7['Pb2Bi2Te5']['S_PBE_low_uV/K'], t7['Pb2Bi2Te5']['S_exp_low_uV/K']),
         '%s；%s' % (t7['Pb2Bi2Te5']['rms_HSE_uV/K'], t7['Pb2Bi2Te5']['S_HSE_low_uV/K']), 'PBE 对，HSE 不对（HSE 导带只剩 Γ 单谷）', 'm*：实验 0.84、PBE 0.83、HSE 0.39'],
    ]
    V = [[str(x).replace('-', '−') if j in (2, 3) else x for j, x in enumerate(row)] for row in V]
    put('verdict', '逐条结论：同一指标在 PBE+SOC 与 HSE06+SOC 下的值、哪个更支持报告的说法、实验能否检验。',
        ['claim', 'metric', 'PBE+SOC', 'HSE06+SOC', 'which supports the claim more', 'experimental check'], V)
    write_xlsx()
    print('→', OUT)


def rng(co, site, fn):
    v = [float(r[f'ICOHP_s-p_within_1eV_below_VBM_{fn}_eV']) for (m, s), r in co.items() if s == site]
    return f'+{min(v):.2f} 到 +{max(v):.2f}'


def dA(co):
    d = [float(r['ICOHP_s-p_within_1eV_below_VBM_HSE_eV']) - float(r['ICOHP_s-p_within_1eV_below_VBM_PBE_eV']) for (m, s), r in co.items() if s == 'A']
    return min(d), max(d)


def write_xlsx():
    import xlsxwriter
    path = os.path.join(OUT, 'HSE_vs_PBE_作图数据.xlsx')
    wb = xlsxwriter.Workbook(path, {'strings_to_formulas': False, 'strings_to_urls': False, 'nan_inf_to_errors': False})
    b = wb.add_format({'bold': True}); t = wb.add_format({'bold': True, 'font_size': 12}); wr = wb.add_format({'text_wrap': True, 'valign': 'top'})
    ix = wb.add_worksheet('说明'); r = 0
    ix.write_string(r, 0, 'HSE06+SOC 与 PBE+SOC 对比：哪些结论在 HSE 下更明显（图 HSE_vs_PBE_summary）', t); r += 1
    ix.write_string(r, 0, '数值全部取自主报告源数据（表 6、7、8，图 7，图 3-3）；(d) 的曲线由两种泛函的 BoltzTraP2 结果直接算出。程序：10_脚本/hse_vs_pbe_summary.py', wr); r += 2
    ix.write_string(r, 0, '工作表', b); ix.write_string(r, 1, '对应子图与内容', b); r += 1
    for name, (note, _, _) in SHEETS.items():
        ix.write_string(r, 0, name); ix.write_string(r, 1, note, wr); r += 1
    ix.set_column(0, 0, 26); ix.set_column(1, 1, 120)
    for name, (note, header, rows) in SHEETS.items():
        ws = wb.add_worksheet(name[:31]); ws.write_string(0, 0, note)
        for j, h in enumerate(header): ws.write_string(1, j, h, b)
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                x = pdi.number(str(v))
                (ws.write_number if isinstance(x, float) else ws.write_string)(i + 2, j, x)
        ws.freeze_panes(2, 0); ws.set_column(0, len(header), 18 if name != 'verdict' else 34)
    wb.close(); pdi.check_book(path)


if __name__ == '__main__':
    main()
