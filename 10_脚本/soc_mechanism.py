"""SOC 机制检验（22_能带基准_1-0/SOC机制检验/）：Bi 为什么把 Γ 处的价带态压下去。

① d1_nosoc_<材料>：PBE、不加 SOC，同一条 c2 路径；与现有的 PBE+SOC（<材料>/c2_pbe）对比。
② d2_<成分>_in_<晶格>_lattice：PBE+SOC，把 B 位 Sb ↔ Bi 换进另一种材料的晶格和坐标，分开化学与体积效应。

Γ 点的态按轨道成分追踪，不按带序号（能量都相对 Σ 谷，即 Γ–M 段内部最高价带的极大，meV）：
- s 态：阳离子 s + Te p_z 的反键态（孤对电子在 Γ 处的态）。无 SOC 时是 Γ 价带顶，有 SOC 时被反转到导带底。
- 阳离子 p 态：无 SOC 时在导带，分成 p_z 单态与 p_xy 二重态（两者以 B 位 p 为主），取 (p_z + 2 p_xy)/3 作为 p 能级中心；
  有 SOC 时其中 j = 1/2 型的组合（B 位 p ≥ 25%，p_z : p_xy ≈ 1 : 2）被压到 s 态之下，成为 Γ 价带顶。
  这就是 Γ 处的带反转；Γ 价带顶在有、无 SOC 时是两个不同的态，不能直接相减。
Sb → Bi 的变化拆成：标量部分（无 SOC 时 p 能级中心的变化）+ SOC 部分（SOC 把 p 态压低的量的变化）。

用法：python 10_脚本/soc_mechanism.py [--fig <输出目录>]      （缺少的计算自动跳过）
"""
import argparse, csv, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import c2_valleys as cv

R = os.path.join(ROOT, '22_能带基准_1-0'); D = os.path.join(R, 'SOC机制检验')
NELECT = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
SIGMA_SEG = range(1, 33)                         # Γ–M 段内部（路径序号，M 在 33）
B, TE, A = slice(0, 2), slice(2, 7), slice(7, 9)  # POSCAR 元素顺序：B 位、Te、A 位


def soc_of(d):
    return any(l.split('=')[0].strip() == 'LSORBIT' and 'T' in l.split('=')[1].upper()
               for l in open(os.path.join(d, 'INCAR')) if '=' in l and not l.strip().startswith('#'))


def gamma_states(d, comp):
    f = os.path.join(d, 'PROCAR_OPT.gz')
    if not os.path.exists(f): f = f[:-3]
    if not os.path.exists(f): return None
    E, W = cv.read_procar(f); soc = soc_of(d)
    nv = NELECT[comp] if soc else NELECT[comp] // 2            # 价带数（SOC 下每条带 1 个电子）
    st = 2 if soc else 1                                       # SOC 下 Kramers 对只取一条
    esig = E[:, nv - 1][list(SIGMA_SEG)].max()
    w = W[0]; tot = w.sum(axis=(1, 2))
    fr = lambda i, sl, o: w[i, sl, o].sum() / tot[i]
    rows = []
    for i in range(nv - 1 - 6 * st, nv + 6 * st, st):
        rows.append(dict(band=i, vb=i < nv, E=1000 * (E[0, i] - esig), cs=fr(i, A, [0]) + fr(i, B, [0]),
                         bp=fr(i, B, [1, 2, 3]), bpz=fr(i, B, [2]), bpxy=fr(i, B, [1, 3]), tep=fr(i, TE, [1, 2, 3]), tepz=fr(i, TE, [2])))
    top = [r for r in rows if r['vb']][-1]; cb = next(r for r in rows if not r['vb'])
    near = [top, cb, next(r for r in rows if r['band'] == cb['band'] + st)]
    s_state = max(near, key=lambda r: r['cs'] + r['tepz'])                # 阳离子 s–Te p_z 反键态
    out = dict(soc=soc, vbm=top['E'], s=s_state['E'], s_in='价带顶' if s_state is top else '导带底', rows=rows)
    if soc:
        p = next((r for r in reversed([r for r in rows if r['vb']]) if r['bp'] >= 0.25 and r['cs'] < 0.10), None)
        out.update(p12=None if p is None else p['E'], p12_bp=None if p is None else p['bp'],
                   p12_pz=None if p is None else p['bpz'] / p['bp'])
    else:
        cbs = [r for r in rows if not r['vb']]
        pz = next((r for r in cbs if r['bpz'] >= 0.20 and r['bpxy'] < 0.05), None)
        pxy = next((r for r in cbs if r['bpxy'] >= 0.30), None)
        out.update(pz=None if pz is None else pz['E'], pxy=None if pxy is None else pxy['E'],
                   pc=None if (pz is None or pxy is None) else (pz['E'] + 2 * pxy['E']) / 3)
    return out


def runs():
    out = {}
    for m in NELECT:
        out[(m, m, True)] = (os.path.join(R, 'c2_pbe') if m == 'Pb2Sb2Te5' else os.path.join(R, m, 'c2_pbe'), m)
        out[(m, m, False)] = (os.path.join(D, f'd1_nosoc_{m}'), m)
    for comp, latt in (('Sn2Bi2Te5', 'Sn2Sb2Te5'), ('Sn2Sb2Te5', 'Sn2Bi2Te5')):
        out[(comp, latt, True)] = (os.path.join(D, f'd2_{comp}_in_{latt}_lattice'), comp)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--fig'); a = ap.parse_args()
    res = {}
    f0 = lambda v: '—' if v is None else f'{v:.0f}'
    print('# Γ 点各态相对 Σ 谷的能量（meV，PBE）')
    print(f"{'成分':10s} {'晶格':10s} {'SOC':3s} {'Γ价带顶':>7s} {'s态':>6s} {'(位置)':6s} {'p_z':>6s} {'p_xy':>6s} {'p中心':>6s} {'p1/2':>6s} {'B-p%':>5s} {'pz比例':>6s}")
    for key, (d, comp) in runs().items():
        r = gamma_states(d, comp) if os.path.isdir(d) else None
        if r is None:
            print(f'{key[0]:10s} {key[1]:10s} {"有" if key[2] else "无":3s}  — 缺少 PROCAR_OPT（{os.path.relpath(d, ROOT)}）'); continue
        res[key] = r
        print(f"{key[0]:10s} {key[1]:10s} {'有' if r['soc'] else '无':3s} {r['vbm']:7.0f} {r['s']:6.0f} {r['s_in']:6s} "
              f"{f0(r.get('pz')):>6s} {f0(r.get('pxy')):>6s} {f0(r.get('pc')):>6s} {f0(r.get('p12')):>6s} "
              f"{'' if r.get('p12_bp') is None else f'{100 * r['p12_bp']:.0f}':>5s} {'' if r.get('p12_pz') is None else f'{r['p12_pz']:.2f}':>6s}")
    print('\n# Sb → Bi 时 Γ 处 p₁/₂ 态（有 SOC 时的 Γ 价带顶）相对 Σ 谷的变化，拆成标量部分与 SOC 部分（meV）')
    dec = {}
    for A_ in ('Pb', 'Sn'):
        sb, bi = f'{A_}2Sb2Te5', f'{A_}2Bi2Te5'
        try:
            lo = {m: res[(m, m, True)]['p12'] - res[(m, m, False)]['pc'] for m in (sb, bi)}   # SOC 把 p 态从能级中心压下的量
            tot = res[(bi, bi, True)]['p12'] - res[(sb, sb, True)]['p12']
            sc = res[(bi, bi, False)]['pc'] - res[(sb, sb, False)]['pc']
        except (KeyError, TypeError):
            print(f'  {A_} 系：数据不全'); continue
        dec[A_] = dict(total=tot, scalar=sc, soc=lo[bi] - lo[sb], low_sb=lo[sb], low_bi=lo[bi])
        print(f"  {A_} 系：总变化 {tot:+.0f} = 标量部分 {sc:+.0f}（无 SOC 时 p 能级中心）+ SOC 部分 {lo[bi] - lo[sb]:+.0f}"
              f"（SOC 把 p 态压下 {-lo[sb]:.0f} → {-lo[bi]:.0f}）")
    s1, s2 = res.get(('Sn2Bi2Te5', 'Sn2Sb2Te5', True)), res.get(('Sn2Sb2Te5', 'Sn2Bi2Te5', True))
    a_, b_ = res.get(('Sn2Sb2Te5', 'Sn2Sb2Te5', True)), res.get(('Sn2Bi2Te5', 'Sn2Bi2Te5', True))
    if s1 and s2 and a_ and b_:
        print(f"\n# 晶格交换（PBE+SOC，p₁/₂ 态）：原生晶格 Sb → Bi {b_['p12'] - a_['p12']:+.0f}；SnSb 晶格中 {s1['p12'] - a_['p12']:+.0f}；"
              f"SnBi 晶格中 {b_['p12'] - s2['p12']:+.0f} meV。同一成分换晶格：SnBi {s1['p12'] - b_['p12']:+.0f}、SnSb {s2['p12'] - a_['p12']:+.0f} meV")
    if a.fig and all((m, m, s) in res for m in ('Sn2Sb2Te5', 'Sn2Bi2Te5') for s in (True, False)):
        figure(res, dec, a.fig)


def figure(res, dec, out):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import paper_figures as pf
    pf.style(); pf.DATA = os.path.join(out, 'source_data')
    C_P, C_S = pf.C_B, pf.C_A                                   # p 态橙、s 态紫，与图 5、图 6 一致
    fig, ax = plt.subplots(figsize=(pf.W1 * 1.45, 3.4))
    cols = [('Sn2Sb2Te5', False, 0.0), ('Sn2Sb2Te5', True, 1.0), ('Sn2Bi2Te5', False, 2.4), ('Sn2Bi2Te5', True, 3.4)]
    hw = 0.28; rows = []
    for m, soc, x in cols:
        r = res[(m, m, soc)]
        lv = [('s', r['s'], C_S, '-')] + ([('p$_{1/2}$', r['p12'], C_P, '-')] if soc else
                                         [('p$_z$', r['pz'], C_P, (0, (2, 1.2))), ('p$_{xy}$', r['pxy'], C_P, (0, (2, 1.2)))])
        for nm, e, c, ls in lv:
            ax.plot([x - hw, x + hw], [e, e], color=c, lw=1.8, ls=ls)
            ax.text(x + hw + 0.04 if soc else x - hw - 0.04, e, nm, color=c, fontsize=6.5, va='center', ha='left' if soc else 'right')
            rows.append([m, 'SOC' if soc else 'no SOC', nm.replace('$', '').replace('_', '').replace('{', '').replace('}', ''), round(e, 1)])
        if not soc:
            ax.plot([x - hw, x + hw], [r['pc']] * 2, color=C_P, lw=0.6, alpha=0.6)
            rows.append([m, 'no SOC', 'p center', round(r['pc'], 1)])
    for m, x0, x1 in (('Sn2Sb2Te5', 0.0, 1.0), ('Sn2Bi2Te5', 2.4, 3.4)):    # p 能级中心 → p₁/₂，s → s
        a0, a1 = res[(m, m, False)], res[(m, m, True)]
        ax.annotate('', (x1 - hw, a1['p12']), (x0 + hw, a0['pc']), arrowprops=dict(arrowstyle='->', color=C_P, lw=0.8))
        ax.plot([x0 + hw, x1 - hw], [a0['s'], a1['s']], color=C_S, lw=0.6, ls=(0, (1, 1)))
        ax.text((x0 + x1) / 2 + 0.12, a1['p12'] + 0.6 * (a0['pc'] - a1['p12']), f"SOC {a1['p12'] - a0['pc']:+.0f}".replace('-', '−'),
                ha='left', va='center', fontsize=6.5, color=C_P)
    ax.axhline(0, color='k', lw=0.6, ls=pf.DASH); ax.text(3.95, 8, 'Σ valley', fontsize=6.5, ha='right', va='bottom')
    ax.set_xticks([0, 1, 2.4, 3.4], ['no SOC', 'SOC', 'no SOC', 'SOC'])
    for x, m in ((0.5, 'Sn2Sb2Te5'), (2.9, 'Sn2Bi2Te5')):
        ax.text(x, -0.09, pf.TEX[m], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=7.5)
    ax.set_xlim(-0.75, 4.0); ax.set_ylabel(r'$E - E_\Sigma$ at $\Gamma$ (meV, PBE)')
    ax.set_title('Γ-point levels: SOC pushes the cation-p$_{1/2}$ state below the s state', fontsize=7.5, pad=3)
    pf.dump('FigS2_Gamma_levels_SOC', 'gamma_levels_PBE.csv', ['material', 'SOC', 'state', 'E_minus_Sigma_meV'], rows)
    if dec:
        pf.dump('FigS2_Gamma_levels_SOC', 'SbBi_decomposition_meV.csv', ['series', 'total', 'scalar_part', 'SOC_part', 'SOC_lowering_Sb', 'SOC_lowering_Bi'],
                [[k, round(v['total']), round(v['scalar']), round(v['soc']), round(-v['low_sb']), round(-v['low_bi'])] for k, v in dec.items()])
    pf.save(fig, out, 'FigS2_Gamma_levels_SOC')


if __name__ == '__main__':
    main()
