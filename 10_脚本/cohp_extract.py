#!/usr/bin/env python3
"""从 COGITO 的交互图（06_COGITO成键分析与COHP/图_交互式COHP/<材料>_COHP交互图.html，PBE、不加 SOC、λ=1）取出逐键、逐轨道的 COHP(E)
（2026-11-05，用于报告的 COHP 图：孤对电子是阳离子 s–Te p 的反键态）。

交互图里每根键一组曲线（customdata 为键名，如 SnTe1），组末一条 'bond COHP' 为该键的总 COHP；轨道对曲线名如 'Sn s - Te p'。
组内长度 435 的曲线是逐项累加的堆叠边界，相邻相减得各轨道对；总 COHP 积分到 E_F 与 bond_info.txt 的 ICOHP 恰差 1.50 倍\n（有的组含两根等长键时为 3.0 倍），这里换算成每根键、每 eV。
符号：COGITO 的 COHP 为负表示成键；本脚本另给 −COHP（成键为正）。
输出：06_COGITO成键分析与COHP/审核与重绘/轨道分辨COHP/<材料>_<键>.csv（E−E_F，−COHP：总、s–p、p–s、p–p、s–s）与 summary.csv。

--run（2026-11-05 补算 f）：
  pbe      原来的 PBE 运行（图_交互式COHP/、结果/），输出 轨道分辨COHP/（默认）
  hse      HSE06 不加 SOC：补算_HSE与SOC/f_cogito_hse/<材料>/cogito_hse/（bond_info.txt 与 COGITOpost 的 COHP 交互图 *.html），输出 轨道分辨COHP_HSE/
  pbe_f    同一批补算里的 PBE 对照：…/cogito_pbe/，输出 轨道分辨COHP_PBE_对照/
交互图曲线单位与 eV/键之比 K 对每次运行重新核对（单键组的 ∫COHP / ICOHP），与 1.5038 相差超过 1% 时报错。
"""
import argparse, base64, csv, glob, json, os, re
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, '06_COGITO成键分析与COHP'); OUT = os.path.join(D, '审核与重绘', '轨道分辨COHP')
MATS = ('Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5')
arr = lambda v: np.frombuffer(base64.b64decode(v['bdata']), dtype=v['dtype']).astype(float) if isinstance(v, dict) else np.array(v, float)


FDIR = os.path.join(D, '补算_HSE与SOC', 'f_cogito_hse')
RUNS = {'pbe': (lambda m: os.path.join(D, '图_交互式COHP', f'{m}_COHP交互图.html'), lambda m: os.path.join(D, '结果', m, 'bond_info.txt'),
                os.path.join(D, '审核与重绘', '轨道分辨COHP')),
        'hse': (lambda m: cohp_html(os.path.join(FDIR, m, 'cogito_hse')), lambda m: os.path.join(FDIR, m, 'cogito_hse', 'bond_info.txt'),
                os.path.join(D, '审核与重绘', '轨道分辨COHP_HSE')),
        'pbe_f': (lambda m: cohp_html(os.path.join(FDIR, m, 'cogito_pbe')), lambda m: os.path.join(FDIR, m, 'cogito_pbe', 'bond_info.txt'),
                  os.path.join(D, '审核与重绘', '轨道分辨COHP_PBE_对照'))}
HTML, BINFO = RUNS['pbe'][0], RUNS['pbe'][1]


def cohp_html(d):
    """目录 d 中 COGITOpost 写出的 COHP 交互图：含 'bond COHP' 曲线的那个 html（晶体键图不含）。"""
    fs = [f for f in sorted(glob.glob(os.path.join(d, '*.html'))) if '"bond COHP"' in open(f, encoding='utf-8', errors='ignore').read()]
    if len(fs) != 1: raise SystemExit(f'[ERROR] {d} 中含 COHP 曲线的 html 有 {len(fs)} 个：{fs}')
    return fs[0]


def groups(m):
    t = open(HTML(m), encoding='utf-8').read()
    j = t.index('[', t.find('Plotly.newPlot')); data, _ = json.JSONDecoder().raw_decode(t[j:])
    out, cur = {}, []
    for d in (d for d in data if d.get('type') == 'scatter'):
        if d.get('name') == 'bond COHP':
            out[d['customdata'][0]] = (cur, d); cur = []
        elif d.get('name'): cur.append(d)
    return out


def bond_info(m):
    rows = []
    for l in open(BINFO(m)).read().splitlines():
        p = l.split()
        try: rows.append((p[0], p[1], float(p[2]), float(p[3])))
        except (ValueError, IndexError): pass
    return rows


K = 1.5038       # 交互图曲线单位 / (eV/键)：B–Te、Sn–Te 单键组的 ∫COHP 与 bond_info 的 ICOHP 之比，四种材料一致


def contributions(g, ny):
    """组内长度为 ny 的曲线是逐项累加的堆叠边界（第 k 条 = 前 k 个轨道对之和，最后一条 = 总 COHP），相邻相减得各轨道对。"""
    cum = [(d['name'], arr(d['x'])) for d in g if len(arr(d['x'])) == ny]
    out, prev = {}, 0.0
    for name, x in cum:
        a, b = [q.split() for q in name.split(' - ')]
        cat, te = (a, b) if b[0] == 'Te' else (b, a)
        out[f'cation_{cat[1]}-Te_{te[1]}'] = x - prev; prev = x
    return out, prev


def main():
    global HTML, BINFO, OUT
    ap = argparse.ArgumentParser(); ap.add_argument('--run', choices=list(RUNS), default='pbe'); a = ap.parse_args()
    HTML, BINFO, OUT = RUNS[a.run]
    os.makedirs(OUT, exist_ok=True); summ = []
    for m in MATS:                                                     # 核对 K：最短的 B–Te 键在 COGITO 里总是单键一组
        G, info = groups(m), bond_info(m); B = m[3:5]
        bt = min((r for r in info if {r[0], r[1]} == {B, 'Te'}), key=lambda r: r[2])
        ks = []
        for name, (g, b) in G.items():
            if re.fullmatch(f'{B}Te\\d+', name):
                y = arr(b['y']); tot = arr(b['x']); ks.append(float(np.trapezoid(tot[y <= 0], y[y <= 0])) / bt[3])
        k = min(ks, key=lambda v: abs(v - K))
        print(f'[K] {m}: 最短 B–Te（{bt[2]:.2f} Å）∫COHP / ICOHP = {k:.4f}（PBE 运行为 {K}）')
        if abs(k / K - 1) > 0.01: raise SystemExit(f'[ERROR] {m} 的 K 与 PBE 运行差 {100 * (k / K - 1):+.1f}%，交互图的单位可能变了，先查清再用')
    for m in MATS:
        A, B = m[:2], m[3:5]; G = groups(m); info = bond_info(m)
        for site, el in (('A', A), ('B', B)):
            ref = [r for r in info if {r[0], r[1]} == {el, 'Te'} and r[3] < -0.3]
            for name, (g, b) in sorted(G.items()):
                if not re.fullmatch(f'{el}Te\\d+', name): continue
                y = arr(b['y']); occ = y <= 0; tot = arr(b['x']); I = float(np.trapezoid(tot[occ], y[occ]))
                if I > -0.3: continue                                           # 远距离的弱相互作用不要
                # 这组含几根键、是哪几根：组内 n 根键（可重复取同一行）的 ICOHP 平均最接近 ∫/(K·n) 的组合。
                # 2026-11-05 审核：Pb 系的 A 位组含两根长度几乎相同的键（3.20 + 3.21、3.21 + 3.22 Å），曲线是两者的平均，不是单根键
                from itertools import combinations_with_replacement as cwr
                best = min(((abs(I / (K * n) - np.mean([q[3] for q in c])), n, c) for n in (1, 2, 3) for c in cwr(ref, n)), key=lambda q: q[0])
                _, n, grp = best; sc = 1.0 / (K * n); r = min(grp, key=lambda q: q[2])
                lens = '+'.join(f'{q[2]:.2f}' for q in sorted(grp, key=lambda q: q[2]))
                parts, last = contributions(g, len(y))
                assert np.allclose(last, tot, atol=1e-6), (m, name)
                cols = {'total': -sc * tot}
                for k in ('cation_s-Te_p', 'cation_p-Te_s', 'cation_p-Te_p', 'cation_s-Te_s', 'cation_d-Te_p', 'cation_d-Te_s'):
                    cols[k] = -sc * parts.get(k, 0 * tot)
                f = os.path.join(OUT, f'{m}_{site}-Te_{r[2]:.2f}A.csv')
                with open(f, 'w', newline='', encoding='utf-8-sig') as fh:
                    # 单位：COHP 对能量（eV）积分得 ICOHP（eV/键），所以 COHP 本身是“每根键”的无量纲量（eV/eV），不是 eV⁻¹
                    w = csv.writer(fh); w.writerow(['E-EF_eV'] + [f'-COHP_{k}_per_bond' for k in cols])
                    for i in range(len(y)): w.writerow([round(float(y[i]), 4)] + [round(float(cols[k][i]), 5) for k in cols])
                integ = lambda v, lo=-99: float(np.trapezoid(v[occ & (y >= lo)], y[occ & (y >= lo)]))
                summ.append([m, site, el, f'{r[2]:.2f}', lens, name, n, round(float(np.mean([q[3] for q in grp])), 4), round(-integ(cols['total']), 4)]
                            + [round(-integ(cols[k]), 4) for k in list(cols)[1:]] + [round(-integ(cols['cation_s-Te_p'], -1.0), 4)])
                print(f'{m} {el}–Te {lens} Å ({name}, {n} bond/group): ICOHP total {-integ(cols["total"]):+.3f} (bond_info {np.mean([q[3] for q in grp]):+.3f}); '
                      f's–p {-integ(cols["cation_s-Te_p"]):+.3f}, p–p {-integ(cols["cation_p-Te_p"]):+.3f}, s–p in top 1 eV {-integ(cols["cation_s-Te_p"], -1.0):+.3f} eV')
    with open(os.path.join(OUT, 'summary.csv'), 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh)
        w.writerow(['material', 'site', 'cation', 'bond_length_A', 'bond_lengths_in_group_A', 'cogito_bond_id', 'bonds_in_group',
                    'ICOHP_bond_info_eV (mean over the bonds in the group)', 'ICOHP_total_eV',
                    'ICOHP_cation_s-Te_p_eV', 'ICOHP_cation_p-Te_s_eV', 'ICOHP_cation_p-Te_p_eV', 'ICOHP_cation_s-Te_s_eV', 'ICOHP_cation_d-Te_p_eV',
                    'ICOHP_cation_d-Te_s_eV', 'ICOHP_cation_s-Te_p_within_1eV_below_EF_eV'])
        w.writerows(summ)


if __name__ == '__main__':
    main()
