"""SOC 机制检验（22_能带基准_1-0/SOC机制检验/）：Bi 为什么把 Γ 处的价带态压下去。

① d1_nosoc_<材料>：PBE、不加 SOC，同一条 c2 路径；与现有的 PBE+SOC（<材料>/c2_pbe）对比。
② d2_<成分>_in_<晶格>_lattice：PBE+SOC，把 B 位 Sb ↔ Bi 换进另一种材料的晶格和坐标，分开化学与体积效应。

每个计算读出（能量都相对 Σ 谷，即 Γ–M 段内部最高价带的极大，meV）：
- Γ 价带顶（Γ 处最高的价带）；
- Γ 处的 “B 位 p 态”：Γ 处 B 位 p 投影 ≥ 25% 的最高价带，即 SOC 下的 p₁/₂ 型态；
  同时给出它的 B 位 p 成分与 p_z /(p_x + p_y + p_z)（j = 1/2 的纯 p 态为 1/3）；
- Γ 处最低导带（反转后的 s 型态）。
然后算 Sb → Bi 的变化，以及“有 SOC 的变化 − 无 SOC 的变化”（双重差分，A 位与 Te 的 SOC 在其中基本抵消）。

用法：python 10_脚本/soc_mechanism.py            （缺少的计算自动跳过）
"""
import os, sys
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


def analyse(d, comp):
    f = os.path.join(d, 'PROCAR_OPT.gz')
    if not os.path.exists(f): f = f[:-3]
    if not os.path.exists(f): return None
    E, W = cv.read_procar(f); soc = soc_of(d)
    nv = NELECT[comp] if soc else NELECT[comp] // 2            # 价带数（SOC 下每条带 1 个电子）
    step = 2 if soc else 1                                     # SOC 下 Kramers 对，取一对中的上面一条
    top = E[:, nv - 1]
    esig = top[list(SIGMA_SEG)].max()
    w = W[0]; tot = w.sum(axis=(1, 2))                          # Γ 点各带的总投影
    bp = w[:, B, 1:4].sum(axis=(1, 2)) / tot; pz = w[:, B, 2].sum(1) / np.maximum(w[:, B, 1:4].sum(axis=(1, 2)), 1e-12)
    ib = next((i for i in range(nv - 1, nv - 1 - 16 * step, -step) if bp[i] >= 0.25), None)
    return dict(soc=soc, gvb=1000 * (top[0] - esig), cb=1000 * (E[0, nv] - esig),
                bp_E=None if ib is None else 1000 * (E[0, ib] - esig), bp_w=None if ib is None else 100 * bp[ib],
                bp_pz=None if ib is None else pz[ib], bp_depth=None if ib is None else (nv - 1 - ib) // step)


def runs():
    out = {}
    for m in NELECT:
        out[(m, m, True)] = (os.path.join(R, 'c2_pbe') if m == 'Pb2Sb2Te5' else os.path.join(R, m, 'c2_pbe'), m)
        out[(m, m, False)] = (os.path.join(D, f'd1_nosoc_{m}'), m)
    for comp, latt in (('Sn2Bi2Te5', 'Sn2Sb2Te5'), ('Sn2Sb2Te5', 'Sn2Bi2Te5')):
        out[(comp, latt, True)] = (os.path.join(D, f'd2_{comp}_in_{latt}_lattice'), comp)
    return out


def main():
    res = {}
    print(f"{'成分':10s} {'晶格':10s} {'SOC':4s} {'Γ价带顶':>8s} {'Γ B位p态':>9s} {'B-p%':>6s} {'pz比例':>7s} {'Γ导带':>7s}   （meV，相对 Σ 谷）")
    for key, (d, comp) in runs().items():
        r = analyse(d, comp) if os.path.isdir(d) else None
        if r is None:
            print(f'{key[0]:10s} {key[1]:10s} {"有" if key[2] else "无":4s}  — 缺少 PROCAR_OPT（{os.path.relpath(d, ROOT)}）'); continue
        res[key] = r
        f = lambda v, n=0: '—' if v is None else f'{v:.{n}f}'
        print(f"{key[0]:10s} {key[1]:10s} {'有' if key[2] else '无':4s} {r['gvb']:8.0f} {f(r['bp_E']):>9s} {f(r['bp_w']):>6s} {f(r['bp_pz'], 2):>7s} {r['cb']:7.0f}"
              + (f"  （B 位 p 态在价带顶下第 {r['bp_depth']} 条）" if r['bp_depth'] else ''))
    print('\n# Sb → Bi：Γ 处 B 位 p 态相对 Σ 谷的变化（meV），以及双重差分')
    for A_ in ('Pb', 'Sn'):
        sb, bi = f'{A_}2Sb2Te5', f'{A_}2Bi2Te5'
        dd = {}
        for soc in (True, False):
            a, b = res.get((sb, sb, soc)), res.get((bi, bi, soc))
            if a and b and a['bp_E'] is not None and b['bp_E'] is not None: dd[soc] = b['bp_E'] - a['bp_E']
        line = '；'.join(f"{'有' if s else '无'} SOC {v:+.0f}" for s, v in dd.items())
        if len(dd) == 2: line += f"；双重差分（Bi 比 Sb 多出的 SOC 的贡献）{dd[True] - dd[False]:+.0f}"
        print(f'  {A_} 系：{line or "数据不全"}')
    s1, s2 = res.get(('Sn2Bi2Te5', 'Sn2Sb2Te5', True)), res.get(('Sn2Sb2Te5', 'Sn2Bi2Te5', True))
    a, b = res.get(('Sn2Sb2Te5', 'Sn2Sb2Te5', True)), res.get(('Sn2Bi2Te5', 'Sn2Bi2Te5', True))
    if s1 and a and s1['bp_E'] is not None:
        print(f"\n# 晶格交换（PBE+SOC）：SnSb 晶格中 Sb → Bi {s1['bp_E'] - a['bp_E']:+.0f} meV" +
              (f"；SnBi 晶格中 Sb → Bi {b['bp_E'] - s2['bp_E']:+.0f} meV" if s2 and b and s2['bp_E'] is not None else '') +
              "。与原生晶格的变化接近，说明主要是化学（B 位元素）效应；差得多则体积效应也重要")


if __name__ == '__main__':
    main()
