"""22_：带边附近几个 k 点的元素-轨道成分（SOC PROCAR，Kramers 对求和）。

用法：python 10_脚本/band_edge_composition.py 22_能带基准_1-0/<材料> [--side c1_HSE_SOC] [--order B,Te,A --counts 2,5,2]
读 <目录>/<side>/PROCAR.gz（或 PROCAR）与同目录 EIGENVAL（取 NELECT）。
口径：SOC 下每条带有 4 个块，只取第 1 块（电荷投影）；价带取第 NELECT−1、NELECT 条（Kramers 对），导带取第 NELECT+1、NELECT+2 条；
份额 = 该类投影 / 球内 s+p+d 总投影（%）。离子顺序默认与 1-0 的 POSCAR 相同：B（Sb/Bi）×2、Te×5、A（Pb/Sn）×2，
脚本先打印 Γ 点各离子的半芯 d 带能量供核对（Bi 5d 最深、Pb 5d 次之、Sn 4d；Sb、Te 的 POTCAR 无 d 半芯）。
"""
import argparse, gzip, os, re, numpy as np

PTS = [((0, 0, 0), 'Γ'), ((0, 0, .5), 'A'), ((1/6, 0, 0), '(1/6,0,0)'), ((1/3, 0, .5), '(1/3,0,1/2)'),
       ((1/12, 0, .5), '(1/12,0,1/2)'), ((1/12, 0, 0), '(1/12,0,0)')]


def read_procar(p):
    L = (gzip.open(p, 'rt') if p.endswith('.gz') else open(p)).read().splitlines()
    nk, nb, ni = map(int, re.findall(r'\d+', L[1])[:3])
    E = np.zeros((nk, nb)); W = np.zeros((nk, nb, ni, 9)); K = np.zeros((nk, 3))
    ik = ib = -1; i = 2
    while i < len(L):
        s = L[i]
        if s.startswith(' k-point'):
            ik += 1; ib = -1
            K[ik] = [float(x) for x in re.findall(r'-?\d+\.\d+', s.split(':')[1])[:3]]
        elif s.startswith('band'):
            ib += 1; E[ik, ib] = float(s.split('# energy')[1].split()[0])
            while not L[i].startswith('ion'): i += 1
            for j in range(ni):                      # 第 1 块：电荷投影
                W[ik, ib, j] = [float(x) for x in L[i + 1 + j].split()[1:10]]
            i += ni + 1
        i += 1
    return E, W, K


ap = argparse.ArgumentParser()
ap.add_argument('dir'); ap.add_argument('--side', default='c1_HSE_SOC')
ap.add_argument('--order', default='B,Te,A'); ap.add_argument('--counts', default='2,5,2')
a = ap.parse_args()
d = os.path.join(a.dir, a.side)
p = os.path.join(d, 'PROCAR.gz') if os.path.exists(os.path.join(d, 'PROCAR.gz')) else os.path.join(d, 'PROCAR')
nel = int(open(os.path.join(d, 'EIGENVAL')).read().split('\n')[5].split()[0])
E, W, K = read_procar(p)
groups, start = {}, 0
for g, c in zip(a.order.split(','), map(int, a.counts.split(','))):
    groups[g] = slice(start, start + c); start += c
assert start == W.shape[2], (start, W.shape[2])

print(f'{a.dir}/{a.side}   NELECT {nel}   离子分组 {a.order} × {a.counts}')
dd = W[0, :, :, 4:9].sum(-1)
for ion in range(W.shape[2]):
    b = np.where(dd[:, ion] > 0.2)[0]
    if len(b): print(f'   核对：离子 {ion + 1} 有 {len(b)} 条半芯 d 带，Γ 处 {E[0, b].min():.1f} … {E[0, b].max():.1f} eV')
v, c = nel - 1, nel                                  # 0-based：价带顶第 NELECT 条，导带底第 NELECT+1 条
vbm, cbm = E[:, v].max(), E[:, c].min()
print(f'   带隙 {cbm - vbm:.4f} eV；份额 = 投影 / 球内 s+p+d 总投影（%），Kramers 对求和\n')
cols = [(g, o) for g in groups if g != 'Te' for o in 'sp'] + [('Te', 's'), ('Te', 'p')]
print(f'   {"k 点":14s}{"带":4s}{"E−边(meV)":>10s}  ' + '  '.join(f'{g}-{o}' for g, o in cols) + '    d')
for kk, name in PTS:
    ik = next((i for i, k in enumerate(K) if np.abs(np.mod(k - np.array(kk) + .5, 1) - .5).max() < 2e-3), None)
    if ik is None: continue
    for tag, bands, ref in (('VB', [v - 1, v], vbm), ('CB', [c, c + 1], cbm)):
        w = W[ik, bands].sum(0); t = w.sum()
        f = lambda g, o: 100 * (w[groups[g], 0].sum() if o == 's' else w[groups[g], 1:4].sum()) / t
        print(f'   {name:14s}{tag:4s}{1000 * (E[ik, bands].mean() - ref):10.1f}  ' +
              '  '.join(f'{f(g, o):4.0f}' for g, o in cols) + f'  {100 * w[:, 4:9].sum() / t:3.0f}')
