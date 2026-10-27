"""c2 路径能带：沿路径找价带/导带的局部极值（谷），给出距带边能量、沿路径有效质量，和（有 PROCAR_OPT 时）成分。

数据：22_能带基准_1-0/<材料>/<side>/path_eigen.txt（k_index band_index energy_eV，k 从 0 起）、可选 PROCAR_OPT(.gz)；
路径：22_能带基准_1-0/c2_path/KPOINTS_OPT；晶格：08_介电复核/源数据/<材料>/POSCAR；NELECT 取 <side>/EIGENVAL 或 --nelect。
Pb2Sb2Te5 的数据在 22_能带基准_1-0/ 根目录下。

口径：
- 价带 = 第 NELECT 条，导带 = 第 NELECT+1 条（SOC，Kramers 简并）。只看这两条带本身的局部极值，窗口 --win（meV）。
- 有效质量 m*/m_e = 7.62 / (d²E/dk²，eV·Å²)，只沿路径方向：段内极值用两侧共 5 点抛物线拟合；
  段端点（对称点）上分别向两侧各取 3 点按 E0 + c·dk² 拟合。不是张量，垂直于路径的方向未知。
- 成分：PROCAR 每条带第 1 块（电荷投影），对 Kramers 对求和；份额 = 投影 / 球内 s+p+d 总投影。
  离子顺序 B×2、Te×5、A×2（与 1-0 的 POSCAR 相同），脚本打印半芯 d 带供核对。

用法：python 10_脚本/c2_valleys.py Sn2Sb2Te5 --side c2_hse [--nelect 68] [--win 150]
"""
import argparse, gzip, os, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HB = 7.62
BREAKS = [(187, 188), (195, 196)]
SEGS = [(0, 33, 'Γ–M'), (33, 52, 'M–K'), (52, 90, 'K–Γ'), (90, 97, 'Γ–A'), (97, 130, 'A–L'), (130, 149, 'L–H'),
        (149, 187, 'H–A'), (188, 195, 'L–M'), (196, 203, 'H–K')]


def read_lat(p):
    L = open(p).read().splitlines(); s = float(L[1])
    return s * np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])


def read_kpoints_opt(p):
    L = open(p).read().splitlines(); n = int(L[1])
    return np.array([[float(x) for x in l.split()[:3]] for l in L[3:3 + n]])


def read_path_eigen(p):
    d = np.loadtxt(p, comments='#'); nk = int(d[:, 0].max()) + 1; nb = int(d[:, 1].max())
    E = np.full((nk, nb), np.nan); E[d[:, 0].astype(int), d[:, 1].astype(int) - 1] = d[:, 2]
    return E


def read_procar(p):
    L = (gzip.open(p, 'rt') if p.endswith('.gz') else open(p)).read().splitlines()
    nk, nb, ni = map(int, re.findall(r'\d+', L[1])[:3])
    E = np.zeros((nk, nb)); W = np.zeros((nk, nb, ni, 9))
    ik = ib = -1; i = 2
    while i < len(L):
        s = L[i]
        if s.startswith(' k-point'):
            ik += 1; ib = -1
        elif s.startswith('band'):
            ib += 1; E[ik, ib] = float(s.split('# energy')[1].split()[0])
            while not L[i].startswith('ion'): i += 1
            for j in range(ni):                      # 第 1 块：电荷投影
                W[ik, ib, j] = [float(x) for x in L[i + 1 + j].split()[1:10]]
            i += ni + 1
        i += 1
    return E, W


def kdist(K, B):
    x = [0.0]
    for i in range(1, len(K)):
        x.append(x[-1] + (0.0 if (i - 1, i) in BREAKS else np.linalg.norm((K[i] - K[i - 1]) @ B)))
    return np.array(x)


def side_mass(x, e, i, d):
    idx = [i + d * n for n in (1, 2, 3)]
    dx = np.array([x[j] - x[i] for j in idx]); de = np.array([e[j] - e[i] for j in idx])
    c = (dx ** 2 @ de) / (dx ** 4).sum()
    return HB / (2 * c) if abs(c) > 1e-9 else np.inf


def extrema(x, e, kind, win):
    ref = np.nanmax(e) if kind == 'max' else np.nanmin(e); out = []
    for a, b, name in SEGS:
        for i in range(a, b + 1):
            nb = [j for j in (i - 1, i + 1) if a <= j <= b]
            if kind == 'max' and all(e[i] >= e[j] for j in nb) or kind == 'min' and all(e[i] <= e[j] for j in nb):
                dE = 1000 * (e[i] - ref)
                if abs(dE) > win: continue
                if a < i - 1 and i + 1 < b:
                    p = np.polyfit(x[i - 2:i + 3] - x[i], e[i - 2:i + 3], 2); m = [(name, HB / (2 * p[0]))]
                else:                                # 对称点：两侧分别给
                    m = [(name if d > 0 else name[::-1], side_mass(x, e, i, d)) for d in (1, -1) if a <= i + 3 * d <= b]
                out.append((i, dE, m))
    merged = {}                                      # 同一 k 序号出现在相邻两段端点：合并质量
    for i, dE, m in out:
        merged.setdefault(i, [dE, []])[1].extend(m)
    return sorted(merged.items(), key=lambda t: -t[1][0] if kind == 'max' else t[1][0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('material'); ap.add_argument('--side', default='c2_hse')
    ap.add_argument('--nelect', type=int); ap.add_argument('--win', type=float, default=150)
    a = ap.parse_args()
    d = os.path.join(ROOT, '22_能带基准_1-0', '' if a.material == 'Pb2Sb2Te5' else a.material, a.side)
    nel = a.nelect or int(open(os.path.join(d, 'EIGENVAL')).read().split('\n')[5].split()[0])
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    B = 2 * np.pi * np.linalg.inv(read_lat(os.path.join(ROOT, '08_介电复核', '源数据', a.material, 'POSCAR'))).T
    x = kdist(K, B)
    E = read_path_eigen(os.path.join(d, 'path_eigen.txt'))
    if E.shape[0] != len(K): sys.exit(f'路径点数 {E.shape[0]} ≠ KPOINTS_OPT {len(K)}')
    v, c = E[:, nel - 1], E[:, nel]

    W = None
    pc = next((os.path.join(d, f) for f in ('PROCAR_OPT.gz', 'PROCAR_OPT') if os.path.exists(os.path.join(d, f))), None)
    A, Bn = a.material[:2], a.material[3:5]
    if pc:
        Ep, W = read_procar(pc)
        if Ep.shape[0] != len(K): sys.exit(f'{pc}: k 点数 {Ep.shape[0]} ≠ {len(K)}')
        dev = np.abs(Ep[:, nel - 1:nel + 1] - E[:, nel - 1:nel + 1]).max()
        dd = W[0, :, :, 4:9].sum(-1)
        semi = [f'离子{j + 1}:{Ep[0, dd[:, j] > 0.2].min():.1f}' for j in range(W.shape[2]) if (dd[:, j] > 0.2).any()]
        print(f'# PROCAR_OPT 与 path_eigen 带边两条带最大差 {1000 * dev:.2f} meV；半芯 d（Γ 最低能量，eV）：{"，".join(semi)}')
        G = {Bn: slice(0, 2), 'Te': slice(2, 7), A: slice(7, 9)}

    gap = np.nanmin(c) - np.nanmax(v)
    print(f'# {a.material} {a.side}：NELECT {nel}，路径带隙 {gap:.4f} eV，VBM #{np.nanargmax(v)}，CBM #{np.nanargmin(c)}')
    hdr = f'{"":4s}{"k#":>4s}  {"k（分数坐标）":22s}{"距带边/meV":>10s}  有效质量（沿路径方向，m_e）'
    if W is not None: hdr += f'{"":18s}' + ' '.join(f'{g}-{o:1s}' for g in (A, Bn) for o in 'sp') + ' Te-s Te-p   d'
    print(hdr)
    for kind, e, tag, bb in (('max', v, '价带', [nel - 2, nel - 1]), ('min', c, '导带', [nel, nel + 1])):
        for i, (dE, m) in extrema(x, e, kind, a.win):
            ms = '，'.join(f'{s} {abs(q):.3f}' for s, q in m)
            line = f'{tag} {i:4d}  {str(np.round(K[i], 3).tolist()):22s}{dE:+10.1f}  {ms:34s}'
            if W is not None:
                w = W[i, bb].sum(0); t = w.sum()
                f = lambda g, o: 100 * (w[G[g], 0].sum() if o == 's' else w[G[g], 1:4].sum()) / t
                line += ' '.join(f'{f(g, o):4.0f}' for g in (A, Bn) for o in 'sp') + f' {f("Te", "s"):4.0f} {f("Te", "p"):4.0f} {100 * w[:, 4:9].sum() / t:3.0f}'
            print(line)


if __name__ == '__main__':
    main()
