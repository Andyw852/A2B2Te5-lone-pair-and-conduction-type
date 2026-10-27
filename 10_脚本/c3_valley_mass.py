"""c3：离轴能谷的完整有效质量张量（HSE，KPOINTS_OPT 小星形）。

c2 路径只给出沿路径方向的质量。对落在 Γ–M、A–L 线上的谷，镜面（含该线与 c 轴）是对称操作：
面内垂直方向是主轴；沿线方向与 c 方向可以耦合，要补 (线, c) 平面的 4 个对角点。所以每个谷取 11 个点：
中心、沿线 ±δ、面内垂直 ±δ、沿 c ±δz、(线, c) 对角 (±δ, ±δz)。另加若干 c2 路径点作校验。

用法：
  生成：python 10_脚本/c3_valley_mass.py make Sn2Bi2Te5 --band 88 --kind v --valleys 11,101,119 --out 22_能带基准_1-0/Sn2Bi2Te5/c3_valley
  拟合：python 10_脚本/c3_valley_mass.py fit 22_能带基准_1-0/Sn2Bi2Te5/c3_valley
谷心由 c2 HSE 路径上该序号附近 5 点的抛物线拟合定出（价带取极大、导带取极小）。
拟合优先读 <目录>/PROCAR_OPT(.gz) 的能量（8 位小数）；没有时读 c3_eigen.txt（三列 k_index band_index energy_eV，
只有 4 位小数，即 0.1 meV，沿 c 的重质量会有 5% 量级的舍入误差）。
"""
import argparse, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HB = 7.62                                    # ħ²/m_e，eV·Å²
BREAKS = [(187, 188), (195, 196)]
SEG_ENDS = [0, 33, 52, 90, 97, 130, 149, 187, 188, 195, 196, 203]


def read_lat(p):
    L = open(p).read().splitlines(); s = float(L[1])
    return s * np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])


def read_kpoints_opt(p):
    L = open(p).read().splitlines(); n = int(L[1])
    return np.array([[float(x) for x in l.split()[:3]] for l in L[3:3 + n]])


def read_eig(p):
    d = np.loadtxt(p, comments='#'); nk = int(d[:, 0].max()) + 1; nb = int(d[:, 1].max())
    E = np.full((nk, nb), np.nan); E[d[:, 0].astype(int), d[:, 1].astype(int) - 1] = d[:, 2]
    return E


def mat_dir(m):
    return os.path.join(ROOT, '22_能带基准_1-0', '' if m == 'Pb2Sb2Te5' else m)


def make(a):
    B = 2 * np.pi * np.linalg.inv(read_lat(os.path.join(ROOT, '08_介电复核', '源数据', a.material, 'POSCAR'))).T
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    E = read_eig(os.path.join(mat_dir(a.material), 'c2_hse', 'path_eigen.txt'))[:, a.band - 1]
    z = np.array([0, 0, 1.0])
    pts, meta = [], []
    for idx in map(int, a.valleys.split(',')):
        if idx in SEG_ENDS or idx - 2 < 0 or idx + 2 >= len(K):
            sys.exit(f'#{idx} 是线段端点，这里只处理线段内部的谷')
        t = np.arange(-2, 3)
        c = np.polyfit(t, E[idx - 2:idx + 3], 2); t0 = -c[1] / (2 * c[0])
        if (c[0] < 0) != (a.kind == 'v'): sys.exit(f'#{idx} 附近不是{"极大" if a.kind == "v" else "极小"}')
        if abs(t0) > 1: sys.exit(f'#{idx} 附近的极值不在 ±1 步内（t0 = {t0:.2f}），先换序号')
        step = K[idx + 1] - K[idx]
        k0 = K[idx] + t0 * step
        eline = step @ B; eline /= np.linalg.norm(eline)
        if abs(eline[2]) > 1e-8: sys.exit(f'#{idx} 所在线段不在面内')
        eperp = np.cross(z, eline)
        to_frac = lambda v: v @ np.linalg.inv(B)
        d, dz = a.delta, a.delta_z
        stencil = [('0', 0, 0, 0), ('l+', d, 0, 0), ('l-', -d, 0, 0), ('p+', 0, d, 0), ('p-', 0, -d, 0),
                   ('z+', 0, 0, dz), ('z-', 0, 0, -dz), ('l+z+', d, 0, dz), ('l-z-', -d, 0, -dz),
                   ('l+z-', d, 0, -dz), ('l-z+', -d, 0, dz)]
        start = len(pts)
        for name, dl, dp, dzz in stencil:
            pts.append(k0 + to_frac(dl * eline + dp * eperp + dzz * z))
        meta.append({'path_index': idx, 'center_frac': k0.tolist(), 't0_steps': t0, 'e_line': eline.tolist(),
                     'e_perp': eperp.tolist(), 'delta': d, 'delta_z': dz, 'first_k': start,
                     'stencil': [s[0] for s in stencil]})
    checks = sorted({0, *map(int, a.valleys.split(','))})
    chk = []
    for i in checks:
        chk.append({'path_index': i, 'k_index': len(pts)}); pts.append(K[i])
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'KPOINTS_OPT'), 'w') as f:
        f.write(f'c3 valley stencils for {a.material} (band {a.band}); centers from c2 HSE path; see c3_points.json\n')
        f.write(f'{len(pts)}\nReciprocal\n')
        for p in pts: f.write(''.join(f'{x:+.14f}  ' for x in p) + '     1\n')
    json.dump({'material': a.material, 'band': a.band, 'kind': a.kind, 'valleys': meta, 'checks': chk, 'n_kpoints': len(pts),
               'B_rows_cart_invA': B.tolist()}, open(os.path.join(a.out, 'c3_points.json'), 'w'), indent=1, ensure_ascii=False)
    print(f'{len(pts)} 个 k 点 → {a.out}/KPOINTS_OPT')
    for m in meta:
        print(f"  谷 #{m['path_index']}：中心 {np.round(m['center_frac'], 5).tolist()}（偏离路径点 {m['t0_steps']:+.2f} 步），δ = {m['delta']}、δz = {m['delta_z']} Å⁻¹")


def fit(a):
    meta = json.load(open(os.path.join(a.dir, 'c3_points.json')))
    pc = next((os.path.join(a.dir, f) for f in ('PROCAR_OPT.gz', 'PROCAR_OPT') if os.path.exists(os.path.join(a.dir, f))), None)
    if pc:
        import gzip
        E = []
        for line in (gzip.open(pc, 'rt') if pc.endswith('.gz') else open(pc)):
            if line.startswith(' k-point'): E.append([])
            elif line.startswith('band'): E[-1].append(float(line.split('# energy')[1].split()[0]))
        E = np.array(E); print(f'# 能量取自 {os.path.basename(pc)}（8 位小数）')
    else:
        E = read_eig(os.path.join(a.dir, 'c3_eigen.txt')); print('# 能量取自 c3_eigen.txt（4 位小数）')
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"c3_eigen.txt 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    b = meta['band'] - 1
    P = read_eig(os.path.join(mat_dir(meta['material']), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], b - 3:b + 5] - P[c['path_index'], b - 3:b + 5]).max() for c in meta['checks'])
    print(f"# {meta['material']} 第 {meta['band']} 条带；校验点与 c2 路径（带边 ±4 条带）最大差 {1000 * dev:.3f} meV")
    for v in meta['valleys']:
        e = {s: E[v['first_k'] + i, b] for i, s in enumerate(v['stencil'])}
        d, dz = v['delta'], v['delta_z']
        H = np.zeros((3, 3))                       # 基：线、面内垂直、c
        H[0, 0] = (e['l+'] + e['l-'] - 2 * e['0']) / d ** 2
        H[1, 1] = (e['p+'] + e['p-'] - 2 * e['0']) / d ** 2
        H[2, 2] = (e['z+'] + e['z-'] - 2 * e['0']) / dz ** 2
        H[0, 2] = H[2, 0] = (e['l+z+'] + e['l-z-'] - e['l+z-'] - e['l-z+']) / (4 * d * dz)
        g = [(e['l+'] - e['l-']) / (2 * d), (e['p+'] - e['p-']) / (2 * d), (e['z+'] - e['z-']) / (2 * dz)]
        w, U = np.linalg.eigh(H)
        m = (-1 if meta['kind'] == 'v' else 1) * HB / w
        md = np.prod(np.abs(m)) ** (1 / 3)
        print(f"谷 #{v['path_index']}：中心 {np.round(v['center_frac'], 4).tolist()}，E0 = {e['0']:.5f} eV")
        print(f"   梯度（eV·Å，应≈0）：线 {g[0]:+.4f}、面内垂直 {g[1]:+.4f}（镜面要求严格为 0）、c {g[2]:+.4f}")
        print(f"   曲率（eV·Å²）：线 {H[0, 0]:+.3f}、面内垂直 {H[1, 1]:+.3f}、c {H[2, 2]:+.3f}、线–c 交叉 {H[0, 2]:+.3f}")
        print(f"   主质量（m_e，正 = 该带的正常载流子质量）：{', '.join(f'{x:.3f}' for x in m)}；"
              f"{'三个主曲率同号，是真极值' if np.all(np.sign(w) == np.sign(w[0])) else '主曲率异号：鞍点，不是谷'}；态密度质量 {md:.3f}")
        print('   主轴（在 线、面内垂直、c 基下）：' + '；'.join(str(np.round(U[:, i], 3).tolist()) for i in range(3)))


ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd', required=True)
s = sub.add_parser('make'); s.add_argument('material'); s.add_argument('--band', type=int, required=True)
s.add_argument('--kind', choices=('v', 'c'), required=True, help='v：价带谷（找极大），c：导带谷（找极小）')
s.add_argument('--valleys', required=True, help='c2 路径上谷的序号，逗号分隔')
s.add_argument('--delta', type=float, default=0.02, help='面内步长，Å⁻¹')
s.add_argument('--delta-z', type=float, default=0.03, help='c 方向步长，Å⁻¹')
s.add_argument('--out', required=True)
s = sub.add_parser('fit'); s.add_argument('dir')
a = ap.parse_args()
make(a) if a.cmd == 'make' else fit(a)
