"""c4 / c5：离轴能谷在镜面内的真实极值与完整质量张量（HSE，KPOINTS_OPT）。

背景：c3 表明 SnBi 在 Γ–M、A–L 线上的价带谷，在路径上的点都不是三维极值。这两条线在镜面内、不在二重轴上，
对称性不要求沿 c 的梯度为零，真实极值在镜面内的一般 kz 处（22_能带基准_1-0/Sn2Bi2Te5/c3_valley/README.md）。
K–Γ、H–A 线在二重轴上，那里的谷沿 c 的梯度由对称性为零，不需要重算。

镜面 = 含 b1（Γ–M 方向）与 c 轴的平面，倒格矢分数坐标 (u, 0, w)。(u, w) 与 (−u, −w) 由时间反演加反演等价，
所以 u ≤ 0 的半平面加完整的 w 周期就覆盖了整个镜面。

三个子命令：
  grid      c4：镜面内均匀网格 u ∈ [--umin, 0]（步长 --du）× w ∈ [−1/2, 1/2)（--nw 等分），另加 c2 路径校验点。
  analyze   读 c4 的 PROCAR_OPT(.gz)：找价带（第 NELECT 条）的局部极大、导带（第 NELECT+1 条）的局部极小，
            每个用 3×3 邻域的二次曲面拟合定出驻点与能量；给出真实价带顶、导带底和带隙；
            为离带边 --win-v 内的价带极大、--win-c 内的导带极小写出 c5 的星形 KPOINTS_OPT（与 c3 相同的 11 点星形）。
  fit       读 c5 的 PROCAR_OPT(.gz)：每个谷的梯度、曲率、主质量；检查驻点是否落在采样半径内（同 multivalley_seebeck.c3_verdict）。

用法（Sn2Bi2Te5）：
  python 10_脚本/c4_mirror_scan.py grid Sn2Bi2Te5 --out 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror
  python 10_脚本/c4_mirror_scan.py analyze 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror      # 同时写出 ../c5_star/KPOINTS_OPT
  python 10_脚本/c4_mirror_scan.py fit 22_能带基准_1-0/Sn2Bi2Te5/c5_star
只依赖 numpy。
"""
import argparse, gzip, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HB = 7.62                                     # ħ²/m_e，eV·Å²
NEL = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
CHECKS = (0, 11, 101, 119)                    # c2 路径校验点：Γ 与三个镜面内的离轴谷
CATION_IONS = (0, 1, 7, 8)                    # 离子 1–2 为 B（Sb/Bi），8–9 为 A（Pb/Sn）
STENCIL = [('0', 0, 0, 0), ('l+', 1, 0, 0), ('l-', -1, 0, 0), ('p+', 0, 1, 0), ('p-', 0, -1, 0),
           ('z+', 0, 0, 1), ('z-', 0, 0, -1), ('l+z+', 1, 0, 1), ('l-z-', -1, 0, -1), ('l+z-', 1, 0, -1), ('l-z+', -1, 0, 1)]


def read_lat(p):
    L = open(p).read().splitlines(); s = float(L[1])
    return s * np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])


def recip(material):
    return 2 * np.pi * np.linalg.inv(read_lat(os.path.join(ROOT, '08_介电复核', '源数据', material, 'POSCAR'))).T


def mat_dir(m):
    return os.path.join(ROOT, '22_能带基准_1-0', '' if m == 'Pb2Sb2Te5' else m)


def read_kpoints_opt(p):
    L = open(p).read().splitlines(); n = int(L[1])
    return np.array([[float(x) for x in l.split()[:3]] for l in L[3:3 + n]])


def read_path_eigen(p):
    d = np.loadtxt(p, comments='#'); nk = int(d[:, 0].max()) + 1; nb = int(d[:, 1].max())
    E = np.full((nk, nb), np.nan); E[d[:, 0].astype(int), d[:, 1].astype(int) - 1] = d[:, 2]
    return E


def read_procar(p):
    """返回能量 E[k, band] 与阳离子 s 成分 fs[k, band]（电荷投影第 1 块，占球内总投影的比例）。"""
    L = (gzip.open(p, 'rt') if p.endswith('.gz') else open(p)).read().splitlines()
    nk, nb, ni = (int(x) for x in L[1].replace(':', ' ').split() if x.isdigit())
    E = np.zeros((nk, nb)); fs = np.zeros((nk, nb))
    ik = ib = -1; i = 2
    while i < len(L):
        s = L[i]
        if s.startswith(' k-point'):
            ik += 1; ib = -1
        elif s.startswith('band'):
            ib += 1; E[ik, ib] = float(s.split('# energy')[1].split()[0])
            while not L[i].startswith('ion'): i += 1
            W = np.array([[float(x) for x in L[i + 1 + j].split()[1:10]] for j in range(ni)])
            fs[ik, ib] = W[list(CATION_IONS), 0].sum() / max(W.sum(), 1e-9)
            i += ni
        i += 1
    if ik + 1 != nk: sys.exit(f'PROCAR 头部写 {nk} 个 k 点，实际读到 {ik + 1}')
    return E, fs


def find_procar(d):
    p = next((os.path.join(d, f) for f in ('PROCAR_OPT.gz', 'PROCAR_OPT') if os.path.exists(os.path.join(d, f))), None)
    if p is None: sys.exit(f'{d} 里没有 PROCAR_OPT(.gz)')
    return p


def frames(B):
    """镜面内的局部直角基：e_line 沿 −b1（与 c2 路径上 Γ→M 段内 #11 处的前进方向一致，也与 c3 相同），e_perp = z × e_line，z。"""
    eline = -B[0] / np.linalg.norm(B[0]); z = np.array([0, 0, 1.0])
    return eline, np.cross(z, eline), z


def write_kpoints(path, title, pts):
    with open(path, 'w') as f:
        f.write(title + '\n' + f'{len(pts)}\nReciprocal\n')
        for p in pts: f.write(''.join(f'{x:+.14f}  ' for x in p) + '     1\n')


# ------------------------------------------------------------------ c4：镜面网格
def cmd_grid(a):
    B = recip(a.material)
    us = np.round(np.arange(0, a.umin - 1e-9, -a.du), 10)
    ws = -0.5 + np.arange(a.nw) / a.nw
    pts = [np.array([u, 0.0, w]) for u in us for w in ws]
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    chk = []
    for i in CHECKS:
        chk.append({'path_index': i, 'k_index': len(pts)}); pts.append(K[i])
    os.makedirs(a.out, exist_ok=True)
    write_kpoints(os.path.join(a.out, 'KPOINTS_OPT'),
                  f'c4 mirror-plane grid for {a.material}: (u,0,w), u in [{a.umin},0] step {a.du}, w in [-1/2,1/2) /{a.nw}; see c4_grid.json', pts)
    meta = {'material': a.material, 'nelect': NEL[a.material], 'us': us.tolist(), 'ws': ws.tolist(),
            'n_grid': len(us) * len(ws), 'checks': chk, 'n_kpoints': len(pts), 'B_rows_cart_invA': B.tolist()}
    json.dump(meta, open(os.path.join(a.out, 'c4_grid.json'), 'w'), indent=1, ensure_ascii=False)
    nb1, cz = np.linalg.norm(B[0]), B[2, 2]
    print(f'{len(pts)} 个 k 点（网格 {len(us)}×{len(ws)} + 校验 {len(chk)}）→ {a.out}/KPOINTS_OPT；'
          f'网格间距：沿线 {a.du * nb1:.4f} Å⁻¹，沿 c {cz / a.nw:.4f} Å⁻¹')


# ------------------------------------------------------------------ c4 分析
def quad_fit(xs, zs, es):
    """E = c0 + gx x + gz z + ½Hxx x² + ½Hzz z² + Hxz x z（最小二乘）。返回 c0, g(2), H(2×2)。"""
    A = np.column_stack([np.ones_like(xs), xs, zs, 0.5 * xs ** 2, 0.5 * zs ** 2, xs * zs])
    c, *_ = np.linalg.lstsq(A, es, rcond=None)
    return c[0], np.array([c[1], c[2]]), np.array([[c[3], c[5]], [c[5], c[4]]])


def extrema(G, us, ws, kind, nb1, cz):
    """G[iu, jw]：网格上一条带的能量。u = 0 一侧的邻居用 (u, w) ≡ (−u, −w) 补；u = umin 一侧没有邻居（记为边界）。"""
    nu, nw = G.shape
    du = us[0] - us[1]
    jneg = lambda j: (nw - j) % nw                # w_j = −1/2 + j/nw，−w_j 对应序号 nw − j（模 nw）

    def val(i, j):
        j %= nw
        if i < 0: return G[-i, jneg(j)]          # u > 0：用 (−u, −w)
        if i >= nu: return None
        return G[i, j]

    sgn = 1 if kind == 'max' else -1
    out = []
    for i in range(nu):
        for j in range(nw):
            nbr = [val(i + di, j + dj) for di in (-1, 0, 1) for dj in (-1, 0, 1) if (di, dj) != (0, 0)]
            if any(sgn * (v - G[i, j]) > 0 for v in nbr if v is not None): continue
            pts = [(di, dj, val(i + di, j + dj)) for di in (-1, 0, 1) for dj in (-1, 0, 1)]
            pts = [p for p in pts if p[2] is not None]
            # 局部坐标（Å⁻¹）：网格序号 i 每加 1，u 减小 du；x 沿 e_line = −b1，所以 x = di·du·|b1|
            xs = np.array([p[0] * du * nb1 for p in pts])
            zs = np.array([p[1] / nw * cz for p in pts]); es = np.array([p[2] for p in pts])
            c0, g, H = quad_fit(xs, zs, es)
            ok = len(pts) == 9 and np.all(np.linalg.eigvalsh(H) * sgn < 0)
            s = -np.linalg.solve(H, g) if ok else np.zeros(2)
            inside = ok and abs(s[0]) <= du * nb1 and abs(s[1]) <= cz / nw
            E = c0 + 0.5 * g @ s if inside else G[i, j]
            if not inside: s = np.zeros(2)
            u = us[i] - s[0] / nb1; w = ws[j] + s[1] / cz
            w = (w + 0.5) % 1 - 0.5
            out.append({'grid': [i, j], 'u': float(u), 'w': float(w), 'E': float(E), 'E_grid': float(G[i, j]),
                        'refined': bool(inside), 'boundary': i == nu - 1,
                        'H_line_c_eVA2': H.tolist(), 'm_line_c_me': [float(-sgn * HB / H[0, 0]) if H[0, 0] else None,
                                                                    float(-sgn * HB / H[1, 1]) if H[1, 1] else None]})
    return sorted(out, key=lambda e: -sgn * e['E'])


def star_points(B, center, d, dz):
    eline, eperp, z = frames(B)
    to_frac = lambda v: v @ np.linalg.inv(B)
    return [center + to_frac(dl * d * eline + dp * d * eperp + dzz * dz * z) for _, dl, dp, dzz in STENCIL]


def cmd_analyze(a):
    meta = json.load(open(os.path.join(a.dir, 'c4_grid.json')))
    m, nel = meta['material'], meta['nelect']
    B = np.array(meta['B_rows_cart_invA']); nb1, cz = np.linalg.norm(B[0]), B[2, 2]
    E, fs = read_procar(find_procar(a.dir))
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"PROCAR 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    P = read_path_eigen(os.path.join(mat_dir(m), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], nel - 5:nel + 4] - P[c['path_index'], nel - 5:nel + 4]).max() for c in meta['checks'])
    print(f'# {m}：校验点与 c2 路径（带边 ±4 条带）最大差 {1000 * dev:.3f} meV' + ('' if dev < 5e-4 else '  ← 超过 0.5 meV，先查输入是否与 c2 相同'))
    us, ws = np.array(meta['us']), np.array(meta['ws']); nu, nw = len(us), len(ws)
    vb_path, cb_path = np.nanmax(P[:, nel - 1]), np.nanmin(P[:, nel])
    res = {}
    for kind, b, ref in (('max', nel - 1, vb_path), ('min', nel, cb_path)):
        G = E[:meta['n_grid'], b].reshape(nu, nw); F = fs[:meta['n_grid'], b].reshape(nu, nw)
        ex = extrema(G, us, ws, kind, nb1, cz)
        for e in ex: e['cation_s'] = float(F[e['grid'][0], e['grid'][1]])
        res[kind] = ex
    vbm = max(vb_path, res['max'][0]['E']); cbm = min(cb_path, res['min'][0]['E'])
    print(f'# 路径：价带顶 {vb_path:.5f}、导带底 {cb_path:.5f} eV，路径带隙 {1000 * (cb_path - vb_path):.1f} meV')
    print(f'# 镜面网格最高价带极大 {res["max"][0]["E"]:.5f} eV（比路径价带顶 {1000 * (res["max"][0]["E"] - vb_path):+.1f} meV）；'
          f'最低导带极小 {res["min"][0]["E"]:.5f} eV（比路径导带底 {1000 * (res["min"][0]["E"] - cb_path):+.1f} meV）')
    print(f'# 真实价带顶 {vbm:.5f}、导带底 {cbm:.5f} eV → 带隙 {1000 * (cbm - vbm):.1f} meV')
    for kind, ref, lab in (('max', vbm, '价带极大'), ('min', cbm, '导带极小')):
        print(f'\n{lab}（按能量排序；E 相对{"价带顶" if kind == "max" else "导带底"}）：')
        for e in res[kind][:a.show]:
            mm = e['m_line_c_me']
            print(f"  (u, w) = ({e['u']:+.4f}, {e['w']:+.4f})  E {1000 * (e['E'] - ref):+7.1f} meV  阳离子 s {100 * e['cation_s']:4.1f}%  "
                  f"质量（线, c）≈ {mm[0]:.3f}, {mm[1]:.3f}  {'二次拟合' if e['refined'] else '网格点（未细化）'}{'；在 umin 边界上' if e['boundary'] else ''}")
    out = {'material': m, 'nelect': nel, 'check_dev_meV': 1000 * dev, 'path_vbm': vb_path, 'path_cbm': cb_path,
           'vbm': vbm, 'cbm': cbm, 'gap_meV': 1000 * (cbm - vbm), 'vb_maxima': res['max'], 'cb_minima': res['min']}
    json.dump(out, open(os.path.join(a.dir, 'c4_extrema.json'), 'w'), indent=1, ensure_ascii=False)
    # c5：星形
    sel = [('v', e) for e in res['max'] if vbm - e['E'] <= a.win_v and not e['boundary']]
    sel += [('c', e) for e in res['min'] if e['E'] - cbm <= a.win_c and not e['boundary']]
    # 只给经二次拟合确认的极值做星形；u = 0 是 Γ–A 线（三重轴），那里的谷已由对称性定好
    sel = [(k, e) for k, e in sel if e['refined'] and abs(e['u']) > 1e-6]
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    pts, vals = [], []
    for n, (kind, e) in enumerate(sel):
        c = np.array([e['u'], 0.0, e['w']])
        vals.append({'label': f'{kind}{n}', 'kind': kind, 'band': nel if kind == 'v' else nel + 1, 'center_frac': c.tolist(),
                     'E_c4_eV': e['E'], 'first_k': len(pts), 'delta': a.delta, 'delta_z': a.delta_z, 'stencil': [s[0] for s in STENCIL]})
        pts += star_points(B, c, a.delta, a.delta_z)
    chk = []
    for i in CHECKS:
        chk.append({'path_index': i, 'k_index': len(pts)}); pts.append(K[i])
    sd = os.path.join(os.path.dirname(os.path.abspath(a.dir)), 'c5_star')
    os.makedirs(sd, exist_ok=True)
    write_kpoints(os.path.join(sd, 'KPOINTS_OPT'), f'c5 stars at the c4 mirror-plane extrema of {m}; see c5_points.json', pts)
    eline, eperp, _ = frames(B)
    json.dump({'material': m, 'nelect': nel, 'valleys': vals, 'checks': chk, 'n_kpoints': len(pts), 'vbm_c4': vbm, 'cbm_c4': cbm,
               'path_vbm': vb_path, 'path_cbm': cb_path,
               'e_line': eline.tolist(), 'e_perp': eperp.tolist(), 'B_rows_cart_invA': B.tolist()},
              open(os.path.join(sd, 'c5_points.json'), 'w'), indent=1, ensure_ascii=False)
    print(f'\n# c5：{len(sel)} 个谷（价带 {sum(k == "v" for k, _ in sel)}、导带 {sum(k == "c" for k, _ in sel)}）× 11 点 + {len(chk)} 个校验点 = {len(pts)} 个 k 点 → {sd}/KPOINTS_OPT')


# ------------------------------------------------------------------ c5 拟合
def cmd_fit(a):
    meta = json.load(open(os.path.join(a.dir, 'c5_points.json')))
    m, nel = meta['material'], meta['nelect']
    E, fs = read_procar(find_procar(a.dir))
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"PROCAR 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    P = read_path_eigen(os.path.join(mat_dir(m), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], nel - 5:nel + 4] - P[c['path_index'], nel - 5:nel + 4]).max() for c in meta['checks'])
    print(f'# {m}：校验点与 c2 路径最大差 {1000 * dev:.3f} meV')
    out = {}
    for v in meta['valleys']:
        b = v['band'] - 1; sg = -1 if v['kind'] == 'v' else 1
        e = {s: E[v['first_k'] + i, b] for i, s in enumerate(v['stencil'])}
        d, dz = v['delta'], v['delta_z']
        H = np.zeros((3, 3))
        H[0, 0] = (e['l+'] + e['l-'] - 2 * e['0']) / d ** 2
        H[1, 1] = (e['p+'] + e['p-'] - 2 * e['0']) / d ** 2
        H[2, 2] = (e['z+'] + e['z-'] - 2 * e['0']) / dz ** 2
        H[0, 2] = H[2, 0] = (e['l+z+'] + e['l-z-'] - e['l+z-'] - e['l-z+']) / (4 * d * dz)
        g = np.array([(e['l+'] - e['l-']) / (2 * d), (e['p+'] - e['p-']) / (2 * d), (e['z+'] - e['z-']) / (2 * dz)])
        w = np.linalg.eigvalsh(H); mp = sg * HB / w
        extremum = np.all(np.sign(w) == sg)
        s = -np.linalg.solve(H, g); dE = -0.5 * g @ np.linalg.solve(H, g)
        inside = abs(s[0]) <= d and abs(s[1]) <= d and abs(s[2]) <= dz
        e0 = e['0'] + dE if (extremum and inside) else e['0']
        ml, mperp = sg * HB / H[0, 0], sg * HB / H[1, 1]
        out[v['label']] = {'kind': v['kind'], 'center_frac': v['center_frac'], 'E0': float(e['0']), 'E_extremum': float(e0),
                           'cation_s': float(fs[v['first_k'], b]),
                           'm': [float(ml), float(mperp), float(np.prod(mp) / (ml * mperp))], 'principal': mp.tolist(),
                           'H_eVA2_line_perp_c': H.tolist(), 'gradient': g.tolist(), 'stationary_shift_invA': s.tolist(),
                           'true_extremum_within_star': bool(extremum and inside)}
        print(f"{v['label']}（{'价带' if v['kind'] == 'v' else '导带'}，中心 {np.round(v['center_frac'], 4).tolist()}）："
              f"极值能量 {e0:.5f} eV，阳离子 s {100 * fs[v['first_k'], b]:.1f}%")
        print(f"   梯度（线、面内垂直、c；eV·Å）{np.round(g, 4).tolist()}；主质量 {np.round(mp, 3).tolist()}；"
              f"驻点位移 {np.round(s, 4).tolist()} Å⁻¹ → {'在采样半径内，是真极值' if extremum and inside else '不满足（' + ('鞍点' if not extremum else '驻点超出采样半径') + '）'}")
    good = lambda k: [o['E_extremum'] for o in out.values() if o['kind'] == k and o['true_extremum_within_star']]
    vbm = max([meta['path_vbm']] + good('v')); cbm = min([meta['path_cbm']] + good('c'))
    for o in out.values():
        o['E_rel_edge_meV'] = float(1000 * ((vbm - o['E_extremum']) if o['kind'] == 'v' else (o['E_extremum'] - cbm)))
    print(f"# 真实价带顶 {vbm:.5f} eV（路径 {meta['path_vbm']:.5f}，高 {1000 * (vbm - meta['path_vbm']):.1f} meV）；"
          f"导带底 {cbm:.5f} eV（路径 {meta['path_cbm']:.5f}）；带隙 {1000 * (cbm - vbm):.1f} meV（路径 {1000 * (meta['path_cbm'] - meta['path_vbm']):.1f} meV）")
    for k, o in out.items():
        print(f"   {k}：离{'价带顶' if o['kind'] == 'v' else '导带底'} {o['E_rel_edge_meV']:.1f} meV" + ('' if o['true_extremum_within_star'] else '（未确认为真极值）'))
    json.dump({'vbm': vbm, 'cbm': cbm, 'gap_meV': 1000 * (cbm - vbm), 'path_vbm': meta['path_vbm'], 'path_cbm': meta['path_cbm'], 'valleys': out},
              open(os.path.join(a.dir, 'c5_masses.json'), 'w'), indent=1, ensure_ascii=False)
    print(f"# 写出 {os.path.join(a.dir, 'c5_masses.json')}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('grid'); s.add_argument('material'); s.add_argument('--out', required=True)
    s.add_argument('--umin', type=float, default=-0.40); s.add_argument('--du', type=float, default=0.025)
    s.add_argument('--nw', type=int, default=12)
    s = sub.add_parser('analyze'); s.add_argument('dir')
    s.add_argument('--win-v', type=float, default=0.060, help='价带极大离价带顶多少 eV 以内做星形')
    s.add_argument('--win-c', type=float, default=0.100, help='导带极小离导带底多少 eV 以内做星形')
    s.add_argument('--delta', type=float, default=0.02); s.add_argument('--delta-z', type=float, default=0.03)
    s.add_argument('--show', type=int, default=8)
    s = sub.add_parser('fit'); s.add_argument('dir')
    a = ap.parse_args()
    {'grid': cmd_grid, 'analyze': cmd_analyze, 'fit': cmd_fit}[a.cmd](a)


if __name__ == '__main__':
    main()
