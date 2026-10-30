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
  local     c6：c4 中未能细化的价带极大（网格二次拟合的驻点超出一格），在它周围另取一块细网格：u 步长 --du（默认 1/60，
            12u = n/5 永远不是半整数，避开 FOCK_FCC 的奇数/24 点）、w 步长 --dw（默认 1/36），各 ±--nu、±--nw 格；
            另在中心加垂直镜面的 ±δ 两点（面内垂直方向的曲率；镜面对称使该方向梯度为零）。
  local-fit 读 c6：每块细网格上找极大，3×3 二次拟合定驻点与（线, c）曲率，加上垂直方向曲率，给出能量与主质量；
            极大落在细网格边界上时报告“未收敛”，不写入。结果 c6_masses.json 与 c5_masses.json 同一格式，
            multivalley_seebeck.true_load 会自动用它替换对应的 c4 网格谷。

KPOINTS_OPT 的排布（2026-10-29，针对 VASP 在 NKRED 2 2 1 下的 FOCK_FCC 内部错误“number of k-points incorrect”）：
  原网格（全部 ky = 0，208 点）第 1 批通过、第 2 批崩溃；第 1 批含落在 12×12×2 自洽网格上的点（Γ、A），第 2 批没有。
  为同时排除“批内点共面”和“批内没有自洽网格点”两种可能，现在：
  - 网格点轮流写成三个等价表示 (u,0,w)、(−u,u,w)、(0,−u,w)（绕 c 轴转 0、120°、240°，能量严格相同），批内不再共面；
  - 每批（KPOINTS_OPT_NKBATCH 个点，默认 36）的第 1 个点是一个“锚点”：c2 路径上同时在自洽网格与镜面内的点
    （#0 Γ、#11、#22、#33 M、#97 A、#108、#119、#130 L），轮流使用；它们也是与 c2 对比的校验点；
  - 另加 4 个对称校验点（两个网格点的另两个等价表示），直接检验三种表示的能量相同。
  INCAR 必须设 KPOINTS_OPT_NKBATCH 等于 --nkbatch，批的边界才与锚点对齐。

  2026-10-29 晚，根因（据 VASP 6.6.0 源码，见 write_kpoints 前的说明）：崩溃与批的排法无关，只取决于有没有分数坐标为
  奇数/24 的点（NKRED = 2 的 x、y 方向）。原网格步长 0.025 时 u = −0.125、−0.375 正是这种点。现在 write_kpoints 遇到这种点直接报错；
  网格改用步长 1/36（u = 0 … −15/36，含 −1/6、−1/3 两个 c2 谷所在的 u），NKRED 2 2 1 不变。
  另有三个选项：--split 把每批另写成 KPOINTS_OPT_b01 …，可各开一个作业并行算；analyze、fit 用 --parts 按批号顺序读各作业的
  PROCAR_OPT；--wgroup、--c5-group 按 w 或按谷排批（当时针对“批内 kz 取值数”的猜测，已不需要，保留备用）。

用法（Sn2Bi2Te5）：
  python 10_脚本/c4_mirror_scan.py grid Sn2Bi2Te5 --out 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_v3 --du 1/36 --umin=-15/36 --split
  python 10_脚本/c4_mirror_scan.py analyze 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_v3 --parts <b01…b06 的目录> --split   # 同时写出 ../c5_star/
  python 10_脚本/c4_mirror_scan.py fit 22_能带基准_1-0/Sn2Bi2Te5/c5_star --parts <c5 各批目录>
  python 10_脚本/c4_mirror_scan.py local 22_能带基准_1-0/Sn2Bi2Te5/c4_mirror_v3 --out 22_能带基准_1-0/Sn2Bi2Te5/c6_local --split   # c6
  python 10_脚本/c4_mirror_scan.py local-fit 22_能带基准_1-0/Sn2Bi2Te5/c6_local --parts <c6 各批目录>
Sn2Sb2Te5（孤对电子谷在价带顶下 130–150 meV，analyze 要放宽窗口）：
  python 10_脚本/c4_mirror_scan.py grid Sn2Sb2Te5 --out 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3 --du 1/36 --umin=-15/36 --split
  python 10_脚本/c4_mirror_scan.py analyze 22_能带基准_1-0/Sn2Sb2Te5/c4_mirror_v3 --parts <b01…b06> --split --win-v 0.20
单个作业顺序跑全部批时，去掉 --parts，PROCAR_OPT 放在目录里即可。
只依赖 numpy。
"""
import argparse, gzip, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HB = 7.62                                     # ħ²/m_e，eV·Å²
NEL = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
CHECKS = (0, 11, 101, 119)                    # c2 路径校验点：Γ 与三个镜面内的离轴谷
ANCHORS = (0, 11, 22, 33, 97, 108, 119, 130)  # c2 路径上同时落在 12×12×2 自洽网格与镜面内的点，用作每批的锚点
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


# FOCK_FCC 的判据（VASP 6.6.0 源码 fock_dbl.F:1407、fock.F 的 SKIP_THIS_KPOINT_IN_FOCK）：对每个 KPOINTS_OPT 点 k，
# 在 12×12×2 全网格上数 q = k_full − k 中 MOD(FLOOR((q+32)·NKP + 0.5), NKRED) = 0 的个数，必须恰好 288/(NKREDX·NKREDY·NKREDZ)。
# 每个方向上 12 个网格坐标减去 k 再乘 12，是 12 个间隔为 1 的数，FLOOR 后模 2 恰好留下一半，与 k 无关；
# 唯一的例外是 12·k 恰为半整数（k = 奇数/24）：这 12 个数全都正好落在整数上，FLOOR 取哪边由浮点舍入决定，个数就不再是一半。
# c4 原网格（步长 0.025）中 u = −0.125、−0.375 正是这种点，崩溃的批都含它们，通过的批都不含。
NKP_SCF, NKRED = (12, 12, 2), (2, 2, 1)


def fock_unsafe(pts, tol=1e-3):
    """返回会触发 FOCK_FCC 判据边界的点：NKRED ≥ 2 的方向上 NKP·k 离半整数不到 tol。"""
    bad = []
    for n, p in enumerate(pts):
        for d in range(3):
            if NKRED[d] >= 2 and abs((NKP_SCF[d] * p[d]) % 1 - 0.5) < tol: bad.append((n, d, p[d])); break
    return bad


def write_kpoints(path, title, pts):
    bad = fock_unsafe(pts)
    if bad: sys.exit(f'{path}：{len(bad)} 个点的分数坐标是奇数/24（NKRED = 2 的方向），会触发 FOCK_FCC 内部错误，'
                     f'例如第 {bad[0][0] + 1} 点 {np.round(pts[bad[0][0]], 6).tolist()}。换网格步长或星形步长')
    with open(path, 'w') as f:
        f.write(title + '\n' + f'{len(pts)}\nReciprocal\n')
        for p in pts: f.write(''.join(f'{x:+.14f}  ' for x in p) + '     1\n')


# ------------------------------------------------------------------ c4：镜面网格
def c3_images(f):
    """绕 c 轴转 0、120°、240° 的等价点（六方 120° 设置下倒格矢变换 b1 → −b1 + b2 → −b2）：(u,v,w) → (−u−v, u, w) → (v, −u−v, w)。"""
    u, v, w = f
    return [np.array([u, v, w]), np.array([-u - v, u, w]), np.array([v, -u - v, w])]


def layout(items, K, nkbatch):
    """每 nkbatch 个点为一批，每批第 1 个是锚点（轮流取 ANCHORS），其后是 nkbatch − 1 个 items。返回 k 点列表与逐点记录。"""
    pts, recs, n, a = [], [], 0, 0
    while n < len(items):
        i = ANCHORS[a % len(ANCHORS)]; a += 1
        recs.append({'kind': 'anchor', 'path_index': i, 'k_index': len(pts)}); pts.append(np.array(K[i]))
        for it in items[n:n + nkbatch - 1]:
            it = dict(it); k = it.pop('k'); it['k_index'] = len(pts); pts.append(np.array(k)); recs.append(it)
        n += nkbatch - 1
    return pts, recs


def fold_kz(w):
    return round((w + 0.5) % 1 - 0.5, 6)


def layout_groups(groups, extras, K, nkbatch, pad):
    """按组排批：每组一批，批首是锚点（优先取 kz 与该组相同的），其后是该组的点；extras 放进有空位、且 kz 相同的批（没有就放进任一有空位的批）。
    除最后一批外，每批用 pad(批内的点, 缺几个) 补足 nkbatch 个点，这样单个作业按 KPOINTS_OPT_NKBATCH 切出的批与组一一对应。
    返回 k 点列表、逐点记录（含 batch）与每批的 [起, 止) 序号。"""
    anc = {}
    for i in ANCHORS: anc.setdefault(fold_kz(K[i][2]), []).append(i)
    bat, na = [], 0
    for g in groups:
        if len(g) + 1 > nkbatch: sys.exit(f'一组 {len(g)} 个点加锚点超过每批 {nkbatch} 个')
        pool = next((anc[z] for z in sorted({fold_kz(it['k'][2]) for it in g}) if z in anc), ANCHORS)
        i = pool[na % len(pool)]; na += 1
        bat.append([{'kind': 'anchor', 'path_index': i, 'k': list(K[i])}] + list(g))
    for it in extras:
        free = [b for b in bat if len(b) < nkbatch]
        if not free: sys.exit('extras 放不下，增大 --nkbatch')
        z = fold_kz(it['k'][2])
        next((b for b in free if z in {fold_kz(x['k'][2]) for x in b}), free[0]).append(it)
    bat.sort(key=lambda b: -len(b))                    # 不满的批放到最后
    for b in bat[:-1]:
        b += pad(b, nkbatch - len(b))
        if len(b) != nkbatch: sys.exit('补点后批的大小仍不等于 nkbatch')
    pts, recs, batches = [], [], []
    for n, b in enumerate(bat):
        s = len(pts)
        for it in b:
            it = dict(it); k = it.pop('k'); it['k_index'] = len(pts); it['batch'] = n; pts.append(np.array(k)); recs.append(it)
        batches.append([s, len(pts)])
    return pts, recs, batches


def batch_summary(pts, batches):
    """每批的点数、kz 取值数、Γ–A 线（kx = ky = 0）上的点数，供对照 FOCK_FCC 崩溃的经验规律。"""
    out = []
    for s, e in batches:
        P = np.array(pts[s:e]); f = P - np.round(P)
        out.append({'k_index': [s, e], 'n': e - s, 'n_kz': len({fold_kz(p[2]) for p in P}),
                    'n_Gamma_A': int(np.sum(np.all(np.isclose(f[:, :2], 0, atol=1e-8), axis=1)))})
    return out


def write_batches(d, stem, title, pts, batches):
    """每批另写一个 KPOINTS_OPT（<stem>_b01 …），可各开一个作业并行算；结果按批号顺序交给 analyze/fit 的 --parts。"""
    for n, (s, e) in enumerate(batches):
        write_kpoints(os.path.join(d, f'{stem}_b{n + 1:02d}'), f'{title} | batch {n + 1}/{len(batches)}, global k_index {s}-{e - 1}', pts[s:e])


def read_parts(d, parts, meta):
    """parts 为空时读 d 里的 PROCAR_OPT；否则按顺序读各批目录的 PROCAR_OPT 并接起来，逐批核对点数。"""
    if not parts: return read_procar(find_procar(d))
    got = [read_procar(find_procar(p)) for p in parts]
    want = [e - s for s, e in meta.get('batches', [])]
    if want and [E.shape[0] for E, _ in got] != want: sys.exit(f'各批 PROCAR 的 k 点数 {[E.shape[0] for E, _ in got]} 与记录 {want} 不符')
    nb = min(E.shape[1] for E, _ in got)
    return np.vstack([E[:, :nb] for E, _ in got]), np.vstack([f[:, :nb] for _, f in got])


def cmd_grid(a):
    B = recip(a.material)
    us = np.round(np.arange(0, a.umin - 1e-9, -a.du), 10)
    ws = -0.5 + np.arange(a.nw) / a.nw
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    items, n = [], 0
    for iu, u in enumerate(us):
        for jw, w in enumerate(ws):
            r = n % 3 if a.rotate else 0; n += 1
            items.append({'kind': 'grid', 'grid': [iu, jw], 'rep': r, 'k': c3_images((u, 0.0, w))[r].tolist()})
    check = {'kind': 'check', 'path_index': 101, 'k': K[101].tolist()}               # 不在自洽网格上的 c2 校验点
    symcheck = lambda iu, jw, r: {'kind': 'symcheck', 'grid': [iu, jw], 'rep': r, 'k': c3_images((us[iu], 0.0, ws[jw]))[r].tolist()}
    batches = None
    if a.wgroup:
        # 按 w 分组：每批只含 a.wgroup 个 w（间隔 nw/wgroup，第一组是 w = −1/2 与 0 两个自洽网格面）× 全部 u，批内 kz 取值少
        if a.nw % a.wgroup: sys.exit('--nw 必须能被 --wgroup 整除')
        st = a.nw // a.wgroup
        groups = [[it for it in items if it['grid'][1] % st == j] for j in range(st)]

        def pad(b, n):                                 # 用批内网格点的另一个旋转表示补位（kz 不变），同时就是对称校验
            g = [it for it in b if it['kind'] == 'grid']
            mid = sorted(g, key=lambda it: (abs(us[it['grid'][0]] - a.umin / 2), it['grid'][1]))
            return [symcheck(*it['grid'], (it['rep'] + 1) % 3) for it in mid[:n]]
        pts, recs, batches = layout_groups(groups, [check], K, a.nkbatch, pad)
    else:
        items.append(check)
        for u0, w0 in ((-0.2, 0.25), (-0.35, -1 / 3)):
            iu, jw = int(np.argmin(abs(us - u0))), int(np.argmin(abs(ws - w0)))
            used = next(it['rep'] for it in items if it['kind'] == 'grid' and it['grid'] == [iu, jw])
            items += [symcheck(iu, jw, r) for r in range(3) if r != used]
        pts, recs = layout(items, K, a.nkbatch)
    os.makedirs(a.out, exist_ok=True)
    title = (f'c4 mirror-plane grid for {a.material}: u in [{a.umin},0] step {a.du}, w in [-1/2,1/2) /{a.nw}; '
             + (f'{"C3-rotated reps" if a.rotate else "(u,0,w) reps"}, batches grouped by w ({a.wgroup} w values each), '
                if a.wgroup else 'C3-rotated reps, ')
             + f'one c2 anchor per batch of {a.nkbatch}; set KPOINTS_OPT_NKBATCH = {a.nkbatch}; see c4_grid.json')
    write_kpoints(os.path.join(a.out, 'KPOINTS_OPT'), title, pts)
    if batches is None: batches = [[s, min(s + a.nkbatch, len(pts))] for s in range(0, len(pts), a.nkbatch)]
    if a.split: write_batches(a.out, 'KPOINTS_OPT', title, pts, batches)
    gk = [[None] * len(ws) for _ in us]
    for r_ in recs:
        if r_['kind'] == 'grid': gk[r_['grid'][0]][r_['grid'][1]] = r_['k_index']
    meta = {'material': a.material, 'nelect': NEL[a.material], 'us': us.tolist(), 'ws': ws.tolist(), 'nkbatch': a.nkbatch,
            'n_grid': len(us) * len(ws), 'grid_k_index': gk,
            'checks': [{'path_index': r_['path_index'], 'k_index': r_['k_index']} for r_ in recs if r_['kind'] in ('anchor', 'check')],
            'symchecks': [{'grid': r_['grid'], 'rep': r_['rep'], 'k_index': r_['k_index']} for r_ in recs if r_['kind'] == 'symcheck'],
            'records': recs, 'n_kpoints': len(pts), 'B_rows_cart_invA': B.tolist(),
            'batches': batches, 'batch_summary': batch_summary(pts, batches)}
    json.dump(meta, open(os.path.join(a.out, 'c4_grid.json'), 'w'), indent=1, ensure_ascii=False)
    nb1, cz = np.linalg.norm(B[0]), B[2, 2]
    na = sum(r_['kind'] == 'anchor' for r_ in recs)
    print(f'{len(pts)} 个 k 点 = 网格 {len(us)}×{len(ws)} + 锚点 {na}（每批 1 个）+ 校验 1 + 对称校验 {len(meta["symchecks"])}'
          f' → {a.out}/KPOINTS_OPT；INCAR 设 KPOINTS_OPT_NKBATCH = {a.nkbatch}（共 {len(batches)} 批）')
    print(f'网格间距：沿线 {a.du * nb1:.4f} Å⁻¹，沿 c {cz / a.nw:.4f} Å⁻¹')
    for n, s in enumerate(meta['batch_summary']):
        print(f'  批 {n + 1}：k_index {s["k_index"][0]}–{s["k_index"][1] - 1}，{s["n"]} 点，kz 取值 {s["n_kz"]} 个，Γ–A 线上 {s["n_Gamma_A"]} 点'
              + (f'  → KPOINTS_OPT_b{n + 1:02d}' if a.split else ''))


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
    E, fs = read_parts(a.dir, a.parts, meta)
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"PROCAR 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    P = read_path_eigen(os.path.join(mat_dir(m), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], nel - 5:nel + 4] - P[c['path_index'], nel - 5:nel + 4]).max() for c in meta['checks'])
    print(f'# {m}：{len(meta["checks"])} 个校验点与 c2 路径（带边 ±4 条带）最大差 {1000 * dev:.3f} meV'
          + ('' if dev < 5e-4 else ('  ← 超过 0.5 meV' + ('；按 --internal-ref 用本次计算的校验点作参考' if a.internal_ref else '，先查输入是否与 c2 相同'))))
    us, ws = np.array(meta['us']), np.array(meta['ws']); nu, nw = len(us), len(ws)
    gk = np.array(meta.get('grid_k_index') or np.arange(meta['n_grid']).reshape(nu, nw))
    if meta.get('symchecks'):
        sdev = max(abs(E[c['k_index'], nel - 5:nel + 4] - E[gk[c['grid'][0], c['grid'][1]], nel - 5:nel + 4]).max() for c in meta['symchecks'])
        print(f'# 对称校验（同一点的不同旋转表示）最大差 {1000 * sdev:.4f} meV' + ('' if sdev < 1e-4 else '  ← 超过 0.1 meV，等价表示不成立，停下检查'))
    if a.internal_ref:      # NKRED 与 c2 不同时：带边参考取本次计算中的 c2 校验点（SnBi 的路径价带顶在 #11、导带底在 Γ）
        vb_path = max(E[c['k_index'], nel - 1] for c in meta['checks']); cb_path = min(E[c['k_index'], nel] for c in meta['checks'])
    else:
        vb_path, cb_path = np.nanmax(P[:, nel - 1]), np.nanmin(P[:, nel])
    res = {}
    for kind, b, ref in (('max', nel - 1, vb_path), ('min', nel, cb_path)):
        G = E[gk, b]; F = fs[gk, b]
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
    items, vals = [], []
    for n, (kind, e) in enumerate(sel):
        c = np.array([e['u'], 0.0, e['w']]); lab = f'{kind}{n}'
        vals.append({'label': lab, 'kind': kind, 'band': nel if kind == 'v' else nel + 1, 'center_frac': c.tolist(),
                     'E_c4_eV': e['E'], 'delta': a.delta, 'delta_z': a.delta_z})
        for (st, *_), kp in zip(STENCIL, star_points(B, c, a.delta, a.delta_z)):
            items.append({'kind': 'star', 'label': lab, 'stencil': st, 'k': kp.tolist()})
    check = {'kind': 'check', 'path_index': 101, 'k': K[101].tolist()}
    nk5, batches = a.nkbatch, None
    if a.c5_group and vals:
        # 每个谷一批（锚点 + 11 点星形），批内 kz 只有 3–4 个；空位用星形中心的另两个旋转表示补（kz 不变，fit 不用它们）
        groups = [[it for it in items if it['label'] == v['label']] for v in vals]
        nk5 = max(len(g) for g in groups) + 2

        def pad(b, n):
            c = next(it for it in b if it.get('stencil') == '0')
            return [{'kind': 'symcheck', 'label': c['label'], 'rep': r, 'k': c3_images(c['k'])[r].tolist()} for r in (1, 2)][:n]
        pts, recs, batches = layout_groups(groups, [check], K, nk5, pad)
    else:
        pts, recs = layout(items + [check], K, nk5)
    if batches is None: batches = [[s, min(s + nk5, len(pts))] for s in range(0, len(pts), nk5)]
    for v in vals:
        v['k_indices'] = {r_['stencil']: r_['k_index'] for r_ in recs if r_['kind'] == 'star' and r_['label'] == v['label']}
    chk = [{'path_index': r_['path_index'], 'k_index': r_['k_index']} for r_ in recs if r_['kind'] in ('anchor', 'check')]
    sd = os.path.join(os.path.dirname(os.path.abspath(a.dir)), 'c5_star')
    os.makedirs(sd, exist_ok=True)
    title = (f'c5 stars at the c4 mirror-plane extrema of {m}; ' + ('one valley per batch, ' if a.c5_group else '')
             + f'one c2 anchor per batch of {nk5}; set KPOINTS_OPT_NKBATCH = {nk5}; see c5_points.json')
    write_kpoints(os.path.join(sd, 'KPOINTS_OPT'), title, pts)
    if a.split: write_batches(sd, 'KPOINTS_OPT', title, pts, batches)
    eline, eperp, _ = frames(B)
    json.dump({'material': m, 'nelect': nel, 'valleys': vals, 'checks': chk, 'n_kpoints': len(pts), 'nkbatch': nk5,
               'batches': batches, 'batch_summary': batch_summary(pts, batches),
               'vbm_c4': vbm, 'cbm_c4': cbm, 'path_vbm': vb_path, 'path_cbm': cb_path, 'internal_ref': bool(a.internal_ref),
               'e_line': eline.tolist(), 'e_perp': eperp.tolist(), 'B_rows_cart_invA': B.tolist(), 'records': recs},
              open(os.path.join(sd, 'c5_points.json'), 'w'), indent=1, ensure_ascii=False)
    print(f'\n# c5：{len(sel)} 个谷（价带 {sum(k == "v" for k, _ in sel)}、导带 {sum(k == "c" for k, _ in sel)}）× 11 点 + 锚点与校验 {len(chk)} 个'
          f' = {len(pts)} 个 k 点 → {sd}/KPOINTS_OPT；INCAR 设 KPOINTS_OPT_NKBATCH = {nk5}（共 {len(batches)} 批）')
    for n, s in enumerate(batch_summary(pts, batches)):
        print(f'  批 {n + 1}：{s["n"]} 点，kz 取值 {s["n_kz"]} 个，Γ–A 线上 {s["n_Gamma_A"]} 点' + (f'  → KPOINTS_OPT_b{n + 1:02d}' if a.split else ''))


# ------------------------------------------------------------------ c5 拟合
def cmd_fit(a):
    meta = json.load(open(os.path.join(a.dir, 'c5_points.json')))
    m, nel = meta['material'], meta['nelect']
    E, fs = read_parts(a.dir, a.parts, meta)
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"PROCAR 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    P = read_path_eigen(os.path.join(mat_dir(m), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], nel - 5:nel + 4] - P[c['path_index'], nel - 5:nel + 4]).max() for c in meta['checks'])
    print(f'# {m}：校验点与 c2 路径最大差 {1000 * dev:.3f} meV' + ('（内部参考模式：只看相对能量）' if meta.get('internal_ref') else ''))
    out = {}
    for v in meta['valleys']:
        b = v['band'] - 1; sg = -1 if v['kind'] == 'v' else 1
        e = {st: E[k, b] for st, k in v['k_indices'].items()}; k0 = v['k_indices']['0']
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
                           'cation_s': float(fs[k0, b]),
                           'm': [float(ml), float(mperp), float(np.prod(mp) / (ml * mperp))], 'principal': mp.tolist(),
                           'H_eVA2_line_perp_c': H.tolist(), 'gradient': g.tolist(), 'stationary_shift_invA': s.tolist(),
                           'true_extremum_within_star': bool(extremum and inside)}
        print(f"{v['label']}（{'价带' if v['kind'] == 'v' else '导带'}，中心 {np.round(v['center_frac'], 4).tolist()}）："
              f"极值能量 {e0:.5f} eV，阳离子 s {100 * fs[k0, b]:.1f}%")
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


# ------------------------------------------------------------------ c6：未细化极大的局部细网格
def cmd_local(a):
    ex = json.load(open(os.path.join(a.c4dir, 'c4_extrema.json')))
    m, nel = ex['material'], ex['nelect']; B = recip(m); nb1, cz = np.linalg.norm(B[0]), B[2, 2]
    eline, eperp, z = frames(B); to_frac = lambda v: v @ np.linalg.inv(B)
    cen = [e for e in ex['vb_maxima'] if not e['refined'] and abs(e['u']) > 1e-6 and ex['vbm'] - e['E'] <= a.win_v]
    if not cen: sys.exit('c4 中没有需要细化的价带极大')
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    items, vals, n = [], [], 0
    for iv, e in enumerate(cen):
        lab = f'g{iv}'; u0, w0 = e['u'], e['w']
        us = [u0 + i * a.du for i in range(-a.nu, a.nu + 1)]; ws = [w0 + j * a.dw for j in range(-a.nw, a.nw + 1)]
        for i, u in enumerate(us):
            for j, w in enumerate(ws):
                r = n % 3; n += 1
                items.append({'kind': 'grid', 'label': lab, 'grid': [i, j], 'k': c3_images((u, 0.0, (w + 0.5) % 1 - 0.5))[r].tolist()})
        c = np.array([u0, 0.0, w0])
        for st, sg in (('p+', 1), ('p-', -1)):
            items.append({'kind': 'perp', 'label': lab, 'stencil': st, 'k': (c + to_frac(sg * a.delta * eperp)).tolist()})
        vals.append({'label': lab, 'from_c4_uw': [u0, w0], 'E_c4_eV': e['E'], 'us': us, 'ws': ws, 'du': a.du, 'dw': a.dw, 'delta': a.delta})
    check = {'kind': 'check', 'path_index': 101, 'k': K[101].tolist()}
    pts, recs = layout(items + [check], K, a.nkbatch)
    batches = [[s_, min(s_ + a.nkbatch, len(pts))] for s_ in range(0, len(pts), a.nkbatch)]
    for v in vals:
        v['grid_k_index'] = [[None] * len(v['ws']) for _ in v['us']]
        for r_ in recs:
            if r_.get('label') != v['label']: continue
            if r_['kind'] == 'grid': v['grid_k_index'][r_['grid'][0]][r_['grid'][1]] = r_['k_index']
            if r_['kind'] == 'perp': v.setdefault('perp_k_index', {})[r_['stencil']] = r_['k_index']
    os.makedirs(a.out, exist_ok=True)
    title = (f'c6 local mirror-plane grids around unrefined c4 maxima of {m}: du {a.du:.6f}, dw {a.dw:.6f}, +perp; '
             f'one c2 anchor per batch of {a.nkbatch}; set KPOINTS_OPT_NKBATCH = {a.nkbatch}; see c6_points.json')
    write_kpoints(os.path.join(a.out, 'KPOINTS_OPT'), title, pts)
    if a.split: write_batches(a.out, 'KPOINTS_OPT', title, pts, batches)
    chk = [{'path_index': r_['path_index'], 'k_index': r_['k_index']} for r_ in recs if r_['kind'] in ('anchor', 'check')]
    json.dump({'material': m, 'nelect': nel, 'valleys': vals, 'checks': chk, 'n_kpoints': len(pts), 'nkbatch': a.nkbatch, 'batches': batches,
               'batch_summary': batch_summary(pts, batches), 'vbm_c4': ex['vbm'], 'path_vbm': ex['path_vbm'], 'path_cbm': ex['path_cbm'],
               'c5_masses': os.path.relpath(os.path.join(os.path.dirname(os.path.abspath(a.c4dir)), 'c5_star', 'c5_masses.json'), a.out),
               'B_rows_cart_invA': B.tolist(), 'records': recs}, open(os.path.join(a.out, 'c6_points.json'), 'w'), indent=1, ensure_ascii=False)
    print(f'{len(cen)} 个未细化的极大 × ({2 * a.nu + 1}×{2 * a.nw + 1} 网格 + 2 个垂直点) + 锚点与校验 {len(chk)} 个 = {len(pts)} 个 k 点'
          f' → {a.out}/KPOINTS_OPT；INCAR 设 KPOINTS_OPT_NKBATCH = {a.nkbatch}（共 {len(batches)} 批）')
    print(f'细网格范围：沿线 ±{a.nu * a.du * nb1:.3f} Å⁻¹、沿 c ±{a.nw * a.dw * cz:.3f} Å⁻¹（c4 一格：沿线 {nb1 / 36:.3f}、沿 c {cz / 12:.3f} Å⁻¹）')
    for v in vals: print(f"  {v['label']}：c4 中心 (u, w) = ({v['from_c4_uw'][0]:+.4f}, {v['from_c4_uw'][1]:+.4f})，E_c4 = {v['E_c4_eV']:.5f} eV")
    for k, s_ in enumerate(batch_summary(pts, batches)):
        print(f'  批 {k + 1}：{s_["n"]} 点，kz 取值 {s_["n_kz"]} 个' + (f'  → KPOINTS_OPT_b{k + 1:02d}' if a.split else ''))


def cmd_local_fit(a):
    meta = json.load(open(os.path.join(a.dir, 'c6_points.json')))
    m, nel = meta['material'], meta['nelect']; B = np.array(meta['B_rows_cart_invA']); nb1, cz = np.linalg.norm(B[0]), B[2, 2]
    E, fs = read_parts(a.dir, a.parts, meta)
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"PROCAR 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    P = read_path_eigen(os.path.join(mat_dir(m), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], nel - 5:nel + 4] - P[c['path_index'], nel - 5:nel + 4]).max() for c in meta['checks'])
    print(f'# {m}：{len(meta["checks"])} 个校验点与 c2 路径最大差 {1000 * dev:.3f} meV' + ('' if dev < 5e-4 else '  ← 超过 0.5 meV，先查输入'))
    c5 = json.load(open(os.path.normpath(os.path.join(a.dir, meta['c5_masses']))))
    vbm = c5['vbm']; b = nel - 1; out = {}
    for v in meta['valleys']:
        gk = np.array(v['grid_k_index']); G = E[gk, b]; nu, nw = G.shape
        i, j = np.unravel_index(np.argmax(G), G.shape)
        if i in (0, nu - 1) or j in (0, nw - 1):
            print(f"{v['label']}：细网格上的最高点在边界 (i, j) = ({i}, {j})，极大不在范围内，需要扩大网格；不写入"); continue
        pts = [(di, dj, G[i + di, j + dj]) for di in (-1, 0, 1) for dj in (-1, 0, 1)]
        xs = np.array([-p[0] * v['du'] * nb1 for p in pts]); zs = np.array([p[1] * v['dw'] * cz for p in pts]); es = np.array([p[2] for p in pts])
        c0, g, H2 = quad_fit(xs, zs, es)
        neg = np.all(np.linalg.eigvalsh(H2) < 0); sft = -np.linalg.solve(H2, g) if neg else np.zeros(2)
        inside = neg and abs(sft[0]) <= v['du'] * nb1 and abs(sft[1]) <= v['dw'] * cz
        e0 = c0 + 0.5 * g @ sft if inside else G[i, j]
        u = v['us'][i] - sft[0] / nb1; w = (v['ws'][j] + sft[1] / cz + 0.5) % 1 - 0.5
        pk = v['perp_k_index']; ec = E[gk[nu // 2, nw // 2], b]
        hp = (E[pk['p+'], b] + E[pk['p-'], b] - 2 * ec) / v['delta'] ** 2
        H = np.array([[H2[0, 0], 0, H2[0, 1]], [0, hp, 0], [H2[0, 1], 0, H2[1, 1]]]); wv = np.linalg.eigvalsh(H)
        ok = bool(inside and np.all(wv < 0))
        out[v['label']] = {'kind': 'v', 'from_c4_uw': v['from_c4_uw'], 'center_frac': [float(u), 0.0, float(w)], 'E_extremum': float(e0),
                           'cation_s': float(fs[gk[i, j], b]), 'm': [float(HB / -H[0, 0]), float(HB / -H[1, 1]), float(HB / -H[2, 2])],
                           'principal': (HB / -wv).tolist(), 'H_eVA2_line_perp_c': H.tolist(), 'true_extremum_within_star': ok,
                           'E_rel_edge_meV': float(1000 * (vbm - e0))}
        print(f"{v['label']}（c4 中心 {np.round(v['from_c4_uw'], 4).tolist()}）：极大 (u, w) = ({u:+.4f}, {w:+.4f})，E = {e0:.5f} eV，"
              f"离价带顶 {1000 * (vbm - e0):.1f} meV（c4 网格点 {1000 * (vbm - v['E_c4_eV']):.1f}），阳离子 s {100 * fs[gk[i, j], b]:.1f}%")
        print(f"   主质量 {np.round(HB / -wv, 3).tolist()} mₑ；{'确认是极大' if ok else '未确认（驻点超出一格或曲率不全为负）'}")
        if e0 > vbm + 1e-4: print('   ← 比 c5 的价带顶还高：价带顶换了位置，停下汇报')
    json.dump({'vbm': vbm, 'path_vbm': c5['path_vbm'], 'valleys': out}, open(os.path.join(a.dir, 'c6_masses.json'), 'w'), indent=1, ensure_ascii=False)
    print(f"# 写出 {os.path.join(a.dir, 'c6_masses.json')}")


def frac(s):
    from fractions import Fraction
    return float(Fraction(s))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('grid'); s.add_argument('material'); s.add_argument('--out', required=True)
    s.add_argument('--umin', type=frac, default=-0.40); s.add_argument('--du', type=frac, default=0.025, help='可写分数，如 1/36')
    s.add_argument('--nw', type=int, default=12)
    s.add_argument('--nkbatch', type=int, default=36, help='每批点数，INCAR 的 KPOINTS_OPT_NKBATCH 必须相同')
    s.add_argument('--no-rotate', dest='rotate', action='store_false', help='不做三重旋转轮换（全部写成 (u,0,w)）')
    s.add_argument('--wgroup', type=int, default=0, help='按 w 分组排批：每批只含这么多个 w（如 2），批内 kz 取值少；0 为原排布')
    s.add_argument('--split', action='store_true', help='另把每批写成 KPOINTS_OPT_b01 …，可各开一个作业')
    s = sub.add_parser('analyze'); s.add_argument('dir')
    s.add_argument('--parts', nargs='+', help='逐批作业的目录（按批号顺序），代替 dir 里的单个 PROCAR_OPT')
    s.add_argument('--c5-group', action='store_true', help='c5 每个谷一批（锚点 + 11 点星形 + 补位），批内 kz 只有 3–4 个')
    s.add_argument('--split', action='store_true', help='另把 c5 的每批写成 KPOINTS_OPT_b01 …')
    s.add_argument('--nkbatch', type=int, default=36, help='c5 每批点数，INCAR 的 KPOINTS_OPT_NKBATCH 必须相同')
    s.add_argument('--internal-ref', action='store_true', help='NKRED 与 c2 不同时：带边参考取本次计算中的 c2 校验点，而不是 c2 的数值')
    s.add_argument('--win-v', type=float, default=0.060, help='价带极大离价带顶多少 eV 以内做星形')
    s.add_argument('--win-c', type=float, default=0.100, help='导带极小离导带底多少 eV 以内做星形')
    s.add_argument('--delta', type=float, default=0.02); s.add_argument('--delta-z', type=float, default=0.03)
    s.add_argument('--show', type=int, default=8)
    s = sub.add_parser('fit'); s.add_argument('dir')
    s.add_argument('--parts', nargs='+', help='逐批作业的目录（按批号顺序），代替 dir 里的单个 PROCAR_OPT')
    s = sub.add_parser('local'); s.add_argument('c4dir'); s.add_argument('--out', required=True)
    s.add_argument('--du', type=frac, default=1 / 60); s.add_argument('--dw', type=frac, default=1 / 36)
    s.add_argument('--nu', type=int, default=3); s.add_argument('--nw', type=int, default=4)
    s.add_argument('--delta', type=float, default=0.02, help='垂直镜面方向的步长（Å⁻¹）')
    s.add_argument('--win-v', type=float, default=0.060, help='只细化离价带顶这么多 eV 以内的极大')
    s.add_argument('--nkbatch', type=int, default=36); s.add_argument('--split', action='store_true')
    s = sub.add_parser('local-fit'); s.add_argument('dir')
    s.add_argument('--parts', nargs='+', help='逐批作业的目录（按批号顺序）')
    a = ap.parse_args()
    {'grid': cmd_grid, 'analyze': cmd_analyze, 'fit': cmd_fit, 'local': cmd_local, 'local-fit': cmd_local_fit}[a.cmd](a)


if __name__ == '__main__':
    main()
