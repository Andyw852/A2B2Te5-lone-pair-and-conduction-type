"""HSE 能谷的多谷模型：同浓度的 Seebeck 系数（p、n 分开，含双极），能谷直接取自 c2 路径。

建谷（自动，窗口 --win，默认 300 meV）：沿 c2 HSE 路径找价带、导带的极值（10_脚本/c2_valleys.py 的口径），然后
- 对称点 Γ、A、L、M、K、H：三个主方向（两个面内、一个沿 c）在路径上都有，直接用路径上的真实色散 ε(k)（非抛物）；
  三个方向都是极值才算谷，鞍点舍去（例如 PbBi 的 Γ 价带：面内极大、沿 c 极小）。等价谷数 Γ、A 为 1，M、L 为 3，K、H 为 2。
- Γ–A 线上的驼峰（kz = ±k0）：沿 c 用路径色散（朝 Γ 一侧到 Γ 截止），面内借用 Γ 处的色散形状；g = 2。
- 线段内部的离轴谷（Γ–M、K–Γ、A–L、H–A，g = 6；M–K、L–H 在区界上，g = 3）：沿线用路径色散（两侧平均），
  面内垂直方向取沿线的抛物质量，沿 c 取余弦估计 m_z = 2·7.62/(W c²)，W 为同一面内位置 kz = 0 与 1/2 的能量差；
  W ≤ 0 说明沿 c 是鞍点，舍去。c3_valley 拟合过的谷，只有在采样点处确认为真极值时（三个主曲率同号，且由梯度与曲率
  外推的驻点落在采样半径以内）才用拟合质量（抛物）整体替换；否则保留路径估计，并在 --list 中注明原因。
  注意：Γ–M、A–L 线在镜面内而不在二重轴上，对称性不要求沿 c 的梯度为零，离轴谷的真实极大可以在镜面内的一般 kz 处，
  路径给出的离轴谷能量只是下限（SnBi 的 c3 即如此）。
- 真实谷（c4_mirror_v3 + c5_star，目前只有 SnBi 的 HSE 价带）：有这两步的结果时，镜面网格覆盖范围内（Γ–M、A–L 线上
  |u| ≤ 网格边界）的路径离轴谷全部换成镜面内的真实极大：c5 星形拟合过的谷用拟合质量（抛物，主轴质量算态数、线与面内垂直
  方向的质量算面内速度）与极值能量；只在 c4 网格上找到、未细化的谷用网格点能量（下限）和网格二次拟合的（线, c）质量，
  面内垂直质量取 c5 两个谷的平均。价带顶改为真实价带顶，其余路径谷的离带边能量、带隙随之改为相对真实价带顶。g = 6。
每个谷当作“三个半轴 k_i(ε) 取自真实色散”的椭球：单个口袋的态数 N(ε) = 2·(4π/3)k1k2k3/(2π)³，
面内 ⟨v_x²⟩ ≈ (v1² + v2²)/6，v_i = dε/dk_i。

输运：两种散射口径。--tau crta（默认）常数弛豫时间，σ(ε) ∝ g(ε)⟨v_x²⟩，与 BoltzTraP 的口径相同，可以和 2.4 节旧的
PBE 数值直接比；--tau adp 谷内声学形变势，τ ∝ 1/g(ε)，σ(ε) ∝ ⟨v_x²⟩（各谷同一形变势，不含谷间散射）。
S = Σ∫σ(ε)(±(ε − μ))(−∂f/∂ε)dε / (T Σ∫σ(ε)(−∂f/∂ε)dε)。给定净掺杂 p − n 解费米能级，少数载流子同样计入。
注意：PbSb 的 Γ 带几乎是线性（类 Dirac）色散，声学形变势下输运函数近乎常数，简并时 S 很小；早先 AMSET 显示 POP 主导，
所以两种口径给出的是这一类带的上下范围。

这是半定量模型：椭球近似、单带序号（不含更高的导带/更深的价带）、散射简化。用于比较材料之间、p 与 n 之间的相对大小，不替代 BoltzTraP/AMSET 的绝对值。

用法：python 10_脚本/multivalley_seebeck.py [--T 300] [--conc 1e19 3e19 1e20] [--materials ...] [--list]
"""
import argparse, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
KB = 8.617333262e-5
HB = 7.62                                  # ħ²/m_e，eV·Å²
BREAKS = [(187, 188), (195, 196)]
SEGS = [(0, 33), (33, 52), (52, 90), (90, 97), (97, 130), (130, 149), (149, 187), (188, 195), (196, 203)]
# 对称点：(路径序号, [(起点, 方向) × 2 面内, (起点, 方向) 沿 c], g)
POINTS = {'Γ': ((0, 90), [(0, 1), (90, -1)], [(90, 1)], 1), 'A': ((97, 187), [(97, 1), (187, -1)], [(97, -1)], 1),
          'L': ((130, 188), [(130, -1), (130, 1)], [(188, 1)], 3), 'M': ((33, 195), [(33, -1), (33, 1)], [(195, -1)], 3),
          'K': ((52, 203), [(52, -1), (52, 1)], [(203, -1)], 2), 'H': ((149, 196), [(149, -1), (149, 1)], [(196, 1)], 2)}
EPS = np.linspace(0, 0.8, 3201)            # 进入带内的能量网格，eV
STAR = (0.02, 0.03)                        # c3 星形采样的半径：沿线与面内垂直 δ、沿 c δz（Å⁻¹）
MZ_SCALE = 1.0                             # 余弦估计的 m_z 乘的系数（敏感性检验用）


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


def kdist(K, B):
    x = [0.0]
    for i in range(1, len(K)):
        x.append(x[-1] + (0.0 if (i - 1, i) in BREAKS else np.linalg.norm((K[i] - K[i - 1]) @ B)))
    return np.array(x)


def seg_of(i):
    return next(s for s in SEGS if s[0] <= i <= s[1])


def ray(x, e, sgn, i0, d):
    """从 i0 沿 d 走到线段端点：返回 (k, ε)，ε = sgn·(e − e[i0]) ≥ 0，只保留单调上升的部分。"""
    a, b = seg_of(i0 + d if i0 in (0, 90, 97, 187, 130, 188, 33, 195, 52, 203, 149, 196) else i0)
    idx = list(range(i0, b + 1)) if d > 0 else list(range(i0, a - 1, -1))
    k = np.abs(x[idx] - x[i0]); eps = sgn * (e[idx] - e[i0])
    keep = [0]
    for j in range(1, len(idx)):
        if eps[j] <= eps[keep[-1]]: break
        keep.append(j)
    return k[keep], eps[keep]


def kane_fit(kr, er, emax=0.30):
    """按 Kane 两带关系 ε(1 + ε/E_K) = (ħ²/2m) k² 拟合射线（取 ε ≤ emax 的单调部分）：返回 (a = ħ²/2m, b = 1/E_K)。
    路径步长约 0.025 Å⁻¹，轻带在 0.25 eV 内只有 2–3 个点，逐段线性插值会使速度失真，所以用这个两参数式。
    b 限制为 ≥ 0（拟合给出超抛物时退回抛物）。"""
    sel = (er > 0) & (er <= emax)
    if sel.sum() < 2: sel = er > 0; sel[3:] = False
    k2, e = kr[sel] ** 2, er[sel]
    if len(e) == 1:
        return e[0] / k2[0], 0.0
    A = np.vstack([k2, -e ** 2]).T                       # ε = a k² − b ε²
    (a, bb), *_ = np.linalg.lstsq(A, e, rcond=None)
    if bb < 0 or a <= 0:
        return float(np.sum(k2 * e) / np.sum(k2 * k2)), 0.0
    return a, bb


def k_and_v2_from_ray(kr, er, extend=True):
    """Kane 拟合给出 k(ε)、v²(ε)。extend=False（沿 c）时 k 不超过射线终点（布里渊区边界或驼峰的另一侧），超出部分速度记 0。"""
    if len(kr) < 2:
        return np.zeros_like(EPS), np.zeros_like(EPS)
    a, b = kane_fit(kr, er)
    k = np.sqrt(EPS * (1 + b * EPS) / a)
    v2 = 4 * a * EPS * (1 + b * EPS) / (1 + 2 * b * EPS) ** 2
    if not extend:
        cap = k > kr[-1]
        k = np.where(cap, kr[-1], k); v2 = np.where(cap, 0.0, v2)
    return k, v2


def parabolic(m):
    return np.sqrt(2 * m * EPS / HB), 2 * HB * EPS / m


def line_mass(x, e, i):
    p = np.polyfit(x[i - 2:i + 3] - x[i], e[i - 2:i + 3], 2)
    return abs(HB / (2 * p[0]))


def c3_load(material):
    f = os.path.join(ROOT, '22_能带基准_1-0', material, 'c3_valley', 'c3_masses.json')
    return json.load(open(f, encoding='utf-8')) if os.path.exists(f) else None


def c3_verdict(entry, kind='v'):
    """c3 拟合的谷能否直接用：三个主曲率同号（价带为极大、导带为极小），且驻点 Δk = −H⁻¹g 在采样半径以内。"""
    H = np.array(entry['H_eVA2_line_perp_c']); g = np.array(entry['gradient'])
    w = np.linalg.eigvalsh(H)
    if not (np.all(w < 0) if kind == 'v' else np.all(w > 0)):
        return False, '主曲率异号，是鞍点'
    dk = -np.linalg.solve(H, g); dE = -0.5 * g @ np.linalg.solve(H, g)
    if abs(dk[0]) > STAR[0] or abs(dk[1]) > STAR[0] or abs(dk[2]) > STAR[1]:
        return False, f'驻点偏离采样点 {np.round(dk, 3).tolist()} Å⁻¹（线、面内垂直、c），超出采样半径；二次外推高 {1000 * abs(dE):.0f} meV'
    return True, ''


def c3_usable(material, kind='v'):
    """返回 (可直接用的谷 {序号: 拟合}, 不能用的谷 {序号: 原因})；没有 c3 时返回 (None, {})。"""
    c3 = c3_load(material)
    if c3 is None: return None, {}
    ok, bad = {}, {}
    for k, v in c3.items():
        good, why = c3_verdict(v, kind)
        (ok.__setitem__(k, v) if good else bad.__setitem__(k, why))
    return ok, bad


def c3_provisional(material):
    """没有 c3，或 c3 中有谷不能直接用：该材料 p 型的多谷估计仍为暂定。已有 c4/c5 的真实谷时不再暂定。"""
    if USE_TRUE and true_load(material): return False
    ok, bad = c3_usable(material)
    return ok is None or bool(bad)


USE_TRUE = True                            # 有 c4/c5 结果时用镜面内的真实价带谷（--path-only 关闭）


def true_load(material):
    """c4_mirror_v3 的镜面极大 + c5_star 的星形拟合 → 真实价带谷。没有这两步结果时返回 None。"""
    d = os.path.join(ROOT, '22_能带基准_1-0', material)
    f4, f5 = os.path.join(d, 'c4_mirror_v3', 'c4_extrema.json'), os.path.join(d, 'c5_star', 'c5_masses.json')
    if not (os.path.exists(f4) and os.path.exists(f5)): return None
    c4, c5 = json.load(open(f4, encoding='utf-8')), json.load(open(f5, encoding='utf-8'))
    c4g = json.load(open(os.path.join(d, 'c4_mirror_v3', 'c4_grid.json'), encoding='utf-8'))
    v5 = [v for v in c5['valleys'].values() if v['kind'] == 'v' and v['true_extremum_within_star']]
    mperp = float(np.mean([v['m'][1] for v in v5]))
    out = []
    for vm in c4['vb_maxima']:
        hit = next((v for v in v5 if abs(v['center_frac'][0] - vm['u']) < 2e-3 and abs(v['center_frac'][2] - vm['w']) < 2e-3), None)
        if hit:
            out.append(dict(u=hit['center_frac'][0], w=hit['center_frac'][2], E=hit['E_extremum'], m=hit['m'], principal=hit['principal'],
                            s=hit['cation_s'], src='c5 星形拟合（真极值）'))
        else:
            H = np.array(vm['H_line_c_eVA2']); ev = np.linalg.eigvalsh(H)
            if np.any(ev >= 0): continue
            out.append(dict(u=vm['u'], w=vm['w'], E=vm['E'], m=[HB / abs(H[0, 0]), mperp, HB / abs(H[1, 1])],
                            principal=[HB / abs(ev[0]), mperp, HB / abs(ev[1])], s=vm['cation_s'],
                            src=f'c4 网格点（未细化，能量为下限）；面内垂直质量取 c5 平均 {mperp:.3f}'))
    return dict(vbm=c5['vbm'], path_vbm=c5['path_vbm'], umin=min(c4g['us']), valleys=out)


def build(material, kind, win, c3=None, side='c2_hse'):
    nel = 88 if material[3:5] == 'Bi' else 68
    base = os.path.join(ROOT, '22_能带基准_1-0', '' if material == 'Pb2Sb2Te5' else material)
    E = read_path_eigen(os.path.join(base, side, 'path_eigen.txt'))
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    A = read_lat(os.path.join(ROOT, '08_介电复核', '源数据', material, 'POSCAR')); B = 2 * np.pi * np.linalg.inv(A).T
    x = kdist(K, B); c = A[2, 2]
    e = E[:, nel - 1] if kind == 'v' else E[:, nel]
    sgn = -1 if kind == 'v' else 1
    edge = e.max() if kind == 'v' else e.min()
    gap = E[:, nel].min() - E[:, nel - 1].max()
    T5 = true_load(material) if (USE_TRUE and kind == 'v' and side == 'c2_hse') else None
    shift = 0.0
    if T5:                                                       # 价带顶改为镜面内的真实价带顶
        shift = T5['vbm'] - e.max(); gap -= shift
        ucut = -T5['umin'] * 66                                  # Γ–M、A–L 线上 |u| ≤ 网格边界的路径序号跨度（M、L 在 66·½ = 33 步处）
    from c2_valleys import extrema
    valleys, seen = [], set()
    for i, (dE, _) in extrema(x, e, 'max' if kind == 'v' else 'min', win):
        off = abs(dE) / 1000 + shift
        if T5 and (0 < i <= ucut or 97 < i <= 97 + ucut): continue   # 由镜面网格的真实谷替代
        name = next((n for n, (ids, *_r) in POINTS.items() if i in ids), None)
        if name:
            if name in seen: continue
            seen.add(name)
            ids, inpl, alongc, g = POINTS[name]
            rays = [ray(x, e, sgn, a0, d) for a0, d in inpl + alongc]
            if any(len(r[0]) < 2 for r in rays): continue          # 某个方向不是极值：鞍点
            (k1, v1), (k2, v2), (k3, _) = [k_and_v2_from_ray(*r) for r in rays]
            valleys.append(dict(name=name, off=off, g=g, k=(k1, k2, k3), kcap=rays[2][0][-1], v2=(v1, v2), src='路径色散'))
        elif 90 < i < 97:                                        # Γ–A 驼峰
            up, dn = ray(x, e, sgn, i, 1), ray(x, e, sgn, i, -1)
            kup, _ = k_and_v2_from_ray(*up, extend=False); kdn, _ = k_and_v2_from_ray(*dn, extend=False)
            g0 = [ray(x, e, sgn, 0, 1), ray(x, e, sgn, 90, -1)]
            if any(len(r[0]) < 2 for r in g0): continue
            (k1, v1), (k2, v2) = [k_and_v2_from_ray(*r) for r in g0]
            valleys.append(dict(name=f'Γ–A kz=±{abs(K[i][2]):.3f}', off=off, g=2, k=(k1, k2, (kup + kdn) / 2), kcap=np.inf, v2=(v1, v2),
                                src='沿 c 路径色散；面内借用 Γ 的色散'))
        else:
            a, b = seg_of(i)
            if i - 2 < a or i + 2 > b: continue
            g = 3 if (33 < i < 52 or 130 < i < 149) else 6
            j = i + 97 if i <= 90 else i - 97
            W = sgn * (e[j] - e[i])
            lab = f'#{i} {np.round(K[i], 3).tolist()}'
            if c3 and str(i) in c3:
                m = c3[str(i)]['m']
                (k1, v1), (k2, v2), (k3, _) = parabolic(m[0]), parabolic(m[1]), parabolic(m[2])
                valleys.append(dict(name=lab, off=off, g=g, k=(k1, k2, k3), kcap=np.pi / c, v2=(v1, v2), src='c3 拟合质量（抛物）'))
                continue
            if W <= 0: continue                                     # 沿 c 是鞍点
            mz = MZ_SCALE * 2 * HB / (W * c * c)
            r1, r2 = ray(x, e, sgn, i, 1), ray(x, e, sgn, i, -1)
            (ka, va), (kb, vb) = k_and_v2_from_ray(*r1), k_and_v2_from_ray(*r2)
            kl = (ka + kb) / 2; vl = (va + vb) / 2
            mp = line_mass(x, e, i)
            (kp, vp), (kz, _) = parabolic(mp), parabolic(mz)
            valleys.append(dict(name=lab, off=off, g=g, k=(kl, kp, kz), kcap=np.pi / c, v2=(vl, vp),
                                src=f'沿线路径色散；面内垂直取沿线质量 {mp:.3f}、沿 c 余弦估计 {mz:.2f}'))
    if T5:
        for t in T5['valleys']:
            off = T5['vbm'] - t['E']
            if off > win / 1000: continue
            p, m = t['principal'], t['m']
            (k1, _), (k2, _), (k3, _) = parabolic(p[0]), parabolic(p[1]), parabolic(p[2])
            (_, v1), (_, v2) = parabolic(m[0]), parabolic(m[1])          # 面内速度：沿线、面内垂直两个方向的质量
            valleys.append(dict(name=f'mirror ({t["u"]:.4f},0,{t["w"]:.4f})', off=off, g=6, k=(k1, k2, k3), kcap=np.pi / c, v2=(v1, v2),
                                src=t['src'], kfrac=(t['u'], 0.0, t['w']), s=t['s']))
    for v in valleys:
        k1, k2, k3 = v['k']; K = v['kcap']
        # 椭球沿 c 被布里渊区边界 |kz| ≤ K 截断：体积 = π k1 k2 ∫(1 − z²/k3²)dz，k3 → ∞ 时趋于柱面 2πk1k2K
        if np.isinf(K):
            vol = 4 * np.pi / 3 * k1 * k2 * k3
        else:
            vol = np.where(k3 <= K, 4 * np.pi / 3 * k1 * k2 * k3, np.pi * k1 * k2 * (2 * K - 2 * K ** 3 / (3 * np.maximum(k3, 1e-12) ** 2)))
        v['N'] = 2 * vol / (2 * np.pi) ** 3 * 1e24                                              # cm⁻³，单个口袋
        v['vx2'] = (v['v2'][0] + v['v2'][1]) / 6
        v['dos'] = np.gradient(v['N'], EPS)
    return valleys, gap


def _occ(T):
    kT = KB * T

    def occ(eta):            # 每个 ε 网格点的占据与 −∂f/∂ε
        f = 1 / (1 + np.exp(np.clip(EPS / kT - eta, -700, 700)))
        return f, f * (1 - f) / kT
    return occ


def fermi_level(V, C, gap, N, T):
    """由电中性 p − n = N 解费米能级（相对价带顶，eV），并返回此时的空穴、电子浓度。"""
    kT = KB * T; occ = _occ(T)

    def carriers(EF):
        p = sum(v['g'] * np.trapezoid(v['N'] * occ((-v['off'] - EF) / kT)[1], EPS) for v in V)
        n = sum(c['g'] * np.trapezoid(c['N'] * occ((EF - gap - c['off']) / kT)[1], EPS) for c in C)
        return p, n
    lo, hi = -1.0, gap + 1.0
    for _ in range(70):
        mid = (lo + hi) / 2; p, n = carriers(mid)
        lo, hi = (mid, hi) if p - n > N else (lo, mid)
    EF = (lo + hi) / 2
    return EF, carriers(EF)


def solve(V, C, gap, N, T, tau='crta'):
    """tau='crta'：常数弛豫时间，σ(ε) ∝ g(ε)⟨v_x²⟩（与 BoltzTraP 的口径相同）；tau='adp'：声学形变势，τ ∝ 1/g(ε)，σ(ε) ∝ ⟨v_x²⟩。"""
    kT = KB * T; occ = _occ(T)
    EF, _ = fermi_level(V, C, gap, N, T)
    num = den = 0.0; parts = []
    for sgn, L in ((1, V), (-1, C)):
        for v in L:
            eta = ((-v['off'] - EF) if sgn > 0 else (EF - gap - v['off'])) / kT
            w = occ(eta)[1]
            sig = v['vx2'] * (v['dos'] if tau == 'crta' else 1.0)
            s0 = v['g'] * np.trapezoid(sig * w, EPS)
            s1 = v['g'] * np.trapezoid(sig * w * (EPS - eta * kT), EPS)
            num += sgn * s1; den += s0; parts.append(('v' if sgn > 0 else 'c', v['name'], s0))
    return num / den / T * 1e6, EF, [(s, nm, x / den) for s, nm, x in parts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--T', type=float, default=300)
    ap.add_argument('--conc', type=float, nargs='+', default=[1e19, 3e19, 1e20])
    ap.add_argument('--materials', nargs='+', default=['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5'])
    ap.add_argument('--win', type=float, default=300)
    ap.add_argument('--list', action='store_true', help='列出每种材料建出的谷')
    ap.add_argument('--tau', choices=('crta', 'adp'), default='crta', help='散射：常数弛豫时间（默认）或声学形变势')
    ap.add_argument('--side', default='c2_hse', help='c2_hse（默认）或 c2_pbe（用 PBE 路径，与旧 BoltzTraP2 结果对照）')
    ap.add_argument('--path-only', action='store_true', help='不用 c4/c5 的真实谷，只用路径（与旧结果对照）')
    ap.add_argument('--mz-scale', type=float, default=1.0, help='余弦估计的沿 c 质量乘以该系数（检验敏感性；SnSb Γ 价带的检验显示余弦估计约偏轻一半）')
    a = ap.parse_args()
    global MZ_SCALE, USE_TRUE
    MZ_SCALE = a.mz_scale; USE_TRUE = not a.path_only
    c3all, c3bad = {}, {}
    for m in a.materials:
        ok, bad = c3_usable(m)
        if ok: c3all[m] = ok
        if bad: c3bad[m] = bad
    models = {}
    for m in a.materials:
        V, gap = build(m, 'v', a.win, c3all.get(m) if a.side == 'c2_hse' else None, a.side); C, _ = build(m, 'c', a.win, side=a.side)
        models[m] = (V, C, gap)
        if a.list:
            print(f'== {m}：带隙 {1000 * gap:.1f} meV' + ('（含 c3）' if m in c3all else ''))
            for k, why in c3bad.get(m, {}).items():
                print(f'   c3 谷 #{k} 未采用：{why}')
            for side, L in (('价带', V), ('导带', C)):
                for v in L:
                    print(f'   {side} {v["name"]:28s} 离带边 {1000 * v["off"]:6.1f} meV  g={v["g"]}  {v["src"]}')
    print(f'# {a.side}；T = {a.T:.0f} K；S 单位 µV/K；面内；多谷模型（路径色散椭球，Kane 拟合）；散射 {"常数弛豫时间" if a.tau == "crta" else "声学形变势 τ∝1/g"}；窗口 {a.win:.0f} meV')
    print(f'{"材料":10s}' + ''.join(f'  p {c:.0e}' for c in a.conc) + ''.join(f'  n {c:.0e}' for c in a.conc) + '   p 型 S 最大 / n 型 |S| 最大（扫 1e15–1e21，含双极）')
    for m in a.materials:
        V, C, gap = models[m]
        Sp = [solve(V, C, gap, c, a.T, a.tau)[0] for c in a.conc]
        Sn = [solve(V, C, gap, -c, a.T, a.tau)[0] for c in a.conc]
        grid = np.logspace(15, 21, 37)
        smax = max(solve(V, C, gap, c, a.T, a.tau)[0] for c in grid); nmax = min(solve(V, C, gap, -c, a.T, a.tau)[0] for c in grid)
        print(f'{m:10s}' + ''.join(f'{s:9.0f}' for s in Sp) + ''.join(f'{s:9.0f}' for s in Sn) + f'   {smax:5.0f} / {nmax:5.0f}')
    print('\n# 3e19 cm⁻³ 时各谷对电导的份额（前 4 个）')
    for m in a.materials:
        V, C, gap = models[m]
        for N, lab in ((3e19, 'p'), (-3e19, 'n')):
            S, EF, parts = solve(V, C, gap, N, a.T, a.tau)
            top = sorted(parts, key=lambda t: -t[2])[:4]
            print(f'{m:10s} {lab}: S {S:5.0f}  E_F−VBM {1000 * EF:+6.0f} meV  ' + '；'.join(f'{s}:{nm} {100 * w:.0f}%' for s, nm, w in top))


if __name__ == '__main__':
    main()
