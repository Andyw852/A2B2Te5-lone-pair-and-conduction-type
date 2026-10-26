"""22_：TRIM 点宇称（不依赖 irrep，直接读 VASP WAVECAR 的平面波系数；用于与 irrep 交叉核对）。

用法（在 c_hse_A 目录里运行，需要 numpy）：
  python trim_parity.py WAVECAR --kpoints 1,7,20,26 --nelect 68 [--bands 65-72]
口径：反演中心在原点（1-0 的四种材料都是，中心 Te 在 (0,0,0)）。TRIM 点 k 满足 2k 为倒格矢，
  ψ_k(−r) 的平面波系数为 c(−G−2k)，宇称 ξ = Σ_G c*(G)·c(−G−2k) / Σ_G |c(G)|²（旋量两分量相加）。
  简并能级（能量差 < --degtol，如 SOC 下的 Kramers 对）在子空间内求 O⁻¹P（O 为平面波重叠矩阵）的本征值，应为 ±1。
  用的是赝波函数（PAW 平面波部分），与 irrep、irvsp 的做法相同。
输出：每个 TRIM 点上 --bands 范围内各能级的能量（eV，可与 EIGENVAL 核对）与宇称，
  以及第 1…NELECT 条带（每个 Kramers 对计一次）的宇称乘积 δ（--first 可改起点）。六方 P-3m1 的强 Z2：(−1)^ν = δ_Γ·δ_A·δ_M³·δ_L³。
  一对关于反演中心对称的原子，其满壳层在每个 TRIM 点上贡献相同的因子（d 壳层为 (−1)^5 = −1），
  所以半芯层只会让 8 个 TRIM 点的 δ 一起翻号，不改变 Z2，也不改变“哪个 TRIM 点与其余不同”。
自查：每个 k 点生成的 G 矢量个数必须等于 WAVECAR 记录的平面波数（不等就报错停止）；宇称必须接近 ±1（偏离 > 0.05 会标出）。
"""
import argparse, numpy as np

C = 0.262465831            # 2m/ħ²，单位 1/(eV·Å²)


def read_header(f):
    recl, nspin, rtag = np.fromfile(f, np.float64, 3)
    f.seek(int(recl))
    h = np.fromfile(f, np.float64, 12)
    nk, nb, encut = int(h[0]), int(h[1]), h[2]
    return int(recl), int(nspin), int(rtag), nk, nb, encut, h[3:12].reshape(3, 3)


def gvectors(k, a, encut, nplw_half):
    """按 VASP 的顺序生成 |k+G|²/C < ENCUT 的 G（x 最快变化；各分量 0,1,…,n,−n,…,−1）。"""
    b = 2 * np.pi * np.linalg.inv(a).T
    n = (np.ceil(np.sqrt(encut * C) / np.linalg.norm(b, axis=1)) + 2).astype(int)
    rng = [np.r_[0:n[i] + 1, -n[i]:0] for i in range(3)]
    G = []
    for g3 in rng[2]:
        for g2 in rng[1]:
            for g1 in rng[0]:
                v = (k + np.array([g1, g2, g3])) @ b
                if v @ v / C < encut:
                    G.append((g1, g2, g3))
    G = np.array(G, dtype=int)
    if len(G) != nplw_half:
        raise SystemExit(f'G 矢量个数 {len(G)} ≠ WAVECAR 中的 {nplw_half}：读法与这份 WAVECAR 不符，停止')
    return G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('wavecar'); ap.add_argument('--kpoints', required=True, help='1 起的 k 点序号（与 IBZKPT 相同），逗号分隔')
    ap.add_argument('--nelect', type=int, required=True); ap.add_argument('--bands', default=None, help='如 65-72（1 起）')
    ap.add_argument('--first', type=int, default=1, help='宇称乘积从第几条带算起（默认 1；如 21 可去掉 Bi 5d 半芯的 20 条带）')
    ap.add_argument('--collinear', action='store_true', help='非 SOC 计算时加（默认按旋量读）'); ap.add_argument('--degtol', type=float, default=1e-3)
    a = ap.parse_args()
    a.spinor = not a.collinear
    f = open(a.wavecar, 'rb')
    recl, nspin, rtag, nk, nb, encut, A = read_header(f)
    cdt = {45200: np.complex64, 45210: np.complex128}[rtag]
    assert nspin == 1, '非共线计算 ISPIN 应为 1'
    b0, b1 = (map(int, a.bands.split('-')) if a.bands else (max(1, a.nelect - 3), min(nb, a.nelect + 4)))
    print(f'WAVECAR: nk={nk} nbands={nb} ENCUT={encut} eV  精度={cdt.__name__}  NELECT={a.nelect}')
    deltas = {}
    for ik in [int(x) - 1 for x in a.kpoints.split(',')]:
        rec = 2 + ik * (nb + 1)
        f.seek(rec * recl)
        hdr = np.fromfile(f, np.float64, 4 + 3 * nb)
        nplw = int(hdr[0]); k = hdr[1:4]; E = hdr[4::3][:nb]
        if not np.allclose(2 * k, np.round(2 * k), atol=1e-6):
            raise SystemExit(f'k 点 {ik + 1} = {k} 不是 TRIM 点')
        nG = nplw // 2 if a.spinor else nplw
        G = gvectors(k, A, encut, nG)
        idx = {tuple(g): i for i, g in enumerate(G)}
        K0 = np.round(2 * k).astype(int)
        perm = np.array([idx[tuple(-g - K0)] for g in G])      # c(−G−2k) 在原列表中的位置
        def coeff(ib):
            f.seek((rec + 1 + ib) * recl)
            c = np.fromfile(f, cdt, nplw).astype(np.complex128)
            return c.reshape(2, nG) if a.spinor else c.reshape(1, nG)
        # 按能量分组（简并子空间）
        groups, cur = [], [0]
        for ib in range(1, max(b1, a.nelect)):
            if abs(E[ib] - E[cur[-1]]) < a.degtol: cur.append(ib)
            else: groups.append(cur); cur = [ib]
        groups.append(cur)
        delta, report = 1.0, []
        for g in groups:
            cs = [coeff(ib) for ib in g]
            # 子空间内：O_mn = <c_m|c_n>，P_mn = <c_m|P c_n>；赝波函数在平面波内积下不严格正交（PAW 的正交性在 S 度规下），
            # 所以取 O⁻¹P 的本征值；子空间在反演下封闭时它们严格为 ±1
            O = np.array([[np.vdot(cm, cn) for cn in cs] for cm in cs])
            P = np.array([[np.vdot(cm, cn[:, perm]) for cn in cs] for cm in cs])
            M = np.linalg.solve(O, P)
            ev = np.sort(np.linalg.eigvals(M).real)
            xi = np.trace(M).real / len(g)
            nneg = int(np.sum(ev < 0))
            if a.first - 1 <= g[0] and g[-1] < a.nelect:        # 占据带：每个 Kramers 对计一次 → (−1)^(负本征值个数/2)
                delta *= (-1) ** (nneg // 2) if len(g) % 2 == 0 else (-1) ** nneg
            elif g[0] < a.nelect:
                print(f'   注意：带 {g[0] + 1}–{g[-1] + 1} 跨过第 NELECT 条，价带顶与导带底简并')
            if g[-1] + 1 >= b0 and g[0] + 1 <= b1:
                flag = '' if np.all(np.abs(np.abs(ev) - 1) < 0.05) else '  ← 偏离 ±1'
                report.append(f'   带 {g[0] + 1:>3}–{g[-1] + 1:<3} E = {E[g[0]]:9.4f} eV  宇称本征值 {np.round(ev, 3).tolist()}  ξ = {xi:+.3f}{flag}')
        deltas[ik + 1] = int(np.round(delta))
        print(f'\nk 点 {ik + 1}  k = {np.round(k, 4).tolist()}  平面波 {nG}×{2 if a.spinor else 1}')
        print('\n'.join(report))
        print(f'   第 {a.first}…{a.nelect} 条带（每个 Kramers 对计一次）的宇称乘积 δ = {deltas[ik + 1]:+d}')
    print('\nδ:', deltas)


if __name__ == '__main__':
    main()
