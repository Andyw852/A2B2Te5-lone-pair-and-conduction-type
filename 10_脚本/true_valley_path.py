"""穿过镜面内真实价带极大的附加能带路径（HSE06+SOC，KPOINTS_OPT）：Γ → Σ′ → R₂′ → L。

c2 的高对称路径只在 k_z = 0（Γ–M）与 k_z = ½（A–L）两处穿过孤对电子谷的“脊”，看不到真实极大（图 S1）。
这条附加路径全部在含 c 轴的镜面 (u, 0, w) 内，折点正好取在 c5 拟合的 Σ′、R₂′ 上，
所以能带图上可以直接画出真实价带顶（Sn₂Bi₂Te₅ 的 Σ′）和孤对电子谷的真实高度。

  make  <材料> --out <目录> [--dk 0.015] [--nkbatch 36] [--split]
        生成 KPOINTS_OPT（及 _bNN 分批文件）与 true_path.json。每批开头放一个 c2 路径上的点作锚点，用来核对。
        所有点经过 FOCK_FCC 奇数/24 检查；落在判据边界上的点沿路径挪开 0.2 步。
  check <目录> [--parts b01 b02 …]
        读 PROCAR_OPT：锚点与 c2 路径的差、折点处能量与 c5 极值的差，写出 true_path_eigen.csv（各带能量，相对真实价带顶）。

例：
  python 10_脚本/true_valley_path.py make Sn2Bi2Te5 --out 22_能带基准_1-0/Sn2Bi2Te5/c7_true_path --split
  python 10_脚本/true_valley_path.py check 22_能带基准_1-0/Sn2Bi2Te5/c7_true_path --parts <b01 b02>
"""
import argparse, csv, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import c4_mirror_scan as c4
import multivalley_seebeck as mv

L_PT = np.array([-0.5, 0.0, 0.5])
ANCHORS = (0, 11, 119, 130)                      # c2 路径上的 Γ、Σ、R₂、L，轮流作每批的锚点


def true_nodes(m):
    """Σ′、R₂′：c5 拟合的真实极大中，离路径 Σ (−1/6, 0, 0)、R₂ (−1/3, 0, ½) 最近的（笛卡尔距离）。"""
    T = mv.true_load(m)
    if T is None: sys.exit(f'{m} 没有 c4/c5 结果')
    B = c4.recip(m); nb1, cz = np.linalg.norm(B[0]), abs(B[2, 2])
    def dist(v, u0, w0):
        dw = abs(v['w'] - w0) % 1; dw = min(dw, 1 - dw)
        return np.hypot((v['u'] - u0) * nb1, dw * cz)
    c5 = [v for v in T['valleys'] if v['src'].startswith('c5')]
    sig = min(c5, key=lambda v: dist(v, -1 / 6, 0.0)); r2 = min(c5, key=lambda v: dist(v, -1 / 3, 0.5))
    return T, B, [('Γ', np.zeros(3)), ('Σ′', np.array([sig['u'], 0, sig['w']])), ('R₂′', np.array([r2['u'], 0, r2['w']])), ('L', L_PT)], (sig, r2)


def cmd_make(a):
    T, B, nodes, _ = true_nodes(a.material)
    pts, labels, seg = [], {}, []
    for (n0, p0), (n1, p1) in zip(nodes[:-1], nodes[1:]):
        ln = np.linalg.norm((p1 - p0) @ B); n = max(2, int(np.ceil(ln / a.dk)))
        for j in range(n + (1 if n1 == 'L' else 0)):
            t = j / n
            for shift in (0, 0.2, -0.2, 0.4, -0.4):                  # 折点（t = 0、1）不挪，它们本身不在判据边界上
                q = p0 + (t + (shift / n if 0 < j < n else 0)) * (p1 - p0)
                if not c4.fock_unsafe([q]): break
            else: sys.exit(f'{n0}→{n1} 第 {j} 点挪不开判据边界')
            if j == 0: labels[len(pts)] = n0
            pts.append(q)
        seg.append(ln)
    labels[len(pts) - 1] = 'L'
    if c4.fock_unsafe([nodes[1][1], nodes[2][1]]): sys.exit('折点 Σ′ 或 R₂′ 落在判据边界上')
    K2 = c4.read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    per = a.nkbatch - 1; allp, rec, checks, batches = [], [], [], []
    for b, s in enumerate(range(0, len(pts), per)):
        pi = ANCHORS[b % len(ANCHORS)]; start = len(allp)
        checks.append(dict(path_index=pi, k_index=len(allp))); allp.append(np.array(K2[pi]))
        for i in range(s, min(s + per, len(pts))): rec.append(len(allp)); allp.append(pts[i])
        batches.append([start, len(allp)])
    os.makedirs(a.out, exist_ok=True)
    title = f'{a.material} true-valley path Γ–Σ′–R₂′–L in the mirror plane (u,0,w); one c2 anchor per batch; KPOINTS_OPT_NKBATCH = {a.nkbatch}'
    c4.write_kpoints(os.path.join(a.out, 'KPOINTS_OPT'), title, allp)
    if a.split: c4.write_batches(a.out, 'KPOINTS_OPT', title, allp, batches)
    x = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(np.array(pts) @ B, axis=0), axis=1))])
    json.dump(dict(material=a.material, nelect=88 if a.material[3:5] == 'Bi' else 68, n_kpoints=len(allp), nkbatch=a.nkbatch,
                   batches=batches, checks=checks, path_k_index=rec, x_invA=x.tolist(), labels={str(k): v for k, v in labels.items()},
                   nodes={n: p.tolist() for n, p in nodes}, segment_lengths_invA=seg, vbm=T['vbm'], path_vbm=T['path_vbm']),
              open(os.path.join(a.out, 'true_path.json'), 'w'), indent=1, ensure_ascii=False)
    print(f'{len(pts)} 个路径点 + {len(checks)} 个锚点 = {len(allp)} 个 k 点，{len(batches)} 批 → {a.out}')
    print('折点：' + '；'.join(f'{n} {np.round(p, 4).tolist()}' for n, p in nodes) + f'；三段长度 {np.round(seg, 3).tolist()} Å⁻¹，步长 ≤ {a.dk} Å⁻¹')


def cmd_check(a):
    meta = json.load(open(os.path.join(a.dir, 'true_path.json'))); m, nel = meta['material'], meta['nelect']
    E, fs = c4.read_parts(a.dir, a.parts, meta)
    if E.shape[0] != meta['n_kpoints']: sys.exit(f"PROCAR 有 {E.shape[0]} 个 k 点，应为 {meta['n_kpoints']}")
    P = c4.read_path_eigen(os.path.join(c4.mat_dir(m), 'c2_hse', 'path_eigen.txt'))
    dev = max(abs(E[c['k_index'], nel - 5:nel + 4] - P[c['path_index'], nel - 5:nel + 4]).max() for c in meta['checks'])
    print(f'# {m}：{len(meta["checks"])} 个锚点与 c2 路径最大差 {1000 * dev:.3f} meV' + ('' if dev < 5e-4 else '  ← 超过 0.5 meV，停下检查'))
    _, _, nodes, (sig, r2) = true_nodes(m)
    k = np.array(meta['path_k_index']); top = E[k, nel - 1]; lab = {int(i): n for i, n in meta['labels'].items()}
    for i, n in lab.items():
        if n in ('Σ′', 'R₂′'):
            ref = sig if n == 'Σ′' else r2
            print(f'{n}：本路径上 {top[i]:.5f} eV，c5 极值 {ref["E"]:.5f} eV，差 {1000 * (top[i] - ref["E"]):+.2f} meV')
    print(f"本路径上价带最高 {top.max():.5f} eV（第 {int(top.argmax())} 点），真实价带顶 {meta['vbm']:.5f} eV，差 {1000 * (top.max() - meta['vbm']):+.2f} meV")
    with open(os.path.join(a.dir, 'true_path_eigen.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); bands = range(nel - 8, nel + 6)
        w.writerow(['i', 'x_invA', 'label', 'k1', 'k2', 'k3'] + [f'band_{b + 1}' for b in bands] + ['cation_s_topVB'])
        K = c4.read_kpoints_opt(os.path.join(a.dir, 'KPOINTS_OPT'))
        for j, (ki, xx) in enumerate(zip(k, meta['x_invA'])):
            w.writerow([j, round(xx, 5), lab.get(j, '')] + [round(v, 6) for v in K[ki]] + [round(E[ki, b] - meta['vbm'], 5) for b in bands]
                       + [round(100 * fs[ki, nel - 1], 2)])
    print(f"# 写出 {os.path.join(a.dir, 'true_path_eigen.csv')}（能量相对真实价带顶，eV）")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter); sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('make'); s.add_argument('material'); s.add_argument('--out', required=True)
    s.add_argument('--dk', type=float, default=0.015); s.add_argument('--nkbatch', type=int, default=36); s.add_argument('--split', action='store_true')
    s = sp.add_parser('check'); s.add_argument('dir'); s.add_argument('--parts', nargs='*')
    a = ap.parse_args(); {'make': cmd_make, 'check': cmd_check}[a.cmd](a)


if __name__ == '__main__':
    main()
