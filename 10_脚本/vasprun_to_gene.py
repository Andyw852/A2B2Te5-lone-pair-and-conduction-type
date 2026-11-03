#!/usr/bin/env python3
"""从 VASP 的 vasprun.xml（可达数 GB）生成 BoltzTraP2 的 GENE 输入，可抽稀到较疏的子网格（2026-11-03）。

两种用途：
1. 网格收敛检验（PBE+SOC）：把 28×28×7/6 密网格抽成 14×14×7 等子网格，与完整密网格比较 BoltzTraP2 的 S，
   找出够用的最小网格，再按它算 HSE。
2. HSE06+SOC 的 KPOINTS_OPT 非自洽结果（--opt）：读 <eigenvalues_kpoints_opt> 里的能带。分批计算时一次给全部批次的
   vasprun.xml，合并、去掉重复的锚点（核对各批在锚点上一致），再核对恰好是 --mesh 的一套不可约点（每类一个）。

流式读取（不把整个文件读进内存），只取：参数（LHFCALC、LSORBIT、NELECT）、自洽网格 divisions、
k 点表、自洽或 KPOINTS_OPT 的本征值、末态结构（finalpos）与元素。GENE 的写法与 eigenval_to_gene.write_gene 相同
（每对 Kramers 简并能带写一条、只写费米能级附近、扣除深能级电子），所以 run_bt2_hse_grid.py 可以直接跑。

抽稀（--sub nx ny nz）：只保留 k·n 为整数的点。Γ 中心网格的子格在对称操作下不变，所以不可约点抽稀后仍是子网格的一套
不可约点。
全网格输入（ISYM = −1，例如 A2B2Te5_SOC_dense 的 PBE 密网格，vasprun.xml 里是全部 nx·ny·nz 个点）：先抽稀，再用 spglib
（含时间反演）约化到不可约点，与 HSE 的 KPOINTS_OPT 输入口径一致；约化时核对对称等价点的本征值相同（费米能级附近的能带，
差 ≤ EQUIV_TOL），不满足即报错。最后用 spglib 数一遍不可约点核对。需要 spglib。

用法：
  python 10_脚本/vasprun_to_gene.py --vasprun <.../vasprun.xml> --material Sn2Sb2Te5 --functional pbe \
      --sub 14 14 7 --out 22_能带基准_1-0/bt2_gridtest/Sn2Sb2Te5/pbe_14x14x7
  不给 --sub 即用完整网格；--opt 读 KPOINTS_OPT 的能带，此时用 --mesh 写明这些点属于哪个均匀网格，例如
  python 10_脚本/vasprun_to_gene.py --opt --mesh 28 28 6 --vasprun b01/vasprun.xml … b13/vasprun.xml \\
      --material Sn2Bi2Te5 --functional hse --out 22_能带基准_1-0/Sn2Bi2Te5/c8_hse_dense/gene
"""
import argparse, hashlib, os, sys
import xml.etree.ElementTree as ET
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eigenval_to_gene as g

FUNC = {'hse': 'HSE06+SOC', 'pbe': 'PBE+SOC'}
NELECT = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
ANCHOR_TOL = 1e-3                          # eV，不同 KPOINTS_OPT 批次在同一点（锚点）上本征值允许的差（c4 验收为 0.5 meV）
EQUIV_TOL = 5e-3                           # eV，对称等价点本征值允许的差（vasprun.xml 只有 4 位小数）
FLUSH = {'set', 'varray', 'array', 'projected', 'projected_kpoints_opt', 'dos', 'scstep', 'dynmat', 'partial', 'total'}


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def parse(path, opt):
    """返回 dict(params, divisions, kpts, eig (nk, nb) eV, basis, positions, symbols)。"""
    stack = []                                   # [(tag, attrib)]
    out = {'params': {}, 'divisions': None, 'kmain': [], 'kopt': [], 'basis': [], 'pos': [], 'symbols': []}
    emain, eopt = {}, {}
    for ev, el in ET.iterparse(path, events=('start', 'end')):
        if ev == 'start':
            stack.append((el.tag, el.attrib)); continue
        tag, at = stack.pop(); tags = [t for t, _ in stack]
        name = at.get('name', '')
        if tag == 'i' and 'parameters' in tags and name in ('LHFCALC', 'LSORBIT', 'NELECT') and name not in out['params']:
            out['params'][name] = (el.text or '').strip()
        elif tag == 'v' and name == 'divisions' and 'calculation' not in tags:
            out['divisions'] = [int(x) for x in el.text.split()]
        elif tag == 'v' and stack and stack[-1][1].get('name') == 'kpointlist':
            if 'eigenvalues_kpoints_opt' in tags:
                out['kopt'].append([float(x) for x in el.text.split()])
            elif 'calculation' not in tags:
                out['kmain'].append([float(x) for x in el.text.split()])
        elif tag == 'r' and 'eigenvalues' in tags and 'projected' not in tags and 'projected_kpoints_opt' not in tags:
            sets = [a.get('comment', '') for t, a in stack if t == 'set']
            spin = next(c for c in sets if c.startswith('spin'))
            if spin.replace(' ', '') != 'spin1':
                raise SystemExit('[ERROR] 有第二个自旋通道（%s），不是 SOC 计算' % spin)
            k = int(next(c for c in sets if c.startswith('kpoint')).split()[1])
            (eopt if 'eigenvalues_kpoints_opt' in tags else emain).setdefault(k, []).append(float(el.text.split()[0]))
        elif tag == 'v' and stack and 'structure' in tags:
            st = next(a for t, a in stack if t == 'structure')
            if st.get('name') == 'finalpos':
                vn = stack[-1][1].get('name')
                if vn == 'basis': out['basis'].append([float(x) for x in el.text.split()])
                elif vn == 'positions': out['pos'].append([float(x) for x in el.text.split()])
        elif tag == 'rc' and any(t == 'array' and a_.get('name') == 'atoms' for t, a_ in stack):
            out['symbols'].append(el[0].text.strip())
        if tag in FLUSH:
            el.clear()
    E = eopt if opt else emain
    ks = out['kopt'] if opt else out['kmain']
    if not E or len(E) != len(ks):
        raise SystemExit('[ERROR] %s：%s 本征值 %d 组、k 点 %d 个，对不上' % (path, 'KPOINTS_OPT' if opt else '自洽', len(E), len(ks)))
    out['eig'] = np.array([E[i + 1] for i in range(len(ks))]); out['kpts'] = np.array(ks)
    return out


def ir_mesh(lat, frac, symbols, mesh):
    import spglib
    nums = {s: i + 1 for i, s in enumerate(sorted(set(symbols)))}
    return spglib.get_ir_reciprocal_mesh(np.array(mesh), (lat, frac, [nums[s] for s in symbols]), is_shift=[0, 0, 0], is_time_reversal=True)


def ibz_count(lat, frac, symbols, mesh):
    try:
        return len(np.unique(ir_mesh(lat, frac, symbols, mesh)[0]))
    except ImportError:
        return None


def reduce_full(ks, E, mesh, lat, frac, symbols, nel):
    """全网格（每个网格点一个 k）→ 不可约点。返回 (ks, E, 对称等价点本征值的最大差 eV)。"""
    mesh = np.array(mesh); x = ks * mesh
    if np.abs(x - np.round(x)).max() > 1e-6:
        raise SystemExit('[ERROR] k 点不在 %s 网格上' % mesh.tolist())
    key = lambda a: tuple(np.mod(np.asarray(a, int), mesh))
    row = {key(a): i for i, a in enumerate(np.round(x).astype(int))}
    if len(row) != len(ks) or len(ks) != mesh.prod():
        raise SystemExit('[ERROR] 不是完整的全网格：%d 个 k 点（去重后 %d），网格应有 %d 个' % (len(ks), len(row), mesh.prod()))
    mapping, grid = ir_mesh(lat, frac, symbols, mesh)
    vrow = np.array([row[key(a)] for a in grid])               # spglib 第 j 个网格点 → vasprun 中的行
    rrow = vrow[mapping]                                        # 它所属的不可约代表点 → vasprun 中的行
    bands = slice(max(nel - 24, 0), min(nel + 20, E.shape[1]))   # 费米能级附近（SOC 序号）
    spread = float(np.abs(E[vrow, bands] - E[rrow, bands]).max())
    if spread > EQUIV_TOL:
        raise SystemExit('[ERROR] 对称等价点的本征值最多差 %.4f eV，对称性与数据不符' % spread)
    reps = vrow[np.unique(mapping)]
    return ks[reps], E[reps], spread


def check_ibz(ks, mesh, lat, frac, symbols):
    """每个 k 点都在网格上，且恰好每个不可约类一个。返回不可约点数；没有 spglib 时返回 None（只查点数以外的不做）。"""
    mesh = np.array(mesh); x = ks * mesh
    if np.abs(x - np.round(x)).max() > 1e-6:
        raise SystemExit('[ERROR] 有 k 点不在 %s 网格上' % mesh.tolist())
    try:
        mapping, grid = ir_mesh(lat, frac, symbols, mesh)
    except ImportError:
        return None
    gi = {tuple(np.mod(a, mesh)): j for j, a in enumerate(grid)}
    cls = [mapping[gi[tuple(np.mod(np.round(k * mesh).astype(int), mesh))]] for k in ks]
    n = len(np.unique(mapping))
    if len(set(cls)) != len(cls) or len(cls) != n:
        raise SystemExit('[ERROR] 网格 %s 的不可约点应为 %d 个、每类一个；这里 %d 个点、%d 类' % (mesh.tolist(), n, len(cls), len(set(cls))))
    return n


def merge(parts, nel):
    """多个 KPOINTS_OPT 批次：去掉重复的 k 点（锚点），并核对重复点上的本征值一致（各批自洽的一致性）。"""
    ks, E, dup = [], [], 0.0
    bands = slice(max(nel - 24, 0), nel + 20)
    for d in parts:
        for k, e in zip(d['kpts'], d['eig']):
            j = next((i for i, q in enumerate(ks) if np.abs(q - k).max() < 1e-6), None)
            if j is None:
                ks.append(k); E.append(e)
            else:
                dup = max(dup, float(np.abs(E[j][bands] - e[bands]).max()))
    if dup > ANCHOR_TOL:
        raise SystemExit('[ERROR] 不同批次在同一 k 点（锚点）上的本征值最多差 %.4f eV，各批自洽不一致' % dup)
    return np.array(ks), np.array(E), dup


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vasprun', required=True, nargs='+', help='一个或多个 vasprun.xml（多个时须同时给 --opt，各批合并）')
    ap.add_argument('--material', required=True, choices=sorted(NELECT))
    ap.add_argument('--functional', required=True, choices=('hse', 'pbe'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--sub', type=int, nargs=3, help='抽稀到的子网格，例如 14 14 7（须整除原网格）')
    ap.add_argument('--opt', action='store_true', help='读 KPOINTS_OPT 的能带（非自洽）')
    ap.add_argument('--mesh', type=int, nargs=3, help='--opt 时这些 k 点所属的均匀网格')
    a = ap.parse_args()
    if len(a.vasprun) > 1 and not a.opt:
        raise SystemExit('[ERROR] 多个 vasprun.xml 只用于 --opt（KPOINTS_OPT 分批）')

    parts = []
    for path in a.vasprun:
        print('[..] 读 %s' % path)
        d = parse(path, a.opt)
        hf = d['params'].get('LHFCALC', 'F').upper().startswith('T')
        if FUNC['hse' if hf else 'pbe'] != FUNC[a.functional] or not d['params'].get('LSORBIT', 'F').upper().startswith('T'):
            raise SystemExit('[ERROR] %s 是 %s（LSORBIT=%s），与 --functional %s 不符' %
                             (path, FUNC['hse' if hf else 'pbe'], d['params'].get('LSORBIT'), a.functional))
        nel = int(round(float(d['params']['NELECT'])))
        if nel != NELECT[a.material]:
            raise SystemExit('[ERROR] %s：NELECT = %d，%s 应为 %d' % (path, nel, a.material, NELECT[a.material]))
        if parts and (np.abs(np.array(d['basis']) - np.array(parts[0]['basis'])).max() > 1e-5 or d['symbols'] != parts[0]['symbols']):
            raise SystemExit('[ERROR] %s 的结构与第一个文件不同' % path)
        parts.append(d)
    d = parts[0]
    mesh = a.mesh if a.opt else d['divisions']
    if mesh is None:
        raise SystemExit('[ERROR] 不知道 k 点属于哪个网格：--opt 时请给 --mesh')
    ks, E, dup = merge(parts, nel) if len(parts) > 1 else (d['kpts'], d['eig'], None)
    if dup is not None:
        print('[OK] %d 批合并为 %d 个 k 点；重复点（锚点）上各批本征值最多差 %.4f eV' % (len(parts), len(ks), dup))
    if a.sub:
        if any(m % s for m, s in zip(mesh, a.sub)):
            raise SystemExit('[ERROR] 子网格 %s 不整除原网格 %s' % (a.sub, mesh))
        x = ks * np.array(a.sub)
        keep = np.all(np.abs(x - np.round(x)) < 1e-6, axis=1)
        ks, E, mesh = ks[keep], E[keep], list(a.sub)
    lat, frac, sym = np.array(d['basis']), np.array(d['pos']), d['symbols']
    full, spread = len(ks) == int(np.prod(mesh)) and len(ks) > 1, None
    n_ibz = ibz_count(lat, frac, sym, mesh)
    if full and n_ibz is None:
        raise SystemExit('[ERROR] 输入是全网格（ISYM = −1），约化到不可约点需要 spglib')
    if full and n_ibz != len(ks):
        nfull = len(ks)
        ks, E, spread = reduce_full(ks, E, mesh, lat, frac, sym, nel)
        print('[OK] 全网格 %d 个点 → 不可约点 %d 个；对称等价点本征值最多差 %.4f eV' % (nfull, len(ks), spread))
    if n_ibz is not None and n_ibz != len(ks):
        raise SystemExit('[ERROR] 网格 %s 的不可约点应为 %d 个，这里有 %d 个' % (mesh, n_ibz, len(ks)))
    check_ibz(ks, mesh, lat, frac, sym)
    grid = 'Gamma %s (%d irreducible k-points)%s%s%s' % (' x '.join(map(str, mesh)), len(ks),
                                                        ', subsampled from %s' % 'x'.join(map(str, d['divisions'])) if a.sub and not a.opt else '',
                                                        ', reduced from full grid (ISYM=-1)' if spread is not None else '',
                                                        ', KPOINTS_OPT (non-SCF%s)' % (', %d batches' % len(parts) if len(parts) > 1 else '') if a.opt else '')
    src = [os.path.abspath(p) for p in a.vasprun]; sha = [sha256(p) for p in a.vasprun]
    info = g.write_gene(a.out, a.material, FUNC[a.functional], hf, lat, sym, frac, ks, E, nel, grid,
                        src[0] if len(src) == 1 else src, sha[0] if len(sha) == 1 else sha,
                        heads=('%s %s from vasprun.xml finalpos, Bohr' % (a.material, FUNC[a.functional]),
                               '%s %s %s, one band per Kramers pair, Ry' % (a.material, FUNC[a.functional], grid)))
    print('[DONE] %s %s：%s；带隙（网格上）%.4f eV；写 %d 条能带，nelect_bt2 = %d；Kramers 最大差 %.4f eV%s → %s'
          % (a.material, info['functional'], grid, info['gap_eV'], info['nbands_written'], info['nelect_bt2'],
             info['kramers_max_split_eV'], '' if n_ibz is None else '；不可约点已用 spglib 核对（每类一个）', a.out))


if __name__ == '__main__':
    main()
