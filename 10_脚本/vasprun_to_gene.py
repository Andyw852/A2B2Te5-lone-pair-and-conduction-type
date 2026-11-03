#!/usr/bin/env python3
"""从 VASP 的 vasprun.xml（可达数 GB）生成 BoltzTraP2 的 GENE 输入，可抽稀到较疏的子网格（2026-11-03）。

两种用途：
1. 网格收敛检验（PBE+SOC）：把 28×28×7/6 密网格抽成 14×14×7 等子网格，与完整密网格比较 BoltzTraP2 的 S，
   找出够用的最小网格，再按它算 HSE。
2. 以后 HSE06+SOC 的 KPOINTS_OPT 非自洽结果（--opt）：读 <eigenvalues_kpoints_opt> 里的能带。

流式读取（不把整个文件读进内存），只取：参数（LHFCALC、LSORBIT、NELECT）、自洽网格 divisions、
k 点表、自洽或 KPOINTS_OPT 的本征值、末态结构（finalpos）与元素。GENE 的写法与 eigenval_to_gene.write_gene 相同
（每对 Kramers 简并能带写一条、只写费米能级附近、扣除深能级电子），所以 run_bt2_hse_grid.py 可以直接跑。

抽稀（--sub nx ny nz）：只保留 k·n 为整数的不可约点。Γ 中心网格的子格在对称操作下不变，保留下来的正好是子网格的
一套不可约点；装了 spglib 时再用它数一遍核对。

用法：
  python 10_脚本/vasprun_to_gene.py --vasprun <.../vasprun.xml> --material Sn2Sb2Te5 --functional pbe \
      --sub 14 14 7 --out 22_能带基准_1-0/bt2_gridtest/Sn2Sb2Te5/pbe_14x14x7
  不给 --sub 即用完整网格；--opt 读 KPOINTS_OPT 的能带，此时用 --mesh 写明这些点属于哪个均匀网格。
"""
import argparse, hashlib, os, sys
import xml.etree.ElementTree as ET
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eigenval_to_gene as g

FUNC = {'hse': 'HSE06+SOC', 'pbe': 'PBE+SOC'}
NELECT = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
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


def ibz_count(lat, frac, symbols, mesh):
    try:
        import spglib
    except ImportError:
        return None
    nums = {s: i + 1 for i, s in enumerate(sorted(set(symbols)))}
    mapping, _ = spglib.get_ir_reciprocal_mesh(mesh, (lat, frac, [nums[s] for s in symbols]), is_shift=[0, 0, 0], is_time_reversal=True)
    return len(np.unique(mapping))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vasprun', required=True)
    ap.add_argument('--material', required=True, choices=sorted(NELECT))
    ap.add_argument('--functional', required=True, choices=('hse', 'pbe'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--sub', type=int, nargs=3, help='抽稀到的子网格，例如 14 14 7（须整除原网格）')
    ap.add_argument('--opt', action='store_true', help='读 KPOINTS_OPT 的能带（非自洽）')
    ap.add_argument('--mesh', type=int, nargs=3, help='--opt 时这些 k 点所属的均匀网格')
    a = ap.parse_args()

    print('[..] 读 %s' % a.vasprun)
    d = parse(a.vasprun, a.opt)
    hf = d['params'].get('LHFCALC', 'F').upper().startswith('T')
    if FUNC['hse' if hf else 'pbe'] != FUNC[a.functional] or not d['params'].get('LSORBIT', 'F').upper().startswith('T'):
        raise SystemExit('[ERROR] 这个 vasprun.xml 是 %s（LSORBIT=%s），与 --functional %s 不符' %
                         (FUNC['hse' if hf else 'pbe'], d['params'].get('LSORBIT'), a.functional))
    nel = int(round(float(d['params']['NELECT'])))
    if nel != NELECT[a.material]:
        raise SystemExit('[ERROR] NELECT = %d，%s 应为 %d' % (nel, a.material, NELECT[a.material]))
    mesh = a.mesh if a.opt else d['divisions']
    if mesh is None:
        raise SystemExit('[ERROR] 不知道 k 点属于哪个网格：--opt 时请给 --mesh')
    ks, E = d['kpts'], d['eig']
    if a.sub:
        if any(m % s for m, s in zip(mesh, a.sub)):
            raise SystemExit('[ERROR] 子网格 %s 不整除原网格 %s' % (a.sub, mesh))
        x = ks * np.array(a.sub)
        keep = np.all(np.abs(x - np.round(x)) < 1e-6, axis=1)
        ks, E, mesh = ks[keep], E[keep], list(a.sub)
    lat, frac, sym = np.array(d['basis']), np.array(d['pos']), d['symbols']
    n_ibz = ibz_count(lat, frac, sym, mesh)
    if n_ibz is not None and n_ibz != len(ks):
        raise SystemExit('[ERROR] 网格 %s 的不可约点应为 %d 个，这里有 %d 个' % (mesh, n_ibz, len(ks)))
    grid = 'Gamma %s (%d irreducible k-points)%s%s' % (' x '.join(map(str, mesh)), len(ks),
                                                      ', subsampled from %s' % 'x'.join(map(str, d['divisions'])) if a.sub and not a.opt else '',
                                                      ', KPOINTS_OPT (non-SCF)' if a.opt else '')
    info = g.write_gene(a.out, a.material, FUNC[a.functional], hf, lat, sym, frac, ks, E, nel, grid,
                        os.path.abspath(a.vasprun), sha256(a.vasprun),
                        heads=('%s %s from vasprun.xml finalpos, Bohr' % (a.material, FUNC[a.functional]),
                               '%s %s %s, one band per Kramers pair, Ry' % (a.material, FUNC[a.functional], grid)))
    print('[DONE] %s %s：%s；带隙（网格上）%.4f eV；写 %d 条能带，nelect_bt2 = %d；Kramers 最大差 %.4f eV%s → %s'
          % (a.material, info['functional'], grid, info['gap_eV'], info['nbands_written'], info['nelect_bt2'],
             info['kramers_max_split_eV'], '' if n_ibz is None else '；不可约点数已用 spglib 核对', a.out))


if __name__ == '__main__':
    main()
