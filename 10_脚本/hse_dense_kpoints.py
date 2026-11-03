#!/usr/bin/env python3
"""HSE06+SOC 加密网格（c8_hse_dense）的 KPOINTS_OPT 分批文件（2026-11-03）。

网格与 PBE 密网格相同：Sb 系 28×28×7（515 个不可约点），Bi 系 28×28×6（450 个）。依据：网格收敛检验
（21_阶段报告与交接/交接_HSE直接BoltzTraP2_20261102.md，10_脚本/bt2_gridtest_compare.py）——面内 14、c 方向 3 都不够；
21×21×7/6 使 Sn₂Sb₂Te₅ 在霍尔浓度下的 S 偏低 17%，Sb/Bi 之比由 1.77 变为 2.14，影响问题一的结论。

不可约点由 spglib（含时间反演）从 POSCAR 求出，坐标为倒格子分数坐标，取 (−1/2, 1/2]。
所有点都检查过不在 FOCK_FCC 判据边界上（12·k_x、12·k_y 不为半整数；见 Sn2Bi2Te5/c4_mirror_v3/README.md）。
每批 36 个点：Γ 作锚点放在每批第一个（用来核对各批自洽一致），其余 35 个为新点。

用法：python 10_脚本/hse_dense_kpoints.py      写到 22_能带基准_1-0/<材料>/c8_hse_dense/（Pb₂Sb₂Te₅ 在 22_能带基准_1-0/c8_hse_dense/）
"""
import json, os, re, sys
import numpy as np
import spglib

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import eigenval_to_gene as g

MESH = {'Pb2Sb2Te5': (28, 28, 7), 'Sn2Sb2Te5': (28, 28, 7), 'Pb2Bi2Te5': (28, 28, 6), 'Sn2Bi2Te5': (28, 28, 6)}
NBANDS = {'Pb2Sb2Te5': 96, 'Sn2Sb2Te5': 96, 'Pb2Bi2Te5': 108, 'Sn2Bi2Te5': 108}   # Sb 系 84 不够（窗口上端只有 1.5 eV）
BATCH = 36


def ibz(m):
    lat, sym, frac = g.read_poscar(os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'POSCAR'))
    nums = {s: i + 1 for i, s in enumerate(sorted(set(sym)))}
    mesh = np.array(MESH[m])
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, (lat, frac, [nums[s] for s in sym]), is_shift=[0, 0, 0], is_time_reversal=True)
    k = grid[np.unique(mapping)] / mesh
    k = k - np.round(k)                                      # 取 (−1/2, 1/2]
    k[np.isclose(k, -0.5)] = 0.5
    k = k[np.lexsort((k[:, 1], k[:, 0], np.abs(k[:, 2])))]    # kz 由小到大
    x = (12 * k[:, :2]) % 1
    assert not np.any(np.isclose(x, 0.5)), 'FOCK_FCC 边界上的点'
    assert np.allclose(k[0], 0)
    return k


def incar(m, src):
    s = open(src).read()
    s = re.sub(r'^SYSTEM.*$', 'SYSTEM = %s c8_hse_dense HSE06+SOC KPOINTS_OPT %s IBZ' % (m, 'x'.join(map(str, MESH[m]))), s, flags=re.M)
    s = re.sub(r'^NBANDS\s*=.*$', 'NBANDS = %d' % NBANDS[m], s, flags=re.M)
    return s.rstrip('\n') + '\nKPOINTS_OPT_NKBATCH = %d\n' % BATCH


def write_kp(path, title, pts):
    with open(path, 'w') as f:
        f.write(title + '\n%d\nReciprocal\n' % len(pts))
        for p in pts:
            f.write('%+.14f  %+.14f  %+.14f       1\n' % tuple(p))


def main():
    for m in MESH:
        bdir = g.bench_dir(m); out = os.path.join(bdir, 'c8_hse_dense'); os.makedirs(out, exist_ok=True)
        k = ibz(m); rest = k[1:]; nb = -(-len(rest) // (BATCH - 1))
        batches = [np.vstack([k[:1], rest[i * (BATCH - 1):(i + 1) * (BATCH - 1)]]) for i in range(nb)]
        mesh = 'x'.join(map(str, MESH[m]))
        for i, b in enumerate(batches):
            write_kp(os.path.join(out, 'KPOINTS_OPT_b%02d' % (i + 1)),
                     '%s c8_hse_dense: IBZ of Gamma %s (%d points); batch %d/%d, first point = Gamma anchor; KPOINTS_OPT_NKBATCH = %d'
                     % (m, mesh, len(k), i + 1, nb, BATCH), b)
        write_kp(os.path.join(out, 'KPOINTS_OPT'), '%s c8_hse_dense: IBZ of Gamma %s, all %d points (single sequential job)' % (m, mesh, len(k)), k)
        open(os.path.join(out, 'INCAR'), 'w').write(incar(m, os.path.join(bdir, 'c2_hse', 'INCAR')))
        open(os.path.join(out, 'KPOINTS'), 'w').write(open(os.path.join(bdir, 'c2_hse', 'KPOINTS')).read())
        json.dump({'material': m, 'mesh': MESH[m], 'n_irreducible': len(k), 'batch_size': BATCH, 'n_batches': nb,
                   'anchor': 'Gamma is the first point of every batch', 'nbands': NBANDS[m],
                   'kpoints': [[round(float(x), 14) for x in p] for p in k],
                   'batch_of_point': [0] + [1 + i // (BATCH - 1) for i in range(len(rest))]},
                  open(os.path.join(out, 'c8_grid.json'), 'w'), indent=0)
        print('%s：%s，%d 个不可约点 → %d 批（每批 Γ + 35，最后一批 %d 个点）；NBANDS %d → %s'
              % (m, mesh, len(k), nb, len(batches[-1]), NBANDS[m], os.path.relpath(out, ROOT)))


if __name__ == '__main__':
    main()
