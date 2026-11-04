#!/usr/bin/env python3
"""把 1-0 能带基准的 12×12×2 自洽 EIGENVAL（HSE06+SOC、PBE+SOC）转成 BoltzTraP2 的 GENE 输入（2026-11-03）。

为什么：BoltzTraP2 的 VASP 读取器要 vasprun.xml，但 Pb₂Sb₂Te₅ 的 HSE vasprun.xml 在已到期的集群上；仓库里四种材料的
HSE、PBE EIGENVAL 都在。GENE 是 BoltzTraP2 自带的通用格式（dft.GENELoader），AutoZT 的 run_boltztrap_crta 传目录即可读。

SOC 的处理：P-3m1 有反演对称，非磁，SOC 能带处处 Kramers 二重简并（第 2j−1 与 2j 条相同）。GENE 只写每对中的一条、
nspin = 1（BoltzTraP2 据此取每态占据 2），与“全部 SOC 能带、每态占据 1”完全等价。
只写费米能级附近的能带（WIN；与 BoltzTraP2 自带的 bandana 相同的做法）：NBANDS 顶端几条高能带没收敛好、简并被破坏，
不能进插值，最顶上 SKIP_TOP 对一律不写，并检查第一条没写入的导带离费米能级足够远（GAP_ABOVE）；去掉的深能级电子从 nelect 中扣除（source.json 的 nelect_bt2）。写入的能带在能量窗口内逐 k 点检查二重简并（容差 2 meV），不满足即报错；窗口外的差只记录（2026-11-04）。

输出（每个目录）：<材料>.structure（晶格、原子笛卡尔坐标，Bohr）、<材料>.energy（Ry）、source.json（泛函、网格、来源文件
SHA-256、带边、简并检查），供 run_bt2_hse_grid.py 核对泛函后再跑。

用法：python 10_脚本/eigenval_to_gene.py      写到 22_能带基准_1-0/bt2_gene/<材料>/{hse,pbe}/
"""
import hashlib, json, os, re
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCH = os.path.join(ROOT, '22_能带基准_1-0')
OUT = os.path.join(BENCH, 'bt2_gene')
RY, BOHR = 13.605693122994, 0.529177210903          # eV/Ry，Å/Bohr
MATS = ('Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5')
KRAMERS_TOL = 2e-3                                   # eV；EIGENVAL 只有 4 位小数
WIN = (3.0, 2.5)                                     # eV，写入费米能级下 3.0、上 2.5 eV 内的能带（化学势只扫 ±1.2 eV）
SKIP_TOP = 2                                         # NBANDS 最顶上的 2 对能带一律不写（迭代对角化时最不收敛）
GAP_ABOVE = 1.8                                      # 第一条没写入的导带，最低点须高于费米能级这么多（eV）


def bench_dir(m):
    return BENCH if m == 'Pb2Sb2Te5' else os.path.join(BENCH, m)    # PbSb 的数据在基准目录顶层


def sha256(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def read_poscar(p):
    L = open(p).read().splitlines(); s = float(L[1])
    lat = s * np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])
    el, n = L[5].split(), [int(x) for x in L[6].split()]
    frac = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(sum(n))])
    assert L[7].strip()[0] in 'Dd', 'POSCAR 应为 Direct 坐标'
    return lat, [e for e, k in zip(el, n) for _ in range(k)], frac


def read_eigenval(p):
    L = open(p).read().splitlines()
    nel, nk, nb = (int(x) for x in L[5].split()[:3])
    i, ks, E = 7, [], []
    for _ in range(nk):
        while not L[i].strip(): i += 1
        ks.append([float(x) for x in L[i].split()[:3]]); i += 1
        E.append([float(L[i + b].split()[1]) for b in range(nb)]); i += nb
    return nel, np.array(ks), np.array(E)                   # E: (nk, nb)，eV


def incar_flag(p, key):
    m = re.search(r'^\s*%s\s*=\s*\.?(\w+)' % key, open(p).read(), re.I | re.M)
    return bool(m) and m.group(1).upper().startswith('T')


def convert(m, kind):
    src = os.path.join(bench_dir(m), {'hse': 'c1_HSE_SOC', 'pbe': 'b_PBE_SOC'}[kind])
    ev, inc, kp = (os.path.join(src, f) for f in ('EIGENVAL', 'INCAR', 'KPOINTS'))
    hf = incar_flag(inc, 'LHFCALC')
    if hf != (kind == 'hse') or not incar_flag(inc, 'LSORBIT'):
        raise SystemExit('[ERROR] %s：INCAR 的 LHFCALC=%s、LSORBIT 与 %s 不符' % (src, hf, kind))
    mesh = open(kp).read().splitlines()[3].split()[:3]
    lat, sym, frac = read_poscar(os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'POSCAR'))
    a_e, c_e = (float(x) * 1e10 for x in open(ev).read().splitlines()[1].split()[1:4:2])
    if abs(np.linalg.norm(lat[0]) - a_e) > 1e-4 or abs(np.linalg.norm(lat[2]) - c_e) > 1e-4:
        raise SystemExit('[ERROR] %s：POSCAR 晶格与 EIGENVAL 不符' % m)
    nel, ks, E = read_eigenval(ev)
    return write_gene(os.path.join(OUT, m, kind), m, {'hse': 'HSE06+SOC', 'pbe': 'PBE+SOC'}[kind], hf, lat, sym, frac, ks, E, nel,
                      'Gamma %s (%d irreducible k-points), 1-0 benchmark SCF' % (' x '.join(mesh), len(ks)),
                      os.path.relpath(ev, ROOT), sha256(ev),
                      heads=('%s %s+SOC 1-0 benchmark (POSCAR of A2B2Te5_SOC_dense), Bohr' % (m, kind.upper()),
                             '%s %s+SOC Gamma %s, one band per Kramers pair, Ry' % (m, kind.upper(), 'x'.join(mesh))))


def write_gene(out, m, functional, hf, lat, sym, frac, ks, E, nel, grid, source, source_sha, heads):
    """写 GENE 目录：每对 Kramers 简并能带写一条（nspin = 1），只写费米能级附近、收敛的能带。
    lat (3,3) Å；frac (natom,3) 分数坐标；ks (nk,3) 倒格子分数坐标；E (nk, nband) eV，nband 为全部 SOC 能带。
    heads：.structure、.energy 两个文件的标题行。"""
    E = np.asarray(E, float); ks = np.asarray(ks, float)
    assert E.shape[1] % 2 == 0 and nel % 2 == 0
    vbm, cbm = E[:, nel - 1].max(), E[:, nel].min()
    if cbm <= vbm:
        raise SystemExit('[ERROR] %s：带隙为负（%.4f eV），不能用半满中点作费米能级' % (m, cbm - vbm))
    ef = (vbm + cbm) / 2
    P = E[:, 0::2]                                           # 每对取一条，(nk, npair)
    keep = [j for j in range(P.shape[1] - SKIP_TOP) if P[:, j].max() >= ef - WIN[0] and P[:, j].min() <= ef + WIN[1]]
    lo, hi = keep[0], keep[-1]
    assert keep == list(range(lo, hi + 1)) and lo < nel // 2 <= hi
    above = P[:, hi + 1].min() - ef                           # 第一条没写入的导带
    if above < GAP_ABOVE:
        raise SystemExit('[ERROR] %s：第一条没写入的导带只高出费米能级 %.2f eV' % (m, above))
    # Kramers 二重简并只在能量窗口内检查（2026-11-04）：能带按整条取舍，一条带只要有一部分进窗口就整条写入，
    # 它在窗口外（远离费米能级）的部分不影响输运，那里的本征值又最接近 NBANDS 顶端、迭代对角化最不收敛。
    # 容差 KRAMERS_TOL 不变；窗口外的最大差照样记进 source.json（kramers_max_split_all_eV 及其位置）。
    A_, B_ = E[:, 2 * lo:2 * hi + 2:2], E[:, 2 * lo + 1:2 * hi + 2:2]
    dk = np.abs(A_ - B_); inwin = (np.minimum(A_, B_) >= ef - WIN[0]) & (np.maximum(A_, B_) <= ef + WIN[1])
    split = dk[inwin].max(); split_all = dk.max()
    ik, jb = np.unravel_index(int(np.argmax(dk)), dk.shape)
    where_all = {'k_index_0based': int(ik), 'soc_bands_1based': [int(2 * (lo + jb) + 1), int(2 * (lo + jb) + 2)],
                 'E_minus_EF_eV': round(float((A_[ik, jb] + B_[ik, jb]) / 2 - ef), 3)}
    if split > KRAMERS_TOL:
        raise SystemExit('[ERROR] %s：窗口内的能带 Kramers 二重简并不满足（最大差 %.4f eV）' % (m, split))
    if split_all > KRAMERS_TOL:
        print('[NOTE] %s：窗口外有 Kramers 差 %.4f eV（k #%d，SOC 能带 %d/%d，E − E_F = %+.2f eV），不影响输运，照常写入'
              % (m, split_all, where_all['k_index_0based'], *where_all['soc_bands_1based'], where_all['E_minus_EF_eV']))
    Eh = P[:, lo:hi + 1]
    os.makedirs(out, exist_ok=True)
    cart = frac @ lat
    with open(os.path.join(out, m + '.structure'), 'w') as f:
        f.write(heads[0] + '\n')
        for v in lat: f.write('%.10f %.10f %.10f\n' % tuple(v / BOHR))
        f.write('%d\n' % len(sym))
        for s_, r in zip(sym, cart): f.write('%s %.10f %.10f %.10f\n' % ((s_,) + tuple(r / BOHR)))
    with open(os.path.join(out, m + '.energy'), 'w') as f:
        f.write(heads[1] + '\n')
        f.write('%d 1 %.10f\n' % (len(ks), ef / RY))
        for k, e in zip(ks, Eh):
            f.write('%.12f %.12f %.12f %d\n' % (k[0], k[1], k[2], len(e)))
            for x in e: f.write('%.10f\n' % (x / RY))
    info = {'material': m, 'functional': functional, 'LHFCALC': hf, 'grid': grid,
            'nelect': nel, 'nelect_bt2': nel - 2 * lo, 'nbands_soc': int(E.shape[1]), 'nbands_written': int(Eh.shape[1]),
            'soc_bands_written_1based': [2 * lo + 1, 2 * hi + 2], 'window_below_above_eV': list(WIN),
            'first_unwritten_cb_min_above_EF_eV': round(float(above), 3),
            'vbm_band_written_1based': nel // 2 - lo, 'kramers_max_split_eV': round(float(split), 5),
            'kramers_check': 'within the energy window only (EF - %.1f to EF + %.1f eV), tolerance %.1f meV' % (WIN[0], WIN[1], 1000 * KRAMERS_TOL),
            'kramers_max_split_all_eV': round(float(split_all), 5), 'kramers_max_split_all_at': where_all,
            'vbm_eV': float(vbm), 'cbm_eV': float(cbm), 'gap_eV': round(float(cbm - vbm), 4),
            'source_file': source, 'source_sha256': source_sha,
            'note': 'GENE nspin=1（每态占据 2）只写每对 Kramers 简并能带中的一条，等价于全部 SOC 能带、每态占据 1。'
                    'BoltzTraP2 的电子数用 nelect_bt2（已扣除未写入的深能级）；价带顶是写入的第 vbm_band_written_1based 条'}
    json.dump(info, open(os.path.join(out, 'source.json'), 'w'), ensure_ascii=False, indent=1)
    return info


def main():
    for m in MATS:
        for kind in ('hse', 'pbe'):
            d = convert(m, kind)
            print('%-10s %-9s 带隙 %.4f eV，%s；SOC 能带 %d 条，写第 %d–%d 条中每对一条（%d 条），nelect %d → %d，Kramers 最大差 %.4f eV' %
                  (m, d['functional'], d['gap_eV'], d['grid'].split(', ')[0], d['nbands_soc'], *d['soc_bands_written_1based'],
                   d['nbands_written'], d['nelect'], d['nelect_bt2'], d['kramers_max_split_eV']))


if __name__ == '__main__':
    main()
