#!/usr/bin/env python3
"""在 1-0 基准的 12×12×2 自洽网格上跑 BoltzTraP2（CRTA），HSE06+SOC 直接用，PBE+SOC 作同网格对照（2026-11-02）。

不改 AutoZT：直接调用 AutoZT step8.1 的 gen_step11_boltztrap.run_boltztrap_crta（与密网格结果同一函数），
只是换输入目录、温度、插值倍数，并把结果写成精简 JSON（paper_figures.bt2_signed 能直接读）。
交接说明见 21_阶段报告与交接/交接_HSE直接BoltzTraP2_20261102.md。

用法（在集群上，需要 BoltzTraP2 环境）：
  python run_bt2_hse_grid.py --autozt <含 gen_step11_boltztrap.py 的目录> --material Sn2Sb2Te5 \
      --functional hse --run-dir <含 vasprun.xml 的作业目录> --mult 20 --out <输出目录> \
      [--kcheck <仓库>/22_能带基准_1-0/Sn2Sb2Te5/bt2_check_kpoints.txt]

输出（--out 目录下）：
  boltztrap_crta.json   T、化学势、带符号载流子浓度（+ 空穴）、面内 S 与 σ/τ 的 xx、yy 分量
  bands_check.txt       （给了 --kcheck 时）插值能带在检验 k 点上的能量，格式同 c2_hse/path_eigen.txt
"""
import argparse, hashlib, json, sys
from pathlib import Path

HA = 27.211386245988                     # eV/Hartree
# PBE 同网格对照用的剪刀目标带隙：与密网格 step8.1（run_bt2_soc.py 的 GAPS）相同，使两者只差 k 网格
PBE_GAPS = {'Pb2Bi2Te5': 0.1141, 'Pb2Sb2Te5': 0.0654, 'Sn2Bi2Te5': 0.1366, 'Sn2Sb2Te5': 0.1352}
NELECT = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
TEMPS = list(range(300, 651, 25))        # 实验点在 300–620 K
KEEP = ('seebeck_xx_V_per_K', 'seebeck_yy_V_per_K', 'sigma_over_tau_xx', 'sigma_over_tau_yy', 'carrier_conc_cm-3')


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def slim(x):
    """嵌套列表保留 6 位有效数字，控制文件大小；NaN 写成 null。"""
    if isinstance(x, list):
        return [slim(v) for v in x]
    return None if x is None or x != x else float('%.6g' % x)


def read_kcheck(p):
    L = Path(p).read_text().splitlines(); n = int(L[1].split()[0])
    return [[float(v) for v in l.split()[:3]] for l in L[3:3 + n]]


def band_check(run_dir, scissor, mult, kp, nelect, bt, out):
    """用同样的插值（同一剪刀、同一倍数）求检验 k 点上的能带，写价带顶下 4 条、导带底上 3 条。"""
    import numpy as np
    from BoltzTraP2 import dft, sphere, fite
    data = dft.DFTData(str(run_dir))
    bt._apply_scissor(data, scissor, HA)
    eq = sphere.get_equivalences(data.atoms, data.magmom, len(data.kpoints) * mult)
    coeffs = fite.fitde3D(data, eq)
    eg = fite.getBands(np.array(kp), eq, data.get_lattvec(), coeffs)[0] * HA     # (nband, nk)，eV
    lines = ['# k_index band_index energy_eV（BoltzTraP2 插值；band_index 从 1 数，%d 为价带顶那条）' % nelect]
    for i in range(len(kp)):
        for b in range(nelect - 4, nelect + 3):
            lines.append('%d %d %.6f' % (i, b + 1, eg[b, i]))
    (out / 'bands_check.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--autozt', required=True, help='含 gen_step11_boltztrap.py 的目录（AutoZT step8.1）')
    ap.add_argument('--material', required=True, choices=sorted(PBE_GAPS))
    ap.add_argument('--functional', required=True, choices=('hse', 'pbe'))
    ap.add_argument('--run-dir', required=True, help='12×12×2 自洽作业目录（含 vasprun.xml）')
    ap.add_argument('--mult', type=int, required=True, help='插值倍数 kmesh_mult（5、20、60 各跑一次）')
    ap.add_argument('--out', required=True)
    ap.add_argument('--kcheck', help='检验 k 点文件（KPOINTS_OPT 格式）')
    a = ap.parse_args()

    run = Path(a.run_dir).resolve(); out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    vr = run / 'vasprun.xml'
    if not vr.is_file():
        raise SystemExit('[ERROR] %s 没有 vasprun.xml' % run)
    sys.path.insert(0, str(Path(a.autozt).resolve()))
    import gen_step11_boltztrap as bt
    bt.PAPER_SCAN = False                 # 与 run_bt2_soc.py 相同：只跑 CRTA，不并入 DPT
    scissor = None if a.functional == 'hse' else PBE_GAPS[a.material]

    res = bt.run_boltztrap_crta(run, temperatures=TEMPS, mu_window_ev=1.2, mu_npts=600, dos_npts=6000,
                                kmesh_mult=a.mult, nelect=None, is_2d=False, scissor_gap_eV=scissor)
    if res.get('nelect') is not None and int(round(res['nelect'])) != NELECT[a.material]:
        raise SystemExit('[ERROR] NELECT = %s，应为 %d：作业目录不对' % (res['nelect'], NELECT[a.material]))
    doc = {k: res[k] for k in ('temperatures_K', 'nelect', 'cell_volume_cm3', 'scissor')}
    doc['mu_rel_efermi_eV'] = res['mu_rel_efermi_eV']
    doc.update({k: slim(res[k]) for k in KEEP})
    doc.update({'material': a.material, 'functional': {'hse': 'HSE06+SOC', 'pbe': 'PBE+SOC'}[a.functional],
                'grid': 'Gamma 12x12x2 (38 irreducible k-points), 1-0 benchmark SCF', 'kmesh_mult': a.mult,
                'input_dir': str(run), 'input_sha256': {'vasprun.xml': sha256(vr)},
                'protocol': 'BoltzTraP2 CRTA via AutoZT run_boltztrap_crta; '
                            + ('no scissor (HSE gap as computed)' if scissor is None else 'scissor to dense-grid target gap %.4f eV' % scissor),
                'note': 'S 与 σ/τ 取面内 xx、yy；载流子浓度 + 为空穴、− 为电子（cm⁻³）'})
    (out / 'boltztrap_crta.json').write_text(json.dumps(doc, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    print('[DONE] %s %s mult=%d → %s' % (a.material, a.functional, a.mult, out / 'boltztrap_crta.json'))

    if a.kcheck:
        try:
            band_check(run, scissor, a.mult, read_kcheck(a.kcheck), NELECT[a.material], bt, out)
            print('[DONE] bands_check.txt')
        except Exception as e:                # 检验失败不影响主结果
            print('[WARN] 能带检验没做成：%r' % (e,))


if __name__ == '__main__':
    main()
