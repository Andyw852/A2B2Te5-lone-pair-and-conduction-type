#!/usr/bin/env python3
"""用 BoltzTraP2（CRTA）直接算 HSE06+SOC 能带的输运；PBE+SOC 同网格作对照（2026-11-02，11-03 改为读 GENE、核对泛函）。

不改 AutoZT：直接调用 AutoZT step8.1 的 gen_step11_boltztrap.run_boltztrap_crta（与密网格结果同一函数），
只是换输入目录、温度、插值倍数，并把结果写成精简 JSON（paper_figures.bt2_signed 能直接读）。
交接说明见 21_阶段报告与交接/交接_HSE直接BoltzTraP2_20261102.md。

输入目录（--run-dir）两种：
  - GENE 目录（推荐）：22_能带基准_1-0/bt2_gene/<材料>/{hse,pbe}/，由 10_脚本/eigenval_to_gene.py 从仓库里 1-0 基准的
    12×12×2 自洽 EIGENVAL 生成；泛函、网格、电子数写在同目录的 source.json。
  - VASP 作业目录（含 vasprun.xml）：泛函从 vasprun.xml 的 LHFCALC 读，网格从 KPOINTS 生成信息读。
两种都先核对：实际泛函必须与 --functional 一致，否则报错退出（11-02 的第一轮误把 PBE 密网格 vasprun.xml 当成 HSE）。

用法（在集群上，需要 BoltzTraP2 环境）：
  python 10_脚本/run_bt2_hse_grid.py --autozt <含 gen_step11_boltztrap.py 的目录> --material Sn2Sb2Te5 \
      --functional hse --run-dir 22_能带基准_1-0/bt2_gene/Sn2Sb2Te5/hse --mult 20 \
      --out 22_能带基准_1-0/bt2_gene/Sn2Sb2Te5/hse/bt2_mult20 \
      --kcheck 22_能带基准_1-0/Sn2Sb2Te5/bt2_check_kpoints.txt

输出（--out 目录下）：
  boltztrap_crta.json   T、化学势、带符号载流子浓度（+ 空穴）、面内 S 与 σ/τ 的 xx、yy 分量，以及核对过的泛函、网格、输入指纹
  bands_check.txt       （给了 --kcheck 时）插值能带在检验 k 点上的能量，格式同 c2_hse/path_eigen.txt（SOC 能带序号）
"""
import argparse, hashlib, json, re, sys
from pathlib import Path

HA = 27.211386245988                     # eV/Hartree
# PBE 同网格对照用的剪刀目标带隙：与密网格 step8.1（run_bt2_soc.py 的 GAPS）相同，使两者只差 k 网格
PBE_GAPS = {'Pb2Bi2Te5': 0.1141, 'Pb2Sb2Te5': 0.0654, 'Sn2Bi2Te5': 0.1366, 'Sn2Sb2Te5': 0.1352}
NELECT = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}
FUNC = {'hse': 'HSE06+SOC', 'pbe': 'PBE+SOC'}
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


def inspect_input(run, material):
    """识别输入并核对：返回 dict(functional, grid, nelect_bt2, soc_row, fingerprint)。
    soc_row(b) 把 SOC 能带序号 b（从 1 数）换成 BoltzTraP2 能带数组的行号。"""
    src = run / 'source.json'
    if src.is_file():                                          # GENE
        s = json.loads(src.read_text(encoding='utf-8'))
        if s['material'] != material:
            raise SystemExit('[ERROR] %s 是 %s 的数据，不是 %s' % (run, s['material'], material))
        energy = next(run.glob('*.energy'))
        srcf = s.get('source_file', s.get('source_eigenval'))                  # 来源文件（11-03 前的 source.json 键名为 source_eigenval）
        eig = Path(__file__).resolve().parents[1] / srcf                        # 在仓库里运行时，核对来源文件没变（绝对路径则原样）
        if eig.is_file() and sha256(eig) != s['source_sha256']:
            raise SystemExit('[ERROR] %s 与生成 GENE 时的指纹不同，先重跑 eigenval_to_gene.py' % eig)
        lo = (s['soc_bands_written_1based'][0] - 1) // 2
        return {'functional': s['functional'], 'grid': s['grid'], 'nelect_bt2': s['nelect_bt2'], 'format': 'GENE',
                'soc_row': lambda b: (b - 1) // 2 - lo,
                'fingerprint': {energy.name: sha256(energy), 'source': srcf, 'source_sha256': s['source_sha256'],
                                'gap_on_grid_eV': s['gap_eV']}}
    vr = run / 'vasprun.xml'
    if not vr.is_file():
        raise SystemExit('[ERROR] %s 既不是 GENE 目录（没有 source.json），也没有 vasprun.xml' % run)
    head = ''
    with open(vr, encoding='utf-8', errors='ignore') as f:                  # incar、kpoints、parameters 都在文件开头
        while '</parameters>' not in head:
            chunk = f.read(1 << 22)
            if not chunk:
                raise SystemExit('[ERROR] %s 里找不到 <parameters>' % vr)
            head += chunk
    hf = re.search(r'name="LHFCALC"\s*>\s*T', head) is not None
    div = re.search(r'name="divisions"\s*>\s*([\d\s]+)<', head)
    return {'functional': FUNC['hse' if hf else 'pbe'], 'nelect_bt2': NELECT[material], 'format': 'VASP',
            'grid': 'Gamma %s (vasprun.xml)' % ('x'.join(div.group(1).split()) if div else '?'),
            'soc_row': lambda b: b - 1, 'fingerprint': {'vasprun.xml': sha256(vr)}}


def band_check(run_dir, scissor, mult, kp, nelect, soc_row, bt, out):
    """用同样的插值（同一剪刀、同一倍数）求检验 k 点上的能带：SOC 序号 nelect−3 … nelect+4（价带顶为 nelect）。"""
    import numpy as np
    from BoltzTraP2 import dft, sphere, fite
    data = dft.DFTData(str(run_dir))
    bt._apply_scissor(data, scissor, HA)
    eq = sphere.get_equivalences(data.atoms, data.magmom, len(data.kpoints) * mult)
    coeffs = fite.fitde3D(data, eq)
    eg = fite.getBands(np.array(kp), eq, data.get_lattvec(), coeffs)[0] * HA     # (nband, nk)，eV
    lines = ['# k_index band_index energy_eV（BoltzTraP2 插值；band_index 为 SOC 能带序号，从 1 数，%d 为价带顶）' % nelect]
    for i in range(len(kp)):
        for b in range(nelect - 3, nelect + 5):
            lines.append('%d %d %.6f' % (i, b, eg[soc_row(b), i]))
    (out / 'bands_check.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--autozt', required=True, help='含 gen_step11_boltztrap.py 的目录（AutoZT step8.1）')
    ap.add_argument('--material', required=True, choices=sorted(PBE_GAPS))
    ap.add_argument('--functional', required=True, choices=('hse', 'pbe'))
    ap.add_argument('--run-dir', required=True, help='GENE 目录（22_能带基准_1-0/bt2_gene/<材料>/<hse|pbe>）或含 vasprun.xml 的作业目录')
    ap.add_argument('--mult', type=int, required=True, help='插值倍数 kmesh_mult（5、20、60 各跑一次）')
    ap.add_argument('--out', required=True)
    ap.add_argument('--kcheck', help='检验 k 点文件（KPOINTS_OPT 格式）')
    a = ap.parse_args()

    run = Path(a.run_dir).resolve(); out = Path(a.out).resolve()
    info = inspect_input(run, a.material)
    if info['functional'] != FUNC[a.functional]:
        raise SystemExit('[ERROR] --functional %s，但 %s 实际是 %s。停止，不要改标签继续跑。' % (a.functional, run, info['functional']))
    print('[OK] 输入核对：%s，%s，%s，BoltzTraP2 电子数 %d' % (info['format'], info['functional'], info['grid'], info['nelect_bt2']))
    out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(a.autozt).resolve()))
    import gen_step11_boltztrap as bt
    bt.PAPER_SCAN = False                 # 与 run_bt2_soc.py 相同：只跑 CRTA，不并入 DPT
    scissor = None if a.functional == 'hse' else PBE_GAPS[a.material]

    res = bt.run_boltztrap_crta(run, temperatures=TEMPS, mu_window_ev=1.2, mu_npts=600, dos_npts=6000,
                                kmesh_mult=a.mult, nelect=info['nelect_bt2'], is_2d=False, scissor_gap_eV=scissor)
    doc = {k: res[k] for k in ('temperatures_K', 'nelect', 'cell_volume_cm3', 'scissor')}
    doc['mu_rel_efermi_eV'] = res['mu_rel_efermi_eV']
    doc.update({k: slim(res[k]) for k in KEEP})
    doc.update({'material': a.material, 'functional': info['functional'], 'functional_checked': True,
                'grid': info['grid'], 'input_format': info['format'], 'kmesh_mult': a.mult,
                'input_dir': str(run), 'input_fingerprint': info['fingerprint'],
                'protocol': 'BoltzTraP2 CRTA via AutoZT run_boltztrap_crta; '
                            + ('no scissor (HSE gap as computed)' if scissor is None else 'scissor to dense-grid target gap %.4f eV' % scissor),
                'note': 'S 与 σ/τ 取面内 xx、yy；载流子浓度 + 为空穴、− 为电子（cm⁻³）；GENE 输入时 nelect 已扣除未写入的深能级'})
    (out / 'boltztrap_crta.json').write_text(json.dumps(doc, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    print('[DONE] %s %s mult=%d → %s' % (a.material, info['functional'], a.mult, out / 'boltztrap_crta.json'))

    if a.kcheck:
        try:
            band_check(run, scissor, a.mult, read_kcheck(a.kcheck), NELECT[a.material], info['soc_row'], bt, out)
            print('[DONE] bands_check.txt')
        except Exception as e:                # 检验失败不影响主结果
            print('[WARN] 能带检验没做成：%r' % (e,))


if __name__ == '__main__':
    main()
