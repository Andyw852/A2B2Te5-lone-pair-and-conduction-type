"""HSE06+SOC 直接 BoltzTraP2（12×12×2 GENE 输入，分支 bt2-hse-gene）的三项验收（2026-11-03）。

1. 插值倍数收敛：mult 5、20、60 的 S（300 K，±1e20、±1e21 cm⁻³）。
2. 网格够不够：PBE+SOC 12×12×2（剪刀同密网格）对 PBE+SOC 28×28×7/6 密网格（A2B2Te5_SOC_dense step8.1），300 K。
3. 插值能带：bands_check.txt 对直接计算的 HSE 路径（c2_hse/path_eigen.txt）与 Sn 系 c5 真实价带谷。

用法：python 10_脚本/bt2_hse_validate.py      只读，不写文件。
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
BENCH = os.path.join(ROOT, '22_能带基准_1-0')
MATS = ('Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5')
NELECT = {'Pb2Sb2Te5': 68, 'Sn2Sb2Te5': 68, 'Pb2Bi2Te5': 88, 'Sn2Bi2Te5': 88}


def bench_dir(m):
    return BENCH if m == 'Pb2Sb2Te5' else os.path.join(BENCH, m)


def gene_json(m, kind, mult=60):
    return os.path.join(BENCH, 'bt2_gene', m, kind, 'bt2_mult%d' % mult, 'boltztrap_crta.json')


def dense_json(m):
    return os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'ke-dft-cpu', 'result', 'step8.1_boltztrap', 'boltztrap_crta.json')


def s_at(path, N, T):
    """面内 S（µV/K）在带符号浓度 N（+ 空穴）、温度 T（须为表中温度）处，按 log|N| 线性插值。"""
    d = json.load(open(path)); it = list(d['temperatures_K']).index(float(T))
    c = np.array(d['carrier_conc_cm-3'][it], float)
    S = (np.array(d['seebeck_xx_V_per_K'][it], float) + np.array(d['seebeck_yy_V_per_K'][it], float)) / 2 * 1e6
    sg = 1 if N > 0 else -1; ok = (sg * c > 0) & np.isfinite(S); x = np.log10(sg * c[ok]); o = np.argsort(x)
    if not x.min() <= np.log10(abs(N)) <= x.max():
        raise ValueError('浓度 %.1e 超出范围' % N)
    return float(np.interp(np.log10(abs(N)), x[o], S[ok][o]))


def read_path(p):
    d = np.loadtxt(p, comments='#'); E = {}
    for k, b, e in d: E[(int(k), int(b))] = e
    return E


def main():
    print('# 1 插值倍数收敛（300 K，面内 S，µV/K；mult 5 / 20 / 60，最大相对差 |20−60|）')
    worst = 0
    for m in MATS:
        for kind in ('hse', 'pbe'):
            cells = []
            for N in (1e20, 1e21, -1e20, -1e21):
                v = [s_at(gene_json(m, kind, k), N, 300) for k in (5, 20, 60)]
                r = abs(v[1] - v[2]) / abs(v[2]); worst = max(worst, r)
                cells.append('%+.0e: %.1f/%.1f/%.1f' % (N, *v))
            print('  %-10s %s  %s' % (m, kind, '  '.join(cells)))
    print('  |mult20 − mult60| 最大 %.2f%%' % (worst * 100))

    print('\n# 2 网格：PBE+SOC 12×12×2（剪刀同密网格）对 28×28×7/6 密网格，300 K，面内 S（µV/K）')
    worst = 0
    for m in MATS:
        row = []
        for N in (3e19, 1e20, 3e20, 1e21, 2e21, -3e19, -1e20, -3e20, -1e21, -2e21):
            a, b = s_at(gene_json(m, 'pbe'), N, 300), s_at(dense_json(m), N, 300)
            rel = (a - b) / abs(b); worst = max(worst, abs(rel)); row.append('%+.0e %.0f/%.0f (%+.0f%%)' % (N, a, b, rel * 100))
        print('  %s:\n    %s\n    %s' % (m, '  '.join(row[:5]), '  '.join(row[5:])))
    print('  最大相对差 %.1f%%' % (worst * 100))

    print('\n# 3 插值能带（HSE，mult 60）对直接计算：路径 204 点，价带顶两条、导带底两条；能量都相对各自的价带顶')
    for m in MATS:
        nel = NELECT[m]
        P = read_path(os.path.join(bench_dir(m), 'c2_hse', 'path_eigen.txt'))
        B = read_path(os.path.join(BENCH, 'bt2_gene', m, 'hse', 'bt2_mult60', 'bands_check.txt'))
        vb_p = max(P[(k, nel)] for k in range(204)); vb_b = max(B[(k, nel)] for k in range(204))
        cb_p = min(P[(k, nel + 1)] for k in range(204)); cb_b = min(B[(k, nel + 1)] for k in range(204))
        out = []
        for b, nm in ((nel - 2, 'VB−1'), (nel, 'VB'), (nel + 1, 'CB'), (nel + 3, 'CB+1')):
            dif = np.array([(B[(k, b)] - vb_b) - (P[(k, b)] - vb_p) for k in range(204)])
            out.append('%s 均方根 %.0f、最大 %.0f meV' % (nm, np.sqrt(np.mean(dif ** 2)) * 1e3, np.abs(dif).max() * 1e3))
        kvb_p = max(range(204), key=lambda k: P[(k, nel)]); kvb_b = max(range(204), key=lambda k: B[(k, nel)])
        print('  %s：带隙（路径上）直接 %.0f、插值 %.0f meV；价带顶 k 序号 直接 %d、插值 %d；%s'
              % (m, (cb_p - vb_p) * 1e3, (cb_b - vb_b) * 1e3, kvb_p, kvb_b, '；'.join(out)))
        if m.startswith('Sn'):
            c5 = json.load(open(os.path.join(bench_dir(m), 'c5_star', 'c5_masses.json')))
            vs = [(k, v) for k, v in c5['valleys'].items() if v['kind'] == 'v']
            for j, (k, v) in enumerate(vs):
                e_dir = v['E_extremum'] - c5['path_vbm']; e_int = B[(204 + j, nel)] - vb_b
                print('     c5 %s（镜面内真实价带谷）：相对路径价带顶 直接 %+.0f、插值 %+.0f meV' % (k, e_dir * 1e3, e_int * 1e3))


if __name__ == '__main__':
    main()
