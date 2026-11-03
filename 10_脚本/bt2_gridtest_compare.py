"""PBE+SOC 网格收敛检验（分支 bt2-gridtest）的比较，用来选 HSE 的 k 网格（2026-11-03）。

1. GENE 路线本身：完整密网格（28×28×7/6，全网格约化到不可约点）经 GENE 跑的 BoltzTraP2，对 AutoZT 原有的密网格结果
   （同一 vasprun.xml、同一剪刀，直接读全网格），300、500 K。
2. 插值倍数：各网格 mult 5 对 20。
3. 子网格对完整密网格（都取 mult 20），300、600 K，p、n 各 3×10¹⁹–2×10²¹ cm⁻³，以及四个样品的霍尔浓度。

用法：python 10_脚本/bt2_gridtest_compare.py      只读，不写文件。
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from bt2_hse_validate import s_at, dense_json

GT = os.path.join(ROOT, '22_能带基准_1-0', 'bt2_gridtest')
FULL = {'Pb2Sb2Te5': '28x28x7', 'Sn2Sb2Te5': '28x28x7', 'Pb2Bi2Te5': '28x28x6', 'Sn2Bi2Te5': '28x28x6'}
SUBS = {'Pb2Sb2Te5': ['14x14x7'], 'Sn2Sb2Te5': ['14x14x7'], 'Pb2Bi2Te5': ['14x14x6', '28x28x3', '14x14x3'], 'Sn2Bi2Te5': ['14x14x6', '28x28x3', '14x14x3']}
NS = (3e19, 1e20, 3e20, 1e21, 2e21, -3e19, -1e20, -3e20, -1e21, -2e21)
HALL = {'Pb2Bi2Te5': -8.11e19, 'Sn2Bi2Te5': 9.07e20, 'Pb2Sb2Te5': 1.04e21, 'Sn2Sb2Te5': 1.72e21}


def gj(m, grid, mult=20):
    return os.path.join(GT, m, 'pbe_' + grid, 'bt2_mult%d' % mult, 'boltztrap_crta.json')


def rel(a, b):
    return (a - b) / abs(b)


def main():
    print('# 0 来源：GENE 的 vasprun.xml 指纹与 AutoZT 密网格结果的输入一致')
    for m in FULL:
        src = json.load(open(os.path.join(GT, m, 'pbe_' + FULL[m], 'source.json')))
        ref = json.load(open(dense_json(m)))['input_files']['step3_uniform/vasprun.xml']
        print('  %s：%s；%s' % (m, '一致' if src['source_sha256'] == ref else '不一致！', src['grid']))

    print('\n# 1 GENE 路线（完整密网格，不可约点，mult 20）对 AutoZT 原结果（全网格，mult 2），面内 S')
    for m in FULL:
        w = 0; cells = []
        for T in (300, 500):
            for N in NS:
                a, b = s_at(gj(m, FULL[m]), N, T), s_at(dense_json(m), N, T); w = max(w, abs(rel(a, b)))
                if T == 300 and abs(N) in (1e20, 1e21): cells.append('%+.0e %.1f/%.1f' % (N, a, b))
        print('  %s：最大相对差 %.1f%%（300、500 K，±3e19–2e21）；例 300 K：%s' % (m, w * 100, '  '.join(cells)))

    print('\n# 2 插值倍数 mult 5 对 20（300 K，±3e19–2e21 的最大相对差）')
    for m in FULL:
        print('  %s：' % m + '；'.join('%s %.1f%%' % (g, 100 * max(abs(rel(s_at(gj(m, g, 5), N, 300), s_at(gj(m, g, 20), N, 300))) for N in NS))
                                    for g in [FULL[m]] + SUBS[m]))

    print('\n# 3 子网格对完整密网格（mult 20），面内 S 的相对差；括号内为最大差出现处')
    for m in FULL:
        for g in SUBS[m]:
            out = []
            for T in (300, 600):
                d = [(abs(rel(s_at(gj(m, g), N, T), s_at(gj(m, FULL[m]), N, T))), N) for N in NS]
                w, Nw = max(d)
                out.append('%d K 最大 %.1f%%（%+.0e）、均值 %.1f%%' % (T, w * 100, Nw, 100 * np.mean([x for x, _ in d])))
            h = HALL[m]; hs = [(s_at(gj(m, g), h, T), s_at(gj(m, FULL[m]), h, T)) for T in (300, 600)]
            out.append('霍尔浓度 %+.2e：300 K %.1f/%.1f、600 K %.1f/%.1f' % (h, hs[0][0], hs[0][1], hs[1][0], hs[1][1]))
            print('  %s %s：%s' % (m, g, '；'.join(out)))


if __name__ == '__main__':
    main()
