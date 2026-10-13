#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
汇总DI与ES.py —— 把【折叠 DI】整理成 Raty 2021 表 2 的形式（ES = 2×DI，按接触重数列出）

为什么用“折叠 DI / 重数”而不是逐周期像分解：
    周期体系里每个原子的“盆”是周期集合；把积分限制在胞内时，任意周期像的盆在该胞内的部分
    都是同一块区域（周期性的直接后果）。因此逐周期像分解在 Bloch 表象下是退化的，
    要真正分出“哪一个像”必须先把波函数局域化（MLWF，文献用的正是这条路）。
    但【折叠 DI】(A,B) = Σ_R DI(A, B+R) 本身是良定义的、且满足求和规则 Σ_B DI(A,B) = 2N_A。
    当结构对称性使最近邻壳层里每个原子对恰好有 m 个等价像时，
        单接触 ES = 2 × 折叠 DI(A,B) / m，
    这正是文献表 2 “1.36 (3×) intra” 的读法（我们用它反推过 PbO/Te 的求和规则，数量级自洽）。
"""
import json, os, sys
import numpy as np

BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
PROJ = os.path.join(BASE, 'Raty2021_ES_ET_补算')
OUT = os.path.join(PROJ, '05_DI_ES', '结果')
SOC_ROOT = os.path.join(BASE, 'band_and_elf', 'huang')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
from pymatgen.core import Structure

RS = 3
TOL = 0.35          # 同一壳层的距离容差 (A)


def shells(structure, a, b, cut=6.0):
    """返回原子 b 相对原子 a 的各周期像距离列表。"""
    L = structure.lattice.matrix
    fr = np.array(structure.frac_coords)
    out = []
    for i in range(-RS, RS + 1):
        for j in range(-RS, RS + 1):
            for k in range(-RS, RS + 1):
                d = np.linalg.norm((fr[b] + np.array([i, j, k]) - fr[a]) @ L)
                if 0.1 < d < cut:
                    out.append(((i, j, k), d))
    return sorted(out, key=lambda x: x[1])


def classify(structure, a, b, R):
    """按“是否跨过 vdW 间隙平面”判层内/层间。

    间隙在 z≈0.5（Te0.428 与 Te0.572 之间，周期等价于 z≈-0.5），
    因此只要两原子之间的 z 区间跨过某个 (n+0.5) 平面，就算层间接触。
    """
    fr = np.array(structure.frac_coords)
    za = fr[a][2]
    dz = fr[b][2] + R[2] - za
    crossed = abs(np.floor(za + dz - 0.5) - np.floor(za - 0.5))
    return '层间' if crossed >= 1 else '层内'


def main():
    res = json.load(open(os.path.join(OUT, 'DI_results_lowdin.json')))
    md = ['# ES = 2×DI（Raty 2021 表 2 口径）—— 逐接触结果', '',
          'DI 用 **KS 波函数（SOC，ISYM=0 的 38 个 IBZ 点 + 空间群轨道展开）+ 全密度 Bader 盆** 计算：',
          '',
          '    DI(A,B) = 2·g_spin·Σ_k w_k (1/m_k) Σ_g Tr[ N(k, π_g⁻¹A) · N(k, π_g⁻¹B)† ]',
          '',
          '其中 N 是盆内自旋求和的一阶密度矩阵块，占据轨道先做 **Löwdin 正交化**',
          '（VASP 存的是 PAW 赝波函数，Σ|C|² 在 0.35–1.30 之间，本身不正交；不正交化时求和规则会偏 21%）。',
          '上面的 DI 是【折叠值】= Σ_R DI(A, B+R)；单接触 ES = 2×折叠 DI / 重数 m（m = 最近邻壳层里该原子对的等价周期像数），',
          '这正是文献 “1.36 (3×) intra” 的读法。','']
    rows_all = {}
    for m in MATS:
        if m not in res:
            continue
        structure = Structure.from_file(os.path.join(SOC_ROOT, m, 'step3_PBE_WAVECAR', 'POSCAR'))
        DI = np.array(res[m]['DI'])
        syms = [str(s) for s in structure.species]
        nat = len(syms)
        chk = res[m]['checks']
        rows = []
        for a in range(nat):
            for b in range(nat):
                if a == b:
                    continue
                sh = shells(structure, a, b)
                if not sh:
                    continue
                d0 = sh[0][1]
                near = [x for x in sh if x[1] < d0 + TOL]
                far = [x for x in sh if x[1] >= d0 + TOL]
                # 只保留“盆对盆”有意义的：折叠 DI 里近壳层贡献为主
                if DI[a, b] < 5e-3:
                    continue
                lab = classify(structure, a, b, near[0][0])
                rows.append(dict(a=a, b=b, pair='-'.join(sorted([syms[a], syms[b]])),
                                 layer=lab, d=round(d0, 3), m=len(near),
                                 DI_folded=float(DI[a, b]), ES_contact=float(2 * DI[a, b] / len(near)),
                                 far_shells=[(list(r), round(dd, 2)) for r, dd in far[:3]],
                                 self_DI=float(DI[a, a])))
        # 汇总到“接触类”：元素对 + 归属 + 距离
        classes = {}
        for r in rows:
            key = (r['pair'], r['layer'], round(r['d'], 2))
            classes.setdefault(key, []).append(r)
        table = []
        for (pair, lab, d), v in sorted(classes.items(), key=lambda kv: (kv[0][1] != '层间', kv[0][2])):
            table.append(dict(pair=pair, layer=lab, d=d, m=int(np.mean([x['m'] for x in v])),
                              n_pairs=len(v),
                              ES=float(np.mean([x['ES_contact'] for x in v])),
                              ES_std=float(np.std([x['ES_contact'] for x in v])),
                              DI_folded=float(np.mean([x['DI_folded'] for x in v]))))
        rows_all[m] = dict(classes=table, pairs=rows, checks=chk)
        md += ['## %s' % m, '',
               '求和规则校验：Σ_B DI(A,B) 与 2N_A 的最大相对误差 **%.3f%%**；DI 反对称残差 %.1e；k 轨道自检 %s（期望 %s）' % (
                   100 * chk['sum_rule_max_rel_err'], chk['asymmetry'],
                   chk['orbit_sizes'][:6], chk['orbit_expect'][:6]), '',
               '| 接触 | 归属 | d / A | 重数 m | 单接触 ES = 2×DI/m / e | 折叠 2×DI(A,B) / e | 类内标准差 |',
               '|---|---|---:|---:|---:|---:|---:|']
        for t in table:
            md.append('| %s | %s | %.2f | %d | **%.3f** | %.3f | %.3f |' % (
                t['pair'], t['layer'], t['d'], t['m'], t['ES'], 2 * t['DI_folded'], t['ES_std']))
        intra = [t for t in table if t['layer'] == '层内']
        inter = [t for t in table if t['layer'] == '层间']
        md += ['', '- **层内**：' + '；'.join('%s %.2f A  %.2f e (%d×)' % (t['pair'], t['d'], t['ES'], t['m']) for t in intra),
               '- **层间**：' + '；'.join('%s %.2f A  %.2f e (%d×)' % (t['pair'], t['d'], t['ES'], t['m']) for t in inter), '']
    md += ['## 与文献对照（Raty 2021 表 2 / 正文）', '',
           '| 体系 | 层内 ES / e | 层间 ES / e | TET / e | Gap / A |', '|---|---|---|---|---|',
           '| Te（文献） | 1.81 (2×) | 0.56 (4×) | – | 3.51 |',
           '| Sb（文献） | 1.45 (3×) | 0.52 (3×) | – | 3.45 |',
           '| a-PbO（文献） | 0.95 (4×) | 0.03 (4×) | 1.16 | 4.09 |',
           '| GST225（文献） | 1.36 (3×) Ge–Te | 0.36 (3×) Te–Te | 0.3 | 3.65 |',
           '| 本工作 | 见上表 | 见上表 | 见 04_Bader全密度_ET/结果 | 3.58–3.61 |', '',
           '注：文献用 DGrid 或 MLWF+CRITIC2 得到逐接触 DI；我们用 Bloch 波函数算折叠 DI 后按对称重数拆分，',
           '两者口径相同（HF 型 pair density、Bader 盆、ES = 2×DI），但数值路径不同；',
           '另外 VASP 存的是 PAW 赝波函数，本结果属“赝波函数口径”，会带约 4% 的原子球内增补项缺失。', '']
    open(os.path.join(OUT, 'ES_逐接触汇总.md'), 'w').write('\n'.join(md))
    json.dump(rows_all, open(os.path.join(OUT, 'ES_逐接触汇总.json'), 'w'), ensure_ascii=False, indent=2)
    print('\n'.join(md))


if __name__ == '__main__':
    main()
