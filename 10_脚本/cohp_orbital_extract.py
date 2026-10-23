# -*- coding: utf-8 -*-
"""从 COGITO 的 all_bonds.json 提取轨道分辨 ICOHP（按实际键长分组）。

用法: python cohp_orbital_extract.py COHP_ROOT OUT.csv
  COHP_ROOT 下需有 <材料>_lam<1|0>/cohp-cogito/step1_cohp/all_bonds.json

口径:
  - COHP 约定：负 = 成键，正 = 反键。
  - all_bonds.json 中每根键按 (X,Te) 与 (Te,X) 两个方向各存一条，每条给该方向的一半 ICOHP；
    本脚本把两个方向相加，得到该键的 ICOHP（eV/胞口径同 bond_info.txt）。
  - 按实际键长分组（四舍五入到 0.01 A），因此 B-Te 的短键/长键会各自成组，而不是用 3.2 A 一刀切。
"""
import json, os, csv, sys
from collections import defaultdict
MATS = ['Pb2Sb2Te5','Sn2Sb2Te5','Pb2Bi2Te5','Sn2Bi2Te5']
A_OF = {'Pb2Sb2Te5':'Pb','Sn2Sb2Te5':'Sn','Pb2Bi2Te5':'Pb','Sn2Bi2Te5':'Sn'}
B_OF = {'Pb2Sb2Te5':'Sb','Sn2Sb2Te5':'Sb','Pb2Bi2Te5':'Bi','Sn2Bi2Te5':'Bi'}
ORBS = ['s-s','s-p','p-s','p-p']
def norm(e):
    return str(e).replace('_d','').replace('_sv','').replace('_pv','').strip()
def main(root, out):
    rows = []
    for m in MATS:
        A, B = A_OF[m], B_OF[m]
        for lam in ('1','0'):
            p = os.path.join(root, '%s_lam%s' % (m, lam), 'cohp-cogito/step1_cohp/all_bonds.json')
            if not os.path.exists(p):
                continue
            d = json.load(open(p))
            agg = defaultdict(float)   # (pair, dkey, orb) -> ICOHP
            ncnt = defaultdict(int)
            dsum = defaultdict(list)
            for k, b in d.items():
                if not isinstance(b, dict):
                    continue
                els = [norm(x) for x in b.get('elements', [])]
                if len(els) != 2:
                    continue
                bl = b.get('bond_length', 99)
                if bl > 4.0:
                    continue
                dkey = round(bl, 2)
                for orb, v in (b.get('orb-orb dict') or {}).items():
                    cm = v.get('cohp_matrix')
                    if not cm:
                        continue
                    s = sum(x for a in cm for bb in a for x in bb)
                    agg[(tuple(els), dkey, orb)] += s
                    ncnt[(tuple(els), dkey, orb)] += 1
                dsum[(tuple(els), dkey)].append(bl)
            for pair, lab in (((A,'Te'),'A-Te'), ((B,'Te'),'B-Te')):
                keys = sorted(set((dkey) for (p2, dkey) in dsum if p2 == tuple(pair)))
                for dkey in keys:
                    vals = {o: agg.get((tuple(pair), dkey, o), 0.0) for o in ORBS}
                    tot = sum(vals.values())
                    rows.append({'material':m, 'lambda':lam, 'bond':lab,
                                 'bond_length_A': dkey,
                                 's-s': round(vals['s-s'],4),
                                 'Bs-Tep' if lab.startswith('B') else 'As-Tep': round(vals['s-p'],4),
                                 'Bp-Tes' if lab.startswith('B') else 'Ap-Tes': round(vals['p-s'],4),
                                 'p-p': round(vals['p-p'],4),
                                 'total': round(tot,4),
                                 'n_directed': ncnt.get((tuple(pair), dkey, 's-s'), 0)})
    cols = ['material','lambda','bond','bond_length_A','s-s','Bs-Tep','As-Tep','Bp-Tes','Ap-Tes','p-p','total','n_directed']
    with open(out,'w',newline='',encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    print('wrote', out, len(rows), 'rows')
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])