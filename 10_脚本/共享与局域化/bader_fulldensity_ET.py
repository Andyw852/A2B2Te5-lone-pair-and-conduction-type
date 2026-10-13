#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bader_fulldensity_ET.py —— 用【全密度（芯+价）】定义 Bader 盆，在价密度上积分，给出文献口径的 ET

Raty 2021 (Phys. Status Solidi RRL 15, 2000534) 原文：
    "The Bader charge is then obtained from the integration of the full charge
     density over an atomic basin Omega. By comparing the Bader charges with
     those of the neutral atoms, one can obtain the number of ET between atoms."
    "In those calculations, the full density (charge + valence) was analyzed to
     define the Bader volumes."

做法：
    CHGCAR_sum = AECCAR0(芯) + AECCAR2(价)  -> pybader 用它划分盆
    在价 CHGCAR 上对每个盆积分 -> N_basin
    q = ZVAL - N_basin                      -> 逐原子电荷转移（单位 e，正=失电子）
"""
import json, os, re, sys
import numpy as np
import pandas as pd

if not hasattr(pd.Series, 'append'):
    pd.Series.append = lambda self, other, **kw: pd.concat([self, other], **kw)

from pybader.interface import Bader, python_config
from pybader import __config__ as PYBADER_CONFIG

BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
SRC = os.path.join(BASE, 'band_and_elf', 'Bader_TET')
OUT = os.path.join(BASE, 'Raty2021_ES_ET_补算', '04_Bader全密度_ET', '结果')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']


def grid_info(path):
    lines = open(path, errors='replace').read().splitlines()
    n = sum(int(x) for x in lines[6].split())
    i = 8 + n
    while lines[i].strip() == '':
        i += 1
    return [int(x) for x in lines[i].split()[:3]]


def parse_zval(outcar):
    txt = open(outcar, errors='replace').read().splitlines()
    vals = {}
    for i, line in enumerate(txt):
        m = re.search(r'TITEL\s*=\s*PAW_\S+\s+([A-Z][a-z]?)', line)
        if m:
            for nxt in txt[i:i + 12]:
                z = re.search(r'ZVAL\s*=\s*([0-9.]+)', nxt)
                if z:
                    vals[m.group(1)] = float(z.group(1))
                    break
    return vals


def species_list(poscar):
    L = open(poscar).read().splitlines()
    sp = L[5].split(); cnt = [int(x) for x in L[6].split()]
    return [e for e, c in zip(sp, cnt) for _ in range(c)]


def run_one(mat):
    d = os.path.join(SRC, mat, 'ES')
    cfg = python_config(PYBADER_CONFIG, 'DEFAULT')
    cfg['threads'] = 8
    print('  reading', os.path.join(d, 'CHGCAR'), flush=True)
    bader = Bader.from_file(os.path.join(d, 'CHGCAR'), **cfg)
    ref = np.zeros(bader.density.shape, dtype=np.float64)
    for f in ('AECCAR0', 'AECCAR2'):
        r = Bader.from_file(os.path.join(d, f), **cfg)
        ref += r.density
    bader.reference = ref
    bader()
    q = np.asarray(bader.atoms_charge, dtype=float)
    vol = np.asarray(bader.atoms_volume, dtype=float)
    elems = species_list(os.path.join(d, 'POSCAR'))
    zval = parse_zval(os.path.join(d, 'OUTCAR'))
    ncell = float(q.sum())
    rows = []
    for i, el in enumerate(elems):
        rows.append({'index': i, 'element': el, 'zval': zval[el],
                     'N_basin_valence': round(float(q[i]), 5),
                     'q_ET': round(float(zval[el] - q[i]), 5),
                     'basin_volume_A3': round(float(vol[i]), 4)})
    return {'material': mat, 'grid': grid_info(os.path.join(d, 'CHGCAR')),
            'reference': 'AECCAR0+AECCAR2 (full core+valence density) for basins; valence CHGCAR for integration',
            'n_electrons_valence_sum': round(ncell, 5),
            'atoms': rows}


def main():
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for m in MATS:
        if len(sys.argv) > 1 and m not in sys.argv[1:]:
            continue
        print('=== ', m, flush=True)
        res[m] = run_one(m)
        json.dump(res[m], open(os.path.join(OUT, m + '_bader_fulldensity.json'), 'w'), ensure_ascii=False, indent=2)
        print('   N_e(sum) =', res[m]['n_electrons_valence_sum'], flush=True)
    json.dump(res, open(os.path.join(OUT, 'Bader_fulldensity_ET.json'), 'w'), ensure_ascii=False, indent=2)
    print('done')


if __name__ == '__main__':
    main()
