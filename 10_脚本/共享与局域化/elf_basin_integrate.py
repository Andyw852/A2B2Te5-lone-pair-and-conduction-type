#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
elf_basin_integrate.py —— ELF f-localization 盆积分（Raty 2021 图 1/图 2 口径）

文献口径：
    "the number of electrons involved in those localization regions (so-called
     f-localization domains) depends critically upon the threshold f used"
    "the numbers correspond to the fractional number of electrons obtained by
     integrating an f-localization basin defined by an isosurface ELF value
     equal to 0.85"  (Raty 2021, Fig.1/Fig.2 caption)

做法（与文献同等定义，纯后处理，不需要新 DFT）：
    1. 读 VASP ELFCAR（价电子 ELF）与 CHGCAR（价电荷密度）。
    2. 在 ELF 网格上取阈值 f，用**周期性** 26 邻居连通域标记得到 f-localization 盆。
    3. 把价电荷按体积权重映射到 ELF 网格，在每个盆内积分 -> 该盆的电子数。
    4. 按盆质心归属原子 / 判定是否位于层间 vdW 间隙。

单位：电子数 e、体积 A^3、分数坐标。
CHGCAR 的数值约定：写入的是 rho(r) * V_cell，所以
    电子数(细网格某体素) = data / N_grid_fine
"""
import json, os, sys
import numpy as np

BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
SRC = os.path.join(BASE, 'band_and_elf', 'Bader_TET')
OUT = os.path.join(BASE, 'Raty2021_ES_ET_补算', '03_ELF盆积分', '结果')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
THRESHOLDS = [0.80, 0.85, 0.90]


def parse_vasp_grid(path, want_grid=False):
    """返回 (lattice(3x3), species, frac positions, grid dims, data array)."""
    lines = open(path, errors='replace').read().splitlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    sp = lines[5].split()
    cnt = [int(x) for x in lines[6].split()]
    n = sum(cnt)
    sp_atoms = [e for e, c in zip(sp, cnt) for _ in range(c)]
    pos = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(8, 8 + n)])
    if lines[7].strip().lower().startswith(('c', 'k')):
        pos = pos @ np.linalg.inv(lat)
    i = 8 + n
    while lines[i].strip() == '':
        i += 1
    gd = [int(x) for x in lines[i].split()[:3]]
    i += 1
    vals = []
    need = gd[0] * gd[1] * gd[2]
    while len(vals) < need:
        vals.extend(float(x) for x in lines[i].split())
        i += 1
    data = np.array(vals[:need]).reshape((gd[2], gd[1], gd[0]))  # 索引 [z,y,x]，x 最快
    return lat, sp_atoms, pos, gd, data


def downsample_charge_elec(charge, elf_grid):
    """粗体素内电子数 = 细体素 data 之和 / 细网格总点数（data 已含 V_cell 因子）。"""
    nz, ny, nx = charge.shape
    gnx, gny, gnz = elf_grid          # VASP 里网格写的是 NGX NGY NGZ
    mz, my, mx = gnz, gny, gnx        # 数据数组索引顺序是 [z,y,x]
    fz, fy, fx = nz // mz, ny // my, nx // mx
    block = charge.reshape(mz, fz, my, fy, mx, fx).sum(axis=(1, 3, 5))
    return block / float(nz * ny * nx)


def label_periodic(mask):
    """周期性 26 邻居连通域标记，返回 (labels, nlab)。"""
    nx, ny, nz = mask.shape[2], mask.shape[1], mask.shape[0]
    lab = np.zeros(mask.shape, dtype=np.int32)
    nlab = 0
    offs = [(dz, dy, dx) for dz in (-1, 0, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dz, dy, dx) != (0, 0, 0)]
    for z0 in range(nz):
        for y0 in range(ny):
            for x0 in range(nx):
                if mask[z0, y0, x0] and lab[z0, y0, x0] == 0:
                    nlab += 1
                    stack = [(z0, y0, x0)]
                    lab[z0, y0, x0] = nlab
                    while stack:
                        z, y, x = stack.pop()
                        for dz, dy, dx in offs:
                            zz, yy, xx = (z + dz) % nz, (y + dy) % ny, (x + dx) % nx
                            if mask[zz, yy, xx] and lab[zz, yy, xx] == 0:
                                lab[zz, yy, xx] = nlab
                                stack.append((zz, yy, xx))
    return lab, nlab


def main():
    os.makedirs(OUT, exist_ok=True)
    summary = {}
    for m in MATS:
        d = os.path.join(SRC, m, 'ES')
        lat, sp_atoms, pos, gd_e, elf = parse_vasp_grid(os.path.join(d, 'ELFCAR'))
        _, _, _, gd_c, chg = parse_vasp_grid(os.path.join(d, 'CHGCAR'))
        elec = downsample_charge_elec(chg, gd_e)
        cellvol = abs(np.linalg.det(lat))
        voxvol = cellvol / elf.size
        # 粗网格坐标（分数）
        mz, my, mx = gd_e[2], gd_e[1], gd_e[0]
        zz, yy, xx = np.meshgrid(np.arange(mz), np.arange(my), np.arange(mx), indexing='ij')
        frac = np.stack([(xx + 0.5) / mx, (yy + 0.5) / my, (zz + 0.5) / mz], axis=-1)
        # 每个体素到最近原子的（周期性）分数距离（用笛卡尔距离，考虑近邻像）
        imgs = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
        dist = np.full(elf.shape + (len(sp_atoms),), 1e9)
        for a, p in enumerate(pos):
            for i, j, k in imgs:
                dv = frac - (p + np.array([i, j, k]))
                cart = dv @ lat
                dist[..., a] = np.minimum(dist[..., a], np.linalg.norm(cart, axis=-1))
        nearest = np.argmin(dist, axis=-1)
        dmin = np.min(dist, axis=-1)

        rec = {'material': m, 'elf_grid': [int(x) for x in gd_e], 'chg_grid': [int(x) for x in gd_c],
               'cell_volume_A3': round(float(cellvol), 3), 'voxel_volume_A3': round(float(voxvol), 6),
               'n_electrons_valence': round(float(elec.sum()), 4),
               'atoms': [{'i': i, 'el': sp_atoms[i], 'frac': [round(float(x), 6) for x in pos[i]]} for i in range(len(sp_atoms))],
               'thresholds': {}}
        for f in THRESHOLDS:
            mask = elf >= f
            lab, nlab = label_periodic(mask)
            doms = []
            for L in range(1, nlab + 1):
                sel = lab == L
                ne = float(elec[sel].sum())
                nvox = int(sel.sum())
                cfrac = frac[sel].mean(axis=0)
                # 质心归属原子（周期性）
                best, bestd = None, 1e9
                for a, p in enumerate(pos):
                    for i, j, k in imgs:
                        dv = cfrac - (p + np.array([i, j, k]))
                        dd = np.linalg.norm(dv @ lat)
                        if dd < bestd:
                            bestd, best = dd, a
                elems = {}
                for a in np.unique(nearest[sel]):
                    elems[sp_atoms[a]] = elems.get(sp_atoms[a], 0) + int((nearest[sel] == a).sum())
                doms.append({'id': L, 'n_voxels': nvox, 'volume_A3': round(nvox * voxvol, 4),
                             'electrons': round(ne, 5), 'centroid_frac': [round(float(x) % 1.0, 5) for x in cfrac],
                             'nearest_atom_index': int(best), 'nearest_atom_el': sp_atoms[best],
                             'centroid_atom_distance_A': round(float(bestd), 3),
                             'nearest_atom_element_histogram': elems})
            doms.sort(key=lambda r: -r['electrons'])
            rec['thresholds'][str(f)] = {'n_domains': nlab, 'total_electrons_in_domains': round(sum(r['electrons'] for r in doms), 5), 'domains': doms}
        summary[m] = rec
        print(m, 'done: f=0.85 domains =', rec['thresholds']['0.85']['n_domains'],
              'electrons in domains =', rec['thresholds']['0.85']['total_electrons_in_domains'])
    json.dump(summary, open(os.path.join(OUT, 'ELF_f_localization_domains.json'), 'w'), ensure_ascii=False, indent=2)
    print('written', os.path.join(OUT, 'ELF_f_localization_domains.json'))


if __name__ == '__main__':
    main()
