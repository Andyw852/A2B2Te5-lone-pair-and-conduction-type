#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export_plot_source_data.py —— 打包“能带绘图源数据”，用于上传 GitHub
（供 band_refine.py / orbital_fatband.py 复现；不复制 60–90 MB 的原始 PROCAR，
  只放压缩后的解析数组 + 原始文件 SHA256 追溯）

输出结构（默认 staging/）：
  源数据/<材料>/{band.dat, band_klabels.txt, band_summary.json, POSCAR,
                INCAR_band, INCAR_fatband, KPOINTS, KPOINTS_OPT,
                投影_解析.npz, 投影_解析.json, 原始来源.json}
  摘要/{band_summary_refined.csv,md,json}
  校验/SHA256.json
"""
import hashlib, json, os, shutil
import numpy as np

BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
SRC = os.path.join(BASE, 'band_and_elf', 'huang')
PROJ = os.path.join(BASE, 'band_refine_20261021')
STAGE = os.path.join(PROJ, 'staging')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def main():
    if os.path.isdir(STAGE):
        shutil.rmtree(STAGE)
    os.makedirs(os.path.join(STAGE, '源数据'))
    os.makedirs(os.path.join(STAGE, '摘要'))
    os.makedirs(os.path.join(STAGE, '校验'))
    manifest = {}

    # 复用同目录 orbital_fatband.py 的解析函数（同一套口径）
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'of', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'orbital_fatband.py'))
    of = importlib.util.module_from_spec(spec); spec.loader.exec_module(of)

    for m in MATS:
        d = os.path.join(SRC, m)
        out = os.path.join(STAGE, '源数据', m)
        os.makedirs(out)
        bdir = os.path.join(d, 'step4_band_plot')
        fdir = os.path.join(d, 'step4_fatband_plot')
        # 1) 原样小文件
        for src, dst in [
            (os.path.join(bdir, 'band.dat'), 'band.dat'),
            (os.path.join(bdir, 'band_klabels.txt'), 'band_klabels.txt'),
            (os.path.join(bdir, 'band_summary.json'), 'band_summary.json'),
            (os.path.join(bdir, 'POSCAR'), 'POSCAR'),
            (os.path.join(fdir, 'INCAR'), 'INCAR_fatband'),
            (os.path.join(bdir, 'KPOINTS'), 'KPOINTS'),
        ]:
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(out, dst))
        for fn in ('KPOINTS_OPT',):
            p = os.path.join(bdir, fn)
            if os.path.exists(p):
                shutil.copy2(p, os.path.join(out, fn))
        # band 计算 INCAR（在 step4_HSE_band）
        inc_b = os.path.join(d, 'step4_HSE_band', 'INCAR')
        if os.path.exists(inc_b):
            shutil.copy2(inc_b, os.path.join(out, 'INCAR_band'))
        # 2) 解析投影 + 能量（与绘图所用完全一致）
        r = of.load_material(m, SRC)              # 显式读原工作区的 PROCAR
        np.savez_compressed(
            os.path.join(out, '投影_解析.npz'),
            kdist=r['kdist'].astype(np.float32),
            E=r['E'].astype(np.float32),
            frac=r['frac'].astype(np.float32),
            tot=r['tot'].astype(np.float32),
            nv=np.array([r['nv']]),
            cats=np.array(['%s-%s' % (c[0], c[2]) for c in r['cats']]),
        )
        # 3) 原始大文件（PROCAR/PROCAR_OPT）的来源与校验
        raw = []
        for fn in ('PROCAR_OPT', 'PROCAR'):
            p = os.path.join(fdir, fn)
            if os.path.exists(p):
                raw.append(dict(file='step4_fatband_plot/' + fn, size_bytes=os.path.getsize(p),
                                sha256=sha256(p)))
        prov = dict(
            material=m,
            note='投影取 PROCAR 每带第 1 个块（电荷投影）；后 3 个块是自旋分量，不参与轨道成分。'
                 '9 条 ion 行的 s/py/pz/px/d 权重按元素归并为 s/p。',
            raw_files=raw,
            parsed=dict(kdist='路径累积 k 距离 (1/Ang)', E='E - E_VBM (eV)，与 band.dat 相同',
                        frac='(nk, nb, ncat) 各 (元素-轨道) 类别的相对份额（按该带投影总和归一）',
                        tot='(nk, nb) 该带投影总权重（原子球内份额，约 0.44–0.58）',
                        cats='类别顺序 = 元素顺序 × (s,p)；元素顺序见 POSCAR'),
        )
        open(os.path.join(out, '原始来源.json'), 'w').write(json.dumps(prov, ensure_ascii=False, indent=2))
        manifest[m] = dict(files=sorted(os.listdir(out)), raw=raw)
        print(m, '->', sorted(os.listdir(out)))

    # 摘要（绘图脚本生成的数值表）
    for fn in ('band_summary_refined.csv', 'band_summary_refined.md', 'band_summary_refined.json'):
        p = os.path.join(PROJ, fn)
        if os.path.exists(p):
            shutil.copy2(p, os.path.join(STAGE, '摘要', fn))
    # 校验清单
    rows = {}
    for root, _, files in os.walk(STAGE):
        for fn in files:
            p = os.path.join(root, fn)
            rel = os.path.relpath(p, STAGE)
            if rel.startswith('校验'):
                continue
            rows[rel] = dict(size=os.path.getsize(p), sha256=sha256(p))
    json.dump(rows, open(os.path.join(STAGE, '校验', 'SHA256.json'), 'w'), ensure_ascii=False, indent=2)
    open(os.path.join(STAGE, '校验', '原始大文件_SHA256.json'), 'w').write(
        json.dumps({m: manifest[m]['raw'] for m in manifest}, ensure_ascii=False, indent=2))
    total = sum(v['size'] for v in rows.values())
    print('staged %d files, %.1f MB -> %s' % (len(rows), total / 1e6, STAGE))


if __name__ == '__main__':
    main()
