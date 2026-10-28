#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orbital_fatband.py —— 一张图：四材料 HSE+SOC 能带，用【颜色+线宽】表示不同元素轨道的参与程度

数据：step4_fatband_plot/{PROCAR_OPT 或 PROCAR}（路径 k 点的 lm 投影）+ step4_band_plot/band.dat（能量）
      PROCAR 的 SOC 每带 4 个块：第 1 块是“电荷”投影（后 3 块是自旋分量，不用于轨道成分）。
      每个块的 9 条 ion 行 = 每个离子的 s/py/pz/px/dxy/dyz/dz2/dxz/x2-y2 权重。

画法（参考 α/β-PbO 的 fat-band 风格）：
  对 6 个“元素-轨道”类别（A-s, A-p, B-s, B-p, Te-s, Te-p）分别画同一条能带，
  线宽 ∝ 该类别在该 (k, band) 的权重 —— 于是“哪条带、在哪些 k、被哪个元素的哪个轨道占据”直接可见。
  颜色：Te 蓝、A 阳离子红、B 阳离子绿；s 用浅色、p 用深色。

输出：一个文件 orbital_fatband_4mats.png/pdf（2×2 面板 + 共享图例）
"""
import json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
SRC = os.path.join(BASE, 'band_and_elf', 'huang')
OUT = os.path.join(BASE, 'band_refine_20261021')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
EROLE = {'Pb': 'A', 'Sn': 'A', 'Sb': 'B', 'Bi': 'B', 'Te': 'Te'}
COL = {('A', 's'): '#e6a0a0', ('A', 'p'): '#c0392b',
       ('B', 's'): '#a3d9a5', ('B', 'p'): '#1e8449',
       ('Te', 's'): '#9ecae1', ('Te', 'p'): '#1f4e9c'}
LW_MAX = 5.2
LW_MIN = 0.15

plt.rcParams.update({'font.size': 10.5, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8})


def read_poscar_elements(path):
    L = open(path).read().splitlines()
    return L[5].split(), [int(x) for x in L[6].split()]


def parse_procar(path, elements, counts):
    """返回 kvec(nk,3), kweight(nk), E(nk,nb), W(nk,nb,ncat)，ncat = 3*n_el（每元素 s,p,d）。"""
    lines = open(path, errors='replace').read().splitlines()
    h = lines[1]
    nk = int(re.search(r'# of k-points:\s*(\d+)', h).group(1))
    nb = int(re.search(r'# of bands:\s*(\d+)', h).group(1))
    ni = int(re.search(r'# of ions:\s*(\d+)', h).group(1))
    ne = len(elements)
    assert ni == sum(counts), (ni, counts)
    ion_el = [e for e, c in zip(elements, counts) for _ in range(c)]
    ncat = 3 * ne
    K = np.zeros((nk, 3)); KW = np.zeros(nk)
    E = np.full((nk, nb), np.nan)
    W = np.zeros((nk, nb, ncat))
    ik = -1; ib = -1
    i = 3
    while i < len(lines):
        s = lines[i]
        if s.lstrip().startswith('k-point'):
            ik += 1
            ib = -1
            m = re.search(r'k-point\s+\d+\s*:\s*(.*?)\s+weight\s*=\s*([\d.Ee+-]+)', s)
            if m:
                K[ik] = [float(x) for x in re.findall(r'[-+]?\d+\.\d+', m.group(1))][:3]
                KW[ik] = float(m.group(2))
            i += 1; continue
        if s.lstrip().startswith('band '):
            ib += 1
            m = re.search(r'# energy\s*([-0-9.]+)', s)
            E[ik, ib] = float(m.group(1))
            i += 1
            while i < len(lines) and (lines[i].strip() == '' or lines[i].lstrip().startswith('ion ')):
                i += 1
            # 第一个块 = 电荷投影（9 条 ion 行）
            for j in range(ni):
                if i >= len(lines):
                    break
                row = lines[i].split()
                if len(row) < 10 or not row[0].isdigit():
                    break
                vals = [float(x) for x in row[1:10]]
                el = ion_el[j]
                e_idx = elements.index(el)
                W[ik, ib, 3 * e_idx + 0] += vals[0]                      # s
                W[ik, ib, 3 * e_idx + 1] += vals[1] + vals[2] + vals[3]  # p
                W[ik, ib, 3 * e_idx + 2] += sum(vals[4:9])               # d
                i += 1
            # 跳过本带剩余内容（tot 行 + 后 3 个自旋块），直到下一个 band/k-point
            while i < len(lines) and not (lines[i].lstrip().startswith('band ') or lines[i].lstrip().startswith('k-point')):
                i += 1
            continue
        i += 1
    if ik + 1 < nk:
        K, KW, E, W = K[:ik + 1], KW[:ik + 1], E[:ik + 1], W[:ik + 1]
    return K, KW, E, W, ion_el


def load_material(m):
    d = os.path.join(SRC, m)
    elements, counts = read_poscar_elements(os.path.join(d, 'step4_band_plot', 'POSCAR'))
    # 路径投影：优先 PROCAR_OPT（只有路径点），否则从 PROCAR 里挑零权重点
    p_opt = os.path.join(d, 'step4_fatband_plot', 'PROCAR_OPT')
    p_all = os.path.join(d, 'step4_fatband_plot', 'PROCAR')
    if os.path.exists(p_opt):
        K, KW, E, W, ion_el = parse_procar(p_opt, elements, counts)
    else:
        K, KW, E, W, ion_el = parse_procar(p_all, elements, counts)
        keep = np.where(np.abs(KW) <= 1e-12)[0]
        K, KW, E, W = K[keep], KW[keep], E[keep], W[keep]
    # 能量表（band.dat 已相对 VBM，且已验证与 PROCAR 一致）
    raw = np.loadtxt(os.path.join(d, 'step4_band_plot', 'band.dat'))
    kdist = raw[:, 0]; Edat = raw[:, 1:]
    if Edat.shape[0] != E.shape[0]:
        raise RuntimeError('%s: band.dat 行数 %d 与 PROCAR 路径点数 %d 不一致' % (m, Edat.shape[0], E.shape[0]))
    # 用 band.dat 的能量（相对 VBM），投影按同一 (k, band) 索引对齐
    # 校验：两者带宽一致
    de = np.abs((Edat - Edat.min()) - (E - E.min()))
    print('  %s: %d k × %d bands, 投影-能量对齐最大偏差 %.3f eV' % (m, E.shape[0], E.shape[1], de.max()))
    labels = []
    for line in open(os.path.join(d, 'step4_band_plot', 'band_klabels.txt')):
        p = line.split()
        if len(p) == 2 and not p[0].startswith('#'):
            labels.append((p[0], float(p[1])))
    sm = json.load(open(os.path.join(d, 'step4_band_plot', 'band_summary.json')))
    nv = int(sm['nelect'])
    # 类别：把 (元素, s/p/d) 合并成 6 类（d 归入 other）
    cats = []
    for el in elements:
        role = EROLE[el]
        cats.append((el, role, 's'))
        cats.append((el, role, 'p'))
    cat_idx = []
    for el, role, orb in cats:
        e = elements.index(el)
        cat_idx.append(3 * e + (0 if orb == 's' else 1))
    Wc = W[:, :, cat_idx]                      # (nk, nb, 6)
    tot = Wc.sum(axis=-1)
    with np.errstate(invalid='ignore', divide='ignore'):
        frac = np.where(tot[..., None] > 1e-6, Wc / tot[..., None], 0.0)
    return dict(m=m, kdist=kdist, E=Edat, Wc=Wc, frac=frac, tot=tot,
                labels=labels, sm=sm, nv=nv, elements=elements, cats=cats)


def panel(ax, r):
    m = r['m']; k = r['kdist']; E = r['E']; frac = r['frac']; tot = r['tot']
    sm = r['sm']; nv = r['nv']
    vbm = 0.0
    nk, nb = E.shape
    cats = r['cats']
    # 背景：所有带，细灰线
    for b in range(nb):
        ax.plot(k, E[:, b], color='0.82', lw=0.45, zorder=1)
    # fat bands：6 个类别各用 LineCollection 画（线段宽度随该点的权重变化）
    from matplotlib.collections import LineCollection
    order = [5, 4, 3, 2, 1, 0] if len(cats) == 6 else list(range(len(cats)))
    for zz, ci in enumerate(order):
        el, role, orb = cats[ci]
        colr = COL[(role, orb)]
        w = frac[:, :, ci]                      # (nk, nb)
        segs, lws = [], []
        for b in range(nb):
            e = E[:, b]
            if e.max() < -3.2 or e.min() > 3.2:
                continue
            wb = w[:, b]
            if wb.max() < 0.04:
                continue
            for j in range(nk - 1):
                segs.append([(k[j], e[j]), (k[j + 1], e[j + 1])])
                lws.append(LW_MIN + LW_MAX * 0.5 * (wb[j] + wb[j + 1]))
        if not segs:
            continue
        ax.add_collection(LineCollection(segs, colors=colr, linewidths=lws,
                                         capstyle='round', zorder=2 + zz))
    # 高对称点
    for name, x in r['labels']:
        ax.axvline(x, color='0.35', lw=0.7, zorder=0)
        ax.text(x, 1.005, name, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=9)
    ax.axhline(0, color='0.3', lw=0.6, ls='--', zorder=0)
    iv = int(np.argmax(E[:, :nv][:, -1] if False else E[:, nv - 1]))
    ax.set_xlim(k[0], k[-1]); ax.set_ylim(-3.2, 3.2)
    ax.set_ylabel('E - E$_{VBM}$ (eV)', fontsize=9.5)
    ax.set_title('%s   %s gap %.3f eV   (A=%s, B=%s)' % (
        m, sm['gap_type'], sm['gap_eV'],
        [e for e in r['elements'] if EROLE[e] == 'A'][0],
        [e for e in r['elements'] if EROLE[e] == 'B'][0]), loc='left', fontsize=10.5)


def main():
    rs = [load_material(m) for m in MATS]
    fig, axs = plt.subplots(2, 2, figsize=(13.2, 10.2), layout='constrained')
    for ax, r in zip(axs.flat, rs):
        panel(ax, r)
    # 共享图例
    handles = []
    for role, orb, lab in [('Te', 'p', 'Te-p'), ('Te', 's', 'Te-s'),
                           ('A', 'p', 'A-p (Pb/Sn)'), ('A', 's', 'A-s'),
                           ('B', 'p', 'B-p (Sb/Bi)'), ('B', 's', 'B-s')]:
        handles.append(Line2D([], [], color=COL[(role, orb)], lw=5.0, label=lab))
    fig.legend(handles=handles, loc='lower center', ncol=6, frameon=False, fontsize=10.5,
               title='line width ∝ orbital contribution of that (element, orbital) category',
               title_fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle('HSE+SOC bands of A$_2$B$_2$Te$_5$: color and width = element-orbital participation',
                 fontsize=13.5, y=1.0)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, 'orbital_fatband_4mats.' + ext), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('written', os.path.join(OUT, 'orbital_fatband_4mats.png'))


if __name__ == '__main__':
    main()
