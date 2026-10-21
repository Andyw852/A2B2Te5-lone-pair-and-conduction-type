#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orbital_fatband.py —— 一张图：四材料 HSE+SOC 能带，用【颜色+线宽】表示不同元素轨道的参与程度

数据：step4_fatband_plot/{PROCAR_OPT 或 PROCAR}（路径 k 点的 lm 投影）+ step4_band_plot/band.dat（能量）
      PROCAR 的 SOC 每带 4 个块：第 1 块是“电荷”投影（后 3 块是自旋分量，不用于轨道成分）。
      每个块的 9 条 ion 行 = 每个离子的 s/py/pz/px/dxy/dyz/dz2/dxz/x2-y2 权重。

画法（分层 fat-band，“线缆”式）：
  每条能带画成总线宽固定的一根“线缆”，由内到外依次是 A-s, B-s, Te-s, A-p, B-p, Te-p 六层，
  每层的厚度 ∝ 该“元素-轨道”在该 (k, band) 的份额 —— 于是“哪条带、在哪些 k、被哪个元素的哪个轨道占据”
  直接可见，且小份额（如 Sn-s 孤对）不会被大份额的层盖住（逐层累积线宽，外层先画、内层后画）。
  颜色：Te 蓝、A 阳离子红、B 阳离子绿；s 用浅色、p 用深色。路径断点（A|L、M|H）处不连线。

输出：一个文件 orbital_fatband_4mats.png/pdf（2×2 面板 + 共享图例）

数据位置（自动识别两种目录结构）：
  * clone 下来的仓库：默认读 12_能带绘图_20261021/源数据/<材料>/（band.dat + 投影_解析.npz），
    输出到 12_能带绘图_20261021/复现输出/；
  * 原工作区：<root>/<材料>/step4_band_plot + step4_fatband_plot（直接解析 PROCAR）。
  --procar-root 07_带边能窗/源数据 可改为从原始 PROCAR(.gz) 重新解析，并与 npz 逐元素比对。

用法：
  python 10_脚本/orbital_fatband.py
  python 10_脚本/orbital_fatband.py --materials Sn2Bi2Te5          # 只画单材料
  python 10_脚本/orbital_fatband.py --procar-root 07_带边能窗/源数据   # 从原始 PROCAR 复算
  python 10_脚本/orbital_fatband.py --data-root <root> --out <dir>
"""
import argparse, gzip, json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_PKG = os.path.normpath(os.path.join(HERE, '..', '12_能带绘图_20261021'))
BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
if os.path.isdir(os.path.join(REPO_PKG, '源数据')):          # 仓库内：读打包好的源数据
    SRC = os.path.join(REPO_PKG, '源数据')
    OUT = os.path.join(REPO_PKG, '复现输出')
else:                                                          # 原工作区
    SRC = os.path.join(BASE, 'band_and_elf', 'huang')
    OUT = os.path.join(BASE, 'band_refine_20261021')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
EROLE = {'Pb': 'A', 'Sn': 'A', 'Sb': 'B', 'Bi': 'B', 'Te': 'Te'}
COL = {('A', 's'): '#e6a0a0', ('A', 'p'): '#c0392b',
       ('B', 's'): '#a3d9a5', ('B', 'p'): '#1e8449',
       ('Te', 's'): '#9ecae1', ('Te', 'p'): '#1f4e9c'}
LAYERS = [('A', 's'), ('B', 's'), ('Te', 's'), ('A', 'p'), ('B', 'p'), ('Te', 'p')]   # 由内到外
LW_TOT = 4.0           # 每条带的总线宽 (pt)
W_MIN = 0.01           # 份额低于此值的层不画
EWIN = (-3.2, 3.2)     # 能窗 (eV, 相对 VBM)

plt.rcParams.update({'font.size': 10.5, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8})


def read_poscar_elements(path):
    L = open(path).read().splitlines()
    return L[5].split(), [int(x) for x in L[6].split()]


def mat_dirs(root, m):
    """(能带目录, 投影目录)。原工作区：<root>/<材料>/step4_band_plot | step4_fatband_plot；
    仓库打包（12_能带绘图_20261021/源数据、07_带边能窗/源数据）：文件直接在 <root>/<材料>/ 下。"""
    d = os.path.join(root, m)
    if os.path.isdir(os.path.join(d, 'step4_band_plot')):
        return os.path.join(d, 'step4_band_plot'), os.path.join(d, 'step4_fatband_plot')
    return d, d


def find_procar(d):
    """路径投影：优先 PROCAR_OPT（只有路径点），否则 PROCAR（需挑零权重点）；支持 .gz。"""
    for fn in ('PROCAR_OPT', 'PROCAR_OPT.gz', 'PROCAR', 'PROCAR.gz'):
        p = os.path.join(d, fn)
        if os.path.exists(p):
            return p, fn.startswith('PROCAR_OPT')
    return None, False


def parse_procar(path, elements, counts):
    """返回 kvec(nk,3), kweight(nk), E(nk,nb), W(nk,nb,ncat)，ncat = 3*n_el（每元素 s,p,d）。"""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', errors='replace') as f:
        lines = f.read().splitlines()
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


def load_material(m, data_root=None, procar_root=None):
    bdir, fdir = mat_dirs(data_root or SRC, m)
    elements, counts = read_poscar_elements(os.path.join(bdir, 'POSCAR'))
    # 能量表（band.dat 已相对 VBM，且已验证与 PROCAR 一致）
    raw = np.loadtxt(os.path.join(bdir, 'band.dat'))
    kdist = raw[:, 0]; Edat = raw[:, 1:]
    # 类别：把 (元素, s/p/d) 合并成 6 类（d 归入 other）
    cats = []
    for el in elements:
        role = EROLE[el]
        cats.append((el, role, 's'))
        cats.append((el, role, 'p'))
    npz = os.path.join(fdir, '投影_解析.npz')
    p_procar, is_opt = find_procar(mat_dirs(procar_root, m)[1] if procar_root else fdir)
    if procar_root and p_procar is None:
        raise FileNotFoundError('%s: %s 下没有 PROCAR_OPT/PROCAR(.gz)' % (m, mat_dirs(procar_root, m)[1]))
    if p_procar is not None:
        K, KW, E, W, ion_el = parse_procar(p_procar, elements, counts)
        if not is_opt:
            keep = np.where(np.abs(KW) <= 1e-12)[0]
            K, KW, E, W = K[keep], KW[keep], E[keep], W[keep]
        if Edat.shape[0] != E.shape[0]:
            raise RuntimeError('%s: band.dat 行数 %d 与 PROCAR 路径点数 %d 不一致' % (m, Edat.shape[0], E.shape[0]))
        # 用 band.dat 的能量（相对 VBM），投影按同一 (k, band) 索引对齐
        # 校验：扣掉 E_VBM 后能窗内逐点一致（能窗外最高几条空带在投影计算里未收敛，偏差大但不画）
        de = np.abs(E - np.median(E - Edat) - Edat)
        inw = (Edat >= EWIN[0]) & (Edat <= EWIN[1])
        print('  %s: %d k × %d bands, 投影-能量对齐偏差 能窗内 %.4f eV（全部 %.3f eV）  [%s]' % (
            m, E.shape[0], E.shape[1], de[inw].max(), de.max(), os.path.relpath(p_procar)))
        if de[inw].max() > 0.05:
            raise RuntimeError('%s: 能窗内 PROCAR 与 band.dat 能量对不上（%.3f eV），投影可能错位' % (m, de[inw].max()))
        cat_idx = []
        for el, role, orb in cats:
            e = elements.index(el)
            cat_idx.append(3 * e + (0 if orb == 's' else 1))
        Wc = W[:, :, cat_idx]                      # (nk, nb, 6)
        tot = Wc.sum(axis=-1)
        with np.errstate(invalid='ignore', divide='ignore'):
            frac = np.where(tot[..., None] > 1e-6, Wc / tot[..., None], 0.0)
        if os.path.exists(npz):                    # 与仓库里打包的解析结果比对
            z = np.load(npz)
            print('     与 %s 比对：frac 最大差 %.2e，tot 最大差 %.2e' % (
                os.path.relpath(npz), np.abs(frac - z['frac']).max(), np.abs(tot - z['tot']).max()))
    elif os.path.exists(npz):
        # 仓库打包的解析数组（与原工作区解析 PROCAR 的结果相同，float32）
        z = np.load(npz)
        if [str(c) for c in z['cats']] != ['%s-%s' % (c[0], c[2]) for c in cats]:
            raise RuntimeError('%s: npz 类别 %s 与 POSCAR 元素顺序不符' % (m, list(z['cats'])))
        frac = z['frac'].astype(float); tot = z['tot'].astype(float)
        Wc = frac * tot[..., None]
        if z['E'].shape != Edat.shape or np.abs(z['E'] - Edat).max() > 1e-3:
            raise RuntimeError('%s: npz 能量与 band.dat 不一致' % m)
        print('  %s: %d k × %d bands  [%s]' % (m, Edat.shape[0], Edat.shape[1], os.path.relpath(npz)))
    else:
        raise FileNotFoundError('%s: %s 下既没有 投影_解析.npz 也没有 PROCAR' % (m, fdir))
    labels = []
    for line in open(os.path.join(bdir, 'band_klabels.txt')):
        p = line.split()
        if len(p) == 2 and not p[0].startswith('#'):
            labels.append((p[0], float(p[1])))
    sm = json.load(open(os.path.join(bdir, 'band_summary.json')))
    nv = int(sm['nelect'])
    return dict(m=m, kdist=kdist, E=Edat, Wc=Wc, frac=frac, tot=tot,
                labels=labels, sm=sm, nv=nv, elements=elements, cats=cats)


def panel(ax, r):
    m = r['m']; k = r['kdist']; E = r['E']; frac = r['frac']
    sm = r['sm']
    nk, nb = E.shape
    # 线段：路径断点（A|L、M|H，相邻两点 k 距离相同）处不连线；只取与能窗相交的带
    jj = np.where(np.diff(k) > 1e-6)[0]
    bands = [b for b in range(nb) if E[:, b].max() >= EWIN[0] and E[:, b].min() <= EWIN[1]]
    Eb = E[:, bands]
    kb = np.broadcast_to(k[:, None], Eb.shape)
    seg = np.stack([np.stack([kb[jj], Eb[jj]], -1), np.stack([kb[jj + 1], Eb[jj + 1]], -1)], axis=2).reshape(-1, 2, 2)
    w = 0.5 * (frac[jj][:, bands] + frac[jj + 1][:, bands]).reshape(len(seg), -1)   # (nseg, ncat)
    # 背景：细灰线（只在投影缺失时露出来）
    ax.add_collection(LineCollection(seg, colors='0.8', linewidths=0.4, zorder=1))
    # 分层：第 i 层线宽 = LW_TOT × (第 0..i 层份额之和)，外层先画、内层后画
    ci = [[c[1:] for c in r['cats']].index(lay) for lay in LAYERS]
    cum = np.cumsum(w[:, ci], axis=1)
    for i in range(len(LAYERS) - 1, -1, -1):
        keep = w[:, ci[i]] >= W_MIN
        if keep.any():
            ax.add_collection(LineCollection(seg[keep], colors=COL[LAYERS[i]], linewidths=LW_TOT * cum[keep, i],
                                             capstyle='round', joinstyle='round', zorder=2 + len(LAYERS) - i))
    # 高对称点
    xs = [x for _, x in r['labels']]
    for x in xs:
        ax.axvline(x, color='0.55', lw=0.6, zorder=0)
    ax.set_xticks(xs, [n for n, _ in r['labels']])
    ax.tick_params(axis='x', length=0, labelsize=9)
    ax.axhline(0, color='0.3', lw=0.6, ls='--', zorder=0)
    ax.set_xlim(k[0], k[-1]); ax.set_ylim(*EWIN)
    ax.set_ylabel('E - E$_{VBM}$ (eV)', fontsize=9.5)
    ax.set_title('%s   %s gap %.3f eV   (A=%s, B=%s)' % (
        re.sub(r'(\d)', r'$_\1$', m), sm['gap_type'], sm['gap_eV'],
        [e for e in r['elements'] if EROLE[e] == 'A'][0],
        [e for e in r['elements'] if EROLE[e] == 'B'][0]), loc='left', fontsize=10.5)


def main(argv=None):
    ap = argparse.ArgumentParser(description='四材料 HSE+SOC 轨道分辨 fat-band（颜色=元素-轨道，线宽∝参与度）')
    ap.add_argument('--data-root', default=SRC,
                    help='含 <材料>/band.dat… 的目录（默认 %(default)s）')
    ap.add_argument('--procar-root', default=None,
                    help='从原始 PROCAR(.gz) 重新解析投影，如 07_带边能窗/源数据；默认用 投影_解析.npz')
    ap.add_argument('--out', default=OUT, help='输出目录（默认 %(default)s）')
    ap.add_argument('--materials', nargs='+', default=MATS, choices=MATS, metavar='MAT',
                    help='要画的材料（默认四个：%s）' % ' '.join(MATS))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    mats = [m for m in MATS if m in a.materials]
    rs = [load_material(m, a.data_root, a.procar_root) for m in mats]
    nrow, ncol = {1: (1, 1), 2: (1, 2)}.get(len(rs), (2, 2))
    fig, axs = plt.subplots(nrow, ncol, figsize=(6.6 * ncol, 5.1 * nrow + 0.4 * (nrow == 1)),
                            layout='constrained', squeeze=False)
    for ax, r in zip(axs.flat, rs):
        panel(ax, r)
    for ax in list(axs.flat)[len(rs):]:
        ax.set_visible(False)
    # 共享图例（按由内到外的层序）
    name = {'A': 'A-%s (Pb/Sn)', 'B': 'B-%s (Sb/Bi)', 'Te': 'Te-%s'}
    handles = [Line2D([], [], color=COL[lay], lw=6.0, label=name[lay[0]] % lay[1]) for lay in LAYERS]
    fig.legend(handles=handles, loc='outside lower center', ncol=6 if ncol > 1 else 3, frameon=False,
               fontsize=10.5, title_fontsize=9.5,
               title='each band = stacked layers (inner → outer as listed); layer thickness ∝ share of that element-orbital')
    fig.suptitle('HSE+SOC bands of A$_2$B$_2$Te$_5$: color = element-orbital, layer thickness = its share',
                 fontsize=13.5 if ncol > 1 else 11)
    stem = 'orbital_fatband_4mats' if mats == MATS else 'orbital_fatband_' + '_'.join(mats)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(a.out, stem + '.' + ext), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('written', os.path.join(a.out, stem + '.png'))


if __name__ == '__main__':
    main()
