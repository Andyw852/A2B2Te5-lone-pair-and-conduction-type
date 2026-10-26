#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
orbital_fatband.py —— 轨道分辨能带（fat-band，文献常用的气泡画法）：每个材料单独一张图

数据：step4_fatband_plot/{PROCAR_OPT 或 PROCAR}（路径 k 点的 lm 投影）+ step4_band_plot/band.dat（能量）
      PROCAR 的 SOC 每带 4 个块：第 1 块是“电荷”投影（后 3 块是自旋分量，不用于轨道成分）。
      每个块的 9 条 ion 行 = 每个离子的 s/py/pz/px/dxy/dyz/dz2/dxz/x2-y2 权重。

画法：
  能带画成细灰线；在每个 (k, 能带) 上按 5 个“元素-轨道”类别画气泡，气泡面积 ∝ 该类别份额：
      A-s 红、A-p 琥珀（圆点）；B-s 绿、B-p 青（菱形）；Te（s+p 合并）蓝（圆点）
  —— A=Pb/Sn，B=Sb/Bi；Te-s 在 VBM±3 eV 内平均仅 2–9%，并入 Te。
  配色取自经色觉校验的分类色板（全两两对比：正常视觉 ΔE ≥ 15.6）；红-绿、红-青两对对色弱者偏近，
  因此 A/B 再用圆点/菱形区分（二次编码）。
  所有气泡按面积从大到小统一排序后绘制，小份额（如 Sn-s 孤对）始终压在大气泡上面，不会被遮住。
  SOC + 反演对称 ⇒ Kramers 二重简并：每对只画一组气泡（能量、份额取两条的平均），细线仍画全部能带。
  路径断点（A|L、M|H）处不连线。
  Bi 系（零权重 HSE 路径方案）默认剔除落在 SCF 网格上的路径点（见 mesh_dup_points），能量零点、带隙按剩余点重算；
  --keep-dup-kpts 保留它们。

输出：每个材料一张 orbital_fatband_<材料>.png/pdf

数据位置（自动识别两种目录结构）：
  * clone 下来的仓库：默认读 12_能带绘图_20261021/源数据/<材料>/（band.dat + 投影_解析.npz），
    输出到 12_能带绘图_20261021/复现输出/；
  * 原工作区：<root>/<材料>/step4_band_plot + step4_fatband_plot（直接解析 PROCAR）。
  --procar-root 07_带边能窗/源数据 可改为从原始 PROCAR(.gz) 重新解析，并与 npz 逐元素比对。

用法：
  python 10_脚本/orbital_fatband.py
  python 10_脚本/orbital_fatband.py --materials Sn2Bi2Te5          # 只画指定材料
  python 10_脚本/orbital_fatband.py --procar-root 07_带边能窗/源数据   # 从原始 PROCAR 复算
  python 10_脚本/orbital_fatband.py --data-root <root> --out <dir>
  python 10_脚本/orbital_fatband.py --emin -2 --emax 2              # 改能窗
"""
import argparse, gzip, json, os, re, sys
from fractions import Fraction
from math import lcm
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle

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
GROUPS = [('A', 's'), ('A', 'p'), ('B', 's'), ('B', 'p'), ('Te', 'sp')]      # 图例顺序
COL = {('A', 's'): '#e34948', ('A', 'p'): '#eda100',        # 红 / 琥珀
       ('B', 's'): '#008300', ('B', 'p'): '#1baf7a',        # 绿 / 青
       ('Te', 'sp'): '#2a78d6'}                             # 蓝
MARK = {'A': 'o', 'B': 'D', 'Te': 'o'}
S_MAX = 46.0           # 份额 = 1 时的气泡面积 (pt²)
W_MIN = 0.03           # 份额低于此值不画气泡
K_STEP = 2             # 每隔几个 k 点画一个气泡（细线仍用全部点）
EWIN = (-3.0, 3.0)     # 能窗 (eV, 相对 VBM)
KEEP_DUP = False       # True: 不剔除落在 SCF 网格上的零权重路径点

plt.rcParams.update({'font.size': 11, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42,
                     'axes.linewidth': 0.9, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'ytick.right': True, 'ytick.minor.visible': True})


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


def mesh_dup_points(d, nk):
    """零权重 HSE 路径方案（Bi 系）：KPOINTS = 带权重的 SCF 网格 + nk 个零权重路径点。
    落在 SCF 网格上的路径点，Fock 项 q→0 的处理与其它路径点不同，价带被整体抬高 0.03–0.47 eV
    （能带上的单点尖峰）。返回 (这些路径点的序号 = band.dat 行号, 路径点分数坐标)；
    KPOINTS_OPT 方案（Sb 系）返回 ([], None)。"""
    p = os.path.join(d, 'KPOINTS')
    if not os.path.exists(p):
        return [], None
    rows = []
    for line in open(p).read().splitlines()[3:]:
        q = line.split()
        try:
            rows.append([float(x) for x in q[:4]])
        except ValueError:
            continue
    rows = np.array([r for r in rows if len(r) == 4])
    if rows.ndim != 2:
        return [], None
    path, grid = rows[rows[:, 3] == 0, :3], rows[rows[:, 3] > 0, :3]
    if len(path) != nk or not len(grid):
        return [], None
    N = np.array([lcm(*(Fraction(x).limit_denominator(48).denominator for x in grid[:, i])) for i in range(3)])
    on = np.all(np.abs(path * N - np.round(path * N)) < 1e-4, axis=1)
    return [int(i) for i in np.where(on)[0]], path


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
        # 校验：扣掉 E_VBM 后作图能窗（±3.2 eV）内逐点一致（能窗外最高几条空带在投影计算里未收敛，偏差大但不画）
        de = np.abs(E - np.median(E - Edat) - Edat)
        inw = np.abs(Edat) <= 3.2
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
    # 剔除落在 SCF 网格上的路径点（Bi 系），能量零点移到剩余点的 VBM，带隙重算
    seg = np.concatenate([[0], np.cumsum(np.diff(kdist) < 1e-6)])   # 路径段号，断点处不连线
    drop, path = ([], None) if KEEP_DUP else mesh_dup_points(bdir, len(kdist))
    if drop:
        keep = np.setdiff1d(np.arange(len(kdist)), drop)
        kdist, Edat, Wc, frac, tot, seg, path = (x[keep] for x in (kdist, Edat, Wc, frac, tot, seg, path))
        vb, cb = Edat[:, nv - 1], Edat[:, nv]
        iv, ic = int(np.argmax(vb)), int(np.argmin(cb))
        Edat = Edat - vb[iv]
        sm = dict(sm, gap_eV=float(cb[ic] - vb[iv]),
                  gap_type='direct' if np.allclose(np.mod(path[iv] + 1e-9, 1), np.mod(path[ic] + 1e-9, 1), atol=1e-4)
                  else 'indirect')
        print('     剔除 %d 个落在 SCF 网格上的路径点；VBM 移到 k=%.3f %s（原零点下 %.3f eV），带隙 %.3f eV (%s)' % (
            len(drop), kdist[iv], np.round(path[iv], 3), vb[iv], sm['gap_eV'], sm['gap_type']))
    return dict(m=m, kdist=kdist, E=Edat, Wc=Wc, frac=frac, tot=tot, seg=seg, dropped=drop,
                labels=labels, sm=sm, nv=nv, elements=elements, cats=cats)


def elem_name(r, role):
    return [e for e in r['elements'] if EROLE[e] == role][0]


def group_frac(r):
    """(nk, nb, 5)：按 GROUPS 归并的份额（Te-s + Te-p → Te）。"""
    cats = [c[1:] for c in r['cats']]
    out = []
    for role, orb in GROUPS:
        idx = [i for i, (ro, ob) in enumerate(cats) if ro == role and ob in orb]
        out.append(r['frac'][:, :, idx].sum(axis=-1))
    return np.stack(out, axis=-1)


def marker_path(m):
    ms = MarkerStyle(m)
    return ms.get_path().transformed(ms.get_transform())


def plot_material(r, path):
    m = r['m']; k = r['kdist']; E = r['E']; sm = r['sm']
    nk, nb = E.shape
    fig, ax = plt.subplots(figsize=(5.4, 4.6), layout='constrained')
    # 细线：全部能带；路径断点（相邻两点 k 距离相同）处不连线
    jj = np.where(np.diff(r['seg']) == 0)[0]
    bands = [b for b in range(nb) if E[:, b].max() >= EWIN[0] - 0.5 and E[:, b].min() <= EWIN[1] + 0.5]
    Eb = E[:, bands]
    kb = np.broadcast_to(k[:, None], Eb.shape)
    seg = np.stack([np.stack([kb[jj], Eb[jj]], -1), np.stack([kb[jj + 1], Eb[jj + 1]], -1)], axis=2).reshape(-1, 2, 2)
    ax.add_collection(LineCollection(seg, colors='0.45', linewidths=0.55, zorder=1))
    # 气泡：每个 Kramers 对一组（能量、份额取平均），每隔 K_STEP 个 k 点
    G = group_frac(r)
    npair = nb // 2
    Ep = 0.5 * (E[:, 0:2 * npair:2] + E[:, 1:2 * npair:2])
    Gp = 0.5 * (G[:, 0:2 * npair:2] + G[:, 1:2 * npair:2])
    ik = np.arange(0, nk, K_STEP)
    xs, ys, ss, cs, ps = [], [], [], [], []
    for g, (role, orb) in enumerate(GROUPS):
        w = Gp[ik, :, g]
        y = Ep[ik]
        sel = (w >= W_MIN) & (y >= EWIN[0] - 0.2) & (y <= EWIN[1] + 0.2)
        ii, bb = np.nonzero(sel)
        xs.append(k[ik][ii]); ys.append(y[ii, bb]); ss.append(S_MAX * w[ii, bb])
        cs += [COL[(role, orb)]] * len(ii); ps += [MARK[role]] * len(ii)
    xs, ys, ss = map(np.concatenate, (xs, ys, ss))
    order = np.argsort(-ss, kind='stable')          # 大气泡先画，小气泡压在上面
    cs = [cs[i] for i in order]
    edge = [tuple(0.62 * np.array(matplotlib.colors.to_rgb(c))) for c in cs]
    paths = {mk: marker_path(mk) for mk in set(MARK.values())}
    sc = ax.scatter(xs[order], ys[order], s=ss[order], c=cs, edgecolors=edge, linewidths=0.35,
                    alpha=0.9, zorder=2)
    sc.set_paths([paths[ps[i]] for i in order])
    # 高对称点、VBM
    xt = [x for _, x in r['labels']]
    names = [n for n, _ in r['labels']]
    # 末段 A|L、M|H… 间距很短：相邻标签按字符数估计会重叠时，后一个下移一行
    span = k[-1] - k[0]
    for i in range(1, len(xt)):
        if not names[i - 1].startswith('\n') and \
                xt[i] - xt[i - 1] < 0.02 * span * (len(names[i - 1]) + len(names[i])) / 2:
            names[i] = '\n' + names[i]
    for x in xt:
        ax.axvline(x, color='0.6', lw=0.6, zorder=0)
    ax.set_xticks(xt, names)
    ax.tick_params(axis='x', length=0, labelsize=10, pad=5)
    ax.axhline(0, color='0.25', lw=0.7, ls='--', zorder=0)
    ax.set_xlim(xt[0], xt[-1]); ax.set_ylim(*EWIN)
    ax.set_ylabel('E $-$ E$_{\\rm VBM}$ (eV)')
    ax.set_title(re.sub(r'(\d)', r'$_\1$', m), loc='left', fontsize=12, fontweight='bold')
    ax.set_title('HSE+SOC, %s gap %.3f eV' % (sm['gap_type'], sm['gap_eV']), loc='right', fontsize=9.5,
                 color='0.3')
    # 图例：真实元素名；标记大小固定（面积 ∝ 份额见说明）
    A, B = elem_name(r, 'A'), elem_name(r, 'B')
    lab = {('A', 's'): A + '-s', ('A', 'p'): A + '-p', ('B', 's'): B + '-s', ('B', 'p'): B + '-p',
           ('Te', 'sp'): 'Te (s+p)'}
    handles = [Line2D([], [], ls='none', marker=MARK[g[0]], ms=8, mfc=COL[g],
                      mec=tuple(0.62 * np.array(matplotlib.colors.to_rgb(COL[g]))), mew=0.6, label=lab[g])
               for g in GROUPS]
    fig.legend(handles=handles, loc='outside lower center', ncol=5, frameon=False, fontsize=10,
               handletextpad=0.2, columnspacing=1.1,
               title='marker area ∝ share of each element-orbital (Kramers pairs drawn once)'
                     + ('\n%d path k-points lying on the SCF k-mesh removed (HSE zero-weight artifact)' % len(r['dropped'])
                        if r['dropped'] else ''),
               title_fontsize=8.5)
    for ext in ('png', 'pdf'):
        fig.savefig(path + '.' + ext, dpi=300)
    plt.close(fig)
    print('written', path + '.png')


def main(argv=None):
    global EWIN, KEEP_DUP
    ap = argparse.ArgumentParser(description='HSE+SOC 轨道分辨 fat-band（气泡面积 ∝ 元素-轨道份额），每个材料一张图')
    ap.add_argument('--data-root', default=SRC,
                    help='含 <材料>/band.dat… 的目录（默认 %(default)s）')
    ap.add_argument('--procar-root', default=None,
                    help='从原始 PROCAR(.gz) 重新解析投影，如 07_带边能窗/源数据；默认用 投影_解析.npz')
    ap.add_argument('--out', default=OUT, help='输出目录（默认 %(default)s）')
    ap.add_argument('--materials', nargs='+', default=MATS, choices=MATS, metavar='MAT',
                    help='要画的材料（默认四个：%s）' % ' '.join(MATS))
    ap.add_argument('--emin', type=float, default=EWIN[0], help='能窗下限 (eV，相对 VBM，默认 %(default)s)')
    ap.add_argument('--emax', type=float, default=EWIN[1], help='能窗上限 (eV，相对 VBM，默认 %(default)s)')
    ap.add_argument('--keep-dup-kpts', action='store_true',
                    help='不剔除 Bi 系落在 SCF 网格上的零权重路径点')
    a = ap.parse_args(argv)
    EWIN, KEEP_DUP = (a.emin, a.emax), a.keep_dup_kpts
    os.makedirs(a.out, exist_ok=True)
    for m in MATS:
        if m in a.materials:
            plot_material(load_material(m, a.data_root, a.procar_root), os.path.join(a.out, 'orbital_fatband_' + m))


if __name__ == '__main__':
    main()
