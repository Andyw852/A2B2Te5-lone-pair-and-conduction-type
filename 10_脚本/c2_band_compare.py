"""c2 路径能带：HSE06+SOC 与 PBE+SOC 叠画（各自以 VBM 为零点）。

数据：22_能带基准_1-0/<材料>/c2_{pbe,hse}/path_eigen.txt（k_index band_index energy_eV，k 从 0 起），
路径：22_能带基准_1-0/c2_path/KPOINTS_OPT（四种材料共用的 204 点），晶格：08_介电复核/源数据/<材料>/POSCAR。
Pb2Sb2Te5 的数据在 22_能带基准_1-0/ 根目录下（c2_pbe/、c2_hse_v2/ 或 c2_hse/）。

用法：python 10_脚本/c2_band_compare.py Sn2Sb2Te5 [--nelect 68] [--emin -0.8 --emax 0.8] [--out 图.png]
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAB = {0: 'Γ', 33: 'M', 52: 'K', 90: 'Γ', 97: 'A', 130: 'L', 149: 'H', 187: 'A', 188: 'L', 195: 'M', 196: 'H', 203: 'K'}
BREAKS = [(187, 188), (195, 196)]              # 相邻两点之间不连续（A|L、M|H）
SUB = {m: m.replace('2', '$_2$').replace('5', '$_5$') for m in ('Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5')}


def read_lat(p):
    L = open(p).read().splitlines(); s = float(L[1])
    return s * np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)])


def read_kpoints_opt(p):
    L = open(p).read().splitlines(); n = int(L[1])
    return np.array([[float(x) for x in l.split()[:3]] for l in L[3:3 + n]])


def read_path_eigen(p):
    d = np.loadtxt(p, comments='#'); nk = int(d[:, 0].max()) + 1; nb = int(d[:, 1].max())
    E = np.full((nk, nb), np.nan); E[d[:, 0].astype(int), d[:, 1].astype(int) - 1] = d[:, 2]
    return E


def kdist(K, B):
    x = [0.0]
    for i in range(1, len(K)):
        x.append(x[-1] + (0.0 if (i - 1, i) in BREAKS else np.linalg.norm((K[i] - K[i - 1]) @ B)))
    return np.array(x)


def segments(n):
    """按 BREAKS 切开的连续段 [(起, 止)]，止为闭区间端点。"""
    cuts = [0] + [b for _, b in BREAKS] + [n]
    return [(cuts[i], cuts[i + 1] - 1) for i in range(len(cuts) - 1)]


def nelect_from_eigenval(p):
    return int(open(p).read().split('\n')[5].split()[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('material')
    ap.add_argument('--nelect', type=int)
    ap.add_argument('--hse-dir', help='默认依次找 c2_hse_v2、c2_hse')
    ap.add_argument('--emin', type=float, default=-0.8)
    ap.add_argument('--emax', type=float, default=0.8)
    ap.add_argument('--out')
    a = ap.parse_args()

    base = os.path.join(ROOT, '22_能带基准_1-0', '' if a.material == 'Pb2Sb2Te5' else a.material)
    pbe = os.path.join(base, 'c2_pbe')
    hse = a.hse_dir or next((os.path.join(base, d) for d in ('c2_hse_v2', 'c2_hse')
                             if os.path.exists(os.path.join(base, d, 'path_eigen.txt'))), None)
    if hse is None:
        sys.exit(f'找不到 {a.material} 的 c2_hse*/path_eigen.txt')
    nel = a.nelect or nelect_from_eigenval(os.path.join(pbe, 'EIGENVAL'))
    K = read_kpoints_opt(os.path.join(ROOT, '22_能带基准_1-0', 'c2_path', 'KPOINTS_OPT'))
    B = 2 * np.pi * np.linalg.inv(read_lat(os.path.join(ROOT, '08_介电复核', '源数据', a.material, 'POSCAR'))).T
    x = kdist(K, B)

    for f in font_manager.findSystemFonts():
        if 'wqy' in f.lower() or 'notosanscjk' in f.lower().replace('-', ''):
            font_manager.fontManager.addfont(f)
            plt.rcParams['font.family'] = font_manager.FontProperties(fname=f).get_name(); break
    plt.rcParams['mathtext.fontset'] = 'dejavusans'
    plt.rcParams['axes.unicode_minus'] = False

    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
    style = {'PBE+SOC': dict(color='#E08A2C', lw=1.1, ls=(0, (4, 2)), zorder=2),
             'HSE06+SOC': dict(color='#1F4E9C', lw=1.5, ls='-', zorder=3)}
    summary = []
    for name, d in (('PBE+SOC', pbe), ('HSE06+SOC', hse)):
        E = read_path_eigen(os.path.join(d, 'path_eigen.txt'))
        if E.shape[0] != len(K):
            sys.exit(f'{d}: 路径点数 {E.shape[0]} ≠ KPOINTS_OPT {len(K)}')
        v, c = E[:, nel - 1], E[:, nel]
        ev = np.nanmax(v); iv, ic = np.nanargmax(v), np.nanargmin(c); gap = np.nanmin(c) - ev
        for b in range(max(0, nel - 12), min(E.shape[1], nel + 12)):
            for s, t in segments(len(K)):
                ax.plot(x[s:t + 1], E[s:t + 1, b] - ev, **style[name],
                        label=f'{name}（$E_g$ = {1000 * gap:.0f} meV）' if (b, s) == (nel - 1, 0) else None)
        mk = dict(ms=6, mfc='white', mew=1.4, mec=style[name]['color'], zorder=5)
        ax.plot(x[iv], 0, 'o', **mk); ax.plot(x[ic], c[ic] - ev, 's', **mk)
        summary.append(f'{name}: VBM #{iv} {np.round(K[iv], 3).tolist()}  CBM #{ic} {np.round(K[ic], 3).tolist()}  E_g {gap:.4f} eV')

    ticks, labels = [], []
    for i in sorted(LAB):
        if (i - 1, i) in BREAKS:
            labels[-1] = f'{labels[-1]}|{LAB[i]}'; continue
        ticks.append(x[i]); labels.append(LAB[i])
    for t in ticks[1:-1]:
        ax.axvline(t, color='#B0B0B0', lw=0.6, zorder=1)
    ax.axhline(0, color='#808080', lw=0.6, ls=':', zorder=1)
    ax.set_xticks(ticks, labels, fontsize=9); ax.set_xlim(x[0], x[-1]); ax.set_ylim(a.emin, a.emax)
    ax.set_ylabel(r'$E - E_\mathrm{VBM}$ (eV)')
    ax.set_title(f'{SUB.get(a.material, a.material)}：c2 路径能带（各自以 VBM 为零点；○ VBM，□ CBM）', fontsize=10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07), ncol=2, fontsize=8.5, frameon=False)
    fig.tight_layout()
    out = a.out or os.path.join(base, f'c2_能带_HSE_vs_PBE_{a.material}.png')
    fig.savefig(out); print('\n'.join(summary)); print('→', out)


if __name__ == '__main__':
    main()
