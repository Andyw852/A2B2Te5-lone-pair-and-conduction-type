"""交接样稿（2026-10-29）：价带口袋在布里渊区中的位置，以及“Γ 单谷 → 孤对电子多谷”转变的能量。

只读已有数据（22_能带基准_1-0/<材料>/c2_{pbe,hse}/path_eigen.txt、PROCAR_OPT.gz，A2B2Te5_SOC_dense/materials/<材料>/POSCAR），
不做新计算。样式沿用 paper_figures.py。这是给接手者的样稿，正式的图 6 由接手者并入 paper_figures.py。

输出到 21_阶段报告与交接/交接_口袋示意图_样稿/：
  样稿A_口袋_布里渊区.png/.pdf     (a) Sn2Sb2Te5、(b) Sn2Bi2Te5 沿 c 俯视；(c) Sn2Bi2Te5 三维透视（k_z 放大）
  样稿B_转变能量.png/.pdf          ΔE = E(最高孤对电子谷) − E(Γ–A 轴上最高价带态)，四种材料 × PBE/HSE
  pocket_positions.csv             每个口袋一行（Sn2Sb2Te5、Sn2Bi2Te5，HSE06+SOC）
  valley_offsets.csv               转变能量表（四种材料 × 两种泛函）

用法：python3 10_脚本/pocket_preview_20261029.py
"""
import csv, os, sys
import numpy as np
import spglib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import paper_figures as pf
import report_figures_20261028 as rf
import c2_valleys as cv
import multivalley_seebeck as mv

OUT = os.path.join(ROOT, '21_阶段报告与交接', '交接_口袋示意图_样稿')
KT = 1000 * mv.KB * 300                     # 300 K 的 k_BT（meV）
AXIS = range(90, 98)                        # 路径上的 Γ–A 段（第二个 Γ 到 A）
C_A, C_B = pf.C_A, pf.C_B


# ------------------------------------------------------------------ 数据
def lattice(m):
    lat, frac, sp = pf.read_poscar(os.path.join(ROOT, 'A2B2Te5_SOC_dense', 'materials', m, 'POSCAR'))
    num = [sorted(set(sp)).index(e) for e in sp]
    rots = spglib.get_symmetry((lat, frac, num), symprec=1e-3)['rotations']
    rots = np.unique(rots.reshape(len(rots), -1), axis=0).reshape(-1, 3, 3)
    assert len(rots) == 12, len(rots)                     # D3d
    return 2 * np.pi * np.linalg.inv(lat).T, rots


def star(k, rots):
    out = []
    for R in rots:
        for s in (1, -1):                                 # 含时间反演（有反演时不增加新点）
            q = s * (np.asarray(k, float) @ R); q -= np.round(q); q[np.isclose(q, -0.5)] = 0.5
            if not any(np.allclose(q, p, atol=1e-6) for p in out): out.append(q)
    return out


def ws(q, B):
    """把分数坐标 q 换到第一布里渊区（面内取最短的等价像），返回笛卡尔坐标（Å⁻¹）。"""
    c = [(q + np.array([i, j, 0])) @ B for i in (-1, 0, 1) for j in (-1, 0, 1)]
    return min(c, key=lambda v: np.hypot(v[0], v[1]))


def top_vb(m, side):
    """路径上最高价带（Kramers 对取平均，与图 5、表 5 一致），相对价带顶，meV。"""
    E, _ = pf.bands(m, side); nel = pf.NEL[m]
    return 1000 * E[:, nel - 2:nel].mean(1)


def offsets():
    """每种材料、每种泛函：Γ 处、Γ–A 轴上最高、最高孤对电子谷（与图 6 同一判据）的能量与位置。"""
    K = pf.kpath(); rows = []
    for m in pf.MATS:
        _, rots = lattice(m)
        for side, fn in (('c2_pbe', 'PBE+SOC'), ('c2_hse', 'HSE06+SOC')):
            v = top_vb(m, side)
            ia = max(AXIS, key=lambda i: v[i])                       # Γ–A 轴上最高的价带态
            lp = rf.lone_pair_offset(m, side); il = lp[1]
            dE = v[il] - v[ia]
            rows.append(dict(material=m, functional=fn, E_Gamma=v[0], E_axis=v[ia], i_axis=ia, k_axis=K[ia],
                             n_axis=len(star(K[ia], rots)), E_LP=v[il], i_LP=il, k_LP=K[il], n_LP=len(star(K[il], rots)),
                             s_LP=lp[2], dE=dE, dE_Gamma=v[il] - v[0]))
    return rows


# ------------------------------------------------------------------ 样稿 A：口袋
def pockets(m):
    """HSE06+SOC：Γ–A 轴上最高态与三组孤对电子谷，展开成全部等价口袋。"""
    B, rots = lattice(m); K = pf.kpath(); v = top_vb(m, 'c2_hse')
    ia = max(AXIS, key=lambda i: v[i])
    groups = [('Γ–A' if ia != 90 else 'Γ', ia, C_B)] + [(nm, i, C_A) for nm, i in pf.SNVAL.items() if nm != 'Γ']
    out = []
    for nm, i, col in groups:
        pts = star(K[i], rots)
        out.append(dict(name=nm, i=i, E=v[i], col=col, frac=pts, cart=[ws(p, B) for p in pts]))
    return B, out


def hexagon(B):
    K0 = np.array([-2 / 3, 1 / 3, 0]) @ B; ang = np.arctan2(K0[1], K0[0]); r0 = np.hypot(*K0[:2])
    a = ang + np.arange(7) * np.pi / 3
    return r0 * np.cos(a), r0 * np.sin(a)


def near(E):
    return E >= -2 * KT                                   # 价带顶 2k_BT 以内：画大、画实


def top_view(ax, m, B, P, letter):
    hx, hy = hexagon(B)
    ax.plot(hx, hy, color='k', lw=0.8, zorder=1)
    r0 = 0.5 * np.linalg.norm(B[0][:2]); a0 = np.arctan2(B[0][1], B[0][0])
    M = [(r0 * np.cos(a0 + t * np.pi / 3), r0 * np.sin(a0 + t * np.pi / 3)) for t in range(6)]   # 6 个 M，在边中点
    for q in M: ax.plot([0, q[0]], [0, q[1]], color='#BDBDBD', lw=0.5, ls=pf.DASH, zorder=0)
    for g in P:
        big = near(g['E']); ms = 5.6 if big else 3.4; al = 1.0 if big else 0.45
        for c, f in zip(g['cart'], g['frac']):
            top = np.isclose(abs(f[2]), 0.5)
            ax.plot(c[0], c[1], 'o', ms=ms, mew=0.9, color=g['col'], mfc='white' if top else g['col'], alpha=al, zorder=3)
    Mq = max(M, key=lambda c: c[0] + 0.1 * c[1])
    ax.text(Mq[0] + 0.04, Mq[1], 'M (L)', ha='left', va='center', fontsize=6.5)
    kv = max(zip(hx, hy), key=lambda t: (round(t[1], 3), t[0]))
    ax.text(kv[0] + 0.03, kv[1] + 0.02, 'K (H)', ha='left', va='bottom', fontsize=6.5)
    ax.annotate('Γ (A)', (0, 0), (-0.62, -0.50), fontsize=6.5, ha='center', va='center',
                arrowprops=dict(arrowstyle='-', lw=0.5, color='#555555', shrinkA=1, shrinkB=3))
    # 右下角小表：每组的能量（meV，相对价带顶）与等价口袋数
    y = -0.66
    for g in P:
        top = all(np.isclose(abs(f[2]), 0.5) for f in g['frac'])
        al = 1.0 if near(g['E']) else 0.45
        ax.plot(0.74, y, 'o', ms=3.4, mew=0.8, color=g['col'], mfc='white' if top else g['col'], alpha=al)
        e = f"{g['E']:.0f}".replace('-', '−')
        ax.text(0.81, y, f"{g['name']}", ha='left', va='center', fontsize=6.0)
        ax.text(1.66, y, f"{e}  ×{len(g['frac'])}", ha='right', va='center', fontsize=6.0)
        y -= 0.115
    n = sum(len(g['frac']) for g in P if near(g['E']))
    ax.set_title(f"{pf.TEX[m]}: {n} pocket{'s' if n > 1 else ''}", fontsize=8, pad=2)
    ax.set_aspect('equal'); ax.set_xlim(-1.0, 1.68); ax.set_ylim(-1.12, 0.98); ax.axis('off')
    ax.text(0.0, 1.06, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', ha='left', va='top')


def persp(ax, m, B, P, zs):
    hx, hy = hexagon(B); h = 0.5 * np.linalg.norm(B[2]) * zs
    for z in (-h, h): ax.plot(hx, hy, [z] * 7, color='k', lw=0.6)
    for x, y in zip(hx[:-1], hy[:-1]): ax.plot([x, x], [y, y], [-h, h], color='k', lw=0.6)
    ax.plot([0, 0], [0, 0], [-h, h], color='#BDBDBD', lw=0.5, ls=pf.DASH)
    for g in P:
        big = near(g['E']); ms = 4.6 if big else 3.0; al = 1.0 if big else 0.45
        for c in g['cart']:
            if np.isclose(abs(c[2]), h / zs, atol=1e-6):
                zz = (h, -h)                                # 区界面上的口袋：顶面、底面是同一个，两处都画
            else:
                zz = (c[2] * zs, -c[2] * zs) if g['name'] == 'Γ–A' else (c[2] * zs,)
            for z in zz:
                ax.plot([c[0]], [c[1]], [z], 'o', ms=ms, color=g['col'], mfc=g['col'] if abs(z) < h - 1e-9 else 'white',
                        mew=0.9, alpha=al)
    ax.text(-0.10, 0, h + 0.06, 'A', fontsize=6.5, ha='center', va='bottom'); ax.text(0.12, -0.10, -0.08, 'Γ', fontsize=6.5, ha='left', va='top')
    ax.view_init(elev=22, azim=-62); ax.set_box_aspect((1, 1, 2 * h / 1.9))
    ax.set_axis_off()
    ax.set_title(f'{pf.TEX[m]}: 3D view ($k_z$ ×{zs:g})', fontsize=8, pad=-14)


def figA():
    data = {m: pockets(m) for m in ('Sn2Sb2Te5', 'Sn2Bi2Te5')}
    fig = plt.figure(figsize=(pf.W2, 2.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05], wspace=0.04)
    for k, m in enumerate(data):
        top_view(fig.add_subplot(gs[0, k]), m, *data[m], f'({"ab"[k]})')
    ax = fig.add_subplot(gs[0, 2], projection='3d'); persp(ax, 'Sn2Bi2Te5', *data['Sn2Bi2Te5'], zs=4)
    ax.text2D(0.0, 1.0, '(c)', transform=ax.transAxes, fontsize=9, fontweight='bold', ha='left', va='top')
    h = [Line2D([], [], ls='', marker='o', ms=5, color=C_B, label=r'B-site p, on $\Gamma$–A axis'),
         Line2D([], [], ls='', marker='o', ms=5, color=C_A, label='Sn 5s lone pair'),
         Line2D([], [], ls='', marker='o', ms=5, color='#555555', label=r'filled: $k_z$ = 0 plane (or inside)'),
         Line2D([], [], ls='', marker='o', ms=5, color='#555555', mfc='white', label=r'open: $k_z$ = ½ zone face')]
    h.append(Line2D([], [], ls='', marker='o', ms=3.4, color=C_A, alpha=0.45, label=r'faded: > $2k_\mathrm{B}T$ below VBM'))
    fig.legend(handles=h, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.13), handletextpad=0.2, columnspacing=0.9)
    for e in ('png', 'pdf'): fig.savefig(os.path.join(OUT, f'样稿A_口袋_布里渊区.{e}'))
    plt.close(fig)
    rows = []
    for m, (B, P) in data.items():
        for g in P:
            for f, c in zip(g['frac'], g['cart']):
                rows.append([m, g['name'], g['i'], *[pf.r(x, 4) for x in f], *[pf.r(x, 4) for x in c],
                             'kz=1/2 face' if np.isclose(abs(f[2]), 0.5) else ('kz=0' if np.isclose(f[2], 0) else 'Γ–A axis'),
                             pf.r(g['E'], 1), len(g['frac'])])
    write('pocket_positions.csv', ['material', 'valley', 'path_k_index', 'k_frac_1', 'k_frac_2', 'k_frac_3',
                                   'kx_1/A', 'ky_1/A', 'kz_1/A', 'plane', 'E-EVBM_meV(HSE06+SOC)', 'n_equivalent'], rows)
    return {m: [(g['name'], len(g['frac']), round(g['E'], 1)) for g in P] for m, (B, P) in data.items()}


# ------------------------------------------------------------------ 样稿 B：转变能量
def figB(rows):
    fig, ax = plt.subplots(figsize=(pf.W1, 2.3))
    ax.axvspan(-KT, KT, color='#EEF3F8', zorder=0, lw=0)
    ax.axvline(0, color='k', lw=0.6, zorder=1)
    by = {(d['material'], d['functional']): d for d in rows}
    for j, m in enumerate(pf.MATS[::-1]):
        p, h = by[(m, 'PBE+SOC')]['dE'], by[(m, 'HSE06+SOC')]['dE']
        ax.plot([p, h], [j, j], color='#BDBDBD', lw=1.2, zorder=1)
        ax.plot(p, j + 0.08, 'o', ms=4.5, color=pf.C_PBE, mec='white', mew=0.6, zorder=3, label='PBE+SOC' if j == 0 else None)
        ax.plot(h, j - 0.08, 'o', ms=4.5, color=pf.C_HSE, mec='white', mew=0.6, zorder=3, label='HSE06+SOC' if j == 0 else None)
        for v, dy, col in ((h, -0.30, pf.C_HSE), (p, 0.30, pf.C_PBE)):
            ax.text(v, j + dy, f'{v:+.0f}'.replace('-', '−'), ha='center', va='center', fontsize=6.0, color=col)
    ax.set_yticks(range(4), [pf.TEX[m] for m in pf.MATS[::-1]]); ax.set_ylim(-0.6, 4.2)
    ax.set_xlim(-262, 135); ax.tick_params(axis='y', length=0, right=False)
    ax.set_xlabel(r'$\Delta E = E_{\mathrm{lone\ pair}} - E_{\Gamma\mathrm{-A}}$ (meV)')
    ax.text(-130, 3.70, r'$\Gamma$–A valley on top' '\n' r'($N_\mathrm{v}$ = 1–2)', ha='center', va='center', fontsize=6.3)
    ax.text(80, 3.70, 'lone-pair valleys' '\n' 'on top (×6 each)', ha='center', va='center', fontsize=6.3)
    ax.text(0, -0.48, r'$\pm k_\mathrm{B}T$', ha='center', va='center', fontsize=6.0)
    ax.legend(loc='lower left', handletextpad=0.2, borderaxespad=0.2, fontsize=6.5)
    for e in ('png', 'pdf'): fig.savefig(os.path.join(OUT, f'样稿B_转变能量.{e}'))
    plt.close(fig)


def write(name, header, rows):
    with open(os.path.join(OUT, name), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def main():
    os.makedirs(OUT, exist_ok=True); pf.style()
    rows = offsets()
    write('valley_offsets.csv',
          ['material', 'functional', 'E_Gamma_meV', 'E_GammaA_axis_max_meV', 'axis_path_k_index', 'axis_k_frac_3', 'axis_n_equivalent',
           'E_lone_pair_max_meV', 'LP_path_k_index', 'LP_k_frac_1', 'LP_k_frac_3', 'LP_n_equivalent', 'LP_cation_s_percent',
           'dE_LP-axis_meV', 'dE_LP-Gamma_meV'],
          [[d['material'], d['functional'], pf.r(d['E_Gamma'], 1), pf.r(d['E_axis'], 1), d['i_axis'], pf.r(d['k_axis'][2], 4), d['n_axis'],
            pf.r(d['E_LP'], 1), d['i_LP'], pf.r(d['k_LP'][0], 4), pf.r(d['k_LP'][2], 4), d['n_LP'], pf.r(100 * d['s_LP'], 1),
            pf.r(d['dE'], 1), pf.r(d['dE_Gamma'], 1)] for d in rows])
    for d in rows:
        print(f"{d['material']:10s} {d['functional']:10s} Γ {d['E_Gamma']:7.1f}  axis {d['E_axis']:7.1f} (#{d['i_axis']}, ×{d['n_axis']})  "
              f"LP {d['E_LP']:7.1f} (#{d['i_LP']}, ×{d['n_LP']})  ΔE {d['dE']:+7.1f}  ΔE(Γ) {d['dE_Gamma']:+7.1f}")
    figB(rows)
    print(figA())
    print('→', OUT)


if __name__ == '__main__':
    main()
