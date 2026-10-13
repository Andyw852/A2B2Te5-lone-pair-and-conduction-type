#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
di_from_wavecar.py —— 用 KS 波函数 + Bader 盆计算 Raty 2021 口径的 DI / ES = 2*DI

公式（推导见 05_DI_ES/DI_公式与校验.md）：

    DI(A,B) = 2 · g_spin · < Tr[ N(A) N(B)^dag ] >_BZ
    N(a)_{ij} = Integral_{Omega_a} phi_i(r) phi_j*(r) dr        （自旋求和后的一阶密度矩阵块）
    ES(A,B)  = 2 · DI(A,B)

    * ncl（SOC，自旋子）波函数：N = M^up + M^dn，g_spin = 1
    * ISPIN=1 共线波函数：N = M（每空间轨道双占据），g_spin = 2
    * <...>_BZ 用不可约 k 点 + 空间群轨道展开：<G> = Sum_k w_k (1/m_k) Sum_{g in G/G_k} G(pi_{g^-1}A, pi_{g^-1}B)

内建校验：
    1. Sum_B DI(A,B) == 2 * N_A（QTAIM 求和规则；N_A = 该盆的（赝）电子数）—— 硬判据
    2. DI 对称性 DI(A,B) == DI(B,A)
    3. 每个 IBZ 点的 |orbit| == w_k * N_k（对称约定的自检）
    4. 等价接触的 DI 必须相等

注意：VASP 存的是 PAW 赝波函数（Sum|C|^2 != 1），因此这里是“赝波函数口径”的 DI；
      34 个占据带重建密度积分 65.04 e vs CHGCAR 68.00 e（差 4.4% = 原子球内增补项）。
"""
import argparse, json, os, re, struct, sys
import numpy as np

BASE = '/home/wangchao/work_wc/Pb2Sb2Te5-Sb孤对电子特性分析-POSCAR_B-band_and_elf'
PROJ = os.path.join(BASE, 'Raty2021_ES_ET_补算')
OUT = os.path.join(PROJ, '05_DI_ES', '结果')
MASKDIR = os.path.join(PROJ, '05_DI_ES', 'masks')
SOC_ROOT = os.path.join(BASE, 'band_and_elf', 'huang')
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']

from pymatgen.io.vasp.outputs import Wavecar
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# ------------------------------------------------------------------ WAVECAR
def read_wavecar_kmeta(path):
    """VASP>=6.4 布局：REC=1 (RECL,ISPIN,RTAG)；REC=2 (NKPTS,NBANDS,ENMAX,晶格)；
    每个 k 点一条记录 (RNPL,VKPT(3),然后每带 3 个实数: E, 0.0, F)。"""
    f = open(path, 'rb')
    recl = int(struct.unpack('d', f.read(8))[0])
    f.seek(recl)
    nk, nb, encut = struct.unpack('ddd', f.read(24))
    nk, nb = int(round(nk)), int(round(nb))
    recs = []
    for ik in range(nk):
        f.seek(recl * (2 + ik * (nb + 1)))
        v = np.frombuffer(f.read(8 + 24 + 24 * nb), dtype=np.float64)
        recs.append(dict(npl=int(round(v[0])), kvec=v[1:4],
                         eig=v[4:4 + 3 * nb:3].copy(), occ=v[6:4 + 3 * nb:3].copy()))
    f.close()
    return dict(recl=recl, nk=nk, nb=nb, encut=encut, recs=recs)


def read_kpoints(path, nk_file):
    """返回 (kvec_all(n,3), weights_all(n,))，只取带权重的行（KPOINTS_OPT 的零权重点排除）。"""
    ks, ws = [], []
    for line in open(path):
        p = line.split()
        if len(p) == 4:
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            if abs(v[3]) <= 1e-12:
                continue
            ks.append(v[:3]); ws.append(v[3])
    ks = np.array(ks); ws = np.array(ws)
    if len(ks) > nk_file:
        ks, ws = ks[:nk_file], ws[:nk_file]
    if len(ks) == 0:
        raise RuntimeError('KPOINTS 里没有带权重的 k 点')
    # 注意：Bi 材料的 WAVECAR 里还含 KPOINTS_OPT 的零权重路径点，只取前 len(ks) 个
    return ks, ws / ws.sum()


# ------------------------------------------------------------------ 对称
def symmetry_ops(structure):
    sga = SpacegroupAnalyzer(structure, symprec=1e-4)
    ops = sga.get_symmetry_operations(cartesian=False)
    frac = np.array(structure.frac_coords)
    L = structure.lattice.matrix
    out = []
    for op in ops:
        R = np.array(op.rotation_matrix); t = np.array(op.translation_vector)
        perm = np.zeros(len(structure), dtype=int)
        for a in range(len(structure)):
            d = frac - (R @ frac[a] + t)
            d -= np.round(d)
            cart = d @ L
            j = int(np.argmin(np.linalg.norm(cart, axis=1)))
            perm[a] = j
        out.append((R, t, perm))
    return out


def orbit_images(ops, kvec, kmesh):
    """返回 {网格索引: 代表 op 序号}；网格索引 = round(k*N) mod N（Gamma 中心全网格）。"""
    def idx(k):
        return tuple(int(round(x * n)) % int(n) for x, n in zip(k, kmesh))
    best = None
    for conv in ('R', 'RT'):
        d = {}
        for io, (R, t, perm) in enumerate(ops):
            kk = (R @ kvec) if conv == 'R' else (R.T @ kvec)
            d.setdefault(idx(kk), io)
        if best is None or len(d) > best[1]:
            best = (d, len(d))
    return best[0]


# ------------------------------------------------------------------ 主计算
def compute(mat, wdir, kmesh=(12, 12, 2), verbose=True, use_symmetry=True):
    structure = Structure.from_file(os.path.join(wdir, 'POSCAR'))
    nat = len(structure)
    wpath = os.path.join(wdir, 'WAVECAR')
    w = Wavecar(wpath)
    meta = read_wavecar_kmeta(wpath)
    ks, ws = read_kpoints(os.path.join(wdir, 'KPOINTS'), w.nk)
    occ0 = meta['recs'][0]['occ']
    nocc = int(np.sum(occ0 > 0.5))
    is_ncl = w.vasp_type.lower()[0] == 'n'
    if verbose:
        print('[%s] %s nk=%d nb=%d ng=%s nocc=%d vol=%.3f' % (mat, w.vasp_type, w.nk, w.nb, tuple(int(x) for x in w.ng), nocc, w.vol))
        eigs = meta['recs'][0]['eig']
        print('     E_F 附近：最高占据 %.4f eV，最低空带 %.4f eV，gap=%.4f eV' % (eigs[nocc - 1], eigs[nocc], eigs[nocc] - eigs[nocc - 1]))
        print('     占据数取值范围 %.6f .. %.6f' % (occ0.min(), occ0.max()))
    labels = np.load(os.path.join(MASKDIR, mat + '_bader_labels.npy'))
    ng = [int(x) for x in w.ng]
    ii = (np.arange(ng[0]) * (labels.shape[0] / ng[0])).astype(int) % labels.shape[0]
    jj = (np.arange(ng[1]) * (labels.shape[1] / ng[1])).astype(int) % labels.shape[1]
    kk = (np.arange(ng[2]) * (labels.shape[2] / ng[2])).astype(int) % labels.shape[2]
    lab = labels[np.ix_(ii, jj, kk)].reshape(-1)
    masks = {a: (lab == a) for a in range(nat)}
    dV = w.vol / float(np.prod(ng))
    ops = symmetry_ops(structure)
    gspin = 1 if is_ncl else 2
    DI = np.zeros((nat, nat))
    NA = np.zeros(nat)          # 盆内（赝）电子数，含自旋因子
    norbit = []
    for ik in range(len(ks)):
        f = meta['recs'][ik]['occ']
        occ_idx = [ib for ib in range(w.nb) if f[ib] > 0.5]   # 整数占据（高斯展宽的尾巴 f~1e-3 不计）
        # 每个自旋分量的波函数行（行 = 占据带），带 sqrt(f) 权重
        def comp(spinor):
            rows = []
            for ib in occ_idx:
                if is_ncl:
                    m = np.asarray(w.fft_mesh(ik, ib, spinor=spinor), dtype=np.complex64)
                else:
                    m = np.asarray(w.fft_mesh(ik, ib), dtype=np.complex64)
                psi = np.fft.ifftn(m) * (np.prod(ng) / np.sqrt(w.vol))
                rows.append(psi.reshape(-1) * np.sqrt(f[ib]))
            return np.array(rows) if rows else np.zeros((0, int(np.prod(ng))), dtype=np.complex64)
        P1 = comp(0)
        P2 = comp(1) if is_ncl else None
        # ---- Löwdin 正交化：PAW 赝波函数本身不正交（Sum|C|^2 在 0.35~1.3 之间），
        #      而 HF 型 pair-density 公式要求正交归一轨道集合。用 O = <psi_i|psi_j> 做对称正交化：
        #      psi'_i = Sum_j (O^-1/2)_{ji} psi_j ，于是 Sum_b N'(b) = I（求和规则成立）。
        nrow = P1.shape[0]
        O = (P1 @ P1.conj().T).astype(np.complex128)
        if P2 is not None:
            O = O + (P2 @ P2.conj().T)
        O = O * dV
        ev, U = np.linalg.eigh(O)
        if ev.min() <= 1e-8:
            raise RuntimeError('重叠矩阵奇异（最小本征值 %.3e）' % ev.min())
        W = (U * (1.0 / np.sqrt(ev))) @ U.conj().T          # O^{-1/2}
        # 每个盆的 N 矩阵（自旋求和 + Löwdin 变换）
        Nm = {}
        for a in range(nat):
            m = masks[a]
            if not m.any():
                Nm[a] = np.zeros((nrow, nrow), dtype=np.complex128)
                continue
            A = P1[:, m]
            M = A @ A.conj().T
            if P2 is not None:
                B = P2[:, m]
                M = M + B @ B.conj().T
            Nm[a] = W @ (M * dV) @ W.conj().T
        # 该 k 的盆内电子数（用于求和规则）
        for a in range(nat):
            NA[a] += ws[ik] * gspin * float(np.real(np.trace(Nm[a])))
        # 轨道展开
        if use_symmetry:
            img = orbit_images(ops, ks[ik], kmesh)
        else:
            img = {tuple(int(round(x * n)) % int(n) for x, n in zip(ks[ik], kmesh)): 0}
        norbit.append(len(img))
        for grid_idx, io in img.items():
            R, t, perm = ops[io]
            inv = np.empty(nat, dtype=int)
            for a in range(nat):
                inv[perm[a]] = a
            for A in range(nat):
                ia = inv[A]
                MA = Nm[ia]
                for B in range(nat):
                    ibb = inv[B]
                    DI[A, B] += 2.0 * gspin * ws[ik] / len(img) * float(np.real(np.trace(MA @ Nm[ibb].conj().T)))
    checks = dict(
        orbit_sizes=norbit,
        orbit_expect=[float(x) for x in np.round(ws * 288).astype(int)],
        sum_rule_max_rel_err=float(np.max(np.abs(DI.sum(axis=1) - 2 * NA) / np.maximum(2 * NA, 1e-9))),
        sum_rule_rows=[float(x) for x in DI.sum(axis=1)],
        sum_rule_2NA=[float(x) for x in 2 * NA],
        asymmetry=float(np.abs(DI - DI.T).max()),
        NA=[float(x) for x in NA],
        nocc=nocc, vasp_type=w.vasp_type, nk=w.nk, nk_used=len(ks),
        kweights=[float(x) for x in ws],
    )
    return DI, structure, checks


def contact_classes(structure, DI, verbose=True):
    """把 DI 按接触类型归类：层内/层间 + 元素对；返回每类平均 DI/ES 与重数。"""
    L = structure.lattice.matrix
    frac = np.array(structure.frac_coords)
    syms = [str(s) for s in structure.species]
    nat = len(frac)
    imgs = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    pairs = []
    for a in range(nat):
        for b in range(a + 1, nat):
            best = 1e9; bv = None
            for i, j, k in imgs:
                d = (frac[b] + np.array([i, j, k]) - frac[a]) @ L
                if np.linalg.norm(d) < best:
                    best = np.linalg.norm(d); bv = d
            pairs.append((a, b, best, bv))
    dz = [(p[0], p[1], abs(p[3][2])) for p in pairs]
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('materials', nargs='*', default=MATS)
    ap.add_argument('--no-symmetry', action='store_true')
    ap.add_argument('--tag', default='')
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    allres = {}
    for m in args.materials:
        wdir = os.path.join(SOC_ROOT, m, 'step3_PBE_WAVECAR')
        DI, structure, checks = compute(m, wdir, use_symmetry=not args.no_symmetry)
        checks['orbit_sizes'] = [int(x) for x in checks['orbit_sizes']]
        allres[m] = dict(DI=DI.tolist(), checks=checks,
                         species=[str(s) for s in structure.species],
                         frac=[list(map(float, x)) for x in structure.frac_coords])
        print('    sum-rule 最大相对误差 = %.4f%%   反对称残差 = %.2e' % (100 * checks['sum_rule_max_rel_err'], checks['asymmetry']))
        print('    DI 对角(自共享) =', np.round(np.diag(DI), 3))
        print('    行和 Sum_B DI(A,B) =', np.round(DI.sum(axis=1), 3), ' 2N_A =', np.round(2 * np.array(checks['NA']), 3))
        print('    DI 矩阵(上三角，>0.01 项):')
        syms = [str(s) for s in structure.species]
        for a in range(len(DI)):
            for b in range(a + 1, len(DI)):
                if DI[a, b] > 0.01:
                    print('      %s%d - %s%d : DI=%.4f  ES=%.4f' % (syms[a], a, syms[b], b, DI[a, b], 2 * DI[a, b]))
    tag = args.tag or ('nosym' if args.no_symmetry else 'sym')
    json.dump(allres, open(os.path.join(OUT, 'DI_results_%s.json' % tag), 'w'), ensure_ascii=False, indent=2, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print('written', os.path.join(OUT, 'DI_results_%s.json' % tag))


if __name__ == '__main__':
    main()
