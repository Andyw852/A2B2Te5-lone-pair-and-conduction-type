"""22_：同一 k 网格上逐点比较 PBE+SOC 与 HSE+SOC 的价带顶、导带底附近能量。

用法：python 10_脚本/kpoint_compare_pbe_hse.py 22_能带基准_1-0
读 <目录>/b_PBE_SOC/EIGENVAL、c1_HSE_SOC/EIGENVAL、b_PBE_SOC/IBZKPT（两者 k 点相同）；SOC 下第 NELECT 条为价带、NELECT+1 条为导带。
"""
import sys, numpy as np
def read_eig(p):
    L = open(p).read().split('\n')
    nel, nk, nb = map(int, L[5].split())
    ks, E, i = [], [], 7
    for _ in range(nk):
        while not L[i].strip(): i += 1
        ks.append([float(x) for x in L[i].split()[:4]]); i += 1
        E.append([float(L[i + j].split()[1]) for j in range(nb)]); i += nb
    return nel, np.array(ks), np.array(E)
def read_mult(p):
    L = open(p).read().split('\n'); n = int(L[1])
    return np.array([int(L[3 + j].split()[3]) for j in range(n)])
d = sys.argv[1]
nel, kb, Eb = read_eig(d + '/b_PBE_SOC/EIGENVAL')
_, kh, Eh = read_eig(d + '/c1_HSE_SOC/EIGENVAL')
mult = read_mult(d + '/b_PBE_SOC/IBZKPT')
assert np.allclose(kb[:, :3], kh[:, :3])
v, c = nel - 1, nel  # 0-based: band 68, 69
out = []
for tag, E in (('PBE', Eb), ('HSE', Eh)):
    vbm, cbm = E[:, v].max(), E[:, c].min()
    out.append((E[:, v] - vbm, E[:, c] - cbm, vbm, cbm))
    print(f'{tag}: VBM {vbm:.4f} at k#{E[:, v].argmax()+1}  CBM {cbm:.4f} at k#{E[:, c].argmin()+1}  gap {cbm-vbm:.4f}')
names = {(0,0,0):'Γ',(0.5,0,0):'M',(0,0,0.5):'A',(0.5,0,0.5):'L'}
def name(k):
    for kk, n in names.items():
        if np.allclose(k, kk, atol=1e-6): return n
    return ''
print('\n k#  kx     ky     kz    mult | VB-VBM PBE  HSE (meV) | CB-CBM PBE  HSE (meV) | gap@k PBE HSE (meV)')
order = np.argsort(out[1][0])[::-1]
for j in order:
    print(f'{j+1:3d} {kb[j,0]:.3f} {kb[j,1]:.3f} {kb[j,2]:.3f} {mult[j]:3d} {name(kb[j,:3]):1s}| {1000*out[0][0][j]:8.1f} {1000*out[1][0][j]:8.1f} | {1000*out[0][1][j]:8.1f} {1000*out[1][1][j]:8.1f} | {1000*(Eb[j,c]-Eb[j,v]):7.1f} {1000*(Eh[j,c]-Eh[j,v]):7.1f}')
for w in (25, 50, 100):
    for t, o in zip(('PBE', 'HSE'), out):
        print(f'{t}: 价带距 VBM ≤{w} meV 的 k 点（含对称等价）数 {mult[o[0] >= -w/1000].sum():3d}；导带 ≤{w} meV {mult[o[1] <= w/1000].sum():3d}')
