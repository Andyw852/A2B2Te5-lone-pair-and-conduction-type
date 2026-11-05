#!/usr/bin/env python3
"""SOC 与泛函对孤对电子（阳离子 s）分布的影响（2026-11-05）：同一条 c2 路径上，PBE 不加 SOC（SOC机制检验/d1_nosoc_<材料>）、
PBE+SOC（c2_pbe）、HSE06+SOC（c2_hse）三种计算的阳离子 s 成分。COHP 做不了 SOC（COGITO 要求 LSORBIT = False，LOBSTER 5.1.1
读不了旋量波函数），用它说明 SOC 不改变“价带顶附近是填满的阳离子 s–Te p 反键态”这一结论。
每个路径 k 点上价带顶以下 1 eV、3 eV 内的阳离子 s 电子数（不加 SOC 时每条带 2 个电子），以及 Σ 谷（路径第 11 点）最高价带的阳离子 s 成分。
输出：22_能带基准_1-0/SOC机制检验/阳离子s成分_PBE_PBE+SOC_HSE+SOC.csv
用法：python 10_脚本/soc_s_weight.py
"""
import gzip, re, numpy as np
import os, csv
B=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '22_能带基准_1-0')
M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']; NEL={m:88 if 'Bi' in m else 68 for m in M}
base=lambda m: B if m=='Pb2Sb2Te5' else f'{B}/{m}'
def procar(p):
    E=[];W=[]; k=-1; blk=0; rows=[]
    for l in gzip.open(p,'rt'):
        if l.startswith(' k-point'): k+=1; E.append([]); W.append([])
        elif l.startswith('band'): E[k].append(float(l.split()[4])); blk=0; rows=[]
        elif l.startswith('ion'): blk=1; rows=[]
        elif blk==1 and re.match(r'\s+\d+\s',l): rows.append([float(x) for x in l.split()[1:10]])
        elif l.startswith('tot') and blk==1: W[k].append(np.array(rows)); blk=2
    return np.array(E), np.array(W)   # E[k,b], W[k,b,ion,orb]
rows=[]
print(f"{'material':10s} {'calc':9s} {'A-s':>6s} {'B-s':>6s} {'cat-s':>6s}   {'A-s':>6s} {'B-s':>6s}  | Σ valley top VB: cat-s %")
for m in M:
    nel=NEL[m]
    for lab,p,soc in (('PBE',f'{B}/SOC机制检验/d1_nosoc_{m}/PROCAR_OPT.gz',False),('PBE+SOC',f'{base(m)}/c2_pbe/PROCAR_OPT.gz',True),('HSE+SOC',f'{base(m)}/c2_hse/PROCAR_OPT.gz',True)):
        E,W=procar(p); nv=nel if soc else nel//2; occ=1 if soc else 2
        vbm=E[:,nv-1].max(); t=np.maximum(W.sum((-1,-2)),1e-9)
        As=W[:,:,7:9,0].sum(-1)/t; Bs=W[:,:,0:2,0].sum(-1)/t
        out=[]
        for lo in (1.0,3.0):
            sel=(E<=vbm+1e-4)&(E>=vbm-lo); sel[:,nv:]=False
            out.append((occ*(As*sel).sum()/len(E), occ*(Bs*sel).sum()/len(E)))
        sv=100*((As+Bs)[11,nv-1])
        rows.append([m, lab, round(out[0][0],3), round(out[0][1],3), round(out[1][0],3), round(out[1][1],3), round(sv,1)])
        print(f"{m:10s} {lab:9s} {out[0][0]:6.3f} {out[0][1]:6.3f} {out[0][0]+out[0][1]:6.3f}   {out[1][0]:6.3f} {out[1][1]:6.3f}  | {sv:5.1f}")
print('前三列：价带顶以下 1 eV 内、每个路径 k 点的阳离子 s 电子数（A 位、B 位、合计）；后两列：3 eV 内')
with open(os.path.join(B, 'SOC机制检验', '阳离子s成分_PBE_PBE+SOC_HSE+SOC.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w=csv.writer(f); w.writerow(['material','calc','A_s_electrons_per_k_top1eV','B_s_electrons_per_k_top1eV','A_s_electrons_per_k_top3eV','B_s_electrons_per_k_top3eV','Sigma_valley_top_VB_cation_s_percent']); w.writerows(rows)
