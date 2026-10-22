"""Controlled distortion (15_): band gap, Gamma levels, valleys and BoltzTraP2 S at fixed n for each lambda structure.
Run: python lambda_bands.py DELIVERY MATERIAL LAMBDA      (LAMBDA 0.00 / 0.50 / 1.00; about 5 min and 10 GB for lambda < 1)
lambda < 1: 15_受控畸变/<m>/lam<LAMBDA>/EIGENVAL.gz with POSCAR_lam<LAMBDA>; lambda = 1: the baseline step3_uniform EIGENVAL,
whose S, Lorenz and n_i are taken from the baseline rows of 12_/剪刀敏感性_300K.csv (same code and settings).
Scissor: the fixed shift of lambda = 1 (target - PBE gap, 12_ section 3.1), so a gap change at lambda < 1 is kept.
"""
from pathlib import Path
import sys,csv,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm
import s_vs_n as sv
import bt2_local as bl

OUT='15_受控畸变/λ扫描_能带与S.csv'

def main(D,m,lam):
 base=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform'
 eig,pos=(base/'EIGENVAL',base/'POSCAR') if lam=='1.00' else (D/'15_受控畸变'/m/f'lam{lam}'/'EIGENVAL.gz',D/'15_受控畸变'/m/f'POSCAR_lam{lam}')
 nel,K,E=sm.read_eigenval(eig);mesh=np.array(json.loads((base/'dense_grid_manifest.json').read_text())['mesh']);idx=sm.grid(K,mesh)
 vbm,cbm=E[:,:nel].max(),E[:,nel:].min();ig=int(np.argmin(np.linalg.norm(np.mod(K+0.5,1)-0.5,axis=1)))
 j=sm.bt2(D,m);shift=j['scissor']['shift_eV'];gap=cbm-vbm+shift
 r={'material':m,'lambda':lam,'PBE_gap_eV':round(float(cbm-vbm),4),'scissored_gap_eV':round(float(gap),4),
    'VBM_k':str([round(float(x),3) for x in K[E[:,:nel].max(1).argmax()]]),'CBM_k':str([round(float(x),3) for x in K[E[:,nel:].min(1).argmin()]]),
    'Gamma_VB_below_VBM_meV':round(float(1000*(vbm-E[ig,nel-1])),1),'Gamma_CB_above_CBM_meV':round(float(1000*(E[ig,nel]-cbm)),1),'Gamma_direct_gap_PBE_eV':round(float(E[ig,nel]-E[ig,nel-1]),4)}
 for edge,bands,ref in (('VB',list(range(nel-6,nel)),vbm),('CB',list(range(nel,nel+6)),cbm)):
  st,neff=sm.stars(sm.pockets(idx,E,mesh,bands,'VBM' if edge=='VB' else 'CBM',ref,0.1),ref)
  r[f'{edge}_Nv_eff']=round(float(neff),2);r[f'{edge}_stars']='; '.join(f"{x['mult']}@{x['dE_meV']}meV" for x in st[:4])
 if lam=='1.00':
  b=next(x for x in csv.DictReader(open(D/sm.OUT/'剪刀敏感性_300K.csv',encoding='utf-8-sig')) if x['material']==m and x['is_baseline']=='True')
  for k in ('n_i_300K_cm3','S_net_zero','electron_S_negative_above_cm3','hole_S_max','hole_S_max_at_cm3','S_p_1e19','S_p_3e19','S_p_1e20','S_n_1e19','S_n_3e19','S_n_1e20','L_p_1e19_1e-8','L_n_1e19_1e-8'):r[k]=b[k]
 else:
  nw,v,c,eb,vv,vuc=bl.fit(D,m,eig=eig,pos=pos);t=bl.transport(nw,v,c,eb,vv,vuc,gap)
  f=sv.curve_from(t['mu_eV'],t['n'][0],t['S'][0]);ts=np.logspace(16,20.3,600)
  se=np.array([f(-x) for x in ts]);sh=np.array([f(x) for x in ts]);k=int(np.argmax(se<0)) if (se<0).any() else None
  N=t['n'][0];o=np.argsort(N);ni=float(np.sqrt(max(np.interp(0,N[o],t['p_only'][0][o])*np.interp(0,N[o],t['n_only'][0][o]),0)))
  r.update({'n_i_300K_cm3':float('%.3g'%ni),'S_net_zero':round(f(0.),1),'electron_S_negative_above_cm3':float('%.2g'%ts[k]) if k is not None else None,
            'hole_S_max':round(float(sh.max()),1),'hole_S_max_at_cm3':float('%.2g'%ts[sh.argmax()])})
  for car,sg in (('p',1),('n',-1)):
   for lv in (1e19,3e19,1e20):r[f'S_{car}_{lv:.0e}'.replace('+','')]=round(bl.at(t,0,sg*lv,'S'),1)
   r[f'L_{car}_1e19_1e-8']=round(1e8*bl.at(t,0,sg*1e19,'kap')/(bl.at(t,0,sg*1e19,'sig')*300),3)
 p=D/OUT;old=[x for x in csv.DictReader(open(p,encoding='utf-8-sig')) if (x['material'],x['lambda'])!=(m,lam)] if p.exists() else []
 rows=sorted(old+[r],key=lambda x:(sm.M.index(x['material']),x['lambda']));keys=list(r)
 with open(p,'w',newline='',encoding='utf-8-sig') as fh:
  w=csv.DictWriter(fh,fieldnames=keys);w.writeheader();w.writerows([{k:x.get(k) for k in keys} for x in rows])
 print(r)

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),sys.argv[2],sys.argv[3])
