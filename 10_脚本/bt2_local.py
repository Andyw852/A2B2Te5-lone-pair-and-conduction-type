"""Local BoltzTraP2 rerun from the dense EIGENVAL: band-gap (scissor) sensitivity and bipolar decomposition of kappa_e.
Run: python bt2_local.py DELIVERY [material ...]      (needs `pip install BoltzTraP2 ase`; about 1-3 min per material)
Same settings as the cluster step8.1 (all 120/144 bands, k-mesh multiplier 5, 4000 DOS bins, mu = mid-gap +-2 eV, 500 points,
CRTA). The fit is linear in the band energies, so a rigid conduction-band shift after the fit equals a scissor before it:
the bands are fitted once and every gap variant reuses them.
Bipolar term: valence and conduction bands are integrated separately at the common mu;
kappa_e = kappa_p + kappa_n + kappa_bp with kappa_bp = sigma_p sigma_n / (sigma_p + sigma_n) (S_p - S_n)^2 T (per direction).
"""
from pathlib import Path
import sys,csv,json
import numpy as np
import ase.io
from BoltzTraP2 import sphere,fite,bandlib as BL,units as U
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm
import s_vs_n as sv

HA=27.211386245988;BOHR_CM=0.529177210903e-8;TR=np.array([300.,500.])
LEVELS=[3e18,1e19,3e19,1e20]
DGAP=[-0.035,0.0,0.05,0.10]          # added to the HSE target gap; -0.035 keeps Pb2Sb2Te5 at 0.030 eV

class Data:
 def __init__(s,atoms,K,E):s.atoms,s.kpoints,s.ebands,s.mommat,s.magmom=atoms,K,E.T/HA,None,None
 def get_lattvec(s):return s.atoms.get_cell().T*U.Angstrom

def fit(D,m):
 u=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform';nel,K,E=sm.read_eigenval(u/'EIGENVAL')
 d=Data(ase.io.read(u/'POSCAR',format='vasp'),K,E);eq=sphere.get_equivalences(d.atoms,None,len(K)*5)
 eb,vv,_=fite.getBTPbands(eq,fite.fitde3D(d,eq),d.get_lattvec())
 return nel,E[:,:nel].max(),E[:,nel:].min(),eb,vv,abs(np.linalg.det(d.get_lattvec()))

def onsager(eb,vv,erange,mur,vuc):
 eps,dos,vvdos,_=BL.BTPDOS(eb,vv,erange=erange,npts=4000)
 N,L0,L1,L2,_=BL.fermiintegrals(eps,dos,vvdos,mur=mur,Tr=TR,dosweight=1.0)
 sig,see,kap,_=BL.calc_Onsager_coefficients(L0,L1,L2,mur,TR,vuc)
 return N,sig,see,kap

def transport(nel,vbm,cbm,eb,vv,vuc,gap):
 """All-band and per-carrier transport at a given scissored gap (eV); mu grid as on the cluster (mid-gap +-2 eV)."""
 eb=eb.copy();eb[nel:]+=(gap-(cbm-vbm))/HA;mid=(vbm+gap/2)/HA;mur=mid+np.linspace(-2,2,500)/HA
 er=(eb.min()-1e-4,eb.max()+1e-4);V=vuc*BOHR_CM**3
 N,sig,see,kap=onsager(eb,vv,er,mur,vuc);Nv,sv_,ssv,kv=onsager(eb[:nel],vv[:nel],er,mur,vuc);Nc,sc,ssc,kc=onsager(eb[nel:],vv[nel:],er,mur,vuc)
 ip=lambda t:(t[...,0,0]+t[...,1,1])/2
 sp,sn=sv_[...,[0,1],[0,1]],sc[...,[0,1],[0,1]];Sp,Sn=ssv[...,[0,1],[0,1]],ssc[...,[0,1],[0,1]]
 kbp=(sp*sn/np.where(sp+sn>0,sp+sn,np.nan)*(Sp-Sn)**2*TR[:,None,None]).mean(-1)
 return {'mu_eV':(mur-mid)*HA,'n':(N+nel)/V,'p_only':(Nv+nel)/V,'n_only':-Nc/V,'S':ip(see)*1e6,'sig':ip(sig),'kap':ip(kap),
         'sig_p':ip(sv_),'sig_n':ip(sc),'S_p':ip(ssv)*1e6,'S_n':ip(ssc)*1e6,'kap_p':ip(kv),'kap_n':ip(kc),'kap_bp':kbp}

def at(r,ti,n,key):
 """Value of r[key] at net concentration n (signed, holes > 0), interpolated in mu like s_vs_n."""
 mu,N=r['mu_eV'],r['n'][ti];i0=int(np.where(np.diff(np.sign(N))!=0)[0][0]);w=np.abs(mu-mu[i0])<0.35;o=np.argsort(N[w])
 return float(np.interp(np.interp(n,N[w][o],mu[w][o]),mu,r[key][ti]))

def main(D,mats):
 sens,bip,check=[],[],[]
 for m in mats:
  nel,vbm,cbm,eb,vv,vuc=fit(D,m);j=sm.bt2(D,m);target=j['scissor']['target_gap_eV']
  for dg in DGAP:
   g=target+dg;r=transport(nel,vbm,cbm,eb,vv,vuc,g)
   if dg==0:
    # Reproduction check against the cluster JSON on the same mu grid (300 K, 1e18 <= |n| <= 1e21).
    Sj=(np.array(j['seebeck_xx_V_per_K'][0])+np.array(j['seebeck_yy_V_per_K'][0]))/2*1e6;Nj=np.array(j['carrier_conc_cm-3'][0])
    sel=(np.abs(Nj)>=1e18)&(np.abs(Nj)<=1e21)
    check.append({'material':m,'max_dS_uV_K':round(float(np.abs(r['S'][0][sel]-Sj[sel]).max()),2),
                  'max_rel_dn':round(float((np.abs(r['n'][0][sel]-Nj[sel])/np.abs(Nj[sel])).max()),4)})
    for ti,T in enumerate(TR):
     for car,sg in (('hole',1),('electron',-1)):
      for lv in [3e18]+LEVELS[1:]:
       x={k:at(r,ti,sg*lv,k) for k in ('S','sig','kap','sig_p','sig_n','S_p','S_n','kap_p','kap_n','kap_bp')}
       bip.append({'material':m,'T_K':int(T),'carrier':car,'net_n_cm3':lv,'S_total_uV_K':round(x['S'],1),'S_p_uV_K':round(x['S_p'],1),'S_n_uV_K':round(x['S_n'],1),
        'minority_sigma_share':round((x['sig_n'] if sg>0 else x['sig_p'])/x['sig'],4),
        'L_total_1e-8':round(1e8*x['kap']/(x['sig']*T),3),'L_unipolar_1e-8':round(1e8*(x['kap_p']+x['kap_n'])/(x['sig']*T),3),
        'kappa_bp_share':round(x['kap_bp']/x['kap'],4),'sum_rule_err':round(abs(x['kap_p']+x['kap_n']+x['kap_bp']-x['kap'])/x['kap'],5)})
   f=sv.curve_from(r['mu_eV'],r['n'][0],r['S'][0]);ts=np.logspace(16,20.3,600)
   se=np.array([f(-t) for t in ts]);sh=np.array([f(t) for t in ts]);k=int(np.argmax(se<0)) if (se<0).any() else None
   ni=float(np.sqrt(np.clip(np.interp(0,np.sort(r['n'][0]),r['p_only'][0][np.argsort(r['n'][0])])*np.interp(0,np.sort(r['n'][0]),r['n_only'][0][np.argsort(r['n'][0])]),0,None)))
   row={'material':m,'gap_eV':round(g,4),'is_baseline':dg==0,'n_i_300K_cm3':float('%.3g'%ni),'S_net_zero':round(f(0.),1),
        'electron_S_negative_above_cm3':float('%.2g'%ts[k]) if k is not None else None,'hole_S_max':round(float(sh.max()),1),'hole_S_max_at_cm3':float('%.2g'%ts[sh.argmax()])}
   for car,sg in (('p',1),('n',-1)):
    for lv in LEVELS:row[f'S_{car}_{lv:.0e}'.replace('+','')]=round(at(r,0,sg*lv,'S'),1)
    for lv in (1e19,1e20):row[f'L_{car}_{lv:.0e}_1e-8'.replace('+','')]=round(1e8*at(r,0,sg*lv,'kap')/(at(r,0,sg*lv,'sig')*300),3)
   sens.append(row);print(row)
 for name,data,out in (('剪刀敏感性_300K.csv',sens,sm.OUT),('双极热导分解.csv',bip,'18_ZT参数与孤对电子'),('本地BoltzTraP2复现核对.csv',check,sm.OUT)):
  with open(D/out/name,'w',newline='',encoding='utf-8-sig') as fh:
   w=csv.DictWriter(fh,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 print(check)

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),sys.argv[2:] or sm.M)
