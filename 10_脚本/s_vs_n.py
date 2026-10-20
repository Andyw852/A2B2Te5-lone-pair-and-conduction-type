"""In-plane S against net carrier concentration down to the intrinsic region, plus the numbers that describe that region.
Run: python s_vs_n.py DELIVERY      (no new DFT; dense BoltzTraP2 CRTA + dense EIGENVAL)
S is interpolated in the chemical potential, not in log|n|: near the neutral point |n| drops below the smallest sampled
value of the 8 meV mu grid, and an interpolation in log|n| then clamps to the last point.
Concentrations are net values (p - n); 'electron' rows are net electron concentrations.
"""
from pathlib import Path
import sys,csv
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm

LEVELS=[1e17,3e17,1e18,3e18,1e19,3e19,1e20]
TEMPS=[300,500]

def curve_from(mu,N,S,full=False):
 """S as a function of net concentration n (holes > 0), interpolated in mu within 0.35 eV of the neutral point."""
 i0=int(np.where(np.diff(np.sign(N))!=0)[0][0]);w=np.abs(mu-mu[i0])<0.35;mw,Nw,Sw=mu[w],N[w],S[w];o=np.argsort(Nw)
 assert np.all(np.diff(Nw[o])>0) and np.all(np.diff(mw[o])<0),'n(mu) not monotonic in the window'
 f=lambda n:float(np.interp(np.interp(n,Nw[o],mw[o]),mw,Sw))
 return (f,(Nw.min(),Nw.max())) if full else f

def curves(D,m,T):
 j=sm.bt2(D,m);ti=j['temperatures_K'].index(T);mu=np.array(j['mu_rel_efermi_eV'])
 N=np.array(j['carrier_conc_cm-3'][ti]);S=(np.array(j['seebeck_xx_V_per_K'][ti])+np.array(j['seebeck_yy_V_per_K'][ti]))/2*1e6
 return curve_from(mu,N,S,full=True)

def intrinsic(D,m,T):
 u=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform';nel,K,E=sm.read_eigenval(u/'EIGENVAL')
 j=sm.bt2(D,m);sc=j['scissor'];Es=E.copy();Es[:,nel:]+=sc['shift_eV'];vbm=E[:,:nel].max()
 z=brentq(lambda z:(1-sm.fermi(Es[:,:nel],z,T)).sum()-sm.fermi(Es[:,nel:],z,T).sum(),vbm-0.5,vbm+0.8)
 return sm.fermi(Es[:,nel:],z,T).sum()/len(K)/j['cell_volume_cm3'],sc['target_gap_eV']

def main(D):
 O=D/sm.OUT;rows,feat=[],[]
 for T in TEMPS:
  for m in sm.M:
   f,(lo,hi)=curves(D,m,T);ni,gap=intrinsic(D,m,T)
   for car,sg in (('hole',1),('electron',-1)):
    for n in LEVELS:
     if lo<=sg*n<=hi:rows.append({'material':m,'T_K':T,'carrier':car,'net_n_cm3':n,'S_in_plane_uV_K':round(f(sg*n),1)})
   ts=np.logspace(16,20.3,600)
   se=np.array([f(-t) for t in ts]);sh=np.array([f(t) for t in ts]);k=int(np.argmax(se<0)) if (se<0).any() else None
   feat.append({'material':m,'T_K':T,'gap_eV':gap,'n_i_cm3':float('%.3g'%ni),'S_at_net_zero_uV_K':round(f(0.),1),
    'electron_S_negative_above_cm3':float('%.2g'%ts[k]) if k is not None else None,
    'hole_S_max_uV_K':round(float(sh.max()),1),'hole_S_max_at_cm3':float('%.2g'%ts[sh.argmax()]),
    'electron_S_min_uV_K':round(float(se.min()),1),'electron_S_min_at_cm3':float('%.2g'%ts[se.argmin()])})
 for name,data in (('S随净浓度_300K_500K.csv',rows),('低浓度特征量.csv',feat)):
  with open(O/name,'w',newline='',encoding='utf-8-sig') as fh:
   w=csv.DictWriter(fh,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 for r in feat:print(r)

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
