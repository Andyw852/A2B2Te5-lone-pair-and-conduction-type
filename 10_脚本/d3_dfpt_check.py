"""What the missing D3 term does to the 08 DFPT phonon quantities (eps_ion, eps0, pop_frequency, B-site stiffness).
Run: python d3_dfpt_check.py DELIVERY      (needs `pip install dftd3`; no new DFT, ~10 s)
D3 is a function of atomic positions only, so it leaves eps_inf and Z* at fixed geometry exactly unchanged; it only adds
a Hessian. Here that Hessian is taken by finite differences of D3(BJ) gradients (VASP IVDW=12 parameters, no ATM term,
VASP cutoffs 50.2 / 20 A) and added to the force constants rebuilt from the 08 Gamma eigenpairs (PBE only, D3 geometry).
pop_frequency follows amset 0.5.1 tools/phonon_frequency.py (no-NAC path): weight_m ~ |sum_k Z*_k e_mk / sqrt(M_k f_m)|.
"""
from pathlib import Path
import sys,json,csv
import numpy as np
from dftd3.interface import RationalDampingParam,DispersionModel
sys.path.insert(0,str(Path(__file__).resolve().parent))
import zt_parameters as zp

Z_NUM={'Pb':82,'Sn':50,'Sb':51,'Bi':83,'Te':52}
BOHR=0.529177210903;HA=27.211386245988;EV_A2=16.02176634   # eV/A^2 -> J/m^2
E_CH=1.602176634e-19;EPS0=8.8541878128e-12

def d3_hessian(lat,sym,frac,h=0.01):
 """D3(BJ)-PBE Hessian in eV/A^2 by central differences (h in bohr)."""
 pos=frac@lat/BOHR;cell=lat/BOHR
 mod=DispersionModel(np.array([Z_NUM[s] for s in sym]),pos,cell,periodic=np.array([True]*3))
 mod.set_realspace_cutoff(50.2/BOHR,50.2/BOHR,20.0/BOHR)
 par=RationalDampingParam(s6=1.0,a1=0.4289,s8=0.7875,a2=4.4407,s9=0.0)
 n=len(sym);H=np.zeros((3*n,3*n))
 for i in range(3*n):
  g=[]
  for s in (1,-1):
   p=pos.copy();p.reshape(-1)[i]+=s*h;mod.update(p,cell);g.append(mod.get_dispersion(par,grad=True)['gradient'].reshape(-1))
  H[i]=(g[0]-g[1])/(2*h)
 mod.update(pos,cell)
 return 0.5*(H+H.T)*HA/BOHR**2

def response(Phi,mm,Z,vol,eps_inf):
 """Modes, eps_ion, pop_frequency (amset weighting) and the highest in-plane LO from force constants Phi (eV/A^2)."""
 w2,V=np.linalg.eigh(Phi*EV_A2/np.sqrt(np.outer(mm,mm)))
 f=np.sign(w2)*np.sqrt(np.abs(w2))/(2*np.pi*1e12);opt=np.argsort(np.abs(f))[3:]
 n=len(mm)//3;p=np.array([sum(Z[k]@(V[3*k:3*k+3,m]/np.sqrt(mm[3*k])) for k in range(n)) for m in range(3*n)])
 eps_ion=E_CH**2/(EPS0*vol)*sum(np.outer(p[m],p[m])/w2[m] for m in opt)
 # amset: weight = mean over directions of |q.(sum_k Z*_k e_k/sqrt(M_k[amu] f))| ~ |sum_k Z*_k e_k/sqrt(M_k)| / sqrt(f)
 wt=np.array([np.linalg.norm(p[m]*np.sqrt(zp.AMU))/np.sqrt(f[m]) if f[m]>0 and m in opt else 0. for m in range(3*n)])
 pop=float((wt*f).sum()/wt.sum())
 q=np.array([1.,0,0]);qz=np.concatenate([q@Z[k] for k in range(n)])/np.sqrt(mm)
 lo=np.linalg.eigvalsh((V*w2)@V.T + E_CH**2/(EPS0*vol)*np.outer(qz,qz)/(q@eps_inf@q))
 return f[opt],eps_ion,pop,float(np.sqrt(lo.max())/(2*np.pi*1e12))

def main(D):
 diel=json.loads((D/'08_介电复核/介电复核.json').read_text());rows=[]
 for m in zp.M:
  src=D/'08_介电复核/源数据'/m;lat,sym=zp.poscar(src/'POSCAR');n=len(sym);vol=abs(np.linalg.det(lat))*1e-30
  frac=np.array([list(map(float,l.split()[:3])) for l in (src/'POSCAR').read_text().splitlines()[8:8+n]])
  Z,modes=zp.dfpt(src/'OUTCAR.gz',n);mm=np.repeat([zp.MASS[x]*zp.AMU for x in sym],3);ei=np.array(diel[m]['eps_inf'])
  w=np.array([2*np.pi*f*1e12 for f,_ in modes]);E=np.array([e.reshape(-1) for _,e in modes]);E/=np.linalg.norm(E,axis=1)[:,None]
  P0=np.sqrt(np.outer(mm,mm))*((E.T*np.sign(w)*w**2)@E)/EV_A2;P0=0.5*(P0+P0.T);Hd=d3_hessian(lat,sym,frac)
  fr=lambda name:np.array([list(map(float,l.split()[:3])) for l in (D/'15_受控畸变'/m/name).read_text().splitlines()[8:8+n]])
  B='Sb' if 'Sb' in sym else 'Bi';u=(((((fr('POSCAR_lam0.00')-fr('POSCAR_lam1.00'))+0.5)%1-0.5)@lat)*(np.array(sym)[:,None]==B)).reshape(-1)
  r={'material':m,'eps_inf_xx':ei[0,0],'eps_ion_xx_VASP':diel[m]['eps_ion'][0][0],'eps_ion_zz_VASP':diel[m]['eps_ion'][2][2],
     'pop_frequency_settings_THz':diel[m]['pop_frequency']}
  for tag,P in (('PBE',P0),('PBE+D3',P0+Hd)):
   f,ion,pop,lo=response(P,mm,Z,vol,ei);e0=ei[0,0]+ion[0,0]
   r.update({f'eps_ion_xx_{tag}':ion[0,0],f'eps_ion_zz_{tag}':ion[2,2],f'eps0_xx_{tag}':e0,f'pop_coupling_xx_{tag}':1/ei[0,0]-1/e0,
             f'pop_frequency_THz_{tag}':pop,f'top_LO_inplane_THz_{tag}':lo,f'lowest_optical_THz_{tag}':float(f.min()),
             f'B_stiffness_eV_A2_{tag}':float(u@P@u/(u@u)),f'B_harmonic_recentre_meV_{tag}':1000*0.5*float(u@P@u)})
  r['D3_hessian_max_eV_A2']=float(np.abs(Hd).max());r['PBE_force_const_max_eV_A2']=float(np.abs(P0).max())
  for k in ('eps0_xx','pop_coupling_xx','pop_frequency_THz','lowest_optical_THz','B_stiffness_eV_A2'):
   r[f'change_{k}_pct']=100*(r[f'{k}_PBE+D3']/r[f'{k}_PBE']-1)
  rows.append(r);print({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()})
 with open(D/zp.OUT/'D3对DFPT声子量的影响.csv','w',newline='',encoding='utf-8-sig') as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
