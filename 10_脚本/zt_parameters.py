"""Every ZT ingredient that existing outputs can give, per material, next to the bonding descriptors the lone-pair argument uses.
Run: python zt_parameters.py DELIVERY      (no new DFT)
Lattice side : Born effective charges and Gamma phonons (08 SOC DFPT OUTCAR.gz), elastic tensor (05 AMSET settings)
               -> VRH moduli, sound velocities, Debye temperature, Cahill kappa_min, Poisson-ratio Grueneisen, Slack kappa_L (empirical).
Electron side: in-plane S, m*_D, m*_S (12_), m*_c (14_), Lorenz number (BoltzTraP2); mobility, PF and weighted mobility only
               from the dense-grid AMSET in 14_散射机制分项迁移率/密网格AMSET (left empty until those runs are there).
"""
from pathlib import Path
import json,csv,sys,gzip,re
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parent))
import scattering_mobility as smob

M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
MASS={'Pb':207.2,'Sn':118.71,'Sb':121.76,'Bi':208.98,'Te':127.6}
OUT='18_ZT参数与孤对电子'
KB=1.380649e-23;H=6.62607015e-34;AMU=1.66053907e-27

def poscar(p):
 L=p.read_text().splitlines();lat=np.array([list(map(float,l.split()[:3])) for l in L[2:5]])*float(L[1].split()[0])
 return lat,[s for s,c in zip(L[5].split(),map(int,L[6].split())) for _ in range(c)]

def dfpt(p,n):
 """Last Born-charge block and the Gamma eigenvectors (mass-weighted) from a VASP IBRION=7/8 OUTCAR."""
 t=gzip.open(p,'rt').read()
 b=t.split('BORN EFFECTIVE CHARGES')[-1].split('\n')
 Z=[];i=0
 while len(Z)<n:
  if b[i].strip().startswith('ion'):Z.append([list(map(float,b[i+k].split()[1:4])) for k in (1,2,3)]);i+=4
  else:i+=1
 blk=t.split('Eigenvectors and eigenvalues of the dynamical matrix')[-1]
 modes=[]
 for m in re.finditer(r'\n\s*\d+ (f|f/i)\s*=\s*([\d.]+) THz.*?\n.*?dz\s*\n((?:.*\n){%d})'%n,blk):
  f=float(m.group(2))*(-1 if m.group(1)=='f/i' else 1);e=np.array([list(map(float,l.split()[3:6])) for l in m.group(3).splitlines()])
  modes.append((f,e))
 return np.array(Z),modes

def elastic(C,rho,lat,sym):
 S=np.linalg.inv(C)
 KV=(C[0,0]+C[1,1]+C[2,2]+2*(C[0,1]+C[0,2]+C[1,2]))/9;GV=(C[0,0]+C[1,1]+C[2,2]-(C[0,1]+C[0,2]+C[1,2])+3*(C[3,3]+C[4,4]+C[5,5]))/15
 KR=1/(S[0,0]+S[1,1]+S[2,2]+2*(S[0,1]+S[0,2]+S[1,2]));GR=15/(4*(S[0,0]+S[1,1]+S[2,2])-4*(S[0,1]+S[0,2]+S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))
 K,G=(KV+KR)/2*1e9,(GV+GR)/2*1e9
 vl,vt=np.sqrt((K+4*G/3)/rho),np.sqrt(G/rho);vm=(((2/vt**3)+(1/vl**3))/3)**(-1/3)
 n=len(sym);V=abs(np.linalg.det(lat))*1e-30;N=n/V
 thD=H/KB*(3*N/(4*np.pi))**(1/3)*vm;nu=(3*K-2*G)/(2*(3*K+G));gam=3*(1+nu)/(2*(2-3*nu))
 kmin=(np.pi/6)**(1/3)*KB*N**(2/3)*(vl+2*vt)
 Mbar=np.mean([MASS[s] for s in sym]);delta=(abs(np.linalg.det(lat))/n)**(1/3);tha=thD*n**(-1/3)
 A=2.43e-6/(1-0.514/gam+0.228/gam**2)
 slack=A*Mbar*tha**3*delta/(gam**2*n**(2/3)*300)
 return {'K_GPa':K/1e9,'G_GPa':G/1e9,'poisson':nu,'v_l_m_s':vl,'v_t_m_s':vt,'v_m_m_s':vm,'theta_D_K':thD,'gamma_poisson':gam,
         'kappa_min_Cahill_W_mK':kmin,'kappa_L_Slack_300K_W_mK':slack}

def interp_n(N,val,n):
 sel=np.sign(N)==np.sign(n);o=np.argsort(np.abs(N[sel]));nn,vv=np.abs(N[sel])[o],val[sel][o];keep=np.concatenate([[True],np.diff(nn)>0])
 return float(np.interp(np.log(abs(n)),np.log(nn[keep]),vv[keep]))

def main(D):
 O=D/OUT;O.mkdir(exist_ok=True);lat_rows,el_rows=[],[]
 mass_tab={(r['material'],r['carrier'],float(r['n_cm3'])):r for r in csv.DictReader(open(D/'12_s占比与Seebeck有效质量/有效质量_Seebeck.csv',encoding='utf-8-sig')) if float(r['T_K'])==300}
 mp=D/'14_散射机制分项迁移率/分项迁移率_全温度全浓度.csv'
 mob={(r['material'],r['carrier'],abs(float(r['doping_cm3']))):r for r in csv.DictReader(open(mp,encoding='utf-8-sig')) if float(r['T_K'])==300} if mp.exists() else {}
 par={r['material']:r for r in csv.DictReader(open(D/'14_散射机制分项迁移率/散射输入参数.csv',encoding='utf-8-sig'))}
 for m in M:
  src=D/'08_介电复核/源数据'/m;lat,sym=poscar(src/'POSCAR');n=len(sym)
  Z,modes=dfpt(src/'OUTCAR.gz',n);A='Pb' if 'Pb' in sym else 'Sn';B='Sb' if 'Sb' in sym else 'Bi'
  rho=sum(MASS[s] for s in sym)*AMU/(abs(np.linalg.det(lat))*1e-30)
  s=yaml.safe_load((D/'05_SOC_AMSET补算'/m/'04_AMSET输运/settings.yaml').read_text());C=np.array(s['elastic_constant'],float)
  r={'material':m,'density_g_cm3':rho/1000}
  for el,nom in ((A,2),(B,3),('Te',-2)):
   zz=Z[[i for i,x in enumerate(sym) if x==el]]
   r[f'Zstar_{el if el=="Te" else ("A" if el==A else "B")}_xx']=float(zz[:,0,0].mean());r[f'Zstar_{el if el=="Te" else ("A" if el==A else "B")}_zz']=float(zz[:,2,2].mean())
  r['Zstar_A_nominal'],r['Zstar_B_nominal']=2,3
  opt=[(f,e) for f,e in modes if abs(f)>0.1]
  r['lowest_optical_THz']=min(f for f,_ in opt);r['mean_optical_THz']=float(np.mean([f for f,_ in opt]))
  # Force constants rebuilt from the Gamma eigenpairs: Phi = M^1/2 (sum w^2 e e^T) M^1/2, in eV/A^2.
  w=np.array([2*np.pi*f*1e12*abs(f)/f for f,_ in modes]);E=np.array([e.reshape(-1) for _,e in modes]);E/=np.linalg.norm(E,axis=1)[:,None]
  mm=np.repeat([MASS[x]*AMU for x in sym],3);Phi=np.sqrt(np.outer(mm,mm))*((E.T*np.sign(w)*w**2)@E)/16.0218
  # Re-centring displacement exactly as in the 15_ frozen scan (POSCARs identical to this DFPT cell), split by site.
  fr=lambda name:np.array([list(map(float,l.split()[:3])) for l in (D/'15_受控畸变'/m/name).read_text().splitlines()[8:8+n]])
  u0=(((fr('POSCAR_lam0.00')-fr('POSCAR_lam1.00'))+0.5)%1-0.5)@lat;site=np.array(sym)[:,None]
  for tag,u in (('B',u0*(site==B)),('A',u0*(site==A))):
   v=u.reshape(-1);r[f'{tag}_offcentre_A']=float(np.abs(u[:,2]).max())
   r[f'{tag}_offcentre_k_eV_A2']=float(v@Phi@v/(v@v));r[f'{tag}_harmonic_recentre_meV_per_cell']=1000*0.5*float(v@Phi@v)
  r.update(elastic(C,rho,lat,sym));r.update({k:float(par[m][k]) for k in ('eps_inf_xx','eps0_xx','omega_LO_THz')})
  lat_rows.append(r)
  j=json.loads((D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step8.1_boltztrap/boltztrap_crta.json').read_text());ti=j['temperatures_K'].index(300)
  N=np.array(j['carrier_conc_cm-3'][ti]);Lz=np.array(j['lorenz_WOhm_per_K2'][ti])
  t=smob.dense(D,m)
  for car,sg in (('hole',1),('electron',-1)):
   for lv in (1e19,1e20):
    x=mass_tab.get((m,car,lv));y=mob.get((m,'p' if sg>0 else 'n',lv))
    pf=None
    if t and 300 in t['temperatures']:
     dop=np.array(t['doping']);tt=t['temperatures'].index(300);di=int(np.argmin(np.abs(dop-sg*lv)))
     if np.isclose(dop[di],sg*lv,rtol=1e-3):
      sig=np.array(t['conductivity'])[di,tt];see=np.array(t['seebeck'])[di,tt];pf=np.mean([sig[a,a]*(see[a,a]*1e-6)**2 for a in (0,1)])*1e3
    mu=float(y['mu_overall_inplane']) if y else None;mD=float(x['m_DOS']) if x else None
    el_rows.append({'material':m,'carrier':car,'n_cm3':lv,'S_in_plane_uV_K':float(x['S_in_plane_uV_K']) if x else None,'m_DOS':mD,
     'm_S_CRTA_in_plane':float(x['m_S_CRTA_in_plane']) if x and x['m_S_CRTA_in_plane'] else None,
     'm_c_in_plane':float(par[m][f"m_c_inplane_{car}_{lv:.0e}".replace('+','')]),'mu_AMSET_in_plane':mu,
     'mu_weighted':mu*mD**1.5 if mu and mD else None,'Lorenz_1e-8_WOhm_K2':interp_n(N,Lz,sg*lv)*1e8,'PF_AMSET_in_plane_mW_mK2':pf,
     'audit_grade':x['audit_grade'] if x else ''})
 for name,rows in (('晶格与成键参数.csv',lat_rows),('电子参数_300K.csv',el_rows)):
  with open(O/name,'w',newline='',encoding='utf-8-sig') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 for r in lat_rows:print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()})
 for r in el_rows:print({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()})

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
