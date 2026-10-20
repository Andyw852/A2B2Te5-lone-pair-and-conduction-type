"""Controlled-distortion structures: move the B (Sb/Bi) and A (Pb/Sn) cations along c between their two Te planes.
Run: python controlled_distortion.py DELIVERY
lambda = 1 is the relaxed PBE-D3 cell, lambda = 0 puts each cation at its octahedron centre (3+3 bonds equal);
lattice and Te positions stay at the relaxed values (frozen scan), so only the off-centring changes.
Writes 15_受控畸变/<material>/POSCAR_lam<value> (0, 0.5, 1) and a bond-length table; no DFT is run.
The A site is nearly centred (0.007-0.018 A), so lambda = 0 is in practice a B-site re-centring.
"""
from pathlib import Path
import csv,sys,itertools
import numpy as np

M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
LAMBDAS=[0.0,0.5,1.0]
OUT='15_受控畸变'

def read(p):
 L=p.read_text().splitlines();lat=np.array([list(map(float,l.split()[:3])) for l in L[2:5]])*float(L[1].split()[0])
 sp,cn=L[5].split(),list(map(int,L[6].split()));sym=[s for s,c in zip(sp,cn) for _ in range(c)]
 F=np.array([list(map(float,l.split()[:3])) for l in L[8:8+len(sym)]])
 return lat,sp,cn,sym,F

def bonds(lat,sym,F,i):
 d=sorted(np.linalg.norm((F[j]-F[i]+np.array(o))@lat) for j,t in enumerate(sym) if t=='Te' for o in itertools.product((-1,0,1),repeat=3))
 return d[:3],d[3:6]

def octa_centre_z(lat,sym,F,i):
 """Mean z of the six nearest Te (unwrapped along c) = centre of the cation's octahedron."""
 c=[]
 for j,t in enumerate(sym):
  if t!='Te':continue
  for o in itertools.product((-1,0,1),repeat=3):
   r=F[j]+np.array(o);c.append((np.linalg.norm((r-F[i])@lat),r[2]))
 return float(np.mean([z for _,z in sorted(c)[:6]]))

def write(p,title,lat,sp,cn,F):
 p.write_text('\n'.join([title,'1.0',*(' '.join(f'{x:.12f}' for x in v) for v in lat),' '.join(sp),' '.join(map(str,cn)),'Direct',*(' '.join(f'{x:.12f}' for x in f) for f in np.mod(F,1))])+'\n')

def main(D):
 O=D/OUT;O.mkdir(exist_ok=True);rows=[]
 for m in M:
  lat,sp,cn,sym,F0=read(D/'A2B2Te5_SOC_dense/materials'/m/'POSCAR');(O/m).mkdir(exist_ok=True)
  A='Pb' if m.startswith('Pb') else 'Sn';B='Sb' if 'Sb' in m else 'Bi'
  zc={i:octa_centre_z(lat,sym,F0,i) for i,s in enumerate(sym) if s in (A,B)}
  for lam in LAMBDAS:
   F=F0.copy()
   for i,z in zc.items():F[i,2]=z+lam*(F0[i,2]-z)
   tag=f'lam{lam:.2f}'
   write(O/m/f'POSCAR_{tag}',f'{m} controlled distortion {tag} (lattice and Te frozen at relaxed PBE-D3)',lat,sp,cn,F)
   for el in (A,B):
    i=next(k for k,s in enumerate(sym) if s==el);s3,l3=bonds(lat,sym,F,i)
    rows.append({'material':m,'structure':tag,'cation':el,'short_A':round(s3[0],4),'long_A':round(l3[0],4),'delta_A':round(l3[0]-s3[0],4),
                 'offcentre_A':round((F[i,2]-zc[i])*lat[2,2],4)})
 with open(O/'键长与偏心.csv','w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 for r in rows:
  if r['structure'] in ('lam1.00','lam0.00'):print(r)

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
