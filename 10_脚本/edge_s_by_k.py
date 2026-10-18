"""Per-k band-edge s character: separate 'which k the edge sits at' from 'which functional'.
extract (on the cluster, next to the 2 GB PROCAR):
  python edge_s_by_k.py extract STEP3_UNIFORM_DIR OUT.npz      # bands nel-4..nel+3, all k; ~2 MB
analyze (locally, after copying the four npz into 12_s占比与Seebeck有效质量/逐k投影/<material>.npz):
  python edge_s_by_k.py analyze DELIVERY
Per valley star (from seebeck_mass.py) the s fractions at the pocket k, and PBE s at the HSE path edge k vs the HSE value there.
"""
from pathlib import Path
import sys,json,csv
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm

def extract(u,out):
 u=Path(u);nel,K,E=sm.read_eigenval(u/'EIGENVAL');bands=list(range(nel-4,nel+4))
 sp=sum(([e]*int(c) for e,c in zip(*(l.split() for l in (u/'POSCAR').read_text().splitlines()[5:7]))),[])
 A=next(e for e in sp if e in ('Pb','Sn'));B=next(e for e in sp if e in ('Sb','Bi'))
 pro=next(p for p in (u/'PROCAR',u/'PROCAR.gz') if p.exists())
 Kp,W,T=sm.read_procar_s(pro,bands,[{i for i,e in enumerate(sp) if e==x} for x in (A,B,'Te')])
 assert np.allclose(Kp,K,atol=1e-5),'PROCAR k list differs from EIGENVAL'
 np.savez_compressed(out,k=K,E=E[:,bands],bands=np.array(bands),nel=nel,s_A=W[...,0],s_B=W[...,1],s_Te=W[...,2],tot=T,
  mesh=np.array(json.loads((u/'dense_grid_manifest.json').read_text())['mesh']),A=A,B=B)

def frac(z,i,b):
 t=z['tot'][i,b]
 return {k:100*float(z[f's_{k}'][i,b])/t for k in ('A','B','Te')}

def analyze(D):
 D=Path(D);src=D/sm.OUT/'逐k投影';rows=[];hse=[]
 info=json.loads((D/sm.OUT/'有效质量_能谷_汇总.json').read_text())['materials']
 for m in sm.M:
  p=src/f'{m}.npz'
  if not p.exists():print('missing',p);continue
  z=np.load(p);mesh=z['mesh'];idx=np.mod(np.rint(z['k']*mesh).astype(int),mesh);nel=int(z['nel']);bands=list(z['bands'])
  find=lambda kf:int(np.argmin(np.abs(((idx/mesh)-np.mod(np.array(kf),1)+0.5)%1-0.5).sum(1)))
  for edge,b in (('VBM',bands.index(nel-1)),('CBM',bands.index(nel))):
   for st in info[m]['valleys'][edge]['stars']:
    i=find(st['k_frac']);f=frac(z,i,b)
    rows.append({'material':m,'edge':edge,'star_dE_meV':st['dE_meV'],'multiplicity':st['mult'],'k_frac':st['k_frac'],
                 **{f's_{k}_pct':round(v,2) for k,v in f.items()},'s_total_pct':round(sum(f.values()),2)})
  # PBE s at the k where the HSE path put the edge (07 cache), against the HSE value at that k.
  st=list(csv.DictReader(open(D/'07_带边能窗/入选态清单.csv',encoding='utf-8-sig')))
  pct=list(csv.DictReader(open(D/'07_带边能窗/带边能窗统计.csv',encoding='utf-8-sig')))
  for edge,b in (('VBM',bands.index(nel-1)),('CBM',bands.index(nel))):
   # 1 meV window = the exact HSE path edge (k location from the state list, percentages from the window table).
   rec=[r for r in st if r['material']==m and r['edge']==edge and float(r['window_meV'])==1.0]
   win=[r for r in pct if r['material']==m and r['edge']==edge and float(r['window_meV'])==1.0]
   if not rec or not win:continue
   kf=[float(rec[0][x]) for x in ('kx','ky','kz')];i=find(kf);f=frac(z,i,b);A=str(z['A']);B=str(z['B']);w=win[0]
   hse.append({'material':m,'edge':edge,'hse_edge_k':kf,'nearest_dense_k':np.round(idx[i]/mesh,4).tolist(),
               'hse_A_s':round(float(w[f'{A}-s']),2),'pbe_A_s_same_k':round(f['A'],2),
               'hse_total_s':round(sum(float(w[f'{x}-s']) for x in (A,B,'Te')),2),'pbe_total_s_same_k':round(sum(f.values()),2)})
 out=D/sm.OUT
 for name,data in (('能谷分辨s占比.csv',rows),('同k点_PBE对HSE.csv',hse)):
  if data:
   with open(out/name,'w',newline='',encoding='utf-8-sig') as fh:
    w=csv.DictWriter(fh,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 for r in hse:print(r)

if __name__=='__main__':
 {'extract':lambda:extract(sys.argv[2],sys.argv[3]),'analyze':lambda:analyze(sys.argv[2])}[sys.argv[1]]()
