"""Is a rigid scissor enough? Compare the old HSE06+SOC path (07_) with the PBE-D3+SOC dense grid at the k-points both contain.
Run: python hse_pbe_bandshape.py DELIVERY
Matched points are the HSE path points that fall exactly on the dense grid (up to k -> -k): the TRIM points, and for the Sb compounds also the
Gamma-A line at kz = j/7. Energies are relative to each functional's own band edges (HSE: path; PBE: dense grid). A rigid scissor is valid only
if Eg_HSE(k) - Eg_PBE(k) is the same at every k. PbBi's HSE run used a slightly different cell (11_ audit section 5.3).
"""
from pathlib import Path
import sys,csv,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm

OUT='12_s占比与Seebeck有效质量/HSE与PBE同k点能带对照.csv'

def main(D):
 rows=[]
 for m in sm.M:
  h=np.load(D/'07_带边能窗/源数据'/m/'解析投影.npz');Kh,Eh=h['K'],h['E']
  base=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform'
  nel,K,E=sm.read_eigenval(base/'EIGENVAL');mesh=np.array(json.loads((base/'dense_grid_manifest.json').read_text())['mesh'])
  hv,hc,pv,pc=Eh[:,nel-1].max(),Eh[:,nel].min(),E[:,nel-1].max(),E[:,nel].min()
  key={tuple(np.mod(np.round(k*mesh),mesh).astype(int)):i for i,k in enumerate(K)};seen=set()
  for i,k in enumerate(Kh):
   g=k*mesh
   if not np.allclose(g,np.round(g),atol=1e-4):continue
   j=next((key[t] for t in (tuple(np.mod(np.round(s*g),mesh).astype(int)) for s in (1,-1)) if t in key),None)
   if j is None or j in seen:continue
   seen.add(j);egh,egp=Eh[i,nel]-Eh[i,nel-1],E[j,nel]-E[j,nel-1]
   rows.append({'material':m,'k':str([round(float(x),3) for x in np.abs(k)]),
     'HSE_VB_below_VBM_meV':round(1000*(hv-Eh[i,nel-1]),1),'HSE_CB_above_CBM_meV':round(1000*(Eh[i,nel]-hc),1),
     'PBE_VB_below_VBM_meV':round(1000*(pv-E[j,nel-1]),1),'PBE_CB_above_CBM_meV':round(1000*(E[j,nel]-pc),1),
     'Eg_HSE_k_meV':round(1000*egh,1),'Eg_PBE_k_meV':round(1000*egp,1),'HSE_minus_PBE_meV':round(1000*(egh-egp),1)})
  r=[x for x in rows if x['material']==m];d=[x['HSE_minus_PBE_meV'] for x in r]
  print(f"{m}: {len(r)} common k; HSE gap {1000*(hc-hv):.0f} meV, PBE gap {1000*(pc-pv):.0f} meV; Eg_HSE(k)-Eg_PBE(k) {min(d):+.0f} .. {max(d):+.0f} meV")
 with open(D/OUT,'w',newline='',encoding='utf-8-sig') as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
