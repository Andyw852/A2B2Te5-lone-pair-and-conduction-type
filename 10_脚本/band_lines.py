"""Plain line band structures: the old HSE06+SOC path (07_) and the PBE-D3+SOC dense grid with the rigid scissor used for transport.
Run: python band_lines.py DELIVERY      (one BoltzTraP2 fit per material for the PBE lines, about 3 min and a few GB each)
Bands are assigned by electron count (SOC: one electron per band), not by a Fermi level, so the gap is max(VB) < min(CB) over every k.
Energy zero: the intrinsic Fermi level, taken at mid-gap (undoped, T -> 0); doping moves E_F toward a band edge, and the transport scans it.
Left: HSE along the whole path. Right: zoom around Gamma-A, HSE against PBE + scissor (the band structure the transport used),
each with its own mid-gap at zero. PBE lines are BoltzTraP2 interpolations of the 28x28x7 (Bi: 28x28x6) grid evaluated on the HSE path.
The Bi compounds' HSE path used zero-weight k-points in the SCF; there the top valence band at k off the 12x12x2 SCF grid lies 0.1-0.5 eV
below the on-grid values (spikes at every grid point), so only the on-grid points (dots) are trusted. The Sb compounds used KPOINTS_OPT: smooth.
"""
from pathlib import Path
import sys,json
import numpy as np
import ase.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from BoltzTraP2 import sphere,fite
sys.path.insert(0,str(Path(__file__).resolve().parent))
import seebeck_mass as sm
import bt2_local as bl

OUT='12_s占比与Seebeck有效质量/能带线图_HSE与PBE'
HSE_C,PBE_C,GAP_C,INK,INK2,GRID='#2a78d6','#eb6834','#e9e8e4','#0b0b0b','#52514e','#d9d8d4'

def label(k):
 f=np.mod(np.round(np.asarray(k,float)*6,6),6)/6;ab=tuple(np.round(f[:2],4));z=round(f[2],4)
 base={(0.,0.):'Γ',(0.5,0.):'M',(0.,0.5):'M',(0.5,0.5):'M',(0.3333,0.3333):'K',(0.6667,0.6667):'K'}.get(ab)
 if base is None or z not in (0.,0.5):return None
 return base if z==0 else {'Γ':'A','M':'L','K':'H'}[base]

def path_x(K,cell):
 kc=K@cell.reciprocal();d=np.r_[0,np.linalg.norm(np.diff(kc,axis=0),axis=1)];brk=d>4*np.median(d[1:]);d[brk]=0
 return np.cumsum(d),np.where(brk)[0]

def pbe_on_path(D,m,K,window=1.5):
 u=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform';nel,Kd,E=sm.read_eigenval(u/'EIGENVAL')
 vbm,cbm=E[:,:nel].max(),E[:,nel:].min();keep=np.where(((E>vbm-window)&(E<cbm+window)).any(0))[0];E=E[:,keep];nv=nel-int(keep[0])
 d=bl.Data(ase.io.read(u/'POSCAR',format='vasp'),Kd,E);eq=sphere.get_equivalences(d.atoms,None,len(Kd)*5)
 eb,_=fite.getBands(K,eq,d.get_lattvec(),fite.fitde3D(d,eq))
 shift=sm.bt2(D,m)['scissor']['shift_eV'];eb=eb.T*bl.HA-vbm;eb[:,nv:]+=shift    # BoltzTraP2 works in Ha; (nk, nband) relative to the PBE VBM, CB scissored
 return eb,nv,cbm-vbm,cbm-vbm+shift

def hse_path(D,m):
 """HSE energies in path order: 解析投影.npz stores each distinct k once and not in path order, so map the KPOINTS(_OPT) path onto it."""
 s=D/'07_带边能窗/源数据'/m;h=np.load(s/'解析投影.npz');f=s/'KPOINTS_OPT' if (s/'KPOINTS_OPT').exists() else s/'KPOINTS'
 L=f.read_text().split('\n');P=np.array([[float(x) for x in l.split()[:4]] for l in L[3:3+int(L[1])]])
 P=P[:,:3] if f.name=='KPOINTS_OPT' else P[P[:,3]==0,:3]                       # Bi compounds: 38 weighted SCF points, then the zero-weight path
 idx=[int(np.argmin(np.linalg.norm(h['K']-p,axis=1))) for p in P];assert max(np.linalg.norm(h['K'][idx]-P,axis=1))<1e-6
 return P,h['E'][idx],f.name=='KPOINTS'

def main(D):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,
                      'axes.edgecolor':INK2,'xtick.color':INK2,'ytick.color':INK2,'axes.labelcolor':INK})
 fig,axs=plt.subplots(4,2,figsize=(12.5,15.5),gridspec_kw={'width_ratios':[2.3,1]});summary={}
 for row,m in enumerate(sm.M):
  K,Eh,zw=hse_path(D,m)
  ongrid=np.array([np.allclose(np.mod(k[:2]*12+1e-6,1),0,atol=1e-4) and np.allclose(np.mod(k[2]*2+1e-6,1),0,atol=1e-4) for k in K])
  nel=sm.read_eigenval(D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step3_uniform/EIGENVAL')[0]
  hv=Eh[:,nel-1].max();Eh=Eh-hv;hg=Eh[:,nel].min()
  cell=ase.io.read(D/'07_带边能窗/源数据'/m/'POSCAR',format='vasp').cell;x,brk=path_x(K,cell)
  lab=[(x[i],label(k)) for i,k in enumerate(K) if label(k)];ep,nv,pg,pgs=pbe_on_path(D,m,K)
  ig=[i for i,k in enumerate(K) if label(k)=='Γ'];ia=[i for i,k in enumerate(K) if label(k)=='A']
  summary[m]={'HSE_gap_path_eV':round(float(hg),4),'HSE_VBM_k':K[Eh[:,nel-1].argmax()].round(3).tolist(),'HSE_CBM_k':K[Eh[:,nel].argmin()].round(3).tolist(),
              'PBE_gap_dense_eV':round(float(pg),4),'PBE_gap_after_scissor_eV':round(float(pgs),4),
              'HSE_Gamma_VB_CB_eV':[round(float(Eh[ig[0],nel-1]),3),round(float(Eh[ig[0],nel]),3)],'HSE_A_VB_CB_eV':[round(float(Eh[ia[0],nel-1]),3),round(float(Eh[ia[0],nel]),3)]}
  Eh=Eh-hg/2;ep=ep-pgs/2                                                        # plot with E_F at mid-gap; the summary above is relative to the VBM
  segs=np.split(np.arange(len(K)),brk)
  for col,ax in enumerate(axs[row]):
   ax.axhspan(-hg/2,hg/2,color=GAP_C,lw=0,zorder=0)
   for s in segs:
    for b in range(Eh.shape[1]):ax.plot(x[s],Eh[s,b],color=HSE_C,lw=1.1,zorder=3,alpha=0.55 if zw else 1)
    if col==1:
     for b in range(ep.shape[1]):ax.plot(x[s],ep[s,b],color=PBE_C,lw=1.1,ls=(0,(4,2)),zorder=2)
   if zw:
    for b in range(Eh.shape[1]):ax.plot(x[ongrid],Eh[ongrid,b],'o',ms=4.2,color=HSE_C,mec='#fcfcfb',mew=0.8,zorder=4)
   for xv,t in lab:ax.axvline(xv,color=GRID,lw=0.8,zorder=1)
   for i in brk:ax.axvline(x[i],color=INK2,lw=1.0,zorder=1)
   ax.axhline(0,color=INK2,lw=0.8,ls=(0,(2,2)),zorder=1)
   seen={};[seen.setdefault(round(xv,4),t) for xv,t in lab]
   ax.set_xticks(list(seen));ax.set_xticklabels(list(seen.values()))
   ax.tick_params(length=0,axis='x');ax.set_ylabel('E − E_F (eV)')
  ax0,ax1=axs[row];ax0.set_xlim(x[0],x[-1]);ax0.set_ylim(-1.5,1.5)
  ax0.set_title(f'{m}   HSE06+SOC   E_g = {1000*hg:.0f} meV (path)',loc='left',fontsize=11,color=INK)
  if zw:ax0.text(0.995,0.97,'zero-weight path: dots = SCF-grid k (trusted),\nlines between them biased low in the VB',transform=ax0.transAxes,ha='right',va='top',fontsize=8.5,color=INK2,zorder=5,bbox=dict(facecolor='#fcfcfb',edgecolor='none',pad=1.5))
  ax0.text(x[-1],0,f' E_F (mid-gap), gap {1000*hg:.0f} meV',va='center',ha='right',fontsize=8.5,color=INK2,zorder=4,
           bbox=dict(facecolor='#fcfcfb',edgecolor='none',pad=1))
  xg,xa=x[ig[1] if len(ig)>1 else ig[0]],x[ia[0]];w=xa-xg;ax1.set_xlim(xg-2.2*w,xa+2.2*w);ax1.set_ylim(-0.45,0.45)
  ax1.set_title(f'near Γ–A;  PBE gap {1000*pg:.0f} → {1000*pgs:.0f} meV after scissor',loc='left',fontsize=9.5,color=INK2)
  print(m,summary[m])
 h1,=axs[0,1].plot([],[],color=HSE_C,lw=1.6,label='HSE06+SOC (07_ path)');h2,=axs[0,1].plot([],[],color=PBE_C,lw=1.6,ls=(0,(4,2)),label='PBE-D3+SOC + rigid scissor (used for transport)')
 fig.legend(handles=[h1,h2],loc='upper center',ncol=2,frameon=False,fontsize=10,bbox_to_anchor=(0.5,0.995))
 fig.text(0.5,0.005,'Dashed line: E_F at mid-gap (undoped). Shaded: HSE gap. Bands are counted by electrons (one per band with SOC), so a positive gap means a semiconductor; each line is a Kramers pair.',
          ha='center',fontsize=9,color=INK2)
 fig.tight_layout(rect=(0,0.012,1,0.975));fig.savefig(D/(OUT+'.png'),dpi=170);fig.savefig(D/(OUT+'.svg'))
 (D/(OUT+'.json')).write_text(json.dumps(summary,ensure_ascii=False,indent=1))

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
