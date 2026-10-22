"""Band-edge s fraction vs Seebeck: Seebeck mass m*_S, DOS mass m*_D and valley stars from the SOC dense grid.
Run: python seebeck_mass.py DELIVERY
Inputs reused as-is: step3_uniform EIGENVAL/OUTCAR (full 28x28x6/7 grid, ISYM=-1), step8.1 BoltzTraP2 CRTA,
07_带边能窗 HSE+SOC path windows, 11_续算审计_DeepSeek A/B/C grades of the fixed-concentration points.
BoltzTraP2 mu is measured from the scissored mid-gap: E - E_VBM = mu + target_gap/2 (see audit 4.2).
Optional: step3_uniform/PROCAR[.gz] of the same dense run (LORBIT=11) -> DOS- and sigma-weighted s fractions on the transport Hamiltonian.
Optional: --dense-s CSV, the s_* columns of an earlier run made with the dense PROCAR (e.g. 17_ on the cluster), merged by (material, T, carrier, n).
"""
from pathlib import Path
import json,csv,sys,gzip,math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
COL={'Pb2Sb2Te5':'#2a78d6','Pb2Bi2Te5':'#eb6834','Sn2Sb2Te5':'#1baf7a','Sn2Bi2Te5':'#eda100'}
MK={'Pb2Sb2Te5':'o','Pb2Bi2Te5':'s','Sn2Sb2Te5':'D','Sn2Bi2Te5':'^'}
KB=8.617333262e-5;H=6.62607015e-34;ME=9.1093837015e-31;KBJ=1.380649e-23
LEVELS=[1e19,3e19,1e20]
OUT='12_s占比与Seebeck有效质量'

def F(j,eta):
 if eta<-30:return math.exp(eta)*math.gamma(j+1)
 return quad(lambda x:x**j/(1+math.exp(x-eta)) if x-eta<700 else 0.,0,max(eta,0)+60,limit=200)[0]
def spb_S(eta,r):
 """|S| in V/K for tau ~ E^r (r=0 CRTA, r=-1/2 acoustic phonon)."""
 return KB*((r+2.5)*F(r+1.5,eta)/((r+1.5)*F(r+0.5,eta))-eta)
def eta_from_S(S,r):return brentq(lambda e:spb_S(e,r)-abs(S),-25,60)
def mass_from(n_cm3,eta,T):return (H**2/(2*KBJ*T))*(abs(n_cm3)*1e6/(4*np.pi*F(0.5,eta)))**(2/3)/ME
def spb_n(m,eta,T):return 4*np.pi*(2*m*ME*KBJ*T/H**2)**1.5*F(0.5,eta)/1e6
def fermi(E,mu,T):return 0.5*(1-np.tanh((E-mu)/(2*KB*T)))
def dfermi(E,mu,T):return 1/(4*KB*T*np.cosh((E-mu)/(2*KB*T))**2)
def S_mass(S,n,T,r):
 try:return mass_from(n,eta_from_S(S,r),T)
 except ValueError:return None

def read_eigenval(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:
  for _ in range(5):next(f)
  nel,nk,nb=map(int,next(f).split());K=np.zeros((nk,3));E=np.zeros((nk,nb))
  for k in range(nk):
   line=next(f)
   while not line.strip():line=next(f)
   K[k]=list(map(float,line.split()[:3]))
   for b in range(nb):E[k,b]=float(next(f).split()[1])
 return nel,K,E
def lattice(p):return np.array([list(map(float,l.split()[:3])) for l in p.read_text().splitlines()[2:5]])

def grid(K,mesh):
 idx=np.mod(np.rint(K*mesh).astype(int),mesh);assert len({tuple(i) for i in idx})==len(K)==mesh.prod(),'EIGENVAL is not a full uniform grid'
 return idx
def pockets(idx,E,mesh,band_ids,edge,ref,window):
 """Local extrema of the edge bands on the periodic grid; Kramers partners and plateau duplicates share one k cluster."""
 G=np.empty(tuple(mesh)+(len(band_ids),));G[idx[:,0],idx[:,1],idx[:,2]]=E[:,band_ids]
 s=1 if edge=='VBM' else -1;X=s*G;hits=[]
 for bi in range(len(band_ids)):
  x=X[...,bi];ok=np.ones(x.shape,bool)
  for d in [(a,b,c) for a in (-1,0,1) for b in (-1,0,1) for c in (-1,0,1) if (a,b,c)!=(0,0,0)]:ok&=x>=np.roll(x,d,axis=(0,1,2))-1e-6
  hits+=[(tuple(ijk),G[tuple(ijk)+(bi,)]) for ijk in np.argwhere(ok) if abs(G[tuple(ijk)+(bi,)]-ref)<=window]
 cl=[]
 for ijk,e in sorted(hits,key=lambda h:-s*h[1]):
  for c in cl:
   if all(min(abs(a-b),m-abs(a-b))<=1 for a,b,m in zip(ijk,c['ijk'],mesh)):
    if all(abs(e-x)>2e-3 for x in c['E']):c['E'].append(e)
    break
  else:cl.append({'ijk':ijk,'E':[e],'k_frac':[round(i/m if i<=m//2 else i/m-1,4) for i,m in zip(ijk,mesh)]})
 return cl
def stars(cl,ref,T=300.):
 """Group pockets by depth (2 meV) -> symmetry stars; N_v,eff = sum mult*exp(-dE/kT)."""
 st=[]
 for c in sorted(cl,key=lambda c:abs(c['E'][0]-ref)):
  dE=1000*abs(c['E'][0]-ref)
  if st and abs(st[-1]['dE_meV']-dE)<2:st[-1]['mult']+=1
  else:st.append({'dE_meV':round(dE,1),'mult':1,'k_frac':c['k_frac']})
 return st,sum(x['mult']*math.exp(-x['dE_meV']/1000/(KB*T)) for x in st)

def bt2(D,m):return json.loads((D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step8.1_boltztrap/boltztrap_crta.json').read_text())
def grades(D):
 p=D/'11_续算审计_DeepSeek/csv/bt2_interpolated_S_sigma.csv'
 if not p.exists():return {}
 return {(r['material'],float(r['temperature_K']),float(r['target_carrier_conc_cm-3'])):r['reliability'] for r in csv.DictReader(open(p,encoding='utf-8-sig'))}

def read_procar_s(p,bands,groups):
 """First (charge) block only; s weight of each ion group and total projected weight for the selected bands."""
 op=gzip.open if p.suffix=='.gz' else open;bset={b:i for i,b in enumerate(bands)}
 with op(p,'rt') as f:
  next(f);h=next(f).replace(':',' ').split();nk,nb,ni=int(h[3]),int(h[7]),int(h[11])
  W=np.zeros((nk,len(bands),len(groups)));T=np.zeros((nk,len(bands)));K=np.zeros((nk,3));k=b=-1;want=False;done=False
  for line in f:
   t=line.lstrip()
   if t.startswith('k-point'):
    k+=1;s=line.split(':')[1].split('weight')[0].replace('-',' -');K[k]=list(map(float,s.split()[:3]))
   elif t.startswith('band'):b=int(line.split()[1])-1;want=b in bset;done=False
   elif t.startswith('ion') and want and not done:
    for j in range(ni):
     row=next(f).split();v=np.array(list(map(float,row[1:10])));T[k,bset[b]]+=v.sum()
     for g,ions in enumerate(groups):
      if j in ions:W[k,bset[b],g]+=v[0]
    done=True
 assert k==nk-1
 return K,W,T

def analyse(D,m,grade):
 base=D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result';u=base/'step3_uniform'
 nel,K,E=read_eigenval(u/'EIGENVAL');mesh=np.array(json.loads((u/'dense_grid_manifest.json').read_text())['mesh'])
 j=bt2(D,m);sc=j['scissor'];assert sc['n_valence_bands']==nel;V=j['cell_volume_cm3'];idx=grid(K,mesh)
 vbm=E[:,:nel].max();cbm_dft=E[:,nel:].min();gap_dft=cbm_dft-vbm;assert abs(gap_dft-sc['dft_gap_eV'])<2e-3,(m,gap_dft,sc)
 Es=E.copy();Es[:,nel:]+=sc['shift_eV'];cbm=cbm_dft+sc['shift_eV'];mid=vbm+sc['target_gap_eV']/2
 out={'material':m,'VBM_eV':vbm,'CBM_scissored_eV':cbm,'midgap_eV':mid,'gap_dft_eV':gap_dft,'gap_target_eV':sc['target_gap_eV'],'mesh':mesh.tolist(),'nk':len(K)}
 nv={}
 for edge,bands,ref in [('VBM',list(range(nel-6,nel)),vbm),('CBM',list(range(nel,nel+6)),cbm_dft)]:
  cl=pockets(idx,E,mesh,bands,edge,ref,0.1);st,neff=stars(cl,ref)
  nv[edge]={'pockets_26meV':sum(x['mult'] for x in st if x['dE_meV']<=26),'pockets_52meV':sum(x['mult'] for x in st if x['dE_meV']<=52),
   'pockets_100meV':sum(x['mult'] for x in st),'Nv_eff_300K':round(neff,2),'stars':st}
 out['valleys']=nv
 mu=np.array(j['mu_rel_efermi_eV'])+mid
 # Neutral point: BoltzTraP2 c(mu)=0 vs the same condition on the raw grid; both sit off mid-gap when the band DOS differ.
 cn={}
 for ti,T in enumerate(j['temperatures_K']):
  N=np.array(j['carrier_conc_cm-3'][ti]);i=int(np.where(np.diff(np.sign(N))!=0)[0][0]);x=mu[i]-N[i]*(mu[i+1]-mu[i])/(N[i+1]-N[i])
  g=brentq(lambda z:(1-fermi(Es[:,:nel],z,T)).sum()-fermi(Es[:,nel:],z,T).sum(),vbm-0.5,cbm+0.5)
  cn[str(int(T))]={'bt2_c0_minus_midgap_meV':round(1000*(x-mid),1),'grid_neutral_minus_midgap_meV':round(1000*(g-mid),1),'bt2_c_at_midgap_cm3':float(np.interp(mid,mu,N))}
 out['neutral_point']=cn
 pro=next((p for p in (u/'PROCAR',u/'PROCAR.gz') if p.exists()),None);ws=None
 if pro:
  el=(D/'A2B2Te5_SOC_dense/materials'/m/'POSCAR').read_text().splitlines()[5:7];sp=sum(([e]*int(c) for e,c in zip(el[0].split(),el[1].split())),[])
  An='Pb' if m.startswith('Pb') else 'Sn';Bn='Sb' if 'Sb' in m else 'Bi'
  bands=list(range(nel-8,nel+8));Kp,W,Tp=read_procar_s(pro,bands,[{i for i,e in enumerate(sp) if e==x} for x in (An,Bn,'Te')])
  assert np.allclose(Kp,K,atol=1e-5),'PROCAR k list differs from EIGENVAL'
  rec=2*np.pi*np.linalg.inv(lattice(u/'POSCAR')).T;G=np.empty(tuple(mesh)+(len(bands),));G[idx[:,0],idx[:,1],idx[:,2]]=Es[:,bands]
  dEdf=np.stack([(np.roll(G,-1,a)-np.roll(G,1,a))*mesh[a]/2 for a in range(3)],-1);v2=((dEdf@np.linalg.inv(rec).T)**2).sum(-1)[idx[:,0],idx[:,1],idx[:,2]]
  ws=(bands,Es[:,bands],W,Tp,v2);out['dense_procar']=str(pro.relative_to(D))
 else:out['dense_procar']=None
 rows=[]
 for ti,T in enumerate(j['temperatures_K']):
  N=np.array(j['carrier_conc_cm-3'][ti]);S=np.array(j['seebeck_V_per_K'][ti]);Sin=(np.array(j['seebeck_xx_V_per_K'][ti])+np.array(j['seebeck_yy_V_per_K'][ti]))/2
  Szz=np.array(j['seebeck_zz_V_per_K'][ti]) if 'seebeck_zz_V_per_K' in j else 3*S-2*Sin
  for carrier,sgn in [('hole',1),('electron',-1)]:
   sel=(sgn*N>0)&((mu<=vbm+0.3)&(mu>=vbm-0.6) if sgn>0 else (mu>=cbm-0.3)&(mu<=cbm+0.6))
   o=np.argsort(np.abs(N[sel]));n_s=np.abs(N[sel])[o];pick=lambda a:a[sel][o]
   for lv in LEVELS:
    if not n_s.min()<=lv<=n_s.max():continue
    i=np.searchsorted(n_s,lv);w=(np.log(lv)-np.log(n_s[i-1]))/(np.log(n_s[i])-np.log(n_s[i-1]));lin=lambda a:a[i-1]+w*(a[i]-a[i-1])
    Sv,Siv,muv=lin(pick(S)),lin(pick(Sin)),lin(pick(mu));edge_e=vbm if sgn>0 else cbm;eta=sgn*(edge_e-muv)/(KB*T)
    p_h=(1-fermi(Es[:,:nel],muv,T)).sum()/len(K)/V;n_e=fermi(Es[:,nel:],muv,T).sum()/len(K)/V
    r={'material':m,'T_K':T,'carrier':carrier,'n_cm3':lv,'audit_grade':grade.get((m,T,sgn*lv),'not_graded'),
     'S_scalar_uV_K':1e6*Sv,'S_in_plane_uV_K':1e6*Siv,'S_zz_uV_K':1e6*lin(pick(Szz)),'mu_minus_edge_meV':1000*(muv-edge_e),'eta':eta,
     'minority_over_majority_grid':(n_e/p_h if sgn>0 else p_h/n_e),'m_DOS':mass_from(lv,eta,T),'sign_ok':bool(np.sign(Sv)==sgn and np.sign(Siv)==sgn)}
    for tag,Sx in (('scalar',Sv),('in_plane',Siv)):
     r[f'm_S_CRTA_{tag}']=S_mass(Sx,lv,T,0.) if r['sign_ok'] else None
     r[f'm_S_APS_{tag}']=S_mass(Sx,lv,T,-0.5) if r['sign_ok'] else None
    r['mS_over_mD_in_plane']=r['m_S_CRTA_in_plane']/r['m_DOS'] if r['m_S_CRTA_in_plane'] else None
    if ws:
     bands,Eb,W,Tp,v2=ws;side=(np.array(bands)<nel) if sgn>0 else (np.array(bands)>=nel)
     for tag,wt in (('dos',dfermi(Eb,muv,T)),('sigma',dfermi(Eb,muv,T)*v2)):
      wt=wt*side;den=(wt*Tp).sum()
      for g,name in enumerate(('A','B','Te')):r[f's_{name}_{tag}_pct']=100*(wt*W[...,g]).sum()/den
      r[f's_cation_{tag}_pct']=r[f's_A_{tag}_pct']+r[f's_B_{tag}_pct'];r[f's_total_{tag}_pct']=r[f's_cation_{tag}_pct']+r[f's_Te_{tag}_pct']
    rows.append(r)
 return out,rows

def windows(D):
 res={}
 for r in csv.DictReader(open(D/'07_带边能窗/带边能窗统计.csv',encoding='utf-8-sig')):
  m=r['material'];A='Pb' if m.startswith('Pb') else 'Sn';B='Sb' if 'Sb' in m else 'Bi';g=lambda k:float(r[k] or 0)
  res[(m,r['edge'],float(r['window_meV']))]={'A_s':g(A+'-s'),'B_s':g(B+'-s'),'Te_s':g('Te-s'),'cation_s':g(A+'-s')+g(B+'-s'),'total_s':g(A+'-s')+g(B+'-s')+g('Te-s'),'Te_p':g('Te-p')}
 return res

def merge_dense(rows,path):
 """Copy s_* columns from a dense-PROCAR run into rows that lack them (same script, same (material, T, carrier, n) keys)."""
 if not path or any('s_total_dos_pct' in r for r in rows):return
 src={(r['material'],float(r['T_K']),r['carrier'],float(r['n_cm3'])):r for r in csv.DictReader(open(path,encoding='utf-8-sig'))}
 for r in rows:
  x=src.get((r['material'],float(r['T_K']),r['carrier'],float(r['n_cm3'])))
  if x:r.update({k:float(v) for k,v in x.items() if k.startswith('s_') and v not in ('',None)})

def s_at(rows,win,m,edge,carrier,lv):
 """Edge s fractions: dense-grid DOS-weighted at (300 K, lv) when available, else the HSE path 50 meV window."""
 r=[x for x in rows if x['material']==m and x['T_K']==300 and x['carrier']==carrier and x['n_cm3']==lv and 's_total_dos_pct' in x]
 if r:return {'A_s':r[0]['s_A_dos_pct'],'B_s':r[0]['s_B_dos_pct'],'Te_s':r[0]['s_Te_dos_pct'],'cation_s':r[0]['s_cation_dos_pct'],'total_s':r[0]['s_total_dos_pct'],'source':'dense'}
 return {**win[(m,edge,50.)],'source':'hse'}

def plot(D,O,rows,info,win):
 plt.rcParams.update({'font.size':9,'axes.edgecolor':'#52514e','axes.labelcolor':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e'})
 fig,axs=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
 ax=axs[0];nn=np.logspace(18,20.7,60)
 for mref in (0.2,0.5,1.0,2.0):
  Sref=[1e6*spb_S(brentq(lambda e:spb_n(mref,e,300)-x,-25,60),0.) for x in nn];ax.plot(nn,Sref,color='#c3c2b7',lw=1,ls='--',zorder=1)
  k=np.argmin(np.abs(np.array(Sref)-60)) if Sref[-1]<60 else -1;ax.text(nn[k],Sref[k],f'{mref}',color='#52514e',fontsize=7,ha='left',va='bottom')
 for m in M:
  j=bt2(D,m);N=np.array(j['carrier_conc_cm-3'][0]);S=1e6*(np.array(j['seebeck_xx_V_per_K'][0])+np.array(j['seebeck_yy_V_per_K'][0]))/2;mu=np.array(j['mu_rel_efermi_eV'])+info[m]['midgap_eV']
  sel=(N>0)&(mu<=info[m]['VBM_eV']+0.3)&(mu>=info[m]['VBM_eV']-0.6);o=np.argsort(N[sel]);ax.plot(N[sel][o],S[sel][o],color=COL[m],lw=2,label=m,zorder=2)
 ax.axvspan(1e18,1e19,color='#c3c2b7',alpha=.18,lw=0);ax.text(1.1e18,12,'bipolar: audit grade C at 1e18',fontsize=7,color='#52514e')
 ax.set_xscale('log');ax.set_xlim(1e18,5e20);ax.set_ylim(0,420);ax.grid(alpha=.25,lw=.6);ax.set_xlabel('hole concentration p (cm$^{-3}$)');ax.set_ylabel('in-plane S (µV/K)')
 ax.set_title('a  Pisarenko, 300 K, BoltzTraP2 CRTA (SOC), in-plane\ndashed: single parabolic band, m* = 0.2 / 0.5 / 1 / 2 m$_e$',loc='left',fontsize=9);ax.legend(fontsize=8,frameon=False,loc='upper right')
 off={'hole':{'Pb2Sb2Te5':(10,-4),'Pb2Bi2Te5':(-12,-6),'Sn2Sb2Te5':(10,10),'Sn2Bi2Te5':(-10,14)},
  'electron':{'Pb2Sb2Te5':(10,-26),'Pb2Bi2Te5':(-12,-26),'Sn2Sb2Te5':(-12,6),'Sn2Bi2Te5':(-12,0)}}
 for ax,edge,carrier,tag in [(axs[1],'VBM','hole','b'),(axs[2],'CBM','electron','c')]:
  top=0
  for m in M:
   r=[x for x in rows if x['material']==m and x['T_K']==300 and x['carrier']==carrier and x['n_cm3']==3e19]
   if not r or r[0]['m_S_CRTA_in_plane'] is None:continue
   s=s_at(rows,win,m,edge,carrier,3e19)['total_s'];r=r[0];mS=r['m_S_CRTA_in_plane'];hi=max(mS,r['m_DOS']);top=max(top,hi)
   ax.plot([s,s],[r['m_DOS'],mS],color='#c3c2b7',lw=1,zorder=1)
   ax.plot(s,mS,MK[m],color=COL[m],ms=10,mec='white',mew=2,zorder=3);ax.plot(s,r['m_DOS'],MK[m],mfc='white',mec=COL[m],mew=2,ms=9,zorder=2)
   dx,dy=off[carrier][m]
   ax.annotate(f"{m}\nN$_{{v,eff}}$={info[m]['valleys'][edge]['Nv_eff_300K']:.1f}",(s,hi if dy>=0 else min(mS,r['m_DOS'])),xytext=(dx,dy),textcoords='offset points',fontsize=8,color='#0b0b0b',ha='left' if dx>0 else 'right',va='bottom' if dy>=0 else 'top')
  xs=[s_at(rows,win,m,edge,carrier,3e19) for m in M];src=xs[0]['source']
  ax.grid(alpha=.25,lw=.6);ax.set_ylabel('effective mass (m$_e$)');ax.set_ylim(0,top*1.25)
  ax.set_xlabel(f'{edge} total s (A-s + B-s + Te-s, %)\n'+('PBE+SOC dense grid, $-\\partial f/\\partial E$ weighted at the same (n, T)' if src=='dense' else 'HSE+SOC path PROCAR, 50 meV inner window'))
  ax.set_xlim(min(x['total_s'] for x in xs)*0.6,max(x['total_s'] for x in xs)*1.25)
  ax.set_title(f'{tag}  {carrier}s, 300 K, |n| = 3e19 cm$^{{-3}}$\nfilled: m*$_S$ from in-plane S (CRTA)   open: m*$_D$ from n(µ)',loc='left',fontsize=9)
 for ext in ('png','svg'):fig.savefig(O/f's占比_Seebeck有效质量.{ext}',dpi=180)
 plt.close(fig)

def main(D,dense=None):
 O=D/OUT;O.mkdir(exist_ok=True)
 info={};rows=[];grade=grades(D)
 for m in M:
  i,r=analyse(D,m,grade);info[m]=i;rows+=r
 merge_dense(rows,dense)
 win=windows(D);keys=list(dict.fromkeys(k for r in rows for k in r))
 with open(O/'有效质量_Seebeck.csv','w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 tab=[]
 for m in M:
  for edge,carrier in [('VBM','hole'),('CBM','electron')]:
   r={x['n_cm3']:x for x in rows if x['material']==m and x['T_K']==300 and x['carrier']==carrier};v=info[m]['valleys'][edge]
   t={'material':m,'edge':edge,**{f'{k}_{int(wd)}meV':round(win[(m,edge,wd)][k],2) for wd in (25.,50.,100.) for k in ('A_s','B_s','Te_s','total_s')},
      'pockets_26meV':v['pockets_26meV'],'pockets_52meV':v['pockets_52meV'],'Nv_eff_300K':v['Nv_eff_300K']}
   sd=s_at(rows,win,m,edge,carrier,3e19)
   if sd['source']=='dense':t.update({f'dense_{k}_3e19':round(sd[k],2) for k in ('A_s','B_s','Te_s','total_s')})
   for lv in LEVELS:
    x=r.get(lv,{});tag=f'{lv:.0e}'.replace('+','')
    t.update({f'grade_{tag}':x.get('audit_grade'),f'S_in_plane_{tag}':x.get('S_in_plane_uV_K'),f'mS_{tag}':x.get('m_S_CRTA_in_plane'),f'mD_{tag}':x.get('m_DOS'),f'mS_APS_{tag}':x.get('m_S_APS_in_plane')})
   tab.append(t)
 with open(O/'s占比_有效质量对照表.csv','w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(tab[0].keys()));w.writeheader();w.writerows(tab)
 corr={}
 for edge,carrier in [('VBM','hole'),('CBM','electron')]:
  for k in ('A_s','cation_s','total_s'):
   for q in ('m_S_CRTA_in_plane','m_DOS','absS_in_plane'):
    for lv in (3e19,1e20):
     xs,ys=[],[]
     for m in M:
      r=[x for x in rows if x['material']==m and x['T_K']==300 and x['carrier']==carrier and x['n_cm3']==lv]
      if r and (q.startswith('abs') or r[0][q] is not None):sd=s_at(rows,win,m,edge,carrier,lv);src=sd['source'];xs.append(sd[k]);ys.append(abs(r[0]['S_in_plane_uV_K']) if q.startswith('abs') else r[0][q])
     if len(xs)>=3:
      rx,ry=np.argsort(np.argsort(xs)),np.argsort(np.argsort(ys))
      corr[f'{carrier}:{k}_{src}~{q}@{lv:.0e}']={'n':len(xs),'pearson':round(float(np.corrcoef(xs,ys)[0,1]),3),'spearman':round(float(np.corrcoef(rx,ry)[0,1]),3)}
 (O/'有效质量_能谷_汇总.json').write_text(json.dumps({'materials':info,'correlations_n4':corr,
  'notes':['mu alignment: BoltzTraP2 mu_rel is measured from the scissored mid-gap (VBM + target_gap/2), not from the OUTCAR Fermi level',
   'm_S_CRTA: SPB inversion of BoltzTraP2 S(n) with energy-independent tau (same assumption as the input); m_S_APS: tau~E^-1/2, the convention for experimental Seebeck mass',
   'in_plane: (S_xx+S_yy)/2, the direction for layered films; scalar: trace/3 as stored in the JSON',
   'm_DOS: parabolic mass reproducing BoltzTraP2 n(mu) at the reduced Fermi level measured from the EIGENVAL band edge (scissored CBM)',
   'audit_grade: A/B/C grade of the same fixed-concentration point in 11_续算审计_DeepSeek (3e19 is not graded there)',
   'minority_over_majority_grid: Fermi occupation on the raw 28x28xNz grid at the same mu; the raw grid overweights light pockets, so it is an upper-bound bipolar indicator',
   'neutral_point: offset of c(mu)=0 from mid-gap in BoltzTraP2 and on the raw grid; a non-zero offset follows from unequal VB/CB DOS',
   'valleys: local extrema of six edge bands on the full grid within 100 meV; Kramers partners and plateaus merged; stars grouped within 2 meV; grid-resolution energies, ring-like maxima split into grid points',
   's fractions: dense = PBE+SOC 28x28xNz PROCAR, -df/dE weighted at the same (n, T) (from the PROCAR or --dense-s); hse = 07 HSE+SOC path PROCAR inner windows, a different Hamiltonian and edge k; figure, table dense_* columns and correlations use dense when present',
   'four materials only: correlation coefficients are descriptive, not statistical evidence']},ensure_ascii=False,indent=1,default=float),encoding='utf-8')
 plot(D,O,rows,info,win)
 for m in M:print(m,'neutral point',info[m]['neutral_point']['300'])
 for t in tab:print(t['material'],t['edge'],'s_total50=',t['total_s_50meV'],'Nv_eff=',t['Nv_eff_300K'],*(f"{k}={t[k]:.3g}" if isinstance(t[k],float) else f"{k}={t[k]}" for k in ('grade_3e19','S_in_plane_3e19','mS_3e19','mD_3e19','grade_1e20','S_in_plane_1e20','mS_1e20','mD_1e20')))
 print(json.dumps({k:v for k,v in corr.items() if 'total' in k},indent=0))

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('delivery');ap.add_argument('--dense-s',help='CSV with s_* columns from a dense-PROCAR run')
 a=ap.parse_args();main(Path(a.delivery).resolve(),a.dense_s)
