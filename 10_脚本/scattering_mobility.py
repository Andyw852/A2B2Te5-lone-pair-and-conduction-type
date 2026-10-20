"""Mechanism-resolved AMSET mobilities (ADP, IMP, POP) for the four materials, with the inputs that set each one.
Run: python scattering_mobility.py DELIVERY
Mobilities: dense-grid AMSET (28x28x7 / 28x28x6 DFT input, interpolation_factor 1), one folder per material:
  14_散射机制分项迁移率/密网格AMSET/<m>/transport.json   (materials without it are skipped; the 12x12x2 results are void)
Inputs table: 05_SOC_AMSET补算/<m>/04_AMSET输运/settings.yaml (same eps, omega_LO, C_ij, gap as the dense runs),
03_SOC形变势/band_edges.json and the step8.1 BoltzTraP2 JSON (m*_c = n e^2 / (sigma/tau), in-plane, mid-gap mu reference).
"""
from pathlib import Path
import json,csv,sys
import numpy as np
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
MECH=['ADP','IMP','POP']
COL={'ADP':'#2a78d6','IMP':'#eb6834','POP':'#1baf7a','overall':'#0b0b0b'}
MK={'ADP':'o','IMP':'s','POP':'D','overall':'|'}
LABEL={'ADP':'ADP  acoustic deformation potential','IMP':'IMP  ionized impurity','POP':'POP  polar optical phonon','overall':'total (all three, rates summed per state)'}
OUT='14_散射机制分项迁移率'
CAP='AMSET on the 28x28x7 (Bi: 28x28x6) SOC bands, interpolation factor 1; deformation potentials from 12x12x2. Absolute mobilities may be 5-25% high (mesh convergence, see README).'
E=1.602176634e-19;ME=9.1093837015e-31

def inplane(a):return (a[0][0]+a[1][1])/2
def inputs(D,m):
 base=D/'05_SOC_AMSET补算'/m
 return yaml.safe_load((base/'04_AMSET输运/settings.yaml').read_text()),json.loads((base/'03_SOC形变势/band_edges.json').read_text())
def dense(D,m):
 p=D/OUT/'密网格AMSET'/m/'transport.json'
 return json.loads(p.read_text()) if p.exists() else None

def mass_c(D,m,n_target,T=300.):
 j=json.loads((D/'A2B2Te5_SOC_dense/materials'/m/'ke-dft-cpu/result/step8.1_boltztrap/boltztrap_crta.json').read_text())
 ti=j['temperatures_K'].index(T);N=np.array(j['carrier_conc_cm-3'][ti]);sig=(np.array(j['sigma_over_tau_xx'][ti])+np.array(j['sigma_over_tau_yy'][ti]))/2
 sel=np.sign(N)==np.sign(n_target);o=np.argsort(np.abs(N[sel]))
 # Majority branch only: |N| grows monotonically away from the neutral point on this side.
 nn,ss=np.abs(N[sel])[o],sig[sel][o];keep=np.concatenate([[True],np.diff(nn)>0]);nn,ss=nn[keep],ss[keep]
 sv=np.exp(np.interp(np.log(abs(n_target)),np.log(nn),np.log(ss)))
 return abs(n_target)*1e6*E**2/sv/ME

def main(D):
 O=D/OUT;O.mkdir(exist_ok=True)
 rows,params=[],[]
 for m in M:
  s,be=inputs(D,m);t=dense(D,m)
  for di,d in (enumerate(np.array(t['doping'])) if t else ()):
   Ts=t['temperatures']
   if not any(np.isclose(abs(d),10.**k,rtol=1e-6) for k in range(17,21)):continue
   for ti,T in enumerate(Ts):
    r={'material':m,'carrier':'p' if d>0 else 'n','doping_cm3':float(f'{d:.3e}'),'T_K':T}
    for k in ['overall']+MECH:
     a=np.array(t['mobility'][k])[di,ti];r[f'mu_{k}_inplane']=inplane(a);r[f'mu_{k}_zz']=a[2][2]
    rows.append(r)
  # PyYAML reads exponents without a dot (e.g. -1e-05) as strings, hence the float casts.
  ei,es=np.array(s['high_frequency_dielectric'],float),np.array(s['static_dielectric'],float);C=np.array(s['elastic_constant'],float)
  p={'material':m,'eps_inf_xx':ei[0,0],'eps_inf_zz':ei[2,2],'eps0_xx':es[0,0],'eps0_zz':es[2,2],
     'POP_coupling_xx_1/einf-1/e0':1/ei[0,0]-1/es[0,0],'omega_LO_THz':s['pop_frequency'],'C11_GPa':C[0,0],'C33_GPa':C[2,2],'C44_GPa':C[3,3],
     'E1_hole_xx_eV':be['hole']['E1_xx_eV'],'E1_hole_iso_eV':be['hole']['E1_iso_eV'],'E1_electron_xx_eV':be['electron']['E1_xx_eV'],'E1_electron_iso_eV':be['electron']['E1_iso_eV'],
     'gap_eV':s['bandgap']}
  for lv in (1e19,1e20):
   for c,sg in (('hole',1),('electron',-1)):p[f'm_c_inplane_{c}_{lv:.0e}'.replace('+','')]=mass_c(D,m,sg*lv)
  params.append(p)
 with open(O/'散射输入参数.csv','w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(params[0]));w.writeheader();w.writerows(params)
 for p in params:print({k:(round(v,4) if isinstance(v,float) else v) for k,v in p.items()})
 if not rows:print('no dense-grid transport.json under',O/'密网格AMSET','- mobility table and figures not written');return
 with open(O/'分项迁移率_全温度全浓度.csv','w',newline='',encoding='utf-8-sig') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 Md=[m for m in M if any(x['material']==m for x in rows)]
 plt.rcParams.update({'font.size':9,'axes.edgecolor':'#52514e','axes.labelcolor':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e'})
 # Figure 1: 300 K dot plot, rows = doping, columns = carrier; log mobility axis.
 fig,axs=plt.subplots(2,2,figsize=(12,6.4),layout='constrained',sharex=True,squeeze=False)
 for ri,lv in enumerate((1e19,1e20)):
  for ci,car in enumerate(('p','n')):
   ax=axs[ri,ci]
   for yi,m in enumerate(Md):
    r=next((x for x in rows if x['material']==m and x['carrier']==car and np.isclose(abs(x['doping_cm3']),lv,rtol=1e-3) and x['T_K']==300),None)
    if r is None:continue
    y=len(Md)-1-yi;vals=[r[f'mu_{k}_inplane'] for k in MECH+['overall']]
    ax.plot([min(vals),max(vals)],[y,y],color='#e1e0d9',lw=1,zorder=1)
    for k in MECH:ax.plot(r[f'mu_{k}_inplane'],y,MK[k],color=COL[k],ms=9,mec='white',mew=1.5,zorder=3,label=LABEL[k] if (ri,ci,yi)==(0,0,0) else None)
    ax.plot(r['mu_overall_inplane'],y,MK['overall'],color=COL['overall'],ms=16,mew=2.5,zorder=4,label=LABEL['overall'] if (ri,ci,yi)==(0,0,0) else None)
    ax.annotate(f"{r['mu_overall_inplane']:.0f}",(r['mu_overall_inplane'],y),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8,color='#0b0b0b')
   ax.set_yticks(range(len(Md)));ax.set_yticklabels(Md[::-1]);ax.set_xscale('log');ax.set_xlim(200,2e5);ax.set_ylim(-0.6,len(Md)-0.4)
   ax.grid(axis='x',which='major',alpha=.25,lw=.6)
   ax.set_title(f"{'holes (p)' if car=='p' else 'electrons (n)'}, |n| = {lv:.0e} cm$^{{-3}}$, 300 K".replace('e+',' × 10^').replace('1 × 10^19','10$^{19}$').replace('1 × 10^20','10$^{20}$'),loc='left',fontsize=9)
   if ri==1:ax.set_xlabel('in-plane mobility limited by each mechanism alone (cm$^2$ V$^{-1}$ s$^{-1}$, log scale)')
 fig.legend(*axs[0,0].get_legend_handles_labels(),loc='outside upper center',ncol=4,frameon=False,fontsize=8)
 fig.text(0.01,-0.01,CAP,fontsize=8,color='#52514e',ha='left',va='top')
 for ext in ('png','svg'):fig.savefig(O/f'分项迁移率_四材料_300K.{ext}',dpi=180,bbox_inches='tight')
 plt.close(fig)
 # Figure 2: temperature dependence at |n| = 1e20 (single-carrier-like at all T per the audit), small multiples.
 fig,axs=plt.subplots(2,len(Md),figsize=(3.5*len(Md),6),layout='constrained',sharey=True,squeeze=False)
 for ri,car in enumerate(('p','n')):
  for ci,m in enumerate(Md):
   ax=axs[ri,ci];rr=sorted([x for x in rows if x['material']==m and x['carrier']==car and np.isclose(abs(x['doping_cm3']),1e20,rtol=1e-3)],key=lambda x:x['T_K'])
   T=[x['T_K'] for x in rr]
   for k in MECH+['overall']:
    ax.plot(T,[x[f'mu_{k}_inplane'] for x in rr],color=COL[k],lw=2.4 if k=='overall' else 1.8,ls='-' if k!='overall' else '--',label=LABEL[k] if (ri,ci)==(0,0) else None)
   ax.set_yscale('log');ax.set_ylim(100,5e5);ax.grid(alpha=.25,lw=.6);ax.set_title(f"{m}, {'holes' if car=='p' else 'electrons'}",loc='left',fontsize=9)
   if ri==1:ax.set_xlabel('T (K)')
   if ci==0:ax.set_ylabel('in-plane mobility (cm$^2$ V$^{-1}$ s$^{-1}$)')
 fig.legend(*axs[0,0].get_legend_handles_labels(),loc='outside upper center',ncol=4,frameon=False,fontsize=8)
 fig.text(0.01,-0.01,'|n| = 1e20 cm^-3 (single-carrier-like at all T). '+CAP,fontsize=8,color='#52514e',ha='left',va='top')
 for ext in ('png','svg'):fig.savefig(O/f'分项迁移率_温度依赖_1e20.{ext}',dpi=180,bbox_inches='tight')
 plt.close(fig)
 for r in rows:
  if r['T_K']==300 and np.isclose(abs(r['doping_cm3']),1e19,rtol=1e-3) or (r['T_K']==300 and np.isclose(abs(r['doping_cm3']),1e20,rtol=1e-3)):
   print(r['material'],r['carrier'],f"{r['doping_cm3']:.0e}",*(f"{k}={r[f'mu_{k}_inplane']:.0f}" for k in ['overall']+MECH),f"zz/in={r['mu_overall_zz']/r['mu_overall_inplane']:.2f}")

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
