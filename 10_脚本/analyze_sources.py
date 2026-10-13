from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).parent/'deps'))
import json,csv,hashlib,gzip,re,zipfile,argparse
import numpy as np
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
def dump(p,x):p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(D,project=None):
 out=D/'09_审查';out.mkdir(exist_ok=True);F=D/'02_图表';checks=[];bond={};diel={};trans={}
 for m in M:
  folder=D/'06_COGITO成键分析与COHP/结果'/m;atoms=json.loads((folder/'all_atoms.json').read_text());unique=json.loads((folder/'all_unique_bonds.json').read_text())
  entries=[]
  for line in (folder/'bond_info.txt').read_text().splitlines()[2:]:
   w=line.split()
   if len(w)!=7 or not w[0].isalpha() or not w[1].isalpha():continue
   a,b=w[:2]
   try:v=list(map(float,w[2:]))
   except ValueError:continue
   d,cohp,coop,n,total=v
   candidates=[x for k,x in unique.items() if sorted(k.split()[:2])==sorted([a,b]) and abs(x['bond length']-d)<.00501 and abs(x['cohp']-cohp)<.000051 and abs(x['coop']-coop)<.000051 and x['degeneracy']==n]
   assert candidates,(m,line)
   entries.append(dict(pair=a+'-'+b,d_A=d,ICOHP_eV_bond=cohp,ICOOP_e_bond=coop,count=n,ICOHP_eV_cell=total))
  q=[v['mulli_partial_chg'] for v in atoms.values()];sums={key:sum(sum(sum(z) for z in v[key].values()) for v in atoms.values()) for key in ['onsite_occup','mulli_occup']}
  meanq={el:float(np.mean([v['mulli_partial_chg'] for v in atoms.values() if v['elem']==el])) for el in sorted(set(v['elem'] for v in atoms.values()))}
  bond[m]=dict(bonds=entries,electron_population_sums=sums,mulliken_charge_sum=sum(q),Q_abs_half_e_cell=sum(abs(x) for x in q)/2,element_mean_mulliken_q=meanq,quality=json.loads((folder/'cogito_summary.json').read_text())['quality'])
  checks.append(dict(check='bond_table_to_unique_cohp_coop',material=m,rows=len(entries),pass_check=True))
  assert abs(sum(q))<.001
  files=[]
  for variant,p in [('noSOC',D/'01_报告/源数据/AMSET_results/分类输运/00_合并原始结果'/m/'transport.json'),('SOC',D/'05_SOC_AMSET补算'/m/'04_AMSET输运/transport.json')]:
   j=json.loads(p.read_text());s=yaml.safe_load(p.with_name('settings.yaml').read_text());shape=(len(j['doping']),len(j['temperatures']),3,3)
   ar={k:np.array(j[k],float) for k in ['conductivity','seebeck','electronic_thermal_conductivity']};ar['mobility']=np.array(j['mobility']['overall'])
   for a in ar.values():assert a.shape==shape and np.isfinite(a).all()
   diag={k:np.diagonal(v,axis1=-2,axis2=-1) for k,v in ar.items()};diag['PF']=diag['conductivity']*diag['seebeck']**2*1e-9
   if project:
    orig=project/('AMSET_results/分类输运/00_合并原始结果' if variant=='noSOC' else 'SOC_AMSET_补算_结果汇总')/m
    if variant=='SOC':orig=orig/'04_AMSET输运'
    for fn in ['transport.json','settings.yaml']:assert sha(orig/fn)==sha(p.with_name(fn)),str(p)
   for mir in D.rglob('transport.json'):
    if m not in str(mir):continue
    is_soc='SOC_AMSET_补算_结果汇总' in str(mir) or (mir.parent.name=='04_AMSET输运' and '源数据' not in str(mir))
    if is_soc==(variant=='SOC'):assert sha(mir)==sha(p),str(mir)
   trans[m+'_'+variant]=dict(source=str(p.relative_to(D)),doping=j['doping'],temperatures=j['temperatures'],settings=s,mean={k:v.mean(axis=-1).tolist() for k,v in diag.items()},diagonal={k:v.tolist() for k,v in diag.items()},sha256=sha(p))
   checks.append(dict(check='transport_finite_shape_and_mirrors',material=m,variant=variant,shape=shape,pass_check=True))
  folder=D/'08_介电复核/源数据'/m
  if (folder/'OUTCAR.gz').exists():
   txt=gzip.open(folder/'OUTCAR.gz','rt',errors='replace').read()
   def grab(tag):
    i=txt.rfind(tag);assert i>=0;rows=[]
    for l in txt[i:].splitlines()[1:]:
     w=l.split()
     if len(w)==3:
      try:rows.append(list(map(float,w)))
      except ValueError:continue
     if len(rows)==3:return np.array(rows)
    raise ValueError(tag)
   inf=grab('MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)');ion=grab('MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION');sta=inf+ion
   assert 'General timing and accounting' in txt
   s=trans[m+'_SOC']['settings'];s0=trans[m+'_noSOC']['settings']
   for key,ar in [('high_frequency_dielectric',inf),('static_dielectric',sta)]:
    assert np.max(abs(ar-np.array(s[key],float)))<1e-5,(m,key,ar,s[key]);assert np.allclose(np.array(s[key],float),np.array(s0[key],float),rtol=0,atol=1e-10)
   diel[m]=dict(eps_inf=inf.tolist(),eps_ion=ion.tolist(),eps_static=sta.tolist(),nkpts=int(re.search(r'NKPTS\s*=\s*(\d+)',txt)[1]),kpoints=(folder/'KPOINTS').read_text(),finished=True,pop_frequency=s.get('pop_frequency'),pop_same_both=s.get('pop_frequency')==s0.get('pop_frequency'),sha256_outcar=hashlib.sha256(txt.encode()).hexdigest())
   checks.append(dict(check='DFPT_OUTCAR_to_both_AMSET_settings',material=m,pass_check=True))
 # Independent recalculation of both comparison CSVs, including all directions and properties.
 for name in ['全温度_全浓度_方向及平均对比.csv','300K_正负1e19_汇总.csv']:
  path=D/'05_SOC_AMSET补算/06_SOC_无SOC对比'/name;rows=list(csv.DictReader(path.open(encoding='utf-8-sig')));nval=0;maxerr=0
  for r in rows:
   m=r['material'];d=float(r['doping_cm3']);t=float(r['temperature_K']);axis=r['direction']
   for key in ['mobility','conductivity','seebeck','kappa_e','PF']:
    k='electronic_thermal_conductivity' if key=='kappa_e' else key;vals=[]
    for variant in ['noSOC','SOC']:
     a=trans[m+'_'+variant];di=next(i for i,x in enumerate(a['doping']) if np.isclose(x,d,rtol=1e-10,atol=0));ti=a['temperatures'].index(t)
     v=a['mean'][k][di][ti] if axis=='mean' else a['diagonal'][k][di][ti]['xyz'.index(axis)];vals.append(v)
     assert np.isclose(float(r[key+'_'+variant]),v,rtol=1e-10,atol=1e-10);maxerr=max(maxerr,abs(float(r[key+'_'+variant])-v));nval+=1
    diff=vals[1]-vals[0];assert np.isclose(float(r[key+'_difference']),diff,rtol=1e-10,atol=1e-10)
    if r[key+'_relative_pct']:assert np.isclose(float(r[key+'_relative_pct']),100*diff/abs(vals[0]),rtol=1e-10,atol=1e-10)
  checks.append(dict(check='recomputed_SOC_comparison_CSV',file=name,rows=len(rows),checked_values=nval,max_abs_error=maxerr,pass_check=True))
 dump(out/'数值审查.json',checks);dump(out/'输运复算.json',trans);dump(out/'成键复算.json',bond);dump(D/'08_介电复核/介电复核.json',diel)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained');x=np.arange(4)
 for ax,key,title in [(axs[0,0],'ICOHP_eV_bond','(a) Bond-energy contrast'),(axs[0,1],'ICOOP_e_bond','(b) Overlap populations (not ICOBI)')]:
  for j,(lab,col) in enumerate([('Short B-Te','#257e9a'),('Long B-Te','#d39a4c')]):
   vals=[]
   for m in M:
    b='Sb' if 'Sb' in m else 'Bi';rows=sorted([r for r in bond[m]['bonds'] if sorted(r['pair'].split('-'))==sorted(['Te',b])],key=lambda r:r['d_A']);vals.append(rows[j][key]*(-1 if key.startswith('ICOHP') else 1))
   ax.bar(x+(j-.5)*.32,vals,width=.32,label=lab,color=col)
   for xi,v in zip(x+(j-.5)*.32,vals):ax.text(xi,v,f'{v:.3f}',ha='center',va='bottom',fontsize=8)
  ax.set_title(title,loc='left');ax.set_ylabel('-ICOHP (eV/bond)' if key.startswith('ICOHP') else 'ICOOP (e/bond)');ax.set_ylim(top=ax.get_ylim()[1]*1.22);ax.legend(frameon=False,ncol=2,loc='upper center')
 ax=axs[1,0]
 for j,(site,col) in enumerate([('A','#267ab0'),('B','#9a6eaa'),('Te','#359873')]):
  vals=[bond[m]['element_mean_mulliken_q'][('Pb' if m.startswith('Pb') else 'Sn') if site=='A' else ('Sb' if 'Sb' in m else 'Bi') if site=='B' else 'Te'] for m in M]
  ax.bar(x+(j-1)*.25,vals,width=.25,label=site,color=col)
 ax.axhline(0,color='#777',lw=.7);ax.set_ylabel('Mean Mulliken charge (e/atom)');ax.set_title('(c) Atomic charge partition',loc='left');ax.legend(ncol=3,frameon=False)
 ax=axs[1,1];vs=[bond[m]['Q_abs_half_e_cell'] for m in M];ax.bar(x,vs,color='#447a96')
 for xi,v in zip(x,vs):ax.text(xi,v,f'{v:.3f}',ha='center',va='bottom')
 ax.set_ylabel(r'$\frac{1}{2}\sum_i |q_i^{\rm Mulliken}|$ (e / 9-atom cell)');ax.set_title('(d) Defined cell charge-transfer descriptor',loc='left')
 for ax in axs.flat:ax.set_xticks(x,M,rotation=15,ha='right');ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
 fig.suptitle('Bond disproportionation and charge redistribution | PBE, no SOC',fontsize=14)
 for ext in ['png','pdf','svg']:fig.savefig(F/('成键_ICOHP_ICOOP_电荷转移.'+ext),dpi=300)
 plt.close(fig)
 print(json.dumps(checks,ensure_ascii=False,indent=1));print('bond charge',[(m,bond[m]['Q_abs_half_e_cell']) for m in M]);print('dielectric',diel)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('delivery',type=Path);p.add_argument('--project',type=Path);a=p.parse_args();run(a.delivery,a.project)
