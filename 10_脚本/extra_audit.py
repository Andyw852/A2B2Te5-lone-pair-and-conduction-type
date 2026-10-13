from pathlib import Path
import sys,json,re,zipfile,csv,gzip,base64
sys.path.insert(0,str(Path(__file__).parent/'deps'))
import numpy as np
from bs4 import BeautifulSoup
from analyze_sources import dump,M
def run(D):
 old=BeautifulSoup((D/'01_报告/完整分析报告.html').read_text(),'html.parser');defects={};checks=[]
 for m,idx in [('Pb2Sb2Te5',1),('Sn2Sb2Te5',2),('Pb2Bi2Te5',3),('Sn2Bi2Te5',4)]:
  with zipfile.ZipFile(D/'03_缺陷形成能与PN型数据'/f'{m}.zip') as z:j=json.loads(z.read('计算结果/formation_energy_results.json'))
  vertices=[];trs=old.find(id=f'm{idx}').find('table').find_all('tr')[1:]
  for v,tr in zip(j['vertex_summary'],trs,strict=True):
   rank=[];ref=j['vertex_summary'][0]['mu']
   for did,d in j['defects'].items():
    shift=sum(n*(ref[e]-v['mu'][e]) for e,n in d['dn'].items());en,q=min((en+shift+int(q)*v['E_F_eq'],int(q)) for q,en in d['Ef0'].items());rank.append((en,did,q))
   rank.sort();cells=[x.get_text(' ',strip=True) for x in tr.find_all('td')]
   assert rank[0][1].startswith('def-'+re.search(r'\[(\d+)\]',cells[1])[1]);assert rank[1][1].startswith('def-'+re.search(r'\[(\d+)\]',cells[3])[1]);assert abs(rank[0][0]-float(cells[2]))<1.1e-6;assert abs((rank[1][0]-rank[0][0])*1000-float(cells[4]))<.0011
   assert rank[0][2]==int(cells[1].split('q=')[-1]);assert np.isclose(float(cells[6]),v['p0_cm3']/v['n0_cm3'],rtol=5e-4)
   vertices.append(dict(vertex=v['vertex'],minimum=rank[0],second=rank[1],delta_meV=(rank[1][0]-rank[0][0])*1000,mu=v['mu'],EF_eq=v['E_F_eq'],n=v['n0_cm3'],p=v['p0_cm3'],net_p_minus_n=v['p0_cm3']-v['n0_cm3']))
  defects[m]=vertices
 checks.append(dict(check='defect_rankings_and_display',vertices=18,pass_check=True))
 # Main SnBi table: precise input doping vs old defect net values, no nearest arbitrary point.
 T=json.loads((D/'09_审查/输运复算.json').read_text());sn=[]
 for v in defects['Sn2Bi2Te5']:
  entry=dict(vertex=v['vertex'],net_defect=v['net_p_minus_n'])
  for variant in ['noSOC','SOC']:
   t=T['Sn2Bi2Te5_'+variant];di=next(i for i,d in enumerate(t['doping']) if np.isclose(d,v['net_p_minus_n'],rtol=1e-5,atol=0));ti=t['temperatures'].index(300)
   entry[variant]=dict(input_doping=t['doping'][di],doping_rounding_relative=t['doping'][di]/v['net_p_minus_n']-1,S300=t['mean']['seebeck'][di][ti],sigma300=t['mean']['conductivity'][di][ti],PF300=t['mean']['PF'][di][ti],S_all=t['mean']['seebeck'][di],sigma_all=t['mean']['conductivity'][di],PF_all=t['mean']['PF'][di],diagonal_S=t['diagonal']['seebeck'][di])
  sn.append(entry)
 checks.append(dict(check='defect_net_doping_matched_both_transport',rows=4,pass_check=True))
 # Raw source totals vs published CSV; then redraw total COHP from raw traces.
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axs=plt.subplots(1,4,figsize=(13,6),sharey=True,layout='constrained')
 for ax,m in zip(axs,M):
  p=D/'06_COGITO成键分析与COHP/图_交互式COHP'/f'{m}_COHP交互图.html';text=p.read_text();mt=re.search(r'Plotly\.newPlot\(\s*"plotDiv",\s*',text);traces,_=json.JSONDecoder().raw_decode(text,mt.end());tr=next(x for x in traces if x.get('name','').strip()=='total COHP')
  def arr(v):return np.frombuffer(base64.b64decode(v['bdata']),dtype=v['dtype']) if isinstance(v,dict) else np.array(v)
  x,y=arr(tr['x']),arr(tr['y']);ax.plot(x,y,color='#485869',lw=.8);ax.fill_betweenx(y,0,x,where=x<=0,interpolate=True,color='#4c9dbc',alpha=.8);ax.fill_betweenx(y,0,x,where=x>=0,interpolate=True,color='#cc685c',alpha=.8);ax.axhline(0,color='#bd403d',ls='--',lw=.8);ax.axvline(0,color='#777',lw=.6);ax.set_title(m);ax.set_xlim(-7,30);ax.set_ylim(-24.5,6);ax.set_xlabel('COHP (source scale)')
  dest=D/'06_COGITO成键分析与COHP/审核与重绘'/f'{m}_total_COHP.csv'
  np.savetxt(dest,np.column_stack([y,x]),delimiter=',',header='E_minus_EF_eV,COHP_source_scale',comments='')
  checks.append(dict(check='total_COHP_extracted_from_raw_HTML',material=m,points=len(x),pass_check=True))
 axs[0].set_ylabel(r'$E-E_F$ (eV)');fig.suptitle('Occupied bonding and antibonding contributions | total COHP')
 for ext in ['png','pdf','svg']:fig.savefig(D/'02_图表'/('COHP_四材料总曲线_审核重绘.'+ext),dpi=300)
 plt.close(fig)
 # Compare denser grid from original OUTCAR; no extrapolation to other compounds.
 def tensors(folder):
  t=gzip.open(folder/'OUTCAR.gz','rt').read();out={}
  for key,tag in [('eps_inf','MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)'),('eps_ion','MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION')]:
   i=t.rfind(tag);assert i>=0;rows=[]
   for line in t[i:].splitlines()[1:]:
    w=line.split()
    if len(w)==3:
     try:rows.append(list(map(float,w)))
     except ValueError:continue
    if len(rows)==3:break
   out[key]=rows
  out['eps_static']=(np.array(out['eps_inf'])+out['eps_ion']).tolist();out['finished']='General timing and accounting' in t;return out
 dump(D/'09_审查/缺陷复算.json',defects);dump(D/'09_审查/SnBi缺陷输运对应.json',sn);dump(D/'09_审查/补充审查.json',checks)
 print('SnBi',[(x['vertex'],x['noSOC']['S300'],x['SOC']['S300'],x['noSOC']['PF300'],x['SOC']['PF300']) for x in sn])
if __name__=='__main__':run(Path(sys.argv[1]))
