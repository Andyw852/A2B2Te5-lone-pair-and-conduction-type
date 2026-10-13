"""Parse SOC charge projections, deduplicate periodic k, and average edge windows.
Run: python band_windows.py DELIVERY [--project PROJECT] ; subsequent runs use gzip sources.
"""
from pathlib import Path
import argparse,gzip,shutil,json,re,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
def parse(p):
 opener=gzip.open if p.suffix=='.gz' else open
 with opener(p,'rt') as f:
  next(f); header=next(f); nk,nb,ni=map(int,re.findall(r':\s*(\d+)',header))
  E=np.full((nk,nb),np.nan);P=np.full((nk,nb,ni,9),np.nan);K=np.zeros((nk,3));KW=np.zeros(nk);occ=np.zeros((nk,nb));ki=bi=None
  for line in f:
   if line.lstrip().startswith('k-point'):
    mt=re.search(r'k-point\s+(\d+)\s*:\s*([\d .+Ee-]+)\s+weight\s*=\s*([\d.Ee+-]+)',line)
    ki=int(mt[1])-1;K[ki]=list(map(float,re.findall(r'[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?',mt[2])));KW[ki]=float(mt[3])
   elif line.startswith('band '):
    mt=re.search(r'band\s+(\d+)\s+# energy\s+([\d.Ee+-]+)\s+# occ.\s+([\d.Ee+-]+)',line)
    bi=int(mt[1])-1;E[ki,bi]=float(mt[2]);occ[ki,bi]=float(mt[3]);read=False
   elif line.lstrip().startswith('ion ') and not read:
    for j in range(ni):
     row=next(f).split();assert int(row[0])==j+1;P[ki,bi,j]=list(map(float,row[1:10]))
    read=True
 assert np.isfinite(E).all() and np.isfinite(P).all() and P.min()>=-1e-8
 return K,KW,E,P,occ

def run(D,project=None):
 src=D/'07_带边能窗/源数据';out=D/'07_带边能窗';src.mkdir(parents=True,exist_ok=True);figdir=D/'02_图表';figdir.mkdir(exist_ok=True)
 all_results={};flat=[];states=[]
 for m in M:
  folder=src/m;folder.mkdir(exist_ok=True);name='PROCAR_OPT' if 'Sb' in m else 'PROCAR'
  if project:
   origin=project/'band_and_elf/huang'/m/'step4_HSE_band'
   for fn in ['POSCAR','INCAR','KPOINTS','KPOINTS_OPT','kpath.json']:
    if (origin/fn).exists():shutil.copy2(origin/fn,folder/fn)
   if not (folder/(name+'.gz')).exists():
    with (origin/name).open('rb') as f,gzip.open(folder/(name+'.gz'),'wb',compresslevel=6) as g:shutil.copyfileobj(f,g)
   (folder/'原始来源.json').write_text(json.dumps({'source':str(origin/name),'sha256_uncompressed':hashlib.sha256((origin/name).read_bytes()).hexdigest()},ensure_ascii=False,indent=2))
  p=folder/(name+'.gz');K,KW,E,P,occ=parse(p)
  pos=(folder/'POSCAR').read_text().splitlines();elems=pos[5].split();counts=list(map(int,pos[6].split()));groups=[e+'-'+l for e in elems for l in ['s','p','d']]
  keep=np.arange(len(K)) if name=='PROCAR_OPT' else np.where(np.abs(KW)<1e-12)[0]
  assert len(keep)>100
  K,E,P=K[keep],E[keep],P[keep]
  grouped=[];start=0
  for n in counts:
   t=P[:,:,start:start+n].sum(axis=2);grouped.extend([t[:,:,0],t[:,:,1:4].sum(axis=2),t[:,:,4:].sum(axis=2)]);start+=n
  G=np.stack(grouped,axis=-1);nv=68 if 'Sb' in m else 88
  # Periodic duplicates are one k location, not additional valleys. Retain all distinct bands.
  keys=np.round(np.mod(K+1e-8,1),7);_,inverse=np.unique(keys,axis=0,return_inverse=True)
  raw_vbm=E[:,:nv].max();raw_cbm=E[:,nv:].min()
  eu=[];gu=[];ku=[];original=[];dup_de=0;dup_dp=0;dup_edge_de=0
  for i in range(max(inverse)+1):
   ids=np.where(inverse==i)[0];eu.append(E[ids].mean(axis=0));gu.append(G[ids].mean(axis=0));ku.append(K[ids[0]]);original.append((keep[ids]+1).tolist())
   dup_de=max(dup_de,float(np.ptp(E[ids],axis=0).max()));dup_dp=max(dup_dp,float(np.ptp(G[ids],axis=0).max()))
   near=((E[ids]>=raw_vbm-.101)&(E[ids]<=raw_cbm+.101)).any(axis=0)
   if near.any():dup_edge_de=max(dup_edge_de,float(np.ptp(E[ids][:,near],axis=0).max()))
  assert dup_edge_de<1e-4, 'Repeated k states differ in the analyzed energy window'
  E,G,K=np.array(eu),np.array(gu),np.array(ku)
  print(m,'duplicate max all-band spread',dup_de,flush=True)
  ev=float(E[:,:nv].max());ec=float(E[:,nv:].min());gap=ec-ev;direct=float((E[:,nv:].min(axis=1)-E[:,:nv].max(axis=1)).min())
  res=dict(gap_eV=gap,min_direct_gap_eV=direct,type='direct' if abs(direct-gap)<1e-5 else 'indirect',VBM_eV=ev,CBM_eV=ec,n_path_records=len(keep),n_unique_k=len(K),duplicate_energy_spread_eV=dup_de,duplicate_edge_energy_spread_eV=dup_edge_de,duplicate_group_projection_spread=dup_dp,groups=groups,windows=[],nv=nv)
  for edge,reference in [('VBM',ev),('CBM',ec)]:
   delta=reference-E if edge=='VBM' else E-reference
   valmask=np.arange(E.shape[1])<nv;side=valmask if edge=='VBM' else ~valmask
   for mev in [0.1,1,5,10,25,50,100]:
    ids=np.argwhere((delta>=-1e-9)&(delta<=mev/1000+1e-9)&side[None,:]);weights=G[ids[:,0],ids[:,1]];raw=weights.sum(axis=1)
    assert len(ids)>0 and raw.min()>0
    meanraw=weights.mean(axis=0);pct=100*meanraw/meanraw.sum();alt=(100*weights/raw[:,None]).mean(axis=0)
    row=dict(edge=edge,window_meV=mev,states=len(ids),unique_k=len(np.unique(ids[:,0])),raw_sum_mean=float(raw.mean()),raw_sum_min=float(raw.min()),raw_sum_max=float(raw.max()),percent=dict(zip(groups,pct.tolist())),state_normalized_percent=dict(zip(groups,alt.tolist())),normalization_difference_max_pp=float(np.max(abs(pct-alt))))
    res['windows'].append(row)
    flat.append(dict(material=m,edge=edge,window_meV=mev,states=len(ids),unique_k=row['unique_k'],raw_sum_mean=row['raw_sum_mean'],**row['percent']))
    for ki,bi in ids:
     states.append(dict(material=m,edge=edge,window_meV=mev,unique_k_id=int(ki),source_k_ids=';'.join(map(str,original[ki])),kx=K[ki,0],ky=K[ki,1],kz=K[ki,2],band=int(bi+1),energy_eV=E[ki,bi],distance_from_edge_meV=delta[ki,bi]*1000,raw_projection_sum=G[ki,bi].sum()))
  all_results[m]=res
  np.savez_compressed(folder/'解析投影.npz',K=K,E=E,G=G,groups=groups)
  print(m,'gap',gap,'k',len(keep),len(K),'duplicates',dup_de,flush=True)
  for r in res['windows']:
   if r['window_meV'] in [1,25,50,100]:print(r['edge'],r['window_meV'],r['states'],{k:round(v,2) for k,v in r['percent'].items() if k.endswith('-s') or k=='Te-p'},flush=True)
 (out/'带边能窗统计.json').write_text(json.dumps(all_results,ensure_ascii=False,indent=2))
 for name,rows in [('带边能窗统计.csv',flat),('入选态清单.csv',states)]:
  fields=list(dict.fromkeys(k for r in rows for k in r))
  with (out/name).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 plt.rcParams.update({'font.size':10,'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 colors=['#399875','#347eb1','#9675ae','#e5a14b','#ce6865','#dce0e5']
 fig,axs=plt.subplots(2,2,figsize=(12,10),layout='constrained')
 for ax,m in zip(axs.flat,M):
  a='Pb' if m.startswith('Pb') else 'Sn';b='Sb' if 'Sb' in m else 'Bi';keys=['Te-p',a+'-s',b+'-s',a+'-p',b+'-p'];rr=[r for r in all_results[m]['windows'] if r['window_meV'] in [1,25,50,100]]
  for y,r in enumerate(rr):
   vs=[r['percent'][k] for k in keys];vs.append(100-sum(vs));left=0
   for j,v in enumerate(vs):
    ax.barh(y,v,left=left,height=.67,color=colors[j],label=(keys+['Other'])[j] if y==0 else None)
    if v>=8:ax.text(left+v/2,y,f'{v:.1f}',ha='center',va='center',fontsize=8,color='white' if j<5 else '#243140')
    left+=v
   ax.text(101,y,f"{r['states']}/{r['unique_k']}",va='center',fontsize=8)
  ax.set_yticks(range(8),[r['edge']+f" {int(r['window_meV'])} meV" for r in rr]);ax.invert_yaxis();ax.set_xlim(0,118);ax.set_xticks([0,25,50,75,100]);ax.axhline(3.5,color='#777',lw=.6);ax.set_xlabel('Normalized charge projection (%)');ax.set_title(f"{m}  |  {all_results[m]['gap_eV']:.4f} eV",loc='left',fontweight='bold');ax.legend(ncol=3,loc='lower center',bbox_to_anchor=(.5,-.28),frameon=False,fontsize=9)
  ax.text(101,-.65,'states/k',fontsize=8)
 fig.suptitle('HSE + SOC band-edge windows | unique path k points',fontsize=16)
 for ext in ['png','pdf','svg']:fig.savefig(figdir/('带边_简并态与能窗.'+ext),dpi=300)
 plt.close(fig)
 return all_results
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('delivery',type=Path);p.add_argument('--project',type=Path);a=p.parse_args();run(a.delivery,a.project)
