from pathlib import Path
import json,csv,hashlib,base64,html
import numpy as np
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
R=Path(sys.argv[1]).resolve(); O=R/'05_SOC_AMSET补算/06_SOC_无SOC对比';O.mkdir(exist_ok=True)
M=['Pb2Bi2Te5','Pb2Sb2Te5','Sn2Bi2Te5','Sn2Sb2Te5']
keys=['mobility','conductivity','seebeck','kappa_e','PF']
labels=['Mobility (cm²/V/s)','Conductivity (S/m)','Seebeck (µV/K)','Electronic κ (W/m/K)','PF (mW/m/K²)']
rows=[]; audit=[]; imgs=[]; summaries=[]
for m in M:
 paths=[R/'01_报告/源数据/AMSET_results/分类输运/00_合并原始结果'/m/'transport.json',R/'05_SOC_AMSET补算'/m/'04_AMSET输运/transport.json']
 data=[json.loads(p.read_text()) for p in paths]; old,new=data
 assert np.allclose(old['doping'],new['doping'],rtol=1e-10,atol=0)
 assert old['temperatures']==new['temperatures']
 nd,nt=len(new['doping']),len(new['temperatures']); arrays=[]
 for j in data:
  a={k:np.array(j[src] if k!='mobility' else j['mobility']['overall'],float) for k,src in [('mobility','mobility'),('conductivity','conductivity'),('seebeck','seebeck'),('kappa_e','electronic_thermal_conductivity')]}
  for v in a.values(): assert v.shape==(nd,nt,3,3) and np.isfinite(v).all()
  a={k:np.diagonal(v,axis1=-2,axis2=-1) for k,v in a.items()}
  a['PF']=a['conductivity']*a['seebeck']**2*1e-9
  arrays.append(a)
 settings=[yaml.safe_load(p.with_name('settings.yaml').read_text()) for p in paths]
 diff={k:[settings[0].get(k),settings[1].get(k)] for k in set(settings[0])|set(settings[1]) if settings[0].get(k)!=settings[1].get(k)}
 audit.append({'material':m,'inputs':[{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths],'settings_differences':diff,'points':nd*nt,'finite_and_shapes_verified':True})
 for di,d in enumerate(new['doping']):
  for ti,t in enumerate(new['temperatures']):
   for ai,axis in enumerate(['x','y','z','mean']):
    row=dict(material=m,temperature_K=t,doping_cm3=d,direction=axis)
    for k in keys:
     a,b=[v[k][di,ti,ai] if ai<3 else v[k][di,ti].mean() for v in arrays]
     row.update({k+'_noSOC':a,k+'_SOC':b,k+'_difference':b-a,k+'_relative_pct':100*(b-a)/abs(a) if abs(a)>1e-12 else ''})
    rows.append(row)
 ti=new['temperatures'].index(300)
 # Every standard concentration at 300K; two carrier signs and five properties.
 fig,axs=plt.subplots(2,5,figsize=(18,7),layout='constrained')
 for ri,sign in enumerate([-1,1]):
  ids=sorted([i for i,d in enumerate(new['doping']) if np.sign(d)==sign and any(np.isclose(abs(d),10.**e,rtol=1e-10) for e in range(17,21))],key=lambda i:abs(new['doping'][i]))
  for ci,k in enumerate(keys):
   ax=axs[ri,ci]
   for ar,style,name in zip(arrays,['s--','o-'],['no SOC','SOC']):ax.plot([abs(new['doping'][i]) for i in ids],[ar[k][i,ti].mean() for i in ids],style,label=name)
   ax.set_xscale('log');ax.set_xlabel('|net doping| (cm⁻³)');ax.set_ylabel(labels[ci]);ax.grid(alpha=.2)
   if k!='seebeck':ax.set_yscale('log')
   else:ax.axhline(0,color='gray',lw=.6)
   if ci==0:ax.set_title('Electron doping' if sign<0 else 'Hole doping');ax.legend()
 fig.suptitle(m+' | 300 K | directional arithmetic means')
 stem=m+'_300K_doping_comparison'
 for ext in ['png','svg']:fig.savefig(O/(stem+'.'+ext),dpi=180)
 plt.close(fig);imgs.append((stem,m+'：300 K默认正负浓度扫描。实线为SOC，虚线为无SOC；各量为三方向算术平均，PF先逐方向计算再平均。'))
 fig,axs=plt.subplots(2,5,figsize=(18,7),layout='constrained')
 for ri,sign in enumerate([-1,1]):
  di=next(i for i,d in enumerate(new['doping']) if np.isclose(d,sign*1e19,rtol=1e-10,atol=0))
  for ci,k in enumerate(keys):
   ax=axs[ri,ci]
   for ar,style,name in zip(arrays,['s--','o-'],['no SOC','SOC']):ax.plot(new['temperatures'],ar[k][di].mean(axis=-1),style,label=name,markersize=3)
   ax.set_xlabel('Temperature (K)');ax.set_ylabel(labels[ci]);ax.grid(alpha=.2)
   if k=='seebeck':ax.axhline(0,color='gray',lw=.6)
   if ci==0:ax.set_title(('n' if sign<0 else 'p')+' = 10¹⁹ cm⁻³');ax.legend()
 fig.suptitle(m+' | fixed net doping | temperature dependence')
 stem=m+'_temperature_comparison'
 for ext in ['png','svg']:fig.savefig(O/(stem+'.'+ext),dpi=180)
 plt.close(fig);imgs.append((stem,m+'：固定净浓度±10¹⁹ cm⁻³的全温区对比。Seebeck保留正负号；高温变号不等于缺陷类型改变。'))
 summary=[r for r in rows if r['material']==m and r['temperature_K']==300 and r['direction']=='mean' and np.isclose(abs(r['doping_cm3']),1e19,rtol=1e-10)]
 summaries+=summary
for name,rr in [('全温度_全浓度_方向及平均对比.csv',rows),('300K_正负1e19_汇总.csv',summaries)]:
 with (O/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rr[0]);w.writeheader();w.writerows(rr)
(O/'输入与校验.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
intro='所有输入清单路径以交付目录为根；两套方案共用 9×9×2 SOC DFPT 介电。四份新SOC AMSET计算已结束。本报告与原分类输运的合并原始结果逐点匹配，覆盖全部温度及浓度（包括缺陷净浓度）。总迁移率使用mobility.overall。σ、S、κe、μ的平均为迹/3；PF_i=σ_ii(S_ii×10⁻⁶)²×1000，单位mW/(m·K²)，PF平均为三个方向PF平均，不能用平均σ乘平均S平方。'
body='<h1>SOC与无SOC输运对比</h1><p>'+intro+'</p><p>本次同时更新SOC波函数和SOC形变势，并非只改变带隙；这是两套输运方案的对比，不能把全部差异单独归于孤对电子或某一个SOC散射效应。原缺陷浓度未重算，沿用同一净掺杂。网格及插值收敛尚未由本报告证明。</p>'
body+='<h2>300 K、净浓度±10¹⁹ cm⁻³</h2><table><tr><th>材料/类型</th>'+''.join('<th>'+k+' 无SOC → SOC</th>' for k in keys)+'</tr>'
for r in summaries:body+='<tr><td>'+r['material']+(' n' if r['doping_cm3']<0 else ' p')+'</td>'+''.join('<td>%.4g → %.4g</td>'%(r[k+'_noSOC'],r[k+'_SOC']) for k in keys)+'</tr>'
body+='</table><p>μ：cm²/(V·s)，σ：S/m，S：μV/K，κe：W/(m·K)，PF：mW/(m·K²)。差值与相对变化、x/y/z方向值及全部缺陷条件保存在完整CSV。</p>'
for stem,cap in imgs:body+='<figure><img src="data:image/png;base64,'+base64.b64encode((O/(stem+'.png')).read_bytes()).decode()+'"><figcaption>'+cap+'</figcaption></figure>'
body+='<h2>设置差异审计</h2><pre>'+html.escape(json.dumps(audit,ensure_ascii=False,indent=2))+'</pre>'
(O/'SOC_无SOC对比报告.html').write_text('<!doctype html><meta charset="utf-8"><style>body{max-width:1500px;margin:30px auto;font:16px/1.8 sans-serif;padding:25px}img{width:100%}td,th{border:1px solid #ccc;padding:8px}table{border-collapse:collapse}pre{white-space:pre-wrap}figure{margin:35px 0}</style>'+body)
print('Verified:',len(rows),'directional comparison rows;',len(imgs),'figures')
for r in summaries:print(r['material'],r['doping_cm3'],'mu',r['mobility_noSOC'],r['mobility_SOC'],'PF',r['PF_noSOC'],r['PF_SOC'])
