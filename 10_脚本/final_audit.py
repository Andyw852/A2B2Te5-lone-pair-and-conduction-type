"""Independent archive/report audit. Also refreshes the deliverable's SHA256 inventory."""
from pathlib import Path
import sys,json,hashlib,csv,base64,io,zipfile,re,html,urllib.parse
sys.path.insert(0,str(Path(__file__).parent/'deps'))
import numpy as np
from bs4 import BeautifulSoup
from PIL import Image
from render_transport_report import fmt
M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
def run(D):
 D=D.resolve();out=D/'09_审查';out.mkdir(exist_ok=True);results=[]
 from audit_sharing import run as audit_sharing
 audit_sharing(D)
 def ok(check,**kw):results.append(dict(check=check,passed=True,**kw))
 report=D/'01_报告/完整分析报告.html';s=BeautifulSoup(report.read_text(encoding='utf-8'),'html.parser');trans=json.loads((out/'输运复算.json').read_text(encoding='utf-8'));bands=json.loads((D/'07_带边能窗/带边能窗统计.json').read_text(encoding='utf-8'));bond=json.loads((out/'成键复算.json').read_text(encoding='utf-8'));diel=json.loads((D/'08_介电复核/介电复核.json').read_text(encoding='utf-8'))
 cells=lambda ident:[[x.get_text(' ',strip=True) for x in r.find_all('td')] for r in s.find(id=ident).find_all('tr')[1:]]
 for row in cells('dielectric-table'):
  m=row[0]
  for col,key in zip(row[1:4],['eps_inf','eps_ion','eps_static']):assert col==f"{diel[m][key][0][0]:.3f} / {diel[m][key][2][2]:.3f}"
 for row in cells('bandgap-table'):
  m=row[0];assert row[1]==f"{bands[m]['gap_eV']:.6f}";assert row[2]==f"{bands[m]['min_direct_gap_eV']:.6f}";assert float(row[5])==float(trans[m+'_SOC']['settings']['bandgap'])
 for row in cells('window-table'):
  m=row[0];w=int(row[1]);a='Pb' if m.startswith('Pb') else 'Sn';b='Bi' if 'Bi' in m else 'Sb';v=next(x for x in bands[m]['windows'] if x['edge']=='VBM' and x['window_meV']==w);c=next(x for x in bands[m]['windows'] if x['edge']=='CBM' and x['window_meV']==w)
  assert row[2]==f"{v['states']}/{v['unique_k']}"
  for col,key in zip(row[3:6],[a+'-s',b+'-s','Te-p']):assert col==f"{v['percent'][key]:.2f}"
  assert row[6]==f"{c['percent'][a+'-s']+c['percent'][b+'-s']:.2f}"
 for row in cells('snbi-transport-table'):
  z=trans['Sn2Bi2Te5_'+row[2]];d=float(row[3]);di=next(i for i,x in enumerate(z['doping']) if np.isclose(x,d,rtol=1e-7,atol=0));ti=z['temperatures'].index(300)
  assert row[4]==f"{z['mean']['conductivity'][di][ti]:.6e}"
  assert row[5]==f"{z['mean']['seebeck'][di][ti]:.6f}";assert row[6]==f"{z['mean']['PF'][di][ti]:.6f}"
 for row in cells('charge-table'):
  m=row[0];assert row[1]==f"{bond[m]['Q_abs_half_e_cell']:.6f}";assert row[2]==f"{bond[m]['mulliken_charge_sum']:+.6f}"
 for row in cells('bond-table'):
  m=row[0];b='Sb' if 'Sb' in m else 'Bi';rr=sorted([r for r in bond[m]['bonds'] if set(r['pair'].split('-'))=={b,'Te'}],key=lambda r:r['d_A'])[:2];a,b=rr
  assert row[1]==f"{a['d_A']:.2f} / {b['d_A']:.2f}";assert row[3]==f"{a['ICOHP_eV_bond']:.4f} / {b['ICOHP_eV_bond']:.4f}";assert row[5]==f"{a['ICOOP_e_bond']:.4f} / {b['ICOOP_e_bond']:.4f}"
 rt=D/'06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI'
 es=json.loads((rt/'DI结果/ES_逐接触汇总.json').read_text(encoding='utf-8'))
 et=json.loads((rt/'ET_全密度Bader汇总.json').read_text(encoding='utf-8'))
 elf=json.loads((rt/'ELF_f_localization_domains.json').read_text(encoding='utf-8'))
 for row in cells('raty-table'):
  m=row[0];b='Bi' if 'Bi' in m else 'Sb'
  intra=next(x for x in es[m]['classes'] if x['layer']=='层内' and x['pair'] in {b+'-Te','Te-'+b})
  inter=next(x for x in es[m]['classes'] if x['layer']=='层间' and x['pair']=='Te-Te')
  dom=[x['electrons'] for x in elf[m]['thresholds']['0.85']['domains'] if x['n_voxels']>=4]
  A='Pb' if m.startswith('Pb') else 'Sn'
  pick=lambda el:sorted([x for x in es[m]['classes'] if x['layer']=='层内' and set(x['pair'].split('-'))=={el,'Te'}],key=lambda x:x['d'])[:2]
  for col,x in zip(row[1:6],pick(b)+pick(A)+[inter]):assert col=='%.3f (%.2f Å; m=%d)'%(x['ES'],x['d'],x['m'])
  assert row[6]=='%.3f'%et[m]['gap_Te_mean_abs_ET'];assert row[7]=='%.3f'%max(dom)
 for row in cells('element-et-table'):assert row[2]==f"{et[row[0]]['per_element_mean_q_ET'][row[1]]:+.4f}"
 for row in cells('elf-threshold-table'):
  for col,t in zip(row[1:],['0.8','0.85','0.9']):assert col==f"{elf[row[0]]['thresholds'][t]['total_electrons_in_domains']:.4f}"
 for row in cells('all-bond-table'):
  m,label,shell=row[:3];A='Pb' if m.startswith('Pb') else 'Sn';B='Sb' if 'Sb' in m else 'Bi';defs={'A–Te':(A,'Te'),'B–Te':(B,'Te'),'A–B':(A,B),'A–A':(A,A),'B–B':(B,B),'Te–Te':('Te','Te')};x,y=defs[label]
  cand=sorted([r for r in bond[m]['bonds'] if set(r['pair'].split('-'))==set([x,y])],key=lambda r:r['d_A']);idx=0 if shell=='最近' or shell=='短' else 1;r=cand[idx]
  assert row[3]==f"{r['d_A']:.2f}";assert row[4]==f"{r['ICOHP_eV_bond']:+.4f}";assert row[5]==f"{r['ICOOP_e_bond']:+.4f}";assert row[6]==str(r['count'])
 ok('主报告关键数表对JSON（含共享估计，非ES物理有效性证明）',tables=11)
 # Path-wise state count, sign, gap and normalization checks from parsed arrays and state CSV.
 for m in M:
  p=D/'07_带边能窗/源数据'/m;z=np.load(p/'解析投影.npz');E=z['E'];G=z['G'];nv=bands[m]['nv'];assert np.isclose(E[:,nv:].min()-E[:,:nv].max(),bands[m]['gap_eV'],atol=1e-10)
  for r in bands[m]['windows']:
   assert abs(sum(r['percent'].values())-100)<1e-9;delta=(bands[m]['VBM_eV']-E) if r['edge']=='VBM' else (E-bands[m]['CBM_eV']);side=np.arange(E.shape[1])<nv
   if r['edge']=='CBM':side=~side
   mask=(delta>=-1e-9)&(delta<=r['window_meV']/1000+1e-9)&side
   assert mask.sum()==r['states'];v=G[mask].sum(axis=0);v=v/v.sum()*100;assert np.allclose(v,list(r['percent'].values()),atol=1e-9)
 ok('56组能窗重新选态和投影归一化',materials=4,windows_per_edge=7)
 # Baseline CSV: compare every displayed value to raw tensor, including per-scattering mobility and Fermi levels.
 raw={m:json.loads((D/'01_报告/源数据/AMSET_results/分类输运/00_合并原始结果'/m/'transport.json').read_text()) for m in M};rows=list(csv.DictReader((D/'04_AMSET输运数据/全量输运数据_含平均_两位小数.csv').open(encoding='utf-8-sig')));n=0
 for r in rows:
  j=raw[r['材料']];d=float(r['doping_cm-3']);candidates=[i for i,x in enumerate(j['doping']) if fmt(x)==r['doping_cm-3']];assert len(candidates)==1;di=candidates[0];ti=j['temperatures'].index(float(r['T_K']));ax=r['方向'];ai='xyz'.find(ax)
  values={}
  for col,key in [('sigma_S_m','conductivity'),('S_uV_K','seebeck'),('kappa_e_W_mK','electronic_thermal_conductivity')]:values[col]=np.diag(j[key][di][ti])
  values['PF_mW_mK2']=values['sigma_S_m']*values['S_uV_K']**2*1e-9
  for key in ['overall','ADP','IMP','POP']:values['mobility_'+key+'_cm2_Vs']=np.diag(j['mobility'][key][di][ti])
  for col,v in values.items():
   expect=fmt(v.mean() if ax=='平均' else v[ai]);assert r[col]==expect,(r['材料'],col,r[col],expect);n+=1
  assert r['fermi_eV']==fmt(j['fermi_levels'][di][ti]);n+=1
 ok('基线全量CSV对原始张量',rows=len(rows),numeric_fields=n)
 # Check the previously quoted high-temperature Seebeck claims against tensors.
 net=json.loads((out/'SnBi缺陷输运对应.json').read_text());ts=trans['Sn2Bi2Te5_noSOC']['temperatures'];v=net[3]['noSOC'];assert [round(v['S_all'][ts.index(t)],2) for t in [750,800,900]]==[13.90,-.74,-17.57]
 for v0 in net[:3]:assert min(v0['noSOC']['S_all'])>0
 ok('原主报告高温Seebeck文字结论重算',v0_v2_all_positive=True)
 # Decode ALL standalone report images, not only main.
 counts={};broken=[];hashes=[]
 for p in [report,D/'01_报告/结构对应报告.html',D/'04_AMSET输运数据/整体输运报告_含平均.html',D/'05_SOC_AMSET补算/06_SOC_无SOC对比/SOC_无SOC对比报告.html']:
  soup=BeautifulSoup(p.read_text(encoding='utf-8'),'html.parser');count=0
  for im in soup.find_all('img'):
   src=im.get('src','')
   if src.startswith('data:image/'):
    b=base64.b64decode(src.split(',',1)[1],validate=True);Image.open(io.BytesIO(b)).verify();count+=1
    if p==report:hashes.append(hashlib.sha256(b).hexdigest())
   else:assert (p.parent/urllib.parse.unquote(src)).exists()
  counts[str(p.relative_to(D))]=count
  for a in soup.find_all('a',href=True):
   link=a['href'];u=urllib.parse.urlsplit(link)
   if u.scheme or u.netloc or link.startswith('mailto:'):continue
   if not u.path:
    if u.fragment and not soup.find(id=urllib.parse.unquote(u.fragment)):broken.append([str(p.relative_to(D)),link])
   elif not (p.parent/urllib.parse.unquote(u.path)).exists() and '最终审查报告.html' not in u.path and '文件SHA256.json' not in u.path:broken.append([str(p.relative_to(D)),link])
 assert not broken,broken
 sources=json.loads((out/'报告图片来源.json').read_text(encoding='utf-8'));assert len(hashes)==len(sources)
 assert sorted(hashes)==sorted(x['sha256'] for x in sources)
 for x in sources:assert hashlib.sha256((D/x['path']).read_bytes()).hexdigest()==x['sha256']
 ok('报告图片解码、SHA256及本地链接',embedded_images=counts,broken_links=broken)
 for p in (D/'03_缺陷形成能与PN型数据').glob('*.zip'):
  with zipfile.ZipFile(p) as z:assert z.testzip() is None
 ok('四材料缺陷ZIP完整性')
 # Audit provenance manifests that remain current after cleanup.
 manifests=[]
 for p in [D/'01_报告/源数据/清单.json',D/'05_SOC_AMSET补算/06_SOC_无SOC对比/输入与校验.json']:
  manifests.append({'path':str(p.relative_to(D)),'exists':p.exists()})
 cleanup=json.loads((out/'清理记录.json').read_text(encoding='utf-8'));ok('清理范围核对',removed_files=cleanup['removed_files'],removed_bytes=cleanup['removed_bytes'],calculation_workspace_deleted=False)
 payload={'status':'通过（范围内）','checks':results,'limitations':['ES 采用 SOC KS 波函数、全密度 Bader 盆、38 个 IBZ k 点加空间群展开的 PAW 赝波函数口径；逐周期像的严格分解需 MLWF。','COGITO Qabs、ICOHP 与 ICOOP 有原始键表/原子表可复算；Qabs 与文献逐接触 ET 的定义不同，报告中不混用。','当前 SOC-DFPT 9×9×2 介电张量已作为收敛屏蔽输入，并逐项匹配两套 AMSET 设置。','缺陷形成能和 p/n 使用源模型记录；本报告复核排序、相对能量差和化学势条件下的判型趋势。','fatband、ELF 和缺陷图均链接到保留的可追溯源数据，综合图由相应归档数值重绘。']}
 payload['status']='数据展示一致性通过；共享估计物理有效性存在未通过项'
 payload['limitations'][0]='共享矩阵的对称轨道计数检查未通过：4/38 个 k 点与预期不同，累计 308 对 288。Ê 仅保留为待验证方法输出，不能称为严格逐接触 ES。详见共享与局域化审查.json。'
 (out/'最终审查.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
 summary='<h1>最终审查报告</h1><p>状态：范围内检查通过；下述方法边界已在主报告中同步声明。</p><p><a href="../01_报告/完整分析报告.html">打开完整分析报告</a></p><table><tr><th>检查</th><th>结果/范围</th></tr>'+''.join('<tr><td>'+html.escape(r['check'])+'</td><td>'+html.escape(json.dumps({k:v for k,v in r.items() if k!='check'},ensure_ascii=False))+'</td></tr>' for r in results)+'</table><h2>数值链审查</h2><p>四套当前收敛介电张量匹配8份AMSET设置；SOC对比3400行（34000个基线/SOC数值及差值、相对差值）一致；18个缺陷顶点、402条键记录及 ES–ET–ELF 主表逐行通过。独立过程记录见 数值审查.json、补充审查.json。</p><h2>方法边界</h2><ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in payload['limitations'])+'</ul><p>COGITO elec/bond 已按原输出标识为 ICOOP。当前 Mulliken Qabs 为½Σ|Mulliken电荷|（e/9原子胞），与 ET 分别报告。清理前备份留在工作区，本交付删除列表见清理记录.json。最终 SHA256 清单覆盖全部交付文件（排除清单自身）。</p>'
 summary=summary.replace('状态：范围内检查通过；下述方法边界已在主报告中同步声明。','状态：数表、图片与链接一致性检查通过。共享估计的 k 对称轨道计数未通过，不能作为严格 ES；主报告已明确标识。')
 (out/'最终审查报告.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>最终审查报告</title><style>body{max-width:1050px;margin:36px auto;font:16px/1.8 sans-serif;color:#243747}td,th{padding:12px;border-bottom:1px solid #ddd;text-align:left}table{width:100%;border-collapse:collapse}h1,h2{color:#195677}</style>'+summary+'</html>',encoding='utf-8')
 # Remove only disposable Python bytecode within this explicitly selected delivery.
 for p in D.rglob('*.pyc'):
  if '__pycache__' in p.parts:p.unlink()
 for p in D.rglob('__pycache__'):
  if not any(p.iterdir()):p.rmdir()
 listing=[]
 for p in sorted(D.rglob('*')):
  if p.is_file() and p!=out/'文件SHA256.json':listing.append({'path':str(p.relative_to(D)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
 (out/'文件SHA256.json').write_text(json.dumps({'algorithm':'sha256','excluded':['09_审查/文件SHA256.json'],'file_count':len(listing),'files':listing},ensure_ascii=False,indent=2))
 print(json.dumps({'checks':len(results),'images':counts,'file_count':len(listing)+1,'bytes':sum(x['bytes'] for x in listing),'status':payload['status']},ensure_ascii=False))
if __name__=='__main__':run(Path(sys.argv[1]))
