"""Independent arithmetic audit; does not certify the DI estimator as physical ES."""
from pathlib import Path
import json,sys,gzip,hashlib,shutil
import numpy as np

def run(D):
 rt=D/'06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI'
 es=json.loads((rt/'DI结果/ES_逐接触汇总.json').read_text())
 matrices=json.loads((rt/'DI结果/DI_results_lowdin.json').read_text())
 et=json.loads((rt/'ET_全密度Bader汇总.json').read_text())
 elf=json.loads((rt/'ELF_f_localization_domains.json').read_text())
 checks={}
 for m,v in es.items():
  d=np.array(matrices[m]['DI']);assert np.isfinite(d).all();assert np.allclose(d,d.T,atol=1e-10)
  for p in v['pairs']:
   assert np.isclose(p['DI_folded'],d[p['a'],p['b']],atol=1e-12)
   assert np.isclose(p['ES_contact'],2*p['DI_folded']/p['m'],atol=1e-12)
  for c in v['classes']:
   ps=[p for p in v['pairs'] if p['pair']==c['pair'] and p['layer']==c['layer'] and float(np.round(p['d'],2))==c['d']]
   assert len(ps)==c['n_pairs'];assert np.isclose(np.mean([p['ES_contact'] for p in ps]),c['ES'],atol=1e-12)
  atoms=et[m]['atoms']
  for a in atoms:assert abs(a['zval']-a['N_basin_valence']-a['q_ET'])<2e-5
  gap=[a for a in atoms if a['index'] in et[m]['gap_Te_indices']]
  assert abs(np.mean([abs(a['q_ET']) for a in gap])-et[m]['gap_Te_mean_abs_ET'])<6e-5
  for el,q in et[m]['per_element_mean_q_ET'].items():assert abs(np.mean([a['q_ET'] for a in atoms if a['element']==el])-q)<6e-5
  for t,vv in elf[m]['thresholds'].items():assert abs(sum(x['electrons'] for x in vv['domains'])-vv['total_electrons_in_domains'])<3e-4
  c=v['checks'];err=np.max(np.abs(d.sum(axis=1)-2*np.array(c['NA']))/(2*np.array(c['NA'])))
  assert np.isclose(err,c['sum_rule_max_rel_err'],atol=1e-12)
  checks[m]={'matrix_pairs':len(v['pairs']),'classes':len(v['classes']),'occupied_spinors':c['nocc'],'symmetry_residual':float(abs(d-d.T).max()),'internal_sum_rule_relative_error':float(err),'orbit_sizes_match':c['orbit_sizes']==c['orbit_expect'],'max_class_estimator_std':max(x['ES_std'] for x in v['classes']),'ET_charge_sum':sum(a['q_ET'] for a in atoms)}
 out={'status':'源矩阵到汇总、ET原子电荷、ELF域求和算术一致','checks':checks,'not_validated':['不能将求和闭合用作文献 ES 的因子/自旋约定验证。','同 k 折叠矩阵的最近邻分摊不构成逐周期像离域积分。','未重跑 WAVECAR 积分、未验证 PAW 全电子重构、占据截断与网格收敛。','COGITO/DI/ET 几何与方法不完全相同，跨数据链只作定性比较。']}
 (D/'09_审查/共享与局域化审查.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps(out,ensure_ascii=False))
if __name__=='__main__':run(Path(sys.argv[1]))
