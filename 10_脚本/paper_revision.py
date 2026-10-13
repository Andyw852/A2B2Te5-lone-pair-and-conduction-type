"""Paper structure, explicit estimator boundaries, and data-driven interpretations."""
import json,re
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

def revise(document,D,used):
 s=BeautifulSoup(document,'html.parser')
 M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
 rt=D/'06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI'
 es=json.loads((rt/'DI结果/ES_逐接触汇总.json').read_text())
 et=json.loads((rt/'ET_全密度Bader汇总.json').read_text())
 elf=json.loads((rt/'ELF_f_localization_domains.json').read_text())
 bond=json.loads((D/'09_审查/成键复算.json').read_text())
 def para(old,new):
  for p in s.find_all('p'):
   if old in p.get_text():p.clear();p.append(BeautifulSoup(new,'html.parser'))
 def table(headers,rows,ident):
  t=s.new_tag('table',id=ident);head=s.new_tag('tr')
  for x in headers:
   td=s.new_tag('th');td.string=x;head.append(td)
  t.append(head)
  for row in rows:
   tr=s.new_tag('tr')
   for x in row:
    td=s.new_tag('td');td.string=str(x);tr.append(td)
   t.append(tr)
  wrap=s.new_tag('div',attrs={'class':'scroll'});wrap.append(t);return wrap
 def addafter(node,txt):
  p=s.new_tag('p');p.append(BeautifulSoup(txt,'html.parser'));node.insert_after(p);return p
 lead=s.select_one('p.lead');lead.clear();lead.append(BeautifulSoup('组分替换如何改变孤对相关的局域成键及近带边电子态，是本报告的核心问题。证据按<b>元素与结构 → 全元素成键 → 电荷与 ELF 局域化 → 带边杂化 → 缺陷平衡 → 输运响应</b>组织。各环节分别给出数值、规律和方法边界；宿主成键的变化不能直接替代导电类型的缺陷证据。','html.parser'))
 intro=s.find(id='s1').find('p')
 addafter(intro,'<b>符号定义：</b>A₂B₂Te₅ 中 A=Pb 或 Sn，B=Sb 或 Bi，Te 为碲。A、B 表示化学式中的元素位置，<b>不是平均原子</b>。Pb₂Sb₂Te₅ 对应 A=Pb、B=Sb；Pb₂Bi₂Te₅ 对应 A=Pb、B=Bi；Sn₂Sb₂Te₅ 对应 A=Sn、B=Sb；Sn₂Bi₂Te₅ 对应 A=Sn、B=Bi。只有标为“平均电荷”或“方向平均”的量才做相应平均。')
 para('本轮纠正了上一版','<b>成键指标定义：</b>COGITO 的 elec/bond 对应 ICOOP，原始 JSON 字段为 coop；ICOHP 是哈密顿加权布居积分。两者分别回答轨道重叠和占据成键能贡献，不能把 2×ICOOP 解释成电子共享 ES。')
 para('为避免把 B–Te','全元素对表覆盖 A–Te、B–Te、A–B、A–A、B–B 与 Te–Te。A–Te 和 B–Te 只将 4 Å 内的近邻列为配位壳层；其余组合列最近接触。<b>Sn₂Sb₂Te₅ 的 Sn–Te 在源表合并为 3.13 Å、每胞 12 个接触，5.36 Å 属远邻，不能标为配位长键。</b>重数表示每胞接触数；正 ICOHP 只表示所用分解下占据反键贡献占优，不等于总排斥能。')
 para('口径边界（必须','<b>共享估计的证据级别：</b>这里将已保存的矩阵记为 D̂，定义 Ê=2D̂/m，保留源脚本的数值归一化。D̂ 由同 k 点、Löwdin 正交化的 PAW 赝轨道盆内重叠构造；Ê 把折叠矩阵按最近邻重数 m 分摊。它是<b>近似共享描述量，不是已验证的文献逐周期像 ES</b>。更远周期像可能贡献于同一个矩阵元，除以 m 并不能证明这些贡献全部来自最近邻。正交化后的求和规则主要检验内部闭合，不能验证 PAW 全电子重构、因子约定或周期像分解。运行日志实际使用 Sb 系 68、Bi 系 88 个占据 SOC 自旋子态；34 个双占据空间轨道的密度示例不能用作全部 SOC 结果的误差估计。4.4% 的总电子数差也不是 ES 的误差条。逐周期像可通过保留交叉 k 信息或局域轨道实现，MLWF 不是唯一数学路线。')
 para('Raty 等人的表格','Raty 等人的文献使用 ES=2×DI、Bader 电荷转移和固定 ELF 等值面内的电子积分区分共享、转移与局域化。本表保留这种并列组织方法，但将本项目的 Ê 明确标为近似共享描述量。ET 使用全密度定盆、价电子积分相对 ZVAL 的电荷；只有芯电子局域在相应盆的假设成立时才等价于完整密度相对中性原子的电荷。ELF 列为 f=0.85 连通域的积分，不能自动解释为完整孤对的电子数。')
 para('表中每一列都对应','B–Te 短接触的 Ê=1.29–1.32 e，长接触为 0.76–0.87 e；A–Te 两组为 0.82–0.90 e，层间 Te–Te 为 0.42–0.43 e。在同一估计定义下，短 B–Te/长 B–Te 为 1.5–1.7，短 B–Te/层间 Te–Te 为 3.0–3.2。A–Te 因而不能从层内成键讨论中省略。Pb→Sn 对短 B–Te 影响小，Sb→Bi 主要降低长 B–Te 的共享估计；这些四材料差异是描述性趋势，尚无收敛误差条支持统计显著性。层间非零 Ê 提示存在电子耦合，应结合 ELF、几何和电荷研究；<b>不能凭分摊估计直接排除范德华作用或证明严格单接触共享大小</b>。')
 para('全元素对比较显示','六类元素对中，B–Te 短接触具有最大的负 ICOHP 与正 ICOOP；A–Te 的 ICOHP 约 −1.05 至 −1.29 eV/接触，也明显参与层内骨架。B–Te 长接触约 −0.50 至 −0.70 eV/接触，弱于 A–Te 近邻。A–A、B–B 和 A–B 在所列近邻中的积分贡献较小；层间 Te–Te 的 ICOHP 为 +0.36 至 +0.39 eV/接触、ICOOP 为 −0.0094 至 −0.0053 e/接触，表明占据态的净反键贡献。非零共享估计与反键积分可以共存：共享衡量离域，ICOHP 衡量所选投影的能量贡献，它们不应给出相同符号。')
 para('上述投影诊断较小','投影诊断较小支持在当前基组内比较 ICOHP/ICOOP，但不替代基组和 k 网格收敛。Mulliken 总人口恢复约 68/88 个价电子；on-site 人口与完整电荷不是同一分区，s 人口大于 2 不能解释为多出自由电子。A–Te、B–Te、Te–Te 应与各自电荷和 ELF 空间分布共同解释。现有数据支持层内成键不均匀及层间耦合的候选图景，不足以单独量化孤对引起的总能量或指定 p/n 类型。')
 para('适合论文的当前表述','<b>综合结论：</b>A–Te 与 B–Te 共同构成层内成键骨架，B–Te 的短长分化比 A–Te 更明显。Pb/Sn 替换主要改变 A–Te 和近带边 A-s 混合，Sb/Bi 替换同时改变 B–Te、电荷转移与 ELF 阈值域。层间 Te–Te 的共享估计非零，但其严格逐周期像 ES 尚未建立。孤对相关机制须由局域化形状、全元素成键、电荷与带边 s/p 混合共同支持；导电类型则由化学势条件下的缺陷竞争和电荷补偿判断。')
 para('文献口径的 ES =','共享矩阵、逐原子 ET 和 ELF=0.80/0.85/0.90 连通域积分均有对应数值记录。Ê 是按重数分摊的赝轨道共享估计；内部求和闭合不构成严格 ES、PAW 全电子重构或周期像分解的独立验证。完整三维谷简并度及缺陷—完美结构同方法能量分解不在现有证据范围内。')
 # Correct misleading definitions anywhere in tables.
 for td in s.find_all('td'):
  txt=td.get_text()
  if '严格按文献逐接触 ET' in txt:td.string='ET 是逐原子电荷；层间列取两侧 Te 的 |q| 平均。ET 不需要 pair-density，不能与胞级 Qabs 混用。'
  if '本轮已用 SOC KS' in txt:td.string='保存矩阵可复现；报告使用 Ê=2D̂/m 的近似共享描述量。其归一化与逐周期像解释未独立验证。'
  if '用更密介电' in txt:td.string='在统一电子结构与输运网格下检验双极输运及其方向性。'
 # Add all relevant A-Te classes to the ES-organized table.
 rows=[];curves=[]
 for m in M:
  A='Pb' if m.startswith('Pb') else 'Sn';B='Sb' if 'Sb' in m else 'Bi'
  classes=es[m]['classes']
  pick=lambda x:sorted([c for c in classes if c['layer']=='层内' and set(c['pair'].split('-'))=={x,'Te'}],key=lambda c:c['d'])[:2]
  bc=pick(B);ac=pick(A);tc=next(c for c in classes if c['layer']=='层间' and c['pair']=='Te-Te')
  vals=bc+ac+[tc];curves.append([c['ES'] for c in vals])
  domains=[d for d in elf[m]['thresholds']['0.85']['domains'] if d['n_voxels']>=4 and d['nearest_atom_index'] in et[m]['gap_Te_indices']]
  rows.append([m,*['%.3f (%.2f Å; m=%d)'%(c['ES'],c['d'],c['m']) for c in vals],f"{et[m]['gap_Te_mean_abs_ET']:.3f}",f"{max(d['electrons'] for d in domains):.3f}"])
 s.find(id='raty-table').parent.replace_with(table(['材料','短 B–Te Ê/e','长 B–Te Ê/e','A–Te 1 Ê/e','A–Te 2 Ê/e','层间 Te–Te Ê/e','层间 Te |ET|/e','ELF 最大域/e'],rows,'raty-table'))
 fig,axs=plt.subplots(1,3,figsize=(15,4.8),layout='constrained');x=np.arange(4)
 for i,(label,col) in enumerate(zip(['Short B–Te','Long B–Te','A–Te 1','A–Te 2','Interlayer Te–Te'],['#216882','#83afc1','#388863','#98bc74','#d29139'])):
  axs[0].bar(x+(i-2)*.15,np.array(curves)[:,i],.15,label=label,color=col)
 axs[0].set_title('(a) E-hat: provisional, symmetry check failed',fontsize=10);axs[0].set_ylabel('e (folded quantity / multiplicity)');axs[0].set_ylim(0,1.85);axs[0].legend(fontsize=7,frameon=False,ncol=2,loc='upper left')
 axs[1].bar(x,[et[m]['gap_Te_mean_abs_ET'] for m in M],color='#388863');axs[1].set_title('(b) Gap-facing Te |ET|');axs[1].set_ylabel('e / atom')
 axs[2].bar(x,[float(r[-1]) for r in rows],color='#ae7098');axs[2].set_title('(c) ELF domain, f=0.85');axs[2].set_ylabel('e / selected connected domain')
 for ax in axs:ax.set_xticks(x,M,rotation=25,ha='right');ax.grid(axis='y',alpha=.15)
 for ext in ['png','pdf','svg']:fig.savefig(D/'02_图表'/('ES_ET_ELF_文献口径.'+ext),dpi=240)
 plt.close(fig)
 (D/'09_审查/共享电荷局域化绘图数据.json').write_text(json.dumps({'materials':M,'sharing_estimator':curves,'sharing_status':'provisional: symmetry orbit check failed','gap_Te_abs_ET':[et[m]['gap_Te_mean_abs_ET'] for m in M],'ELF_selected_domain':[float(r[-1]) for r in rows]},ensure_ascii=False,indent=2),encoding='utf-8')
 # Refresh the inlined plot after replacement.
 import base64
 for im in s.find_all('img'):
  if 'ES–ET–ELF 双面板' in im.get('alt',''):
   im['src']='data:image/png;base64,'+base64.b64encode((D/'02_图表/ES_ET_ELF_文献口径.png').read_bytes()).decode()
   im.parent.find('figcaption').string='共享、电荷与局域化的分量对照。左图五类接触的 Ê 明确包括 A–Te，因对称轨道计数未通过而仅为待验证方法输出；中图是层间 Te 的逐原子 |ET|；右图为 ELF=0.85 最大选定连通域电子数。三种分区和归一化不同，分别讨论，不构造单一活性分数。Bi 体系的 |ET| 和该阈值域积分较高，但不意味着所有元素的孤对均更活跃。'
 # Full element charge and threshold tables, explicitly separate datasets.
 anchor=s.find(id='raty-table').parent
 literature=table(['文献体系','层内 ES/e（接触重数）','层间 ES/e（接触重数）','所选原子 TET/e'],[['GST225（Raty 2021 表2）','1.36（3×）','0.36（3×）','0.30（间隙 Te）'],['α-PbO（同表）','0.95（4×）','0.03（4×）','1.16']], 'literature-table')
 anchor.insert_after(literature)
 addafter(literature,'文献表 2 的数值保留作已发表结果参照。括号是接触重数，不再把 ES 乘一次重数作为单键值。GST225 的 TET 对应间隙处 Te 原子，不能与本胞所有原子的电荷绝对值总和对比。文献的严格 ES 与本表 Ê 的计算路径及验证程度不同，不以数值接近作为方法正确性的证明。<a href="../06_COGITO成键分析与COHP/源数据/文献方法/Raty2021.pdf">Raty 2021 原文（第 8 页，表 2）</a>。')
 charge_rows=[];threshold_rows=[]
 for m in M:
  for el,q in et[m]['per_element_mean_q_ET'].items():charge_rows.append([m,el,f'{q:+.4f}'])
  threshold_rows.append([m,*[f"{elf[m]['thresholds'][t]['total_electrons_in_domains']:.4f}" for t in ['0.8','0.85','0.9']]])
 anchor.insert_after(table(['材料','元素','平均 ET=q / e·atom⁻¹'],charge_rows,'element-et-table'))
 addafter(s.find(id='element-et-table').parent,'逐元素电荷显示 A、B 均向 Te 侧转移电子。A 的平均电荷约 +0.63 至 +0.66 e，B 约 +0.51 至 +0.69 e；Te 为负。Sb→Bi 时 B 的电荷转移和层间 Te 的 |ET| 增大，Pb→Sn 的变化较小。层间 Te 的 |ET| 小于全体 Te 平均值，说明同为 Te 的不同结构位置不可混为一个等效原子。')
 anchor.insert_after(table(['材料','ELF≥0.80 全域积分/e','ELF≥0.85 全域积分/e','ELF≥0.90 全域积分/e'],threshold_rows,'elf-threshold-table'))
 addafter(s.find(id='elf-threshold-table').parent,'ELF 阈值升高后积分明显减少，表明“电子数”依赖所选区域。它不能直接成为不依赖阈值的孤对占据数。Sb 体系在 f=0.80 存在归于 Sb 邻域的较大连通域，而 Bi 体系相同阈值主要为 Te 域；因此不能把层间 Te 在 f=0.85 的趋势推广为 Bi 比 Sb 更强的阳离子孤对。域归属还受网格和连通性影响。')
 addafter(anchor,'<b>结构口径：</b>COGITO 与 DI 的距离分别来自各自源结构。Pb₂Bi₂Te₅ 的 B–Te 为 3.08/3.25 Å（COGITO）和 3.09/3.30 Å（DI）；层间 Te–Te 为 3.61 和 3.63 Å。Sn₂Sb₂Te₅ 的 COGITO Sn–Te 合并壳层为 3.13 Å，而 DI 结构给出 3.12/3.14 Å。故这些数值只作方法内趋势和跨方法定性对照，不能当作完全相同几何的严格逐键交叉验证。')
 warning=s.new_tag('p',attrs={'class':'note'});warning.string='共享矩阵的质量检查：四材料的 38 个 k 点中，第 16、19、35、38 点（从 1 编号）的记录轨道数为 12、6、12、6，预期为 6、2、6、2；累计 308 对预期 288。轨道计数检验未通过，提示空间群展开实现仍需验证。下表 Ê 保留为待验证的方法输出，不承担定量成键排序或排除范德华间隙的结论；ET 和 ELF 列来自独立计算，不因该项自动失效。'
 anchor.insert_before(warning)
 for cap in s.find_all('figcaption'):
  if '蓝柱为 ICOHP（左轴' in cap.get_text():cap.string='六类元素对最近壳层的积分成键指标。左面板 ICOHP，右面板 ICOOP，颜色表示材料；每个面板使用独立且统一的纵轴。B–Te 和 A–Te 的成键贡献显著，阳离子之间及 Te–Te 接触较弱或净反键。不同物理量的柱高不可跨面板比较。'
  if '：归档 ELF' in cap.get_text():
   m=next(m for m in M if m in cap.get_text());A='Pb' if m.startswith('Pb') else 'Sn';B='Sb' if 'Sb' in m else 'Bi'
   trend='Sb 邻域呈现明显高 ELF 区，而 A 位邻域的值较低。' if B=='Sb' else 'Bi 与 A 位邻域以较低 ELF 环带为主，Te 周围保持高 ELF 区。'
   cap.string=f'{m} 的 (110) ELF 切面，色标统一为 0–1，圆内标注为具体元素。{trend}Te 的高值区不是阳离子孤对电子数；切面形状提示元素选择性的局域化，需结合阈值域积分和近带边 {A}-s/{B}-s 与 Te-p 投影判断。'
 # Remove work-log wording from visible text, retain exact source filenames.
 repl={'本轮已算出 ES':'共享估计与电荷分区','本轮':'本研究','重算':'计算','旧图':'图','原 HSE+SOC':'HSE+SOC','归档 ELF':'ELF','原主报告表 10 的复核与 SOC 扩展':'缺陷条件下的 SOC 输运','沿用归档图片，本研究新计算的是图 2 的近边能窗统计。':'能窗统计用于量化近带边投影。','反键/排斥':'净反键','反键或排斥贡献':'净反键贡献','文献逐接触 ET':'逐原子 ET','更红':'局域函数更高','可测共享':'非零共享估计','单接触 ES = 2×折叠 DI / 对称重数 m':'近似 Ê = 2×D̂ / m'}
 for t in list(s.find_all(string=True)):
  if t.parent.name in ['style','script']:continue
  v=str(t)
  for a,b in repl.items():v=v.replace(a,b)
  if v!=str(t):t.replace_with(v)
 for td in s.find_all('td'):
  if '保存矩阵可复现' in td.get_text():td.string='源矩阵到汇总算术一致，但 k 对称轨道计数检验未通过。Ê 仅为待验证共享估计，不作严格 ES 证据。'
 for p in s.find_all('p'):
  if '层间非零 Ê 提示存在电子耦合' in p.get_text():
   p.string=p.get_text().replace('层间非零 Ê 提示存在电子耦合，应结合 ELF、几何和电荷研究；','这些数值趋势仅描述待验证的 Ê 输出，不形成独立物理证据。层间耦合需要结合 ELF、几何和电荷研究；')
 # Reorder complete sections and renumber all cross-references consistently.
 order=['s1','s2','s4','s3','s5','s6','s7','s8'];body=s.body
 blocks=[s.find(id=i) for i in order]
 for b in blocks:b.extract()
 for b in blocks:body.append(b)
 mapping={'4':'3','3':'4'}
 for t in list(s.find_all(string=True)):
  if t.parent.name in ['style','script']:continue
  v=re.sub(r'(?<![\d.])([34])\.(\d+)(?=\s|节|\b)',lambda m:mapping[m[1]]+'.'+m[2],str(t)) if t.parent.name in ['h3','h2'] else str(t)
  if v!=str(t):t.replace_with(v)
 titles=['研究问题与元素定义','计算方法与介电输入','全元素成键、电荷与局域化','带边简并态与能量窗口','缺陷竞争与电荷平衡','载流子条件与输运','孤对相关机制与综合结论','数据审查与复现']
 nav=s.find('nav');nav.clear()
 for i,(b,title) in enumerate(zip(blocks,titles),1):
  b.find('h2').string=f'{i}　{title}';a=s.new_tag('a',href='#'+b['id']);a.string=f'{i} {title}';nav.append(a)
 # Every figure and table has one sequential number and an interpretation nearby.
 for i,cap in enumerate(s.find_all('figcaption'),1):
  txt=re.sub(r'^图\s*[\d]+(?:[-a-z]\d*)?\s*','',cap.get_text(' ',strip=True));cap.string=f'图 {i}　'+txt
 # References are descriptive so they remain valid when figures move.
 for t in list(s.find_all(string=True)):
  if t.parent.name in ['style','script']:continue
  txt=str(t).replace('图 2 的能窗统计','带边能窗图的统计').replace('图 2 的近边能窗统计','带边能窗图的统计')
  txt=re.sub(r'图([1-4])([ab])',lambda m:'该材料的'+('稳定性面板' if m[2]=='a' else '载流子面板'),txt)
  if txt!=str(t):t.replace_with(txt)
 for i,t in enumerate(s.find_all('table'),1):
  c=s.new_tag('caption');c.string=f'表 {i}　'+ ' / '.join(x.get_text(' ',strip=True) for x in t.find_all('th')[:3]);t.insert(0,c)
 return str(s)
