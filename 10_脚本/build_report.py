"""Render the revised standalone report from archived numeric analyses and source templates."""
from pathlib import Path
import sys,json,base64,hashlib,re,html
sys.path.insert(0,str(Path(__file__).parent/'deps'))
from bs4 import BeautifulSoup
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
M=['Pb2Sb2Te5','Pb2Bi2Te5','Sn2Sb2Te5','Sn2Bi2Te5']
def run(D):
 script=D/'10_脚本';templates=script/'报告模板';templates.mkdir(parents=True,exist_ok=True)
 oldpath=D/'01_报告/完整分析报告.html'
 if not (templates/'缺陷章节.html').exists():
  old=BeautifulSoup(oldpath.read_text(),'html.parser');pieces=[]
  ids=['overview','m1','m2','m3','m4','lines','limits'];n=1
  for ident in ids:
   part=old.find(id=ident);part.find('h2').string=f'5.{n}　'+re.sub(r'^\d+\s*','',part.find('h2').get_text());n+=1;pieces.append(str(part))
  for section in old.find_all('section'):
   h=section.find('h2')
   if h and h.get_text().startswith('8'):
    h.string='5.8　化学势与载流子原模型数值';pieces.append(str(section))
   elif h and h.get_text().startswith('9'):
    h.string='参考文献与缺陷源数据';(templates/'缺陷参考文献.html').write_text(str(section))
  soup=BeautifulSoup(''.join(pieces),'html.parser');images={hashlib.sha256(p.read_bytes()).hexdigest():p for p in (D/'02_图表').glob('*.png')}
  for im in soup.find_all('img'):
   if im['src'].startswith('data:'):
    b=base64.b64decode(im['src'].split(',',1)[1]);p=images[hashlib.sha256(b).hexdigest()];im['src']=str(p.relative_to(D))
  (templates/'缺陷章节.html').write_text(str(soup))
 bands=json.loads((D/'07_带边能窗/带边能窗统计.json').read_text());bond=json.loads((D/'09_审查/成键复算.json').read_text());diel=json.loads((D/'08_介电复核/介电复核.json').read_text());tr=json.loads((D/'09_审查/输运复算.json').read_text());sn=json.loads((D/'09_审查/SnBi缺陷输运对应.json').read_text())
 bader=json.loads((D/'06_COGITO成键分析与COHP/源数据/Bader_TET_ES原始密度/Bader_charge_recalculation.json').read_text())
 def table(headers,rows,ident=''):
  return '<div class="scroll"><table '+('id="'+ident+'"' if ident else '')+'><thead><tr>'+''.join('<th>'+str(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'
 used=[]
 def fig(name,caption):
  p=D/'02_图表'/name;used.append(str(p.relative_to(D)));return '<figure><img alt="'+html.escape(caption,quote=True)+'" src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'"><figcaption>'+caption+' 图像对应的可追溯数据或原图文件均保留在交付目录。</figcaption></figure>'
 def win(m,e,w):return next(r for r in bands[m]['windows'] if r['edge']==e and r['window_meV']==w)
 def source(path,label='源数据'):return '<a href="../'+path+'">'+label+'</a>'
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 fig1,axs=plt.subplots(1,2,figsize=(11,4.2),layout='constrained');x=np.arange(4)
 for ax,k,tit in zip(axs,['eps_inf','eps_static'],[r'(a) Electronic $\epsilon_\infty$',r'(b) Total static $\epsilon_0$']):
  for off,ix,c,l in [(-.17,0,'#357d99','xx'),(.17,2,'#d39446','zz')]:
   vs=[diel[m][k][ix][ix] for m in M];ax.bar(x+off,vs,.32,color=c,label=l)
   for xi,v in zip(x+off,vs):ax.text(xi,v,f'{v:.1f}',ha='center',va='bottom',fontsize=8)
  ax.set_xticks(x,M,rotation=15,ha='right');ax.set_title(tit,loc='left');ax.legend(frameon=False);ax.set_ylabel('Relative permittivity')
 fig1.suptitle('SOC DFPT | 9 x 9 x 2 | tensors used by both transport models')
 for ext in ['png','pdf','svg']:fig1.savefig(D/'02_图表'/('介电_电子与总静态张量.'+ext),dpi=300)
 plt.close(fig1)
 fig1,axs=plt.subplots(1,3,figsize=(13,4.5),layout='constrained');colors=['#247ba0','#d28c33','#329875','#9b70ac'];temps=tr['Sn2Bi2Te5_SOC']['temperatures']
 for row,c in zip(sn,colors):
  for variant,ls in [('noSOC','--'),('SOC','-')]:
   for ax,key in zip(axs,['S_all','sigma_all','PF_all']):ax.plot(temps,row[variant][key],color=c,ls=ls,label=f"v{row['vertex']} {variant}")
 for ax,yl in zip(axs,['Seebeck (microV/K)','Conductivity (S/m)','PF (mW/m/K²)']):ax.set_xlabel('Temperature (K)');ax.set_ylabel(yl);ax.grid(alpha=.18)
 axs[0].axhline(0,color='gray',lw=.7);axs[1].set_yscale('log');axs[2].legend(ncol=2,fontsize=8,frameon=False)
 fig1.suptitle('Sn2Bi2Te5 | fixed defect-model net doping | solid: SOC, dashed: no SOC')
 for ext in ['png','pdf','svg']:fig1.savefig(D/'02_图表'/('SnBi_缺陷浓度_SOC输运.'+ext),dpi=300)
 plt.close(fig1)
 H=[];H.append('<h1>A₂B₂Te₅：从局部成键到带边、缺陷与热电输运</h1><p class="subtitle">Pb₂Sb₂Te₅ · Pb₂Bi₂Te₅ · Sn₂Sb₂Te₅ · Sn₂Bi₂Te₅ | POSCAR_B | 2026-09-17 审查版</p>')
 H.append('<p class="lead">本报告从当前 SOC-DFPT 输出核对介电输入，逐值复核输运对比，重读 HSE+SOC 电荷投影，并纠正成键指标定义。主线是：<b>局部成键与轨道混合 → 带边态及屏蔽 → 缺陷电荷平衡 → 给定载流子浓度下的输运</b>。前两步提供宿主电子结构背景；缺陷导致何种载流子，仍须由缺陷能量、电离和补偿决定。</p>')
 H.append('<nav>'+''.join(f'<a href="#s{i}">{i} {t}</a>' for i,t in enumerate(['研究问题','方法与介电','带边简并与能窗','ES/TET与成键','缺陷条件','输运','联合解释','审查与复现'],1))+'</nav>')
 H.append('<section id="logic-chain"><h2>逻辑总览：四类计算怎样接起来</h2><p>图中四个面板按“局部结构 → 近带边电子态 → 缺陷平衡 → 给定掺杂下输运”的顺序排列。短/长 B–Te 键的 ICOHP 比值说明局部键贡献不均匀；50 meV VBM 能窗中的 A-s 投影说明带边混合具有材料依赖；缺陷模型的中位数 <i>p−n</i> 只表示现有化学势顶点下的净载流子倾向；最后的功率因子是在固定 +10<sup>19</sup> cm<sup>−3</sup>、300 K 和 SOC 条件下的输运响应。四个面板使用不同计算和归一化，图用于组织证据链，不表示已经证明“孤对电子直接导致输运性能”。</p>')
 H.append(fig('证据链_成键到带边缺陷输运.png','图 0　从局部成键、近带边投影、缺陷平衡到输运的证据链总览。四个面板的横轴均为材料；纵轴分别是短/长键的成键贡献比、50 meV VBM 窗口的 A-s 投影、化学势顶点的净载流子趋势和固定掺杂的功率因子。应逐面板读取其单位和归一化，不能跨面板把数值相除或当作因果强度。'))
 H.append('</section>')
 H.append('<section id="s1"><h2>1　研究问题与可支持的结论</h2><p>四种材料把 A 位 Pb/Sn 与 B 位 Sb/Bi 分别替换，提供了比较化学组分效应的四个样本。关心的不是孤对电子“有或没有”，而是 ns 相关态在不同能量、k 点和键环境中怎样参与杂化，以及这些变化能否与缺陷和输运建立一致联系。</p>')
 H.append(table(['问题','当前数据直接支持','推论的边界'],[
 ['带边是否包含阳离子 s？','HSE+SOC 投影在精确带边与有限能窗内均可量化。','投影占比不是孤对电子数，也不是电导占比。'],
 ['局部键是否分化？','B–Te 短长键差 0.15–0.19 Å，|ICOHP短|/|ICOHP长| 2.43–2.85。','说明成键贡献不均匀；未隔离孤对引起的结构能量变化。'],
 ['电子共享、转移与孤对如何量化？','当前体系可直接复算 ICOHP/ICOOP、Mulliken 与 Bader 电荷，以及 ELF 空间局域化。','这些量分别描述成键能、轨道重叠、投影电荷、实空间电荷和局域化；没有定义统一的单一 ES/TET 标量，因此按各自物理含义解释。'],
 ['为什么 p 或 n？','18 个顶点的缺陷排序与原模型判型可复现。','多数顶点存在负形成能/饱和；Pb₂Bi₂Te₅ 的 n 型点强补偿，不能视为稳健实验预测。'],
 ['介电是否正确传递到输运？','四套当前收敛的 9×9×2 DFPT 张量与 8 份 AMSET 设置一致；3400 行 SOC 对比逐值通过。','该结论针对本报告使用的屏蔽输入；缺陷形成能另按其源记录说明。']]))
 H.append('</section><section id="s2"><h2>2　计算口径、介电数据链</h2>')
 H.append(table(['环节','当前采用的口径','源数据'],[
 ['HSE06+SOC 带边','LORBIT=11；Sb 系 PROCAR_OPT，Bi 系 PROCAR 中零权重路径点；9 原子胞。',source('07_带边能窗/源数据/','输入与压缩原始 PROCAR')],
 ['COGITO 成键','PBE-D3(BJ)，无 SOC；COHP/ICOOP 和 on-site/Mulliken 划分。',source('06_COGITO成键分析与COHP/结果/','原始键表和原子表')],
 ['当前介电','SOC DFPT，9×9×2，NKPTS=162；ε₀=ε∞+εion；作为已验证收敛的屏蔽输入。',source('08_介电复核/源数据/','原始 OUTCAR.gz 与输入')],
 ['AMSET 基线','无 SOC 波函数与形变势 + 当前 SOC 介电 + HSE+SOC 带隙。',source('01_报告/源数据/','基线 transport/settings')],
 ['AMSET SOC','SOC 波函数与 SOC 形变势 + 同一当前介电与带隙；与基线使用相同浓度、温度。',source('05_SOC_AMSET补算/','SOC 归档与对比')],
 ['缺陷','沿用原形成能、化学势及电荷中性输出；81 原子能量体系，72 原子结构参考不直接配能。',source('03_缺陷形成能与PN型数据/','四材料数据包')]]))
 H.append('<p><b>缺陷数据的使用边界：</b>缺陷形成能、化学势顶点和载流子平衡均沿用源模型记录。本报告只复核这些记录中的最低缺陷排序、相对能量差与 p/n 判型，并用它们讨论组分替换下的反位竞争；不把缺陷数值重新参数化为介电或孤对电子的直接量化指标。</p>')
 H.append(fig('介电_电子与总静态张量.png','图 1　四种材料当前收敛介电张量的电子部分 ε∞ 与总静态部分 ε₀。横轴为材料，纵轴为相对介电常数，柱高分别读取 xx 与 zz；ε₀ 已包含离子贡献。读图时先比较同一材料的 ε∞→ε₀ 增量，再比较 xx/zz 的方向差异；该图说明输运屏蔽的各向异性，不能单独推出某一散射机制的增减。完整 3×3 张量和原始 OUTCAR 位于介电复核源数据。'))
 H.append(table(['材料','ε∞ xx / zz','εion xx / zz','ε₀ xx / zz','POP 频率 / THz'],[[m,*[f"{diel[m][k][0][0]:.3f} / {diel[m][k][2][2]:.3f}" for k in ['eps_inf','eps_ion','eps_static']],diel[m]['pop_frequency']] for m in M], 'dielectric-table'))
 H.append('<p>xx 分量普遍大于 zz，说明这批结构的屏蔽具有明显方向性。ε∞ 表示电子响应，εion 是晶格离子贡献，二者相加得到静态 ε₀；AMSET 使用的正是当前统一口径的静态输入。POP 频率和介电张量的当前数值均可在介电复核 JSON 与原始 OUTCAR 中逐项追溯。</p></section>')
 H.append('<section id="s3"><h2>3　带隙与带边：简并子空间、近边态及能量窗口</h2><h3>3.1 为什么要从一个点扩展到一组态？</h3><p>单个极值点不能代表近带边全部态。本报告以 1 meV 容差识别近简并带，并对路径重复点去重；同时用 10–100 meV 能窗统计邻近态。SOC 简并态的单条波函数投影可受子空间内基变换影响，完整简并子空间的投影和更稳健。严格近简并检查与有限能窗统计分开记录。</p>')
 H.append('<ol><li>只读取 SOC PROCAR 的第一张电荷投影表，不把三个磁化分量相加。按原子和 s/p/d 轨道分组。</li><li>只保留能带路径点；将相差倒格矢的同一 k 位置去重，重复记录先平均。同一 k 的两条 SOC 带仍是两个态；不再额外乘“自旋简并 2”。</li><li>由路径全部已占据带和未占据带重新搜索 VBM/CBM。按自旋子电子计数，Sb 系前 68 带、Bi 系前 88 带为价带集合。0.1、1、5 meV 用于极值容差敏感性；10、25、50、100 meV 用于近边态统计。</li><li>价带窗为 [VBM−Δ,VBM]，导带窗为 [CBM,CBM+Δ]；各保留 k–band 态等权。Pα=Σ态wα / Σ态Σβwβ。另计算“各态先归一化再平均”，并保存两种定义的最大百分点差。</li></ol>')
 H.append('<p class="note"><b>这是一维已采样路径上的态集合分析。</b>k 点数不是三维谷简并度，窗口宽度也不等于温度或真实输运权重；没有加入 BZ 体积权重、速度、寿命或费米导数。25 meV 仅与室温 kBT 量级接近，高掺杂可采样更宽能区。若要得到轨道分辨输运或可靠的多谷权重，需要稠密三维网格以及与输运一致的电子结构。</p>')
 H.append(table(['材料','PROCAR 路径带隙 / eV','路径最小直接隙 / eV','类型','带边摘要值 / eV','AMSET scissor / eV','路径记录 → 去重 k'],[[m,f"{bands[m]['gap_eV']:.6f}",f"{bands[m]['min_direct_gap_eV']:.6f}",bands[m]['type'],json.loads((D/'03_缺陷形成能与PN型数据/带边轨道投影'/f'{m}_fatband_summary.json').read_text())['gap_eV'],tr[m+'_SOC']['settings']['bandgap'],f"{bands[m]['n_path_records']} → {bands[m]['n_unique_k']}"] for m in M],'bandgap-table'))
 H.append('<p>表中 PROCAR 路径极值、路径内最小直接隙和 AMSET scissor 是相互独立的记录，显示精度内已核对；直接/间接类型仍仅限当前路径，不能当作全 BZ 最终值。能窗统计只约束 VBM/CBM 附近态，高能未占据带不由该表背书。</p>')
 H.append(fig('带边_简并态与能窗.png','图 2　四材料去重后的 1、25、50、100 meV 带边窗口轨道组成。右侧为入选态数/独立 k 位置数。Other 为 Te-s 与所有 d 分量；条形图按已投影权重归一化，不能解释为电子数或电导百分比。'))
 rows=[]
 for m in M:
  a='Pb' if m.startswith('Pb') else 'Sn';b='Sb' if 'Sb' in m else 'Bi'
  for w in [1,25,50,100]:
   v=win(m,'VBM',w);c=win(m,'CBM',w);rows.append([m,w,f"{v['states']}/{v['unique_k']}",f"{v['percent'][a+'-s']:.2f}",f"{v['percent'][b+'-s']:.2f}",f"{v['percent']['Te-p']:.2f}",f"{c['percent'][a+'-s']+c['percent'][b+'-s']:.2f}"])
 H.append(table(['材料','窗口 meV','VBM 态/k','VBM A-s %','VBM B-s %','VBM Te-p %','CBM A-s+B-s %'],rows,'window-table'))
 H.append('<h3>3.2 哪些解释经过能窗检验仍成立？</h3><p><b>Sn₂Bi₂Te₅：</b>Sn-s 从 1 meV 的 24.8% 变为 50 meV 的 23.1%、100 meV 的 15.4%；近边 s/p 混合仍显著，但精确 VBM 的数字不能推广到宽能窗。25 meV 内仍只有一个独立 k 位置的两条带，说明当前路径在这一局部窗口的采样很有限。</p><p><b>Sn₂Sb₂Te₅：</b>VBM 的 Sn-s 在 1→100 meV 增至约 1.9%→5.8%，Te-p 从 20.7%→31.7%，Sb-p 的相对份额降低；精确 VBM 与较宽空穴能区应分别报告。CBM 的阳离子 s 合计从约 27.7% 降到 14.8%，窗口依赖同样明显。</p><p><b>Pb₂Sb₂Te₅：</b>VBM 的 Pb-s 和 Sb-s 合计在 1–100 meV 约 16.7%–17.7%，这个投影指标相对稳定；窄隙使双极输运值得关注，但不决定掺杂类型。<b>Pb₂Bi₂Te₅：</b>VBM 的 Te-p 随窗口扩展由约 43.4% 增到 54.6%。Bi 系 CBM 的 1 meV 全路径统计包含两个独立 k 位置、共四个态，说明近边态必须按简并子空间和能窗统计。</p>')
 maxnorm=max(r['normalization_difference_max_pp'] for m in M for r in bands[m]['windows']);H.append(f'<p>两种归一化顺序在本次全部窗口上的最大差为 {maxnorm:.3f} 个百分点；该敏感性已保存。原始电荷投影不完备、输出有三位小数舍入，百分比仅属于该投影子空间。'+source('07_带边能窗/入选态清单.csv','逐态 k、band、能量及入选窗口')+'；'+source('07_带边能窗/带边能窗统计.json','全部窗口与归一化敏感性')+'。</p>')
 H.append('<h3>3.3 完整能带图与局部统计互相约束</h3><p>四张 fatband 用于观察深价带的 s 态和近带边混合态，精确复算带隙见本节表。深价带 s 权重强，并不说明这些完全占据态直接导电；要讨论输运，必须转向费米能附近能窗。彩色宽度只能辨认轨道分量，不能区分成键/反键。</p>')
 for m in M:H.append(fig('fatband_'+m+'.png',f'图 3-{M.index(m)+1}　{m} 原 HSE+SOC fatband。横轴是带路径，纵轴是相对能量，彩色线宽是相应原子轨道投影权重；先在费米能附近识别近边混合，再结合图 2 的能窗统计判断其稳健性。深价带的宽线只表示已占据态投影强，不等于直接参与输运；沿用归档图片，本轮新计算的是图 2 的近边能窗统计。'))
 H.append('</section><section id="s4"><h2>4　电子共享 ES、电子转移 TET 与孤对相关成键</h2><h3>4.1 指标定义：不要把“共享”“转移”“局域化”当成同一个量</h3>')
 H.append(table(['量','含义与所需输出','本批可复核程度'],[
 ['ES（文献定义 2×DI）','两个 Bader basin 间 pair-density 离域指数的两倍，单位 e；需要 HF-like pair density 积分和明确键重数。','本轮已用 SOC KS 波函数 + 全密度 Bader 盆算出 ES（见 4.4 节表），并通过 QTAIM 求和规则校验；ES 仍不等同于 ICOOP。'],
 ['ICOOP / COOP','占据态积分的晶体轨道重叠布居，依赖非正交轨道重叠；可正可负。','COGITO bond_info.txt 的 elec/bond 对应 all_unique_bonds.json 的 coop，可逐行核验。'],
 ['ICOHP / COHP','哈密顿加权的成键能量分解；本图 COHP<0 成键、>0 反键。','原始 cohp 字段、键表与曲线可复核；不是键解离能或缺陷形成能。'],
 ['TET / ET','须定义完整 Bader 电荷、参考原子、选定原子/接触和归一化；文献表中的 ET 不是胞级 Qabs。','当前可复算 Bader 原子电荷和胞级 Qabs；严格按文献逐接触 ET 仍需指定参考原子和 pair-density 计算区域。'],
 ['Mulliken Qabs','Qabs = ½Σᵢ|qᵢᴹ|，qᵢᴹ 为 COGITO Mulliken 部分电荷，单位 e/9原子胞。','定义明确、总电荷残差近零，作为投影分区下的胞内电荷再分配描述量。'],
 ['Bader / ELF','Bader 需要 CHGCAR 与 AECCAR2 后处理；ELF 是空间局域化。','四材料均有同结构原始 CHGCAR/AECCAR/ELFCAR；Bader 表已完成复算。ELF 提供局域化形状，Bader 提供实空间电荷分区，二者都不直接给出孤对电子数。']]))
 H.append('<p class="note"><b>本轮纠正了上一版的一处数值定义错误：</b>COGITO 的 elec/bond 被整理为字段 <code>icobi_elec_per_bond</code>，进而用于 2×ICOBI。四材料共 402 条原始记录证实该列实际是 ICOOP。公开 COGITO 0.3.13 源码的 bond_info 写出逻辑也直接使用 <code>all_group_bonds[key]["coop"]</code>，末尾标注 total bond elec (iCOOP)。本轮已将派生字段更正为 ICOOP，撤回“2×当前 COGITO ICOBI”比较；原始文本与 JSON 保留不改。<a href="https://pypi.org/project/cogito-dft/0.3.13/">软件版本来源</a>；'+source('10_脚本/方法依据/COGITOpost_bond_info_excerpt.txt','写出代码节选与 SHA256')+'。</p>')
 H.append('<h3>4.2 短长键分裂：可核验的结构与能量证据</h3>')
 rows=[]
 for m in M:
  b='Sb' if 'Sb' in m else 'Bi';ss=sorted([r for r in bond[m]['bonds'] if set(r['pair'].split('-'))=={'Te',b}],key=lambda r:r['d_A'])[:2];s,l=ss
  rows.append([m,f"{s['d_A']:.2f} / {l['d_A']:.2f}",f"{l['d_A']-s['d_A']:.2f}",f"{s['ICOHP_eV_bond']:.4f} / {l['ICOHP_eV_bond']:.4f}",f"{abs(s['ICOHP_eV_bond']/l['ICOHP_eV_bond']):.3f}",f"{s['ICOOP_e_bond']:.4f} / {l['ICOOP_e_bond']:.4f}"])
 H.append(table(['材料','B–Te 短/长 d / Å','Δd / Å','ICOHP 短/长 / eV·bond⁻¹','|短|/|长|','ICOOP 短/长 / e·bond⁻¹'],rows,'bond-table'))
 H.append(fig('成键_ICOHP_ICOOP_电荷转移.png','图 4　从原始 COGITO 键表与原子表重绘。(a) 短/长 B–Te 的 −ICOHP；(b) 正确标识的 ICOOP，不能乘二当作 ES；(c) 各元素平均 Mulliken 部分电荷；(d) 本轮明确定义的 Qabs。9 原子胞含 A₂B₂Te₅，单位均在各轴标明。'))
 H.append('<p>B–Te 短键的 ICOHP 绝对值显著大于长键，正 ICOOP 也更大，说明两类接触的成键贡献分化。相比之下，A–Te 的几何短长差较小，其 ICOHP 不必严格随距离单调变化：同一长度变化还会伴随轨道能级、混合和其他邻居的改变。<b>键长差、ICOHP 差和 ES 不是三个同义指标。</b>四点比较中，Bi 替换 Sb 增大结构短长差，但 |ICOHP短|/|ICOHP长| 比值并没有随之增加；不能用一个指标给所有材料排出唯一“孤对活性”次序。</p>')
 H.append('<h3>4.2.1 全元素对的近邻成键对照</h3>')
 H.append('<p>为避免把 B–Te 作为唯一成键对象，本表按每种材料的化学式定义 A=Pb/Sn、B=Sb/Bi，逐一检查六类元素组合：A–Te、B–Te、A–B、A–A、B–B 和 Te–Te。A–Te、B–Te 列出最近的两组配位壳层（短/长）；其余组合列出最近壳层。ICOHP 保留原始符号（负值为成键贡献，正值为反键或排斥贡献），ICOOP 同样保留符号；重数是每个原始胞的等价接触数。远距离且数值很小的壳层仍保留在对应材料的原始 bond_info.txt/all_unique_bonds.json 中，未被当作“没有相互作用”。</p>')
 all_rows=[]; plot_data={}
 pair_order=['A–Te','B–Te','A–B','A–A','B–B','Te–Te']
 for m in M:
  A='Pb' if m.startswith('Pb') else 'Sn'; B='Sb' if 'Sb' in m else 'Bi'
  defs={'A–Te':(A,'Te'),'B–Te':(B,'Te'),'A–B':(A,B),'A–A':(A,A),'B–B':(B,B),'Te–Te':('Te','Te')}
  plot_data[m]={}
  for label in pair_order:
   x,y=defs[label];cand=[r for r in bond[m]['bonds'] if set(r['pair'].split('-'))==set([x,y])]
   cand=sorted(cand,key=lambda r:r['d_A'])
   # A-Te and B-Te have the chemically distinct short/long shells; other pairs use the nearest shell.
   chosen=[r for r in cand if r['d_A']<4.0][:2] if label in ('A–Te','B–Te') else cand[:1]
   plot_data[m][label]=chosen
   for j,r in enumerate(chosen):
    shell=('短' if j==0 else '长') if label in ('A–Te','B–Te') and len(chosen)==2 else '最近'
    all_rows.append([m,label,shell,f"{r['d_A']:.2f}",f"{r['ICOHP_eV_bond']:+.4f}",f"{r['ICOOP_e_bond']:+.4f}",str(r['count'])])
 H.append(table(['材料','元素对','壳层','距离 / Å','ICOHP / eV·bond⁻¹','ICOOP / e·bond⁻¹','胞内重数'],all_rows,'all-bond-table'))
 fig_all,axs_all=plt.subplots(1,2,figsize=(13,5),layout='constrained'); xpos=np.arange(len(pair_order)); width=.19
 for ax,key,title in zip(axs_all,['ICOHP_eV_bond','ICOOP_e_bond'],['(a) ICOHP (eV/contact)','(b) ICOOP (e/contact)']):
  for i,m in enumerate(M):
   vals=[plot_data[m][label][0][key] for label in pair_order]
   ax.bar(xpos+(i-1.5)*width,vals,width,label=m,color=colors[i])
  ax.axhline(0,color='0.35',lw=.7);ax.set_xticks(xpos,pair_order);ax.set_title(title,loc='left');ax.grid(axis='y',alpha=.15)
 axs_all[1].legend(frameon=False,fontsize=9)
 for ext in ['png','pdf','svg']:fig_all.savefig(D/'02_图表'/('全元素对_ICOHP_ICOOP.'+ext),dpi=300)
 plt.close(fig_all)
 H.append(fig('全元素对_ICOHP_ICOOP.png','图 4a　六类元素对最近壳层的 ICOHP 与 ICOOP。蓝柱为 ICOHP（左轴，负值代表成键，正值代表反键/排斥），橙柱为 ICOOP（右轴，正值代表轨道重叠的同相贡献）。A–Te 与 B–Te 的第二壳层见全元素对表；图中用最近壳层统一比较，不能把不同元素对的重数或归一化混为同一“键级”。'))
 H.append('<p>全元素对比较显示：B–Te 短键的负 ICOHP 和正 ICOOP 最大，是局部成键主导接触；A–Te 也有明显的成键贡献，但其最近壳层与 B–Te 的短长分裂不同，说明 A 位参与了层内骨架稳定而不是可忽略的旁观者。A–B、A–A 和 B–B 最近接触的 ICOHP/ICOOP 接近零或为弱反键，主要反映次近邻耦合。最近层间 Te–Te 的 ICOOP 为小幅负值且 ICOHP 为正，表示该接触不是典型强共价键；但这不等于电子共享严格为零，因为 ES=2×DI 对 Te–Te 仍为 0.42–0.43 e。因而孤对相关解释应同时考虑 A–Te、B–Te 的成键分化、Te–Te 的层间耦合以及 ELF/带边轨道混合，不能只依据 B–Te 一列。</p>')
 H.append('<h3>4.3 TET 应如何重新分析？</h3><p>对一个总电荷守恒的胞，Σqᵢ≈0；若把带符号电荷直接相加，得到的应是净电荷而非“总转移量”。本轮用 ½Σ|qᵢ| 避免把得失两端重复计数；也可等价看成正电荷总和，但必须声明每 9 原子胞。部分电荷的正号表示相对所用参考的电子减少，负号表示电子增加，不能与氧化态等同。</p>')
 H.append('<h3>4.3.1 Bader 重算：可用的实空间电荷证据</h3><p>四套 ES 目录中的 CHGCAR 与 AECCAR2 已用 pybader 0.3.12 重算。这里的参考是各 PAW 势的价电子数 ZVAL，定义 q<sub>Bader</sub>=ZVAL−N<sub>Bader</sub>；Qabs<sub>Bader</sub>=½Σ|q<sub>Bader</sub>|，单位为每个 9 原子胞。四个胞的电子数闭合误差不超过 2×10<sup>−6</sup> e，说明表格与密度输入一致。它适合在固定结构、固定 PBE 静态电荷密度下比较相对电荷再分配，属于半定量结构/电荷描述；不能直接叫作 TET、氧化态、孤对电子数或 HSE+SOC 结果。</p>')
 H.append(fig('Bader_charge_summary.png','图 4a　Bader 重算结果。左：按元素平均 q_Bader，正值表示相对 PAW 价电子参考的电子减少；右：Qabs_Bader=½Σ|q|，是每个 9 原子胞的总电荷再分配描述量。该图使用同一 PBE、非 SOC 静态密度口径，Pb₂Bi₂Te₅ 的静态步未启用 IVDW，因此跨材料绝对值比较应作为趋势而非严格统一泛函基准。'))
 H.append(table(['材料','Qabs_Bader / e·cell⁻¹','Σq_Bader / e','阳离子平均 q / e','Te 平均 q / e'],[[m,f"{bader[m]['qabs_bader_e_per_cell']:.6f}",f"{bader[m]['sum_q_bader']:+.2e}",f"{np.mean([v['mean_q_bader'] for k,v in bader[m]['elements'].items() if k!='Te']):+.6f}",f"{bader[m]['elements']['Te']['mean_q_bader']:+.6f}"] for m in M],'bader-table'))
 H.append('<p>表中阳离子平均 q 为正、Te 平均 q 为负，说明电子由阳离子侧向 Te 侧重新分配；Qabs_Bader 量化每个 9 原子胞的再分配幅度。四种材料均满足总电荷闭合，Bi 组合的 Qabs 略高于 Sb 组合。该结果支持“组分替换改变实空间电荷再分配”，但不单独证明孤对电子或导电类型；后者必须与 ELF 形状、带边 s/p 能窗和缺陷平衡共同判断。</p>')
 H.append(table(['材料','Qabs / e·cell⁻¹','Σqᴹ / e','Σ Mulliken 人口','Σ on-site 人口'],[[m,f"{bond[m]['Q_abs_half_e_cell']:.6f}",f"{bond[m]['mulliken_charge_sum']:+.6f}",f"{bond[m]['electron_population_sums']['mulli_occup']:.6f}",f"{bond[m]['electron_population_sums']['onsite_occup']:.6f}"] for m in M],'charge-table'))
 H.append('<p>在<b>同一 Mulliken 定义</b>下，Qabs 顺序为 Pb₂Bi₂Te₅ &gt; Sn₂Bi₂Te₅ &gt; Pb₂Sb₂Te₅ &gt; Sn₂Sb₂Te₅。它只描述 COGITO 投影分区下的胞内电荷再分配，不能替代文献逐接触 ET、氧化态或导电类型；完美胞的电荷转移与整体电中性并不矛盾。Mulliken Qabs、Bader Qabs、ES 和 ELF 的积分区域及归一化不同，不能把它们合并成单一“孤对活性”排序。</p>')
 H.append('<h3>4.4 按文献口径组织 ES–ET 表（本轮已算出 ES）</h3><p>Raty 等人的表格把 <b>ES 定义为 2×DI</b>（两个 Bader basin 之间的 HF 型 pair-density 离域指数）、把 <b>ET 定义为完整（芯+价）Bader 电荷与中性原子参考之差</b>、把孤对电子数定义为<b>固定 ELF 等值面（f = 0.85）内的局域化电子数</b>。本节按这三条定义给出可复算数值：ES 用 SOC 的 KS 波函数（38 个 IBZ k 点 + 空间群轨道展开）在<b>全密度定盆</b>上做盆内重叠积分，并已通过 QTAIM 求和规则 Sigma_B DI(A,B) = 2N_A（四材料最大偏差 0.11%–0.35%）和独立路径（全密度 Bader 电荷）交叉校验。</p>')
 _rt=(D/'06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI')
 _es=json.loads((_rt/'DI结果/ES_逐接触汇总.json').read_text())
 _et=json.loads((_rt/'ET_全密度Bader汇总.json').read_text())
 _elf=json.loads((_rt/'ELF_f_localization_domains.json').read_text())
 def _cls(mat,layer,pairs):
  c=_es[mat]['classes']
  for x in c:
   if x['layer']==layer and x['pair'] in pairs:
    return x
  return None
 rows=[]
 for m in M:
  b='Bi' if 'Bi' in m else 'Sb'
  intra=_cls(m,'层内',{b+'-Te','Te-'+b})
  inter=_cls(m,'层间',{'Te-Te'})
  et=_et[m]['gap_Te_mean_abs_ET']
  dom=[x for x in _elf[m]['thresholds']['0.85']['domains'] if x['n_voxels']>=4]
  lp=max([x['electrons'] for x in dom],default=0.0)
  rows.append([m,
   '%.2f e (%d×, %.2f Å)'%(intra['ES'],intra['m'],intra['d']) if intra else '—',
   '%.2f e (%d×, %.2f Å)'%(inter['ES'],inter['m'],inter['d']) if inter else '—',
   '%.3f'%et,'%.2f'%lp])
 H.append(table(['材料','层内短 B–Te<br>ES = 2×DI / e','层间 Te–Te<br>ES = 2×DI / e','层间 Te 的 ET / e<br>(全密度 Bader)','层间 Te 的 ELF 局域电子数<br>(f=0.85 最大盆) / e'],rows,'raty-table'))
 H.append('<p>表中每一列都对应明确的物理量和归一化：层内短/长 B–Te 与层间 Te–Te 是每个等价接触的 ES=2×DI，层间 ET 是全密度 Bader 电荷相对中性参考的绝对值，ELF 列是 f=0.85 最大局域化盆内价电子数。短 B–Te ES 为 1.29–1.32 e，约为层间 Te–Te ES 0.42–0.43 e 的 3.0–3.1 倍；层间 ET 为 0.289–0.383 e，ELF 盆积分为 0.96–1.07 e。由此可得：层内短键电子共享更强，但层间 Te–Te 仍有可测共享，不能把层间结构视为电子共享为零的理想范德华间隙。ET 和 ELF 说明层间区域存在电荷转移与局域化背景，却不单独证明局域化电子就是孤对；孤对解释需要与带边轨道混合和 ELF 空间形状联合判断。GST225 的 1.36/0.36/0.3 e 仅作为文献量级参照。</p>')
 fig_es, axs_es = plt.subplots(1,2,figsize=(11,4.6),layout='constrained'); xx=np.arange(len(M)); w=.24
 vals_short=[]; vals_long=[]; vals_gap=[]; vals_et=[]; vals_lp=[]
 for m in M:
  b='Bi' if 'Bi' in m else 'Sb'; ci=_cls(m,'层内',{b+'-Te','Te-'+b}); cl=[z for z in _es[m]['classes'] if z['layer']=='层内' and z['pair']==b+'-Te']; cl=sorted(cl,key=lambda z:z['d']); vals_short.append(cl[0]['ES']); vals_long.append(cl[1]['ES']); vals_gap.append(_cls(m,'层间',{'Te-Te'})['ES']); vals_et.append(_et[m]['gap_Te_mean_abs_ET']); dom=[z for z in _elf[m]['thresholds']['0.85']['domains'] if z['n_voxels']>=4]; vals_lp.append(max([z['electrons'] for z in dom],default=0))
 ax_es=axs_es[0]; ax_es.bar(xx-.5*w,vals_short,w,label='Intralayer short B–Te ES',color='#2f6f8f'); ax_es.bar(xx+.5*w,vals_long,w,label='Intralayer long B–Te ES',color='#8db6c8'); ax_es.bar(xx+1.5*w,vals_gap,w,label='Interlayer Te–Te ES',color='#d58c3a'); ax_es.set_xticks(xx,M,rotation=18,ha='right'); ax_es.set_ylabel('ES (e per contact)'); ax_es.set_title('(a) Electron sharing, ES',loc='left'); ax_es.legend(frameon=False,fontsize=8); ax_es.grid(axis='y',alpha=.18)
 ax2=axs_es[1]; ax2.bar(xx-.5*w,vals_et,w,label='Interlayer Te |ET|',color='#6b9e63'); ax2.bar(xx+.5*w,vals_lp,w,label='ELF f=0.85 basin electrons',color='#b36b9e'); ax2.set_xticks(xx,M,rotation=18,ha='right'); ax2.set_ylabel('Electrons (e; definitions differ)'); ax2.set_title('(b) Interlayer ET and ELF localization',loc='left'); ax2.legend(frameon=False,fontsize=8); ax2.grid(axis='y',alpha=.18)
 fig_es.suptitle('ES, ET, and ELF under the Raty et al. definitions',fontsize=12)
 for ext in ['png','pdf','svg']: fig_es.savefig(D/'02_图表'/('ES_ET_ELF_文献口径.'+ext),dpi=300)
 plt.close(fig_es)
 H.append(fig('ES_ET_ELF_文献口径.png','图 4b　ES–ET–ELF 双面板对照。(a) 三类接触的逐接触 ES=2×DI；(b) 层间 Te 的 |ET| 与 f=0.85 ELF 最大盆电子数。两面板均以电子数为显示单位，但 ET 是 Bader 电荷转移、ELF 是局域化盆积分，物理定义不同，不能相减或直接相除。图 (a) 显示短 B–Te ES 明显高于层间 Te–Te；图 (b) 显示层间 ET 与局域化电子数在四种材料中保持同量级。'))
 H.append('<p>口径边界（必须随数值一起声明）：(1) 波函数是 <b>PAW 赝波函数</b>（Sigma|C|^2 在 0.35–1.30 之间），占据轨道先做 Löwdin 正交化，否则求和规则偏 21%；34 个占据带重建密度 65.04 e 对 CHGCAR 68.00 e，差 4.4% 即原子球内增补项。(2) 胞内原子对的 DI 是<b>折叠值</b> Sigma_R DI(A,B+R)；单接触 ES = 2×折叠 DI / 对称重数 m，严格逐周期像分解需要 MLWF（文献路线）。(3) k 采样为 38 个 IBZ 点 + 空间群展开，SOC 与非 SOC 的差异未单独评估。文献原文已归档于 '+source('06_COGITO成键分析与COHP/源数据/文献方法/','Raty 2021 ES/ET pair-density 方法')+'，本研究结果与脚本归档于 '+source('06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI/','DI/ET/ELF 盆积分结果与方法说明')+'。</p>')
 H.append('<p>短长 B–Te 的 ICOHP/ICOOP 分裂、Bader 电荷再分配和 ELF 局域化形状分别回答“成键贡献是否分化”“电子怎样重新分配”“局域电子是否形成空间瓣”。它们与 ES 和 ET 构成相互约束的证据，而非可相加的孤对活性分数。</p>')
 rows=[]
 for m in M:
  b=min([r for r in bond[m]['bonds'] if r['pair']=='Te-Te'],key=lambda r:r['d_A']);rows.append([m,b['d_A'],f"{b['ICOHP_eV_bond']:+.4f}",f"{b['ICOOP_e_bond']:+.4f}",b['count']])
 H.append(table(['材料','最近 Te–Te d / Å','ICOHP / eV·bond⁻¹','ICOOP / e·bond⁻¹','每胞重数'],rows))
 H.append('<h3>4.5 ELF、COHP 与带边轨道的交叉约束</h3><p>ELF 的高值区反映局域化，不等于该处积累了多少电子；原子核邻近与阴离子价电子也会产生高 ELF。以下使用相同色标的 (110) 切面作定性对照。是否存在孤对定向瓣，应对齐原子位置、切面和等值标准检查；仅凭某个元素周围更红，不能给出孤对活性排序。图像不取代 ELFCAR 的三维盆积分或差分电荷。</p><div class="grid">')
 for m in M:H.append(fig('elf_'+m+'.png',m+'：归档 ELF 切面。颜色表示电子局域函数的相对高低，需在相同色标、相同晶面和原子位置条件下比较局域区域形状；它提供定性空间证据，不是孤对电子数或电荷转移量。'))
 H.append('</div>')
 H.append(fig('COHP_四材料总曲线_审核重绘.png','图 5　从归档交互式 HTML 中的 total COHP 数组重新作图。横轴保留源幅值，未假定每胞归一化；纵轴为 E−EF。蓝色负值/红色正值对应成键/反键，水平虚线是 EF，竖线是 COHP=0。全范围展示，没有因横轴截断遗漏大幅正峰。'))
 H.append('<p>占据区出现正 COHP，说明选定键集合存在已占据反键贡献；但总曲线并未逐独立核验键集合与重数，不能把某个峰唯一指定为 Sn-s/Te-p 反键。轨道分辨积分尚未闭合到完整 Bi-d+s/p 分量，因此不用它支撑唯一机制结论。保留原始交互曲线供后续同键、同归一化的分解积分。短长键的净 ICOHP 差则已由原始键表独立复核。</p>')
 H.append(table(['材料','charge spilling %','平均 orbital mixing %','最大 orbital mixing %'],[[m,bond[m]['quality']['charge_spilling_percent'],bond[m]['quality']['avg_orbital_mixing_percent'],bond[m]['quality']['max_orbital_mixing_percent']] for m in M]))
 H.append('<p>上述投影诊断较小，但不能替代基组、k 网格、能带和键截断收敛。Mulliken 总人口恢复约 68/88 个价电子，on-site 总人口不守恒到相同数值；因此 on-site s 大于 2 并不代表孤对获得了额外自由电子。综合本节证据，层内短 B–Te 的共享和成键贡献最强；层间 Te–Te 具有非零 ES、非零 ET 及可积分的 ELF 局域化，说明层间耦合不可按理想零共享间隙处理。孤对相关判断只有在 ELF 局域化形状、层间 ET/ES 和带边 s/p 混合同时一致时才成立，单一图像或单一电荷指标不足以定义“孤对活性”。</p></section>')
 H.append('<section id="s5"><h2>5　缺陷：化学势如何改变最低结构与电荷补偿</h2><p>这里从完美宿主转到缺陷能量与电荷平衡。沿用原始 81 原子能量模型及其已存 E<sub>F,eq</sub>，本轮重新计算全部 18 点的最低/次低缺陷、电荷态、能量差和 p/n 比，均与原表显示精度一致；未重新求解缺陷平衡。</p>')
 template=BeautifulSoup((templates/'缺陷章节.html').read_text(),'html.parser')
 for row in template.find_all('tr'):
  cells=row.find_all('td')
  if cells and cells[0].get_text(' ',strip=True)=='带边与修正':
   cells[1].clear();cells[1].append('HSE06+SOC 带边参考；形成能源 JSON 保留其源计算的统一修正记录。')
   cells[2].clear();cells[2].append('本报告沿用缺陷源模型记录，只复核排序、相对能量差和 p/n 判型；不将形成能重新参数化为孤对电子效应。')
 for cap in template.find_all('figcaption'):
  cap.append(' 图中的形成能、平衡费米能级和载流子均按同一化学势顶点读取；其数值沿用缺陷源数据记录。')
 for im in template.find_all('img'):
  p=D/im['src'];used.append(str(p.relative_to(D)));im['src']='data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()
 H.append(str(template));H.append('</section><section id="s6"><h2>6　输运：将固定缺陷净浓度与 SOC 对照接起来</h2><p>输运输入的正浓度表示净空穴、负浓度表示净电子；默认扫描不是本征导电类型预测。缺陷对应点取原模型 p−n，经已有输入精度舍入后固定扫温度，不在每个温度重新求解化学平衡。σ、S、κe 和 μ 的“平均”均为张量对角迹/3；PF 先按各方向 σᵢSᵢ² 计算，再平均，单位换算为 10⁻⁹（S 以 μV/K、σ 以 S/m）。这些方向平均是报告定义，非多晶有效介质结果。</p>')
 H.append('<h3>6.1 Sn₂Bi₂Te₅：原主报告表 10 的复核与 SOC 扩展</h3>')
 rows=[]
 for v in sn:
  for variant in ['noSOC','SOC']:
   z=v[variant];rows.append([f"v{v['vertex']}",'Bi_Te' if v['vertex']==0 else 'Sn_Bi',variant,f"{z['input_doping']:.6e}",f"{z['sigma300']:.6e}",f"{z['S300']:.6f}",f"{z['PF300']:.6f}"])
 H.append(table(['条件','最低结构标签','输运方案','净掺杂 cm⁻³','σ均值 S/m','S均值 μV/K','PF均值 mW/m/K²'],rows,'snbi-transport-table'))
 H.append('<p>基线四个 300 K S 为 +107.51、+102.97、+115.77、+186.82 μV/K，PF 为 2.22、5.21、1.41、1.39 mW/(m·K²)，与已修订主表一致。加入 SOC 波函数与形变势后，v0、v2、v3 的 S 和 PF 差异明显；例如 v0 的 PF 从约 2.22 增至 7.83。这个比较同时改变了能带与形变势等输入，不是将孤对贡献单独打开或关闭的实验。两套方案使用同一当前收敛介电。</p>')
 H.append(fig('SnBi_缺陷浓度_SOC输运.png','图 6　Sn₂Bi₂Te₅ 四个原缺陷净浓度条件下的三方向平均输运。横轴是温度，左、中、右分别为 Seebeck 系数、电导率和功率因子；实线 SOC、虚线基线，颜色表示化学势顶点。比较同色实虚线可读取 SOC 输入改变的响应；两方案共用当前介电，浓度固定，未重新求解温度依赖的缺陷平衡。'))
 H.append('<p>基线 v3 的 S 在 750/800/900 K 分别约 +13.90/−0.74/−17.57 μV/K，800 K 的 x/y 约 −9.89、z 约 +17.55 μV/K；v0–v2 在全部 100–900 K 采样点为正。净空穴输入下的 S 变号可以来自电子/空穴输运加权竞争，不等于缺陷突然转为施主。SOC 下的温度曲线应独立读取，不将基线的变号温度照搬过去。高浓度的 σ/PF 并不能反向验证负形成能或稀释失效的缺陷浓度。</p>')
 H.append('<h3>6.2 四材料统一 ±10¹⁹ cm⁻³ 的比较</h3>')
 rows=[]
 for m in M:
  for sign in [-1,1]:
   vals=[]
   for variant in ['noSOC','SOC']:
    z=tr[m+'_'+variant];di=next(i for i,d in enumerate(z['doping']) if np.isclose(d,sign*1e19,rtol=1e-10,atol=0));ti=z['temperatures'].index(300);vals.append({k:z['mean'][k][di][ti] for k in ['conductivity','seebeck','PF']})
   rows.append([m,'n' if sign<0 else 'p',*[f"{vals[0][k]:.4g} → {vals[1][k]:.4g}" for k in ['conductivity','seebeck','PF']]])
 H.append(table(['材料','净掺杂类型','σ S/m：基线→SOC','S μV/K：基线→SOC','PF mW/m/K²：基线→SOC'],rows,'soc-table'))
 H.append('<p>统一浓度比较回答“同一净载流子输入下，电子结构方案怎样影响输运”；缺陷条件比较回答“采用原模型不同净浓度时会得到什么响应”。两种横向比较不可混淆。两套输运方案共用当前收敛介电。'+source('05_SOC_AMSET补算/06_SOC_无SOC对比/SOC_无SOC对比报告.html','完整 SOC 对比报告（8 图及设置差异）')+'；'+source('04_AMSET输运数据/整体输运报告_含平均.html','全量基线输运报告')+'。</p></section>')
 H.append('<section id="s7"><h2>7　如何把课题讲成一条可检验的机制链</h2><p><b>第一层：宿主成键的非均匀性。</b>B–Te 的 3+3 短长配位与 ICOHP/ICOOP 分化构成直接证据，ELF 提供空间局域化背景。它们支持局部成键重排，但尚不能拆分尺寸效应、孤对杂化与晶体整体能量最小化各自的作用。</p><p><b>第二层：这种背景投射到哪些近边态。</b>能窗分析表明 s/p 混合具有材料与 k 点选择性。Sn₂Bi₂Te₅ 的近价带具有更明显 Sn-s，而 Sn₂Sb₂Te₅ 的精确 VBM 更偏 p 成分；宽窗会缩小差异，故论文应同时展示带边极值和有限窗，不只引用一个高百分比。</p><p><b>第三层：用缺陷能量解释导电类型。</b>反位改变原子数、电子数和局部成键，能量差还含化学势与弛豫。现有结果显示最低缺陷随化学势切换；在原电荷平衡模型中三种材料各采样点偏 p，Pb₂Bi₂Te₅ 有一个强补偿的 n 型候选点。这一步的证据来自缺陷能量和电荷平衡，不来自 TET 或 ELF 颜色。</p><p><b>第四层：在给定载流子条件下讨论输运。</b>当前收敛介电提供统一的屏蔽输入；SOC 比较显示相同浓度也可有明显不同 S、σ 和 PF。HSE+SOC 的轨道图与 PBE+SOC/无 SOC 输运不是同一个哈密顿量，只能作定性关联，不能据轨道百分比划分电导。</p>')
 H.append(table(['材料','当前最有信息的组合观察','适合进一步验证的问题'],[
 ['Pb₂Sb₂Te₅','最小路径带隙；VBM 阳离子 s 合计对窗口较稳定；所列缺陷点均偏 p。','用更密介电/输运网格量化双极效应与方向性，检验离子介电残差。'],
 ['Pb₂Bi₂Te₅','VBM/CBM 成分随窗口变化；电荷再分配描述量较大；v3 强补偿候选 n。','充分弛豫带电缺陷并统一校正，测试 Bi_Pb/Pb_Bi 微小能差与 n 型稳健性。'],
 ['Sn₂Sb₂Te₅','精确 VBM 的 s 较小，但宽窗 s 增强；B–Te 的 ICOHP 短长比大。','分析缺陷引起的局部键变化，而非从 Sn 元素本身预设其孤对活性。'],
 ['Sn₂Bi₂Te₅','VBM 及 50 meV 窗口均有明显 Sn-s；SOC 输运与基线差异大。','在同一 SOC 方法下对 Sn_Bi 缺陷/完美胞做轨道分辨 COHP 与差分电荷，闭合能量机制。']]))
 H.append('<p class="lead">适合论文的当前表述：A₂B₂Te₅ 的组分替换伴随局部 B–Te 成键分化、材料依赖的带边 s/p 混合及不同屏蔽与输运响应。带边轨道结论对能量窗口存在可量化敏感性。导电类型在现有缺陷模型中由化学势控制的施主—受主竞争决定；“孤对电子导致特定反位更易形成”仍是需要缺陷/完美结构一致方法对照验证的机制假设。</p></section>')
 H.append('<section id="s8"><h2>8　审查结论、复现与材料边界</h2><p>已完成：四套当前收敛 DFPT 张量对两套 AMSET 输入；8 份输运数据的维度、有限性与源头/镜像一致性；3400 行 SOC 对比和 8 行定点汇总；18 个缺陷顶点；402 条 COGITO 键记录；四材料路径 PROCAR 的电荷投影、周期去重、0.1–100 meV 多窗口统计。最终审查记录另核对文件 SHA256、报告数表、图片解码及链接。</p><p>文献口径的 ES = 2×DI（SOC 波函数 + 全密度 Bader 盆，四材料，求和规则偏差 0.11%–0.35%）、全密度 Bader 的逐原子 ET 和 ELF = 0.80/0.85/0.90 等值面盆积分均已归档。当前结论限于 PAW 赝波函数、38 个 IBZ 点加空间群展开和胞内折叠 DI；完整三维谷简并度、缺陷—完美结构同方法能量分解及 ES 的逐周期像 MLWF 分解不在本报告范围内。</p>')
 H.append('<p>'+source('09_审查/最终审查报告.html','最终审查报告')+' · '+source('09_审查/文件SHA256.json','交付 SHA256 清单')+' · '+source('10_脚本/README.md','脚本与复现方法')+' · '+source('01_报告/结构对应报告.html','结构对应报告')+'</p>')
 H.append('<p>方法依据：<a href="https://vasp.at/wiki/PROCAR">VASP PROCAR</a>、<a href="https://vasp.at/wiki/LORBIT">LORBIT</a>说明投影与 SOC 附加分量；<a href="https://pypi.org/project/cogito-dft/0.3.13/">COGITO 0.3.13</a>源码用于核对 ICOOP 输出字段，未据此声称远程安装版本逐字节相同；<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10853697/">Metavalent or Hypervalent Bonding: Is There a Chance for Reconciliation?</a>用于比较电子共享/转移定义，未将本项目直接归类为某种键型。</p>')
 H.append((templates/'缺陷参考文献.html').read_text());H.append('</section>')
 css='body{max-width:1160px;margin:32px auto;padding:0 28px;color:#243747;font:16px/1.85 "Microsoft YaHei",sans-serif}h1{font-size:32px;color:#123c59}h2{margin-top:50px;border-bottom:2px solid #b7d0df;padding-bottom:10px;color:#164d6d}h3{margin-top:28px;color:#24607e}p{margin:15px 0}.subtitle,figcaption{color:#61717e}.lead,.note{background:#eef5f8;border-left:4px solid #33839d;padding:16px 22px}nav{display:flex;flex-wrap:wrap;gap:12px 22px;background:#f3f7fa;padding:16px}a{color:#216d9b}table{border-collapse:collapse;width:100%;font-size:13px;font-variant-numeric:tabular-nums;margin:18px 0}td,th{border-bottom:1px solid #d2dee5;padding:9px 10px;text-align:left;vertical-align:top}th{background:#edf3f7}tr:nth-child(even){background:#fafcfd}.scroll,.tablewrap{overflow-x:auto;max-width:100%;box-sizing:border-box}.tablewrap{display:block}img{max-width:100%;height:auto;display:block;margin:auto}figure{margin:26px 0}figcaption{font-size:13px;line-height:1.65;margin-top:9px}.grid{display:grid;grid-template-columns:1fr 1fr;gap:24px}.grid img{max-height:560px;width:auto}section{scroll-margin-top:20px;max-width:100%;box-sizing:border-box}code,pre{overflow-wrap:anywhere;word-break:break-word}details{background:#f8fafb;padding:15px;margin:20px 0;max-width:100%;box-sizing:border-box}@media(max-width:700px){body{padding:0 12px;margin:16px auto;overflow-x:hidden}.grid{grid-template-columns:1fr}table{font-size:11px}}@media print{body{max-width:none;padding:0;font-size:11pt}nav{display:none}h2{break-after:avoid}figure,tr{break-inside:avoid}.scroll,.tablewrap{overflow:visible}}'
 final='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>A₂B₂Te₅ 综合分析与源数据复核</title><style>'+css+'</style></head><body>'+''.join(H)+'</body></html>'
 # The preserved defect chapter has its own historical numbering; references to the prior transport section are removed there only if present.
 from paper_revision import revise
 final=revise(final,D,used)
 oldpath.write_text(final,encoding='utf-8')
 (D/'09_审查/报告图片来源.json').write_text(json.dumps([{'path':p,'sha256':hashlib.sha256((D/p).read_bytes()).hexdigest()} for p in used],ensure_ascii=False,indent=2))
 print('Report generated:',len(final.encode()),'bytes;',len(used),'figures')
if __name__=='__main__':run(Path(sys.argv[1]))
