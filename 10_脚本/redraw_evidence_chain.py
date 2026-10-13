from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '02_图表'
OUT.mkdir(exist_ok=True)
MATS = ['Pb2Sb2Te5', 'Pb2Bi2Te5', 'Sn2Sb2Te5', 'Sn2Bi2Te5']
COLORS = ['#2878a0', '#d28d32', '#3b9b7b', '#9a6caf']

bands = json.loads((ROOT / '07_带边能窗/带边能窗统计.json').read_text(encoding='utf-8'))
bonds = json.loads((ROOT / '09_审查/成键复算.json').read_text(encoding='utf-8'))
defects = json.loads((ROOT / '09_审查/缺陷复算.json').read_text(encoding='utf-8'))
transport = json.loads((ROOT / '09_审查/输运复算.json').read_text(encoding='utf-8'))

ratios, edge_s, net, pf = [], [], [], []
for m in MATS:
    B = 'Bi' if 'Bi' in m else 'Sb'
    bb = sorted([x for x in bonds[m]['bonds'] if set(x['pair'].split('-')) == {B, 'Te'}], key=lambda x: x['d_A'])[:2]
    ratios.append(abs(bb[0]['ICOHP_eV_bond']) / abs(bb[1]['ICOHP_eV_bond']))
    w = next(x for x in bands[m]['windows'] if x['edge'] == 'VBM' and x['window_meV'] == 50)
    A = 'Sn' if m.startswith('Sn') else 'Pb'
    edge_s.append(w['percent'][A + '-s'])
    vals = np.array([x['net_p_minus_n'] for x in defects[m]], dtype=float)
    net.append(float(np.median(vals)))
    z = transport[m + '_SOC']
    ti = z['temperatures'].index(300)
    di = min(range(len(z['doping'])), key=lambda i: abs(z['doping'][i] - 1e19))
    pf.append(float(z['mean']['PF'][di][ti]))

plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.titleweight': 'bold', 'figure.dpi': 150})
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
fig.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.16, wspace=0.24, hspace=0.34)
x = np.arange(len(MATS))
labels = ['Pb₂Sb₂Te₅', 'Pb₂Bi₂Te₅', 'Sn₂Sb₂Te₅', 'Sn₂Bi₂Te₅']

def bar(a, vals, title, ylabel, fmt='{:.2f}', zero=False):
    bars = a.bar(x, vals, color=COLORS, width=0.66)
    a.set_xticks(x, labels, rotation=18, ha='right')
    a.set_title(title, loc='left')
    a.set_ylabel(ylabel)
    a.grid(axis='y', color='#dfe6ea', linewidth=.8)
    a.set_axisbelow(True)
    if zero: a.axhline(0, color='#555', linewidth=.9)
    for b, v in zip(bars, vals):
        a.text(b.get_x() + b.get_width()/2, v + (0.025 * (max(vals)-min(vals) or 1)), fmt.format(v), ha='center', va='bottom', fontsize=9)

bar(ax[0,0], ratios, '(a) Short/long B–Te bond strength', r'|ICOHP$_{short}$| / |ICOHP$_{long}$|')
bar(ax[0,1], edge_s, '(b) Cation-s character near the VBM', 'A-s projection in 50 meV window (%)', '{:.1f}')
bar(ax[1,0], np.array(net) / 1e19, '(c) Defect-model carrier balance', r'median (p−n) / 10$^{19}$ cm$^{-3}$', '{:.1f}', zero=True)
bar(ax[1,1], pf, '(d) SOC transport at 300 K', r'PF at +10$^{19}$ cm$^{-3}$ (mW m$^{-1}$ K$^{-2}$)', '{:.2f}')
fig.suptitle('Evidence chain for A₂B₂Te₅: bond asymmetry → band-edge mixing → defect balance → transport', fontsize=15)
fig.text(0.5, 0.035, 'Panels summarize separate calculations; they show material-dependent association, not a demonstrated causal law.', ha='center', fontsize=9, color='#53636d')
for ext in ['png', 'svg', 'pdf']:
    fig.savefig(OUT / f'证据链_成键到带边缺陷输运.{ext}', bbox_inches='tight')
plt.close(fig)
print('wrote evidence-chain figure')
