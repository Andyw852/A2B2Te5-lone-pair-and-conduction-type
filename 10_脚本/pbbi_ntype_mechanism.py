#!/usr/bin/env python3
"""Pb2Bi2Te5 n 型倾向的元素来源：长键、反位与阳离子空位（报告一节）。

只读取仓库中已有的计算输出，不运行任何 DFT：
  03_缺陷形成能与PN型数据/<材料>.zip  -> 计算结果/formation_energy_results.json、全部顶点计算结果.json
  06_COGITO成键分析与COHP/结果/<材料>/icohp_short_long.json
  06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI/DI结果/ES_逐接触汇总.json
  06_COGITO成键分析与COHP/源数据/Raty2021_ES_ET_DI/ET_全密度Bader汇总.json
  07_带边能窗/带边能窗统计.csv
  09_审查/共享与局域化审查.json

用法（以交付目录为当前目录）：
  python 10_脚本/pbbi_ntype_mechanism.py .

输出：13_PbBi_n型机制分析_20261125/ 下的 HTML 报告、独立 SVG 图和源数据 CSV/JSON。
只依赖 Python 标准库。
"""
import csv
import html
import json
import math
import sys
import zipfile
from pathlib import Path

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
OUT = ROOT / '13_PbBi_n型机制分析_20261125'
FIG = OUT / '图'
SRC = OUT / '源数据'

MATS = {  # 材料: (A, B)
    'Pb2Bi2Te5': ('Pb', 'Bi'),
    'Pb2Sb2Te5': ('Pb', 'Sb'),
    'Sn2Bi2Te5': ('Sn', 'Bi'),
    'Sn2Sb2Te5': ('Sn', 'Sb'),
}
FOCUS = 'Pb2Bi2Te5'
SUB = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')


def pretty(m):
    return m.translate(SUB)


def esc(s):
    return html.escape(str(s), quote=True)


# ---------------------------------------------------------------- 读取数据
def load_defects(m):
    with zipfile.ZipFile(ROOT / '03_缺陷形成能与PN型数据' / f'{m}.zip') as z:
        fe = json.loads(z.read('计算结果/formation_energy_results.json'))
        vs = json.loads(z.read('计算结果/全部顶点计算结果.json'))
    return fe, vs


def formation(d, q, mu, mu0, EF):
    """E_f(q, E_F) = Ef0[q] - sum_i dn_i (mu_i - mu0_i) + q E_F，Ef0 对应顶点 v0 的化学势、E_F = VBM。"""
    dmu = sum(d['dn'].get(e, 0) * (mu[e] - mu0[e]) for e in mu)
    return d['Ef0'][q] - dmu + int(q) * EF


def lowest_by_name(defs, name, mu, mu0, EF, qsign=None):
    """同名缺陷（所有不等价位点、所有电荷态）在给定 E_F 下的最低形成能。qsign: None/+1/-1 限定电荷符号。"""
    best = None
    for k, d in defs.items():
        if d['name'] != name:
            continue
        for q in d['Ef0']:
            qq = int(q)
            if qsign is not None and (qq == 0 or (qq > 0) != (qsign > 0)):
                continue
            e = formation(d, q, mu, mu0, EF)
            if best is None or e < best[0]:
                best = (e, qq, k)
    return best


def envelope(defs, mu, mu0, EF, sign, exclude=()):
    """带正电（sign>0）或带负电（sign<0）缺陷的最低形成能及其标识。"""
    best = (math.inf, None, None)
    for k, d in defs.items():
        if d['name'] in exclude:
            continue
        for q in d['Ef0']:
            qq = int(q)
            if qq == 0 or (qq > 0) != (sign > 0):
                continue
            e = formation(d, q, mu, mu0, EF)
            if e < best[0]:
                best = (e, qq, d['name'])
    return best


def pinning_level(defs, mu, mu0):
    """最低带正电缺陷 = 最低带负电缺陷 的费米能级（相对 VBM）。差值随 E_F 单调增，二分求解。"""
    f = lambda x: envelope(defs, mu, mu0, x, +1)[0] - envelope(defs, mu, mu0, x, -1)[0]
    lo, hi = -2.0, 2.0
    assert f(lo) < 0 < f(hi)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    x = 0.5 * (lo + hi)
    return x, envelope(defs, mu, mu0, x, +1), envelope(defs, mu, mu0, x, -1)


def main():
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    SRC.mkdir(exist_ok=True)

    bond_dir = ROOT / '06_COGITO成键分析与COHP' / '结果'
    raty = ROOT / '06_COGITO成键分析与COHP' / '源数据' / 'Raty2021_ES_ET_DI'
    es_all = json.loads((raty / 'DI结果' / 'ES_逐接触汇总.json').read_text(encoding='utf-8'))
    et_all = json.loads((raty / 'ET_全密度Bader汇总.json').read_text(encoding='utf-8'))
    share_audit = json.loads((ROOT / '09_审查' / '共享与局域化审查.json').read_text(encoding='utf-8'))
    with open(ROOT / '07_带边能窗' / '带边能窗统计.csv', encoding='utf-8-sig') as fh:
        edge_rows = list(csv.DictReader(fh))

    R = {}  # 每个材料的汇总
    checks = []
    vertex_rows = []
    for m, (A, B) in MATS.items():
        fe, vs = load_defects(m)
        defs = fe['defects']
        Eg = fe['E_gap']
        mu0 = vs[0]['chemical_potentials_eV_per_atom']
        r = {'A': A, 'B': B, 'Eg': Eg, 'vertices': []}

        # 复核：重建的形成能应与源记录的最低缺陷一致
        for vx in vs:
            mu = vx['chemical_potentials_eV_per_atom']
            EFeq = vx['E_F_eq_relative_VBM_eV']
            md = vx['minimum_defect']
            rec = formation(defs[md['defect_id']], str(md['q_min']) if md['q_min'] <= 0 else f"+{md['q_min']}", mu, mu0, EFeq)
            checks.append(abs(rec - md['formation_energy_eV']))

            x, don, acc = pinning_level(defs, mu, mu0)
            mid = Eg / 2
            pin_type = 'n' if x > mid else 'p'
            row = {
                'material': m, 'vertex': vx['vertex'], 'mu': mu,
                'EF_eq_model': EFeq, 'type_model': vx['carrier_type_model'],
                'EF_pin': x, 'EF_pin_minus_CBM': x - Eg, 'type_pin': pin_type,
                'donor': f'{don[2]}({don[1]:+d})', 'acceptor': f'{acc[2]}({acc[1]:+d})',
                'E_at_pin': don[0],
                'min_defect_model': f"{md['name']}(q={md['q_min']:+d})",
                'min_defect_Ef': md['formation_energy_eV'],
            }
            row['n_cm3_model'], row['p_cm3_model'] = vx['n_cm3'], vx['p_cm3']
            row['ionized_model'] = {c['defect']: c['ionized_cm3'] for c in vx['source_vertex_record'].get('top_concentration', [])}
            row['BA_plus_at_eq'] = lowest_by_name(defs, f'{B}_{A}', mu, mu0, EFeq, +1)[0]
            row['AB_minus_at_eq'] = lowest_by_name(defs, f'{A}_{B}', mu, mu0, EFeq, -1)[0]
            # 在 E_F = CBM 处比较三类关键缺陷
            row['BA_plus_at_CBM'] = lowest_by_name(defs, f'{B}_{A}', mu, mu0, Eg, +1)[0]
            row['AB_minus_at_CBM'] = lowest_by_name(defs, f'{A}_{B}', mu, mu0, Eg, -1)[0]
            vA = lowest_by_name(defs, f'v_{A}', mu, mu0, Eg, -1)
            row['vA_minus_at_CBM'] = vA[0]
            row['vA_q_at_CBM'] = vA[1]
            acc_cbm = envelope(defs, mu, mu0, Eg, -1)
            don_cbm = envelope(defs, mu, mu0, Eg, +1)
            row['lowest_acceptor_at_CBM'] = f'{acc_cbm[2]}({acc_cbm[1]:+d})'
            row['lowest_acceptor_E_at_CBM'] = acc_cbm[0]
            row['lowest_donor_at_CBM'] = f'{don_cbm[2]}({don_cbm[1]:+d})'
            row['lowest_donor_E_at_CBM'] = don_cbm[0]
            r['vertices'].append(row)
            vertex_rows.append(row)

        best = max(r['vertices'], key=lambda t: t['EF_pin'])
        r['best_vertex'] = best
        r['n_vertices'] = len(vs)
        r['n_type_model'] = sum(1 for t in r['vertices'] if t['type_model'] == 'n')
        r['agree'] = sum(1 for t in r['vertices'] if t['type_model'] == t['type_pin'])

        # 与化学势无关的反位交换能
        def ef0(name, q):
            c = [d['Ef0'][q] for d in defs.values() if d['name'] == name and q in d['Ef0']]
            return min(c)
        r['exchange'] = ef0(f'{B}_{A}', '0') + ef0(f'{A}_{B}', '0')
        r['ionized_pair'] = ef0(f'{B}_{A}', '+1') + ef0(f'{A}_{B}', '-1')
        pair = [d for d in defs.values() if d['name'] == f'{B}_{A}__{A}_{B}_pair']
        r['bound_pair'] = pair[0]['Ef0']['0'] if pair else None
        tBA = [d['transitions'] for d in defs.values() if d['name'] == f'{B}_{A}'][0]
        tAB = [d['transitions'] for d in defs.values() if d['name'] == f'{A}_{B}'][0]
        r['BA_level'] = tBA.get('(0/+1)')
        r['AB_level'] = tAB.get('(-1/0)')

        # 成键：COGITO 键长与 ICOHP
        sl = json.loads((bond_dir / m / 'icohp_short_long.json').read_text(encoding='utf-8'))
        for p in sl:
            el = set(p['pair'].split('-'))
            if el == {'Te', A}:
                s, l = p['short'], p['long']
                r['A_Te'] = {'d_s': s['dist_ang'], 'd_l': l['dist_ang'], 'icohp_s': s['icohp_ev_per_bond'],
                             'icohp_l': l['icohp_ev_per_bond'], 'icoop_s': s['icoop_elec_per_bond'], 'icoop_l': l['icoop_elec_per_bond']}
                # 长于 4 Å 的“长键”属于远邻（Sn2Sb2Te5 的 5.36 Å），只取短键作为配位壳层
                r['A_Te']['single_shell'] = l['dist_ang'] > 4.0
            elif el == {'Te', B}:
                s, l = p['short'], p['long']
                r['B_Te'] = {'d_s': s['dist_ang'], 'd_l': l['dist_ang'], 'icohp_s': s['icohp_ev_per_bond'],
                             'icohp_l': l['icohp_ev_per_bond'], 'icoop_s': s['icoop_elec_per_bond'], 'icoop_l': l['icoop_elec_per_bond']}
        a = r['A_Te']
        r['dA'] = a['d_s'] if a['single_shell'] else 0.5 * (a['d_s'] + a['d_l'])
        r['dB'] = 0.5 * (r['B_Te']['d_s'] + r['B_Te']['d_l'])
        r['mismatch'] = r['dA'] - r['dB']

        # ES（Raty 2021 口径的待验证估计）与 ET（全密度 Bader）
        cls = [c for c in es_all[m]['classes'] if c['layer'] == '层内' and c['d'] < 3.5]
        r['ES_A'] = [(c['d'], c['ES']) for c in cls if set(c['pair'].split('-')) == {A, 'Te'}]
        bt = sorted([(c['d'], c['ES']) for c in cls if set(c['pair'].split('-')) == {B, 'Te'}])
        r['ES_B_short'], r['ES_B_long'] = bt[0], bt[-1]
        q = et_all[m]['per_element_mean_q_ET']
        r['qA'], r['qB'], r['qTe'] = q[A], q[B], q['Te']

        # 带边轨道组成
        def edge(e, w, col):
            for row in edge_rows:
                if row['material'] == m and row['edge'] == e and row['window_meV'] == w:
                    return float(row[col] or 0)
        r['VBM_As_01'] = edge('VBM', '0.1', f'{A}-s')
        r['VBM_As_50'] = edge('VBM', '50', f'{A}-s')
        r['CBM_Bp_01'] = edge('CBM', '0.1', f'{B}-p')
        r['VBM_Tep_01'] = edge('VBM', '0.1', 'Te-p')
        R[m] = r

    max_check = max(checks)
    assert max_check < 1e-4, max_check
    agree_total = sum(r['agree'] for r in R.values())
    n_total = sum(r['n_vertices'] for r in R.values())

    # v3 形成能曲线（Pb2Bi2Te5 与 Sn2Bi2Te5 的最有利 n 型顶点）
    curves = {}
    for m in ['Pb2Bi2Te5', 'Sn2Bi2Te5']:
        A, B = MATS[m]
        fe, vs = load_defects(m)
        defs = fe['defects']
        mu0 = vs[0]['chemical_potentials_eV_per_atom']
        bv = R[m]['best_vertex']
        mu = bv['mu']
        xs = [round(-0.15 + 0.0025 * i, 4) for i in range(int(0.5 / 0.0025) + 1)]
        names = {'BA': f'{B}_{A}', 'AB': f'{A}_{B}', 'vA': f'v_{A}'}
        cv = {key: [] for key in ['BA', 'AB', 'vA', 'other']}
        for x in xs:
            for key, nm in names.items():
                b = lowest_by_name(defs, nm, mu, mu0, x)
                cv[key].append((x, b[0], b[1]))
            # 其余缺陷（含中性态）的最低包络
            bo = None
            for k, d in defs.items():
                if d['name'] in names.values():
                    continue
                for qq in d['Ef0']:
                    e = formation(d, qq, mu, mu0, x)
                    if bo is None or e < bo[0]:
                        bo = (e, int(qq), d['name'])
            cv['other'].append((x, bo[0], bo[1], bo[2]))
        curves[m] = {'vertex': bv['vertex'], 'Eg': R[m]['Eg'], 'EF_pin': bv['EF_pin'], 'EF_eq': bv['EF_eq_model'], 'cv': cv, 'names': names}

    write_sources(R, vertex_rows, curves, max_check, agree_total, n_total)
    figs = {
        'fig1': fig_pinning(R),
        'fig2': fig_bonds(R),
        'fig3': fig_exchange(R),
        'fig4': fig_vbm_s(R),
        'fig5': fig_cbm_defects(R),
        'fig6': fig_curves(curves),
    }
    for k, (svg, name) in figs.items():
        (FIG / name).write_text(standalone(svg), encoding='utf-8')
    page = build_html(R, vertex_rows, curves, figs, max_check, agree_total, n_total, share_audit)
    (OUT / 'PbBi_n型机制分析.html').write_text(page, encoding='utf-8')
    print('written', OUT)


# ---------------------------------------------------------------- 源数据
def write_sources(R, vertex_rows, curves, max_check, agree_total, n_total):
    with open(SRC / '钉扎费米能级_全部顶点.csv', 'w', newline='', encoding='utf-8-sig') as fh:
        cols = ['material', 'vertex', 'type_model', 'EF_eq_model', 'EF_pin', 'EF_pin_minus_CBM', 'type_pin', 'donor', 'acceptor',
                'E_at_pin', 'min_defect_model', 'min_defect_Ef', 'BA_plus_at_CBM', 'AB_minus_at_CBM', 'vA_minus_at_CBM', 'vA_q_at_CBM',
                'lowest_donor_at_CBM', 'lowest_donor_E_at_CBM', 'lowest_acceptor_at_CBM', 'lowest_acceptor_E_at_CBM', 'mu']
        w = csv.writer(fh)
        w.writerow(cols)
        for r in vertex_rows:
            w.writerow([json.dumps(r[c], ensure_ascii=False) if c == 'mu' else (f'{r[c]:.6f}' if isinstance(r[c], float) else r[c]) for c in cols])
    with open(SRC / '成键与反位交换.csv', 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh)
        w.writerow(['material', 'A', 'B', 'A-Te_d_short', 'A-Te_d_long', 'A-Te_ICOHP_short', 'A-Te_ICOHP_long', 'B-Te_d_short', 'B-Te_d_long',
                    'B-Te_ICOHP_short', 'B-Te_ICOHP_long', 'mean_d_A', 'mean_d_B', 'd_A_minus_d_B', 'ES_A-Te', 'ES_B-Te_short', 'ES_B-Te_long',
                    'q_ET_A', 'q_ET_B', 'q_ET_Te', 'exchange_E_BA0+AB0', 'ionized_pair_E_BA+ + AB-', 'bound_pair_E', 'BA_(0/+1)_level', 'AB_(-1/0)_level', 'E_gap',
                    'VBM_A-s_%_0.1meV', 'VBM_A-s_%_50meV', 'CBM_B-p_%_0.1meV'])
        for m, r in R.items():
            a, b = r['A_Te'], r['B_Te']
            esA = sum(e for _, e in r['ES_A']) / len(r['ES_A'])
            w.writerow([m, r['A'], r['B'], a['d_s'], '' if a['single_shell'] else a['d_l'], a['icohp_s'], '' if a['single_shell'] else a['icohp_l'],
                        b['d_s'], b['d_l'], b['icohp_s'], b['icohp_l'], f"{r['dA']:.4f}", f"{r['dB']:.4f}", f"{r['mismatch']:.4f}",
                        f'{esA:.4f}', f"{r['ES_B_short'][1]:.4f}", f"{r['ES_B_long'][1]:.4f}", r['qA'], r['qB'], r['qTe'],
                        f"{r['exchange']:.5f}", f"{r['ionized_pair']:.5f}", '' if r['bound_pair'] is None else r['bound_pair'],
                        r['BA_level'], r['AB_level'], r['Eg'], f"{r['VBM_As_01']:.2f}", f"{r['VBM_As_50']:.2f}", f"{r['CBM_Bp_01']:.2f}"])
    with open(SRC / 'n型最有利顶点_形成能曲线.csv', 'w', newline='', encoding='utf-8-sig') as fh:
        w = csv.writer(fh)
        w.writerow(['material', 'vertex', 'EF_rel_VBM', 'curve', 'defect', 'Ef_eV', 'q'])
        for m, c in curves.items():
            for key, pts in c['cv'].items():
                nm = c['names'].get(key, '其余缺陷包络')
                for p in pts:
                    w.writerow([m, c['vertex'], p[0], key, p[3] if key == 'other' else nm, f'{p[1]:.6f}', p[2]])
    summary = {
        'generated_by': '10_脚本/pbbi_ntype_mechanism.py',
        'note': '只读取仓库已有输出，未运行新的 DFT；E_F* 为本节后处理判据，不是完整电荷中性求解。',
        'reconstruction_max_abs_error_eV': max_check,
        'pinning_vs_model_type_agreement': f'{agree_total}/{n_total}',
        'materials': {m: {k: v for k, v in r.items() if k not in ('vertices',)} for m, r in R.items()},
    }
    (SRC / '汇总.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding='utf-8')


# ---------------------------------------------------------------- SVG 工具
FONT = 'system-ui,-apple-system,"PingFang SC","Microsoft YaHei","Noto Sans CJK SC","WenQuanYi Zen Hei","Segoe UI",sans-serif'
TOKENS_LIGHT = ('--surface:#ffffff;--ink:#191919;--ink2:#666666;--grid:#e7e7e7;--axis:#c2c2c2;--band:#f1f1f1;'
                '--s1:#257de0;--s2:#f47230;--s3:#00b894;--hl:#6a4ce2;--mute:#a3a3a3;')
TOKENS_DARK = ('--surface:#191919;--ink:#efefef;--ink2:#9c9c9c;--grid:#303030;--axis:#474747;--band:#242424;'
               '--s1:#338ced;--s2:#e76624;--s3:#00a888;--hl:#9887ef;--mute:#6f6f6f;')
SVG_CSS = """
.viz{font-family:%s;}
.viz text{fill:var(--ink2);font-size:12px;}
.viz .ink{fill:var(--ink);}
.viz .t-b{font-weight:600;}
.viz .small{font-size:11px;}
.viz .grid{stroke:var(--grid);stroke-width:1;}
.viz .axis{stroke:var(--axis);stroke-width:1;}
.viz .ref{stroke:var(--ink2);stroke-width:1;}
.viz .band{fill:var(--band);}
.viz .bg{fill:var(--surface);}
.viz .f-s1{fill:var(--s1);} .viz .f-s2{fill:var(--s2);} .viz .f-s3{fill:var(--s3);}
.viz .f-hl{fill:var(--hl);} .viz .f-mute{fill:var(--mute);} .viz .f-surf{fill:var(--surface);}
.viz .k-s1{stroke:var(--s1);} .viz .k-s2{stroke:var(--s2);} .viz .k-s3{stroke:var(--s3);}
.viz .k-hl{stroke:var(--hl);} .viz .k-mute{stroke:var(--mute);} .viz .k-ink2{stroke:var(--ink2);} .viz .k-surf{stroke:var(--surface);}
.viz .line{fill:none;stroke-width:2;stroke-linejoin:round;stroke-linecap:round;}
.viz .ext{stroke-dasharray:5 4;opacity:.75;}
.viz .hit{fill:transparent;stroke:none;}
.viz [data-tip]{cursor:default;}
""" % FONT


def standalone(svg):
    """独立 SVG 文件：嵌入浅色 token 与样式，可单独打开或插入论文草稿。"""
    style = '<style>svg{%s}%s</style>' % (TOKENS_LIGHT, SVG_CSS)
    return svg.replace('<svg ', '<svg xmlns="http://www.w3.org/2000/svg" ', 1).replace('>', '>' + style + '<rect class="bg" width="100%" height="100%"/>', 1)


def nice_ticks(lo, hi, n=5):
    span = hi - lo
    raw = span / n
    mag = 10 ** math.floor(math.log10(raw))
    step = min((s * mag for s in (1, 2, 2.5, 5, 10) if s * mag >= raw), default=10 * mag)
    t0 = math.ceil(lo / step - 1e-9) * step
    out = []
    t = t0
    while t <= hi + 1e-9:
        out.append(round(t, 10))
        t += step
    return out


def fmt_tick(v, step):
    dec = max(0, -int(math.floor(math.log10(step))) if step < 1 else 0)
    if abs(step * 10 ** dec - round(step * 10 ** dec)) > 1e-9:
        dec += 1
    if abs(v) < 1e-12:
        v = 0.0
    s = f'{v:.{dec}f}'
    return s.replace('-', '−')


class Axes:
    def __init__(self, x0, y0, w, h, xd, yd):
        self.x0, self.y0, self.w, self.h = x0, y0, w, h
        self.xd, self.yd = xd, yd

    def X(self, v):
        return self.x0 + (v - self.xd[0]) / (self.xd[1] - self.xd[0]) * self.w

    def Y(self, v):
        return self.y0 + self.h - (v - self.yd[0]) / (self.yd[1] - self.yd[0]) * self.h

    def frame(self, xt=None, yt=None, xlabel='', ylabel='', xgrid=False, ygrid=True):
        s = []
        if yt:
            step = yt[1] - yt[0] if len(yt) > 1 else 1
            for t in yt:
                y = self.Y(t)
                if ygrid:
                    s.append(f'<line class="grid" x1="{self.x0}" x2="{self.x0 + self.w}" y1="{y:.1f}" y2="{y:.1f}"/>')
                s.append(f'<text x="{self.x0 - 8}" y="{y + 4:.1f}" text-anchor="end" style="font-variant-numeric:tabular-nums">{fmt_tick(t, step)}</text>')
        if xt:
            step = xt[1] - xt[0] if len(xt) > 1 else 1
            for t in xt:
                x = self.X(t)
                if xgrid:
                    s.append(f'<line class="grid" x1="{x:.1f}" x2="{x:.1f}" y1="{self.y0}" y2="{self.y0 + self.h}"/>')
                s.append(f'<text x="{x:.1f}" y="{self.y0 + self.h + 18}" text-anchor="middle" style="font-variant-numeric:tabular-nums">{fmt_tick(t, step)}</text>')
        s.append(f'<line class="axis" x1="{self.x0}" x2="{self.x0 + self.w}" y1="{self.y0 + self.h}" y2="{self.y0 + self.h}"/>')
        if xlabel:
            s.append(f'<text x="{self.x0 + self.w / 2:.1f}" y="{self.y0 + self.h + 40}" text-anchor="middle">{xlabel}</text>')
        if ylabel:
            cx, cy = self.x0 - 52, self.y0 + self.h / 2
            s.append(f'<text transform="translate({cx:.1f},{cy:.1f}) rotate(-90)" text-anchor="middle">{ylabel}</text>')
        return ''.join(s)


def svg_open(w, h, label):
    return f'<svg class="viz" viewBox="0 0 {w} {h}" width="100%" role="img" aria-label="{esc(label)}" preserveAspectRatio="xMidYMid meet">'


def dot(x, y, r, cls, tip, ring=True, hollow=False):
    fill = 'f-surf' if hollow else cls.replace('k-', 'f-')
    stroke = cls if hollow else ('k-surf' if ring else '')
    sw = 2
    return (f'<g data-tip="{esc(tip)}"><circle class="hit" cx="{x:.1f}" cy="{y:.1f}" r="{r + 7}"/>'
            f'<circle class="{fill} {stroke}" cx="{x:.1f}" cy="{y:.1f}" r="{r}" stroke-width="{sw}"/></g>')


def legend_item(x, y, kind, cls, text):
    if kind == 'dot':
        g = f'<circle class="{cls}" cx="{x + 6}" cy="{y - 4}" r="5"/>'
    elif kind == 'hollow':
        g = f'<circle class="f-surf {cls}" stroke-width="2" cx="{x + 6}" cy="{y - 4}" r="5"/>'
    elif kind == 'line':
        g = f'<line class="line {cls}" x1="{x}" x2="{x + 18}" y1="{y - 4}" y2="{y - 4}"/>'
    elif kind == 'dash':
        g = f'<line class="line ext {cls}" x1="{x}" x2="{x + 18}" y1="{y - 4}" y2="{y - 4}"/>'
    elif kind == 'band':
        g = f'<rect class="band" x="{x}" y="{y - 10}" width="18" height="12" rx="2"/>'
    else:
        g = f'<rect class="{cls}" x="{x}" y="{y - 10}" width="12" height="12" rx="2"/>'
    off = 26 if kind in ('line', 'dash', 'band') else 18
    return g + f'<text x="{x + off}" y="{y}">{text}</text>'


# ---------------------------------------------------------------- 图 1：钉扎费米能级
def fig_pinning(R):
    W, H = 760, 330
    order = ['Pb2Bi2Te5', 'Pb2Sb2Te5', 'Sn2Bi2Te5', 'Sn2Sb2Te5']
    ax = Axes(130, 56, 590, 200, (-0.8, 0.2), (0, 4))
    s = [svg_open(W, H, '图1 各化学势顶点的钉扎费米能级相对导带底的位置')]
    xt = nice_ticks(-0.8, 0.2, 5)
    s.append(ax.frame(xt=xt, xlabel='E<tspan baseline-shift="sub" font-size="9">F</tspan>* − CBM（eV）；0 = 导带底，左侧负值越大越偏 p 型', ygrid=False, xgrid=True))
    rowh = ax.h / 4
    for i, m in enumerate(order):
        r = R[m]
        yc = ax.y0 + rowh * (i + 0.5)
        cls = 'k-hl' if m == FOCUS else 'k-mute'
        # 带隙区（−Eg..0）
        x1, x2 = ax.X(-r['Eg']), ax.X(0)
        s.append(f'<rect class="band" x="{x1:.1f}" y="{yc - 16:.1f}" width="{x2 - x1:.1f}" height="32" rx="3"/>')
        s.append(f'<line class="grid" x1="{ax.x0}" x2="{ax.x0 + ax.w}" y1="{yc:.1f}" y2="{yc:.1f}"/>')
        tcls = 'ink t-b' if m == FOCUS else ''
        s.append(f'<text class="{tcls}" x="{ax.x0 - 10}" y="{yc + 4:.1f}" text-anchor="end">{pretty(m)}</text>')
        pts = sorted(r['vertices'], key=lambda t: t['EF_pin_minus_CBM'])
        last_x = -1e9
        lift = 0
        for t in pts:
            x = ax.X(t['EF_pin_minus_CBM'])
            lift = (lift + 1) % 2 if x - last_x < 22 else 0
            last_x = x
            tip = (f"{pretty(m)} v{t['vertex']}｜E_F* = {t['EF_pin']:+.3f} eV（相对 VBM），E_F* − CBM = {t['EF_pin_minus_CBM']:+.3f} eV｜"
                   f"决定交点：{t['donor']} 与 {t['acceptor']}｜原模型判型 {t['type_model']}，E_F,eq = {t['EF_eq_model']:+.3f} eV")
            s.append(dot(x, yc, 6, cls, tip, hollow=(t['type_model'] == 'p')))
            ly = yc - 13 - 11 * lift
            if t['EF_pin_minus_CBM'] > 0:
                s.append(f'<text class="ink small t-b" x="{ax.X(0) + 8:.1f}" y="{ly:.1f}" text-anchor="start">v{t["vertex"]}：+{t["EF_pin_minus_CBM"]:.3f} eV</text>')
            else:
                s.append(f'<text class="small" x="{x:.1f}" y="{ly:.1f}" text-anchor="middle">v{t["vertex"]}</text>')
    zx = ax.X(0)
    s.append(f'<line class="ref" x1="{zx:.1f}" x2="{zx:.1f}" y1="{ax.y0 - 6}" y2="{ax.y0 + ax.h}"/>')
    s.append(f'<text class="ink" x="{zx + 6:.1f}" y="{ax.y0 - 8}">CBM</text>')
    ly = 24
    s.append(legend_item(130, ly, 'dot', 'f-hl', '实心 = 原模型判 n 型'))
    s.append(legend_item(300, ly, 'hollow', 'k-hl', '空心 = 原模型判 p 型'))
    s.append(legend_item(470, ly, 'band', '', '灰带 = 该材料带隙（VBM→CBM）'))
    s.append('</svg>')
    return ''.join(s), '图1_钉扎费米能级_全部顶点.svg'


# ---------------------------------------------------------------- 图 2：ES–键长
def fig_bonds(R):
    W, H = 760, 400
    ax = Axes(80, 60, 620, 270, (2.96, 3.34), (0.5, 2.15))
    s = [svg_open(W, H, '图2 层内 A–Te 与 B–Te 接触的共享电子数 ES 与键长')]
    s.append(ax.frame(xt=nice_ticks(2.96, 3.34, 8), yt=[0.6, 1.0, 1.4, 1.8], xlabel='键长 d（Å，ES 计算所用几何）', ylabel='共享电子数 ES'))
    for yv, lab in [(1.0, 'ES ≈ 1：metavalent 键（约 2c–1e，Raty 2021）'), (2.0, 'ES ≈ 2：经典二中心二电子共价键')]:
        y = ax.Y(yv)
        s.append(f'<line class="ref" x1="{ax.x0}" x2="{ax.x0 + ax.w}" y1="{y:.1f}" y2="{y:.1f}"/>')
        s.append(f'<text class="small" x="{ax.x0 + ax.w}" y="{y - 6:.1f}" text-anchor="end">{lab}</text>')
    focus_pts = []
    for m, r in R.items():
        A, B = r['A'], r['B']
        items = [(d, e, 'k-s1', f'{A}–Te') for d, e in r['ES_A']]
        items += [(r['ES_B_short'][0], r['ES_B_short'][1], 'k-s2', f'{B}–Te 短'), (r['ES_B_long'][0], r['ES_B_long'][1], 'k-s3', f'{B}–Te 长')]
        for d, e, cls, lab in items:
            tip = f'{pretty(m)}｜{lab}｜d = {d:.2f} Å，ES = {e:.3f}'
            if m == FOCUS:
                focus_pts.append((d, e, cls, lab, tip))
            else:
                s.append(dot(ax.X(d), ax.Y(e), 4.5, cls, tip))
    for d, e, cls, lab, tip in focus_pts:
        x, y = ax.X(d), ax.Y(e)
        s.append(f'<circle class="f-surf k-hl" cx="{x:.1f}" cy="{y:.1f}" r="10" stroke-width="2"/>')
        s.append(dot(x, y, 6, cls, tip))
    # Pb2Bi2Te5 直接标注
    lab_pos = {'Pb–Te': (0, 26), 'Bi–Te 短': (16, 4), 'Bi–Te 长': (-14, 26)}
    seen = set()
    for d, e, cls, lab, tip in focus_pts:
        if lab in seen:
            continue
        seen.add(lab)
        dx, dy = lab_pos.get(lab, (14, 4))
        anchor = 'start' if dx > 0 else ('end' if dx < 0 else 'middle')
        val = f'{lab} {d:.2f} Å / ES {e:.2f}'
        s.append(f'<text class="ink small" x="{ax.X(d) + dx:.1f}" y="{ax.Y(e) + dy:.1f}" text-anchor="{anchor}">{val}</text>')
    ly = 24
    s.append(legend_item(80, ly, 'dot', 'f-s1', 'A–Te（Pb/Sn 位）'))
    s.append(legend_item(240, ly, 'dot', 'f-s2', 'B–Te 短键'))
    s.append(legend_item(360, ly, 'dot', 'f-s3', 'B–Te 长键'))
    s.append(f'<circle class="f-surf k-hl" cx="{486}" cy="{ly - 4}" r="7" stroke-width="2"/><text x="{500}" y="{ly}">紫圈 = {pretty(FOCUS)}（带数值标注）</text>')
    s.append('</svg>')
    return ''.join(s), '图2_ES与键长.svg'


# ---------------------------------------------------------------- 图 3：交换能–位点失配
def fig_exchange(R):
    W, H = 760, 380
    ax = Axes(80, 50, 600, 260, (0.0, 0.125), (0.15, 0.55))
    s = [svg_open(W, H, '图3 反位交换能与 A、B 位平均键长差')]
    s.append(ax.frame(xt=nice_ticks(0, 0.125, 5), yt=nice_ticks(0.15, 0.55, 4), xlabel='|〈d(A–Te)〉 − 〈d(B–Te)〉|（Å）：A、B 两个阳离子位的尺寸差',
                      ylabel='交换能 E(B_A⁰)+E(A_B⁰)（eV）'))
    s.append('<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
             '<path d="M0,0 L10,5 L0,10 z" style="fill:var(--ink2)"/></marker></defs>')
    P = {m: (ax.X(abs(r['mismatch'])), ax.Y(r['exchange'])) for m, r in R.items()}
    for a, b in [('Pb2Sb2Te5', 'Pb2Bi2Te5'), ('Sn2Sb2Te5', 'Sn2Bi2Te5')]:
        (x1, y1), (x2, y2) = P[a], P[b]
        L = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / L, (y2 - y1) / L
        s.append(f'<line class="k-ink2" stroke-width="1.2" x1="{x1 + ux * 12:.1f}" y1="{y1 + uy * 12:.1f}" x2="{x2 - ux * 12:.1f}" y2="{y2 - uy * 12:.1f}" marker-end="url(#arr)"/>')
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        s.append(f'<text class="small" x="{mx + 10:.1f}" y="{my + 16:.1f}">Sb → Bi</text>')
    offs = {'Pb2Bi2Te5': (12, 22, 'start'), 'Pb2Sb2Te5': (-14, -10, 'end'), 'Sn2Bi2Te5': (14, 16, 'start'), 'Sn2Sb2Te5': (-14, -8, 'end')}
    for m, r in R.items():
        x, y = P[m]
        cls = 'k-hl' if m == FOCUS else 'k-mute'
        pair = '—' if r['bound_pair'] is None else f"{r['bound_pair']:.3f} eV"
        tip = (f"{pretty(m)}｜A 位 〈d〉 = {r['dA']:.3f} Å，B 位 〈d〉 = {r['dB']:.3f} Å，差 {r['mismatch']:+.3f} Å｜"
               f"交换能 {r['exchange']:.3f} eV｜束缚反位对 {pair}")
        s.append(dot(x, y, 7, cls, tip))
        dx, dy, anc = offs[m]
        tcls = 'ink t-b' if m == FOCUS else 'ink'
        s.append(f'<text class="{tcls}" x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anc}">{pretty(m)}  {r["exchange"]:.2f} eV</text>')
    s.append(legend_item(80, 24, 'dot', 'f-hl', f'{pretty(FOCUS)}'))
    s.append(legend_item(220, 24, 'dot', 'f-mute', '对照材料'))
    s.append(f'<text x="340" y="24">箭头：同一 A 位下把 Sb 换成 Bi</text>')
    s.append('</svg>')
    return ''.join(s), '图3_反位交换能与位点失配.svg'


# ---------------------------------------------------------------- 图 4：VBM 中 A-s
def fig_vbm_s(R):
    W, H = 760, 300
    order = ['Pb2Bi2Te5', 'Pb2Sb2Te5', 'Sn2Bi2Te5', 'Sn2Sb2Te5']
    ax = Axes(80, 40, 620, 190, (0, 4), (0, 30))
    s = [svg_open(W, H, '图4 价带顶 0.1 meV 能窗内 A 位 s 轨道的投影占比')]
    s.append(ax.frame(yt=[0, 10, 20, 30], ylabel='VBM 中 A-s 占比（%）'))
    bw = 56
    for i, m in enumerate(order):
        r = R[m]
        xc = ax.X(i + 0.5)
        v = r['VBM_As_01']
        y = ax.Y(v)
        cls = 'f-hl' if m == FOCUS else 'f-mute'
        tip = f"{pretty(m)}｜VBM（0.1 meV 能窗）{r['A']}-s = {v:.1f}%｜50 meV 能窗 = {r['VBM_As_50']:.1f}%"
        h = ax.Y(0) - y
        s.append(f'<g data-tip="{esc(tip)}"><rect class="hit" x="{xc - bw / 2 - 8:.1f}" y="{ax.y0}" width="{bw + 16}" height="{ax.h}"/>'
                 f'<path class="{cls}" d="M{xc - bw / 2:.1f},{ax.Y(0):.1f} V{y + 4:.1f} q0,-4 4,-4 H{xc + bw / 2 - 4:.1f} q4,0 4,4 V{ax.Y(0):.1f} Z"/></g>')
        s.append(f'<text class="ink" x="{xc:.1f}" y="{y - 8:.1f}" text-anchor="middle">{v:.1f}%</text>')
        tcls = 'ink t-b' if m == FOCUS else ''
        s.append(f'<text class="{tcls}" x="{xc:.1f}" y="{ax.y0 + ax.h + 20}" text-anchor="middle">{pretty(m)}</text>')
        s.append(f'<text class="small" x="{xc:.1f}" y="{ax.y0 + ax.h + 36}" text-anchor="middle">{r["A"]}-s</text>')
    s.append(f'<line class="axis" x1="{ax.x0}" x2="{ax.x0 + ax.w}" y1="{ax.Y(0):.1f}" y2="{ax.Y(0):.1f}"/>')
    s.append('</svg>')
    return ''.join(s), '图4_VBM中A位s轨道占比.svg'


# ---------------------------------------------------------------- 图 5：E_F = CBM 时三类缺陷
def fig_cbm_defects(R):
    W, H = 760, 360
    order = ['Pb2Bi2Te5', 'Pb2Sb2Te5', 'Sn2Bi2Te5', 'Sn2Sb2Te5']
    vals = []
    for m in order:
        b = R[m]['best_vertex']
        vals += [b['BA_plus_at_CBM'], b['AB_minus_at_CBM'], b['vA_minus_at_CBM']]
    lo, hi = min(vals + [0]) - 0.05, max(vals) + 0.1
    yt = nice_ticks(lo, hi, 5)
    ax = Axes(80, 50, 620, 220, (0, 4), (min(yt[0], lo), max(yt[-1], hi)))
    s = [svg_open(W, H, '图5 最有利 n 型顶点、E_F 置于导带底时三类缺陷的形成能')]
    s.append(ax.frame(yt=yt, ylabel='E_f（E_F = CBM）／eV'))
    bw, gap = 34, 2
    keys = [('BA_plus_at_CBM', 'f-s1', 'B_A⁺ 施主反位'), ('AB_minus_at_CBM', 'f-s2', 'A_B⁻ 受主反位'), ('vA_minus_at_CBM', 'f-s3', 'A 空位（负电态）')]
    y0 = ax.Y(0)
    for i, m in enumerate(order):
        b = R[m]['best_vertex']
        A, B = R[m]['A'], R[m]['B']
        names = {'BA_plus_at_CBM': f'{B}_{A}⁺', 'AB_minus_at_CBM': f'{A}_{B}⁻', 'vA_minus_at_CBM': f"v_{A}（q={b['vA_q_at_CBM']:+d}）"}
        xc = ax.X(i + 0.5)
        x = xc - (3 * bw + 2 * gap) / 2
        for key, cls, _ in keys:
            v = b[key]
            y = ax.Y(v)
            top, hgt = (y, y0 - y) if v >= 0 else (y0, y - y0)
            tip = f"{pretty(m)} v{b['vertex']}｜E_F = CBM｜{names[key]} = {v:+.3f} eV"
            s.append(f'<g data-tip="{esc(tip)}"><rect class="hit" x="{x - 2:.1f}" y="{ax.y0}" width="{bw + 4}" height="{ax.h}"/>'
                     f'<rect class="{cls}" x="{x:.1f}" y="{top:.1f}" width="{bw}" height="{max(hgt, 1):.1f}" rx="3"/></g>')
            x += bw + gap
        lowest_acc = min(b['AB_minus_at_CBM'], b['vA_minus_at_CBM'])
        verdict = '施主最低 → 可越过 CBM' if b['BA_plus_at_CBM'] < lowest_acc else '受主最低 → 钉在 CBM 下'
        tcls = 'ink t-b' if m == FOCUS else ''
        s.append(f'<text class="{tcls}" x="{xc:.1f}" y="{ax.y0 + ax.h + 20}" text-anchor="middle">{pretty(m)}（v{b["vertex"]}）</text>')
        s.append(f'<text class="small" x="{xc:.1f}" y="{ax.y0 + ax.h + 36}" text-anchor="middle">{verdict}</text>')
    s.append(f'<line class="axis" x1="{ax.x0}" x2="{ax.x0 + ax.w}" y1="{y0:.1f}" y2="{y0:.1f}"/>')
    lx = 80
    for _, cls, lab in keys:
        s.append(legend_item(lx, 24, 'sq', cls, lab))
        lx += 170
    s.append('</svg>')
    return ''.join(s), '图5_E_F在CBM时的施主与受主.svg'


# ---------------------------------------------------------------- 图 6：形成能–E_F 曲线
def fig_curves(curves):
    W, H = 780, 450
    s = [svg_open(W, H, '图6 Pb2Bi2Te5 与 Sn2Bi2Te5 在最有利 n 型顶点的形成能随费米能级变化')]
    panels = ['Pb2Bi2Te5', 'Sn2Bi2Te5']
    yd = (-0.4, 1.0)
    pw = 280
    style = {'BA': ('k-s1', '施主反位 B_A'), 'AB': ('k-s2', '受主反位 A_B'), 'vA': ('k-s3', 'A 空位'), 'other': ('k-mute', '其余缺陷最低包络')}
    for pi, m in enumerate(panels):
        c = curves[m]
        A, B = MATS[m]
        ax = Axes(80 + pi * (pw + 100), 100, pw, 270, (-0.15, 0.35), yd)
        cid = f'clip{pi}'
        s.append(f'<defs><clipPath id="{cid}"><rect x="{ax.x0}" y="{ax.y0}" width="{ax.w}" height="{ax.h}"/></clipPath></defs>')
        s.append(f'<rect class="band" x="{ax.x0}" y="{ax.y0}" width="{ax.X(0) - ax.x0:.1f}" height="{ax.h}"/>')
        s.append(f'<rect class="band" x="{ax.X(c["Eg"]):.1f}" y="{ax.y0}" width="{ax.x0 + ax.w - ax.X(c["Eg"]):.1f}" height="{ax.h}"/>')
        s.append(ax.frame(xt=[-0.1, 0.0, 0.1, 0.2, 0.3], yt=[-0.4, 0.0, 0.4, 0.8], xlabel='E_F − VBM（eV）', ylabel='形成能 E_f（eV）' if pi == 0 else ''))
        s.append(f'<text class="ink t-b" x="{ax.x0}" y="{ax.y0 - 14}">{pretty(m)}，顶点 v{c["vertex"]}</text>')
        s.append(f'<text class="small" x="{ax.X(0) + 3:.1f}" y="{ax.y0 + 12}">VBM</text>')
        s.append(f'<text class="small" x="{ax.X(c["Eg"]) + 3:.1f}" y="{ax.y0 + 12}">CBM</text>')
        lbl = {'BA': f'{B}_{A}', 'AB': f'{A}_{B}', 'vA': f'v_{A}', 'other': '其余'}
        ends = []
        s.append(f'<g clip-path="url(#{cid})">')
        for key in ['other', 'vA', 'AB', 'BA']:
            cls, _ = style[key]
            pts = c['cv'][key]
            inside = [(x, y) for x, y, *_ in pts if 0 <= x <= c['Eg'] + 1e-9]
            left = [(x, y) for x, y, *_ in pts if x <= 0]
            right = [(x, y) for x, y, *_ in pts if x >= c['Eg'] - 1e-9]
            for seg, ext in [(left, True), (inside, False), (right, True)]:
                if len(seg) < 2:
                    continue
                d = ' '.join(f'{ax.X(x):.1f},{ax.Y(y):.1f}' for x, y in seg)
                s.append(f'<polyline class="line {cls}{" ext" if ext else ""}" points="{d}"/>')
            # 右端直接标注（只标注仍在纵轴范围内的曲线）
            xe, ye = pts[-1][0], pts[-1][1]
            if yd[0] <= ye <= yd[1]:
                ends.append([ax.Y(ye) + 4, lbl[key]])
        s.append('</g>')
        ends.sort()
        for i in range(1, len(ends)):
            ends[i][0] = max(ends[i][0], ends[i - 1][0] + 13)
        for yy, t in ends:
            s.append(f'<text class="small" x="{ax.x0 + ax.w + 5:.1f}" y="{yy:.1f}">{t}</text>')
        # E_F* 标记
        xp = c['EF_pin']
        bp = min(c['cv']['BA'], key=lambda p: abs(p[0] - xp))
        s.append(f'<line class="ref" x1="{ax.X(xp):.1f}" x2="{ax.X(xp):.1f}" y1="{ax.y0}" y2="{ax.y0 + ax.h}"/>')
        tip = f'{pretty(m)} v{c["vertex"]}｜E_F* = {xp:+.3f} eV（相对 VBM），带隙 {c["Eg"]:.4f} eV｜E_F* − CBM = {xp - c["Eg"]:+.3f} eV'
        s.append(dot(ax.X(xp), ax.Y(bp[1]), 5, 'k-hl', tip))
        s.append(f'<text class="ink small" x="{ax.X(xp) + 6:.1f}" y="{ax.y0 + ax.h - 8}">E_F* = {xp:+.3f}</text>')
    lx = 80
    for key in ['BA', 'AB', 'vA', 'other']:
        cls, lab = style[key]
        s.append(legend_item(lx, 24, 'line', cls, lab))
        lx += 160
    s.append(legend_item(80, 50, 'dash', 'k-ink2', '虚线 = 带隙外的公式外推'))
    s.append(legend_item(300, 50, 'band', '', '灰底 = 带隙以外'))
    s.append(f'<circle class="f-hl k-surf" stroke-width="2" cx="{486}" cy="{46}" r="5"/><text x="{498}" y="50">紫点 = 钉扎费米能级 E_F*</text>')
    s.append('</svg>')
    return ''.join(s), '图6_n型最有利顶点形成能曲线.svg'


# ---------------------------------------------------------------- HTML
PAGE_CSS = """
:root{color-scheme:light;%(L)s--page:#f9f9f9;--card:#ffffff;--border:rgba(11,11,11,.10);--accentbg:#f1eefc;--okbg:#eaf6ec;--okink:#006c10;--warnbg:#fdf3e3;--warnink:#8a5300;--nobg:#fdecec;--noink:#a32020;}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;%(D)s--page:#0b0b0b;--card:#191919;--border:rgba(255,255,255,.10);--accentbg:#231e3a;--okbg:#14261a;--okink:#53c556;--warnbg:#2b2112;--warnink:#f0b04a;--nobg:#2d1616;--noink:#f08080;}}
:root[data-theme="dark"]{color-scheme:dark;%(D)s--page:#0b0b0b;--card:#191919;--border:rgba(255,255,255,.10);--accentbg:#231e3a;--okbg:#14261a;--okink:#53c556;--warnbg:#2b2112;--warnink:#f0b04a;--nobg:#2d1616;--noink:#f08080;}
*{box-sizing:border-box}
body{margin:0;background:var(--page);color:var(--ink);font-family:%(F)s;line-height:1.75;font-size:16px;}
main{max-width:980px;margin:0 auto;padding:32px 16px 64px;}
h1{font-size:28px;line-height:1.35;margin:0 0 8px;}
h2{font-size:21px;margin:44px 0 10px;padding-top:8px;border-top:1px solid var(--border);}
h3{font-size:17px;margin:24px 0 6px;}
p{margin:10px 0;}
.meta{color:var(--ink2);font-size:14px;margin-bottom:20px;}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:18px 20px;margin:16px 0;}
.claim{background:var(--accentbg);border-left:4px solid var(--hl);border-radius:6px;padding:12px 16px;margin:16px 0;}
.claim b{color:var(--ink)}
.summary li{margin:6px 0;}
figure{margin:22px 0;background:var(--card);border:1px solid var(--border);border-radius:10px;padding:14px 14px 8px;}
figure .sv{overflow-x:auto;-webkit-overflow-scrolling:touch;}
figure .sv svg{min-width:640px;display:block;}
figcaption{color:var(--ink2);font-size:14px;line-height:1.65;padding:6px 4px 6px;}
figcaption b{color:var(--ink);}
.tablewrap{overflow-x:auto;margin:12px 0 18px;border:1px solid var(--border);border-radius:10px;background:var(--card);}
table{border-collapse:collapse;width:100%%;font-size:14px;}
caption{caption-side:top;text-align:left;padding:10px 12px 4px;color:var(--ink2);font-size:13.5px;}
th,td{padding:7px 10px;border-bottom:1px solid var(--grid);text-align:left;vertical-align:top;}
th{font-weight:600;color:var(--ink);background:var(--page);white-space:nowrap;}
td.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap;}
tr.focus td{background:var(--accentbg);}
.tag{display:inline-block;border-radius:999px;padding:0 9px;font-size:12.5px;font-weight:600;white-space:nowrap;}
.tag.ok{background:var(--okbg);color:var(--okink);} .tag.part{background:var(--warnbg);color:var(--warnink);} .tag.no{background:var(--nobg);color:var(--noink);}
.chain{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin:14px 0;}
.node{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:12px 14px;font-size:14.5px;}
.node h4{margin:0 0 4px;font-size:14.5px;}
.node.bi{border-top:3px solid var(--s3);} .node.pb{border-top:3px solid var(--s1);} .node.out{grid-column:1/-1;border-top:3px solid var(--hl);}
.node .ev{color:var(--ink2);font-size:13px;margin-top:4px;}
.arrow{grid-column:1/-1;text-align:center;color:var(--ink2);font-size:13px;margin:-4px 0;}
details{margin:10px 0;} summary{cursor:pointer;color:var(--ink2);font-size:14px;}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;background:var(--page);border:1px solid var(--border);border-radius:4px;padding:0 4px;word-break:break-all;}
.quote{background:var(--card);border:1px solid var(--border);border-left:4px solid var(--hl);border-radius:6px;padding:14px 18px;}
#tip{position:fixed;pointer-events:none;z-index:10;max-width:340px;background:var(--card);color:var(--ink);border:1px solid var(--border);border-radius:8px;padding:8px 10px;font-size:13px;line-height:1.5;box-shadow:0 4px 18px rgba(0,0,0,.18);display:none;}
.swipe{display:none;color:var(--ink2);font-size:12px;text-align:center;padding-top:4px;}
@media (max-width:640px){.swipe{display:block}.chain{grid-template-columns:1fr}h1{font-size:23px}main{padding:20px 16px 48px}}
"""

TIP_JS = """
(function(){var tip=document.getElementById('tip');
function show(e){var t=e.target.closest&&e.target.closest('[data-tip]');if(!t){tip.style.display='none';return;}
tip.innerHTML=t.getAttribute('data-tip').split('｜').map(function(s){var d=document.createElement('div');d.textContent=s;return d.outerHTML;}).join('');
tip.style.display='block';var x=e.clientX+14,y=e.clientY+14,r=tip.getBoundingClientRect();
if(x+r.width>innerWidth-8)x=e.clientX-r.width-14;if(y+r.height>innerHeight-8)y=e.clientY-r.height-14;tip.style.left=x+'px';tip.style.top=y+'px';}
document.addEventListener('pointermove',show);document.addEventListener('pointerleave',function(){tip.style.display='none';});
document.addEventListener('scroll',function(){tip.style.display='none';},{passive:true});})();
"""


def sci(v):
    if v != v:
        return '—'
    e = int(math.floor(math.log10(abs(v))))
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return f'{v / 10 ** e:.1f}×10{str(e).translate(sup)}'


def f3(v, signed=True):
    s = f'{v:+.3f}' if signed else f'{v:.3f}'
    return s.replace('-', '−')


def build_html(R, vertex_rows, curves, figs, max_check, agree_total, n_total, share_audit):
    F = R[FOCUS]
    fb = F['best_vertex']
    S = R['Sn2Bi2Te5']
    sb = S['best_vertex']
    P = R['Pb2Sb2Te5']
    order = ['Pb2Bi2Te5', 'Pb2Sb2Te5', 'Sn2Bi2Te5', 'Sn2Sb2Te5']
    # n 型点在其自身 E_F,eq 处的施主/受主简并程度与补偿（源模型输出）
    deg_meV = abs(fb['BA_plus_at_eq'] - fb['AB_minus_at_eq']) * 1000
    ion_d = fb['ionized_model'].get('Bi_Pb', float('nan'))
    ion_a = fb['ionized_model'].get('Pb_Bi', float('nan'))
    net_n = fb['n_cm3_model'] - fb['p_cm3_model']
    vA = {m: R[m]['best_vertex']['vA_minus_at_CBM'] for m in R}
    Q = R['Sn2Sb2Te5']
    orbit_bad = all(not v.get('orbit_sizes_match', True) for v in share_audit['checks'].values())

    css = PAGE_CSS % {'L': TOKENS_LIGHT, 'D': TOKENS_DARK, 'F': FONT}
    h = []
    h.append('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">')
    h.append('<title>Pb₂Bi₂Te₅ n 型机制</title>')
    h.append(f'<style>{css}{SVG_CSS}</style></head><body><main>')
    h.append('<h1>Pb₂Bi₂Te₅ 的 n 型倾向从哪里来：长键、反位缺陷与阳离子空位</h1>')
    h.append('<div class="meta">报告补充章节 · 2026-11-25 · 只使用仓库中已有的缺陷形成能、COGITO 成键、ES/ET 与带边投影数据，<b>未运行新的 DFT</b>。'
             '生成脚本 <code>10_脚本/pbbi_ntype_mechanism.py</code>，源数据见本目录 <code>源数据/</code>。</div>')

    h.append('<div class="claim">待检验的命题：<b>“PbBiTe 是 n 型的原因与 Pb 和 Bi 元素有关：它们更容易形成长键的共价键，可以去形成反位缺陷。”</b></div>')
    h.append('<div class="card summary"><b>结论速览</b><ul>')
    h.append(f'<li><b>与 Pb、Bi 有关——成立，但需要两者的组合。</b>四种 A₂B₂Te₅ 共 {n_total} 个化学势顶点中，只有 {pretty(FOCUS)} 的 v{fb["vertex"]} 出现 n 型；'
             f'把 Bi 换成 Sb 或把 Pb 换成 Sn，n 型都消失（图 1）。</li>')
    h.append(f'<li><b>“长键”——主要是 Bi 位的性质，而且更接近 metavalent 键而不是经典共价键。</b>Pb–Te 两组键几乎等长（{F["A_Te"]["d_s"]:.2f}/{F["A_Te"]["d_l"]:.2f} Å）；'
             f'长键在 Bi–Te（{F["B_Te"]["d_l"]:.2f} Å），其 ES ≈ {F["ES_B_long"][1]:.2f}，低于 1（图 2）。</li>')
    h.append(f'<li><b>“长键使反位容易”——是 Bi 的作用，不是 Pb 的作用。</b>阳离子交换能随 A、B 两个位点的尺寸差单调变化；Sb→Bi 降低交换能，Sn→Pb 反而升高（图 3）。'
             f'最容易互换的是 Sn₂Bi₂Te₅（{S["exchange"]:.2f} eV），但它是 p 型，说明“反位容易”只决定补偿强弱，不决定导电符号。</li>')
    h.append(f'<li><b>命题中缺失的一环：Pb 空位受主难以形成。</b>在最有利 n 型的条件下把 E_F 放到 CBM，两种 Pb 化合物的 v_Pb 仍为 {f3(vA["Pb2Bi2Te5"])} / {f3(vA["Pb2Sb2Te5"])} eV，'
             f'两种 Sn 化合物的 v_Sn 已降到 {f3(vA["Sn2Bi2Te5"])} / {f3(vA["Sn2Sb2Te5"])} eV，正是 v_Sn 把 Sn₂Bi₂Te₅ 钉回 p 型（图 5、6）。'
             f'这与 Pb 6s² 惰性的图像一致（Pb₂Bi₂Te₅ 的 VBM 中 Pb-s 仅 {F["VBM_As_01"]:.1f}%，Sn₂Bi₂Te₅ 中 Sn-s 为 {S["VBM_As_01"]:.1f}%），但 VBM 组成并不能单独解释四种材料的差别（图 4）。</li>')
    h.append(f'<li><b>可靠性：</b>唯一的 n 型点施主与受主几乎简并、强烈补偿，结果对计算误差非常敏感；见第 7 节。</li></ul></div>')

    # 0 证据链
    h.append('<h2>0　命题拆解与证据链</h2>')
    h.append('<p>把原命题拆成三个可检验的子命题，再补上数据显示必不可少的第四环。下图是本节的论证结构，各框下方的灰字是对应证据。</p>')
    h.append('<div class="chain">')
    h.append(f'<div class="node bi"><h4>① Bi：长而弱的 Bi–Te 键</h4>Bi 八面体的短/长键 {F["B_Te"]["d_s"]:.2f}/{F["B_Te"]["d_l"]:.2f} Å，长键 ES ≈ {F["ES_B_long"][1]:.2f}，属半满 p-σ（metavalent）型长键。'
             f'<div class="ev">图 2、表 2</div></div>')
    h.append(f'<div class="node pb"><h4>③ Pb：6s² 惰性</h4>Pb–Te 八面体近乎规则（无长短键分化）；Pb₂Bi₂Te₅ 的 VBM 中 Pb-s 仅 {F["VBM_As_01"]:.1f}%（Sn₂Bi₂Te₅ 中 Sn-s 为 {S["VBM_As_01"]:.1f}%）。'
             f'<div class="ev">图 4、表 4</div></div>')
    h.append(f'<div class="node bi"><h4>② Bi 位与 Pb 位相似 → 阳离子交换代价低</h4>两个位点平均键长差 {abs(F["mismatch"]):.3f} Å，Bader 电荷 +{F["qA"]:.2f}/+{F["qB"]:.2f}；交换能 {F["exchange"]:.2f} eV，低于 Pb₂Sb₂Te₅ 的 {P["exchange"]:.2f} eV。'
             f'<div class="ev">图 3、表 3</div></div>')
    h.append(f'<div class="node pb"><h4>④ Pb 空位受主难以形成</h4>E_F = CBM 时 v_Pb 为 {f3(vA["Pb2Bi2Te5"])} / {f3(vA["Pb2Sb2Te5"])} eV（两种 Pb 化合物），v_Sn 为 {f3(vA["Sn2Bi2Te5"])} / {f3(vA["Sn2Sb2Te5"])} eV；Pb 空位无法把费米能级拉回带隙。'
             f'<div class="ev">图 5、图 6</div></div>')
    h.append('<div class="arrow">②（反位多、强补偿）＋ ④（缺少其他受主） ↓</div>')
    h.append(f'<div class="node out"><h4>结果：贫 Pb/富 Te 条件下，浅施主 Bi_Pb 与受主 Pb_Bi 几乎持平且略占上风，E_F* 越过 CBM {fb["EF_pin_minus_CBM"]*1000:.0f} meV → n 型倾向</h4>'
             f'这是一个由补偿差额决定的 n 型：原模型中 Bi_Pb 与 Pb_Bi 的电离浓度分别为 {sci(ion_d)} 与 {sci(ion_a)} cm⁻³，净电子只剩 {sci(net_n)} cm⁻³。<div class="ev">图 1、图 6、表 1</div></div>')
    h.append('</div>')

    h.append('<div class="tablewrap"><table><caption>表 0　原命题各部分与现有数据的对应</caption><thead><tr><th>子命题</th><th>检验所用数据</th><th>结论</th></tr></thead><tbody>')
    h.append('<tr><td>n 型与 Pb 和 Bi 元素有关</td><td>四种材料 18 个顶点的钉扎费米能级（图 1）</td><td><span class="tag ok">成立</span> 只有 Pb＋Bi 组合出现 n 型</td></tr>')
    h.append('<tr><td>Pb、Bi 更容易形成长键</td><td>COGITO 键长/ICOHP，ES、ET（图 2，表 2）</td><td><span class="tag part">部分成立</span> 长键在 Bi 位；Pb–Te 不分化；键型更接近 metavalent</td></tr>')
    h.append('<tr><td>长键使反位缺陷容易形成</td><td>与化学势无关的阳离子交换能（图 3，表 3）</td><td><span class="tag part">部分成立</span> Bi 降低交换能，Pb 反而升高</td></tr>')
    h.append('<tr><td>反位缺陷导致 n 型</td><td>E_F = CBM 时的施主/受主竞争（图 5、6）</td><td><span class="tag part">需补条件</span> 只有在 Pb 空位受主被抑制、且处于贫 Pb 条件时成立</td></tr>')
    h.append('</tbody></table></div>')

    # 1 判据
    h.append('<h2>1　先确认现象：哪些材料、哪些条件给出 n 型</h2>')
    h.append('<p>原模型（缺陷数据包中的电荷平衡输出）在多数顶点出现负形成能与浓度饱和，直接比较载流子浓度并不稳健。为此本节使用一个只依赖现有形成能、物理意义更直接的判据：'
             '在每个化学势顶点下，求<b>最低能带正电缺陷与最低能带负电缺陷形成能相等时的费米能级 E_F*</b>，即本征缺陷对费米能级的钉扎位置。'
             f'E_F* 高于带隙中点记为 n 倾向，否则记为 p 倾向。用这一判据得到的类型与原模型在 <b>{agree_total}/{n_total}</b> 个顶点上一致。'
             f'作为数据复核，由各数据包 v0 记录重建的形成能与源记录中每个顶点的最低缺陷全部吻合（最大偏差 {max_check:.1g} eV，即在源数据的有效位数内一致）。</p>')
    h.append(f'<figure><div class="sv">{figs["fig1"][0]}</div><div class="swipe">← 左右滑动查看完整图 →</div><figcaption><b>图 1　各材料、各化学势顶点的钉扎费米能级相对导带底的位置。</b>横轴为 E_F* − CBM（eV），0 为导带底；每行灰带是该材料自身的带隙（左端为 VBM）。'
             f'实心点为原模型判 n 型，空心点为 p 型；紫色行为 {pretty(FOCUS)}。悬停可查看决定交点的缺陷对。'
             f'<b>读法：</b>只有 {pretty(FOCUS)} 的 v{fb["vertex"]} 越过导带底（+{fb["EF_pin_minus_CBM"]:.3f} eV）；与它只差一个元素的 Pb₂Sb₂Te₅、Sn₂Bi₂Te₅ 最高也停在 CBM 以下 '
             f'{-P["best_vertex"]["EF_pin_minus_CBM"]:.3f} 和 {-sb["EF_pin_minus_CBM"]:.3f} eV。</figcaption></figure>')
    h.append('<div class="tablewrap"><table><caption>表 1　各材料最有利于 n 型的顶点（E_F* 最高者）。完整 18 个顶点见下方折叠表与 <code>源数据/钉扎费米能级_全部顶点.csv</code>。</caption>'
             '<thead><tr><th>材料</th><th>n 型顶点数（原模型）</th><th>最有利顶点</th><th>E_F*（相对 VBM）</th><th>带隙</th><th>E_F* − CBM</th><th>决定交点的施主 / 受主</th></tr></thead><tbody>')
    for m in order:
        r = R[m]
        b = r['best_vertex']
        cls = ' class="focus"' if m == FOCUS else ''
        h.append(f'<tr{cls}><td>{pretty(m)}</td><td class="num">{r["n_type_model"]}/{r["n_vertices"]}</td><td>v{b["vertex"]}</td><td class="num">{f3(b["EF_pin"])} eV</td>'
                 f'<td class="num">{r["Eg"]:.4f} eV</td><td class="num"><b>{f3(b["EF_pin_minus_CBM"])} eV</b></td><td>{esc(b["donor"])} / {esc(b["acceptor"])}</td></tr>')
    h.append('</tbody></table></div>')
    h.append('<details><summary>展开：全部 18 个顶点</summary><div class="tablewrap"><table><thead><tr><th>材料</th><th>顶点</th><th>原模型最低缺陷</th><th>原模型 E_F,eq</th><th>原模型判型</th>'
             '<th>E_F*</th><th>E_F* − CBM</th><th>本判据</th><th>决定交点的施主 / 受主</th></tr></thead><tbody>')
    for t in vertex_rows:
        cls = ' class="focus"' if t['material'] == FOCUS and t['type_model'] == 'n' else ''
        h.append(f'<tr{cls}><td>{pretty(t["material"])}</td><td>v{t["vertex"]}</td><td>{esc(t["min_defect_model"])}，{f3(t["min_defect_Ef"])} eV</td><td class="num">{f3(t["EF_eq_model"])}</td>'
                 f'<td>{t["type_model"]}</td><td class="num">{f3(t["EF_pin"])}</td><td class="num">{f3(t["EF_pin_minus_CBM"])}</td><td>{t["type_pin"]}</td><td>{esc(t["donor"])} / {esc(t["acceptor"])}</td></tr>')
    h.append('</tbody></table></div></details>')
    h.append(f'<p>表 1 的最后一列说明，在四种材料中，决定费米能级位置的几乎都是<b>阳离子反位对 B_A⁺ / A_B⁻</b>（个别顶点为 Te_B⁺ 或 B_Te⁻）。'
             '因此原命题把注意力放在反位缺陷上是对的；接下来要回答的是：为什么只有在 Pb₂Bi₂Te₅ 中，这对反位能把费米能级推过导带底。</p>')

    # 2 长键
    h.append('<h2>2　“长键的共价键”：长键在 Bi 位，键型接近 metavalent</h2>')
    h.append(f'<figure><div class="sv">{figs["fig2"][0]}</div><div class="swipe">← 左右滑动查看完整图 →</div><figcaption><b>图 2　层内 A–Te 与 B–Te 接触的共享电子数 ES 与键长。</b>颜色区分键类型，紫圈为 {pretty(FOCUS)} 并带数值标注；两条参考线取自仓库收录的 Raty 等（2021）：'
             'metavalent 键相邻原子共享约 1 个电子，经典共价单键约 2 个。<b>读法：</b>所有 A–Te 键和 B–Te 长键的 ES 都在 0.75–0.9 之间，B–Te 短键约 1.3，没有一个接近 2。'
             '这些八面体配位的键都是“电子欠缺”的长键，不是经典的二中心二电子共价键。Bi–Te 长键是四种材料中最长、ES 最小的一类。'
             f'注意：ES 的 k 点对称轨道计数检查{"未通过" if orbit_bad else "已通过"}，审查将其列为待验证方法输出，此处只作定性比较。</figcaption></figure>')
    h.append('<div class="tablewrap"><table><caption>表 2　八面体配位键的几何与成键描述符。d、ICOHP 取自 COGITO 键表；ES 取自 ES 逐接触汇总（几何略有不同）；q 为全密度 Bader 电荷（ET）。</caption>'
             '<thead><tr><th>材料</th><th>A–Te d（Å）</th><th>A–Te ICOHP（eV/键）</th><th>A–Te ES</th><th>B–Te d 短/长（Å）</th><th>B–Te ICOHP 短/长</th><th>B–Te ES 短/长</th><th>q(A) / q(B) / q(Te)</th></tr></thead><tbody>')
    for m in order:
        r = R[m]
        a, b = r['A_Te'], r['B_Te']
        ad = f'{a["d_s"]:.2f}' if a['single_shell'] else f'{a["d_s"]:.2f} / {a["d_l"]:.2f}'
        ai = f'{a["icohp_s"]:.2f}' if a['single_shell'] else f'{a["icohp_s"]:.2f} / {a["icohp_l"]:.2f}'
        esA = ' / '.join(f'{e:.2f}' for _, e in r['ES_A'])
        cls = ' class="focus"' if m == FOCUS else ''
        h.append(f'<tr{cls}><td>{pretty(m)}</td><td class="num">{ad}</td><td class="num">{ai.replace("-", "−")}</td><td class="num">{esA}</td><td class="num">{b["d_s"]:.2f} / {b["d_l"]:.2f}</td>'
                 f'<td class="num">{b["icohp_s"]:.2f} / {b["icohp_l"]:.2f}</td><td class="num">{r["ES_B_short"][1]:.2f} / {r["ES_B_long"][1]:.2f}</td>'
                 f'<td class="num">+{r["qA"]:.2f} / +{r["qB"]:.2f} / {r["qTe"]:.2f}</td></tr>'.replace('<td class="num">-', '<td class="num">−'))
    h.append('</tbody></table></div>')
    h.append(f'<p><b>对原命题的修正：</b>Pb 自身并不形成“更长的键”。Pb–Te 两组键几乎等长（{F["A_Te"]["d_s"]:.2f}/{F["A_Te"]["d_l"]:.2f} Å，ICOHP 比 {F["A_Te"]["icohp_s"]/F["A_Te"]["icohp_l"]:.2f}），'
             'Pb 处在近乎规则的八面体中，这正是 6s² 孤对不表现立体化学活性的结构特征。长短键分化来自 B 位：Bi 的长键（ES ≈ 0.76）比 Sb 的长键更长、共享电子更少。'
             '因此建议论文中把“长键的共价键”写成“Bi 位的 metavalent 长键”。</p>')

    # 3 交换能
    h.append('<h2>3　“从而容易形成反位”：交换难易由 A、B 位点的相似度决定</h2>')
    h.append('<p>反位缺陷的形成能含有化学势项，跨材料直接比较没有意义。但把一对互补反位相加时，化学势项正好抵消：'
             '<b>交换能 E(B_A⁰) + E(A_B⁰)</b> 表示把一个 A 和一个 B 互换位置的本征代价，可以在四种材料之间直接比较。</p>')
    h.append(f'<figure><div class="sv">{figs["fig3"][0]}</div><div class="swipe">← 左右滑动查看完整图 →</div><figcaption><b>图 3　阳离子交换能与 A、B 两个位点平均键长差的关系。</b>横轴为 A–Te 与 B–Te 平均键长之差的绝对值，纵轴为交换能；箭头表示同一 A 位下把 Sb 换成 Bi。'
             f'<b>读法：</b>两个位点越接近，交换越便宜，四个点单调排列。Sb→Bi 两次都显著降低交换能（Pb 系 {P["exchange"]:.2f}→{F["exchange"]:.2f} eV；Sn 系 {R["Sn2Sb2Te5"]["exchange"]:.2f}→{S["exchange"]:.2f} eV），'
             f'因为 Bi 的长键把 B 位八面体撑大，使之接近 A 位；而在 B 相同时，Pb 系的交换能高于 Sn 系。只有 4 个数据点，此关系应视为趋势提示，而非拟合规律。</figcaption></figure>')
    h.append('<div class="tablewrap"><table><caption>表 3　反位交换的能量与位点相似度。交换能、离化反位对能量与化学势无关；跃迁能级相对 VBM。</caption>'
             '<thead><tr><th>材料</th><th>〈d〉A / 〈d〉B（Å）</th><th>位点尺寸差</th><th>|q(A) − q(B)|</th><th>交换能 B_A⁰+A_B⁰</th><th>离化对 B_A⁺+A_B⁻</th><th>束缚反位对</th><th>B_A (0/+1)</th><th>A_B (−1/0)</th><th>带隙</th></tr></thead><tbody>')
    for m in sorted(order, key=lambda k: abs(R[k]['mismatch'])):
        r = R[m]
        cls = ' class="focus"' if m == FOCUS else ''
        bp = '—' if r['bound_pair'] is None else f'{r["bound_pair"]:.3f}'
        h.append(f'<tr{cls}><td>{pretty(m)}</td><td class="num">{r["dA"]:.3f} / {r["dB"]:.3f}</td><td class="num">{abs(r["mismatch"]):.3f} Å</td><td class="num">{abs(r["qA"] - r["qB"]):.3f} e</td>'
                 f'<td class="num"><b>{r["exchange"]:.3f}</b> eV</td><td class="num">{r["ionized_pair"]:.3f} eV</td><td class="num">{bp} eV</td><td class="num">{r["BA_level"]:.3f}</td><td class="num">{r["AB_level"]:.3f}</td><td class="num">{r["Eg"]:.4f}</td></tr>')
    h.append('</tbody></table></div>')
    h.append(f'<p>表 3 给出两点补充。第一，Bi 系中 A、B 的 Bader 电荷几乎相同（{pretty(FOCUS)}：Pb +{F["qA"]:.2f}、Bi +{F["qB"]:.2f}），化学环境的相似也支持“容易互换”；'
             '但电荷差的排序与交换能并不完全一致，几何尺寸是更好的描述量。第二，四种材料的 B_A 施主 (0/+1) 能级都高于各自带隙，都是浅施主，'
             '所以“施主是否浅”不能区分 n 型和 p 型。</p>')
    h.append(f'<p><b>关键推论：</b>交换能最低的是 Sn₂Bi₂Te₅（{S["exchange"]:.3f} eV），而它在所有顶点都是 p 型。可见反位容易只意味着施主、受主反位都多、补偿强，并不决定净载流子的符号。'
             '原命题从“容易反位”直接推到“n 型”，中间缺了一步。</p>')

    # 4 Pb 的角色
    h.append('<h2>4　缺失的一环：Pb 的 6s² 惰性抑制阳离子空位受主</h2>')
    h.append(f'<figure><div class="sv">{figs["fig4"][0]}</div><div class="swipe">← 左右滑动查看完整图 →</div><figcaption><b>图 4　价带顶（0.1 meV 能窗，HSE06+SOC）中 A 位 s 轨道的投影占比。</b>紫色为 {pretty(FOCUS)}，悬停可查看 50 meV 能窗的值。'
             f'<b>读法：</b>Pb 的 6s 受相对论效应收缩、能量很深，几乎不参与价带顶（Pb₂Bi₂Te₅ 中 {F["VBM_As_01"]:.1f}%）；Sn₂Bi₂Te₅ 的价带顶中 Sn-s 达 {S["VBM_As_01"]:.1f}%，'
             f'这与 Sn/Ge 碲化物中“ns–Te 5p 反键态抬高价带顶、阳离子空位易形成”的常见图像一致。但 Sn₂Sb₂Te₅ 的精确 VBM 中 Sn-s 只有 {Q["VBM_As_01"]:.1f}%（50 meV 能窗 {Q["VBM_As_50"]:.1f}%），'
             '而它的 v_Sn 同样很低（图 5）。因此 VBM 的 s 占比只是一个辅助指标，稳定的规律是 A 元素本身：Pb 空位总是难形成，Sn 空位总是容易形成。</figcaption></figure>')
    h.append(f'<figure><div class="sv">{figs["fig5"][0]}</div><div class="swipe">← 左右滑动查看完整图 →</div><figcaption><b>图 5　各材料在最有利 n 型的顶点上、把 E_F 放在导带底时，三类关键缺陷的形成能。</b>'
             '蓝色为施主反位 B_A⁺，橙色为受主反位 A_B⁻，绿色为 A 位空位的最低负电态。要让 E_F 越过 CBM，蓝柱必须低于橙柱和绿柱中较低者。'
             f'<b>读法：</b>只有 {pretty(FOCUS)} 满足这一条件；在 Sn₂Bi₂Te₅ 中，v_Sn 在 CBM 处已降到 {f3(sb["vA_minus_at_CBM"])} eV，远低于施主 Bi_Sn⁺（{f3(sb["BA_plus_at_CBM"])} eV），费米能级被空位拉回带隙。</figcaption></figure>')
    h.append('<div class="tablewrap"><table><caption>表 4　最有利 n 型顶点、E_F = CBM 时的施主与受主（单位 eV）。“最低施主/受主”取全部缺陷与电荷态。</caption>'
             '<thead><tr><th>材料（顶点）</th><th>VBM A-s 占比</th><th>B_A⁺</th><th>A_B⁻</th><th>A 空位</th><th>最低施主</th><th>最低受主</th><th>结果</th></tr></thead><tbody>')
    for m in order:
        r = R[m]
        b = r['best_vertex']
        cls = ' class="focus"' if m == FOCUS else ''
        ok = b['lowest_donor_E_at_CBM'] < b['lowest_acceptor_E_at_CBM']
        h.append(f'<tr{cls}><td>{pretty(m)}（v{b["vertex"]}）</td><td class="num">{r["VBM_As_01"]:.1f}%</td><td class="num">{f3(b["BA_plus_at_CBM"])}</td><td class="num">{f3(b["AB_minus_at_CBM"])}</td>'
                 f'<td class="num">{f3(b["vA_minus_at_CBM"])}（q={b["vA_q_at_CBM"]:+d}）</td><td>{esc(b["lowest_donor_at_CBM"])} {f3(b["lowest_donor_E_at_CBM"])}</td>'
                 f'<td>{esc(b["lowest_acceptor_at_CBM"])} {f3(b["lowest_acceptor_E_at_CBM"])}</td><td>{"<span class=\"tag ok\">E_F 可越过 CBM</span>" if ok else "<span class=\"tag no\">钉在 CBM 以下</span>"}</td></tr>')
    h.append('</tbody></table></div>')

    h.append(f'<p><b>小结：</b>阳离子空位的形成能按 A 元素分成两组。在 E_F = CBM 时，v_Pb 在两种 Pb 化合物中都为正（{f3(vA["Pb2Bi2Te5"])}、{f3(vA["Pb2Sb2Te5"])} eV），'
             f'v_Sn 在两种 Sn 化合物中都为负（{f3(vA["Sn2Bi2Te5"])}、{f3(vA["Sn2Sb2Te5"])} eV）。所以 Pb 的作用是关掉阳离子空位这条受主通道。'
             '但 Pb 只是必要条件：Pb₂Sb₂Te₅ 同样没有空位受主，却因为 Pb_Sb 受主反位低于 Sb_Pb 施主而保持 p 型，还需要 Bi 让施主反位足够有利。这正是“Pb＋Bi 组合”的含义。</p>')

    # 5 曲线
    h.append(f'<h2>5　直接对比：{pretty(FOCUS)} 与 Sn₂Bi₂Te₅ 的形成能图</h2>')
    h.append('<p>两者 B 位都是 Bi，差别只在 A 位的 Pb 与 Sn，是检验“Pb 的作用”最干净的一组对照。下图取两者各自 E_F* 最高的顶点（都是贫 A、相对富 Te 的条件）。</p>')
    h.append(f'<figure><div class="sv">{figs["fig6"][0]}</div><div class="swipe">← 左右滑动查看完整图 →</div><figcaption><b>图 6　形成能随费米能级的变化。</b>每条线是该缺陷所有位点、所有电荷态中的最低值，斜率即电荷态；白色区域为带隙，灰底和虚线为带隙外的公式外推；'
             '灰线是其余全部缺陷的最低包络；紫点和竖线标出 E_F*。'
             f'<b>读法：</b>左图中，Bi_Pb（施主，蓝）与 Pb_Bi（受主，橙）在带隙内几乎重合，v_Pb（绿）一直处在高位，二者交点 E_F* 落在 CBM 以上 {fb["EF_pin_minus_CBM"]*1000:.0f} meV。'
             f'右图中，v_Sn 以 −2 的斜率迅速下降，在 CBM 附近成为最低的受主，E_F* 被钉在 {f3(sb["EF_pin"])} eV，低于带隙中点。</figcaption></figure>')
    h.append(f'<p>这组对照说明：<b>把 Sn 换成 Pb 并没有让反位更容易（交换能 {S["exchange"]:.2f}→{F["exchange"]:.2f} eV，反而升高），但移走了阳离子空位这条受主通道</b>，'
             '使反位对的平衡在贫 Pb 条件下能够偏向施主一侧。</p>')

    # 6 改写
    h.append('<h2>6　建议的论文表述</h2>')
    h.append(f'<div class="quote">Pb₂Bi₂Te₅ 的 n 型倾向源于 Pb 与 Bi 的组合效应。Bi 形成长的 Bi–Te metavalent 键（{F["B_Te"]["d_l"]:.2f} Å，ES ≈ {F["ES_B_long"][1]:.2f}），'
             f'使 Bi 位八面体在尺寸和电荷上接近 Pb 位（平均键长差 {abs(F["mismatch"]):.3f} Å，Bader 电荷 +{F["qA"]:.2f}/+{F["qB"]:.2f}），阳离子交换能降至 {F["exchange"]:.2f} eV，'
             '因此 Bi_Pb 与 Pb_Bi 反位大量共存并相互补偿；在贫 Pb、富 Te 条件下，浅施主 Bi_Pb 与受主 Pb_Bi 能量相当。与此同时，Pb 空位受主难以形成'
             f'（E_F = CBM 时仍为 {f3(fb["vA_minus_at_CBM"])} eV），这与 Pb 6s² 的惰性一致（VBM 中 Pb-s 仅约 {F["VBM_As_01"]:.0f}%），费米能级因而可以接近并越过导带底。'
             '在 Sn₂Bi₂Te₅ 中，Sn 空位受主占据这一位置，材料保持 p 型。</div>')
    h.append('<p>如果需要保留原句的结构，最小修改是：“……与 Pb 和 Bi 元素的组合有关：<b>Bi 的长键</b>使 Bi、Pb 两个位点相近，容易形成反位缺陷；<b>而 Pb 的惰性 6s² 抑制了阳离子空位</b>，使施主反位 Bi_Pb 得以主导。”</p>')

    # 7 局限
    h.append('<h2>7　局限与待验证事项</h2>')
    h.append('<div class="tablewrap"><table><thead><tr><th>问题</th><th>具体情况</th><th>对本节结论的影响</th></tr></thead><tbody>')
    h.append(f'<tr><td>n 型点高度敏感</td><td>{pretty(FOCUS)} v{fb["vertex"]} 在其 E_F,eq（{f3(fb["EF_eq_model"])} eV）处，Bi_Pb⁺ 与 Pb_Bi⁻ 仅差 {deg_meV:.1f} meV；电离浓度 {sci(ion_d)} / {sci(ion_a)} cm⁻³，净电子 {sci(net_n)} cm⁻³</td>'
             '<td>n 型是补偿后的小差额，带电缺陷弛豫或修正误差即可改变符号；应表述为“n 型倾向”，不是稳健预测</td></tr>')
    h.append('<tr><td>负形成能与浓度饱和</td><td>多数顶点最低缺陷的形成能为负，原模型浓度超过稀释极限</td><td>原模型的载流子浓度不可直接引用；本节的 E_F* 判据只比较能量交点，受影响较小，但仍依赖同一套能量</td></tr>')
    h.append('<tr><td>超胞对应关系</td><td>能量来自 81 原子超胞，归档 POSCAR 为 72 原子旧超胞</td><td>无法核对缺陷弛豫结构，也无法做 Bi_Pb 周围的局部成键分析</td></tr>')
    h.append(f'<tr><td>极小带隙</td><td>带隙 {R["Pb2Sb2Te5"]["Eg"]:.3f}–{max(r["Eg"] for r in R.values()):.3f} eV，多个 E_F,eq 与 E_F* 落在带外</td><td>带外部分是公式外推；带隙修正会移动 E_F* 与 CBM 的相对位置</td></tr>')
    h.append(f'<tr><td>ES 方法尚未验证</td><td>共享估计的 k 对称轨道计数检查{"未通过" if orbit_bad else "已通过"}（见 <code>09_审查/共享与局域化审查.json</code>）</td><td>图 2 只用于定性区分 ES≈1 与 ES≈2，不据此给出精确排序</td></tr>')
    h.append('<tr><td>不同哈密顿量</td><td>带边投影为 HSE06+SOC，缺陷能量来自另一套计算</td><td>图 4 与图 5/6 之间只能作定性关联</td></tr>')
    h.append('<tr><td>小样本趋势</td><td>交换能与位点尺寸差的单调关系只有 4 个点</td><td>是提示性趋势，需要更多替代组成（如 Ge 系）检验</td></tr>')
    h.append('<tr><td>判据本身</td><td>E_F* 只比较最低带电缺陷的交点，不含态密度与多缺陷求和</td><td>与原模型判型在 18/18 个顶点一致，但不能替代完整的电荷中性求解</td></tr>')
    h.append('</tbody></table></div>')
    h.append('<h3>把机制坐实所需的计算</h3><ol>'
             '<li>在一致的超胞中充分弛豫 Bi_Pb⁺、Pb_Bi⁻、v_Pb²⁻ 与 Sn 系对应缺陷，并做有限尺寸修正；检验 v3 的施主优势是否超过误差尺度。</li>'
             '<li>用 HSE 或带隙修正后的带边重新求电荷中性解，确认 E_F* 与 CBM 的相对位置。</li>'
             '<li>对 Bi_Pb 与 Bi_Sn 的局部环境做 COHP/差分电荷分析，直接检验“长键使 Bi 适应 A 位”。</li>'
             '<li>比较 v_Pb 与 v_Sn 形成前后的价带顶 s–p 反键权重，闭合“6s² 惰性 → 空位难形成”这一环。</li></ol>')

    # 8 复现
    h.append('<h2>8　复现与数据文件</h2>')
    h.append('<p>在交付目录下运行：</p><p><code>python 10_脚本/pbbi_ntype_mechanism.py .</code></p>'
             '<p>脚本只依赖 Python 标准库，读取 <code>03_缺陷形成能与PN型数据/*.zip</code>、<code>06_COGITO成键分析与COHP/</code>、<code>07_带边能窗/带边能窗统计.csv</code> 与 <code>09_审查/共享与局域化审查.json</code>，生成：</p><ul>'
             '<li><code>源数据/钉扎费米能级_全部顶点.csv</code>：18 个顶点的 E_F*、决定缺陷对，以及 E_F = CBM 时的施主与受主形成能（图 1、5，表 1、4）</li>'
             '<li><code>源数据/成键与反位交换.csv</code>：键长、ICOHP、ES、Bader 电荷、交换能、跃迁能级、带边投影（图 2–4，表 2–3）</li>'
             '<li><code>源数据/n型最有利顶点_形成能曲线.csv</code>：图 6 的逐点数据</li>'
             '<li><code>源数据/汇总.json</code>：全部中间量与复核结果</li>'
             '<li><code>图/</code>：图 1–6 的独立 SVG（浅色），可直接插入论文草稿</li></ul>'
             '<p>形成能按 E_f(q, E_F) = E_f⁰(q) − Σ Δn_i(μ_i − μ_i⁰) + q·E_F 由各数据包的 v0 记录重建；E_F 相对缺陷计算所用 VBM。</p>')

    h.append('<div id="tip" role="tooltip"></div>')
    h.append(f'<script>{TIP_JS}</script></main></body></html>')
    return '\n'.join(h)


if __name__ == '__main__':
    main()
