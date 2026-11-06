#!/usr/bin/env python3
"""导出“能带粗细与颜色渐变”独立绘图包的数据（2026-11-06，用户要图 10-1 的脚本与数据，放在一个文件夹里）。

用的是图 10-1 自己的取数函数 transport_figures.bands_fn，所以导出的就是图里画的数：
  - 高对称路径 Γ–M–K–Γ–A–L–H–A（到第二个 A 为止）上价带顶以下 14 条、以上 6 条能带（SOC，Kramers 对各自成行），
    能量相对本泛函的真实价带顶，及每个 k 点、每条带的阳离子 s 成分（%）；
  - Sn 系的附加路径 Γ–Σ′–R₂′–L（不在高对称路径上，单独画在右侧小图），k 距离为原始值（作图时横向放大 1.8 倍）；
  - 300 K 的化学势 µ：样品霍尔浓度与同类型 3×10¹⁹ cm⁻³，取自 BoltzTraP2，与能带同一零点；
  - 每种材料的带隙（含附加路径）。
HSE06+SOC 即图 10-1；PBE+SOC 用同一画法，另附。
输出：21_阶段报告与交接/绘图脚本_能带粗细颜色渐变/data/<泛函>/
用法：python 10_脚本/export_fatband_kit.py
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_figures as pf
import transport_figures as tf

KIT = os.path.join(pf.ROOT, '21_阶段报告与交接', '绘图脚本_能带粗细颜色渐变')
FN = {'HSE': 'HSE06+SOC', 'PBE': 'PBE+SOC'}


def w(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        c = csv.writer(f); c.writerow(header); c.writerows(rows)


def main():
    for fn in ('HSE', 'PBE'):
        d = os.path.join(KIT, 'data', FN[fn]); mu_rows, meta = [], []
        for m in pf.MATS:
            B = tf.bands_fn(m, fn); x, E, tp, fs, nel = B['x'], B['E'], B['tp'], B['fs'], B['nel']
            bs = range(nel - 14, nel + 6)
            w(os.path.join(d, f'{m}_main_path.csv'), ['k_index', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV', 'cation_s_percent'],
              [[i, pf.r(x[i], 5), pf.LAB.get(i, ''), pf.bandname(b, nel), pf.r(E[i, b], 5), pf.r(100 * fs[i, b], 2)]
               for b in bs for i in range(pf.KEND + 1)])
            if tp:
                w(os.path.join(d, f'{m}_off_path_G-S-R2-L.csv'), ['i', 'k_distance_1/A', 'label', 'band', 'E-EVBM_eV', 'cation_s_percent'],
                  [[i, pf.r(tp['x'][i] / pf.TP_XSCALE, 5), tp['labels'].get(i, ''), pf.bandname(b, nel), pf.r(tp['E'][i, b], 5), pf.r(100 * tp['fs'][i, b], 2)]
                   for b in bs for i in range(len(tp['x']))])
            NH = pf.HALL[m]['N']
            for case, N in (('Hall', NH), ('3e19', np.sign(NH) * pf.N_REF)):
                mu_rows.append([m, case, f'{N:.3e}', pf.r(B['mu'](N), 5), 'yes' if (case == 'Hall' and pf.HALL[m]['flagged']) else 'no'])
            meta.append([m, pf.r(1000 * B['gap'], 1), 'yes' if tp else 'no', nel])
        w(os.path.join(d, 'chemical_potential_300K.csv'), ['material', 'case', 'carrier_conc_signed_cm-3', 'mu-EVBM_eV', 'Hall_flagged_by_instrument'], mu_rows)
        w(os.path.join(d, 'band_gap.csv'), ['material', 'Eg_meV (incl. off path)', 'has_off_path', 'VB_top_band_index_1based'], meta)
        xlsx(d, FN[fn])
        print('→', d)


def xlsx(d, fn):
    """同一份数据的 Excel（每个 CSV 一页，数值为数字；xlsxwriter 写出，不合并单元格，读回检查），便于 Origin 作图。"""
    import glob, xlsxwriter, plot_data_index as pdi
    path = os.path.join(d, f'作图数据_{fn}.xlsx')
    wb = xlsxwriter.Workbook(path, {'strings_to_formulas': False, 'strings_to_urls': False, 'nan_inf_to_errors': False}); bold = wb.add_format({'bold': True})
    ix = wb.add_worksheet('说明')
    for r, t in enumerate([f'{fn} 能带（粗细与颜色随阳离子 s 成分渐变）的作图数据，与同目录的 CSV 相同。',
                           '*_main_path：高对称路径 Γ–M–K–Γ–A–L–H–A；横轴 k_distance_1/A，纵轴 E-EVBM_eV，按 band 分线；cation_s_percent 决定颜色与线宽。',
                           '*_off_path：Sn 系附加路径 Γ–Σ′–R₂′–L（不在高对称路径上，图中单独成小图、横向 ×1.8）。',
                           'chemical_potential_300K：300 K 的化学势（霍尔浓度与 3×10¹⁹ cm⁻³），画成水平虚线。band_gap：子图标题的带隙。',
                           '能量零点均为本泛函的真实价带顶。生成：10_脚本/export_fatband_kit.py；作图：plot_fatband_gradient.py。']):
        ix.write_string(r, 0, t)
    ix.set_column(0, 0, 120)
    for f in sorted(glob.glob(os.path.join(d, '*.csv'))):
        name = os.path.basename(f)[:-4].replace('Pb2Sb2Te5', 'PbSb').replace('Pb2Bi2Te5', 'PbBi').replace('Sn2Sb2Te5', 'SnSb').replace('Sn2Bi2Te5', 'SnBi')[:31]
        ws = wb.add_worksheet(name)
        for i, row in enumerate(csv.reader(open(f, encoding='utf-8-sig'))):
            for j, v in enumerate(row):
                x = pdi.number(v)
                if i == 0: ws.write_string(i, j, v, bold)
                elif isinstance(x, float): ws.write_number(i, j, x)
                else: ws.write_string(i, j, x)
        ws.freeze_panes(1, 0)
    wb.close(); pdi.check_book(path)


if __name__ == '__main__':
    main()
