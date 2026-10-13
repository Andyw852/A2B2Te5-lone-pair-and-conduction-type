#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""render_transport_report.py —— 由合并原始结果重出 AMSET 输运三件整合产物。

产物（默认写入 <root>/AMSET_results/分类输运/）：
  整体输运报告_含平均.md
  整体输运报告_含平均.html
  全量输运数据_含平均_两位小数.csv

输入：
  <root>/AMSET_results/分类输运/00_合并原始结果/<材料>/transport.json[--data-suffix]
  <root>/AMSET_results/分类输运/缺陷条件索引.csv           （缺陷条件标签与顺序）
  <root>/AMSET_results/分类输运/02_缺陷对应浓度/<缺陷>/<vN>/条件说明.json  （n、p、形成能）

用法：
  python render_transport_report.py                       # 就地覆盖三个产物
  python render_transport_report.py --data-suffix .9x9x1 --out-dir /tmp/x --tag .check
"""
import argparse, csv, html as _html, json, re, sys
from pathlib import Path
import numpy as np

MAT_FALLBACK = ["Sn2Sb2Te5", "Pb2Sb2Te5", "Sn2Bi2Te5", "Pb2Bi2Te5"]
DIRS = ["x", "y", "z", "平均"]
TH = ["T/K", "方向", "σ/(S/m)", "S/(μV/K)", "κe/(W/m·K)", "PF/(mW/m·K²)",
      "μ总/(cm²/V·s)", "μADP", "μIMP", "μPOP"]
MD_HEADER = "|" + "|".join(TH) + "|"
MD_SEP = "|" + "|".join(["---"] * len(TH)) + "|"
CSV_COLS = ["材料", "分类", "缺陷", "条件", "doping_cm-3", "T_K", "方向", "sigma_S_m",
            "S_uV_K", "kappa_e_W_mK", "PF_mW_mK2", "mobility_overall_cm2_Vs",
            "mobility_ADP_cm2_Vs", "mobility_IMP_cm2_Vs", "mobility_POP_cm2_Vs", "fermi_eV"]
CSS = ("body{font:15px/1.6 system-ui,sans-serif;max-width:1500px;margin:auto;padding:28px;"
       "color:#203040}h1,h2{color:#123e60}h2{border-bottom:2px solid #a9c8da;padding-top:24px}"
       "nav{position:sticky;top:0;background:#eef5fa;padding:12px;z-index:2}"
       "nav a{margin-right:24px;color:#12618d}.intro{white-space:pre-wrap;font:inherit}"
       ".scroll{overflow:auto}table{border-collapse:collapse;width:100%;"
       "font-variant-numeric:tabular-nums;font-size:13px}th,td{padding:6px 9px;"
       "border:1px solid #d5e0e8;text-align:right;white-space:nowrap}th{background:#e8f0f6}"
       ".avg{background:#e7f3ef;font-weight:600}section{margin-bottom:28px}details{margin:12px 0}"
       "@media print{nav{position:static}thead{display:table-header-group}.scroll{overflow:visible}"
       "table{font-size:8px}body{padding:0}}")


def fmt(v):
    """报告显示规则：两位小数；|v|>=1e4 或 0<|v|<0.01 用科学计数法（尾数两位，指数不带前导零）。"""
    a = abs(float(v))
    if a >= 1e4 or (0.0 < a < 1e-2):
        mant, exp = ("%.2e" % float(v)).split("e")
        e = int(exp)
        return "%se%s%d" % (mant, "+" if e >= 0 else "-", abs(e))
    return "%.2f" % float(v)


def load_json(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def material_order(idx):
    order = []
    for m in idx:
        if m not in order:
            order.append(m)
    for m in MAT_FALLBACK:
        if m not in order:
            order.append(m)
    return order


def defect_rows(idx_path):
    rows = list(csv.DictReader(idx_path.open(encoding="utf-8-sig")))
    out = {}
    for r in rows:
        out.setdefault(r["材料"], []).append(r)
    return out


def condition_title(conc, kind, meta):
    dtype = "n型" if conc < 0 else "p型"
    if kind == "默认扫描":
        return "默认扫描 · %s · %s cm⁻³" % (dtype, fmt(conc))
    return "缺陷条件 %s / %s · %s · %s cm⁻³" % (meta["defect"], meta["vertex"], dtype, fmt(conc))


def build_prose():
    L = []
    A = L.append
    A("## 1. 计算范围与读表说明"); A("")
    A("- 四材料：Sn₂Sb₂Te₅、Pb₂Sb₂Te₅、Sn₂Bi₂Te₅、Pb₂Bi₂Te₅。四份 AMSET 作业均正常完成，退出码 0:0。")
    A("- 默认扫描每材料8个浓度：−10¹⁷至−10²⁰、+10¹⁷至+10²⁰ cm⁻³，每数量级一点；负为n型，正为p型。")
    A("- 另含18个缺陷主导化学势条件的净浓度p−n，不将这些条件限制在默认扫描范围内。")
    A("- 温度100–900 K，每50 K，共17点。合计50个材料-浓度/条件组、850个温度-浓度点，加入x/y/z/平均后共3400条记录。")
    A("- 每项显示两位小数；绝对值≥10⁴或0<绝对值<0.01时采用科学计数法（尾数两位小数）。显示舍入不改变原始计算数据。温度、浓度和所有输运值均按此规则显示。")
    A("- **2026-09-16 重算**：以下全部数值取自 9×9×2（kz≥2）DFPT 介电重算后重跑的 AMSET 结果；重算前的 9×9×1 版本另存为 transport.json.9x9x1。")
    A("")
    A("## 2. Pb₂Bi₂Te₅ / def-013_Bi_Pb / v3 浓度核查"); A("")
    A("该条件的原始缺陷平衡结果为：电子浓度 n = 5.198625706562293e18 cm⁻³，空穴浓度 p = 2.0716843364895095e15 cm⁻³。因此 AMSET 采用的净掺杂约定为 p−n = −5.196554022225803e18 cm⁻³，显示为 **−5.20e18 cm⁻³**。这不是少一个数量级；原始数据指数是 e18，若为 −5.20e19 才会大十倍。")
    A("")
    A("Pb₂Bi₂Te₅ 的 v0–v2 是 p 型，只有 v3 是 n 型，因为 v3 的平衡费米能级 E_F=+0.18062 eV 高于带隙 0.1141 eV，电子浓度约为空穴浓度的 2509.37 倍。v3 的最低缺陷 Bi_Pb(q=+1) 为施主，形成能 0.089813 eV；第二低缺陷 Pb_Bi(q=−1) 的形成能为 0.090517 eV，仅差 0.000704 eV。两者几乎简并并强烈补偿：电离浓度分别约为 1.987e20 和 1.935e20 cm⁻³，差值约 5.22e18 cm⁻³，与净电子浓度同量级。因此 v3 是**强补偿后轻微偏 n 型**，不能把单个缺陷的约 2e20 cm⁻³ 电离浓度直接当作自由电子浓度，也不能将 AMSET 净掺杂放大十倍。")
    A("")
    A("详细记录：02_缺陷对应浓度/Pb2Bi2Te5/def-013_Bi_Pb/v3/浓度核查说明.md。")
    A("")
    A("## 3. 平均的定义及方法限制"); A("")
    A("对σ、S、κe、各通道迁移率，平均=(xx+yy+zz)/3=Tr(A)/3。PF_i=σ_ii×(S_ii×10⁻⁶)²×1000，单位mW/(m·K²)；PF平均=(PF_x+PF_y+PF_z)/3，不等于平均σ×平均S²。")
    A("")
    A("这是三方向描述性算术平均，不是多晶有效介质输运，也不是完整非对角功率因子张量。S保留正负号，不平均绝对值。ADP/IMP/POP为单机制迁移率，总迁移率不是这三个值的算术和。μADP、μIMP、μPOP单位同μ总。")
    A("")
    A("当前保留无SOC的S3/S4波函数与S7形变势，结合PBE+SOC介电与历史HSE+SOC带隙；使用完整介电张量、有效极性声子频率以及重排后的弹性张量。不是全链SOC。默认IMP不能自动重现原缺陷体系强补偿下的总离化杂质浓度；缺陷条件为固定净浓度扫温度，不是逐温度缺陷平衡或显式缺陷结构输运。负形成能、饱和及稀释失效情形仅作为模型情景。未新增DFT或声称网格收敛。2026-09-16 起，介电张量改用 9×9×2（kz≥2）DFPT 重算值。")
    A("")
    A("## 4. 全量结果"); A("")
    A("下面每一组均列出全部17个温度及x/y/z/平均；不只展示300 K代表点。费米能级另列（eV，按AMSET原始输出参考，不另行对齐VBM），不做方向平均。完整非对角矩阵保留于00_合并原始结果/<材料>/transport.json。")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--data-suffix", default="")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    root = Path(a.root)
    base = root / "AMSET_results" / "分类输运"
    outdir = Path(a.out_dir) if a.out_dir else base
    outdir.mkdir(parents=True, exist_ok=True)

    idx = defect_rows(base / "缺陷条件索引.csv")
    mats = material_order(idx)
    prose = build_prose()
    md = ["# 四材料 AMSET 输运整体报告（全量数值与三方向平均）", "", prose, ""]
    html_parts, nav, rows_out = [], [], []

    for m in mats:
        data = load_json(base / "00_合并原始结果" / m / ("transport.json" + a.data_suffix))
        dop = np.array(data["doping"], float)
        T = np.array(data["temperatures"], float)
        sig = np.array(data["conductivity"], float)
        s = np.array(data["seebeck"], float)
        ke = np.array(data["electronic_thermal_conductivity"], float)
        ferm = np.array(data["fermi_levels"], float)
        mu = {k: np.array(v, float) for k, v in data["mobility"].items()}
        nd, nt = len(dop), len(T)
        n_default = 8

        # 缺陷条件：按索引 CSV 顺序匹配到 doping 尾部
        drows = idx.get(m, [])
        used, meta_for = set(), {}
        for i in range(n_default, nd):
            best, bj = None, None
            for j, r in enumerate(drows):
                if j in used:
                    continue
                d = abs(float(r["实际输入浓度_cm-3"]) - dop[i])
                if best is None or d < best:
                    best, bj = d, j
            if bj is None:
                raise SystemExit("无法为 %s 的第 %d 个掺杂点匹配缺陷条件" % (m, i))
            used.add(bj)
            r = drows[bj]
            cj = load_json(base / "02_缺陷对应浓度" / r["结果目录"].split("/", 1)[1] / "条件说明.json")
            meta_for[i] = {"defect": r["最低形成能缺陷"], "vertex": r["化学势条件"],
                           "n": cj["n_cm3"], "p": cj["p_cm3"], "ef": cj["formation_energy_eV"]}

        md += ["### %s" % m, ""]
        html_parts.append('<h2 id="%s">%s</h2>' % (m, m))
        nav.append('<a href="#%s">%s</a>' % (m, m))

        for i in range(nd):
            kind = "默认扫描" if i < n_default else "缺陷条件"
            meta = meta_for.get(i)
            title = condition_title(dop[i], kind, meta)
            src = None
            if meta is not None:
                src = "来源模型：n=%s、p=%s cm⁻³；最低形成能=%s eV。" % (fmt(meta["n"]), fmt(meta["p"]), fmt(meta["ef"]))
                if meta["ef"] < 0:
                    src += "负形成能，注意稀释平衡模型失效风险。"
            fermvals = "；".join("%.2f K: %s" % (T[t], fmt(ferm[i, t])) for t in range(nt))
            fermline = "费米能级/eV：" + fermvals

            body_rows, md_rows = [], []
            for t in range(nt):
                q = {}
                for k in range(3):
                    sg, sv = sig[i, t, k, k], s[i, t, k, k]
                    q.setdefault("sig", []).append(sg)
                    q.setdefault("s", []).append(sv)
                    q.setdefault("ke", []).append(ke[i, t, k, k])
                    q.setdefault("pf", []).append(sg * (sv * 1e-6) ** 2 * 1000.0)
                    for key in ("overall", "ADP", "IMP", "POP"):
                        q.setdefault(key, []).append(mu[key][i, t, k, k])
                avg = {key: float(np.mean(v)) for key, v in q.items()}
                for di, dn in enumerate(DIRS):
                    vals = {key: (avg[key] if dn == "平均" else q[key][di]) for key in avg}
                    cells = ["%.2f" % T[t], dn] + [fmt(vals[key]) for key in
                             ("sig", "s", "ke", "pf", "overall", "ADP", "IMP", "POP")]
                    md_rows.append("|" + "|".join(cells) + "|")
                    if dn == "平均":
                        body_rows.append('<tr class="avg">' + "".join("<td>%s</td>" % c for c in cells) + "</tr>")
                    else:
                        body_rows.append("<tr>" + "".join("<td>%s</td>" % c for c in cells) + "</tr>")
                    rows_out.append({"材料": m, "分类": kind,
                                     "缺陷": meta["defect"] if meta else "",
                                     "条件": meta["vertex"] if meta else "",
                                     "doping_cm-3": fmt(dop[i]), "T_K": "%.2f" % T[t], "方向": dn,
                                     "sigma_S_m": fmt(vals["sig"]), "S_uV_K": fmt(vals["s"]),
                                     "kappa_e_W_mK": fmt(vals["ke"]), "PF_mW_mK2": fmt(vals["pf"]),
                                     "mobility_overall_cm2_Vs": fmt(vals["overall"]),
                                     "mobility_ADP_cm2_Vs": fmt(vals["ADP"]),
                                     "mobility_IMP_cm2_Vs": fmt(vals["IMP"]),
                                     "mobility_POP_cm2_Vs": fmt(vals["POP"]),
                                     "fermi_eV": fmt(ferm[i, t])})

            md += ["#### " + title, ""]
            hp = ["<section>", "<h3>%s</h3>" % title]
            if src:
                md += [src, ""]
                hp.append("<p>%s</p>" % src)
            md += [fermline, "", MD_HEADER, MD_SEP] + md_rows
            hp += ["<details><summary>各温度费米能级（eV）</summary><p>%s</p></details>" % fermvals,
                   '<div class="scroll"><table><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div></section>'
                   % ("".join("<th>%s</th>" % h for h in TH), "".join(body_rows))]
            html_parts.append("".join(hp))
            md.append("")

    md_text = "\n".join(md[:-1]) + "\n"
    (outdir / ("整体输运报告_含平均%s.md" % a.tag)).write_text(md_text, encoding="utf-8", newline="\n")

    doc = ('<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
           '<meta name="viewport" content="width=device-width, initial-scale=1">'
           '<title>四材料输运整体报告</title><style>%s</style><body>'
           '<header><h1>四材料 AMSET 输运整体报告</h1>'
           '<p>全量数值 · x / y / z / 三方向平均 · 两位小数</p></header>'
           '<section><pre class="intro">\n%s\n</pre></section><nav>%s </nav>%s</body></html>'
           % (CSS, _html.escape(prose), " ".join(nav), "".join(html_parts)))
    (outdir / ("整体输运报告_含平均%s.html" % a.tag)).write_text(doc, encoding="utf-8", newline="\n")

    with (outdir / ("全量输运数据_含平均_两位小数%s.csv" % a.tag)).open(
            "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows_out)

    print("材料顺序: %s" % mats)
    print("条件组: %d | 记录: %d" % (sum(len(idx.get(m, [])) + 8 for m in mats), len(rows_out)))
    print("输出 -> %s" % outdir)


if __name__ == "__main__":
    main()
