"""Reproduce concentration comparisons from existing SOC AMSET outputs only."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "05_SOC_AMSET补算"
if not SOURCE.is_dir():
    SOURCE = ROOT / "孤对电子与导电类型_综合分析_20260918" / "05_SOC_AMSET补算"
OUT = Path(__file__).resolve().parent
MATERIALS = ["Pb2Sb2Te5", "Pb2Bi2Te5", "Sn2Sb2Te5", "Sn2Bi2Te5"]
TEMPERATURES = [300, 500, 700]
CONCENTRATIONS = [1e17, 1e18, 1e19, 1e20]
COLORS = ["#24558A", "#CF7832", "#247A60", "#985486"]
LABELS = {
    "Pb2Sb2Te5": r"Pb$_2$Sb$_2$Te$_5$",
    "Pb2Bi2Te5": r"Pb$_2$Bi$_2$Te$_5$",
    "Sn2Sb2Te5": r"Sn$_2$Sb$_2$Te$_5$",
    "Sn2Bi2Te5": r"Sn$_2$Bi$_2$Te$_5$",
}
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "DejaVu Sans"],
    "font.size": 10,
    "axes.unicode_minus": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scalar_mean(tensor):
    return float(np.trace(tensor) / 3)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    OUT.mkdir(parents=True, exist_ok=True)
    rows, hashes = [], {}
    max_offdiag_ratio = 0.0
    for material in MATERIALS:
        source = SOURCE / material / "04_AMSET输运" / "transport.json"
        hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
        data = json.loads(source.read_text(encoding="utf-8"))
        doping = np.asarray(data["doping"], float)
        temperatures = np.asarray(data["temperatures"], float)
        for branch, sign in [("p_input", 1), ("n_input", -1)]:
            for temperature in TEMPERATURES:
                tis = np.flatnonzero(np.isclose(temperatures, temperature, rtol=0, atol=1e-9))
                assert len(tis) == 1, (material, temperature, tis)
                ti = int(tis[0])
                for concentration in CONCENTRATIONS:
                    dis = np.flatnonzero(np.isclose(doping, sign * concentration, rtol=1e-12, atol=0))
                    assert len(dis) == 1, (material, concentration, dis)
                    di = int(dis[0])
                    st = np.asarray(data["seebeck"][di][ti], float)
                    ct = np.asarray(data["conductivity"][di][ti], float)
                    mt = np.asarray(data["mobility"]["overall"][di][ti], float)
                    kt = np.asarray(data["electronic_thermal_conductivity"][di][ti], float)
                    assert all(a.shape == (3, 3) and np.isfinite(a).all() for a in [st, ct, mt, kt])
                    sd, cd, md = np.diag(st), np.diag(ct), np.diag(mt)
                    # S is in microvolt/K and conductivity in S/m. PF is mW/(m K^2).
                    # Directional arithmetic mean, not mean(S)^2 * mean(sigma).
                    pfd = sd**2 * cd * 1e-9
                    for tensor in [st, ct, mt, kt]:
                        offdiag = tensor - np.diag(np.diag(tensor))
                        ratio = np.max(np.abs(offdiag)) / max(np.max(np.abs(np.diag(tensor))), 1e-30)
                        max_offdiag_ratio = max(max_offdiag_ratio, float(ratio))
                    row = {
                        "material": material,
                        "temperature_K": temperature,
                        "input_branch": branch,
                        "net_doping_p_minus_n_cm3": sign * concentration,
                        "abs_net_doping_cm3": concentration,
                        "seebeck_mean_uV_K": scalar_mean(st),
                        "conductivity_mean_S_m": scalar_mean(ct),
                        "mobility_overall_mean_cm2_V_s": scalar_mean(mt),
                        "PF_direction_mean_mW_m_K2": float(np.mean(pfd)),
                        "kappa_e_mean_W_m_K": scalar_mean(kt),
                        "PF_product_of_means_mW_m_K2_NOT_USED": scalar_mean(st)**2 * scalar_mean(ct) * 1e-9,
                        "seebeck_sign_matches_net_input": bool(np.sign(scalar_mean(st)) == sign),
                    }
                    for axis, ss, cc, mm, pp in zip("xyz", sd, cd, md, pfd):
                        row[f"seebeck_{axis}_uV_K"] = float(ss)
                        row[f"conductivity_{axis}_S_m"] = float(cc)
                        row[f"mobility_overall_{axis}_cm2_V_s"] = float(mm)
                        row[f"PF_{axis}_mW_m_K2"] = float(pp)
                    row["source_doping_index"] = di
                    row["source_temperature_index"] = ti
                    row["source_file"] = str(source)
                    rows.append(row)
    assert len(rows) == 96
    write_csv(OUT / "SOC_四材料_300_500_700K_净浓度对照.csv", rows)

    # Independently compare the existing +/-1e19, 300 K summary with raw JSON.
    reference_path = SOURCE / "06_SOC_无SOC对比" / "300K_正负1e19_汇总.csv"
    reference = list(csv.DictReader(reference_path.open(encoding="utf-8-sig", newline="")))
    keys = {
        "seebeck_SOC": "seebeck_mean_uV_K",
        "conductivity_SOC": "conductivity_mean_S_m",
        "mobility_SOC": "mobility_overall_mean_cm2_V_s",
        "PF_SOC": "PF_direction_mean_mW_m_K2",
        "kappa_e_SOC": "kappa_e_mean_W_m_K",
    }
    comparisons = []
    for ref in reference:
        if ref["direction"] != "mean" or float(ref["temperature_K"]) != 300:
            continue
        matches = [r for r in rows if r["material"] == ref["material"] and r["temperature_K"] == 300
                   and np.isclose(r["net_doping_p_minus_n_cm3"], float(ref["doping_cm3"]), rtol=1e-12, atol=0)]
        assert len(matches) == 1
        row = matches[0]
        for reference_key, computed_key in keys.items():
            raw, old = row[computed_key], float(ref[reference_key])
            passed = bool(np.isclose(raw, old, rtol=1e-10, atol=1e-10))
            comparisons.append({
                "material": row["material"], "net_doping_cm3": row["net_doping_p_minus_n_cm3"],
                "parameter": computed_key, "raw_json_recomputed": raw,
                "existing_summary": old, "absolute_difference": abs(raw-old), "passed": passed,
            })
    assert len(comparisons) == 40
    assert all(r["passed"] for r in comparisons)
    write_csv(OUT / "300K_正负1e19_独立核验.csv", comparisons)
    selected = [r for r in rows if r["temperature_K"] == 300 and r["abs_net_doping_cm3"] == 1e19]
    write_csv(OUT / "300K_正负1e19_重新提取.csv", selected)

    for field, ylabel, title, filename in [
        ("seebeck_mean_uV_K", r"$\overline{S}$ (μV/K)", "现有 SOC 输运：Seebeck 与净载流子输入", "图1_Seebeck_净浓度_三温度"),
        ("PF_direction_mean_mW_m_K2", r"$\overline{PF}$ (mW m$^{-1}$ K$^{-2}$)", "现有 SOC 输运：方向功率因子均值与净载流子输入", "图2_功率因子_净浓度_三温度"),
    ]:
        fig, axes = plt.subplots(2, 3, figsize=(13.7, 8.1), sharex=True, sharey=True)
        for ri, branch in enumerate(["p_input", "n_input"]):
            for ci, temperature in enumerate(TEMPERATURES):
                ax = axes[ri, ci]
                for material, color in zip(MATERIALS, COLORS):
                    vals = [r for r in rows if r["material"] == material and r["temperature_K"] == temperature and r["input_branch"] == branch]
                    vals.sort(key=lambda r: r["abs_net_doping_cm3"])
                    ax.plot([r["abs_net_doping_cm3"] for r in vals], [r[field] for r in vals],
                            color=color, marker="o", lw=1.8, ms=5, label=LABELS[material])
                ax.set_xscale("log")
                ax.set_xticks(CONCENTRATIONS)
                ax.xaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
                ax.xaxis.set_minor_formatter(NullFormatter())
                ax.grid(True, alpha=0.18)
                ax.axhline(0, color="#5B6672", lw=0.8, alpha=0.7)
                label = "正净输入 p−n > 0" if branch == "p_input" else "负净输入 p−n < 0"
                ax.set_title(f"{temperature} K · {label}", fontsize=10.5)
                if ci == 0:
                    ax.set_ylabel(ylabel)
                if ri == 1:
                    ax.set_xlabel(r"$|p-n|$ 输入 (cm$^{-3}$)")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.948), ncol=4, frameon=False)
        fig.suptitle(title, fontsize=15, y=0.99)
        fig.text(0.5, 0.017, "每条线连接 4 个已有计算点；非拟合。p/n 分栏仅按净输入符号；输运符号可能不同。\n方向算术均值用于材料对照，不代表已建模多晶有效输运；未使用实验数据。", ha="center", va="bottom", fontsize=9, color="#46525E")
        fig.tight_layout(rect=(0.01, 0.075, 0.995, 0.895), h_pad=2)
        fig.savefig(OUT / f"{filename}.png", dpi=200)
        fig.savefig(OUT / f"{filename}.pdf")
        plt.close(fig)

    mismatches = [r for r in rows if not r["seebeck_sign_matches_net_input"]]
    write_csv(OUT / "净输入与Seebeck符号不同_筛查.csv", mismatches)
    audit = {
        "source_sha256": hashes,
        "rows_extracted": len(rows),
        "reference_comparisons": len(comparisons),
        "all_reference_comparisons_passed": True,
        "max_absolute_reference_difference": max(r["absolute_difference"] for r in comparisons),
        "max_offdiagonal_to_diagonal_ratio_selected_tensors": max_offdiag_ratio,
        "net_input_and_seebeck_sign_mismatch_rows": len(mismatches),
        "no_new_dft_or_transport_calculation": True,
    }
    (OUT / "核验摘要.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    text = f"""现有 SOC 输运的净浓度对照（2026-10-07）

范围与来源
仅后处理 20260918 归档的四材料 05_SOC_AMSET补算/<材料>/04_AMSET输运/transport.json。
选择已有共同点：净输入 p−n = ±10^17、±10^18、±10^19、±10^20 cm^-3，T = 300、500、700 K，共 96 行。
未重跑 AMSET 或 DFT，未插入新的计算浓度，未进行有效质量拟合，也未使用实验数据。图中直线仅连接已有点，不能据此定位精确零点、极值或转变浓度。

量及单位
S 用 μV/K；sigma 用 S/m；mu 用 cm^2/(V s)；PF 用 mW/(m K^2)；电子热导用 W/(m K)。
S、sigma、AMSET mobility.overall、电子热导的方向均值均为 (xx+yy+zz)/3。
mu 原样读取 AMSET 的 overall 迁移率并作方向平均，不以 sigma/[e|p−n|] 重新推算，不称 Hall 迁移率；双极区二者不必一致。
PF 对各方向先算 S_ii^2 sigma_ii，再取算术平均，即 sum_i(S_ii^2 sigma_ii)/3 ×10^-9。没有用 mean(S)^2 mean(sigma)。主 CSV 额外保留这个未采用的乘积供核查。
PF 这里是 x/y/z 对角方向结果的算术均值；被选数据的最大非对角/对角比为 {max_offdiag_ratio:.3g}，在此坐标系内近似对角。它不是对任意旋转后全张量公式或实际多晶混合规则的替代。
方向平均是可复现的描述统计，不自动代表实验多晶体的有效 S、sigma 或 PF。

独立核验
300 K、±10^19 cm^-3 的 8 组结果，与旧 SOC/无 SOC 汇总中的 SOC 列独立核对了 S、sigma、mu、PF、电子热导，共 40 项，全部通过 rtol=1e-10、atol=1e-10 检查。
最大绝对差：{audit['max_absolute_reference_difference']:.6g}（不同单位混合，仅表示数值复核精度，详见逐项 CSV）。

净输入不等于单一多数载流子浓度或 Seebeck 判型
AMSET 输入有符号净浓度 p−n：正值表示空穴数目多于电子，负值相反；|p−n| 并不是 p+n，也不是 p 或 n 单独的数目。窄带隙、高温或低净输入时，两种载流子都可能显著参与。
在两载流子描述下 S=(sigma_n S_n+sigma_p S_p)/(sigma_n+sigma_p)，因此净输入符号不能直接决定 S 的符号；当前 96 个点中有 {len(mismatches)} 个方向平均 S 与净输入符号不同。其清单单独保存。
明确例子：300 K、p−n=−10^17 cm^-3，四种材料 S 分别为 Pb2Sb2Te5 +6.21、Pb2Bi2Te5 +61.12、Sn2Sb2Te5 +224.70、Sn2Bi2Te5 +256.85 μV/K。负净输入仍可出现正 S。
反向例子：Pb2Sb2Te5、700 K、p−n=+10^19 cm^-3，S=−5.53 μV/K。图中始终保留 S 正负号，未改绘为绝对值。
这些现象与双极输运及两支输运权重竞争相容；transport.json 本身没有逐支 n、p、sigma_n、sigma_p，所以本次没有从总量倒推出唯一的双极分解，也没有宣称获得实验 Hall 类型。

对有效质量补算的意义
这些曲线提供现有统一浓度基线，表明不能直接把 s 百分比或 Seebeck 大小一一对应。低净输入和高温结果尤其不能套单抛物带 Pisarenko 公式强行提取质量。
下一步应从同一 SOC 能带补三维曲率、DOS/输运质量和带谷能差，并在模型适用的单载流子区间比较；非抛物/双极区保留完整多带输运。现有曲线不确定哪种轨道杂化造成了迁移率差别。

文件
图1：Seebeck 与净输入，p/n 净输入分栏、300/500/700 K 分列，PNG 与矢量 PDF。
图2：按方向先计算再平均的功率因子，同样分栏分温度，PNG 与矢量 PDF。
主 CSV：96 行，含均值、三个方向值、源文件及原始索引。
独立核验 CSV：40 项旧汇总比对；重新提取 CSV：8 组 300 K 正负10^19结果。
核验摘要 JSON 保存原始文件 SHA256，build_soc_concentration_plots.py 可重新生成全部本目录结果。
"""
    (OUT / "结果解释与口径.txt").write_text(text, encoding="utf-8-sig")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
