"""Cross-material Spearman rank correlations between s-fraction descriptors and band/transport quantities (n = 4).
Usage: python rank_table.py <结果目录>   (reads chain.json and mech.json, writes rank_table.json)"""
import json, sys, numpy as np
D = sys.argv[1]
c = json.load(open(f"{D}/chain.json")); mch = json.load(open(f"{D}/mech.json"))
M = ["Pb2Sb2Te5", "Pb2Bi2Te5", "Sn2Sb2Te5", "Sn2Bi2Te5"]
A = {"Pb2Sb2Te5": "Pb", "Pb2Bi2Te5": "Pb", "Sn2Sb2Te5": "Sn", "Sn2Bi2Te5": "Sn"}
def top05(m): return [r for r in mch[m]["profile"] if r["E_mid"] == -0.25][0][A[m] + "-s"]
desc = {"VBM点阳离子s": [c[m]["s_VBM_0.1"]["cat"] for m in M],
        "顶部0.5eV A-s": [top05(m) for m in M],
        "CBM点阳离子s": [c[m]["s_CBM_0.1"]["cat"] for m in M]}
out = {"Eg_HSE": [c[m]["Eg_target"] for m in M], "N_V": [c[m]["NV_300"] for m in M],
       "VB谷数52meV": [c[m]["VB_valleys_52meV"] for m in M], "D_hole": [c[m]["D_hole"] for m in M],
       "S_BT_p1e19": [c[m]["BT_p1e19"]["S_uVK"] for m in M], "S_BT_p1e20": [c[m]["BT_p1e20"]["S_uVK"] for m in M],
       "S_AM_p1e19": [c[m]["AM_p1e19"]["S"] for m in M], "mu_h_AM": [c[m]["AM_p1e19"]["mu"] for m in M],
       "PF_AM_p": [c[m]["AM_p1e19"]["PF"] for m in M], "L_BT_p1e19": [c[m]["BT_p1e19"]["L_1e8"] for m in M]}
rk = lambda a: np.argsort(np.argsort(a))
sp = lambda a, b: round(float(np.corrcoef(rk(np.array(a)), rk(np.array(b)))[0, 1]), 2)
tab = {d: {o: sp(dv, ov) for o, ov in out.items()} for d, dv in desc.items()}
kB = 8.617333e-5
pred = {m: kB * 1e6 * (2.5 + np.log(c[m]["NV_300"] / 1e19)) for m in M}  # non-degenerate S at p = 1e19, r = 0
json.dump({"materials": M, "desc": desc, "out": out, "spearman": tab, "S_pred_nondeg_p1e19": pred},
          open(f"{D}/rank_table.json", "w"), ensure_ascii=False, indent=1)
for d, row in tab.items(): print(d, row)
