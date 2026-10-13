"""Within-material test of the s-p antibonding mechanism on HSE06+SOC path projections (07_带边能窗/源数据/*/解析投影.npz)."""
import json, sys, numpy as np
R = sys.argv[1]; OUT = sys.argv[2]
MATS = {"Pb2Sb2Te5": ("Pb", "Sb", 68), "Pb2Bi2Te5": ("Pb", "Bi", 88), "Sn2Sb2Te5": ("Sn", "Sb", 68), "Sn2Bi2Te5": ("Sn", "Bi", 88)}
res = {}
for m, (A, B, nel) in MATS.items():
    z = np.load(f"{R}/07_带边能窗/源数据/{m}/解析投影.npz", allow_pickle=False)
    K, E, G, groups = z["K"], z["E"], z["G"], list(z["groups"])
    gi = {g: i for i, g in enumerate(groups)}
    tot = G.sum(axis=2); tot[tot == 0] = np.nan
    frac = lambda name: G[:, :, gi[name]] / tot
    Ev = E[:, nel - 1].max(); Ec = E[:, nel].min()
    Erel = E - Ev
    r = {"VBM_check": float(Ev), "gap_path": float(Ec - Ev)}
    # (1) energy profile of A-s / B-s / Te-s weight across the valence band (0.5 eV bins)
    vb = Erel[:, :nel]; prof = []
    for lo in np.arange(-13.0, 0.0, 0.5):
        msk = (vb >= lo) & (vb < lo + 0.5)
        if msk.sum() < 4: continue
        row = {"E_mid": round(lo + 0.25, 2), "n": int(msk.sum())}
        for nm in (f"{A}-s", f"{B}-s", "Te-s", f"{A}-p", f"{B}-p", "Te-p"):
            row[nm] = float(np.nanmean(frac(nm)[:, :nel][msk]) * 100)
        prof.append(row)
    r["profile"] = prof
    # (2) top valence band (Kramers pair: nel-1, nel-2 averaged) along the path: E(k) vs cation-s weight
    top = E[:, nel - 1] - Ev
    wA = (frac(f"{A}-s")[:, nel - 1] + frac(f"{A}-s")[:, nel - 2]) / 2
    wB = (frac(f"{B}-s")[:, nel - 1] + frac(f"{B}-s")[:, nel - 2]) / 2
    w = wA + wB
    ok = np.isfinite(w)
    x = (w / (1 - 2 * w))[ok]; y = top[ok]
    pr = float(np.corrcoef(w[ok], y)[0, 1]); prA = float(np.corrcoef(wA[ok], y)[0, 1])
    slope, icpt = np.polyfit(x, y, 1)
    # rank correlation (Spearman)
    rk = lambda a: np.argsort(np.argsort(a))
    sp = float(np.corrcoef(rk(w[ok]), rk(y))[0, 1])
    r.update(top_VB_pearson_cat_s=pr, top_VB_pearson_A_s=prA, top_VB_spearman_cat_s=sp,
             fit_Delta_eff_eV=float(slope), n_k=int(ok.sum()),
             top_VB_E_range=float(y.max() - y.min()), w_range=[float(w[ok].min()), float(w[ok].max())])
    # (3) same test restricted to the top 0.3 eV of the band (near-edge states only)
    near = ok & (top > -0.3)
    if near.sum() > 5:
        r["near_pearson_cat_s"] = float(np.corrcoef(w[near], top[near])[0, 1]); r["near_n"] = int(near.sum())
    # s weight at the HSE VBM k
    ik = int(np.argmax(E[:, nel - 1])); r["k_VBM_HSE"] = K[ik].round(4).tolist(); r["w_VBM"] = float(w[ik]); r["wA_VBM"] = float(wA[ik])
    res[m] = r
json.dump(res, open(OUT, "w"), indent=1, ensure_ascii=False)
for m, r in res.items():
    print("==", m, {k: v for k, v in r.items() if k != "profile"})
    for row in r["profile"]:
        if row["E_mid"] > -4.5 or row["E_mid"] < -7.5: print("   ", row)
