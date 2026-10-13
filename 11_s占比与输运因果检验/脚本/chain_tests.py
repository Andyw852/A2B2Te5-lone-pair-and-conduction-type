"""Link-by-link tests of: band-edge s fraction -> band quantities -> transport (pure computation).
Inputs (repo-relative): dense PBE+SOC EIGENVAL, HSE+SOC path projections (07), BoltzTraP2 CRTA, AMSET SOC, deformation potentials."""
import json, csv, sys, numpy as np
R = sys.argv[1]; OUT = sys.argv[2]
MATS = ["Pb2Sb2Te5", "Pb2Bi2Te5", "Sn2Sb2Te5", "Sn2Bi2Te5"]
kB = 8.617333262e-5
def read_eigenval(path):
    with open(path) as f: lines = f.read().split("\n")
    nel, nk, nb = map(int, lines[5].split())
    E = np.empty((nk, nb)); K = np.empty((nk, 3)); i = 7
    for ik in range(nk):
        K[ik] = list(map(float, lines[i].split()[:3])); i += 1
        for ib in range(nb):
            E[ik, ib] = float(lines[i].split()[1]); i += 1
        i += 1
    return nel, K, E
def vol(poscar):
    L = [list(map(float, l.split()[:3])) for l in open(poscar).read().split("\n")[2:5]]
    return abs(np.linalg.det(np.array(L)))
res = {}
for m in MATS:
    r = {}
    base = f"{R}/A2B2Te5_SOC_dense/materials/{m}/ke-dft-cpu/result"
    nel, K, E = read_eigenval(f"{base}/step3_uniform/EIGENVAL")
    V = vol(f"{base}/step3_uniform/POSCAR") * 1e-24  # cm^3
    nk = len(K); iv = nel - 1; ic = nel           # SOC: one electron per band
    vb, cb = E[:, :nel], E[:, nel:]
    Ev, Ec = vb.max(), cb.min()
    r["Eg_PBE"] = Ec - Ev
    r["k_VBM"] = K[np.argmax(E[:, iv])].round(4).tolist(); r["k_CBM"] = K[np.argmin(E[:, ic])].round(4).tolist()
    for T in (300,):
        r[f"NV_{T}"] = np.exp(-(Ev - vb) / (kB * T)).sum() / (nk * V)
        r[f"NC_{T}"] = np.exp(-(cb - Ec) / (kB * T)).sum() / (nk * V)
    for w in (0.026, 0.05, 0.1):
        r[f"VB_states_{int(w*1000)}meV"] = (vb > Ev - w).sum() / (nk * V)
        r[f"CB_states_{int(w*1000)}meV"] = (cb < Ec + w).sum() / (nk * V)
    # valleys: local maxima of top VB pair (Kramers: take bands nel-1 and nel-3) on periodic grid
    kf = np.round(K * np.array([28, 28, (7 if E.shape[0] == 5488 else 6)])).astype(int)
    nz = 7 if E.shape[0] == 5488 else 6
    grid = {}
    for ik, (a, b, c) in enumerate(kf): grid[(a % 28, b % 28, c % nz)] = ik
    def valleys(band, sign, Eedge, win):
        out = []
        for (a, b, c), ik in grid.items():
            e = sign * E[ik, band]
            if sign * Eedge - e > win: continue
            ok = True
            for da in (-1, 0, 1):
                for db in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        if da == db == dc == 0: continue
                        j = grid.get(((a + da) % 28, (b + db) % 28, (c + dc) % nz))
                        if j is not None and sign * E[j, band] > e + 1e-6: ok = False
            if ok: out.append((round(sign * Eedge - e, 4), K[ik].round(3).tolist()))
        return out
    for win in (0.052, 0.103):
        vv = valleys(iv, +1, Ev, win) + valleys(iv - 2, +1, Ev, win)
        cc = valleys(ic, -1, Ec, win) + valleys(ic + 2, -1, Ec, win)
        r[f"VB_valleys_{int(win*1000)}meV"] = len(vv); r[f"CB_valleys_{int(win*1000)}meV"] = len(cc)
        if win == 0.103: r["VB_valley_list"] = sorted(vv)[:12]; r["CB_valley_list"] = sorted(cc)[:12]
    # BoltzTraP2 at fixed carrier concentration (300 K)
    d = json.load(open(f"{base}/step8.1_boltztrap/boltztrap_crta.json")); it = d["temperatures_K"].index(300.0)
    N = np.array(d["carrier_conc_cm-3"][it]); o = np.argsort(N)
    def at(key, n): return float(np.interp(n, N[o], np.array(d[key][it])[o]))
    r["Eg_target"] = d["scissor"]["target_gap_eV"]
    for lab, n in (("p1e19", 1e19), ("p1e20", 1e20), ("n1e19", -1e19), ("n1e20", -1e20)):
        S = at("seebeck_V_per_K", n); st = at("sigma_over_tau_S_per_m_s", n); L = at("lorenz_WOhm_per_K2", n)
        r[f"BT_{lab}"] = dict(S_uVK=S * 1e6, sigma_tau_per_carrier=st / (abs(n) * 1e6), L_1e8=L * 1e8, PF_tau=st * S * S)
    Sarr = np.array(d["seebeck_V_per_K"][it]) * 1e6; r["BT_Smax_p"] = float(Sarr.max()); r["BT_p_at_Smax"] = float(N[np.argmax(Sarr)])
    # AMSET SOC 300 K +-1e19 (mean)
    for row in csv.DictReader(open(f"{R}/05_SOC_AMSET补算/06_SOC_无SOC对比/300K_正负1e19_汇总.csv", encoding="utf-8-sig")):
        if row["material"] == m:
            lab = "p" if float(row["doping_cm3"]) > 0 else "n"
            r[f"AM_{lab}1e19"] = dict(S=float(row["seebeck_SOC"]), mu=float(row["mobility_SOC"]), sigma=float(row["conductivity_SOC"]),
                                       kappa_e=float(row["kappa_e_SOC"]), PF=float(row["PF_SOC"]),
                                       S_noSOC=float(row["seebeck_noSOC"]), mu_noSOC=float(row["mobility_noSOC"]))
    t = json.load(open(f"{R}/05_SOC_AMSET补算/{m}/04_AMSET输运/transport.json"))
    r["AM_mobility_keys"] = list(t["mobility"].keys()) if isinstance(t["mobility"], dict) else type(t["mobility"]).__name__
    b = json.load(open(f"{R}/05_SOC_AMSET补算/{m}/03_SOC形变势/band_edges.json"))
    r["D_hole"] = b["hole"]["E1_iso_eV"]; r["D_elec"] = b["electron"]["E1_iso_eV"]; r["k_hole_DP"] = b["hole"]["hits"][0]["k_frac"]
    res[m] = r
# s fractions (HSE+SOC path windows)
for row in csv.DictReader(open(f"{R}/07_带边能窗/带边能窗统计.csv", encoding="utf-8-sig")):
    m = row["material"]; A = "Pb" if row["Pb-s"] else "Sn"; B = "Sb" if row["Sb-s"] else "Bi"
    f = lambda k: float(row[k]) if row[k] else 0.0
    res[m][f"s_{row['edge']}_{row['window_meV']}"] = dict(A=f(A + "-s"), B=f(B + "-s"), Te=f("Te-s"), cat=f(A + "-s") + f(B + "-s"))
json.dump(res, open(OUT, "w"), indent=1, ensure_ascii=False, default=float)
for m, r in res.items():
    print("==", m)
    for k, v in r.items(): print("  ", k, v)
