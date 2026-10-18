#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit the BoltzTraP2 CRTA JSON and extract S and sigma/tau at fixed net
carrier concentrations by interpolation (never extrapolation).

Important conventions recovered from the generator:
  * carrier_conc_cm-3 is (N + nelect)/V_cell with negative = n-type electron,
    positive = p-type hole;
  * the scissor shifts the conduction bands and re-centres the Fermi level at
    the middle of the target gap, so mu_rel_efermi_eV = 0 is the scissored
    mid-gap, not the original PBE Fermi level;
  * sigma/tau and kappa_e/tau are per-relaxation-time quantities;
  * the scalar tensor value is the arithmetic diagonal average
    (xx+yy+zz)/3 for 3D, not an in-plane average and not a tensor invariant.

Writes csv/bt2_checks.csv, csv/bt2_mu_grid_summary.csv,
csv/bt2_interpolated_S_sigma.csv, csv/bt2_interpolation_sensitivity.csv,
json/bt2_audit_summary.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

import audit_common as ac

TARGETS = [-1e20, -1e19, -1e18, 1e18, 1e19, 1e20]
TARGET_LABEL = {-1e20: "n_1e20", -1e19: "n_1e19", -1e18: "n_1e18",
                1e18: "p_1e18", 1e19: "p_1e19", 1e20: "p_1e20"}


def finite_or_none(x):
    x = float(x)
    return x if np.isfinite(x) else None


def safe_rel(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 1e-300)
    if not np.any(mask):
        return None
    return float(np.max(np.abs(a[mask] - b[mask]) / np.abs(b[mask])))


def roots_of_sign_change(x, y):
    roots = []
    for i in range(len(x) - 1):
        if y[i] == 0:
            roots.append(float(x[i]))
        elif y[i] * y[i + 1] < 0:
            try:
                roots.append(float(brentq(lambda z: float(PchipInterpolator(x, y)(z)),
                                          float(x[i]), float(x[i + 1]))))
            except ValueError:
                pass
    if y[-1] == 0:
        roots.append(float(x[-1]))
    return roots


def main():
    hse = ac.read_json(ac.JSON / "s_projection_summary.json")
    inv = {}
    import csv as _csv
    with (ac.CSV / "inventory.csv").open(encoding="utf-8-sig") as fh:
        for row in _csv.DictReader(fh):
            inv[row["material"]] = row
    dense_gap = {}
    with (ac.CSV / "band_edge_check.csv").open(encoding="utf-8-sig") as fh:
        for row in _csv.DictReader(fh):
            dense_gap[row["material"]] = float(row["dense_gap_eV"])

    checks = []
    grid_rows = []
    interp_rows = []
    sens_rows = []
    summary = {"generated_by": "audit_bt2.py", "targets_cm-3": TARGETS,
               "interpolation": "PCHIP on mu(c) and on S(mu), sigma(mu); linear interpolation as sensitivity check",
               "materials": {}}

    for mat in ac.MATERIALS:
        res = ac.RESULT / mat / "ke-dft-cpu" / "result"
        bt = ac.read_json(res / "step8.1_boltztrap/boltztrap_crta.json")
        pos = ac.read_poscar(res / "step3_uniform/POSCAR")
        T = np.asarray(bt["temperatures_K"], dtype=float)
        mu = np.asarray(bt["mu_rel_efermi_eV"], dtype=float)
        c = np.asarray(bt["carrier_conc_cm-3"], dtype=float)
        S = np.asarray(bt["seebeck_V_per_K"], dtype=float)
        sig = np.asarray(bt["sigma_over_tau_S_per_m_s"], dtype=float)
        kap = np.asarray(bt["kappa_e_over_tau_W_per_m_K_s"], dtype=float)
        lor = np.asarray(bt["lorenz_WOhm_per_K2"], dtype=float)
        pf = np.asarray(bt["power_factor_over_tau"], dtype=float)
        sxx = np.asarray(bt["sigma_over_tau_xx"], dtype=float)
        syy = np.asarray(bt["sigma_over_tau_yy"], dtype=float)
        Sxx = np.asarray(bt["seebeck_xx_V_per_K"], dtype=float)
        Syy = np.asarray(bt["seebeck_yy_V_per_K"], dtype=float)
        kxx = np.asarray(bt["kappa_e_over_tau_xx"], dtype=float)
        kyy = np.asarray(bt["kappa_e_over_tau_yy"], dtype=float)
        if "sigma_over_tau_zz" in bt:
            szz = np.asarray(bt["sigma_over_tau_zz"], dtype=float)
            Szz = np.asarray(bt["seebeck_zz_V_per_K"], dtype=float)
            kzz = np.asarray(bt["kappa_e_over_tau_zz"], dtype=float)
            zz_source = "stored"
        else:
            szz = 3.0 * sig - sxx - syy
            Szz = 3.0 * S - Sxx - Syy
            kzz = 3.0 * kap - kxx - kyy
            zz_source = "reconstructed_from_diagonal_average"

        dmu = float(mu[1] - mu[0])
        mu_uniform = bool(np.allclose(np.diff(mu), dmu, rtol=0, atol=1e-9))
        shapes_ok = all(a.shape == (len(T), len(mu)) for a in (c, S, sig, kap, lor, pf))
        finite_ok = all(np.isfinite(a).all() for a in (c, S, sig, kap, lor, pf))
        monotonic = bool(np.all(np.diff(c, axis=1) < 0))
        sig_positive = bool(np.all(sig > 0))
        # diagonal-average checks
        avg_sig_err = safe_rel(sig, (sxx + syy + szz) / 3.0)
        avg_S_err = safe_rel(S, (Sxx + Syy + Szz) / 3.0)
        avg_kap_err = safe_rel(kap, (kxx + kyy + kzz) / 3.0)
        inplane_sig = (sxx + syy) / 2.0
        inplane_S = (Sxx + Syy) / 2.0
        inplane_kap = (kxx + kyy) / 2.0
        anis_sig = finite_or_none(np.nanmedian(np.where(np.abs(inplane_sig) > 1e-300, szz / inplane_sig, np.nan)))
        anis_S = finite_or_none(np.nanmedian(np.where(np.abs(inplane_S) > 1e-300, Szz / inplane_S, np.nan)))
        anis_kap = finite_or_none(np.nanmedian(np.where(np.abs(inplane_kap) > 1e-300, kzz / inplane_kap, np.nan)))
        xy_sig_diff = safe_rel(sxx, syy)
        xy_S_diff = safe_rel(Sxx, Syy)
        # lorenz and pf identities
        lorenz_calc = np.where(sig > 0, kap / (sig * T[:, None]), np.nan)
        lor_err = safe_rel(lorenz_calc, lor)
        pf_calc = S ** 2 * sig
        pf_err = safe_rel(pf_calc, pf)
        # scissor
        sc = bt.get("scissor", {})
        target = float(bt.get("target_gap_eV"))
        dft = float(sc.get("dft_gap_eV"))
        shift = float(sc.get("shift_eV"))
        shift_identity = float((target - dft) - shift)
        hse_gap = float(hse["materials"][mat]["HSE_gap_eV"])
        # volume / nelect
        vol_expected = pos["volume_A3"] * 1e-24
        vol_err = abs(float(bt["cell_volume_cm3"]) - vol_expected) / vol_expected
        nelect = float(bt["nelect"])
        # sign behaviour at the scan ends
        sign_at_minus2 = np.sign(c[:, 0])
        sign_at_plus2 = np.sign(c[:, -1])
        sign_ok = bool(np.all(sign_at_minus2 > 0) and np.all(sign_at_plus2 < 0))
        # zero crossings of c and S
        for iT, temp in enumerate(T):
            cc = c[iT]
            SS = S[iT]
            sigg = sig[iT]
            kapp = kap[iT]
            mu_c0 = None
            if cc[0] > 0 > cc[-1]:
                try:
                    mu_c0 = float(PchipInterpolator(cc[::-1], mu[::-1])(0.0))
                except ValueError:
                    mu_c0 = None
            mu_S0 = None
            c_S0 = None
            roots = roots_of_sign_change(mu, SS)
            if roots:
                if mu_c0 is not None:
                    mu_S0 = min(roots, key=lambda z: abs(z - mu_c0))
                else:
                    mu_S0 = roots[0]
                c_S0 = float(PchipInterpolator(mu, cc)(mu_S0))
            grid_rows.append({
                "material": mat, "temperature_K": temp,
                "dmu_eV": dmu, "mu_min_eV": float(mu[0]), "mu_max_eV": float(mu[-1]),
                "mu_at_c0_eV": mu_c0,
                "mu_at_S0_eV": mu_S0,
                "carrier_conc_at_S0_cm-3": c_S0,
                "mu_c0_minus_mu_S0_eV": (None if (mu_c0 is None or mu_S0 is None) else mu_c0 - mu_S0),
                "c_at_mu0_cm-3": float(np.interp(0.0, mu, cc)),
                "S_at_mu0_uV_K": float(np.interp(0.0, mu, SS) * 1e6),
                "c_at_mu_minus2_cm-3": float(cc[0]), "c_at_mu_plus2_cm-3": float(cc[-1]),
                "S_at_mu_minus2_uV_K": float(SS[0] * 1e6), "S_at_mu_plus2_uV_K": float(SS[-1] * 1e6),
                "sigma_at_mu_minus2_S_per_m_s": float(sigg[0]), "sigma_at_mu_plus2_S_per_m_s": float(sigg[-1]),
                "min_Sigma_over_tau_S_per_m_s": float(np.min(sigg)),
                "argmin_Sigma_mu_eV": float(mu[int(np.argmin(sigg))]),
            })

            # interpolation
            mu_of_c = PchipInterpolator(cc[::-1], mu[::-1], extrapolate=False)
            S_of_mu = PchipInterpolator(mu, SS, extrapolate=False)
            sig_of_mu = PchipInterpolator(mu, sigg, extrapolate=False)
            kap_of_mu = PchipInterpolator(mu, kapp, extrapolate=False)
            Sxx_of_mu = PchipInterpolator(mu, Sxx[iT], extrapolate=False)
            Syy_of_mu = PchipInterpolator(mu, Syy[iT], extrapolate=False)
            Szz_of_mu = PchipInterpolator(mu, Szz[iT], extrapolate=False)
            sxx_of_mu = PchipInterpolator(mu, sxx[iT], extrapolate=False)
            syy_of_mu = PchipInterpolator(mu, syy[iT], extrapolate=False)
            szz_of_mu = PchipInterpolator(mu, szz[iT], extrapolate=False)
            kxx_of_mu = PchipInterpolator(mu, kxx[iT], extrapolate=False)
            kyy_of_mu = PchipInterpolator(mu, kyy[iT], extrapolate=False)
            kzz_of_mu = PchipInterpolator(mu, kzz[iT], extrapolate=False)
            for target_c in TARGETS:
                within = bool(cc[-1] <= target_c <= cc[0])
                row = {
                    "material": mat, "temperature_K": temp, "target_label": TARGET_LABEL[target_c],
                    "target_carrier_conc_cm-3": target_c,
                    "nominal_type": "n" if target_c < 0 else "p",
                    "within_scan_range": within,
                }
                if not within:
                    row.update({"reliability": "X_out_of_range", "reasons": "target outside computed mu scan; no extrapolation"})
                    interp_rows.append(row)
                    continue
                # bracket
                k = int(np.searchsorted(-cc, -target_c))
                i = k - 1 if k > 0 else 0
                i = min(max(i, 0), len(cc) - 2)
                c_hi, c_lo = float(cc[i]), float(cc[i + 1])
                mu_at_c_hi, mu_at_c_lo = float(mu[i]), float(mu[i + 1])
                mu_p = float(mu_of_c(target_c))
                mu_l = mu_at_c_hi + (target_c - c_hi) / (c_lo - c_hi) * (mu_at_c_lo - mu_at_c_hi)
                S_p = float(S_of_mu(mu_p))
                S_l = float(SS[i] + (target_c - c_hi) / (c_lo - c_hi) * (SS[i + 1] - SS[i]))
                sig_p = float(sig_of_mu(mu_p))
                sig_l = float(sigg[i] + (target_c - c_hi) / (c_lo - c_hi) * (sigg[i + 1] - sigg[i]))
                kap_p = float(kap_of_mu(mu_p))
                kap_l = float(kapp[i] + (target_c - c_hi) / (c_lo - c_hi) * (kapp[i + 1] - kapp[i]))
                nearest_index = i if abs(target_c - c_hi) <= abs(target_c - c_lo) else i + 1
                nearest_mu = float(mu[nearest_index])
                nearest_S = float(SS[nearest_index])
                nearest_sig = float(sigg[nearest_index])
                # reliability diagnostics
                sign_consistent = bool(np.sign(S_p) == np.sign(target_c)) if S_p != 0 else False
                near_c0 = (mu_c0 is not None and abs(mu_p - mu_c0) < 1.5 * dmu)
                near_S0 = (mu_S0 is not None and abs(mu_p - mu_S0) < 1.5 * dmu)
                ratio_to_S0 = (None if (c_S0 is None or c_S0 == 0) else abs(target_c) / abs(c_S0))
                pchip_linear_S_diff = abs(S_p - S_l) * 1e6
                pchip_linear_sig_rel = abs(sig_p - sig_l) / abs(sig_l) if sig_l != 0 else None
                reasons = []
                if not sign_consistent:
                    reasons.append("S_sign_opposite_to_nominal_carrier")
                if near_c0:
                    reasons.append("within_1.5_grid_steps_of_net_zero")
                if near_S0:
                    reasons.append("within_1.5_grid_steps_of_S_zero")
                if ratio_to_S0 is not None and target_c < 0 and ratio_to_S0 < 2:
                    reasons.append("electron_target_below_2x_concentration_at_S_zero")
                if abs(S_p) < 20e-6:
                    reasons.append("abs_S_below_20uV_K")
                if pchip_linear_sig_rel is not None and pchip_linear_sig_rel > 0.02:
                    reasons.append("sigma_interpolation_sensitive_above_2pct")
                if pchip_linear_S_diff > max(5.0, 0.05 * abs(S_p * 1e6)):
                    reasons.append("S_interpolation_sensitive_above_5pct_or_5uV_K")
                # Tiered reliability.  C = sign/bipolar failure; B = marginal
                # (near the S=0 or n=0 crossing, small |S|, or interpolation
                # sensitive); A = single-carrier-like at this grid resolution.
                min_grid_distance = min(
                    abs(mu_p - mu_c0) / dmu if mu_c0 is not None else np.inf,
                    abs(mu_p - mu_S0) / dmu if mu_S0 is not None else np.inf)
                if (not sign_consistent) or near_c0 or near_S0 or (
                        target_c < 0 and ratio_to_S0 is not None and ratio_to_S0 < 2):
                    reliability = "C_bipolar_or_sign_unreliable"
                else:
                    b_reasons = []
                    if target_c < 0 and ratio_to_S0 is not None and ratio_to_S0 < 3:
                        b_reasons.append("electron_target_below_3x_concentration_at_S_zero")
                    if min_grid_distance < 2.5:
                        b_reasons.append("within_2.5_grid_steps_of_c0_or_S0")
                    if abs(S_p) < 20e-6:
                        b_reasons.append("abs_S_below_20uV_K")
                    if pchip_linear_S_diff > max(5.0, 0.05 * abs(S_p * 1e6)):
                        b_reasons.append("S_interpolation_sensitive_above_5pct_or_5uV_K")
                    if pchip_linear_sig_rel is not None and pchip_linear_sig_rel > 0.02:
                        b_reasons.append("sigma_interpolation_sensitive_above_2pct")
                    if b_reasons:
                        reliability = "B_marginal"
                        reasons.extend(b_reasons)
                        reasons = list(dict.fromkeys(reasons))
                    else:
                        reliability = "A_single_carrier_like"
                row.update({
                    "mu_pchip_rel_midgap_eV": mu_p,
                    "mu_linear_rel_midgap_eV": mu_l,
                    "mu_midgap_from_VBM_eV": mu_p + target / 2.0,
                    "bracket_i": i,
                    "bracket_mu_low_eV": min(mu_at_c_hi, mu_at_c_lo),
                    "bracket_mu_high_eV": max(mu_at_c_hi, mu_at_c_lo),
                    "bracket_c_low_cm-3": min(c_hi, c_lo),
                    "bracket_c_high_cm-3": max(c_hi, c_lo),
                    "bracket_mu_width_eV": abs(mu_at_c_lo - mu_at_c_hi),
                    "bracket_c_width_cm-3": abs(c_hi - c_lo),
                    "nearest_node_mu_eV": nearest_mu,
                    "nearest_node_S_uV_K": nearest_S * 1e6,
                    "nearest_node_sigma_S_per_m_s": nearest_sig,
                    "mu_distance_to_nearest_node_eV": abs(mu_p - nearest_mu),
                    "S_uV_K": S_p * 1e6,
                    "S_linear_uV_K": S_l * 1e6,
                    "sigma_over_tau_S_per_m_s": sig_p,
                    "sigma_over_tau_linear_S_per_m_s": sig_l,
                    "kappa_e_over_tau_W_per_m_K_s": kap_p,
                    "kappa_e_over_tau_linear_W_per_m_K_s": kap_l,
                    "S_in_plane_uV_K": 0.5 * (float(Sxx_of_mu(mu_p)) + float(Syy_of_mu(mu_p))) * 1e6,
                    "S_zz_uV_K": float(Szz_of_mu(mu_p)) * 1e6,
                    "S_scalar_uV_K": S_p * 1e6,
                    "sigma_in_plane_S_per_m_s": 0.5 * (float(sxx_of_mu(mu_p)) + float(syy_of_mu(mu_p))),
                    "sigma_zz_S_per_m_s": float(szz_of_mu(mu_p)),
                    "sigma_scalar_S_per_m_s": sig_p,
                    "kappa_in_plane_W_per_m_K_s": 0.5 * (float(kxx_of_mu(mu_p)) + float(kyy_of_mu(mu_p))),
                    "kappa_zz_W_per_m_K_s": float(kzz_of_mu(mu_p)),
                    "kappa_scalar_W_per_m_K_s": kap_p,
                    "S_in_plane_over_scalar": (0.5 * (float(Sxx_of_mu(mu_p)) + float(Syy_of_mu(mu_p))) / S_p) if S_p != 0 else None,
                    "sigma_in_plane_over_scalar": (0.5 * (float(sxx_of_mu(mu_p)) + float(syy_of_mu(mu_p))) / sig_p) if sig_p != 0 else None,
                    "tensor_sign_heterogeneous": bool(
                        len({np.sign(S_p), np.sign(0.5 * (float(Sxx_of_mu(mu_p)) + float(Syy_of_mu(mu_p)))),
                             np.sign(float(Szz_of_mu(mu_p)))}) > 1),
                    "lorenz_WOhm_per_K2": float(kap_p / (sig_p * temp)) if sig_p > 0 else None,
                    "S_sign": "n" if S_p < 0 else "p",
                    "sign_consistent": sign_consistent,
                    "mu_at_c0_eV": mu_c0,
                    "mu_at_S0_eV": mu_S0,
                    "carrier_conc_at_S0_cm-3": c_S0,
                    "ratio_abs_c_to_abs_c_at_S0": ratio_to_S0,
                    "mu_distance_to_c0_grid_steps": (None if mu_c0 is None else abs(mu_p - mu_c0) / dmu),
                    "mu_distance_to_S0_grid_steps": (None if mu_S0 is None else abs(mu_p - mu_S0) / dmu),
                    "pchip_linear_S_diff_uV_K": pchip_linear_S_diff,
                    "pchip_linear_sigma_rel_diff": pchip_linear_sig_rel,
                    "reliability": reliability,
                    "reasons": ";".join(reasons),
                })
                interp_rows.append(row)

        checks.append({
            "material": mat,
            "temperatures_K": json.dumps(T.tolist()),
            "mu_npts": len(mu), "mu_min_eV": float(mu[0]), "mu_max_eV": float(mu[-1]),
            "mu_uniform": mu_uniform, "dmu_eV": dmu,
            "shapes_ok": shapes_ok, "finite_ok": finite_ok,
            "carrier_conc_strictly_decreasing": monotonic,
            "sigma_over_tau_positive": sig_positive,
            "carrier_sign_at_mu_minus2_positive": bool(np.all(sign_at_minus2 > 0)),
            "carrier_sign_at_mu_plus2_negative": bool(np.all(sign_at_plus2 < 0)),
            "carrier_sign_convention_ok": sign_ok,
            "carrier_conc_at_mu0_cm-3_by_T": json.dumps([float(np.interp(0.0, mu, c[iT])) for iT in range(len(T))]),
            "max_abs_carrier_conc_at_mu0_cm-3": float(np.max(np.abs([np.interp(0.0, mu, c[iT]) for iT in range(len(T))]))),
            "carrier_conc_zero_offset_note": "c(mu=0) is not forced to zero by the generator; this is the empirical offset of the concentration axis at the scissored mid-gap.",
            "nelect": nelect, "cell_volume_cm3": float(bt["cell_volume_cm3"]),
            "cell_volume_from_poscar_cm3": vol_expected, "cell_volume_rel_error": vol_err,
            "sigma_average_vs_diag_max_rel_error": avg_sig_err,
            "seebeck_average_vs_diag_max_rel_error": avg_S_err,
            "kappa_average_vs_diag_max_rel_error": avg_kap_err,
            "zz_source": zz_source,
            "sigma_zz_over_inplane_median": anis_sig,
            "seebeck_zz_over_inplane_median": anis_S,
            "kappa_zz_over_inplane_median": anis_kap,
            "sigma_xx_vs_yy_max_rel_diff": xy_sig_diff,
            "seebeck_xx_vs_yy_max_rel_diff": xy_S_diff,
            "lorenz_identity_max_rel_error": lor_err,
            "power_factor_identity_max_rel_error": pf_err,
            "scissor_applied": bool(sc.get("applied")),
            "target_gap_eV": target,
            "dft_gap_eV": dft,
            "scissor_shift_eV": shift,
            "shift_identity_error_eV": shift_identity,
            "HSE_gap_eV": hse_gap,
            "target_minus_HSE_meV": 1000.0 * (target - hse_gap),
            "dense_gap_EIGENVAL_eV": dense_gap[mat],
            "dft_gap_minus_dense_meV": 1000.0 * (dft - dense_gap[mat]),
            "mu_zero_is_recentered_midgap": True,
            "VBM_rel_mu_zero_eV": -target / 2.0,
            "CBM_rel_mu_zero_eV": target / 2.0,
            "note_mu_label": "mu_rel_efermi_eV is relative to the scissor-recentred mid-gap Fermi level, not the original PBE Fermi level.",
            "note_offdiagonal": "Only diagonal tensor components are stored; off-diagonal components are not auditable from this JSON.",
            "note_absolute_sigma": "sigma/tau and kappa/tau are per-tau; no absolute conductivity or thermal conductivity is reported.",
        })
        print(f"[{mat}] T={T.tolist()} dmu={dmu:.5f} monotonic={monotonic} "
              f"sign_ok={sign_ok} target={target:.4f} HSE={hse_gap:.4f}")

    # add interpolation sensitivity summary rows (compact)
    for row in interp_rows:
        sens_rows.append({
            "material": row["material"], "temperature_K": row["temperature_K"],
            "target_label": row["target_label"],
            "S_pchip_uV_K": row.get("S_uV_K"), "S_linear_uV_K": row.get("S_linear_uV_K"),
            "S_diff_uV_K": row.get("pchip_linear_S_diff_uV_K"),
            "sigma_pchip": row.get("sigma_over_tau_S_per_m_s"),
            "sigma_linear": row.get("sigma_over_tau_linear_S_per_m_s"),
            "sigma_rel_diff": row.get("pchip_linear_sigma_rel_diff"),
            "bracket_mu_width_eV": row.get("bracket_mu_width_eV"),
            "bracket_c_width_cm-3": row.get("bracket_c_width_cm-3"),
            "reliability": row.get("reliability"),
            "reasons": row.get("reasons"),
        })

    reliable_rows = [row for row in interp_rows
                     if str(row.get("reliability", "")).startswith("A")]
    marginal_rows = [row for row in interp_rows
                     if str(row.get("reliability", "")).startswith("B")]
    unreliable_rows = [row for row in interp_rows
                       if str(row.get("reliability", "")).startswith(("C", "X"))]
    ac.write_csv(ac.CSV / "bt2_reliable_S_sigma.csv", reliable_rows)
    ac.write_csv(ac.CSV / "bt2_marginal_S_sigma.csv", marginal_rows)
    ac.write_csv(ac.CSV / "bt2_unreliable_S_sigma.csv", unreliable_rows)
    ac.write_csv(ac.CSV / "bt2_checks.csv", checks)
    ac.write_csv(ac.CSV / "bt2_mu_grid_summary.csv", grid_rows)
    ac.write_csv(ac.CSV / "bt2_interpolated_S_sigma.csv", interp_rows)
    ac.write_csv(ac.CSV / "bt2_interpolation_sensitivity.csv", sens_rows)
    summary["materials"] = {r["material"]: r for r in checks}
    summary["interpolation_rows"] = len(interp_rows)
    summary["reliability_counts"] = {}
    for row in interp_rows:
        key = row.get("reliability", "unknown")
        summary["reliability_counts"][key] = summary["reliability_counts"].get(key, 0) + 1
    ac.write_json(ac.JSON / "bt2_audit_summary.json", summary)
    print("wrote", ac.CSV / "bt2_interpolated_S_sigma.csv")
    print("reliability counts:", summary["reliability_counts"])


if __name__ == "__main__":
    main()
