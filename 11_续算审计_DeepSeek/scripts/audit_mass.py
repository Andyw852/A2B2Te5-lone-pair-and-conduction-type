#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent re-check of the four SOC local effective-mass fits.

Reads only the original step3_uniform / step8.2_mass / step8.3_mass_fit files.
Re-implements the quadratic fit (see audit_common.fit_quadratic) instead of
importing the production fit_effective_mass.py.

Key checks:
  * dense-grid global VBM/CBM versus mass_metadata edges (band, k, energy);
  * KPOINTS fractions <-> Cartesian q <-> mass_metadata offsets;
  * local EIGENVAL versus electron_band.csv / hole_band.csv;
  * band tracking (is the selected band still the edge band at every local k?);
  * full / inner / outer / no-linear / fractional-coordinate fits;
  * gradient, stationary point, curvature sign, inner-window sensitivity;
  * whether the sampled extremum is at the centre or at the mesh boundary.

Writes csv/mass_fit_recheck.csv, csv/mass_fit_variants.csv,
csv/band_edge_check.csv, csv/band_tracking.csv, csv/mass_point_diagnostics.csv
and json/mass_audit_summary.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import audit_common as ac


def dense_edges(kpts, e, occ):
    occupied = occ > 0.5
    vbm = float(np.max(e[occupied]))
    cbm = float(np.min(e[~occupied]))
    vmask = occupied & (np.abs(e - vbm) < 1e-6)
    cmask = (~occupied) & (np.abs(e - cbm) < 1e-6)
    vlist = [{"k_index_0based": int(i), "band_index_1based": int(j + 1),
              "k_fractional": kpts[i].tolist(), "energy_eV": float(e[i, j])}
             for i, j in zip(*np.where(vmask))]
    clist = [{"k_index_0based": int(i), "band_index_1based": int(j + 1),
              "k_fractional": kpts[i].tolist(), "energy_eV": float(e[i, j])}
             for i, j in zip(*np.where(cmask))]
    return {"vbm_eV": vbm, "cbm_eV": cbm, "gap_eV": cbm - vbm,
            "vbm_states": vlist, "cbm_states": clist,
            "n_occupied": int(occupied.sum()), "n_empty": int((~occupied).sum())}


def match_edge(meta_edge, dense_list):
    target_k = np.asarray(meta_edge["center_fractional"], dtype=float)
    best = None
    for item in dense_list:
        dk = np.asarray(item["k_fractional"], dtype=float) - target_k
        # compare modulo reciprocal lattice in fractional space (wrap to [-0.5,0.5))
        dk = dk - np.round(dk)
        dist = float(np.linalg.norm(dk))
        if best is None or dist < best[0]:
            best = (dist, item)
    if best is None:
        return None
    dist, item = best
    return {"fractional_distance_wrapped": dist,
            "dense_k_fractional": item["k_fractional"],
            "dense_band_index_1based": item["band_index_1based"],
            "dense_energy_eV": item["energy_eV"],
            "same_band": bool(int(meta_edge["band_index_1based"]) == int(item["band_index_1based"])),
            "energy_difference_eV": float(meta_edge["energy_eV"] - item["energy_eV"])}


def verify_kpoints(meta, kp, lattice, eig_kpoints=None):
    B = ac.reciprocal_rows(lattice)
    C = np.linalg.inv(B)
    kfrac = kp["kpoints"]
    n_e = len(meta["points_by_edge"]["electron"])
    all_meta = meta["points_by_edge"]["electron"] + meta["points_by_edge"]["hole"]
    frac_meta = np.asarray([p["k_fractional_unwrapped"] for p in all_meta], dtype=float)
    q_meta = np.asarray([p["offset_cartesian_Ainv"] for p in all_meta], dtype=float)
    grid_meta = np.asarray([p["grid_index"] for p in all_meta], dtype=int)
    # KPOINTS must match metadata fractions exactly (within print precision)
    frac_err = float(np.max(np.abs(kfrac - frac_meta))) if kfrac is not None else None
    # q from fractional difference (centre = index 62 in each 125 block)
    q_from_frac = np.zeros_like(q_meta)
    for carrier, start in (("electron", 0), ("hole", n_e)):
        center = np.asarray(meta["edges"][carrier]["center_fractional"], dtype=float)
        q_from_frac[start:start + 125] = (kfrac[start:start + 125] - center) @ B
    q_err = float(np.max(np.abs(q_from_frac - q_meta)))
    # grid_index * step must equal q_meta
    q_grid = grid_meta.astype(float) * float(meta["step_cartesian_Ainv"])
    q_grid_err = float(np.max(np.abs(q_grid - q_meta)))
    # The centre point must be in the point list at grid index (0,0,0)
    center_indices = [i for i, g in enumerate(grid_meta) if np.all(g == 0)]
    center_k_err = None
    if center_indices:
        center_k_err = float(np.max(np.abs(kfrac[center_indices] - frac_meta[center_indices])))
    eig_k_err = None
    if eig_kpoints is not None and kfrac is not None:
        eig_k_err = float(np.max(np.abs(np.asarray(eig_kpoints, dtype=float) - kfrac)))
    return {"n_kpoints": int(len(kfrac)) if kfrac is not None else None,
            "n_metadata": int(len(all_meta)),
            "max_abs_eigenval_kpoints_vs_KPOINTS": eig_k_err,
            "max_abs_kpoint_fraction_error": frac_err,
            "max_abs_q_from_fraction_error_Ainv": q_err,
            "max_abs_q_from_grid_index_error_Ainv": q_grid_err,
            "center_kpoint_indices_0based": center_indices,
            "max_abs_center_fraction_error": center_k_err,
            "reciprocal_rows_Ainv": B.tolist(),
            "inverse_reciprocal_rows": C.tolist()}


def compare_to_old(report, old):
    out = {}
    old_m = old.get("principal_masses_me")
    new_m = report.get("principal_masses_me")
    if old_m and new_m and len(old_m) == len(new_m):
        out["mass_abs_diff_me"] = [None if (a is None or b is None) else abs(a - b)
                                   for a, b in zip(old_m, new_m)]
        out["mass_rel_diff"] = [None if (a is None or b is None or b == 0) else abs(a - b) / abs(b)
                                for a, b in zip(old_m, new_m)]
    out["old_masses_me"] = old_m
    out["old_gradient_norm"] = old.get("gradient_norm_eV_A")
    out["old_rmse_eV"] = old.get("rmse_eV")
    out["old_stationary_shift"] = old.get("quadratic_stationary_shift_Ainv")
    out["old_inside_box"] = old.get("quadratic_stationary_point_inside_sample_box")
    out["old_warnings"] = old.get("warnings")
    return out


def fit_variants(q, energy, lattice, carrier):
    B = ac.reciprocal_rows(lattice)
    C = np.linalg.inv(B)
    qmax = np.max(np.abs(q), axis=0)
    inner_mask = np.all(np.abs(q) <= 0.005 + 1e-12, axis=1)
    outer_mask = np.any(np.abs(q) >= qmax - 1e-12, axis=1)
    variants = {}

    def add(name, mask=None, linear=True, fractional=False):
        qq = q if mask is None else q[mask]
        ee = energy if mask is None else energy[mask]
        if len(qq) < (10 if linear else 7):
            variants[name] = {"available": False, "reason": f"only {len(qq)} points"}
            return
        if fractional:
            qf = qq @ C
            fit = ac.fit_quadratic(qf, ee, linear=linear)
            hf = fit["hessian_eV_A2"]
            h = C @ hf @ C.T
            fit["hessian_eV_A2"] = h
            fit["hessian_eigenvalues"] = np.linalg.eigvalsh(h)
            fit["gradient_eV_A"] = C @ fit["gradient_eV_A"] if linear else np.zeros(3)
        else:
            fit = ac.fit_quadratic(qq, ee, linear=linear)
        eigval, eigvec, masses = ac.principal_masses(fit["hessian_eV_A2"], carrier)
        shift = ac.stationary_shift(fit["hessian_eV_A2"], fit["gradient_eV_A"]) if linear else None
        variants[name] = {
            "available": True, "n_points": int(len(qq)),
            "linear_term": bool(linear), "fractional_fit": bool(fractional),
            "principal_masses_me": masses,
            "hessian_eigenvalues_eV_A2": eigval.tolist(),
            "gradient_norm_eV_A": float(np.linalg.norm(fit["gradient_eV_A"])),
            "gradient_eV_A": np.asarray(fit["gradient_eV_A"]).tolist(),
            "rmse_eV": fit["rmse_eV"],
            "max_abs_residual_eV": fit["max_abs_residual_eV"],
            "energy_span_eV": fit["energy_span_eV"],
            "condition_number": fit["condition_number"],
            "stationary_shift_Ainv": None if shift is None else np.asarray(shift).tolist(),
            "stationary_inside_box": None if shift is None else ac.inside_box(q, shift),
            "curvature_sign_ok": bool(np.all(eigval > 1e-8)),
        }

    add("full_linear")
    add("inner_central_3x3x3", inner_mask)
    add("outer_shell", outer_mask)
    add("full_no_linear", linear=False)
    add("full_linear_fractional_coords", fractional=True)
    # 1D axis fits
    axis = ac.axis_curvatures(q, energy)
    variants["axis_1d_curvature"] = {"available": True, "curvature_eV_A2": axis,
                                     "note": "diagonal d2E/dq_i2 through the centre; compare to H_ii"}
    return variants, inner_mask, outer_mask


def main():
    summary = {"materials": {}, "generated_by": "audit_mass.py"}
    recheck_rows = []
    variant_rows = []
    edge_rows = []
    tracking_rows = []
    point_rows = []
    variant_summary_rows = []

    for mat in ac.MATERIALS:
        res = ac.RESULT / mat / "ke-dft-cpu" / "result"
        s3, s81, s82, s83 = (res / "step3_uniform", res / "step8.1_boltztrap",
                             res / "step8.2_mass", res / "step8.3_mass_fit")
        pos = ac.read_poscar(s3 / "POSCAR")
        lattice = pos["lattice_rows_A"]
        B = ac.reciprocal_rows(lattice)
        C = np.linalg.inv(B)

        # ---------- dense grid ----------
        kd, ed, od, hd = ac.read_eigenval(s3 / "EIGENVAL")
        kp3 = ac.read_kpoints(s3 / "KPOINTS")
        de = dense_edges(kd, ed, od)
        meta = ac.read_json(s82 / "mass_metadata.json")
        bt = ac.read_json(s81 / "boltztrap_crta.json")
        old = ac.read_json(s83 / "mass_tensor.json")
        edge_match = {c: match_edge(meta["edges"][c], de["vbm_states" if c == "hole" else "cbm_states"])
                      for c in ("electron", "hole")}
        gap_check = {
            "material": mat,
            "dense_nk": hd["nk"], "dense_nb": hd["nb"], "dense_ispin": hd["ispin"],
            "kpoints_mode": kp3["mode"], "kpoints_mesh": "x".join(map(str, kp3["mesh"])) if kp3["mesh"] else "",
            "kpoints_nk": kp3["nk"], "kpoints_product": int(np.prod(kp3["mesh"])) if kp3["mesh"] else None,
            "dense_vbm_eV": de["vbm_eV"], "dense_cbm_eV": de["cbm_eV"], "dense_gap_eV": de["gap_eV"],
            "metadata_hole_k": json.dumps(meta["edges"]["hole"]["center_fractional"]),
            "metadata_hole_band": meta["edges"]["hole"]["band_index_1based"],
            "metadata_hole_E_eV": meta["edges"]["hole"]["energy_eV"],
            "metadata_electron_k": json.dumps(meta["edges"]["electron"]["center_fractional"]),
            "metadata_electron_band": meta["edges"]["electron"]["band_index_1based"],
            "metadata_electron_E_eV": meta["edges"]["electron"]["energy_eV"],
            "metadata_hole_vs_dense": json.dumps(edge_match["hole"], ensure_ascii=False),
            "metadata_electron_vs_dense": json.dumps(edge_match["electron"], ensure_ascii=False),
            "bt2_dft_gap_eV": bt.get("scissor", {}).get("dft_gap_eV"),
            "gap_dense_minus_bt2_meV": 1000.0 * (de["gap_eV"] - float(bt.get("scissor", {}).get("dft_gap_eV", np.nan))),
            "bt2_target_gap_eV": bt.get("target_gap_eV"),
            "bt2_shift_eV": bt.get("scissor", {}).get("shift_eV"),
            "metadata_hole_band_matches_dense_max": bool(edge_match["hole"] and edge_match["hole"]["same_band"]),
            "metadata_electron_band_matches_dense_min": bool(edge_match["electron"] and edge_match["electron"]["same_band"]),
        }
        edge_rows.append(gap_check)

        # ---------- local mesh ----------
        kl, el, ol, hl = ac.read_eigenval(s82 / "EIGENVAL")
        kpl = ac.read_kpoints(s82 / "KPOINTS")
        vcheck = verify_kpoints(meta, kpl, lattice, kl)
        n_e = len(meta["points_by_edge"]["electron"])
        local_center_energy = {}
        local_extremum = {}
        tracking = {}
        for carrier in ("electron", "hole"):
            start = 0 if carrier == "electron" else n_e
            q = np.asarray([p["offset_cartesian_Ainv"] for p in meta["points_by_edge"][carrier]], dtype=float)
            points = meta["points_by_edge"][carrier]
            band = int(meta["edges"][carrier]["band_index_1based"]) - 1
            energy_eig = el[start:start + len(points), band]
            energy_csv = np.asarray([
                float(row["energy_eV"]) for row in
                __import__("csv").DictReader((s83 / (carrier + "_band.csv")).open(encoding="utf-8-sig"))])
            csv_vs_eig = float(np.max(np.abs(energy_csv - energy_eig)))
            center_idx = int(np.argmin(np.linalg.norm(q, axis=1)))
            center_energy = float(energy_eig[center_idx])
            local_center_energy[carrier] = center_energy

            # band tracking
            selected_is_edge = []
            edge_band_ids = []
            degenerate_gap = []
            near_bands = []
            for j in range(len(points)):
                ik = start + j
                occ_idx = np.where(ol[ik] > 0.5)[0]
                empty_idx = np.where(ol[ik] <= 0.5)[0]
                if carrier == "electron":
                    edge_band = int(empty_idx[0]) if len(empty_idx) else None
                else:
                    edge_band = int(occ_idx[-1]) if len(occ_idx) else None
                edge_band_ids.append(edge_band + 1 if edge_band is not None else None)
                ok = bool(edge_band == band or (edge_band is not None and abs(el[ik, edge_band] - el[ik, band]) < 1e-5))
                selected_is_edge.append(ok)
                # Kramers partner (adjacent band) degeneracy
                partner = band + 1 if band + 1 < el.shape[1] else band - 1
                degenerate_gap.append(float(el[ik, partner] - el[ik, band]))
                # other bands within 50 meV (excluding selected and partner)
                close = []
                for b in range(el.shape[1]):
                    if b in (band, partner):
                        continue
                    if abs(el[ik, b] - el[ik, band]) < 0.05:
                        close.append({"band_index_1based": b + 1,
                                      "energy_difference_meV": 1000.0 * float(el[ik, b] - el[ik, band])})
                near_bands.append(close)

            # full fit and variants
            report = {}
            eigval, eigvec, masses = ac.principal_masses(
                ac.fit_quadratic(q, energy_eig)["hessian_eV_A2"], carrier)
            full = ac.fit_quadratic(q, energy_eig)
            shift = ac.stationary_shift(full["hessian_eV_A2"], full["gradient_eV_A"])
            inside = ac.inside_box(q, shift) if shift is not None else None
            report = {
                "carrier": carrier,
                "band_index_1based": band + 1,
                "center_fractional": meta["edges"][carrier]["center_fractional"],
                "center_energy_eV": center_energy,
                "metadata_edge_energy_eV": meta["edges"][carrier]["energy_eV"],
                "center_minus_metadata_edge_eV": center_energy - float(meta["edges"][carrier]["energy_eV"]),
                "principal_masses_me": masses,
                "hessian_eigenvalues_eV_A2": eigval.tolist(),
                "hessian_eV_A2": full["hessian_eV_A2"].tolist(),
                "gradient_eV_A": full["gradient_eV_A"].tolist(),
                "gradient_norm_eV_A": float(np.linalg.norm(full["gradient_eV_A"])),
                "rmse_eV": full["rmse_eV"],
                "max_abs_residual_eV": full["max_abs_residual_eV"],
                "energy_span_eV": full["energy_span_eV"],
                "relative_rmse_to_span": full["rmse_eV"] / full["energy_span_eV"] if full["energy_span_eV"] > 0 else None,
                "quadratic_stationary_shift_Ainv": None if shift is None else np.asarray(shift).tolist(),
                "quadratic_stationary_inside_box": inside,
                "curvature_sign_ok": bool(np.all(eigval > 1e-8)),
                "n_points": len(q),
                "condition_number": full["condition_number"],
                "csv_vs_eigenval_max_abs_eV": csv_vs_eig,
            }
            # predicted energy at the quadratic stationary point
            if shift is not None:
                e_stat = full["intercept_eV"] + 0.5 * float(shift @ full["hessian_eV_A2"] @ shift)
            else:
                e_stat = None
            report["quadratic_stationary_energy_eV"] = e_stat

            actual = ac.actual_extremum(q, energy_eig, carrier)
            report["sampled_extremum"] = actual
            report["center_is_sampled_extremum"] = bool(actual["index"] == center_idx)
            report["energy_gain_center_to_extremum_meV"] = 1000.0 * abs(actual["energy_eV"] - center_energy)
            report["actual_extremum_at_boundary"] = actual["on_boundary"]
            report["actual_extremum_distance_from_center_Ainv"] = actual["distance_from_center_Ainv"]
            report["quadratic_stationary_vs_actual_distance_Ainv"] = (
                None if shift is None else float(np.linalg.norm(shift - np.asarray(actual["q_Ainv"]))))
            report["quadratic_stationary_energy_minus_actual_meV"] = (
                None if e_stat is None else 1000.0 * (e_stat - actual["energy_eV"]))

            # variants
            variants, inner_mask, outer_mask = fit_variants(q, energy_eig, lattice, carrier)
            # compact full/inner/outer variant summary
            def var_mass(name, idx):
                v = variants.get(name, {})
                if not v.get("available") or not v.get("principal_masses_me"):
                    return None
                return v["principal_masses_me"][idx]

            def rel_pct(a, b):
                if a is None or b is None or b == 0:
                    return None
                return 100.0 * abs(a - b) / abs(b)

            vr = {"material": mat, "carrier": carrier}
            for idx in range(3):
                f = var_mass("full_linear", idx)
                inn = var_mass("inner_central_3x3x3", idx)
                out = var_mass("outer_shell", idx)
                vr[f"mass_{idx+1}_full_me"] = f
                vr[f"mass_{idx+1}_inner_me"] = inn
                vr[f"mass_{idx+1}_outer_me"] = out
                vr[f"mass_{idx+1}_inner_change_pct"] = rel_pct(inn, f)
                vr[f"mass_{idx+1}_outer_change_pct"] = rel_pct(out, f)
            inner_changes = [vr.get(f"mass_{i}_inner_change_pct") for i in range(1, 4)]
            outer_changes = [vr.get(f"mass_{i}_outer_change_pct") for i in range(1, 4)]
            inner_changes = [x for x in inner_changes if x is not None]
            outer_changes = [x for x in outer_changes if x is not None]
            vr["max_inner_change_pct"] = max(inner_changes) if inner_changes else None
            vr["max_outer_change_pct"] = max(outer_changes) if outer_changes else None
            vr["full_mass_range_me"] = json.dumps([masses[0], masses[1], masses[2]])
            variant_summary_rows.append(vr)

            # compare old JSON
            report["old_vs_new"] = compare_to_old(report, old["reports"][carrier])
            # inner window change
            inner = variants.get("inner_central_3x3x3", {})
            if inner.get("available"):
                old_eig = np.asarray(ac.principal_masses(full["hessian_eV_A2"], carrier)[0])
                new_eig = np.asarray(inner["hessian_eigenvalues_eV_A2"])
                denom = np.maximum(np.abs(old_eig), 1e-8)
                report["inner_window_relative_eigenvalue_change"] = (np.abs(new_eig - old_eig) / denom).tolist()
            else:
                report["inner_window_relative_eigenvalue_change"] = None

            # reliability label
            reasons = []
            if not report["curvature_sign_ok"]:
                reasons.append("wrong_curvature_sign_at_center")
            if report["actual_extremum_at_boundary"]:
                reasons.append("sampled_extremum_on_mesh_boundary")
            if not report["center_is_sampled_extremum"]:
                reasons.append("center_not_sampled_extremum")
            if inside is False:
                reasons.append("quadratic_stationary_point_outside_mesh")
            if report["gradient_norm_eV_A"] > 1e-3:
                reasons.append("gradient_above_1e-3_eV_A")
            if report["relative_rmse_to_span"] is not None and report["relative_rmse_to_span"] > 0.01:
                reasons.append("relative_rmse_above_1pct")
            if report.get("inner_window_relative_eigenvalue_change") and max(report["inner_window_relative_eigenvalue_change"]) > 0.10:
                reasons.append("inner_window_change_above_10pct")
            if not reasons:
                reliability = "A_reliable_local_curvature"
            elif reasons == ["gradient_above_1e-3_eV_A"]:
                reliability = "B_gradient_warning_only"
            elif "wrong_curvature_sign_at_center" in reasons:
                reliability = "D_invalid_curvature"
            elif "sampled_extremum_on_mesh_boundary" in reasons:
                reliability = "C_not_at_band_edge_mesh_boundary"
            else:
                reliability = "C_needs_recenter"
            report["reliability"] = reliability
            report["reliability_reasons"] = reasons
            report["variants"] = variants

            # point diagnostics
            for j, p in enumerate(points):
                point_rows.append({
                    "material": mat, "carrier": carrier, "point_index_0based": j,
                    "grid_index": json.dumps(p["grid_index"]),
                    "q_x_Ainv": q[j, 0], "q_y_Ainv": q[j, 1], "q_z_Ainv": q[j, 2],
                    "k_fractional": json.dumps(p["k_fractional_unwrapped"]),
                    "energy_eV": float(energy_eig[j]),
                    "csv_energy_eV": float(energy_csv[j]),
                    "is_center": bool(j == center_idx),
                    "is_sampled_extremum": bool(j == actual["index"]),
                    "on_mesh_boundary": bool(np.any(np.abs(q[j]) >= np.max(np.abs(q), axis=0) - 1e-12)),
                    "selected_band_is_edge_band": bool(selected_is_edge[j]),
                    "highest_occupied_band": int(np.where(ol[start + j] > 0.5)[0][-1] + 1),
                    "lowest_empty_band": int(np.where(ol[start + j] <= 0.5)[0][0] + 1),
                    "degenerate_partner_gap_eV": float(degenerate_gap[j]),
                    "n_other_bands_within_50meV": len(near_bands[j]),
                    "fit_residual_eV": float(full["residual_eV"][j]),
                })

            tracking[carrier] = {
                "n_points": len(points),
                "n_selected_is_edge_band": int(sum(selected_is_edge)),
                "selected_band_is_edge_band_all_points": bool(all(selected_is_edge)),
                "edge_band_ids_1based": sorted(set(x for x in edge_band_ids if x is not None)),
                "degenerate_partner_gap_min_eV": float(min(degenerate_gap)),
                "degenerate_partner_gap_max_eV": float(max(degenerate_gap)),
                "n_points_with_other_band_within_50meV": int(sum(1 for x in near_bands if x)),
                "near_band_examples": [x for x in near_bands if x][:5],
            }
            tracking_rows.append({"material": mat, "carrier": carrier,
                                  **{k: v for k, v in tracking[carrier].items()
                                     if k not in ("near_band_examples",)},
                                  "near_band_examples": json.dumps(tracking[carrier]["near_band_examples"], ensure_ascii=False)})
            summary["materials"].setdefault(mat, {})[carrier] = report

            # variant rows
            for name, v in variants.items():
                if not v.get("available"):
                    variant_rows.append({"material": mat, "carrier": carrier, "variant": name,
                                         "available": False, "reason": v.get("reason", "")})
                    continue
                variant_rows.append({
                    "material": mat, "carrier": carrier, "variant": name,
                    "available": True, "n_points": v.get("n_points"),
                    "linear_term": v.get("linear_term"), "fractional_fit": v.get("fractional_fit"),
                    "mass_1_me": v["principal_masses_me"][0] if "principal_masses_me" in v else None,
                    "mass_2_me": v["principal_masses_me"][1] if "principal_masses_me" in v else None,
                    "mass_3_me": v["principal_masses_me"][2] if "principal_masses_me" in v else None,
                    "hessian_eigenvalues_eV_A2": json.dumps(v.get("hessian_eigenvalues_eV_A2")),
                    "gradient_norm_eV_A": v.get("gradient_norm_eV_A"),
                    "rmse_eV": v.get("rmse_eV"),
                    "max_abs_residual_eV": v.get("max_abs_residual_eV"),
                    "condition_number": v.get("condition_number"),
                    "stationary_shift_Ainv": json.dumps(v.get("stationary_shift_Ainv")),
                    "stationary_inside_box": v.get("stationary_inside_box"),
                    "curvature_sign_ok": v.get("curvature_sign_ok"),
                    "axis_curvature_x_eV_A2": (v.get("curvature_eV_A2") or {}).get("x") if name == "axis_1d_curvature" else None,
                    "axis_curvature_y_eV_A2": (v.get("curvature_eV_A2") or {}).get("y") if name == "axis_1d_curvature" else None,
                    "axis_curvature_z_eV_A2": (v.get("curvature_eV_A2") or {}).get("z") if name == "axis_1d_curvature" else None,
                    "note": v.get("note", ""),
                })

            # recheck row
            recheck_rows.append({
                "material": mat, "carrier": carrier,
                "band_index_1based": band + 1,
                "center_fractional": json.dumps(meta["edges"][carrier]["center_fractional"]),
                "metadata_edge_energy_eV": meta["edges"][carrier]["energy_eV"],
                "local_center_energy_eV": center_energy,
                "center_minus_metadata_meV": 1000.0 * report["center_minus_metadata_edge_eV"],
                "csv_vs_eigenval_max_abs_eV": csv_vs_eig,
                "mass_1_me": masses[0], "mass_2_me": masses[1], "mass_3_me": masses[2],
                "mass_ordering": "ascending electron-energy Hessian eigenvalues; masses mapped in the same order",
                "gradient_norm_eV_A": report["gradient_norm_eV_A"],
                "rmse_eV": report["rmse_eV"],
                "relative_rmse_to_span": report["relative_rmse_to_span"],
                "curvature_sign_ok": report["curvature_sign_ok"],
                "stationary_shift_Ainv": json.dumps(report["quadratic_stationary_shift_Ainv"]),
                "stationary_inside_box": inside,
                "sampled_extremum_q_Ainv": json.dumps(actual["q_Ainv"]),
                "sampled_extremum_grid_index": json.dumps(points[actual["index"]]["grid_index"]),
                "sampled_extremum_is_on_boundary": actual["on_boundary"],
                "sampled_extremum_distance_from_center_Ainv": actual["distance_from_center_Ainv"],
                "center_is_sampled_extremum": report["center_is_sampled_extremum"],
                "energy_gain_center_to_extremum_meV": report["energy_gain_center_to_extremum_meV"],
                "quadratic_stationary_energy_minus_actual_meV": report["quadratic_stationary_energy_minus_actual_meV"],
                "selected_band_is_edge_band_all_points": tracking[carrier]["selected_band_is_edge_band_all_points"],
                "n_points_selected_not_edge_band": len(points) - tracking[carrier]["n_selected_is_edge_band"],
                "degenerate_partner_gap_min_eV": tracking[carrier]["degenerate_partner_gap_min_eV"],
                "degenerate_partner_gap_max_eV": tracking[carrier]["degenerate_partner_gap_max_eV"],
                "n_points_other_band_within_50meV": tracking[carrier]["n_points_with_other_band_within_50meV"],
                "old_masses_me": json.dumps(old["reports"][carrier].get("principal_masses_me")),
                "old_warnings": json.dumps(old["reports"][carrier].get("warnings"), ensure_ascii=False),
                "reliability": reliability,
                "reliability_reasons": ";".join(reasons),
            })

        summary["materials"][mat]["local_mesh"] = {
            "kpoints_check": vcheck,
            "n_local_kpoints": hl["nk"], "n_local_bands": hl["nb"],
            "mass_metadata_grid": meta.get("grid"),
            "mass_metadata_step_Ainv": meta.get("step_cartesian_Ainv"),
            "mass_tensor_input_hashes": old.get("input_files"),
        }
        print(f"[{mat}] dense gap {de['gap_eV']:.6f} eV; "
              f"hole {summary['materials'][mat]['hole']['reliability']}; "
              f"electron {summary['materials'][mat]['electron']['reliability']}")

    ac.write_csv(ac.CSV / "mass_fit_recheck.csv", recheck_rows)
    ac.write_csv(ac.CSV / "mass_fit_variants.csv", variant_rows)
    ac.write_csv(ac.CSV / "mass_variant_summary.csv", variant_summary_rows)
    ac.write_csv(ac.CSV / "band_edge_check.csv", edge_rows)
    ac.write_csv(ac.CSV / "band_tracking.csv", tracking_rows)
    ac.write_csv(ac.CSV / "mass_point_diagnostics.csv", point_rows)
    ac.write_json(ac.JSON / "mass_audit_summary.json", summary)
    print("wrote", ac.CSV / "mass_fit_recheck.csv")
    print("wrote", ac.CSV / "mass_fit_variants.csv")
    print("wrote", ac.CSV / "mass_point_diagnostics.csv")


if __name__ == "__main__":
    main()
