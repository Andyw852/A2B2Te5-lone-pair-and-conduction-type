#!/usr/bin/env python3
"""Fit a local effective-mass tensor for ONE user-tracked band.

Inputs are the metadata from make_mass_kpoints.py and a CSV containing
kpoint_index_1based,energy_eV. Requires numpy. Does not select or track bands.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HBAR2_OVER_ME_EV_A2 = 7.619964


def read_data(metadata_path: Path, energies_path: Path):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8-sig"))
    records = metadata["points"]
    ids = [int(point["kpoint_index_1based"]) for point in records]
    if not records or len(set(ids)) != len(ids) or min(ids) < 1:
        raise ValueError("Metadata contains duplicate, invalid or missing point indices")
    q = np.asarray([point["offset_cartesian_Ainv"] for point in records], dtype=float)
    if q.shape != (len(records), 3) or not np.isfinite(q).all():
        raise ValueError("Metadata must contain finite Cartesian offsets in Angstrom^-1")
    energy_map = {}
    with energies_path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"kpoint_index_1based", "energy_eV"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("CSV requires columns kpoint_index_1based,energy_eV")
        for row in reader:
            key = int(row["kpoint_index_1based"])
            energy = float(row["energy_eV"])
            if key in energy_map:
                raise ValueError(f"Duplicate CSV k-point index {key}")
            if not np.isfinite(energy):
                raise ValueError(f"Nonfinite energy at CSV k-point index {key}")
            energy_map[key] = energy
    missing = sorted(set(ids) - energy_map.keys())
    extra = sorted(energy_map.keys() - set(ids))
    if missing or extra:
        raise ValueError(f"CSV and metadata indices must match exactly; missing={missing[:10]}, extra={extra[:10]}")
    return metadata, ids, q, np.asarray([energy_map[key] for key in ids])


def quadratic_fit(q: np.ndarray, energy: np.ndarray):
    """Scale coordinates for conditioning, then return physical Cartesian H."""
    if len(q) < 10:
        raise ValueError("At least ten independent points are required for a full 3D quadratic")
    scale = float(np.max(abs(q)))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Sampling offsets have zero extent")
    x, y, z = (q / scale).T
    design = np.column_stack([np.ones(len(q)), x, y, z,
                              x*x/2, y*y/2, z*z/2, x*y, x*z, y*z])
    # Subtract the mean to avoid losing small band-edge differences on a large
    # arbitrary absolute energy reference.
    reference = float(np.mean(energy))
    beta, _, rank, singular = np.linalg.lstsq(design, energy - reference, rcond=None)
    if rank != 10:
        raise ValueError(f"Sampling does not determine a full 3D Hessian (rank={rank}/10)")
    hessian = np.array([[beta[4], beta[7], beta[8]],
                        [beta[7], beta[5], beta[9]],
                        [beta[8], beta[9], beta[6]]]) / scale**2
    predicted = reference + design @ beta
    residual = energy - predicted
    eigenvalues, eigenvectors = np.linalg.eigh(hessian)
    return {
        "intercept_eV": reference + float(beta[0]),
        "gradient_eV_A": beta[1:4] / scale,
        "hessian_eV_A2": hessian,
        "hessian_eigenvalues_eV_A2": eigenvalues,
        "principal_axes_cartesian_columns": eigenvectors,
        "residual_eV": residual,
        "predicted_eV": predicted,
        "rmse_eV": float(np.sqrt(np.mean(residual**2))),
        "max_absolute_residual_eV": float(np.max(abs(residual))),
        "energy_span_eV": float(np.ptp(energy)),
        "scaled_design_condition_number": float(singular[0]/singular[-1]),
        "n_points": len(q), "n_parameters": 10, "degrees_of_freedom": len(q)-10,
    }


def analyze(q, energy, carrier, gradient_threshold=0.001,
            residual_threshold=0.0001, relative_residual_threshold=0.01,
            curvature_change_threshold=0.10):
    fit = quadratic_fit(q, energy)
    sign = 1 if carrier == "electron" else -1
    curvature = sign * fit["hessian_eigenvalues_eV_A2"]
    warnings = []
    if np.any(curvature <= 0):
        warnings.append("WRONG_CURVATURE_SIGN: the fitted point is not a strict minimum for electrons / maximum for holes; do not report a physical mass tensor.")
    small = np.abs(fit["hessian_eigenvalues_eV_A2"]) < 1e-8
    if np.any(small):
        warnings.append("NEAR_ZERO_CURVATURE: at least one mass is numerically ill-conditioned or divergent.")
    masses = [float(HBAR2_OVER_ME_EV_A2 / value) if value > 1e-8 else None for value in curvature]
    gnorm = float(np.linalg.norm(fit["gradient_eV_A"]))
    if gnorm > gradient_threshold:
        warnings.append("NONZERO_GRADIENT: the selected center is not stationary within the configured gradient threshold; refine the extremum.")
    stationary_shift = None
    inside_box = None
    if not np.any(small):
        stationary_shift = -np.linalg.solve(fit["hessian_eV_A2"], fit["gradient_eV_A"])
        inside_box = bool(np.all(stationary_shift >= q.min(axis=0)) and np.all(stationary_shift <= q.max(axis=0)))
        if not inside_box:
            warnings.append("STATIONARY_POINT_OUTSIDE_MESH: the quadratic stationary point is outside the sampled Cartesian box; recenter and resample.")
    span = fit["energy_span_eV"]
    relative_rmse = fit["rmse_eV"] / span if span > 0 else None
    if fit["rmse_eV"] > max(residual_threshold, relative_residual_threshold * span):
        warnings.append("NONPARABOLIC_OR_NOISY: quadratic residual exceeds both configured absolute and relative thresholds; inspect band tracking, numerical noise and nonparabolicity.")
    if fit["degrees_of_freedom"] == 0:
        warnings.append("NO_RESIDUAL_DEGREES_OF_FREEDOM: ten points determine ten parameters but cannot test quadratic fit quality.")
    if fit["scaled_design_condition_number"] > 1e8:
        warnings.append("ILL_CONDITIONED_DESIGN: sample geometry is poorly conditioned for a 3D quadratic fit.")
    # Compare the central Cartesian cube with the full mesh where possible.
    # This is an internal window diagnostic, not k-point or DFT convergence.
    max_axis = np.max(abs(q), axis=0)
    inner_mask = np.all(abs(q) < (max_axis[None, :] - 1e-12), axis=1)
    inner = {"available": False, "reason": "Fewer than ten independent inner points"}
    if int(inner_mask.sum()) >= 10:
        try:
            inner_fit = quadratic_fit(q[inner_mask], energy[inner_mask])
            denominator = np.maximum(np.abs(fit["hessian_eigenvalues_eV_A2"]), 1e-8)
            relative_changes = np.abs(inner_fit["hessian_eigenvalues_eV_A2"] - fit["hessian_eigenvalues_eV_A2"]) / denominator
            inner = {"available": True, "n_points": int(inner_mask.sum()),
                     "hessian_eigenvalues_eV_A2": inner_fit["hessian_eigenvalues_eV_A2"].tolist(),
                     "relative_eigenvalue_changes_vs_full": relative_changes.tolist(),
                     "rmse_eV": inner_fit["rmse_eV"]}
            if np.max(relative_changes) > curvature_change_threshold:
                warnings.append("WINDOW_DEPENDENT_CURVATURE: inner and full meshes give different Hessian eigenvalues; check nonparabolicity / noise and repeat sampling at smaller steps.")
        except ValueError as exc:
            inner["reason"] = str(exc)
    report = {key: value.tolist() if isinstance(value, np.ndarray) else value
              for key, value in fit.items() if key not in ("residual_eV", "predicted_eV")}
    report.update({
        "carrier": carrier,
        "model": "E(q)=E0+g.q+0.5*q.T*H*q; q in Cartesian Angstrom^-1",
        "mass_conversion_hbar2_over_me_eV_A2": HBAR2_OVER_ME_EV_A2,
        "principal_masses_me": masses,
        "all_curvatures_have_required_sign": bool(np.all(curvature > 1e-8)),
        "axis_order": "Ascending eigenvalues of the ELECTRON ENERGY Hessian; columns of principal_axes_cartesian_columns correspond to masses in the same order.",
        "gradient_norm_eV_A": gnorm,
        "quadratic_stationary_shift_Ainv": stationary_shift.tolist() if stationary_shift is not None else None,
        "quadratic_stationary_point_inside_sample_box": inside_box,
        "relative_rmse_to_energy_span": relative_rmse,
        "inner_window_diagnostic": inner,
        "thresholds": {"gradient_eV_A": gradient_threshold,
                       "absolute_rmse_eV": residual_threshold,
                       "relative_rmse": relative_residual_threshold,
                       "relative_inner_window_curvature_change": curvature_change_threshold},
        "warnings": warnings,
        "scope": "Local curvature for one user-tracked band only. No automatic band tracking, global-extremum proof, valley degeneracy, DOS mass, scattering or material validation.",
    })
    return report, fit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--energies", required=True, type=Path)
    parser.add_argument("--carrier", required=True, choices=["electron", "hole"])
    parser.add_argument("--output", required=True, type=Path, help="Output JSON filename")
    parser.add_argument("--gradient-threshold", type=float, default=0.001, help="eV Angstrom")
    parser.add_argument("--residual-threshold", type=float, default=0.0001, help="Absolute RMSE threshold in eV")
    parser.add_argument("--relative-residual-threshold", type=float, default=0.01)
    parser.add_argument("--curvature-change-threshold", type=float, default=0.10)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    try:
        for name in ("gradient_threshold", "residual_threshold", "relative_residual_threshold", "curvature_change_threshold"):
            if not np.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
                raise ValueError(f"--{name.replace('_', '-')} must be finite and positive")
        sources = [args.metadata.resolve(strict=True), args.energies.resolve(strict=True)]
        if args.output.resolve() in sources:
            raise ValueError("Output must not overwrite an input file")
        if args.output.exists() and not args.force:
            raise ValueError("Output already exists; choose a new filename or pass --force")
        metadata, ids, q, energy = read_data(*sources)
        report, fit = analyze(q, energy, args.carrier, args.gradient_threshold,
                              args.residual_threshold, args.relative_residual_threshold,
                              args.curvature_change_threshold)
        report["input_files"] = [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
        report["center_fractional"] = metadata.get("center_fractional")
        report["point_residuals"] = [{"kpoint_index_1based": i, "energy_eV": float(e),
                                       "predicted_energy_eV": float(pred), "residual_eV": float(res)}
                                      for i, e, pred, res in zip(ids, energy, fit["predicted_eV"], fit["residual_eV"])]
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Saved local quadratic fit to {args.output.resolve()}")
        print(f"Principal {args.carrier} curvature masses / me: {report['principal_masses_me']}")
        print(f"RMSE = {report['rmse_eV']:.6g} eV; gradient norm = {report['gradient_norm_eV_A']:.6g} eV A")
        for warning in report["warnings"]:
            print("WARNING: " + warning)
    except (OSError, ValueError, KeyError, TypeError, np.linalg.LinAlgError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
