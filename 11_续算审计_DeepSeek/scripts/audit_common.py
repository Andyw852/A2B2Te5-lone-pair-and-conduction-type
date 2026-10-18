#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared, read-only parsers and independent fitting helpers for the A2B2Te5 audit.

This module deliberately does NOT import the old fit_effective_mass.py.  The
quadratic fit and all unit conversions are re-implemented here so that the
reported numbers are an independent check, not a replay of the production code.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HBAR2_OVER_ME_EV_A2 = 7.619964
MATERIALS = ["Pb2Bi2Te5", "Pb2Sb2Te5", "Sn2Bi2Te5", "Sn2Sb2Te5"]
A_SITE = {"Pb2Bi2Te5": "Pb", "Pb2Sb2Te5": "Pb", "Sn2Bi2Te5": "Sn", "Sn2Sb2Te5": "Sn"}
B_SITE = {"Pb2Bi2Te5": "Bi", "Pb2Sb2Te5": "Sb", "Sn2Bi2Te5": "Bi", "Sn2Sb2Te5": "Sb"}

WORK = Path(__file__).resolve().parents[2]
DENSE = WORK / "A2B2Te5_SOC_dense"
RESULT = DENSE / "materials"
OUT = WORK / "11_续算审计_DeepSeek"
CSV = OUT / "csv"
JSON = OUT / "json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_eigenval(path: Path):
    """Read a VASP EIGENVAL (non-collinear or collinear) without external libs.

    Returns (kpoints_frac, energies, occupations, header) where energies and
    occupations have shape (nk, nb) for the non-collinear/ISPIN=1 case and
    (ispin, nk, nb) otherwise.  The A2B2Te5 files are non-collinear.
    """
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"EIGENVAL too short: {path}")
    first = [int(x) for x in lines[0].split()[:4]]
    ispin = first[3] if len(first) >= 4 else 1
    header_fields = lines[5].split()
    nelect, nk, nb = int(header_fields[0]), int(header_fields[1]), int(header_fields[2])
    kpts = np.zeros((nk, 3), dtype=float)
    eig = np.zeros((ispin, nk, nb), dtype=float)
    occ = np.zeros((ispin, nk, nb), dtype=float)
    i = 6
    for ik in range(nk):
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            raise ValueError(f"EIGENVAL ended before k-point {ik}")
        kpts[ik] = [float(x) for x in lines[i].split()[:3]]
        i += 1
        for ib in range(nb):
            if i >= len(lines):
                raise ValueError(f"EIGENVAL ended before band {ib} of k {ik}")
            f = lines[i].split()
            # non-collinear ISPIN=1: band E occ
            # collinear ISPIN=2: band E_up E_down occ_up occ_down
            if ispin == 1:
                eig[0, ik, ib] = float(f[1])
                occ[0, ik, ib] = float(f[2]) if len(f) > 2 else 1.0
            else:
                eig[:, ik, ib] = [float(x) for x in f[1:1 + ispin]]
                occ[:, ik, ib] = [float(x) for x in f[1 + ispin:1 + 2 * ispin]]
            i += 1
    header = {"nelect": nelect, "nk": nk, "nb": nb, "ispin": ispin,
              "first_line": first, "raw_header": header_fields}
    if ispin == 1:
        return kpts, eig[0], occ[0], header
    return kpts, eig, occ, header


def read_poscar(path: Path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    scale = float(lines[1].split()[0])
    lattice = np.asarray([[float(x) for x in lines[j].split()[:3]] for j in (2, 3, 4)], dtype=float)
    if scale < 0:
        scale = (-scale / abs(np.linalg.det(lattice))) ** (1.0 / 3.0)
    lattice = lattice * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    volume = abs(float(np.linalg.det(lattice)))
    return {"lattice_rows_A": lattice, "species": species, "counts": counts,
            "natoms": int(sum(counts)), "volume_A3": volume,
            "scale": scale, "title": lines[0].strip()}


def read_kpoints(path: Path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    title = lines[0].strip()
    nk = int(lines[1].split()[0])
    mode = lines[2].strip()
    if mode.lower().startswith("g") and nk == 0:
        mesh = [int(x) for x in lines[3].split()[:3]]
        shift = [float(x) for x in lines[4].split()[:3]]
        return {"title": title, "mode": "Gamma", "mesh": mesh, "shift": shift,
                "kpoints": None, "nk": int(np.prod(mesh))}
    kp = np.zeros((nk, 3), dtype=float)
    weights = np.zeros(nk, dtype=float)
    for i in range(nk):
        f = lines[3 + i].split()
        kp[i] = [float(x) for x in f[:3]]
        weights[i] = float(f[3]) if len(f) > 3 else 1.0
    return {"title": title, "mode": mode, "mesh": None, "shift": None,
            "kpoints": kp, "weights": weights, "nk": nk}


def reciprocal_rows(lattice_rows_A: np.ndarray) -> np.ndarray:
    """Rows are b_i with a_i . b_j = 2 pi delta_ij."""
    return 2.0 * math.pi * np.linalg.inv(lattice_rows_A).T


def read_potcar_zval(path: Path):
    zvals = []
    titles = []
    for line in path.read_text(errors="replace").splitlines():
        if "TITEL" in line:
            titles.append(line.split("=", 1)[1].strip() if "=" in line else line.strip())
        if "ZVAL" in line and "=" in line:
            # POMASS = ...; ZVAL = 14.000
            part = line.split("ZVAL", 1)[1]
            val = part.split("=", 1)[1].split()[0]
            zvals.append(float(val))
    return zvals, titles


def read_outcar_nelect(path: Path):
    for line in path.read_text(errors="replace").splitlines():
        if "NELECT" in line and "=" in line:
            try:
                return float(line.split("=", 1)[1].split()[0])
            except ValueError:
                pass
    return None


def read_json(path: Path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_csv(path: Path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def fit_quadratic(q: np.ndarray, energy: np.ndarray, linear: bool = True):
    """Independent E(q)=E0+g.q+0.5 q^T H q fit, q in Cartesian Angstrom^-1.

    Coordinates are scaled only for conditioning; H is converted back to eV A^2.
    """
    q = np.asarray(q, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if q.shape != (len(energy), 3):
        raise ValueError("q must be (N,3)")
    scale = float(np.max(np.abs(q)))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("zero q extent")
    x, y, z = (q / scale).T
    if linear:
        design = np.column_stack([np.ones(len(q)), x, y, z,
                                  x * x / 2, y * y / 2, z * z / 2,
                                  x * y, x * z, y * z])
        n_par = 10
    else:
        design = np.column_stack([np.ones(len(q)),
                                  x * x / 2, y * y / 2, z * z / 2,
                                  x * y, x * z, y * z])
        n_par = 7
    if len(q) < n_par:
        raise ValueError(f"need at least {n_par} points, got {len(q)}")
    reference = float(np.mean(energy))
    beta, residuals, rank, singular = np.linalg.lstsq(design, energy - reference, rcond=None)
    if rank != n_par:
        raise ValueError(f"rank deficient fit: {rank}/{n_par}")
    if linear:
        g = beta[1:4] / scale
        h = np.array([[beta[4], beta[7], beta[8]],
                      [beta[7], beta[5], beta[9]],
                      [beta[8], beta[9], beta[6]]]) / scale ** 2
        intercept = reference + beta[0]
    else:
        g = np.zeros(3)
        h = np.array([[beta[1], beta[4], beta[5]],
                      [beta[4], beta[2], beta[6]],
                      [beta[5], beta[6], beta[3]]]) / scale ** 2
        intercept = reference + beta[0]
    predicted = reference + design @ beta
    residual = energy - predicted
    eigval, eigvec = np.linalg.eigh(h)
    return {
        "intercept_eV": float(intercept),
        "gradient_eV_A": g,
        "hessian_eV_A2": h,
        "hessian_eigenvalues": eigval,
        "hessian_eigenvectors": eigvec,
        "rmse_eV": float(np.sqrt(np.mean(residual ** 2))),
        "max_abs_residual_eV": float(np.max(np.abs(residual))),
        "energy_span_eV": float(np.ptp(energy)),
        "condition_number": float(singular[0] / singular[-1]),
        "n_points": int(len(q)),
        "n_parameters": int(n_par),
        "dof": int(len(q) - n_par),
        "predicted_eV": predicted,
        "residual_eV": residual,
    }


def principal_masses(hessian: np.ndarray, carrier: str):
    """Masses in the production-code ordering.

    The eigenvectors/values are those of the *electron-energy* Hessian in
    ascending eigenvalue order.  For a hole, the physical curvature is -H, so
    a more negative H eigenvalue gives a smaller hole mass.  This deliberately
    reproduces the old mass_tensor.json ordering so that old and new values can
    be compared element by element; it is not an ascending-mass ordering.
    """
    sign = 1.0 if carrier == "electron" else -1.0
    eigval, eigvec = np.linalg.eigh(np.asarray(hessian, dtype=float))
    curvature = sign * eigval
    masses = []
    for value in curvature:
        masses.append(float(HBAR2_OVER_ME_EV_A2 / value) if value > 1e-8 else None)
    return curvature, eigvec, masses


def stationary_shift(hessian: np.ndarray, gradient: np.ndarray):
    try:
        shift = -np.linalg.solve(hessian, gradient)
    except np.linalg.LinAlgError:
        return None
    return shift


def inside_box(q: np.ndarray, point: np.ndarray, tol: float = 1e-12):
    return bool(np.all(point >= q.min(axis=0) - tol) and np.all(point <= q.max(axis=0) + tol))


def actual_extremum(q: np.ndarray, energy: np.ndarray, carrier: str):
    idx = int(np.argmin(energy) if carrier == "electron" else np.argmax(energy))
    qmax = np.max(np.abs(q), axis=0)
    on_boundary = bool(np.any(np.abs(q[idx]) >= qmax - 1e-12))
    return {"index": idx, "energy_eV": float(energy[idx]), "q_Ainv": q[idx].tolist(),
            "on_boundary": on_boundary,
            "max_abs_q_Ainv": float(np.max(np.abs(q[idx]))),
            "distance_from_center_Ainv": float(np.linalg.norm(q[idx]))}


def axis_curvatures(q: np.ndarray, energy: np.ndarray, tol: float = 1e-12):
    """1D quadratic curvature along each Cartesian axis through q=0."""
    out = {}
    for axis, name in enumerate("xyz"):
        others = [i for i in range(3) if i != axis]
        mask = np.all(np.abs(q[:, others]) <= tol, axis=1)
        if int(mask.sum()) < 3:
            out[name] = None
            continue
        x = q[mask, axis]
        order = np.argsort(x)
        x = x[order]
        y = energy[mask][order]
        coef = np.polyfit(x, y, 2)
        out[name] = float(2.0 * coef[0])
    return out


def fmt(value, nd=6):
    if value is None:
        return ""
    if isinstance(value, (bool, np.bool_)):
        return "true" if value else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        if not np.isfinite(value):
            return ""
        return f"{float(value):.{nd}g}"
    return str(value)
