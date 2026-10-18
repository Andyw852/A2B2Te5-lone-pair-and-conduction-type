#!/usr/bin/env python3
"""Generate a Cartesian local k mesh for a candidate band extremum.

Requires numpy. Does not locate extrema, assign bands, run VASP, or fit masses.
The output KPOINTS_OPT requires VASP >= 6.3 and a separate uniform SCF KPOINTS.
Keep the same structure, functional, SOC and pseudopotentials as the reference.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np


def read_lattice(path: Path) -> np.ndarray:
    """Read lattice rows in Angstrom, including all POSCAR scale conventions."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 5:
        raise ValueError("POSCAR must have at least five lines")
    scales = np.array([float(x) for x in lines[1].split()])
    lattice = np.array([[float(x) for x in line.split()[:3]] for line in lines[2:5]])
    if lattice.shape != (3, 3) or not np.isfinite(lattice).all():
        raise ValueError("Invalid lattice vectors")
    volume = abs(float(np.linalg.det(lattice)))
    if volume < 1e-12 or not np.isfinite(scales).all():
        raise ValueError("Singular lattice or invalid scale")
    if len(scales) == 1 and scales[0] != 0:
        scale = scales[0] if scales[0] > 0 else (-scales[0] / volume) ** (1 / 3)
        lattice *= scale
    elif len(scales) == 3 and np.all(scales > 0):
        lattice *= scales[None, :]
    else:
        raise ValueError("Scale must be one nonzero number or three positive numbers")
    return lattice


def make_mesh(lattice: np.ndarray, center: np.ndarray, step: float, grid: int):
    """Reciprocal rows B = 2*pi*inv(A).T; k_cart = k_frac @ B."""
    if grid < 3 or grid % 2 != 1:
        raise ValueError("--grid must be an odd integer >= 3")
    if not np.isfinite(step) or step <= 0:
        raise ValueError("--step must be finite and positive (Angstrom^-1)")
    if np.shape(center) != (3,) or not np.isfinite(center).all():
        raise ValueError("--center must contain three finite reciprocal fractions")
    reciprocal = 2 * np.pi * np.linalg.inv(lattice).T
    half = grid // 2
    indices = np.array(list(itertools.product(range(-half, half + 1), repeat=3)))
    offsets = indices * step
    fractional = center[None, :] + offsets @ np.linalg.inv(reciprocal)
    return reciprocal, indices, offsets, fractional


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poscar", required=True, type=Path)
    parser.add_argument("--center", required=True, nargs=3, type=float, metavar=("KX", "KY", "KZ"),
                        help="Candidate extremum in fractional reciprocal coordinates")
    parser.add_argument("--step", type=float, default=0.005,
                        help="Cartesian spacing in Angstrom^-1 (default: 0.005)")
    parser.add_argument("--grid", type=int, default=5, help="Odd grid length on each Cartesian axis")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--force", action="store_true", help="Overwrite these two output files only")
    args = parser.parse_args()
    try:
        source = args.poscar.resolve(strict=True)
        lattice = read_lattice(source)
        center = np.array(args.center)
        reciprocal, indices, offsets, fractional = make_mesh(lattice, center, args.step, args.grid)
        targets = [args.output / "KPOINTS_OPT", args.output / "mass_kpoints_metadata.json"]
        if not args.force and any(path.exists() for path in targets):
            raise ValueError("Output already exists; choose a new directory or pass --force")
        if any(path.resolve() == source for path in targets):
            raise ValueError("Output must not overwrite the source POSCAR")
        # Do not wrap: the continuous fractional coordinates retain the local
        # Cartesian displacements even when the mesh crosses a zone boundary.
        records = []
        for i, (idx, cart, frac) in enumerate(zip(indices, offsets, fractional), start=1):
            records.append({"kpoint_index_1based": i, "grid_index": idx.tolist(),
                            "offset_cartesian_Ainv": cart.tolist(), "k_fractional_unwrapped": frac.tolist()})
        metadata = {
            "purpose": "Local three-dimensional sampling around a USER-SUPPLIED candidate extremum",
            "poscar": str(source), "poscar_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "lattice_rows_A": lattice.tolist(), "reciprocal_rows_Ainv_including_2pi": reciprocal.tolist(),
            "center_fractional": center.tolist(), "step_cartesian_Ainv": args.step,
            "grid": args.grid, "number_of_points": len(records),
            "center_index_1based": next(r["kpoint_index_1based"] for r in records if r["grid_index"] == [0, 0, 0]),
            "coordinate_convention": "k_cart = k_frac @ reciprocal_rows; no wrapping",
            "cautions": [
                "The candidate center is not automatically a global or local extremum.",
                "KPOINTS_OPT is an additional sampling file, not a replacement for the uniform SCF KPOINTS.",
                "Track band character near crossings; do not blindly fit a fixed band number.",
                "Check step/window convergence, nonparabolicity and residual linear term before reporting masses.",
                "A local mesh does not establish valley degeneracy or a full-BZ DOS mass."
            ],
            "points": records,
        }
        output = [f"Local Cartesian mass mesh: grid={args.grid}, step={args.step:.8g} A^-1",
                  str(len(records)), "Reciprocal"]
        output += [" ".join(f"{x:+.15f}" for x in frac) + " 1.0" for frac in fractional]
        args.output.mkdir(parents=True, exist_ok=True)
        targets[0].write_text("\n".join(output) + "\n", encoding="utf-8")
        targets[1].write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote {len(records)} sampling points to {targets[0].resolve()}")
        print(f"Metadata: {targets[1].resolve()}")
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
