#!/usr/bin/env python3
"""Generate local SOC band-edge sampling inputs from the dense-grid EIGENVAL."""
from pathlib import Path
import argparse
import itertools
import json
import math
import shutil

import numpy as np

CFG = {
    "Pb2Bi2Te5": {"nbands": 144, "encut": 370},
    "Pb2Sb2Te5": {"nbands": 120, "encut": 360},
    "Sn2Bi2Te5": {"nbands": 144, "encut": 370},
    "Sn2Sb2Te5": {"nbands": 120, "encut": 370},
}


def read_eigenval(path):
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError("EIGENVAL is too short")
    ispin = int(lines[0].split()[3])
    nk, nb = (int(x) for x in lines[5].split()[1:3])
    kpts = []
    eig = np.zeros((ispin, nk, nb), dtype=float)
    i = 6
    for k in range(nk):
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            raise ValueError("EIGENVAL ended before all k-points")
        kpts.append([float(x) for x in lines[i].split()[:3]])
        i += 1
        for b in range(nb):
            fields = lines[i].split()
            eig[:, k, b] = [float(x) for x in fields[1:1 + ispin]]
            i += 1
    occ = np.zeros((ispin, nk, nb), dtype=float)
    i = 6
    for k in range(nk):
        while not lines[i].strip():
            i += 1
        i += 1
        for b in range(nb):
            fields = lines[i].split()
            if len(fields) > 1 + ispin:
                occ[:, k, b] = [float(x) for x in fields[1 + ispin:1 + 2 * ispin]]
            i += 1
    return np.asarray(kpts), eig, occ


def read_lattice(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    scale = float(lines[1].split()[0])
    lattice = np.asarray([[float(x) for x in lines[j].split()[:3]] for j in (2, 3, 4)])
    if scale < 0:
        scale = (-scale / abs(np.linalg.det(lattice))) ** (1.0 / 3.0)
    return lattice * scale


def main():
    cwd = Path.cwd()
    parser = argparse.ArgumentParser()
    parser.add_argument("--material", required=True)
    mat = parser.parse_args().material
    if mat not in CFG:
        raise SystemExit("[ERROR] unsupported material: %s" % mat)
    source = cwd / "step3_uniform" / "EIGENVAL"
    if not source.is_file():
        raise SystemExit("[ERROR] dense step3_uniform/EIGENVAL is missing")
    kfrac, eig, occ = read_eigenval(source)
    e = eig[0]
    o = occ[0]
    occupied = o > 0.5
    if not occupied.any() or occupied.all():
        raise SystemExit("[ERROR] cannot locate both occupied and empty states")
    vbm_flat = np.where(occupied, e, -np.inf)
    cbm_flat = np.where(occupied, np.inf, e)
    v_idx = np.unravel_index(int(np.argmax(vbm_flat)), vbm_flat.shape)
    c_idx = np.unravel_index(int(np.argmin(cbm_flat)), cbm_flat.shape)
    edges = {
        "hole": {"k_index": int(v_idx[0]), "band_index_1based": int(v_idx[1] + 1),
                 "energy_eV": float(e[v_idx]), "center_fractional": kfrac[v_idx[0]].tolist()},
        "electron": {"k_index": int(c_idx[0]), "band_index_1based": int(c_idx[1] + 1),
                      "energy_eV": float(e[c_idx]), "center_fractional": kfrac[c_idx[0]].tolist()},
    }
    lattice = read_lattice(cwd / "POSCAR")
    reciprocal = 2.0 * math.pi * np.linalg.inv(lattice).T
    inv_reciprocal = np.linalg.inv(reciprocal)
    step = 0.005
    grid = 5
    half = grid // 2
    indices = np.asarray(list(itertools.product(range(-half, half + 1), repeat=3)), dtype=int)
    points_by_edge = {}
    for carrier in ("electron", "hole"):
        points = []
        center = np.asarray(edges[carrier]["center_fractional"], dtype=float)
        for idx in indices:
            offset = idx.astype(float) * step
            frac = center + offset @ inv_reciprocal
            points.append({"grid_index": idx.tolist(), "offset_cartesian_Ainv": offset.tolist(),
                           "k_fractional_unwrapped": frac.tolist()})
        points_by_edge[carrier] = points
    out = cwd / "step8.2_mass"
    out.mkdir(exist_ok=True)
    all_points = points_by_edge["electron"] + points_by_edge["hole"]
    kp = ["Local SOC mass meshes: electron then hole", str(len(all_points)), "Reciprocal"]
    kp += [" ".join("%+.15f" % x for x in p["k_fractional_unwrapped"]) + " 1.0" for p in all_points]
    (out / "KPOINTS").write_text("\n".join(kp) + "\n", encoding="utf-8")
    meta = {"material": mat, "grid": grid, "step_cartesian_Ainv": step,
            "lattice_rows_A": lattice.tolist(), "reciprocal_rows_Ainv": reciprocal.tolist(),
            "points_by_edge": points_by_edge, "edges": edges,
            "band_tracking": "fixed band index selected from dense EIGENVAL; inspect PROCAR near crossings"}
    (out / "mass_metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    shutil.copy2(cwd / "POSCAR", out / "POSCAR")
    shutil.copy2(cwd / "POTCAR", out / "POTCAR")
    c = CFG[mat]
    (out / "INCAR").write_text(
        """SYSTEM = %s local SOC mass
ISTART = 0
ICHARG = 11
GGA = PE
IVDW = 12
PREC = Accurate
ENCUT = %d
LREAL = .FALSE.
LASPH = .TRUE.
ALGO = Normal
EDIFF = 1E-8
NELM = 100
ISMEAR = 0
SIGMA = 0.01
IBRION = -1
NSW = 0
ISYM = -1
LWAVE = .FALSE.
LCHARG = .FALSE.
LORBIT = 11
NCORE = 4
KPAR = 2
LSORBIT = .TRUE.
LNONCOLLINEAR = .TRUE.
SAXIS = 0 0 1
MAGMOM = 27*0
NBANDS = %d
""" % (mat, c["encut"], c["nbands"]), encoding="utf-8")
    (out / "submit.sh").write_text(
        """#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name=%s-S8.2_mass
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular
cd $SLURM_SUBMIT_DIR
cp ../step3_uniform/CHGCAR CHGCAR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMPI_MCA_btl=^openib
export UCX_TLS=sm,self
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_ncl
""" % mat, encoding="utf-8")
    (out / "submit.sh").chmod(0o755)
    print("[DONE] %s: local mass mesh generated; VBM/CBM bands are in mass_metadata.json" % mat)


if __name__ == "__main__":
    main()
