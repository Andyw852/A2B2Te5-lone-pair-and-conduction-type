#!/usr/bin/env python3
"""Generate the isolated 3-D dense SOC VASP input without touching old data."""
from pathlib import Path
import argparse
import json
import shutil

CFG = {
    "Pb2Bi2Te5": {"mesh": (28, 28, 6), "encut": 370, "nbands": 144},
    "Pb2Sb2Te5": {"mesh": (28, 28, 7), "encut": 360, "nbands": 120},
    "Sn2Bi2Te5": {"mesh": (28, 28, 6), "encut": 370, "nbands": 144},
    "Sn2Sb2Te5": {"mesh": (28, 28, 7), "encut": 370, "nbands": 120},
}

cwd = Path.cwd()
parser = argparse.ArgumentParser()
parser.add_argument("--material", required=True)
mat = parser.parse_args().material
if mat not in CFG:
    raise SystemExit("[ERROR] unsupported material: %s" % mat)
if not (cwd / "POSCAR").is_file() or not (cwd / "POTCAR").is_file():
    raise SystemExit("[ERROR] POSCAR/POTCAR are required in the material root")

cfg = CFG[mat]
out = cwd / "step3_uniform"
out.mkdir(exist_ok=True)
shutil.copy2(cwd / "POSCAR", out / "POSCAR")
shutil.copy2(cwd / "POTCAR", out / "POTCAR")

nx, ny, nz = cfg["mesh"]
(out / "KPOINTS").write_text(
    "Dense SOC Gamma mesh\n0\nGamma\n  %d  %d  %d\n0 0 0\n" % (nx, ny, nz),
    encoding="utf-8",
)

(out / "INCAR").write_text(
    """SYSTEM = %s dense SOC

ISTART = 0
ICHARG = 2
GGA    = PE
IVDW   = 12
PREC   = Accurate
ENCUT  = %d
LREAL  = .FALSE.
LASPH  = .TRUE.
ALGO   = Normal
EDIFF  = 1E-8
NELM   = 200
NELMIN = 6
ISMEAR = 0
SIGMA  = 0.01
IBRION = -1
NSW    = 0
ISYM   = -1
LWAVE  = .TRUE.
LCHARG = .TRUE.
LORBIT = 11
NEDOS  = 5000
LOPTICS = .FALSE.
NCORE  = 4
KPAR   = 2
LSORBIT = .TRUE.
LNONCOLLINEAR = .TRUE.
SAXIS = 0 0 1
MAGMOM = 27*0
NBANDS = %d
""" % (mat, cfg["encut"], cfg["nbands"]),
    encoding="utf-8",
)

(out / "submit.sh").write_text(
    """#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name=%s-S3_dense_SOC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular
cd $SLURM_SUBMIT_DIR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMPI_MCA_btl=^openib
export UCX_TLS=sm,self
export TF_CELL_CONSTRAINT=none
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_ncl
""" % mat,
    encoding="utf-8",
)
(out / "submit.sh").chmod(0o755)
(out / "workflow_method.txt").write_text(
    "DIM=3d\nMETHOD=PBE-D3(BJ)+SOC\nGRID=%dx%dx%d\n" % (nx, ny, nz),
    encoding="utf-8",
)
(out / "dense_grid_manifest.json").write_text(
    json.dumps({"material": mat, "mesh": [nx, ny, nz], "soc": True,
                "lcharg": True, "encut_eV": cfg["encut"],
                "nbands": cfg["nbands"]}, indent=2) + "\n",
    encoding="utf-8",
)
print("[DONE] %s: dense SOC input %dx%dx%d generated" % (mat, nx, ny, nz))
