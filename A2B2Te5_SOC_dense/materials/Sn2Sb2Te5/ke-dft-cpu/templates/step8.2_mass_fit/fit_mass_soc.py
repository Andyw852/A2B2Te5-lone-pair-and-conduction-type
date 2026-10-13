#!/usr/bin/env python3
"""Fit the two SOC band-edge mass tensors from the local VASP output."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fit_effective_mass import analyze  # noqa: E402


def read_eigenval(path):
    lines = path.read_text(errors="replace").splitlines()
    ispin = int(lines[0].split()[3])
    nk, nb = (int(x) for x in lines[5].split()[1:3])
    eig = np.zeros((ispin, nk, nb), dtype=float)
    i = 6
    for k in range(nk):
        while not lines[i].strip():
            i += 1
        i += 1
        for b in range(nb):
            fields = lines[i].split()
            eig[:, k, b] = [float(x) for x in fields[1:1 + ispin]]
            i += 1
    return eig


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


cwd = Path.cwd()
parser = argparse.ArgumentParser()
parser.add_argument("--material", required=True)
mat = parser.parse_args().material
meta_path = cwd / "step8.2_mass" / "mass_metadata.json"
eigen_path = cwd / "step8.2_mass" / "EIGENVAL"
if not meta_path.is_file() or not eigen_path.is_file():
    raise SystemExit("[ERROR] mass metadata or EIGENVAL is missing")
input_files = {
    "step8.2_mass/mass_metadata.json": sha256(meta_path),
    "step8.2_mass/EIGENVAL": sha256(eigen_path),
}
meta = json.loads(meta_path.read_text(encoding="utf-8"))
eig = read_eigenval(eigen_path)[0]
out = cwd / "step8.2_mass_fit"
out.mkdir(exist_ok=True)
reports = {}
for carrier in ("electron", "hole"):
    edge = meta["edges"][carrier]
    points = meta["points_by_edge"][carrier]
    start = 0 if carrier == "electron" else len(meta["points_by_edge"]["electron"])
    band = int(edge["band_index_1based"]) - 1
    energies = eig[start:start + len(points), band]
    if len(energies) != len(points):
        raise SystemExit("[ERROR] local EIGENVAL point count does not match metadata")
    md = {"center_fractional": edge["center_fractional"], "points": [
        {"kpoint_index_1based": i + 1,
         "offset_cartesian_Ainv": p["offset_cartesian_Ainv"]}
        for i, p in enumerate(points)]}
    md_path = out / (carrier + "_metadata.json")
    csv_path = out / (carrier + "_band.csv")
    md_path.write_text(json.dumps(md, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["kpoint_index_1based", "energy_eV"])
        writer.writerows((i + 1, "%.12f" % value) for i, value in enumerate(energies))
    q = np.asarray([p["offset_cartesian_Ainv"] for p in points], dtype=float)
    report, _fit = analyze(q, energies, carrier)
    report.update({"band_index_1based": int(edge["band_index_1based"]),
                   "edge_energy_eV": edge["energy_eV"],
                   "center_fractional": edge["center_fractional"],
                   "metadata": str(md_path), "energies": str(csv_path)})
    reports[carrier] = report

result = {"material": mat, "soc": True, "status": "done",
          "input_files": input_files, "grid": meta["grid"],
          "step_cartesian_Ainv": meta["step_cartesian_Ainv"],
          "band_tracking": meta["band_tracking"], "reports": reports}
(out / "mass_tensor.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
    encoding="utf-8",
)
for carrier, report in reports.items():
    print("%s principal masses (m0): %s" % (carrier, report["principal_masses_me"]))
    for warning in report["warnings"]:
        print("WARNING %s: %s" % (carrier, warning))
print("[DONE] mass_tensor.json written")
