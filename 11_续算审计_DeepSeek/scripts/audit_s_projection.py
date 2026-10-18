#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit the old HSE06+SOC s-projection / 50 meV band-edge window data.

Questions answered with file evidence:
  * Is the 50 meV VBM window below the VBM (as it must be) or above it?
  * Which structure/functional/k-sampling produced the s projection?
  * Is that band edge the same as the PBE+SOC dense-grid edge used by the new
    BoltzTraP2 transport?
  * Do the A-site s percentages in the 2026-10-13 report equal the CSV values?
  * Is the s percentage a normalized orbital projection rather than a count?

Writes csv/s_projection_audit.csv, csv/s_projection_source_check.csv,
json/s_projection_summary.json.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import audit_common as ac

BANDEDGE = ac.WORK / "07_带边能窗"
REPORT_A_S = {"Pb2Sb2Te5": 7.34, "Pb2Bi2Te5": 2.56,
              "Sn2Sb2Te5": 3.63, "Sn2Bi2Te5": 23.06}
ELEM_COLS = ["Sb-s", "Sb-p", "Sb-d", "Te-s", "Te-p", "Te-d",
             "Pb-s", "Pb-p", "Pb-d", "Bi-s", "Bi-p", "Bi-d",
             "Sn-s", "Sn-p", "Sn-d"]


def gzip_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_incar(path: Path):
    tags = {}
    for line in path.read_text(errors="replace").splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            tags[k.strip().upper()] = v.strip()
    return {k: tags.get(k) for k in ["SYSTEM", "LHFCALC", "HFSCREEN", "AEXX",
                                     "LSORBIT", "GGA", "ENCUT", "NBANDS", "ISYM",
                                     "ISPIN", "EDIFF", "SIGMA"] if k in tags}


def read_poscar_full(path: Path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    scale = float(lines[1].split()[0])
    lattice = np.asarray([[float(x) for x in lines[j].split()[:3]] for j in (2, 3, 4)])
    if scale < 0:
        scale = (-scale / abs(np.linalg.det(lattice))) ** (1.0 / 3.0)
    lattice = lattice * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    mode = lines[7].strip().lower()
    start = 8
    if mode.startswith("s"):
        start = 9
    coords = np.asarray([[float(x) for x in lines[start + i].split()[:3]] for i in range(n)])
    if mode.startswith("c") or mode.startswith("k"):
        coords = coords @ np.linalg.inv(lattice)
    return {"lattice": lattice, "species": species, "counts": counts,
            "coords_frac": coords, "volume": abs(float(np.linalg.det(lattice)))}


def main():
    states = pd.read_csv(BANDEDGE / "入选态清单.csv", encoding="utf-8-sig")
    stats = pd.read_csv(BANDEDGE / "带边能窗统计.csv", encoding="utf-8-sig")
    mass_summary = ac.read_json(ac.JSON / "mass_audit_summary.json")
    audit_rows = []
    source_rows = []
    summary = {"generated_by": "audit_s_projection.py", "materials": {}}

    for mat in ac.MATERIALS:
        sub = states[states.material == mat]
        vbm = float(sub[sub.edge == "VBM"].energy_eV.max())
        cbm = float(sub[sub.edge == "CBM"].energy_eV.min())
        hse_gap = cbm - vbm
        edge_k = {}
        for edge, val, fn in (("VBM", vbm, "max"), ("CBM", cbm, "min")):
            d = sub[sub.edge == edge]
            row = d.loc[d.energy_eV.idxmax() if fn == "max" else d.energy_eV.idxmin()]
            edge_k[edge] = [float(row.kx), float(row.ky), float(row.kz)]

        src = json.loads((BANDEDGE / "源数据" / mat / "原始来源.json").read_text())
        incar = parse_incar(BANDEDGE / "源数据" / mat / "INCAR")
        kp_lines = (BANDEDGE / "源数据" / mat / "KPOINTS").read_text(errors="replace").splitlines()
        kp_title = kp_lines[0].strip()
        kp_n = int(kp_lines[1].split()[0])
        proc = sorted((BANDEDGE / "源数据" / mat).glob("PROCAR*.gz"))
        proc_hash = gzip_sha256(proc[0]) if proc else None
        hash_match = bool(proc_hash == src.get("sha256_uncompressed"))

        pos_hse = read_poscar_full(BANDEDGE / "源数据" / mat / "POSCAR")
        pos_dense = read_poscar_full(ac.RESULT / mat / "ke-dft-cpu/result/step3_uniform/POSCAR")
        same_species = pos_hse["species"] == pos_dense["species"] and pos_hse["counts"] == pos_dense["counts"]
        lat_diff = float(np.max(np.abs(pos_hse["lattice"] - pos_dense["lattice"])))
        coord_diff = float(np.max(np.abs(pos_hse["coords_frac"] - pos_dense["coords_frac"]))) if same_species else None
        poscar_same = bool((BANDEDGE / "源数据" / mat / "POSCAR").read_bytes() ==
                           (ac.RESULT / mat / "ke-dft-cpu/result/step3_uniform/POSCAR").read_bytes())

        pbe_hole = mass_summary["materials"][mat]["hole"]
        pbe_elec = mass_summary["materials"][mat]["electron"]
        pbe_edge = {"VBM": pbe_hole["center_fractional"], "CBM": pbe_elec["center_fractional"]}
        pbe_E = {"VBM": pbe_hole["metadata_edge_energy_eV"], "CBM": pbe_elec["metadata_edge_energy_eV"]}
        edge_shift = {}
        for edge in ("VBM", "CBM"):
            dk = np.asarray(edge_k[edge]) - np.asarray(pbe_edge[edge])
            dk = dk - np.round(dk)
            edge_shift[edge] = {"hse_k": edge_k[edge], "pbe_k": pbe_edge[edge],
                                "wrapped_fractional_distance": float(np.linalg.norm(dk)),
                                "delta_k": dk.tolist(),
                                "hse_E_eV": vbm if edge == "VBM" else cbm,
                                "pbe_E_eV": pbe_E[edge]}

        # 50 meV state window checks
        for edge in ("VBM", "CBM"):
            d = sub[(sub.edge == edge) & (np.isclose(sub.window_meV, 50.0))]
            st = stats[(stats.material == mat) & (stats.edge == edge) & (np.isclose(stats.window_meV, 50.0))]
            if d.empty or st.empty:
                continue
            a_col = ac.A_SITE[mat] + "-s"
            b_col = ac.B_SITE[mat] + "-s"
            row = st.iloc[0]
            if edge == "VBM":
                below = bool(np.all(d.energy_eV <= vbm + 1e-9))
                above = False
                max_dist = float(np.max(vbm - d.energy_eV)) * 1000.0
                min_dist = float(np.min(vbm - d.energy_eV)) * 1000.0
            else:
                below = False
                above = bool(np.all(d.energy_eV >= cbm - 1e-9))
                max_dist = float(np.max(d.energy_eV - cbm)) * 1000.0
                min_dist = float(np.min(d.energy_eV - cbm)) * 1000.0
            elem_sum = float(sum(float(row[c]) for c in ELEM_COLS if c in row and pd.notna(row[c])))
            a_s = float(row[a_col])
            b_s = float(row[b_col])
            te_s = float(row["Te-s"])
            total_s = float(sum(float(row[c]) for c in ELEM_COLS if c.endswith("-s") and c in row and pd.notna(row[c])))
            # The 2026-10-13 report table is a VBM-side A-site s value; do not
            # compare it to the CBM-side rows.
            report_a = REPORT_A_S.get(mat) if edge == "VBM" else None
            audit_rows.append({
                "material": mat, "edge": edge, "window_meV": 50.0,
                "n_states": int(row["states"]), "n_unique_k": int(row["unique_k"]),
                "a_site": ac.A_SITE[mat], "b_site": ac.B_SITE[mat],
                "A_site_s_pct": a_s, "B_site_s_pct": b_s, "Te_s_pct": te_s,
                "total_s_pct": total_s, "raw_sum_mean": float(row["raw_sum_mean"]),
                "sum_all_element_orbitals_pct": elem_sum,
                "window_is_below_VBM": below, "window_is_above_CBM": above,
                "min_distance_from_edge_meV": min_dist,
                "max_distance_from_edge_meV": max_dist,
                "state_energy_min_eV": float(d.energy_eV.min()),
                "state_energy_max_eV": float(d.energy_eV.max()),
                "HSE_VBM_eV": vbm, "HSE_CBM_eV": cbm, "HSE_gap_eV": hse_gap,
                "PBE_dense_VBM_eV": pbe_E["VBM"], "PBE_dense_CBM_eV": pbe_E["CBM"],
                "PBE_dense_gap_eV": pbe_E["CBM"] - pbe_E["VBM"],
                "HSE_VBM_k": json.dumps(edge_k["VBM"]), "HSE_CBM_k": json.dumps(edge_k["CBM"]),
                "PBE_VBM_k": json.dumps(pbe_edge["VBM"]), "PBE_CBM_k": json.dumps(pbe_edge["CBM"]),
                "VBM_k_shift": json.dumps(edge_shift["VBM"]),
                "CBM_k_shift": json.dumps(edge_shift["CBM"]),
                "report_A_site_s_pct": report_a,
                "report_minus_computed_A_site_s_pct": (None if report_a is None else report_a - a_s),
                "source_functional": "HSE06+SOC (LHFCALC, AEXX=0.25, HFSCREEN=0.2)",
                "source_kpoints": f"{kp_title} (nk={kp_n})",
                "source_poscar_same_as_dense": poscar_same,
                "structure_lattice_max_abs_diff_A": lat_diff,
                "structure_coord_max_abs_diff_frac": coord_diff,
                "source_procar": str(proc[0]) if proc else "",
                "source_procar_uncompressed_sha256_match": hash_match,
                "report_phrase_VBM_above_50meV_is_correct": (not below) if edge == "VBM" else None,
            })

        source_rows.append({
            "material": mat,
            "source_json": json.dumps(src, ensure_ascii=False),
            "source_procar": str(proc[0]) if proc else "",
            "source_procar_sha256_uncompressed": proc_hash,
            "source_procar_hash_match": hash_match,
            "incar": json.dumps(incar, ensure_ascii=False),
            "functional": "HSE06+SOC" if incar.get("LHFCALC", "").upper() == ".TRUE." else "unknown",
            "kpoints_title": kp_title, "kpoints_nk": kp_n,
            "kpoints_is_dense_28x28xN": False,
            "kpoints_type": "HSE band path (38 uniform SCF k-points + path), not the 28x28xN dense transport grid",
            "HSE_VBM_eV": vbm, "HSE_CBM_eV": cbm, "HSE_gap_eV": hse_gap,
            "PBE_dense_VBM_eV": pbe_E["VBM"], "PBE_dense_CBM_eV": pbe_E["CBM"],
            "PBE_dense_gap_eV": pbe_E["CBM"] - pbe_E["VBM"],
            "HSE_VBM_k": json.dumps(edge_k["VBM"]), "HSE_CBM_k": json.dumps(edge_k["CBM"]),
            "PBE_VBM_k": json.dumps(pbe_edge["VBM"]), "PBE_CBM_k": json.dumps(pbe_edge["CBM"]),
            "VBM_k_shift": json.dumps(edge_shift["VBM"]),
            "CBM_k_shift": json.dumps(edge_shift["CBM"]),
            "absolute_eigenvalue_shift_note": "HSE and PBE absolute eigenvalues use different Hamiltonians/pseudopotentials and are not comparable; only gaps and k-point shifts are used.",
            "poscar_same_as_dense": poscar_same,
            "structure_lattice_max_abs_diff_A": lat_diff,
            "structure_coord_max_abs_diff_frac": coord_diff,
            "species_and_counts_same": same_species,
        })
        summary["materials"][mat] = {
            "HSE_VBM_eV": vbm, "HSE_CBM_eV": cbm, "HSE_gap_eV": hse_gap,
            "HSE_VBM_k": edge_k["VBM"], "HSE_CBM_k": edge_k["CBM"],
            "PBE_dense_VBM_eV": pbe_E["VBM"], "PBE_dense_CBM_eV": pbe_E["CBM"],
            "PBE_dense_gap_eV": pbe_E["CBM"] - pbe_E["VBM"],
            "PBE_dense_VBM_k": pbe_edge["VBM"], "PBE_dense_CBM_k": pbe_edge["CBM"],
            "edge_shift": edge_shift,
            "poscar_same_as_dense": poscar_same,
            "structure_lattice_max_abs_diff_A": lat_diff,
            "source_procar_hash_match": hash_match,
        }
        print(f"[{mat}] HSE gap {hse_gap:.6f} eV; PBE dense gap "
              f"{pbe_E['CBM']-pbe_E['VBM']:.6f} eV; POSCAR same={poscar_same}")

    ac.write_csv(ac.CSV / "s_projection_audit.csv", audit_rows)
    ac.write_csv(ac.CSV / "s_projection_source_check.csv", source_rows)
    ac.write_json(ac.JSON / "s_projection_summary.json", summary)
    print("wrote", ac.CSV / "s_projection_audit.csv")
    print("wrote", ac.CSV / "s_projection_source_check.csv")


if __name__ == "__main__":
    main()
