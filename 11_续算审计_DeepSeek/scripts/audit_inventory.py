#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Inventory and hash audit for the four A2B2Te5 SOC dense materials.

Writes:
  csv/inventory.csv
  json/source_manifest.json
  json/inventory_summary.json
No original file is modified.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

import audit_common as ac

KEY_INCAR = ["ENCUT", "NBANDS", "LSORBIT", "LNONCOLLINEAR", "ISYM", "ICHARG",
             "EDIFF", "NELM", "ISMEAR", "SIGMA", "LORBIT", "LWAVE", "LCHARG",
             "KPAR", "NCORE", "IVDW", "GGA", "PREC", "LREAL", "LASPH"]


def parse_incar(path: Path):
    tags = {}
    for line in path.read_text(errors="replace").splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        tags[k.strip().upper()] = v.strip()
    return {k: tags.get(k) for k in KEY_INCAR if k in tags}


def outcar_summary(path: Path):
    text = path.read_text(errors="replace")
    nelect = None
    m = re.search(r"NELECT\s*=\s*([0-9.]+)", text)
    if m:
        nelect = float(m.group(1))
    done = "General timing and accounting informations for this job" in text
    elapsed = None
    m = re.search(r"Elapsed time \(sec\):\s*([0-9.]+)", text)
    if m:
        elapsed = float(m.group(1))
    cpu = None
    m = re.search(r"Total CPU time used \(sec\):\s*([0-9.]+)", text)
    if m:
        cpu = float(m.group(1))
    return {"nelect": nelect, "completed": done, "elapsed_s": elapsed, "cpu_s": cpu}


def history_summary(path: Path):
    rows = []
    if not path.is_file():
        return rows
    for line in path.read_text().splitlines():
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return rows


def main():
    inv_rows = []
    manifest = {"generated_by": "11_续算审计_DeepSeek/scripts/audit_inventory.py",
                "root": str(ac.WORK), "hash_algorithm": "sha256",
                "files": {}, "notes": []}

    def add_hash(label: str, path: Path):
        if not path.is_file():
            manifest["files"][label] = {"path": str(path), "exists": False}
            return None
        h = ac.sha256(path)
        manifest["files"][label] = {"path": str(path), "exists": True,
                                    "size_bytes": path.stat().st_size, "sha256": h}
        return h

    # Project-level templates and settings
    for rel in [
        "A2B2Te5_SOC_dense/project_setting/templates/step3_uniform/gen_dense_soc.py",
        "A2B2Te5_SOC_dense/project_setting/templates/step8.1_boltztrap/run_bt2_soc.py",
        "A2B2Te5_SOC_dense/project_setting/templates/step8.1_boltztrap/gen_step11_boltztrap.py",
        "A2B2Te5_SOC_dense/project_setting/templates/step8.1_boltztrap/ke_common.py",
        "A2B2Te5_SOC_dense/project_setting/templates/step8.2_mass/gen_mass_soc.py",
        "A2B2Te5_SOC_dense/project_setting/templates/step8.3_mass_fit/fit_mass_soc.py",
        "A2B2Te5_SOC_dense/project_setting/templates/step8.3_mass_fit/fit_effective_mass.py",
        "A2B2Te5_SOC_dense/project_setting/tf_A2B2Te5_SOC_dense.yaml",
        "A2B2Te5_SOC_dense/reports/A2B2Te5_综合研究记录_20261013.md",
    ]:
        add_hash("project/" + Path(rel).name + "@" + rel, ac.WORK / rel)

    # Historical code versions that were actually referenced by provenance
    for rel in [
        "补算实施方案_20261007/fit_effective_mass.py",
        "补算实施方案_20261007/fit_effective_mass_validation.json",
    ]:
        add_hash("historical/" + Path(rel).name, ac.WORK / rel)

    for mat in ac.MATERIALS:
        res = ac.RESULT / mat / "ke-dft-cpu" / "result"
        s3 = res / "step3_uniform"
        s81 = res / "step8.1_boltztrap"
        s82 = res / "step8.2_mass"
        s83 = res / "step8.3_mass_fit"
        pos = ac.read_poscar(s3 / "POSCAR")
        inc = parse_incar(s3 / "INCAR")
        kp = ac.read_kpoints(s3 / "KPOINTS")
        eigk, eige, eigo, eigh = ac.read_eigenval(s3 / "EIGENVAL")
        out = outcar_summary(s3 / "OUTCAR")
        zvals, titles = ac.read_potcar_zval(s3 / "POTCAR")
        zsum = float(np.dot(zvals, pos["counts"])) if len(zvals) == len(pos["counts"]) else None
        meta = ac.read_json(s82 / "mass_metadata.json")
        mt = ac.read_json(s83 / "mass_tensor.json")
        bt = ac.read_json(s81 / "boltztrap_crta.json")
        hist = history_summary(res / "provenance" / "history.jsonl")
        failed_marker = (s81 / ".autozt_gen_failed").is_file()

        # Verify mass_tensor source hashes
        mt_hash_ok = {}
        for rel, expected in mt.get("input_files", {}).items():
            p = res / rel
            actual = ac.sha256(p) if p.is_file() else None
            mt_hash_ok[rel] = {"expected": expected, "actual": actual,
                               "match": bool(actual == expected)}

        # Provenance generator hashes versus current project templates
        prov83 = ac.read_json(res / "provenance" / "step8.3_mass_fit.json")
        prov81 = ac.read_json(res / "provenance" / "step8.1_boltztrap.json")
        prov3 = ac.read_json(res / "provenance" / "step3_uniform.json")
        cur_fit_eff = (ac.DENSE / "project_setting/templates/step8.3_mass_fit/fit_effective_mass.py")
        cur_fit_soc = (ac.DENSE / "project_setting/templates/step8.3_mass_fit/fit_mass_soc.py")
        cur_bt2 = (ac.DENSE / "project_setting/templates/step8.1_boltztrap/gen_step11_boltztrap.py")
        cur_common = (ac.DENSE / "project_setting/templates/step8.1_boltztrap/ke_common.py")
        cur_run = (ac.DENSE / "project_setting/templates/step8.1_boltztrap/run_bt2_soc.py")
        cur_dense = (ac.DENSE / "project_setting/templates/step3_uniform/gen_dense_soc.py")
        gen_match = {
            "step3_gen_dense_soc": {
                "provenance": prov3.get("generator", {}).get("sha256"),
                "current": ac.sha256(cur_dense)},
            "step8.1_gen_step11_boltztrap": {
                "provenance": prov81.get("inputs", {}).get("gen_step11_boltztrap.py", {}).get("sha256"),
                "current": ac.sha256(cur_bt2)},
            "step8.1_ke_common": {
                "provenance": prov81.get("inputs", {}).get("ke_common.py", {}).get("sha256"),
                "current": ac.sha256(cur_common)},
            "step8.1_run_bt2_soc": {
                "provenance": prov81.get("inputs", {}).get("run_bt2_soc.py", {}).get("sha256"),
                "current": ac.sha256(cur_run)},
            "step8.3_fit_effective_mass": {
                "provenance": prov83.get("inputs", {}).get("fit_effective_mass.py", {}).get("sha256"),
                "current": ac.sha256(cur_fit_eff)},
            "step8.3_fit_mass_soc": {
                "provenance": prov83.get("inputs", {}).get("fit_mass_soc.py", {}).get("sha256"),
                "current": ac.sha256(cur_fit_soc)},
        }
        for k, v in gen_match.items():
            v["match"] = bool(v["provenance"] == v["current"])

        # Hash important local files (excluding multi-GB vasprun here)
        for label, p in [
            ("step3_EIGENVAL", s3 / "EIGENVAL"), ("step3_INCAR", s3 / "INCAR"),
            ("step3_KPOINTS", s3 / "KPOINTS"), ("step3_POSCAR", s3 / "POSCAR"),
            ("step3_POTCAR", s3 / "POTCAR"), ("step3_OUTCAR", s3 / "OUTCAR"),
            ("step8.1_boltztrap_crta", s81 / "boltztrap_crta.json"),
            ("step8.2_EIGENVAL", s82 / "EIGENVAL"),
            ("step8.2_mass_metadata", s82 / "mass_metadata.json"),
            ("step8.2_KPOINTS", s82 / "KPOINTS"),
            ("step8.3_mass_tensor", s83 / "mass_tensor.json"),
            ("step8.3_electron_band", s83 / "electron_band.csv"),
            ("step8.3_hole_band", s83 / "hole_band.csv"),
        ]:
            add_hash(f"{mat}/{label}", p)

        prov_retries_81 = [x for x in hist if x.get("step") == "step8.1_boltztrap"]
        row = {
            "material": mat,
            "poscar_title": pos["title"],
            "formula_species": " ".join(pos["species"]),
            "formula_counts": " ".join(str(x) for x in pos["counts"]),
            "natoms": pos["natoms"],
            "volume_A3": pos["volume_A3"],
            "lattice_a_A": float(np.linalg.norm(pos["lattice_rows_A"][0])),
            "lattice_b_A": float(np.linalg.norm(pos["lattice_rows_A"][1])),
            "lattice_c_A": float(np.linalg.norm(pos["lattice_rows_A"][2])),
            "incar": json.dumps(inc, ensure_ascii=False),
            "kpoints_title": kp["title"],
            "kpoints_mode": kp["mode"],
            "kpoints_mesh": "x".join(str(x) for x in kp["mesh"]) if kp["mesh"] else "",
            "kpoints_nk": kp["nk"],
            "eigenval_nelect": eigh["nelect"],
            "eigenval_nk": eigh["nk"],
            "eigenval_nbands": eigh["nb"],
            "eigenval_ispin": eigh["ispin"],
            "outcar_nelect": out["nelect"],
            "outcar_completed": out["completed"],
            "outcar_elapsed_s": out["elapsed_s"],
            "potcar_zval": ";".join(str(x) for x in zvals),
            "potcar_zval_sum": zsum,
            "nelect_consistent": bool(out["nelect"] is not None and eigh["nelect"] == out["nelect"] and
                                      (zsum is None or abs(zsum - eigh["nelect"]) < 1e-6)),
            "mass_metadata_grid": meta.get("grid"),
            "mass_metadata_step_Ainv": meta.get("step_cartesian_Ainv"),
            "mass_metadata_n_electron": len(meta.get("points_by_edge", {}).get("electron", [])),
            "mass_metadata_n_hole": len(meta.get("points_by_edge", {}).get("hole", [])),
            "mass_metadata_hole_k": json.dumps(meta["edges"]["hole"]["center_fractional"]),
            "mass_metadata_hole_band": meta["edges"]["hole"]["band_index_1based"],
            "mass_metadata_hole_E": meta["edges"]["hole"]["energy_eV"],
            "mass_metadata_electron_k": json.dumps(meta["edges"]["electron"]["center_fractional"]),
            "mass_metadata_electron_band": meta["edges"]["electron"]["band_index_1based"],
            "mass_metadata_electron_E": meta["edges"]["electron"]["energy_eV"],
            "mass_tensor_status": mt.get("status"),
            "mass_tensor_inputs_match": all(v["match"] for v in mt_hash_ok.values()),
            "mass_tensor_inputs": json.dumps(mt_hash_ok, ensure_ascii=False),
            "bt2_target_gap_eV": bt.get("target_gap_eV"),
            "bt2_dft_gap_eV": bt.get("scissor", {}).get("dft_gap_eV"),
            "bt2_scissor_shift_eV": bt.get("scissor", {}).get("shift_eV"),
            "bt2_scissor_applied": bt.get("scissor", {}).get("applied"),
            "bt2_nelect": bt.get("nelect"),
            "bt2_cell_volume_cm3": bt.get("cell_volume_cm3"),
            "bt2_temperatures": json.dumps(bt.get("temperatures_K")),
            "bt2_mu_npts": len(bt.get("mu_rel_efermi_eV", [])),
            "bt2_mu_min": min(bt.get("mu_rel_efermi_eV", [None])),
            "bt2_mu_max": max(bt.get("mu_rel_efermi_eV", [None])),
            "bt2_gen_failed_marker": failed_marker,
            "bt2_retry_records": len(prov_retries_81),
            "provenance_81_ts": prov81.get("ts"),
            "provenance_83_ts": prov83.get("ts"),
            "gen_hash_match": json.dumps({k: v["match"] for k, v in gen_match.items()}),
            "gen_hash_detail": json.dumps(gen_match, ensure_ascii=False),
        }
        inv_rows.append(row)

    ac.write_csv(ac.CSV / "inventory.csv", inv_rows)
    manifest["notes"].append(
        "SHA256 of step3_uniform/vasprun.xml (multi-GB) is recorded in logs/step3_vasprun_sha256.txt "
        "and is not re-hashed in this JSON to keep the inventory fast.")
    manifest["inventory_rows"] = len(inv_rows)
    ac.write_json(ac.JSON / "source_manifest.json", manifest)
    ac.write_json(ac.JSON / "inventory_summary.json", {
        "materials": ac.MATERIALS,
        "rows": inv_rows,
        "all_mass_tensor_input_hashes_match": all(r["mass_tensor_inputs_match"] for r in inv_rows),
        "all_nelect_consistent": all(r["nelect_consistent"] for r in inv_rows),
        "current_template_matches_provenance": {
            r["material"]: json.loads(r["gen_hash_match"]) for r in inv_rows},
    })
    print("inventory rows:", len(inv_rows))
    for r in inv_rows:
        print(r["material"], "natoms", r["natoms"], "mesh", r["kpoints_mesh"],
              "nelect", r["eigenval_nelect"], "mass_hash_ok", r["mass_tensor_inputs_match"],
              "gen_match", r["gen_hash_match"])


if __name__ == "__main__":
    main()
