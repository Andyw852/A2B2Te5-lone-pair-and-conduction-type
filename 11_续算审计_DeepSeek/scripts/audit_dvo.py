#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit whether D/V/O (disordered / intermediate / ordered rocksalt models)
have any actual energies, and separate them from the old 81-atom defect runs.

Hard evidence:
  * 补算实施方案_20261007/rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED has geometry
    and bookkeeping only (243 jobs, all NOT_RUN; validation says 0 DFT);
  * 交付核验.json says new_DFT_executed=false and
    actual_new_material_masses_or_defect_energies_available=false;
  * no OUTCAR/vasprun/OSZICAR/CONTCAR/.h5 files exist in that package;
  * the four /mnt/d/tf_data/<material>/defect-dft-cpu projects are 81-atom
    3x3x1 cells with step4 status=partial and are NOT D/V/O;
  * 03_缺陷形成能与PN型数据/*.zip is an even older 72-atom 2x2x2 dataset
    with formation_energy_results.json, also NOT D/V/O.

Writes csv/dvo_status.csv, csv/dvo_old_defects.csv,
json/dvo_audit_summary.json.
"""
from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd

import audit_common as ac

PLAN = ac.WORK / "补算实施方案_20261007"
if not PLAN.is_dir():
    PLAN = ac.WORK.parent / "补算实施方案_20261007"
DVO = PLAN / "rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED"
OLD_TF = Path("/mnt/d/tf_data")
OLD_ZIP = ac.WORK / "03_缺陷形成能与PN型数据"


def read_zip_poscar_counts(zpath: Path):
    with zipfile.ZipFile(zpath) as z:
        names = [n for n in z.namelist() if n.endswith("bulk/POSCAR") or n.endswith("/POSCAR")]
        if not names:
            return None
        text = z.read(names[0]).decode("utf-8", errors="replace").splitlines()
        species = text[5].split()
        counts = [int(x) for x in text[6].split()]
        return {"member": names[0], "species": species, "counts": counts,
                "natoms": sum(counts)}


def read_zip_formation(zpath: Path):
    with zipfile.ZipFile(zpath) as z:
        names = [n for n in z.namelist() if n.endswith("formation_energy_results.json")]
        if not names:
            return None
        d = json.loads(z.read(names[0]).decode("utf-8", errors="replace"))
        return {"member": names[0], "status": d.get("status"),
                "E_gap": d.get("E_gap"), "E_VBM": d.get("E_VBM"),
                "dominant_defect": d.get("dominant_defect"),
                "carrier_type": d.get("carrier_type"),
                "n_defects": len(d.get("defects", {}))}


def main():
    jobs = pd.read_csv(DVO / "jobs.csv", encoding="utf-8-sig")
    manifest = json.loads((DVO / "manifest.json").read_text())
    validation = json.loads((DVO / "validation.json").read_text())
    delivery = json.loads((PLAN / "交付核验.json").read_text())
    independent = json.loads((DVO / "independent_audit.json").read_text())
    review = pd.read_csv(PLAN / "复核_所需数据与可得结果.csv", encoding="utf-8-sig")

    all_files = [p for p in DVO.rglob("*") if p.is_file()]
    energy_like = [p for p in all_files if any(s in p.name.upper() for s in
                                               ["OUTCAR", "VASPRUN", "OSZICAR", "CONTCAR", ".H5", "ENERGY"])]
    extensions = {}
    for p in all_files:
        ext = p.suffix.lower() or "<none>"
        extensions[ext] = extensions.get(ext, 0) + 1

    status_counts = jobs["status"].value_counts().to_dict()
    stage_counts = jobs["stage"].value_counts().to_dict()
    parent_counts = jobs["parent_id"].value_counts().to_dict()
    charge_counts = jobs["intended_charge_state"].value_counts().to_dict()
    defect_counts = jobs["defect_type"].value_counts().to_dict()

    # old 81-atom AutoZT defect projects
    old_rows = []
    for mat in ac.MATERIALS:
        base = OLD_TF / mat / "defect-dft-cpu" / "result"
        pos = base / "step1_bulk/POSCAR"
        counts = None
        species = None
        if pos.is_file():
            lines = pos.read_text(errors="replace").splitlines()
            species = lines[5].split()
            counts = [int(x) for x in lines[6].split()]
        summary = base / "step4_analysis/formation_energy_summary.json"
        status = None
        n_def = None
        if summary.is_file():
            d = json.loads(summary.read_text())
            status = d.get("status")
            n_def = len(d.get("defects", {})) if isinstance(d.get("defects"), dict) else None
        old_rows.append({
            "material": mat,
            "path": str(base),
            "bulk_species": " ".join(species or []),
            "bulk_counts": " ".join(str(x) for x in (counts or [])),
            "bulk_natoms": sum(counts) if counts else None,
            "supercell": "3x3x1 (81 atoms)" if counts and sum(counts) == 81 else "unknown",
            "step4_status": status,
            "step4_n_defects": n_def,
            "is_DVO": False,
            "note": "old 81-atom AutoZT defect project; not D/V/O; partial step4; do not use as D/V/O energy",
        })

    # old 72-atom zips
    zip_rows = []
    for zpath in sorted(OLD_ZIP.glob("*.zip")):
        pos = read_zip_poscar_counts(zpath)
        form = read_zip_formation(zpath)
        zip_rows.append({
            "material": zpath.stem,
            "zip": str(zpath),
            "bulk_member": pos["member"] if pos else None,
            "bulk_species": " ".join(pos["species"]) if pos else None,
            "bulk_counts": " ".join(str(x) for x in pos["counts"]) if pos else None,
            "bulk_natoms": pos["natoms"] if pos else None,
            "supercell": "2x2x2 (72 atoms)" if pos and pos["natoms"] == 72 else "unknown",
            "formation_status": form.get("status") if form else None,
            "E_gap": form.get("E_gap") if form else None,
            "dominant_defect": form.get("dominant_defect") if form else None,
            "carrier_type": form.get("carrier_type") if form else None,
            "n_defects": form.get("n_defects") if form else None,
            "is_DVO": False,
            "note": "old 2x2x2 72-atom dataset; has old formation energies but is not the 144-atom D/V/O model",
        })

    # Top-level D/V/O status rows
    status_rows = [
        {"item": "DVO_parent_geometries", "count": int(len(list((DVO / "parents").glob("*/POSCAR")))),
         "status": "geometry_only", "evidence": "parents/*/POSCAR exist; no parent CONTCAR/OUTCAR",
         "completed": True, "has_energy": False},
        {"item": "DVO_defect_geometries", "count": int(len(list((DVO / "defects").glob("*/*/POSCAR")))),
         "status": "geometry_only", "evidence": "defects/*/*/POSCAR + defect.json; no OUTCAR/vasprun",
         "completed": True, "has_energy": False},
        {"item": "DVO_jobs", "count": int(len(jobs)),
         "status": json.dumps(status_counts), "evidence": "jobs.csv; all rows NOT_RUN",
         "completed": False, "has_energy": False},
        {"item": "DVO_DFT_calculations", "count": 0,
         "status": "0", "evidence": "validation.json df_t_calculations_performed=0; 交付核验 new_DFT_executed=false",
         "completed": False, "has_energy": False},
        {"item": "DVO_energy_files", "count": int(len(energy_like)),
         "status": "none", "evidence": "no OUTCAR/vasprun/OSZICAR/CONTCAR/.h5/energy files anywhere in package",
         "completed": False, "has_energy": False},
        {"item": "DVO_actual_material_energies", "count": 0,
         "status": "absent",
         "evidence": "交付核验 actual_new_material_masses_or_defect_energies_available=false",
         "completed": False, "has_energy": False},
        {"item": "DVO_validation", "count": int(validation.get("serialized_poscar_files_checked", 0)),
         "status": validation.get("result"),
         "evidence": "geometry/bookkeeping only; df_t_calculations_performed=0",
         "completed": True, "has_energy": False},
        {"item": "DVO_independent_audit", "count": int(independent.get("all_serialized_poscars_verified_against_parent_and_mutation", 0)),
         "status": independent.get("result"),
         "evidence": independent.get("scientific_scope"),
         "completed": True, "has_energy": False},
    ]

    ac.write_csv(ac.CSV / "dvo_status.csv", status_rows)
    ac.write_csv(ac.CSV / "dvo_old_81atom.csv", old_rows)
    ac.write_csv(ac.CSV / "dvo_old_72atom_zips.csv", zip_rows)
    old_all = old_rows + zip_rows
    old_fields = sorted({k for row in old_all for k in row})
    ac.write_csv(ac.CSV / "dvo_old_defects.csv", old_all, fieldnames=old_fields)
    summary = {
        "generated_by": "audit_dvo.py",
        "dvo_root": str(DVO),
        "package_status": manifest.get("status"),
        "parents": manifest.get("counts", {}).get("parents"),
        "defects": manifest.get("counts", {}).get("defects"),
        "jobs": manifest.get("counts", {}).get("jobs"),
        "jobs_status_counts": status_counts,
        "jobs_stage_counts": stage_counts,
        "jobs_parent_counts": parent_counts,
        "jobs_charge_state_counts": charge_counts,
        "jobs_defect_type_counts": defect_counts,
        "file_extensions": extensions,
        "n_files": len(all_files),
        "energy_like_files": [str(p) for p in energy_like],
        "validation": validation,
        "delivery": {k: delivery.get(k) for k in
                     ["date", "new_DFT_executed",
                      "actual_new_material_masses_or_defect_energies_available",
                      "completion_scope", "review_scope"]},
        "independent_audit_scope": independent.get("scientific_scope"),
        "review_table_rows": len(review),
        "old_81atom_projects": old_rows,
        "old_72atom_zips": zip_rows,
        "conclusion": {
            "DVO_has_actual_energies": False,
            "reason": "243 geometry/job entries are all NOT_RUN; validation and delivery both state zero DFT; no energy/output files exist.",
            "old_81atom_is_DVO": False,
            "old_72atom_is_DVO": False,
        },
        "minimal_priority_list": [
            {"priority": 1, "status": "missing", "task": "Relax the 9 Pb2Sb2Te5 D/V/O parent cells (PBE-D3, consistent ENCUT/k-mesh; SOC decision recorded) and verify they remain the intended rocksalt/ordered occupation."},
            {"priority": 2, "status": "missing", "task": "Regenerate the 12 neutral defect types from each relaxed parent (not the unrelaxed POSCARs), run static energies, and compute neutral formation energies with a common chemical-potential convention."},
            {"priority": 3, "status": "missing", "task": "For the competitive defects (Sb_Pb, Pb_Sb, Pb_on_Vac, Sb_on_Vac, V_Te, V_Pb, V_Sb), add relevant charge states across the gap, finite-size corrections, and a charge-neutrality solution."},
            {"priority": 4, "status": "missing", "task": "Repeat the same D/V/O protocol for Pb2Bi2Te5 to test the Bi/Sb contrast; do not infer the contrast from the old 81-atom project."},
            {"priority": 5, "status": "missing", "task": "Only after 1-4: connect defect-limited n/p and mobility to transport; static 0 K formation energies alone do not give an annealing temperature or a P/N switch."},
        ],
        "completed_but_not_energies": [
            "9 parent geometries and 234 single-site defect geometries",
            "243-row job table, all NOT_RUN",
            "geometry/bookkeeping validation PASS; independent mutation audit PASS",
            "6 synthetic tests of the effective-mass fitting tool (not real material masses)",
        ],
        "needs_fix": [
            "All 9 parent cells are unrelaxed (a=6.2 A artificial starting guess).",
            "All 234 defect cells are built from the unrelaxed parent and must be regenerated from relaxed CONTCAR.",
            "No POTCAR/INCAR/KPOINTS/submit script in the D/V/O package.",
            "O-site equivalence after relaxation must be re-checked; D/V random cells are not SQS.",
            "The old 81-atom step4 summaries are status=partial; the old 72-atom zips are a different cell and dataset.",
        ],
    }
    ac.write_json(ac.JSON / "dvo_audit_summary.json", summary)
    print("DVO parents/defects/jobs:", summary["parents"], summary["defects"], summary["jobs"])
    print("jobs status:", status_counts)
    print("energy-like files:", len(energy_like))
    print("old 81-atom:", [(r["material"], r["bulk_natoms"], r["step4_status"]) for r in old_rows])
    print("old 72-atom:", [(r["material"], r["bulk_natoms"], r["formation_status"]) for r in zip_rows])
    print("wrote", ac.CSV / "dvo_status.csv")


if __name__ == "__main__":
    main()
