#!/usr/bin/env python3
"""Independent serialized-model audit plus synthetic software mapping checks.

Synthetic inputs are injected in memory only. No fabricated relaxed CONTCAR or
mapped scientific results are saved. The audit result is bookkeeping evidence.
"""

import argparse
import csv
import importlib.util
import itertools
import json
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np


def parse_text(text):
    lines = text.splitlines()
    lattice = np.array([[float(x) for x in line.split()] for line in lines[2:5]])
    symbols = [s for s, n in zip(lines[5].split(), map(int, lines[6].split())) for _ in range(n)]
    frac = np.array([[float(x) for x in line.split()] for line in lines[8:]])
    return lattice, frac, symbols


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", type=Path, required=True)
    args = parser.parse_args()
    root = args.models.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    lattice = np.array(manifest["lattice_vectors_A_rows"])
    a = manifest["initial_cubic_a_A"]
    with (root / "jobs.csv").open(encoding="utf-8-sig", newline="") as handle:
        jobs = list(csv.DictReader(handle))
    assert len(jobs) == manifest["counts"]["jobs"]
    assert all(j["status"] == "NOT_RUN" for j in jobs)
    assert len({j["job_id"] for j in jobs}) == len(jobs)
    expected_types = {"Sb_Pb", "Pb_Sb", "Pb_on_Vac", "Sb_on_Vac", "V_Te", "V_Pb", "V_Sb",
                      "Te_on_Vac", "Pb_Te", "Sb_Te", "Te_Pb", "Te_Sb"}
    if manifest["material"] == "Pb2Bi2Te5":
        expected_types = {s.replace("Sb", "Bi") for s in expected_types}
    audit_count = 0
    vacancies = {}
    for parent in manifest["parents"]:
        pid = parent["parent_id"]
        site_map = json.loads((root / "parents" / pid / "site_map.json").read_text(encoding="utf-8"))
        sites = site_map["sites"]
        assert len({s["site_id"] for s in sites}) == len(sites)
        frac = np.array([s["fractional_coordinates"] for s in sites])
        grid_coords = frac @ lattice / (a / 2)
        assert np.max(np.abs(grid_coords - np.rint(grid_coords))) < 1e-10
        assert len(np.unique(np.round(frac, 10), axis=0)) == len(sites)
        vacancies[pid] = [sum(s["occupancy"] == "Vac" and s["base_plane_index"] == q for s in sites)
                          for q in range(5)]
        if pid.startswith(("O", "V")):
            assert vacancies[pid][:-1] == [0] * 4
        if pid.startswith("O"):
            transformed = np.mod(np.array([3 * a, 0, 0]) @ np.linalg.inv(lattice) - frac, 1)
            for i, x in enumerate(transformed):
                diff = x - frac
                diff -= np.rint(diff)
                j = np.argmin(np.sum(diff ** 2, axis=1))
                assert np.linalg.norm(diff[j]) < 1e-10
                assert sites[i]["occupancy"] == sites[j]["occupancy"]
        parent_jobs = [j for j in jobs if j["parent_id"] == pid]
        assert {j["defect_type"] for j in parent_jobs if j["defect_type"]} == expected_types
        for job in parent_jobs:
            path = root / job["poscar"]
            lat, actual_frac, actual_symbols = parse_text(path.read_text(encoding="ascii"))
            assert np.allclose(lat, lattice, atol=1e-12)
            symbols = [s["occupancy"] for s in sites]
            if job["defect_type"]:
                mutation = json.loads((path.parent / "defect.json").read_text(encoding="utf-8"))
                indices = [i for i, s in enumerate(sites) if s["site_id"] == mutation["changed_site"]["site_id"]]
                assert len(indices) == 1
                i = indices[0]
                assert symbols[i] == mutation["source_occupancy"]
                symbols[i] = mutation["target_occupancy"]
            expected = [(sym, s["fractional_coordinates"]) for sym, s in zip(symbols, sites) if sym != "Vac"]
            element_order = path.read_text(encoding="ascii").splitlines()[5].split()
            expected = [item for element in element_order for item in expected if item[0] == element]
            assert actual_symbols == [e[0] for e in expected]
            assert len(actual_symbols) == int(job["atom_count"])
            assert np.allclose(actual_frac, np.array([e[1] for e in expected]), atol=1e-12)
            audit_count += 1
    # Independently test mapping functionality; in-memory synthetic perturbations
    # are not submitted or saved as relaxed-material evidence.
    utility_path = Path(__file__).with_name("regenerate_defects_from_relaxed_parent.py")
    spec = importlib.util.spec_from_file_location("mapping_utility", utility_path)
    mapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mapper)
    pid = "O_ordered"
    opts = SimpleNamespace(site_map=root / "parents" / pid / "site_map.json",
                           relaxed_poscar=root / "parents" / pid / "POSCAR",
                           templates=root / "defects" / pid,
                           ghost_method="neighbor", max_displacement=1.5, min_distance=1.5)
    real_read = mapper.read_poscar
    original = real_read(opts.relaxed_poscar)
    original_lat, original_frac, species, counts, selective = original
    checks = {}
    for scenario, scale, shift in [("identity", 1.0, np.zeros(3)),
                                   ("affine_cell_change", 1.02, np.zeros(3)),
                                   ("periodic_boundary_translation", 1.0, np.array([0.06, 0, 0]))]:
        if scenario == "periodic_boundary_translation":
            assert np.any(original_frac + shift >= 1.0)
        mapper.read_poscar = lambda path: (original_lat * scale, np.mod(original_frac + shift, 1),
                                            species, counts, selective)
        plans, diagnostics, _, _ = mapper.prepare(opts)
        for defect_id, text, info in plans:
            actual_lat, actual_frac, actual_symbols = parse_text(text)
            reference = (root / "defects" / pid / defect_id / "POSCAR").read_text(encoding="ascii")
            ref_lat, ref_frac, ref_symbols = parse_text(reference)
            assert np.allclose(actual_lat, ref_lat * scale, atol=1e-11)
            assert Counter(actual_symbols) == Counter(ref_symbols)
            # Insertion is appended before element regrouping, so compare by
            # matching each atom of a given species modulo periodic boundaries.
            for symbol in set(actual_symbols):
                out_xyz = actual_frac[np.array(actual_symbols) == symbol]
                ref_xyz = np.mod(ref_frac[np.array(ref_symbols) == symbol] + shift, 1)
                differences = out_xyz[:, None, :] - ref_xyz[None, :, :]
                differences -= np.rint(differences)
                assert np.max(np.min(np.linalg.norm(differences, axis=-1), axis=1)) < 1e-11
            assert info["minimum_periodic_pair_distance_A"] > 3.09
        checks[scenario] = {"result": "PASS", "candidate_count": len(plans),
                            "max_nonaffine_displacement_A": diagnostics["maximum_nonaffine_wrapped_displacement_A"]}
    # Large nonaffine motion must stop mapping rather than silently misassign.
    moved = original_frac.copy()
    moved[0, 0] += 0.2
    mapper.read_poscar = lambda path: (original_lat, np.mod(moved, 1), species, counts, selective)
    try:
        mapper.prepare(opts)
        raise AssertionError("Large-displacement guard did not trigger")
    except ValueError as exc:
        assert "exceeds --max-displacement" in str(exc)
    checks["large_displacement_refused"] = {"result": "PASS"}
    mapper.read_poscar = real_read
    translations = np.array([v for v in itertools.product(range(-2, 3), repeat=3) if any(v)]) @ lattice
    result = {"result": "PASS_INDEPENDENT_GEOMETRY_AND_SOFTWARE_AUDIT",
              "scientific_scope": "No physical relaxation, energies, stability, or electronic predictions verified.",
              "all_serialized_poscars_verified_against_parent_and_mutation": audit_count,
              "all_twelve_defect_types_present_in_every_parent": True,
              "ordered_inversion_occupancy_verified": True,
              "shortest_tested_periodic_translation_A": float(np.linalg.norm(translations, axis=1).min()),
              "vacancies_per_cation_plane": vacancies,
              "regenerator_software_checks": checks,
              "synthetic_check_note": "Geometry perturbations were in-memory software fixtures; no fabricated CONTCAR saved."}
    target = root / "independent_audit.json"
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"result": result["result"], "poscars": audit_count,
                      "software_checks": list(checks), "audit_file": str(target)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
