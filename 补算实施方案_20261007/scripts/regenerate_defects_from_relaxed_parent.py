#!/usr/bin/env python3
"""Map existing one-site defect proposals onto a consistently relaxed parent.

Strictly assumes CONTCAR preserves the original atom identities/order. Element
blocks can be checked; permutations of identical atoms cannot be proven absent.
Ghost positions are initial estimates, never measured empty-site coordinates.
No electronic charge settings, energies, relaxation, or DFT jobs are generated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_poscar(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    scale_fields = lines[1].split()
    if len(scale_fields) != 1 or float(scale_fields[0]) <= 0:
        raise ValueError("Only one positive POSCAR scale factor is supported; convert other conventions first")
    scale = float(scale_fields[0])
    lattice = np.array([[float(x) for x in line.split()[:3]] for line in lines[2:5]]) * scale
    if np.linalg.det(lattice) <= 0:
        raise ValueError("Expected a right-handed nondegenerate lattice")
    species = lines[5].split()
    try:
        counts = [int(x) for x in lines[6].split()]
    except ValueError as exc:
        raise ValueError("Expected a VASP5 POSCAR with element labels") from exc
    if len(species) != len(counts) or any(n <= 0 for n in counts):
        raise ValueError("Invalid element labels/counts")
    index = 7
    selective = lines[index].lower().startswith("s")
    if selective:
        index += 1
    mode = lines[index].strip().lower()
    if not mode.startswith(("d", "c", "k")):
        raise ValueError("Unsupported POSCAR coordinate mode")
    start = index + 1
    coords = np.array([[float(x) for x in line.split()[:3]]
                       for line in lines[start:start + sum(counts)]])
    if coords.shape != (sum(counts), 3):
        raise ValueError("Missing atom coordinates")
    if not mode.startswith("d"):
        coords = coords * scale @ np.linalg.inv(lattice)
    return lattice, np.mod(coords, 1.0), species, counts, selective


def minimum_image(frac_deltas, lattice):
    """Return rigorously bounded nearest images for a moderately conditioned cell."""
    shape = frac_deltas.shape
    reduced = np.asarray(frac_deltas).reshape(-1, 3)
    reduced = reduced - np.rint(reduced)
    best = reduced.copy()
    best_sq = np.sum((best @ lattice) ** 2, axis=1)
    smallest_singular_value = np.linalg.svd(lattice, compute_uv=False)[-1]
    bound = int(np.ceil(np.sqrt(best_sq.max(initial=0)) / smallest_singular_value + 0.5))
    if bound > 4:
        raise ValueError("Cell is too skewed for this utility; reduce the cell basis before mapping")
    for offset in itertools.product(range(-bound, bound + 1), repeat=3):
        trial = reduced + np.array(offset)
        norm_sq = np.sum((trial @ lattice) ** 2, axis=1)
        use = norm_sq < best_sq - 1e-20
        best[use] = trial[use]
        best_sq[use] = norm_sq[use]
    return best.reshape(shape)


def distance_matrix(frac, lattice):
    delta = minimum_image(frac[:, None, :] - frac[None, :, :], lattice)
    distances = np.linalg.norm(delta @ lattice, axis=-1)
    np.fill_diagonal(distances, np.inf)
    return distances


def output_poscar(title, lattice, frac, symbols, element_order):
    order = [i for element in element_order for i, sym in enumerate(symbols) if sym == element]
    counts = [symbols.count(element) for element in element_order]
    if any(n == 0 for n in counts):
        raise ValueError("Unexpected complete removal of an element block")
    lines = [title, "1.0"]
    lines.extend(" ".join(f"{x:.14f}" for x in v) for v in lattice)
    lines.extend([" ".join(element_order), " ".join(map(str, counts)), "Direct"])
    lines.extend(" ".join(f"{x:.14f}" for x in np.mod(frac[i], 1.0)) for i in order)
    return "\n".join(lines) + "\n"


def prepare(args):
    site_map_path = args.site_map.resolve()
    relaxed_path = args.relaxed_poscar.resolve()
    template_dir = args.templates.resolve()
    mapping = json.loads(site_map_path.read_text(encoding="utf-8"))
    parent_id = mapping["parent_id"]
    sites = mapping["sites"]
    original_lattice = np.array(mapping["lattice_vectors_A_rows"])
    occupied = sorted((s for s in sites if s["occupancy"] != "Vac"),
                      key=lambda s: s["poscar_atom_index_1based"])
    assert [s["poscar_atom_index_1based"] for s in occupied] == list(range(1, len(occupied) + 1))
    old_frac = np.array([s["fractional_coordinates"] for s in occupied])
    expected_symbols = [s["occupancy"] for s in occupied]
    lattice, frac, element_order, counts, selective = read_poscar(relaxed_path)
    actual_symbols = [symbol for symbol, count in zip(element_order, counts) for _ in range(count)]
    if actual_symbols != expected_symbols:
        raise ValueError("Relaxed parent species/count/block order differs from site_map; mapping refused")
    if args.min_distance <= 0 or args.max_displacement <= 0:
        raise ValueError("Both distance thresholds must be positive")
    displacement_frac = minimum_image(frac - old_frac, lattice)
    displacement_A = np.linalg.norm(displacement_frac @ lattice, axis=1)
    maximum = float(displacement_A.max())
    if maximum > args.max_displacement:
        raise ValueError(f"Maximum nonaffine wrapped displacement {maximum:.4f} A exceeds "
                         f"--max-displacement {args.max_displacement:g}. Inspect reconstruction/atom identity "
                         "before explicitly changing the threshold.")
    parent_distances = distance_matrix(frac, lattice)
    parent_minimum = float(parent_distances.min())
    if parent_minimum < args.min_distance:
        raise ValueError(f"Parent minimum separation {parent_minimum:.4f} A below safety threshold")
    sites_by_id = {site["site_id"]: site for site in sites}
    atom_by_site = {site["site_id"]: i for i, site in enumerate(occupied)}
    templates = sorted(template_dir.glob("*/defect.json"))
    if not templates:
        raise ValueError("--templates must contain defect subdirectories with defect.json")
    source_hashes = {"site_map_sha256": sha256(site_map_path), "relaxed_poscar_sha256": sha256(relaxed_path)}
    diagnostics = {
        "maximum_nonaffine_wrapped_displacement_A": maximum,
        "rms_nonaffine_wrapped_displacement_A": float(np.sqrt(np.mean(displacement_A ** 2))),
        "original_cell_volume_A3": float(np.linalg.det(original_lattice)),
        "input_parent_cell_volume_A3": float(np.linalg.det(lattice)),
        "input_parent_minimum_pair_distance_A": parent_minimum,
        "max_displacement_threshold_A": args.max_displacement,
        "min_distance_threshold_A": args.min_distance,
        "input_had_selective_dynamics_flags": selective,
        "selective_dynamics_note": "Output has no constraints; review relaxation settings separately.",
    }
    warnings = ["Same-element atom identity/order preservation is assumed and cannot be proven from CONTCAR alone.",
                "Passing the displacement threshold does not prove the parent retained its intended phase or occupation.",
                "No evidence of relaxation completion or force convergence is inferred from an input filename.",
                "Ghost coordinates are initial guesses; all defect candidates still require relaxation."]
    plans = []
    for template_path in templates:
        template = json.loads(template_path.read_text(encoding="utf-8"))
        if template["parent_id"] != parent_id:
            raise ValueError("Template parent_id does not match site_map")
        site_id = template["changed_site"]["site_id"]
        site = sites_by_id[site_id]
        source, target = template["source_occupancy"], template["target_occupancy"]
        if site["occupancy"] != source:
            raise ValueError("Template source occupancy disagrees with original site_map")
        new_frac = frac.copy()
        symbols = actual_symbols.copy()
        ghost_info = None
        if source == "Vac":
            if target == "Vac":
                raise ValueError("No-op defect template")
            ghost_frac = np.array(site["fractional_coordinates"])
            ghost_info = {"method": args.ghost_method,
                          "note": "Estimated initial position in input parent cell; not a known physical empty-site coordinate."}
            if args.ghost_method == "neighbor":
                te_atoms = [i for i, s in enumerate(occupied) if s["occupancy"] == "Te"]
                te_delta = minimum_image(old_frac[te_atoms] - ghost_frac, original_lattice)
                te_dist = np.linalg.norm(te_delta @ original_lattice, axis=1)
                nearest = np.argsort(te_dist)[:6]
                if len(nearest) != 6:
                    raise ValueError("Need six original neighboring Te atoms for neighbor interpolation")
                selected = [te_atoms[int(i)] for i in nearest]
                average_displacement = displacement_frac[selected].mean(axis=0)
                ghost_frac = ghost_frac + average_displacement
                ghost_info.update({"original_neighbor_site_ids": [occupied[i]["site_id"] for i in selected],
                                   "original_neighbor_distances_A": te_dist[nearest].tolist(),
                                   "mean_nonaffine_neighbor_displacement_fractional": average_displacement.tolist()})
            ghost_frac = np.mod(ghost_frac, 1.0)
            ghost_info["estimated_fractional_coordinates"] = ghost_frac.tolist()
            ghost_delta = minimum_image(frac - ghost_frac, lattice)
            ghost_minimum = float(np.linalg.norm(ghost_delta @ lattice, axis=1).min())
            minimum = min(parent_minimum, ghost_minimum)
            ghost_info["nearest_occupied_atom_distance_A"] = ghost_minimum
            new_frac = np.vstack((new_frac, ghost_frac))
            symbols.append(target)
        else:
            atom_index = atom_by_site[site_id]
            if target == "Vac":
                keep = [i for i in range(len(symbols)) if i != atom_index]
                new_frac = new_frac[keep]
                symbols = [symbols[i] for i in keep]
                minimum = float(parent_distances[np.ix_(keep, keep)].min())
            else:
                symbols[atom_index] = target
                minimum = parent_minimum
        if minimum < args.min_distance:
            raise ValueError(f"Collision threshold failed for {template['defect_id']}: {minimum:.4f} A")
        expected_count = Counter(actual_symbols)
        if source != "Vac":
            expected_count[source] -= 1
        if target != "Vac":
            expected_count[target] += 1
        assert Counter(symbols) == +expected_count
        metadata = {"defect_id": template["defect_id"], "parent_id": parent_id,
                    "status": "MAPPED_GEOMETRY_NOT_RUN", "intended_charge_state": 0,
                    "charge_note": "POSCAR does not define electronic charge; charged states remain separate jobs.",
                    "defect_type": template["defect_type"], "site_id": site_id,
                    "source_occupancy": source, "target_occupancy": target,
                    "input_parent": str(relaxed_path), "original_site_map": str(site_map_path),
                    "template": str(template_path), "template_sha256": sha256(template_path),
                    **source_hashes, "ghost_position_estimate": ghost_info,
                    "species_counts": dict(Counter(symbols)), "atom_count": len(symbols),
                    "minimum_periodic_pair_distance_A": minimum,
                    "parent_mapping_diagnostics": diagnostics, "warnings": warnings}
        text = output_poscar(f"NOT_RUN mapped {parent_id} {template['defect_id']}",
                             lattice, new_frac, symbols, element_order)
        plans.append((template["defect_id"], text, metadata))
    if len({p[0] for p in plans}) != len(plans):
        raise ValueError("Repeated candidate directory name")
    return plans, diagnostics, warnings, source_hashes


def main(args):
    output = args.out.resolve() if args.out else None
    if not args.check_only and output is None:
        raise ValueError("--out is required unless --check-only")
    if output is not None and output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty directory: {output}")
    plans, diagnostics, warnings, source_hashes = prepare(args)
    if not args.check_only:
        output.mkdir(parents=True, exist_ok=True)
        jobs = []
        for defect_id, poscar, metadata in plans:
            folder = output / defect_id
            folder.mkdir()
            (folder / "POSCAR").write_text(poscar, encoding="ascii")
            (folder / "mapped_defect.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            jobs.append({"defect_id": defect_id, "atom_count": metadata["atom_count"],
                         "intended_charge_state": 0, "status": "NOT_RUN", "poscar": f"{defect_id}/POSCAR"})
        manifest = {"status": "GEOMETRY_MAPPED_NO_DFT_PERFORMED", "candidate_count": len(plans),
                    "ghost_method": args.ghost_method, "diagnostics": diagnostics,
                    "source_hashes": source_hashes, "warnings": warnings}
        (output / "mapping_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        with (output / "jobs.csv").open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(jobs[0]))
            writer.writeheader()
            writer.writerows(jobs)
    print(json.dumps({"result": "CHECK_ONLY_PASS" if args.check_only else "MAPPED_GEOMETRY_WRITTEN_NOT_RUN",
                      "candidate_count": len(plans), "output": str(output) if output else None,
                      "diagnostics": diagnostics, "warnings": warnings}, ensure_ascii=False))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-map", type=Path, required=True)
    parser.add_argument("--relaxed-poscar", type=Path, required=True,
                        help="same atom identities/order as original parent; convergence must be checked separately")
    parser.add_argument("--templates", type=Path, required=True, help="one parent's directory containing */defect.json")
    parser.add_argument("--ghost-method", choices=("affine", "neighbor"), required=True,
                        help="explicit empty-site starting-position assumption, always requiring defect relaxation")
    parser.add_argument("--min-distance", type=float, default=1.5, help="heuristic collision threshold in A")
    parser.add_argument("--max-displacement", type=float, default=1.5,
                        help="maximum nonaffine wrapped displacement in A before mapping is refused")
    parser.add_argument("--out", type=Path, help="new or empty output folder")
    parser.add_argument("--check-only", action="store_true", help="validate/prepare mapping without writing files")
    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
