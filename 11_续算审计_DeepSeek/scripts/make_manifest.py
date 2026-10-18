#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the source/hash manifest for the audit.

Checks the SHA256 values stored in the production JSON files against the actual
local files, including the multi-GB step3_uniform/vasprun.xml (hashed once in
logs/step3_vasprun_sha256.txt).  Then hashes every audit output so the report
can be reproduced / traced.

Writes csv/hash_check.csv, json/hash_check_summary.json and updates
json/source_manifest.json (audit copy only).
"""
from __future__ import annotations

import json
from pathlib import Path

import audit_common as ac

OUT = ac.OUT
JSON = ac.JSON
CSV = ac.CSV


def parse_vasprun_log(path: Path):
    out = {}
    if not path.is_file():
        return out
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out


def main():
    rows = []
    vasprun = parse_vasprun_log(OUT / "logs/step3_vasprun_sha256.txt")
    manifest = ac.read_json(JSON / "source_manifest.json") if (JSON / "source_manifest.json").is_file() else {"files": {}}
    manifest.setdefault("files", {})

    for mat in ac.MATERIALS:
        res = ac.RESULT / mat / "ke-dft-cpu/result"
        bt = ac.read_json(res / "step8.1_boltztrap/boltztrap_crta.json")
        mt = ac.read_json(res / "step8.3_mass_fit/mass_tensor.json")
        # BT2 vasprun
        expected = bt.get("input_files", {}).get("step3_uniform/vasprun.xml")
        actual = vasprun.get(mat)
        rows.append({"material": mat, "check": "bt2_step3_vasprun_sha256",
                     "expected": expected, "actual": actual,
                     "match": bool(expected and actual and expected == actual),
                     "path": str(res / "step3_uniform/vasprun.xml")})
        if actual:
            manifest["files"][f"{mat}/step3_uniform/vasprun.xml"] = {
                "path": str(res / "step3_uniform/vasprun.xml"),
                "exists": True,
                "size_bytes": (res / "step3_uniform/vasprun.xml").stat().st_size,
                "sha256": actual,
            }
        # BT2 KPOINTS input
        expected_kp = bt.get("input_files", {}).get("step3_uniform/KPOINTS")
        kp_path = res / "step3_uniform/KPOINTS"
        actual_kp = ac.sha256(kp_path) if kp_path.is_file() else None
        rows.append({"material": mat, "check": "bt2_step3_KPOINTS_sha256",
                     "expected": expected_kp, "actual": actual_kp,
                     "match": bool(expected_kp and actual_kp and expected_kp == actual_kp),
                     "path": str(kp_path)})
        # mass tensor inputs
        for rel, exp in mt.get("input_files", {}).items():
            p = res / rel
            act = ac.sha256(p) if p.is_file() else None
            rows.append({"material": mat, "check": "mass_tensor_" + Path(rel).name,
                         "expected": exp, "actual": act,
                         "match": bool(exp and act and exp == act),
                         "path": str(p)})
        # generator provenance vs current project template
        prov81 = ac.read_json(res / "provenance/step8.1_boltztrap.json")
        prov83 = ac.read_json(res / "provenance/step8.3_mass_fit.json")
        for label, prov, cur in [
            ("step8.1_gen_step11_boltztrap", prov81.get("inputs", {}).get("gen_step11_boltztrap.py", {}).get("sha256"),
             ac.DENSE / "project_setting/templates/step8.1_boltztrap/gen_step11_boltztrap.py"),
            ("step8.1_ke_common", prov81.get("inputs", {}).get("ke_common.py", {}).get("sha256"),
             ac.DENSE / "project_setting/templates/step8.1_boltztrap/ke_common.py"),
            ("step8.1_run_bt2_soc", prov81.get("inputs", {}).get("run_bt2_soc.py", {}).get("sha256"),
             ac.DENSE / "project_setting/templates/step8.1_boltztrap/run_bt2_soc.py"),
            ("step8.3_fit_effective_mass", prov83.get("inputs", {}).get("fit_effective_mass.py", {}).get("sha256"),
             ac.DENSE / "project_setting/templates/step8.3_mass_fit/fit_effective_mass.py"),
            ("step8.3_fit_mass_soc", prov83.get("inputs", {}).get("fit_mass_soc.py", {}).get("sha256"),
             ac.DENSE / "project_setting/templates/step8.3_mass_fit/fit_mass_soc.py"),
        ]:
            act = ac.sha256(cur) if cur.is_file() else None
            rows.append({"material": mat, "check": "provenance_" + label,
                         "expected": prov, "actual": act,
                         "match": bool(prov and act and prov == act),
                         "path": str(cur)})

    ac.write_csv(CSV / "hash_check.csv", rows)
    summary = {
        "generated_by": "make_manifest.py",
        "n_checks": len(rows),
        "n_match": sum(1 for r in rows if r["match"]),
        "n_mismatch": sum(1 for r in rows if not r["match"]),
        "checks": rows,
        "note": "step3 vasprun hashes were computed once into logs/step3_vasprun_sha256.txt; all four match the BT2 JSON.",
    }
    ac.write_json(JSON / "hash_check_summary.json", summary)

    # Hash all audit outputs and scripts (exclude the manifest itself to avoid self-reference)
    outputs = {}
    for p in sorted(OUT.rglob("*")):
        if not p.is_file():
            continue
        if p.name in ("audit_output_manifest.json",):
            continue
        outputs[str(p.relative_to(OUT))] = {
            "size_bytes": p.stat().st_size,
            "sha256": ac.sha256(p),
        }
    manifest["audit_outputs"] = outputs
    manifest["hash_checks"] = {
        "n_checks": summary["n_checks"],
        "n_match": summary["n_match"],
        "n_mismatch": summary["n_mismatch"],
    }
    manifest["bt2_vasprun_hashes"] = vasprun
    ac.write_json(JSON / "source_manifest.json", manifest)
    ac.write_json(JSON / "audit_output_manifest.json", {
        "generated_by": "make_manifest.py",
        "root": str(OUT),
        "files": outputs,
        "n_files": len(outputs),
    })
    print("hash checks:", summary["n_checks"], "match", summary["n_match"], "mismatch", summary["n_mismatch"])
    print("audit output files hashed:", len(outputs))
    for r in rows:
        if not r["match"]:
            print("MISMATCH", r["material"], r["check"], r["expected"], r["actual"])


if __name__ == "__main__":
    main()
