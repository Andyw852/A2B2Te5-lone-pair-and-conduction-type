#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compact key-results table for quick reading; all values are copied from the
detailed CSV/JSON outputs, not recomputed independently."""
from __future__ import annotations

import json

import pandas as pd

import audit_common as ac

rows = []


def add(category, item, value, evidence):
    rows.append({"category": category, "item": item, "value": value, "evidence": evidence})


inv = pd.read_csv(ac.CSV / "inventory.csv", encoding="utf-8-sig")
for _, r in inv.iterrows():
    add("inventory", r["material"] + "_mesh", f"{r['kpoints_mesh']} / nk={r['eigenval_nk']}",
        "step3 KPOINTS/EIGENVAL")
    add("inventory", r["material"] + "_nelect", str(r["eigenval_nelect"]), "EIGENVAL/OUTCAR/POTCAR")
    add("inventory", r["material"] + "_gap_eV", f"{r['bt2_dft_gap_eV']:.6f}", "BT2 JSON scissor.dft_gap")

mass = pd.read_csv(ac.CSV / "mass_fit_recheck.csv", encoding="utf-8-sig")
for _, r in mass.iterrows():
    add("mass", f"{r['material']}_{r['carrier']}_reliability", r["reliability"],
        "independent 125-point quadratic fit")
    add("mass", f"{r['material']}_{r['carrier']}_extremum_on_boundary",
        str(r["sampled_extremum_is_on_boundary"]), "energy minimum/maximum over the 5x5x5 mesh")
    add("mass", f"{r['material']}_{r['carrier']}_masses_me",
        f"{r['mass_1_me']}/{r['mass_2_me']}/{r['mass_3_me']}", "electron-energy Hessian ascending order")

bt = pd.read_csv(ac.CSV / "bt2_checks.csv", encoding="utf-8-sig")
for _, r in bt.iterrows():
    add("bt2", r["material"] + "_target_gap_eV", f"{r['target_gap_eV']:.4f}", "BT2 JSON")
    add("bt2", r["material"] + "_HSE_gap_eV", f"{r['HSE_gap_eV']:.6f}", "07 HSE path")
    add("bt2", r["material"] + "_c_at_mu0", r["carrier_conc_at_mu0_cm-3_by_T"],
        "BT2 JSON carrier_conc at mu=0")

for name, label in [("bt2_reliable_S_sigma.csv", "A"),
                    ("bt2_marginal_S_sigma.csv", "B"),
                    ("bt2_unreliable_S_sigma.csv", "C")]:
    d = pd.read_csv(ac.CSV / name, encoding="utf-8-sig")
    add("bt2_reliability", f"rows_{label}", str(len(d)), name)

s = pd.read_csv(ac.CSV / "s_projection_audit.csv", encoding="utf-8-sig")
for _, r in s[s.edge == "VBM"].iterrows():
    add("s_projection", r["material"] + "_VBM_50meV_below", str(r["window_is_below_VBM"]),
        "states in 入选态清单.csv")
    add("s_projection", r["material"] + "_A_site_s_pct", f"{r['A_site_s_pct']:.3f}", "带边能窗统计.csv")
    add("s_projection", r["material"] + "_total_s_pct", f"{r['total_s_pct']:.3f}", "带边能窗统计.csv")

dvo = json.load(open(ac.JSON / "dvo_audit_summary.json"))
add("dvo", "parents", str(dvo["parents"]), "manifest.json")
add("dvo", "defects", str(dvo["defects"]), "manifest.json")
add("dvo", "jobs_status", json.dumps(dvo["jobs_status_counts"]), "jobs.csv")
add("dvo", "DVO_has_actual_energies", str(dvo["conclusion"]["DVO_has_actual_energies"]),
    "validation/交付核验/no energy files")
add("dvo", "old_81atom_is_DVO", str(dvo["conclusion"]["old_81atom_is_DVO"]), "step1_bulk POSCAR 81 atoms")
add("dvo", "old_72atom_is_DVO", str(dvo["conclusion"]["old_72atom_is_DVO"]), "zip bulk POSCAR 72 atoms")

hc = json.load(open(ac.JSON / "hash_check_summary.json"))
add("hash", "checks", f"{hc['n_match']}/{hc['n_checks']} match", "hash_check.csv")
add("hash", "mismatches", str(hc["n_mismatch"]), "generator-version mismatches only")

ac.write_csv(ac.CSV / "key_results_summary.csv", rows)
print("wrote", ac.CSV / "key_results_summary.csv", "rows", len(rows))
