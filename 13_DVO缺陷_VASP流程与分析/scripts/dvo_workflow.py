#!/usr/bin/env python3
"""VASP workflow and analysis for the Pb2Sb2Te5 D/V/O defect comparison.

Builds on 补算实施方案_20261007/rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED (9 parents,
234 defect templates) and its regenerate_defects_from_relaxed_parent.py. Nothing
here runs VASP; it writes inputs, collects finished outputs and analyses them.

  parents  : relax+static job folders for the 9 parent cells
  refs     : chemical-potential reference phases (Pb, Sb, Te, PbTe, Sb2Te3, 9-atom layered Pb2Sb2Te5)
  defects  : map templates onto relaxed parents (via the existing tool) and write defect jobs
  collect  : read finished static/OUTCAR files into energies.csv (+ core potentials)
  analyze  : neutral formation energies, D/V/O distributions, dominant defects,
             electron-counting doping tendency; charged states and frozen-defect
             neutrality when charged jobs and parent DOSCARs exist
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUTS = HERE.parent / "inputs"
PACKAGE = ROOT / "补算实施方案_20261007"
DVO = PACKAGE / "rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED"
REGEN = PACKAGE / "scripts" / "regenerate_defects_from_relaxed_parent.py"
POTCAR_SRC = sorted((ROOT / "A2B2Te5_SOC_dense" / "materials").glob("*/POTCAR"))
LAYERED = ROOT / "A2B2Te5_SOC_dense" / "materials" / "{material}" / "POSCAR"
KB = 8.617333262e-5
# Electron counting with Pb2+, Sb3+, Te2-: + donor, - acceptor (the literature's classification).
NOMINAL_Q = {"Sb_Pb": 1, "Pb_Sb": -1, "Pb_on_Vac": 2, "Sb_on_Vac": 3, "Te_on_Vac": -2,
             "V_Pb": -2, "V_Sb": -3, "V_Te": 2, "Pb_Te": -2, "Sb_Te": -1, "Te_Pb": 2, "Te_Sb": 1}
KEY_TYPES = ("Sb_Pb", "Pb_Sb", "Pb_on_Vac", "Sb_on_Vac", "V_Te")
# Bi sits on the B site like Sb; defect names from the model generator use the B element.
canonical = lambda t: t.replace("Bi", "Sb")
FAMILY_COLOR = {"D": "#2a78d6", "V": "#eb6834", "O": "#1baf7a"}
FAMILY_MARK = {"D": "o", "V": "s", "O": "D"}


# ---------------------------------------------------------------- inputs
def potcar_blocks():
    """Same PAW datasets as the project's SOC runs: Pb_d, Sn_d, Sb, Bi_d, Te."""
    out = {}
    for src in POTCAR_SRC:
        for b in src.read_text().split("End of Dataset\n"):
            if b.strip():
                title = re.search(r"TITEL\s*=\s*PAW_PBE\s+(\S+)", b).group(1)
                out.setdefault(title.split("_")[0], b + "End of Dataset\n")
    return out


def zval():
    return {e: float(re.search(r"ZVAL\s*=\s*([\d.]+)", b).group(1)) for e, b in potcar_blocks().items()}


def read_poscar(path):
    lines = Path(path).read_text().splitlines()
    scale = float(lines[1].split()[0])
    lattice = np.array([[float(x) for x in l.split()[:3]] for l in lines[2:5]]) * scale
    species, counts = lines[5].split(), [int(x) for x in lines[6].split()]
    i = 8 if lines[7].lower().startswith("s") else 7
    n = sum(counts)
    coords = np.array([[float(x) for x in l.split()[:3]] for l in lines[i + 1:i + 1 + n]])
    if not lines[i].strip().lower().startswith("d"):
        coords = coords * scale @ np.linalg.inv(lattice)
    symbols = [s for s, c in zip(species, counts) for _ in range(c)]
    return lattice, np.mod(coords, 1.0), species, counts, symbols


def write_poscar(path, title, lattice, frac, symbols):
    order = list(dict.fromkeys(symbols))
    lines = [title, "1.0"] + [" ".join(f"{x:.12f}" for x in v) for v in lattice]
    lines += [" ".join(order), " ".join(str(symbols.count(e)) for e in order), "Direct"]
    lines += [" ".join(f"{x:.12f}" for x in frac[i]) for e in order for i, s in enumerate(symbols) if s == e]
    Path(path).write_text("\n".join(lines) + "\n")


def fill(template, **kw):
    text = (INPUTS / template).read_text()
    for k, v in kw.items():
        text = text.replace("{" + k + "}", str(v))
    # INCAR comments are dropped so that only ASCII keywords reach VASP.
    keep = lambda l: l.strip() != "" and not (template.startswith("INCAR") and l.lstrip().startswith("#"))
    return "\n".join(l for l in text.splitlines() if keep(l)) + "\n"


def write_job(folder, poscar_src, name, relax_tpl, static_tpl, args, kpoints=True, ispin=1, nelect=None, meta=None):
    folder.mkdir(parents=True, exist_ok=False)
    shutil.copy(poscar_src, folder / "POSCAR")
    species = Path(poscar_src).read_text().splitlines()[5].split()
    if not args.no_potcar:
        blocks = potcar_blocks()
        (folder / "POTCAR").write_text("".join(blocks[e] for e in species))
    nel = f"NELECT  = {nelect}" if nelect is not None else ""
    (folder / "INCAR.relax").write_text(fill(relax_tpl, name=name, ispin=ispin, nelect=nel))
    (folder / "INCAR.static").write_text(fill(static_tpl, name=name, ispin=ispin, nelect=nel))
    if kpoints:
        k = args.kmesh
        (folder / "KPOINTS").write_text(f"D/V/O common mesh\n0\n{'Monkhorst-Pack' if args.shift else 'Gamma'}\n{k[0]} {k[1]} {k[2]}\n0 0 0\n")
    sub = fill("submit.sh.template", partition=args.partition, jobname=name[:60], nodes=args.nodes,
               ntasks=args.ntasks, vasp=args.vasp, restarts=args.restarts)
    (folder / "submit.sh").write_text(sub)
    (folder / "job.json").write_text(json.dumps(meta or {}, indent=1, ensure_ascii=False) + "\n")


def neutral_nelect(poscar):
    _, _, species, counts, _ = read_poscar(poscar)
    z = zval()
    return int(round(sum(z[s] * c for s, c in zip(species, counts))))


def cmd_parents(args):
    out = Path(args.out) / "parents"
    rows = []
    for pdir in sorted((Path(args.models) / "parents").iterdir()):
        meta = json.loads((pdir / "parent.json").read_text())
        if meta["family"] not in args.families:
            continue
        nel = neutral_nelect(pdir / "POSCAR")
        write_job(out / pdir.name, pdir / "POSCAR", f"PbSbTe_{pdir.name}", "INCAR.parent_relax", "INCAR.static", args,
                  meta={"kind": "parent", "parent_id": pdir.name, "family": meta["family"], "nelect": nel,
                        "models": str(Path(args.models).resolve()),
                        "source_poscar_sha256": meta["fingerprint_sha256"], "status": "NOT_RUN"})
        rows.append((pdir.name, meta["family"], nel))
    print(json.dumps({"written": len(rows), "out": str(out), "parents": rows}, ensure_ascii=False))


def ref_structures():
    from ase.spacegroup import crystal
    from ase.build import bulk
    refs = {
        "Pb": bulk("Pb", "fcc", a=4.95),
        "Sb": crystal("Sb", [(0, 0, 0.2336)], spacegroup=166, cellpar=[4.307, 4.307, 11.273, 90, 90, 120], primitive_cell=True),
        "Te": crystal("Te", [(0.2636, 0, 1 / 3)], spacegroup=152, cellpar=[4.457, 4.457, 5.929, 90, 90, 120], primitive_cell=True),
        "PbTe": bulk("PbTe", "rocksalt", a=6.46),
        "Sb2Te3": crystal(["Sb", "Te", "Te"], [(0, 0, 0.3988), (0, 0, 0), (0, 0, 0.7872)], spacegroup=166,
                          cellpar=[4.264, 4.264, 30.458, 90, 90, 120], primitive_cell=True),
    }
    return refs


def cmd_refs(args):
    out = Path(args.out) / "refs"
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "_poscar"
    tmp.mkdir(exist_ok=True)
    made = []
    drop = {"Sb", "Sb2Te3"} if "Bi" in args.material else set()
    for name, atoms in ref_structures().items():
        if name in drop:
            continue
        p = tmp / f"{name}.vasp"
        write_poscar(p, f"{name} literature start (relaxed with ISIF=3)", atoms.cell[:], atoms.get_scaled_positions(), atoms.get_chemical_symbols())
        made.append((name, p))
    made.append((f"{args.material}_layered9", Path(str(LAYERED).format(material=args.material))))
    if "Bi" in args.material:
        from ase.spacegroup import crystal
        bi = crystal("Bi", [(0, 0, 0.2339)], spacegroup=166, cellpar=[4.546, 4.546, 11.862, 90, 90, 120], primitive_cell=True)
        b2 = crystal(["Bi", "Te", "Te"], [(0, 0, 0.4001), (0, 0, 0), (0, 0, 0.7907)], spacegroup=166,
                     cellpar=[4.386, 4.386, 30.497, 90, 90, 120], primitive_cell=True)
        for name, atoms in (("Bi", bi), ("Bi2Te3", b2)):
            p = tmp / f"{name}.vasp"
            write_poscar(p, f"{name} literature start (relaxed with ISIF=3)", atoms.cell[:], atoms.get_scaled_positions(), atoms.get_chemical_symbols())
            made.append((name, p))
    for name, p in made:
        write_job(out / name, p, f"ref_{name}", "INCAR.ref_relax", "INCAR.ref_static", args, kpoints=False,
                  meta={"kind": "ref", "ref_id": name, "formula": dict(Counter(read_poscar(p)[4])), "status": "NOT_RUN"})
    shutil.rmtree(tmp)
    print(json.dumps({"written": [m[0] for m in made], "out": str(out)}, ensure_ascii=False))


def finished(job):
    st = (job / "STATUS")
    return st.exists() and st.read_text().strip() == "DONE"


def cmd_defects(args):
    pdir = Path(args.parents)
    out = Path(args.out) / "defects"
    charges = {}
    if args.charges:
        for line in Path(args.charges).read_text(encoding="utf-8-sig").splitlines()[1:]:
            if line.strip() and not line.startswith("#"):
                pid, did, qs = [x.strip() for x in line.split(",", 2)]
                charges[(pid, did)] = [int(q) for q in qs.replace('"', "").split()]
    made, skipped, existing = [], [], []
    for job in sorted(pdir.iterdir()):
        meta = json.loads((job / "job.json").read_text())
        pid = meta["parent_id"]
        if meta["family"] not in args.families:
            continue
        if not finished(job):
            skipped.append((pid, "parent STATUS is not DONE"))
            continue
        models = Path(meta.get("models", args.models))
        if not models.exists():
            models = Path(args.models)
        mapped = Path(args.out) / "_mapped" / pid
        if not mapped.exists():
            cmd = [sys.executable, str(REGEN), "--site-map", str(models / "parents" / pid / "site_map.json"),
                   "--relaxed-poscar", str(job / "relax" / "CONTCAR"), "--templates", str(models / "defects" / pid),
                   "--ghost-method", "neighbor", "--out", str(mapped)]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode:
                skipped.append((pid, "mapping refused: " + res.stderr.strip().splitlines()[-1]))
                continue
        for ddir in sorted(p for p in mapped.iterdir() if p.is_dir()):
            dm = json.loads((ddir / "mapped_defect.json").read_text(encoding="utf-8"))
            dtype, site_no = dm["defect_type"], int(dm["defect_id"].split("__s")[-1])
            if canonical(dtype) not in {canonical(x) for x in args.types} or site_no > args.max_sites:
                continue
            nel = neutral_nelect(ddir / "POSCAR")
            dn = {e: dm["species_counts"].get(e, 0) - c for e, c in Counter(read_poscar(job / "relax" / "CONTCAR")[4]).items()}
            for q in [0] if args.charges is None else charges.get((pid, dm["defect_id"]), []):
                name = f"{pid}__{dm['defect_id']}__q{q:+d}"
                if (out / name).exists():
                    existing.append(name)
                    continue
                ispin = 2 if (args.spin == "auto" and (nel - q) % 2) else 1
                write_job(out / name, ddir / "POSCAR", name, "INCAR.defect_relax", "INCAR.static", args,
                          ispin=ispin, nelect=(nel - q) if q else None,
                          meta={"kind": "defect", "parent_id": pid, "family": meta["family"], "defect_id": dm["defect_id"],
                                "defect_type": dtype, "site_id": dm["site_id"], "q": q, "nominal_q": NOMINAL_Q[canonical(dtype)],
                                "models": str(models),
                                "delta_n": dn, "nelect": nel - q, "ghost_position_estimate": dm["ghost_position_estimate"],
                                "status": "NOT_RUN"})
                made.append(name)
    print(json.dumps({"written": len(made), "already_present": len(existing), "skipped_parents": skipped, "out": str(out)}, ensure_ascii=False))


# ---------------------------------------------------------------- collect
def outcar_info(path):
    text = Path(path).read_text(errors="replace")
    e = re.findall(r"energy\(sigma->0\)\s*=\s*([-\d.]+)", text)
    nel = re.findall(r"NELECT\s*=\s*([\d.]+)", text)
    blocks = text.split("average (electrostatic) potential at core")
    core = None
    if len(blocks) > 1:
        lines, body = blocks[-1].splitlines(), []
        for l in lines[1:]:
            if re.match(r"^\s*\d+\s*-?\d+\.\d+", l):
                body.append(l)
            elif body:
                break
        core = [float(v) for _, v in re.findall(r"(\d+)\s*(-?\d+\.\d+)", " ".join(body))]
    mag = re.findall(r"number of electron\s+[\d.]+\s+magnetization\s+([-\d.]+)", text)
    return {"energy": float(e[-1]) if e else None, "nelect": float(nel[-1]) if nel else None,
            "finished": "General timing" in text, "core": core, "mag": float(mag[-1]) if mag else 0.0}


def min_image(df, lattice):
    df = df - np.rint(df)
    best = np.full(df.shape[:-1], np.inf)
    for o in np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3):
        best = np.minimum(best, np.linalg.norm((df + o) @ lattice, axis=-1))
    return best


def cmd_collect(args):
    rows, cores = [], {}
    for jf in sorted(Path(args.root).rglob("job.json")):
        job, meta = jf.parent, json.loads(jf.read_text(encoding="utf-8"))
        st = job / "static" / "OUTCAR"
        relax = job / "relax" / "OUTCAR"
        row = {"job": str(job.relative_to(args.root)), **{k: meta.get(k) for k in ("kind", "parent_id", "family", "ref_id", "defect_id", "defect_type", "site_id", "q", "nominal_q", "models")}}
        row["delta_n"] = json.dumps(meta.get("delta_n", {}))
        row["relax_converged"] = relax.exists() and "reached required accuracy" in relax.read_text(errors="replace")
        if not st.exists():
            row["status"] = "NO_STATIC"
            rows.append(row)
            continue
        info = outcar_info(st)
        lat, frac, _, _, sym = read_poscar(job / "static" / "POSCAR")
        ok = info["finished"] and info["energy"] is not None
        row.update({"status": ("DONE" if row["relax_converged"] else "RELAX_NOT_CONVERGED") if ok else "STATIC_INCOMPLETE",
                    "energy_eV": info["energy"], "nelect": info["nelect"], "n_atoms": len(sym),
                    "volume_A3": abs(np.linalg.det(lat)), "magnetization": info["mag"],
                    "formula": json.dumps(dict(Counter(sym)))})
        if meta.get("kind") == "defect" and meta.get("ghost_position_estimate"):
            # Did the atom placed into a structural vacancy stay there?
            g = np.array(meta["ghost_position_estimate"]["estimated_fractional_coordinates"])
            row["nearest_atom_to_filled_site_A"] = float(min_image(frac - g, lat).min())
        if info["core"] is not None:
            cores[row["job"]] = {"core": info["core"], "frac": frac.tolist(), "lattice": lat.tolist(), "symbols": sym}
        rows.append(row)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(args.out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    Path(args.out).with_suffix(".core.json").write_text(json.dumps(cores))
    print(json.dumps(Counter(r["status"] for r in rows), ensure_ascii=False))


# ---------------------------------------------------------------- analyze
def chemical_potentials(refs, override):
    """Absolute mu (eV/atom). Default: Pb2B2Te5 bounded by PbTe + B2Te3 (its pseudo-binary), Te-rich and Te-poor ends."""
    if override:
        return json.loads(Path(override).read_text(encoding="utf-8"))
    B = "Bi" if "Bi2Te3" in refs else "Sb"
    if B == "Bi":
        refs = {k.replace("Bi", "Sb"): {"E": v["E"], "n": {e.replace("Bi", "Sb"): c for e, c in v["n"].items()}} for k, v in refs.items()}
    need = {"Pb", "Sb", "Te", "PbTe", "Sb2Te3"}
    if not need <= set(refs):
        missing = sorted(need - set(refs))
        raise SystemExit(f"reference energies missing for {missing}; finish refs or pass --mu mu.json")
    el = {e: refs[e]["E"] / refs[e]["n"][e] for e in ("Pb", "Sb", "Te")}
    dH = lambda r: refs[r]["E"] - sum(c * el[e] for e, c in refs[r]["n"].items())
    h_pbte, h_sb2te3 = dH("PbTe"), dH("Sb2Te3")
    conds = {}
    for name, dte in (("Te-rich", 0.0), ("Te-poor", max(h_pbte, h_sb2te3 / 3))):
        dpb, dsb = h_pbte - dte, (h_sb2te3 - 3 * dte) / 2
        conds[name] = {"Pb": el["Pb"] + dpb, "Sb": el["Sb"] + dsb, "Te": el["Te"] + dte,
                       "delta": {"Pb": dpb, "Sb": dsb, "Te": dte}}
    info = {"B_element": B, "dH_PbTe_eV_fu": h_pbte, f"dH_{B}2Te3_eV_fu": h_sb2te3}
    if "Pb2Sb2Te5_layered9" in refs:
        info["E_layered_minus_2PbTe_B2Te3_eV_fu"] = refs["Pb2Sb2Te5_layered9"]["E"] - 2 * refs["PbTe"]["E"] / refs["PbTe"]["n"]["Pb"] - refs["Sb2Te3"]["E"] / (refs["Sb2Te3"]["n"]["Sb"] / 2)
    if B == "Bi":
        conds = {k: {**{e.replace("Sb", "Bi"): v for e, v in c.items() if e != "delta"}, "delta": {e.replace("Sb", "Bi"): x for e, x in c["delta"].items()}} for k, c in conds.items()}
    return {"conditions": conds, "info": info}


def sites_per_cell(pid, models=None):
    """Candidate sites of each canonical defect type in one parent cell."""
    models = Path(models) if models and Path(models).exists() else DVO
    occ = Counter(s["occupancy"] for s in json.loads((models / "parents" / pid / "site_map.json").read_text())["sites"])
    nb = occ["Sb"] + occ["Bi"]
    return {"Sb_Pb": occ["Pb"], "Pb_Sb": nb, "Pb_on_Vac": occ["Vac"], "Sb_on_Vac": occ["Vac"], "Te_on_Vac": occ["Vac"],
            "V_Pb": occ["Pb"], "V_Sb": nb, "V_Te": occ["Te"], "Pb_Te": occ["Te"], "Sb_Te": occ["Te"],
            "Te_Pb": occ["Pb"], "Te_Sb": nb}


def alignment(core, host_job, def_job, exclude_A):
    """Mean core-potential shift of atoms far from the defect; atoms matched to the host by position."""
    h, d = core.get(host_job), core.get(def_job)
    if not h or not d:
        return 0.0, 0
    lat = np.array(d["lattice"])
    hf, df = np.array(h["frac"]), np.array(d["frac"])
    dist = min_image(df[:, None, :] - hf[None, :, :], lat)
    match = dist.argmin(1)
    good = dist[np.arange(len(df)), match] < 0.6
    hs = np.array(h["symbols"])[match]
    good &= hs == np.array(d["symbols"])
    # Defect centre: host atoms without a partner, or defect atoms without a host partner.
    lonely = [df[i] for i in range(len(df)) if not good[i]] + [hf[j] for j in set(range(len(hf))) - set(match[good])]
    if not lonely:
        return 0.0, 0
    far = good & (min_image(df[:, None, :] - np.array(lonely)[None], lat).min(1) > exclude_A)
    dv = np.array(d["core"])[far] - np.array(h["core"])[match[far]]
    return float(dv.mean()) if far.any() else 0.0, int(far.sum())


def host_dos(job_dir):
    p = Path(job_dir) / "static" / "DOSCAR"
    if not p.exists():
        return None
    lines = p.read_text().splitlines()
    nedos, efermi = int(float(lines[5].split()[2])), float(lines[5].split()[3])
    a = np.array([[float(x) for x in l.split()[:3]] for l in lines[6:6 + nedos]])
    # Reference electron count = integrated DOS at the neutral host Fermi level (robust to a truncated window).
    return a[:, 0], a[:, 1], efermi, float(np.interp(efermi, a[:, 0], a[:, 2]))


def solve_ef(T, dos, nel_ref, vol_cm3, defect_terms, frozen=None):
    """Charge neutrality: sum_q q c_q + (N0 - N(EF))/V = 0; returns EF (absolute), free-carrier net (n-p) and c_q."""
    E, g = dos[0], dos[1]
    dE = E[1] - E[0]

    def occupation(mu):
        x = np.clip((E - mu) / (KB * T), -200, 200)
        return float((g / (1 + np.exp(x))).sum() * dE)

    def state(mu):
        conc = {}
        for key, (sites, levels) in defect_terms.items():
            w = {q: math.exp(-min(max((e0 + q * mu) / (KB * T), -700), 700)) for q, e0 in levels.items()}
            z = sum(w.values())
            if frozen is None:
                conc[key] = {q: sites * wq for q, wq in w.items()}
            else:
                conc[key] = {q: frozen[key] * wq / z for q, wq in w.items()}
        net_free = (occupation(mu) - nel_ref) / vol_cm3
        charge = sum(q * c for v in conc.values() for q, c in v.items()) - net_free
        return charge, net_free, conc

    lo, hi = E[0] + 0.5, E[-1] - 0.5
    from scipy.optimize import brentq
    mu = brentq(lambda m: state(m)[0], lo, hi, xtol=1e-6)
    _, net, conc = state(mu)
    return mu, net, conc


def cmd_analyze(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = [r for r in csv.DictReader(open(args.energies, encoding="utf-8-sig")) if r.get("status") == "DONE"]
    core = json.loads(Path(args.energies).with_suffix(".core.json").read_text()) if Path(args.energies).with_suffix(".core.json").exists() else {}
    refs = {r["ref_id"]: {"E": float(r["energy_eV"]), "n": json.loads(r["formula"])} for r in rows if r["kind"] == "ref"}
    hosts = {r["parent_id"]: r for r in rows if r["kind"] == "parent"}
    mu = chemical_potentials(refs, args.mu)
    conds = mu["conditions"]
    root = Path(args.energies).parent if args.root is None else Path(args.root)
    defects = []
    for r in rows:
        if r["kind"] != "defect" or r["parent_id"] not in hosts:
            continue
        h = hosts[r["parent_id"]]
        q = int(r["q"])
        dn = json.loads(r["delta_n"])
        dv, nfar = alignment(core, h["job"], r["job"], args.exclude) if q else (0.0, 0)
        L = float(h["volume_A3"]) ** (1 / 3)
        e_mp = q * q * 2.8373 * 14.39965 / (2 * args.eps * L) if q else 0.0
        base = float(r["energy_eV"]) - float(h["energy_eV"]) + q * dv + e_mp
        d = {**r, "q": q, "dE_raw": float(r["energy_eV"]) - float(h["energy_eV"]), "align_dV": dv, "align_atoms": nfar, "E_MP": e_mp}
        for c, m in conds.items():
            # Ef(EF) = base - sum dn*mu + q*EF_abs; EF measured from the host Fermi level of its own static run.
            d[f"Ef0_{c}"] = base - sum(n * m[e] for e, n in dn.items())
        defects.append(d)
    if not defects:
        raise SystemExit("no finished defect rows")
    # Charged rows: Ef(EF) = Ef0 + q*EF_abs. For reading, also give the value at the host's own Fermi level.
    efermi = {pid: (host_dos(root / h["job"]) or (None,) * 4)[2] for pid, h in hosts.items()}
    for d in defects:
        if d["q"] and efermi.get(d["parent_id"]) is not None:
            for c in conds:
                d[f"Ef_at_hostEF_{c}"] = d[f"Ef0_{c}"] + d["q"] * efermi[d["parent_id"]]
    with open(out / "形成能_逐位点.csv", "w", newline="", encoding="utf-8-sig") as f:
        keys = list(dict.fromkeys(k for d in defects for k in d))
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(defects)
    neutral = [d for d in defects if d["q"] == 0]
    summary, dominant = [], {}
    for c in conds:
        for fam in sorted({d["family"] for d in neutral}):
            per_type = defaultdict(list)
            for d in neutral:
                if d["family"] == fam:
                    per_type[d["defect_type"]].append(d)
            for t, ds in sorted(per_type.items()):
                e = np.array([x[f"Ef0_{c}"] for x in ds])
                rec = {"condition": c, "family": fam, "defect_type": t, "nominal_q": NOMINAL_Q[canonical(t)], "n_sites_sampled": len(e),
                       "n_parents": len({x["parent_id"] for x in ds}), "Ef0_min": e.min(), "Ef0_median": float(np.median(e)), "Ef0_max": e.max(),
                       "inserted_atom_moved_gt1A": sum(1 for x in ds if x.get("nearest_atom_to_filled_site_A") and float(x["nearest_atom_to_filled_site_A"]) > 1.0)}
                for T in args.T_anneal:
                    # Configurational average over sampled sites: Ef_eff = -kT ln <exp(-Ef/kT)>.
                    rec[f"Ef0_eff_{int(T)}K"] = float(-KB * T * np.log(np.mean(np.exp(-(e - e.min()) / (KB * T)))) + e.min())
                summary.append(rec)
            for T in args.T_anneal:
                fam_rows = [s for s in summary if s["condition"] == c and s["family"] == fam]
                best = sorted(fam_rows, key=lambda s: s[f"Ef0_eff_{int(T)}K"])[:3]
                vol = np.mean([float(hosts[p]["volume_A3"]) for p in {d["parent_id"] for d in neutral if d["family"] == fam}]) * 1e-24
                d0 = next(d for d in neutral if d["family"] == fam)
                nsite = sites_per_cell(d0["parent_id"], d0.get("models"))
                net = sum(s["nominal_q"] * nsite[canonical(s["defect_type"])] / vol * math.exp(-max(s[f"Ef0_eff_{int(T)}K"], 0) / (KB * T)) for s in fam_rows)
                dominant[f"{c}|{fam}|{int(T)}K"] = {"lowest_three": [(b["defect_type"], round(b[f"Ef0_eff_{int(T)}K"], 3), b["nominal_q"]) for b in best],
                                                     "nominal_net_donor_minus_acceptor_cm3": net,
                                                     "tendency": "n" if net > 0 else "p",
                                                     "negative_Ef_present": any(s["Ef0_min"] < 0 for s in fam_rows)}
    with open(out / "形成能_家族统计.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    # mu-independent comparison of the same defect type between families (same delta_n).
    shifts = []
    c0 = next(iter(conds))
    for t in sorted({s["defect_type"] for s in summary}):
        by = {s["family"]: s for s in summary if s["condition"] == c0 and s["defect_type"] == t}
        if "D" in by:
            for fam in ("V", "O"):
                if fam in by:
                    shifts.append({"defect_type": t, "compare": f"{fam}-D", "dEf_min": by[fam]["Ef0_min"] - by["D"]["Ef0_min"],
                                   "dEf_median": by[fam]["Ef0_median"] - by["D"]["Ef0_median"]})
    hostE = {}
    for fam in sorted({h["family"] for h in hosts.values()}):
        e = [float(h["energy_eV"]) / 16 for h in hosts.values() if h["family"] == fam]
        hostE[fam] = {"E_per_fu_mean": float(np.mean(e)), "E_per_fu_min": float(np.min(e)), "n": len(e)}
    if "O" in hostE:
        for fam in hostE:
            hostE[fam]["minus_O_meV_per_fu"] = 1000 * (hostE[fam]["E_per_fu_mean"] - hostE["O"]["E_per_fu_mean"])
    charged = [d for d in defects if d["q"] != 0]
    neutrality = {}
    if charged and args.dos:
        neutrality = frozen_neutrality(args, root, hosts, defects, conds, efermi)
    (out / "DVO_结论摘要.json").write_text(json.dumps({"chemical_potentials": mu, "dominant_by_condition_family_T": dominant,
        "same_type_shift_between_families": shifts, "host_energy_per_formula_unit": hostE, "charged_neutrality": neutrality,
        "notes": ["Neutral Ef do not depend on EF; dominant = lowest configurational-average Ef0 among sampled sites",
                  "tendency uses electron counting (Pb2+, Sb3+, Te2-) with full ionisation; it is the literature-level classification, not a carrier prediction",
                  "same_type_shift is mu-independent because both families add/remove the same atoms",
                  "negative_Ef_present: the dilute limit fails for that family/condition; report ranking only",
                  "charged_neutrality appears only with charged jobs and --dos; it needs the parent static DOSCAR"]},
        indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    plot_strips(plt, out, neutral, conds, args)
    print(json.dumps({k: v["lowest_three"][:1] + [v["tendency"]] for k, v in dominant.items()}, ensure_ascii=False, indent=0))


def frozen_neutrality(args, root, hosts, defects, conds, efermi):
    """Equilibrate charge states and concentrations at T_anneal, then freeze totals and re-solve at T_room."""
    res = {}
    for fam in sorted({d["family"] for d in defects}):
        # Only parents that have charged states can move EF; neutral-only parents are skipped, not labelled.
        pids = sorted({d["parent_id"] for d in defects if d["family"] == fam and d["q"] != 0})
        for pid in pids:
            dos = host_dos(root / hosts[pid]["job"])
            if dos is None:
                continue
            h = hosts[pid]
            vol = float(h["volume_A3"]) * 1e-24
            nel_ref = dos[3]
            nsite = sites_per_cell(pid, h.get("models"))
            for c in conds:
                terms = {}
                for d in defects:
                    if d["parent_id"] != pid:
                        continue
                    key = d["defect_id"]
                    sites = nsite[canonical(d["defect_type"])] / vol / max(1, sum(1 for x in defects if x["parent_id"] == pid and x["defect_type"] == d["defect_type"] and x["q"] == 0))
                    # e0 + q*EF_abs; the static energies already contain the host-referenced electron reservoir at EF=0.
                    terms.setdefault(key, (sites, {}))[1][d["q"]] = d[f"Ef0_{c}"]
                for Ta in args.T_anneal:
                    mu_a, net_a, conc_a = solve_ef(Ta, dos, nel_ref, vol, terms)
                    frozen = {k: sum(v.values()) for k, v in conc_a.items()}
                    mu_r, net_r, _ = solve_ef(args.T_room, dos, nel_ref, vol, terms, frozen)
                    res[f"{pid}|{c}|{int(Ta)}K"] = {"EF_anneal_minus_host_EF_eV": mu_a - dos[2], "net_free_n_minus_p_anneal_cm3": net_a,
                                                    "EF_room_minus_host_EF_eV": mu_r - dos[2], "net_free_n_minus_p_room_cm3": net_r,
                                                    "type_room": "undetermined" if abs(net_r) < args.min_carriers else ("n" if net_r > 0 else "p"),
                                                    "charged_sites": sorted({d["defect_id"] for d in defects if d["parent_id"] == pid and d["q"] != 0})}
    return res


def plot_strips(plt, out, neutral, conds, args):
    present = {d["defect_type"] for d in neutral}
    types = [x for t in list(KEY_TYPES) + sorted(NOMINAL_Q) for x in present if canonical(x) == t]
    types = list(dict.fromkeys(types))
    fams = [f for f in ("D", "V", "O") if f in {d["family"] for d in neutral}]
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e"})
    fig, axs = plt.subplots(len(conds), 1, figsize=(max(7, 0.85 * len(types) + 2), 3.6 * len(conds)), layout="constrained", squeeze=False)
    for ax, c in zip(axs[:, 0], conds):
        for k, fam in enumerate(fams):
            labelled = False
            for i, t in enumerate(types):
                e = [d[f"Ef0_{c}"] for d in neutral if d["family"] == fam and d["defect_type"] == t]
                if e:
                    x = i + (k - (len(fams) - 1) / 2) * 0.27 + np.linspace(-0.09, 0.09, len(e))
                    ax.plot(x, e, FAMILY_MARK[fam], color=FAMILY_COLOR[fam], ms=6, mec="white", mew=1, alpha=0.9,
                            label=None if labelled else fam)
                    labelled = True
        ax.axhline(0, color="#c3c2b7", lw=0.8)
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels([f"{t}\n({NOMINAL_Q[canonical(t)]:+d})" for t in types], fontsize=8)
        ax.set_ylabel("neutral Ef (eV)")
        ax.grid(axis="y", alpha=0.25, lw=0.6)
        ax.set_title(f"{c}: every sampled site; D = random, V = vacancy layer + mixed cations, O = Sb-Pb-Pb-Sb-Vac", loc="left", fontsize=9)
        ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
    for ext in ("png", "svg"):
        fig.savefig(out / f"DVO_中性形成能分布.{ext}", dpi=180)
    plt.close(fig)


def cmd_plan(args):
    """Phase table from the template job list (counts only)."""
    rows = list(csv.DictReader(open(Path(args.models) / "jobs.csv", encoding="utf-8-sig")))
    phase = []
    for r in rows:
        if r["stage"] == "PARENT_RELAX_REQUIRED":
            p = "1a"
        else:
            fam = r["parent_id"][0]
            site = int(r["job_id"].split("__s")[-1])
            key = canonical(r["defect_type"]) in KEY_TYPES
            p = "1b" if key and fam in "DO" else "2" if key and fam == "V" else "3"
        phase.append({**r, "phase": p})
    cnt = Counter(p["phase"] for p in phase)
    with open(Path(args.out), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["phase", "job_id", "parent_id", "defect_type", "site_id", "atom_count"])
        w.writeheader()
        w.writerows({k: p[k] for k in w.fieldnames} for p in sorted(phase, key=lambda p: (p["phase"], p["job_id"])))
    print(json.dumps(dict(sorted(cnt.items())), ensure_ascii=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("parents", "refs", "defects"):
        s = sub.add_parser(name)
        s.add_argument("--out", required=True, help="work root; jobs go to <out>/parents, <out>/refs, <out>/defects")
        s.add_argument("--kmesh", type=int, nargs=3, default=[2, 2, 2])
        s.add_argument("--shift", action="store_true", help="Monkhorst-Pack shifted mesh as in the paper (default Gamma)")
        s.add_argument("--partition", default="cpu192")
        s.add_argument("--nodes", type=int, default=1)
        s.add_argument("--ntasks", type=int, default=48)
        s.add_argument("--vasp", default="/public/home/wangchao/software/vasp.6.4.3/bin/vasp_std")
        s.add_argument("--restarts", type=int, default=3)
        s.add_argument("--no-potcar", action="store_true", help="do not write POTCAR (copy it on the cluster)")
        s.add_argument("--families", default="DVO")
        s.add_argument("--models", default=str(DVO), help="D/V/O model package (parents/, defects/, jobs.csv)")
        s.add_argument("--material", default="Pb2Sb2Te5", help="refs: which layered 9-atom POSCAR and B2Te3 to include")
        if name == "defects":
            s.add_argument("--parents", required=True, help="the <out>/parents folder with finished parent jobs")
            s.add_argument("--types", nargs="+", default=list(KEY_TYPES))
            s.add_argument("--max-sites", type=int, default=4)
            s.add_argument("--charges", help="CSV parent_id,defect_id,\"q1 q2 ...\" for a charged batch; omit for q=0")
            s.add_argument("--spin", choices=("off", "auto"), default="off", help="auto: ISPIN=2 for odd electron counts")
    s = sub.add_parser("collect")
    s.add_argument("--root", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("analyze")
    s.add_argument("--energies", required=True)
    s.add_argument("--root", help="work root used by collect (default: folder of energies.csv)")
    s.add_argument("--out", required=True)
    s.add_argument("--mu", help="JSON {condition: {Pb:, Sb:, Te:}} absolute mu in eV/atom; overrides refs")
    s.add_argument("--T-anneal", type=float, nargs="+", default=[473.0, 573.0])
    s.add_argument("--T-room", type=float, default=300.0)
    s.add_argument("--eps", type=float, default=236.0, help="isotropic static dielectric for the image-charge term")
    s.add_argument("--exclude", type=float, default=6.0, help="radius (A) around the defect excluded from alignment")
    s.add_argument("--dos", action="store_true", help="solve frozen-defect neutrality with parent static DOSCAR")
    s.add_argument("--min-carriers", type=float, default=1e17, help="|n-p| below this (cm^-3) is reported as undetermined")
    s = sub.add_parser("plan")
    s.add_argument("--out", required=True)
    s.add_argument("--models", default=str(DVO))
    args = ap.parse_args()
    if hasattr(args, "families"):
        args.families = set(args.families)
    {"parents": cmd_parents, "refs": cmd_refs, "defects": cmd_defects, "collect": cmd_collect,
     "analyze": cmd_analyze, "plan": cmd_plan}[args.cmd](args)


if __name__ == "__main__":
    main()
