#!/usr/bin/env python3
"""Run the existing AutoZT BT2 CRTA implementation on the dense SOC grid."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gen_step11_boltztrap as bt  # noqa: E402

GAPS = {
    "Pb2Bi2Te5": 0.1141,
    "Pb2Sb2Te5": 0.0654,
    "Sn2Bi2Te5": 0.1366,
    "Sn2Sb2Te5": 0.1352,
}

cwd = Path.cwd()
parser = argparse.ArgumentParser()
parser.add_argument("--material", required=True)
mat = parser.parse_args().material
if mat not in GAPS:
    raise SystemExit("[ERROR] unsupported material: %s" % mat)
uniform = cwd / "step3_uniform"
if not (uniform / "vasprun.xml").is_file():
    raise SystemExit("[ERROR] dense step3_uniform/vasprun.xml is missing")
if not (uniform / "KPOINTS").is_file():
    raise SystemExit("[ERROR] dense step3_uniform/KPOINTS is missing")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# strict_done_marker requires relative upstream fingerprints so AutoZT can
# invalidate this post-processing result when the dense-grid input changes.
input_files = {
    "step3_uniform/vasprun.xml": sha256(uniform / "vasprun.xml"),
    "step3_uniform/KPOINTS": sha256(uniform / "KPOINTS"),
}

# The paper-specific DPT merge is intentionally disabled here. This step is
# the SOC BT2 CRTA cross-check; AMSET remains the separate scattering result.
bt.PAPER_SCAN = False
result = bt.run_boltztrap_crta(
    uniform,
    temperatures=[300, 500, 700],
    mu_window_ev=2.0,
    # The dense SOC vasprun.xml files are multi-GB. Keep the first CRTA
    # cross-check below the AutoZT run:gen timeout; convergence can be tested
    # separately after the finite result is available.
    mu_npts=500,
    dos_npts=1500,
    kmesh_mult=2,
    nelect=None,
    is_2d=False,
    scissor_gap_eV=GAPS[mat],
)
result.update({
    "material": mat,
    "dim": "3d",
    "soc": True,
    "status": "done",
    "input_files": input_files,
    "source_grid": str(uniform / "KPOINTS"),
    "target_gap_eV": GAPS[mat],
    "protocol": "BoltzTraP2 CRTA, SOC dense uniform grid, same gap target as AMSET",
})
out = cwd / "step8.1_boltztrap"
out.mkdir(exist_ok=True)
(out / "boltztrap_crta.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
    encoding="utf-8",
)
(out / "boltztrap_summary.txt").write_text(
    "material=%s\nmesh=%s\ntarget_gap_eV=%.4f\n" %
    (mat, result.get("s3_grid", {}).get("mesh", "see KPOINTS"), GAPS[mat]),
    encoding="utf-8",
)
print("[DONE] %s: dense SOC BT2 CRTA written" % mat)
