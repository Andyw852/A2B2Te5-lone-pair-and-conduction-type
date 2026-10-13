#!/usr/bin/env python3
"""Mark a fresh isolated SOC project as seeded from an existing structure."""
from pathlib import Path
import hashlib
import json

cwd = Path.cwd()
poscar = cwd / "POSCAR"
potcar = cwd / "POTCAR"
if not poscar.is_file() or poscar.stat().st_size == 0:
    raise SystemExit("[ERROR] missing POSCAR")
if not potcar.is_file() or potcar.stat().st_size == 0:
    raise SystemExit("[ERROR] missing POTCAR")

out = cwd / "step0_seed"
out.mkdir(exist_ok=True)
manifest = {
    "method": "isolated dense SOC grid from verified perfect-cell geometry",
    "material": cwd.name,
    "poscar_sha256": hashlib.sha256(poscar.read_bytes()).hexdigest(),
    "potcar_sha256": hashlib.sha256(potcar.read_bytes()).hexdigest(),
}
(out / "seed_manifest.json").write_text(
    json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print("[DONE] seed manifest written; no VASP job submitted")
