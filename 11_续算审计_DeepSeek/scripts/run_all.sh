#!/bin/bash
# Reproduce the read-only audit in 11_续算审计_DeepSeek.
# No original A2B2Te5_SOC_dense file is written or modified.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ../csv ../json ../logs
python3 audit_inventory.py 2>&1 | tee ../logs/audit_inventory.log
python3 audit_mass.py 2>&1 | tee ../logs/audit_mass.log
python3 audit_s_projection.py 2>&1 | tee ../logs/audit_s_projection.log
python3 audit_bt2.py 2>&1 | tee ../logs/audit_bt2.log
python3 audit_dvo.py 2>&1 | tee ../logs/audit_dvo.log
python3 make_key_summary.py 2>&1 | tee ../logs/make_key_summary.log
python3 make_manifest.py 2>&1 | tee ../logs/make_manifest.log
echo "[DONE] audit outputs in $(cd .. && pwd)"
