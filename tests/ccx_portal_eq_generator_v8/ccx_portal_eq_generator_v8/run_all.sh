#!/usr/bin/env bash
set -euo pipefail

CCX_BIN="${CCX:-../../src/CalculiX}"

echo "[run_all] Using CCX=${CCX_BIN}"

python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out portal_eq.inp --gravity

if [[ ! -x "${CCX_BIN}" ]]; then
  echo "[run_all] ERROR: CCX executable not found/executable at: ${CCX_BIN}"
  echo "         Set env var CCX to your CalculiX binary, e.g.:"
  echo "         CCX=../../../src/CalculiX bash run_all.sh"
  exit 2
fi

"${CCX_BIN}" -i portal_eq

python3 post_hinge_theta.py portal_eq.dat portal_eq.inp
