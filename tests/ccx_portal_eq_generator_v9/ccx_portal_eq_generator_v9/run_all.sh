\
#!/usr/bin/env bash
set -euo pipefail

JOB="${JOB:-portal_eq}"

# CCX autodetect (override with CCX=/path/to/CalculiX)
if [[ -z "${CCX:-}" ]]; then
  for cand in ../../src/CalculiX ../../../src/CalculiX ../../../../src/CalculiX ../src/CalculiX; do
    if [[ -x "$cand" ]]; then
      CCX="$cand"
      break
    fi
  done
fi
: "${CCX:?Set CCX=/path/to/CalculiX (could not auto-detect)}"

echo "[run_all] Using CCX=${CCX}"

python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out "${JOB}.inp" --gravity

echo "[run_all] ENER request sanity (should exist in STEP1 if --gravity):"
grep -n "^\*EL PRINT" "${JOB}.inp" || true
grep -n "^ENER$" "${JOB}.inp" || true

"${CCX}" -i "${JOB}"
python3 post_hinge_theta.py "${JOB}.dat" "${JOB}.inp"
