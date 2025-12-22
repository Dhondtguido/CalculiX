\
#!/usr/bin/env bash
set -euo pipefail
CCX="${CCX:-../../src/CalculiX}"
JOB="${JOB:-portal_eq}"

python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out "${JOB}.inp" --gravity
echo "[run_all] ENERGY PRINT sanity:"
grep -n "ENERGY PRINT" "${JOB}.inp" || true

"${CCX}" -i "${JOB}"
python3 post_hinge_theta.py "${JOB}.dat" "${JOB}.inp"
