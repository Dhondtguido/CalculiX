#!/usr/bin/env bash
set -euo pipefail

# Build ccx if needed (example):
#   cd /workspace/CalculiX/src && make

CCX=${CCX:-ccx}

run_case() {
  local name=$1
  echo "Running ${name}.inp"
  ${CCX} -i "${name}" > "${name}.log"
}

extract_theta() {
  local dat=$1
  awk '
    /displacements \(v\(i\),i=1\.\.ndof\)/ {found=1; next}
    found && NF>=7 {print $7; exit}
  ' "${dat}"
}

run_case test_rot_spring2_linear
linear_theta=$(extract_theta test_rot_spring2_linear.dat)
if [ -z "${linear_theta}" ]; then
  echo "Failed to read rotation from test_rot_spring2_linear.dat" >&2
  exit 1
fi
python3 - <<PY
import math
val=float("${linear_theta}")
expected=0.001
if abs(val-expected) > 1e-6:
    raise SystemExit(f"Linear rotation check failed: {val} vs {expected}")
print(f"Linear rotation OK: {val}")
PY

run_case test_rot_spring2_nonlinear
nonlinear_theta=$(extract_theta test_rot_spring2_nonlinear.dat)
if [ -z "${nonlinear_theta}" ]; then
  echo "Failed to read rotation from test_rot_spring2_nonlinear.dat" >&2
  exit 1
fi
python3 - <<PY
val=float("${nonlinear_theta}")
if val <= 0.001:
    raise SystemExit(f"Nonlinear rotation check failed: {val} <= 0.001")
print(f"Nonlinear rotation OK: {val}")
PY
