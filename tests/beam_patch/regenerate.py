#!/usr/bin/env python3
from pathlib import Path

BASE = Path(__file__).resolve().parent

CANT = """\
** Cantilever beam (B32R), 1 element, end load in +X
*NODE,NSET=NALL
1, 0, 0,   0
2, 0, 0,  50
3, 0, 0, 100
*ELEMENT,TYPE=B32R,ELSET=EALL
1,1,2,3
*NSET,NSET=TIP
1
*NSET,NSET=FIX
3
*BOUNDARY
3,1,6
*MATERIAL,NAME=MAT
*ELASTIC
1.E7, .3
*BEAM SECTION,ELSET=EALL,MATERIAL=MAT,SECTION=RECT
2.,2.
1.,0.,0.
*STEP
*STATIC
*CLOAD
1,1,1.
*NODE PRINT,NSET=TIP
U
*END STEP
"""

SS = """\
** Simply supported beam (B32R), 1 element, midspan load in +X
*NODE,NSET=NALL
1, 0, 0,   0
2, 0, 0,  50
3, 0, 0, 100
*ELEMENT,TYPE=B32R,ELSET=EALL
1,1,2,3
*NSET,NSET=MID
2
*BOUNDARY
1,1,3
1,6,6
3,1,2
*MATERIAL,NAME=MAT
*ELASTIC
1.E7, .3
*BEAM SECTION,ELSET=EALL,MATERIAL=MAT,SECTION=RECT
2.,2.
1.,0.,0.
*STEP
*STATIC
*CLOAD
2,1,1.
*NODE PRINT,NSET=MID
U
*END STEP
"""

CHECK = """\
#!/usr/bin/env python3
import argparse, math, re
from pathlib import Path

def parse_last_displacements_block(text: str):
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if "displacements (vx,vy,vz)" in line.lower():
            start = i
    if start is None:
        raise RuntimeError("No 'displacements (vx,vy,vz)' block found in .dat")

    rows = {}
    for j in range(start + 1, len(lines)):
        s = lines[j].strip()
        if not s:
            continue
        toks = s.split()
        if len(toks) == 4:
            try:
                node = int(toks[0])
                vx = float(toks[1].replace("D", "E"))
                vy = float(toks[2].replace("D", "E"))
                vz = float(toks[3].replace("D", "E"))
            except ValueError:
                continue
            rows[node] = {"U1": vx, "U2": vy, "U3": vz}
    if not rows:
        raise RuntimeError("Found displacements header, but no numeric rows parsed")
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat", type=Path)
    ap.add_argument("--node", type=int, required=True)
    ap.add_argument("--comp", required=True)
    ap.add_argument("--expected", type=float, required=True)
    ap.add_argument("--rtol", type=float, default=0.10)
    ap.add_argument("--atol", type=float, default=0.0)
    args = ap.parse_args()

    rows = parse_last_displacements_block(args.dat.read_text(errors="ignore"))
    comp = args.comp.upper()
    got = rows[args.node][comp]
    exp = args.expected
    tol = max(args.atol, args.rtol * abs(exp))
    ok = math.isfinite(got) and abs(got - exp) <= tol
    print(f"[check_dat] {args.dat.name}: node={args.node} {comp} got={got:.10g} exp={exp:.10g} tol={tol:.10g} -> {'OK' if ok else 'FAIL'}")
    raise SystemExit(0 if ok else 2)

if __name__ == "__main__":
    main()
"""

RUN = """\
#!/usr/bin/env bash
set -euo pipefail
CCX="${CCX:-ccx}"
PY="${PY:-python3}"
cd "$(dirname "$0")"

echo "[beam_patch] Test 1: cantilever"
$CCX -i cantilever_b32r_unitload
$PY check_dat.py cantilever_b32r_unitload.dat --node 1 --comp U1 --expected 0.02250649 --rtol 0.10

echo "[beam_patch] Test 2: simply supported"
$CCX -i simply_supported_b32r_midload
$PY check_dat.py simply_supported_b32r_midload.dat --node 2 --comp U1 --expected 0.00113403 --rtol 0.10

echo "[beam_patch] ALL OK"
"""

def main():
    BASE.mkdir(parents=True, exist_ok=True)
    (BASE/"cantilever_b32r_unitload.inp").write_text(CANT)
    (BASE/"simply_supported_b32r_midload.inp").write_text(SS)
    (BASE/"check_dat.py").write_text(CHECK)
    (BASE/"run_tests.sh").write_text(RUN)
    # permissions
    import os
    os.chmod(BASE/"check_dat.py", 0o755)
    os.chmod(BASE/"run_tests.sh", 0o755)
    print(f"[regenerate] Wrote tests into {BASE}")

if __name__ == "__main__":
    main()
