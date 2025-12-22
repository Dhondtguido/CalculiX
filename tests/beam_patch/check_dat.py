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
