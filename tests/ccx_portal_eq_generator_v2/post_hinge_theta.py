#!/usr/bin/env python3
import re, sys
from pathlib import Path

def parse_last_block(dat_text: str, set_name: str):
    hdr = re.compile(rf"^\s*displacements.*\b{re.escape(set_name)}\b.*\btime\b", re.I)
    lines = dat_text.splitlines()
    idxs = [i for i,l in enumerate(lines) if hdr.search(l)]
    if not idxs:
        raise RuntimeError(f"Missing displacement block for set {set_name}")
    i = idxs[-1] + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    data = {}
    while i < len(lines):
        l = lines[i].strip()
        if not l:
            break
        if l.lower().startswith(("displacements","forces")):
            break
        parts = l.split()
        if len(parts) >= 4 and parts[0].isdigit():
            nid = int(parts[0])
            data[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
        i += 1
    return data

def main(fn: str):
    text = Path(fn).read_text(encoding="utf-8", errors="ignore")
    rot = parse_last_block(text, "ROT_ALL")
    th = {nid: rot[nid][2] for nid in (1101,1102,1103,1104)}
    print("theta_z (rad) from ROT nodes (U3):")
    for nid in (1101,1102,1103,1104):
        print(f"  {nid}: {th[nid]: .6e}")
    print("\nHinge relative rotations:")
    print(f"  Left : dtheta = {(th[1102]-th[1101]):.6e} rad")
    print(f"  Right: dtheta = {(th[1104]-th[1103]):.6e} rad")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 post_hinge_theta.py <jobname.dat>")
        raise SystemExit(2)
    main(sys.argv[1])
