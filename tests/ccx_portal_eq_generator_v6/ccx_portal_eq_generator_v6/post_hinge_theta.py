\
#!/usr/bin/env python3
"""
post_hinge_theta.py (v6)

python3 post_hinge_theta.py <job.dat> <job.inp>
Prints theta_z from ROT nodes (U3) and hinge relative rotations.
"""
import re, sys
from pathlib import Path
from typing import Dict, List, Tuple

def read_nset_ids(inp_text: str, set_name: str) -> List[int]:
    pat = re.compile(rf"^\*NSET\s*,\s*NSET\s*=\s*{re.escape(set_name)}\s*$", re.I)
    lines = inp_text.splitlines()
    for i, l in enumerate(lines):
        if pat.match(l.strip()):
            ids: List[int] = []
            j = i + 1
            while j < len(lines):
                s = lines[j].strip()
                if not s:
                    j += 1
                    continue
                if s.startswith("*"):
                    break
                for p in [p.strip() for p in s.split(",") if p.strip()]:
                    if p.isdigit():
                        ids.append(int(p))
                j += 1
            if not ids:
                raise RuntimeError(f"NSET {set_name} found but empty")
            return ids
    raise RuntimeError(f"NSET {set_name} not found")

def parse_last_disp_block(dat_text: str, set_name: str) -> Dict[int, Tuple[float, float, float]]:
    hdr = re.compile(rf"^\s*displacements.*\b{re.escape(set_name)}\b.*\btime\b", re.I)
    lines = dat_text.splitlines()
    idxs = [i for i, l in enumerate(lines) if hdr.search(l)]
    if not idxs:
        raise RuntimeError(f"Missing displacement block for set {set_name}")
    i = idxs[-1] + 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    data: Dict[int, Tuple[float, float, float]] = {}
    while i < len(lines):
        l = lines[i].strip()
        if not l:
            break
        if l.lower().startswith(("displacements", "forces", "stresses", "strains")):
            break
        parts = l.split()
        if len(parts) >= 4 and parts[0].isdigit():
            data[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
        i += 1
    return data

def main(dat_fn: str, inp_fn: str) -> None:
    dat = Path(dat_fn).read_text(encoding="utf-8", errors="ignore")
    inp = Path(inp_fn).read_text(encoding="utf-8", errors="ignore")

    rot_all = read_nset_ids(inp, "ROT_ALL")
    disp = parse_last_disp_block(dat, "ROT_ALL")

    print("theta_z (rad) from ROT nodes (U3):")
    th = {}
    for nid in rot_all:
        th[nid] = disp.get(nid, (float("nan"),)*3)[2]
        print(f"  {nid}: {th[nid]: .6e}")

    lcol = read_nset_ids(inp, "ROT_L_COL")[0]
    lbem = read_nset_ids(inp, "ROT_L_BEM")[0]
    rcol = read_nset_ids(inp, "ROT_R_COL")[0]
    rbem = read_nset_ids(inp, "ROT_R_BEM")[0]

    print("\nHinge relative rotations:")
    print(f"  Left : dtheta = {(th[lbem] - th[lcol]): .6e} rad")
    print(f"  Right: dtheta = {(th[rbem] - th[rcol]): .6e} rad")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 post_hinge_theta.py <job.dat> <job.inp>")
        raise SystemExit(2)
    main(sys.argv[1], sys.argv[2])
