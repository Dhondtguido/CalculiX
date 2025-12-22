\
#!/usr/bin/env python3
"""
post_hinge_theta.py v9

Usage:
  python3 post_hinge_theta.py portal_eq.dat portal_eq.inp

Reads last printed displacements for:
- ROT_ALL (rigid-body ROT nodes): U3 == theta_z (rad)
- REF_TOP: top joint reference node displacement (U1..U3)

Outputs:
- theta_z at each ROT node
- dtheta at left/right hinge
- top drift proxy: U1 at REF_TOP
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
        raise RuntimeError(f"Missing displacement block for set {set_name} (analysis may have aborted before printing)")
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
    rot = parse_last_disp_block(dat, "ROT_ALL")
    ref_top = read_nset_ids(inp, "REF_TOP")
    ref = parse_last_disp_block(dat, "REF_TOP")

    th = {}
    print("theta_z (rad) from ROT nodes (U3):")
    for nid in rot_all:
        th[nid] = rot[nid][2]
        print(f"  {nid}: {th[nid]: .6e}")

    # Identify left/right by matching set names in input
    lcol = read_nset_ids(inp, "ROT_L_COL")[0]
    lbem = read_nset_ids(inp, "ROT_L_BEM")[0]
    rcol = read_nset_ids(inp, "ROT_R_COL")[0]
    rbem = read_nset_ids(inp, "ROT_R_BEM")[0]

    print("\nHinge relative rotations:")
    print(f"  Left : dtheta = {(th[lbem] - th[lcol]): .6e} rad")
    print(f"  Right: dtheta = {(th[rbem] - th[rcol]): .6e} rad")

    rt = ref_top[0]
    u = ref[rt]
    print("\nTop REF node displacement (REF_TOP):")
    print(f"  node {rt}: U1={u[0]: .6e}  U2={u[1]: .6e}  U3={u[2]: .6e}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 post_hinge_theta.py <job.dat> <job.inp>")
        raise SystemExit(2)
    main(sys.argv[1], sys.argv[2])
