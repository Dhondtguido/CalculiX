#!/usr/bin/env python3
"""
post_hinge_theta.py  (v4)

Usage:
  python3 post_hinge_theta.py <job.dat> <job.inp>

Reads ROT node IDs from the .inp (*NSET, NSET=ROT_ALL)
and then prints U3 (theta_z) and relative hinge rotations.
"""
import re, sys
from pathlib import Path

def read_nset_ids(inp_text: str, set_name: str):
    # very simple parser: finds "*NSET, NSET=<set_name>" then reads subsequent lines until next "*"
    pat = re.compile(rf"^\*NSET\s*,\s*NSET\s*=\s*{re.escape(set_name)}\s*$", re.I)
    lines = inp_text.splitlines()
    for i,l in enumerate(lines):
        if pat.match(l.strip()):
            ids=[]
            j=i+1
            while j < len(lines):
                s=lines[j].strip()
                if not s:
                    j += 1
                    continue
                if s.startswith("*"):
                    break
                parts=[p.strip() for p in s.split(",") if p.strip()]
                for p in parts:
                    if p.isdigit():
                        ids.append(int(p))
                j += 1
            if not ids:
                raise RuntimeError(f"NSET {set_name} found but empty")
            return ids
    raise RuntimeError(f"NSET {set_name} not found in inp")

def parse_last_disp_block(dat_text: str, set_name: str):
    hdr = re.compile(rf"^\s*displacements.*\b{re.escape(set_name)}\b.*\btime\b", re.I)
    lines = dat_text.splitlines()
    idxs = [i for i,l in enumerate(lines) if hdr.search(l)]
    if not idxs:
        raise RuntimeError(f"Missing displacement block for set {set_name} (analysis probably aborted before printing)")
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

def main(dat_fn: str, inp_fn: str):
    dat = Path(dat_fn).read_text(encoding="utf-8", errors="ignore")
    inp = Path(inp_fn).read_text(encoding="utf-8", errors="ignore")

    rot_ids = read_nset_ids(inp, "ROT_ALL")
    # we expect 4
    rot = parse_last_disp_block(dat, "ROT_ALL")

    th = {nid: rot.get(nid, (float("nan"),float("nan"),float("nan")))[2] for nid in rot_ids}
    print("theta_z (rad) from ROT nodes (U3):")
    for nid in rot_ids:
        print(f"  {nid}: {th[nid]: .6e}")

    # hinge pairs: use specific sets if present, else assume [Lcol,Lbem,Rcol,Rbem] order
    try:
        lcol = read_nset_ids(inp, "ROT_L_COL")[0]
        lbem = read_nset_ids(inp, "ROT_L_BEM")[0]
        rcol = read_nset_ids(inp, "ROT_R_COL")[0]
        rbem = read_nset_ids(inp, "ROT_R_BEM")[0]
    except Exception:
        if len(rot_ids) >= 4:
            lcol, lbem, rcol, rbem = rot_ids[:4]
        else:
            print("\nNot enough ROT nodes to compute hinge dtheta")
            return

    print("\nHinge relative rotations:")
    print(f"  Left : dtheta = {(th[lbem]-th[lcol]):.6e} rad   (ROT_BEM - ROT_COL)")
    print(f"  Right: dtheta = {(th[rbem]-th[rcol]):.6e} rad   (ROT_BEM - ROT_COL)")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 post_hinge_theta.py <job.dat> <job.inp>")
        raise SystemExit(2)
    main(sys.argv[1], sys.argv[2])
