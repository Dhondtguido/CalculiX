# ''#!/usr/bin/env python3
# ´Extract hinge relative rotations and compute hinge moments from the input spring table.

# Usage:
#   python3 extract_hinge_mtheta.py portal_frame_solids_hinges_v6.dat

# It prints:
#   - top drift (from REF_COL node 1003)
#   - hinge relative rotations:
#         left  = U3(1102) - U3(1101)
#         right = U3(1104) - U3(1103)
#   - hinge moments evaluated on the piecewise-linear M-θ curve (symmetric)}}''

import re, sys

def parse_last_displacements(dat_text, set_name):
    # Finds last block "displacements (vx,vy,vz) for set <set_name>"
    pattern = re.compile(rf"displacements\s*\(vx,vy,vz\)\s*for\s*set\s+{re.escape(set_name)}\s+and\s+time.*", re.IGNORECASE)
    lines = dat_text.splitlines()
    idxs=[i for i,l in enumerate(lines) if pattern.search(l)]
    if not idxs:
        raise RuntimeError(f"Could not find displacement block for set {set_name}")
    i0=idxs[-1]+1
    data={}
    for l in lines[i0:]:
        if not l.strip():
            break
        parts=l.split()
        if len(parts)>=4 and parts[0].isdigit():
            nid=int(parts[0]); ux=float(parts[1]); uy=float(parts[2]); uz=float(parts[3])
            data[nid]=(ux,uy,uz)
    return data

def piecewise_force(u, pts):
    # pts: list of (u,f) sorted
    if u<=pts[0][0]:
        # extrapolate first segment
        (u0,f0),(u1,f1)=pts[0],pts[1]
        k=(f1-f0)/(u1-u0)
        return f0 + k*(u-u0)
    if u>=pts[-1][0]:
        (u0,f0),(u1,f1)=pts[-2],pts[-1]
        k=(f1-f0)/(u1-u0)
        return f1 + k*(u-u1)
    for (u0,f0),(u1,f1) in zip(pts[:-1],pts[1:]):
        if u0<=u<=u1:
            k=(f1-f0)/(u1-u0)
            return f0 + k*(u-u0)
    raise RuntimeError("unreachable")

def main(fn):
    dat=open(fn,"r",encoding="utf-8",errors="ignore").read()
    rot=parse_last_displacements(dat,"ROT_ALL")
    ref=parse_last_displacements(dat,"REF_COL")

    th1101=rot[1101][2]; th1102=rot[1102][2]
    th1103=rot[1103][2]; th1104=rot[1104][2]
    dL=th1102-th1101
    dR=th1104-th1103

    ux_top=ref[1003][0]
    H=3.0
    drift=ux_top/H

    pts=[(-0.020,-2.2e5),(-0.002,-2.0e5),(0.0,0.0),(0.002,2.0e5),(0.020,2.2e5)]
    ML=piecewise_force(dL,pts)
    MR=piecewise_force(dR,pts)

    print(f"Top drift: ux(top)= {ux_top:.6e} m  -> drift ratio = {drift:.6e} ({drift*100:.3f} %)")
    print(f"Left hinge : dtheta = {dL:.6e} rad  -> Mz = {ML/1000:.3f} kN*m")
    print(f"Right hinge: dtheta = {dR:.6e} rad  -> Mz = {MR/1000:.3f} kN*m")

if __name__ == '__main__':
    if len(sys.argv)<2:
        print("Usage: python3 extract_hinge_mtheta.py <jobname.dat>")
        sys.exit(2)
    main(sys.argv[1])
