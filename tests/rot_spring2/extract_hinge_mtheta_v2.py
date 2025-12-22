#!/usr/bin/env python3
"""
Extract hinge relative rotations (SPRING2 on ROT-node U3) and compute hinge moments
from a piecewise-linear M-θ table.

Usage:
  python3 extract_hinge_mtheta_v2.py portal_frame_solids_hinges_v6.dat

Outputs:
  - top drift (ux at REF_COL node 1003 / H)
  - left hinge  dtheta = U3(1102) - U3(1101)
  - right hinge dtheta = U3(1104) - U3(1103)
  - hinge moments Mz evaluated on the piecewise-linear curve (symmetric)

Notes:
  - ROT nodes: U3 is theta_z [rad]
  - The parser is tolerant to blank lines right after the header (CalculiX does that).
"""

import re
import sys
import pathlib
from typing import Dict, Tuple, List


def parse_last_block(dat_text: str, kind: str, set_name: str) -> Dict[int, Tuple[float, float, float]]:
    """Parse the LAST displacement/force block for a given set."""
    # Flexible header: tolerate "for set" / "for node set" etc.
    hdr = re.compile(
        rf"^\s*{re.escape(kind)}\s*\(.*\)\s*for\s+.*\b{re.escape(set_name)}\b.*\btime\b.*$",
        re.IGNORECASE,
    )

    lines = dat_text.splitlines()
    idxs = [i for i, l in enumerate(lines) if hdr.search(l)]
    if not idxs:
        raise RuntimeError(f"Could not find '{kind}' block for set {set_name} in .dat")

    i = idxs[-1] + 1

    # Skip initial blank lines after header (CalculiX prints one)
    while i < len(lines) and not lines[i].strip():
        i += 1

    data: Dict[int, Tuple[float, float, float]] = {}
    while i < len(lines):
        l = lines[i]
        if not l.strip():
            break
        # Stop if we ran into another header
        if l.lower().lstrip().startswith(("displacements", "forces")):
            break

        parts = l.split()
        if len(parts) >= 4 and parts[0].lstrip("+-").isdigit():
            nid = int(parts[0])
            v1, v2, v3 = float(parts[1]), float(parts[2]), float(parts[3])
            data[nid] = (v1, v2, v3)
        i += 1

    return data


def piecewise_linear(u: float, pts: List[Tuple[float, float]]) -> float:
    """Evaluate piecewise-linear f(u) given sorted points."""
    pts = sorted(pts, key=lambda x: x[0])
    if len(pts) < 2:
        raise ValueError("Need at least two points for a piecewise-linear curve")

    if u <= pts[0][0]:
        (u0, f0), (u1, f1) = pts[0], pts[1]
        k = (f1 - f0) / (u1 - u0)
        return f0 + k * (u - u0)

    if u >= pts[-1][0]:
        (u0, f0), (u1, f1) = pts[-2], pts[-1]
        k = (f1 - f0) / (u1 - u0)
        return f1 + k * (u - u1)

    for (u0, f0), (u1, f1) in zip(pts[:-1], pts[1:]):
        if u0 <= u <= u1:
            k = (f1 - f0) / (u1 - u0)
            return f0 + k * (u - u0)

    raise RuntimeError("Unexpected: u not bracketed")


def main(dat_fn: str) -> None:
    dat = pathlib.Path(dat_fn).read_text(encoding="utf-8", errors="ignore")

    rot = parse_last_block(dat, "displacements", "ROT_ALL")
    ref = parse_last_block(dat, "displacements", "REF_COL")

    th1101 = rot[1101][2]
    th1102 = rot[1102][2]
    th1103 = rot[1103][2]
    th1104 = rot[1104][2]

    dL = th1102 - th1101
    dR = th1104 - th1103

    ux_top = ref[1003][0]
    H = 3.0
    drift = ux_top / H

    # Symmetric M-θ table used in v6 (N*m)
    pts = [
        (-0.020, -2.2e5),
        (-0.002, -2.0e5),
        ( 0.000,  0.0),
        ( 0.002,  2.0e5),
        ( 0.020,  2.2e5),
    ]

    ML = piecewise_linear(dL, pts)
    MR = piecewise_linear(dR, pts)

    print(f"Top drift: ux(top)= {ux_top:.6e} m  -> drift ratio = {drift:.6e} ({drift*100:.3f} %)")
    print(f"Left hinge : dtheta = {dL:.6e} rad  -> Mz = {ML/1000:.3f} kN*m")
    print(f"Right hinge: dtheta = {dR:.6e} rad  -> Mz = {MR/1000:.3f} kN*m")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 extract_hinge_mtheta_v2.py <jobname.dat>")
        raise SystemExit(2)
    main(sys.argv[1])
