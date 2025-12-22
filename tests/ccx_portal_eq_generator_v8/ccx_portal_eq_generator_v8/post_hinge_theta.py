#!/usr/bin/env python3
"""
Postprocess portal_eq.dat:
- reads hinge_map.json (preferred) or parses portal_eq.inp
- extracts U1 of dummy nodes from *NODE PRINT blocks (ROT_DUMMIES)
- computes hinge relative rotation dtheta = U1(dummyB) - U1(dummyA)
- computes hinge moment M(dtheta) by interpolating the *SPRING NONLINEAR table (theta, M)

Outputs:
- hinge_mtheta.csv   time, dtheta_<hinge>, M_<hinge>  (time may be NaN if not found in .dat)
"""
import argparse, json, re, csv
from pathlib import Path

def parse_spring_tables(inp_text: str):
    """
    Returns dict elset -> list[(theta, M)] from:
      *SPRING, ELSET=NAME, NONLINEAR
      6, 6
      theta, M
      ...
    """
    lines = inp_text.splitlines()
    tables = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.upper().startswith("*SPRING") and "ELSET=" in line.upper() and "NONLINEAR" in line.upper():
            m = re.search(r"ELSET=([^,\s]+)", line, flags=re.IGNORECASE)
            if not m:
                i += 1
                continue
            elset = m.group(1)
            i += 1
            # skip dof line
            if i < len(lines) and re.match(r"^\s*\d+\s*,\s*\d+\s*$", lines[i]):
                i += 1
            pts=[]
            while i < len(lines):
                s = lines[i].strip()
                if not s or s.startswith("**"):
                    i += 1
                    continue
                if s.startswith("*"):
                    break
                # theta, M
                parts = [p.strip() for p in s.split(",")]
                if len(parts) >= 2:
                    try:
                        th = float(parts[0])
                        M = float(parts[1])
                        pts.append((th,M))
                    except ValueError:
                        pass
                i += 1
            if pts:
                pts = sorted(pts, key=lambda x: x[0])
                tables[elset] = pts
            continue
        i += 1
    return tables

def interp_piecewise(pts, x):
    """Linear interpolation/extrapolation for sorted pts [(x_i,y_i)]."""
    if not pts:
        return float("nan")
    if x <= pts[0][0]:
        # extrapolate with first segment slope
        x0,y0=pts[0]; x1,y1=pts[1] if len(pts)>1 else (x0+1.0,y0)
        m = (y1-y0)/(x1-x0) if x1!=x0 else 0.0
        return y0 + m*(x-x0)
    if x >= pts[-1][0]:
        x0,y0=pts[-2] if len(pts)>1 else (pts[-1][0]-1.0, pts[-1][1])
        x1,y1=pts[-1]
        m = (y1-y0)/(x1-x0) if x1!=x0 else 0.0
        return y1 + m*(x-x1)
    # inside
    for (x0,y0),(x1,y1) in zip(pts[:-1], pts[1:]):
        if x0 <= x <= x1:
            if x1==x0:
                return y0
            a = (x - x0)/(x1 - x0)
            return y0 + a*(y1-y0)
    return pts[-1][1]

def load_hinge_map(workdir: Path, inp_text: str):
    hm = workdir/"hinge_map.json"
    if hm.exists():
        return json.loads(hm.read_text())
    # fallback: try to parse (very minimal)
    raise RuntimeError("hinge_map.json not found; run generator in the same folder as the input")

def parse_time_from_context(lines, idx, lookback=40):
    # Look backward for something like "time=" or "actual total time="
    for j in range(max(0, idx-lookback), idx)[::-1]:
        s = lines[j]
        m = re.search(r"(?:time\s*=|actual total time=)\s*([0-9Ee+\-\.]+)", s, flags=re.IGNORECASE)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
    return float("nan")

def parse_node_print_blocks(dat_text: str, node_ids):
    """
    Returns list of (time, {node: u1}).
    We scan all "displacements" blocks and collect u1 for nodes in node_ids.
    """
    lines = dat_text.splitlines()
    node_ids = set(node_ids)
    records = []
    i = 0
    while i < len(lines):
        if "displacements" in lines[i].lower():
            t = parse_time_from_context(lines, i)
            # Advance to first numeric node line
            i += 1
            values = {}
            while i < len(lines):
                s = lines[i].strip()
                if not s:
                    # blank ends block
                    break
                if s.startswith("*") or s.lower().startswith("step") or "displacements" in s.lower():
                    break
                # node line: nid u1 u2 u3 [maybe more]
                m = re.match(r"^(\d+)\s+([0-9Ee+\-\.]+)\s+([0-9Ee+\-\.]+)\s+([0-9Ee+\-\.]+)", s)
                if m:
                    nid = int(m.group(1))
                    if nid in node_ids:
                        u1 = float(m.group(2))
                        values[nid] = u1
                i += 1
            # keep only if we captured at least one target node
            if values:
                records.append((t, values))
        i += 1
    return records

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dat", type=Path)
    ap.add_argument("inp", type=Path)
    ap.add_argument("--out", type=Path, default=Path("hinge_mtheta.csv"))
    args = ap.parse_args()

    inp_text = args.inp.read_text()
    tables = parse_spring_tables(inp_text)
    hm = load_hinge_map(args.inp.parent, inp_text)

    hinges = hm["hinges"]
    # collect dummy nodes
    dummy_nodes = []
    for h in hinges:
        dummy_nodes += [h["dummyA"], h["dummyB"]]
    dummy_nodes = sorted(set(dummy_nodes))

    dat_text = args.dat.read_text(errors="ignore")
    recs = parse_node_print_blocks(dat_text, dummy_nodes)
    if not recs:
        raise RuntimeError("No displacement blocks found in .dat (did the run abort early, or NODE PRINT not requested?)")

    # Build time series: we only keep blocks where we have *all* dummy nodes
    rows = []
    for t,vals in recs:
        if not all(n in vals for n in dummy_nodes):
            continue
        row = {"t": t}
        for h in hinges:
            thA = vals[h["dummyA"]]
            thB = vals[h["dummyB"]]
            dth = thB - thA
            row[f"dtheta_{h['name']}"] = dth
            pts = tables.get(h["name"])
            row[f"M_{h['name']}"] = interp_piecewise(pts, dth) if pts else float("nan")
        rows.append(row)

    if not rows:
        raise RuntimeError("Found displacement blocks, but none contained all dummy nodes. Try FREQUENCY=1 for ROT_DUMMIES.")

    # write CSV
    cols = ["t"]
    for h in hinges:
        cols += [f"dtheta_{h['name']}", f"M_{h['name']}"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    # print quick summary
    for h in hinges:
        dths = [r[f"dtheta_{h['name']}"] for r in rows]
        Ms   = [r[f"M_{h['name']}"] for r in rows]
        print(f"{h['name']}: dtheta range [{min(dths):.3e}, {max(dths):.3e}] rad, "
              f"M range [{min(Ms):.3e}, {max(Ms):.3e}] N*m")
    print(f"Wrote {args.out} with {len(rows)} samples")

if __name__ == "__main__":
    main()
