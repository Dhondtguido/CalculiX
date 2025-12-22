#!/usr/bin/env python3
"""
Generate a planar portal frame (B31 beams) for CalculiX with rotational hinges via SPRING2 on DOF 6,
and a horizontal earthquake acceleration history applied as body-acceleration (*DLOAD, GRAV, AMPLITUDE=...).

Key trick:
- CalculiX internally expands B31 beams into solids and only creates rotational "knot" DOFs if they appear in an SPC/MPC
  (or in nodal moments). Therefore, to use SPRING2 on DOF 6 we *activate* UR3 by referencing DOF 6 in a harmless MPC:
      UR3(node) == U1(dummy_node)
  We then output U1(dummy_node) which is robustly printable as a translation DOF.

Units: m, N, kg, s.
"""
import argparse, csv, json
from pathlib import Path

def read_eq_csv(path: Path):
    """CSV with header: t, ax  (t in s, ax in m/s^2)."""
    t=[]
    ax=[]
    with open(path, newline="") as f:
        r=csv.DictReader(f)
        if "t" not in r.fieldnames or "ax" not in r.fieldnames:
            raise RuntimeError(f"{path} must have columns t,ax (got {r.fieldnames})")
        for row in r:
            t.append(float(row["t"]))
            ax.append(float(row["ax"]))
    if len(t) < 2:
        raise RuntimeError("earthquake CSV must have at least 2 rows")
    return t, ax

def fmt(x: float) -> str:
    return f"{x:.8e}"

def write_amplitude(lines, name, t, a):
    lines.append(f"*AMPLITUDE, NAME={name}")
    pairs=[f"{ti:.8e}, {ai:.8e}" for ti,ai in zip(t,a)]
    for i in range(0,len(pairs),4):
        lines.append(", ".join(pairs[i:i+4]))

def write_spring_table(out, elset: str, My: float, ky: float, kpost: float):
    """
    Bilinear symmetric M-θ:
      M = ky*θ for |θ|<=θy, then M = sign(θ)*My + kpost*(|θ|-θy)
    Provided as 5 points (covers +/-).
    """
    th_y = My/ky
    out.append(f"*SPRING, ELSET={elset}, NONLINEAR")
    out.append("6, 6")
    pts = [
        (-10*th_y, -My - kpost*(10*th_y - th_y)),
        (-th_y,    -My),
        (0.0,       0.0),
        (th_y,      My),
        (10*th_y,   My + kpost*(10*th_y - th_y)),
    ]
    for th,M in pts:
        out.append(f"{fmt(th)}, {fmt(M)}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eq", type=Path, required=True, help="CSV with columns t,ax (m/s^2)")
    ap.add_argument("--out", type=Path, default=Path("portal_eq.inp"))
    ap.add_argument("--L", type=float, default=4.0, help="span [m]")
    ap.add_argument("--H", type=float, default=3.0, help="height [m]")
    ap.add_argument("--ncol", type=int, default=10, help="segments per column")
    ap.add_argument("--nbeam", type=int, default=12, help="segments along beam")
    ap.add_argument("--E", type=float, default=210e9)
    ap.add_argument("--nu", type=float, default=0.3)
    ap.add_argument("--rho", type=float, default=7850.0)
    ap.add_argument("--bcol", type=float, default=0.30, help="column rect width [m]")
    ap.add_argument("--hcol", type=float, default=0.30, help="column rect height [m]")
    ap.add_argument("--bbeam", type=float, default=0.25, help="beam rect width [m]")
    ap.add_argument("--hbeam", type=float, default=0.40, help="beam rect height [m]")
    ap.add_argument("--dt", type=float, default=0.002, help="dynamic step dt [s]")
    ap.add_argument("--T", type=float, default=5.0, help="dynamic duration [s]")
    ap.add_argument("--gravity", action="store_true", help="include gravity")
    ap.add_argument("--My", type=float, default=2.0e5, help="hinge yield moment [N*m]")
    ap.add_argument("--ky", type=float, default=1.0e8, help="hinge initial stiffness [N*m/rad]")
    ap.add_argument("--kpost", type=float, default=1.0e6, help="hinge post-yield stiffness [N*m/rad]")
    args = ap.parse_args()

    t,ax = read_eq_csv(args.eq)

    L,H = args.L, args.H
    ncol, nbeam = args.ncol, args.nbeam

    # --- Nodes (planar, z=0) ---
    nodes = {}
    nid = 1
    left_col=[]
    for i in range(ncol+1):
        y = H*i/ncol
        nodes[nid]=(0.0,y,0.0)
        left_col.append(nid); nid += 1

    right_col=[]
    for i in range(ncol+1):
        y = H*i/ncol
        nodes[nid]=(L,y,0.0)
        right_col.append(nid); nid += 1

    beam_mid=[]
    for i in range(1,nbeam):
        x = L*i/nbeam
        nodes[nid]=(x,H,0.0)
        beam_mid.append(nid); nid += 1

    # joints on column side
    n_tl = left_col[-1]
    n_tr = right_col[-1]

    # duplicate joint nodes on beam side
    n_tl_d = nid; nodes[n_tl_d]=nodes[n_tl]; nid += 1
    n_tr_d = nid; nodes[n_tr_d]=nodes[n_tr]; nid += 1

    # dummy nodes for UR3 activation
    dummy={}
    def add_dummy(for_node: int) -> int:
        nonlocal nid
        d=nid
        nodes[d]=nodes[for_node]
        nid += 1
        dummy[for_node]=d
        return d

    d_tl   = add_dummy(n_tl)
    d_tl_d = add_dummy(n_tl_d)
    d_tr   = add_dummy(n_tr)
    d_tr_d = add_dummy(n_tr_d)

    beam_nodes=[n_tl_d] + beam_mid + [n_tr_d]

    # --- Elements ---
    elems=[]
    eid=1
    col_eids=[]
    for a,b in zip(left_col[:-1], left_col[1:]):
        elems.append((eid,"B31",a,b,"E_COL")); col_eids.append(eid); eid += 1
    for a,b in zip(right_col[:-1], right_col[1:]):
        elems.append((eid,"B31",a,b,"E_COL")); col_eids.append(eid); eid += 1

    beam_eids=[]
    for a,b in zip(beam_nodes[:-1], beam_nodes[1:]):
        elems.append((eid,"B31",a,b,"E_BEAM")); beam_eids.append(eid); eid += 1

    hinges=[
        ("H_TL", n_tl, n_tl_d, d_tl, d_tl_d),
        ("H_TR", n_tr, n_tr_d, d_tr, d_tr_d),
    ]
    hinge_eids=[]
    for name,na,nb,_,_ in hinges:
        elems.append((eid,"SPRING2",na,nb,name)); hinge_eids.append(eid); eid += 1

    # --- Sets ---
    nset_fix=[left_col[0], right_col[0]]
    nset_top=[n_tl_d, n_tr_d]
    nset_dummies=[d_tl, d_tl_d, d_tr, d_tr_d]

    out=[]
    out.append("*HEADING")
    out.append("Portal frame B31 + rotational hinges (SPRING2 DOF6) + earthquake body-acceleration")
    out.append("** Units: m, N, kg, s. Rotations: rad. Moments: N*m.")
    out.append("")

    out.append("*NODE")
    for k,(x,y,z) in sorted(nodes.items()):
        out.append(f"{k}, {fmt(x)}, {fmt(y)}, {fmt(z)}")

    out.append("*NSET, NSET=FIXED_BASE")
    out.append(", ".join(str(n) for n in nset_fix))
    out.append("*NSET, NSET=TOP_NODES")
    out.append(", ".join(str(n) for n in nset_top))
    out.append("*NSET, NSET=ROT_DUMMIES")
    out.append(", ".join(str(n) for n in nset_dummies))

    # Elements
    out.append("*ELEMENT, TYPE=B31, ELSET=E_COL")
    for (eid2,et,a,b,tag) in elems:
        if eid2 in col_eids:
            out.append(f"{eid2}, {a}, {b}")

    out.append("*ELEMENT, TYPE=B31, ELSET=E_BEAM")
    for (eid2,et,a,b,tag) in elems:
        if eid2 in beam_eids:
            out.append(f"{eid2}, {a}, {b}")

    for (name,na,nb,_,_),eid_h in zip(hinges, hinge_eids):
        out.append(f"*ELEMENT, TYPE=SPRING2, ELSET={name}")
        out.append(f"{eid_h}, {na}, {nb}")

    # Material + sections
    out.append("*MATERIAL, NAME=STEEL")
    out.append("*ELASTIC")
    out.append(f"{fmt(args.E)}, {fmt(args.nu)}")
    out.append("*DENSITY")
    out.append(f"{fmt(args.rho)}")

    out.append("*BEAM SECTION, MATERIAL=STEEL, ELSET=E_COL, SECTION=RECT")
    out.append(f"{fmt(args.bcol)}, {fmt(args.hcol)}")
    out.append("0., 0., 1.")
    out.append("*BEAM SECTION, MATERIAL=STEEL, ELSET=E_BEAM, SECTION=RECT")
    out.append(f"{fmt(args.bbeam)}, {fmt(args.hbeam)}")
    out.append("0., 0., 1.")

    # Supports
    out.append("*BOUNDARY")
    out.append("FIXED_BASE, 1, 3, 0.")

    # Translation continuity at hinges (release UR3 only).
    # Put dependent (duplicate) first to avoid SPC-on-dependent errors in cascade.
    out.append("** Translation continuity at hinges (release only UR3):")
    for dep,ind in [(n_tl_d, n_tl), (n_tr_d, n_tr)]:
        for dof in (1,2,3):
            out.append("*EQUATION")
            out.append("2")
            out.append(f"{dep}, {dof},  1.0, {ind}, {dof}, -1.0")

    # Activate UR3 (DOF6) on both hinge-side nodes via MPC to a dummy translation DOF:
    out.append("** Activate UR3 (DOF6) at hinge nodes: UR3(node) = U1(dummy).")
    for node in [n_tl, n_tl_d, n_tr, n_tr_d]:
        d = dummy[node]
        out.append("*EQUATION")
        out.append("2")
        out.append(f"{node}, 6,  1.0, {d}, 1, -1.0")
        # prevent dummy drift in U2/U3 (leave U1 free!)
        out.append("*BOUNDARY")
        out.append(f"{d}, 2, 3, 0.")

    # Springs
    for name,_,_,_,_ in hinges:
        write_spring_table(out, name, args.My, args.ky, args.kpost)

    # Earthquake amplitude
    write_amplitude(out, "EQX", t, ax)

    # STEP 1: gravity preload (static)
    out.append("*STEP")
    out.append("*STATIC")
    out.append("0.1, 1.0, 1e-8, 0.1")
    if args.gravity:
        out.append("*DLOAD")
        out.append("E_COL,  GRAV, 9.81, 0., -1., 0.")
        out.append("E_BEAM, GRAV, 9.81, 0., -1., 0.")
    # Energy output request must be present already in step 1 if you want it later in nonlinear analyses:
    out.append("*EL FILE, ELSET=E_COL")
    out.append("ELKE, ELSE")
    out.append("*EL FILE, ELSET=E_BEAM")
    out.append("ELKE, ELSE")
    out.append("*NODE PRINT, NSET=TOP_NODES")
    out.append("U")
    out.append("*END STEP")

    # STEP 2: dynamic earthquake
    out.append("*STEP, NLGEOM")
    out.append("*DYNAMIC")
    out.append(f"{fmt(args.dt)}, {fmt(args.T)}")
    if args.gravity:
        out.append("*DLOAD")
        out.append("E_COL,  GRAV, 9.81, 0., -1., 0.")
        out.append("E_BEAM, GRAV, 9.81, 0., -1., 0.")
    out.append("*DLOAD, AMPLITUDE=EQX")
    out.append("E_COL,  GRAV, 1.0, 1., 0., 0.")
    out.append("E_BEAM, GRAV, 1.0, 1., 0., 0.")
    out.append("*NODE PRINT, NSET=ROT_DUMMIES, FREQUENCY=1")
    out.append("U")
    out.append("*NODE PRINT, NSET=TOP_NODES, FREQUENCY=10")
    out.append("U")
    # repeat (ok; first appearance was step 1)
    out.append("*EL FILE, ELSET=E_COL")
    out.append("ELKE, ELSE")
    out.append("*EL FILE, ELSET=E_BEAM")
    out.append("ELKE, ELSE")
    out.append("*END STEP")

    args.out.write_text("\n".join(out) + "\n")
    mapping = {
        "rot_dummy_meaning": "U1(dummy) == UR3(node) [rad]",
        "hinges": [
            {"name":"H_TL","nodeA":n_tl,"nodeB":n_tl_d,"dummyA":d_tl,"dummyB":d_tl_d},
            {"name":"H_TR","nodeA":n_tr,"nodeB":n_tr_d,"dummyA":d_tr,"dummyB":d_tr_d},
        ],
    }
    Path("hinge_map.json").write_text(json.dumps(mapping, indent=2))
    print(f"Wrote {args.out}")
    print(f"Nodes: {len(nodes)} (incl. dummies)")
    print(f"Elements: {len(col_eids)+len(beam_eids)+len(hinge_eids)} (B31 + SPRING2)")
    print("Hinge map written to hinge_map.json")

if __name__ == "__main__":
    main()
