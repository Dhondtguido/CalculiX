#!/usr/bin/env python3
"""
Generate CalculiX input for portal frame with solid elements and plastic hinges.
Updated for Problema 4 (Tarea3_DC_2015) specifications.

Default geometry:
- Frame: H=3.0m, L=5.0m (Problema 4)
- Columns: 40cm × 60cm (S1 from Problema 2)
- Beam: 25cm × 50cm (S2 from Problema 2)

Hinges:
- Bilinear M-θ curve (simplified for now, can be updated with RC section analysis)
- Default: θ_y=0.002 rad, M_y=200 kN·m (placeholder, update with real values)

Usage:
    python3 generate_portal_eq_ccx.py --eq earthquake_A0.1g.csv --out portal_A0.1g.inp
"""

import argparse
import math, csv, json
from pathlib import Path

G = 9.81

def read_eq_csv(path: Path):
    """Read CSV with (t, acc). If |acc|max <= 5 -> assume g, else m/s^2."""
    rows=[]
    with open(path, "r", newline="") as f:
        r = csv.reader(f)
        _ = next(r, None)  # header (optional)
        for row in r:
            if not row:
                continue
            if row[0].strip().startswith("#"):
                continue
            rows.append((float(row[0]), float(row[1])))
    maxabs = max(abs(v) for _,v in rows) if rows else 0.0
    is_g = (maxabs <= 5.0)
    return rows, is_g

def write_amplitude(f, name, tv_pairs):
    f.write(f"*AMPLITUDE, NAME={name}\n")
    per_line = 4
    buf=[]
    for t,v in tv_pairs:
        buf.append(f"{t:.6f}, {v:.8e}")
        if len(buf) == per_line:
            f.write(", ".join(buf) + "\n")
            buf=[]
    if buf:
        f.write(", ".join(buf) + "\n")

def structured_block_nodes(x0,x1,y0,y1,z0,z1,nx,ny,nz,start_id):
    """Return nodes dict {id:(x,y,z)} and grid mapping (i,j,k)->node_id."""
    nodes={}
    grid={}
    dx=(x1-x0)/nx
    dy=(y1-y0)/ny
    dz=(z1-z0)/nz
    nid=start_id
    for k in range(nz+1):
        z=z0+dz*k
        for j in range(ny+1):
            y=y0+dy*j
            for i in range(nx+1):
                x=x0+dx*i
                nodes[nid]=(x,y,z)
                grid[(i,j,k)]=nid
                nid += 1
    return nodes, grid, nid

def structured_block_hex_elems(grid,nx,ny,nz,start_eid):
    """Return list of (eid,n1..n8) for C3D8R with standard ordering."""
    elems=[]
    eid=start_eid
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n1=grid[(i,  j,  k)]
                n2=grid[(i+1,j,  k)]
                n3=grid[(i+1,j+1,k)]
                n4=grid[(i,  j+1,k)]
                n5=grid[(i,  j,  k+1)]
                n6=grid[(i+1,j,  k+1)]
                n7=grid[(i+1,j+1,k+1)]
                n8=grid[(i,  j+1,k+1)]
                elems.append((eid,n1,n2,n3,n4,n5,n6,n7,n8))
                eid += 1
    return elems, eid

def select_face_nodes(nodes, *, x=None, y=None, z=None, tol=1e-9, y_range=None, x_range=None):
    out=[]
    for nid,(xx,yy,zz) in nodes.items():
        ok=True
        if x is not None and abs(xx-x) > tol: ok=False
        if y is not None and abs(yy-y) > tol: ok=False
        if z is not None and abs(zz-z) > tol: ok=False
        if y_range is not None:
            y0,y1=y_range
            if not (yy >= y0 - tol and yy <= y1 + tol): ok=False
        if x_range is not None:
            x0,x1=x_range
            if not (xx >= x0 - tol and xx <= x1 + tol): ok=False
        if ok:
            out.append(nid)
    return sorted(out)

def write_nset(f, name, ids):
    f.write(f"*NSET, NSET={name}\n")
    per=16
    for i in range(0, len(ids), per):
        f.write(", ".join(str(v) for v in ids[i:i+per]) + "\n")

def main():
    ap=argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eq", type=str, required=True, help="CSV with t, acc (g or m/s^2)")
    ap.add_argument("--out", type=str, default="portal_eq.inp")

    # Geometry (Problema 4 defaults)
    ap.add_argument("--H", type=float, default=3.0, help="Column height (m)")
    ap.add_argument("--L", type=float, default=5.0, help="Frame span (m)")

    # Column section (S1: 40cm × 60cm from Problema 2)
    ap.add_argument("--col_b", type=float, default=0.40, help="Column width in-plane (m)")
    ap.add_argument("--col_h", type=float, default=0.60, help="Column height/depth (m)")
    ap.add_argument("--col_t", type=float, default=0.20, help="Column out-of-plane thickness (m)")

    # Beam section (S2: 25cm × 50cm from Problema 2)
    ap.add_argument("--beam_b", type=float, default=0.25, help="Beam width (m)")
    ap.add_argument("--beam_h", type=float, default=0.50, help="Beam height/depth (m)")
    ap.add_argument("--beam_t", type=float, default=0.20, help="Beam out-of-plane thickness (m)")

    # Mesh density
    ap.add_argument("--nx_col", type=int, default=2, help="Elems across column width")
    ap.add_argument("--ny_col", type=int, default=16, help="Elems along column height")
    ap.add_argument("--nz", type=int, default=1, help="Elems through thickness (all members)")
    ap.add_argument("--nx_beam", type=int, default=20, help="Elems along beam span")
    ap.add_argument("--ny_beam", type=int, default=2, help="Elems across beam height")

    # Material (concrete)
    ap.add_argument("--E", type=float, default=25e9, help="Young's modulus (Pa)")
    ap.add_argument("--nu", type=float, default=0.20, help="Poisson's ratio")
    ap.add_argument("--rho", type=float, default=2400.0, help="Density (kg/m³)")

    # Analysis parameters
    ap.add_argument("--dt", type=float, default=None, help="Time step (auto if None)")
    ap.add_argument("--damp_alpha", type=float, default=0.0)
    ap.add_argument("--damp_beta", type=float, default=2.0e-4)

    # Hinge parameters (simplified, can be updated with RC section analysis)
    ap.add_argument("--hinge_theta_y", type=float, default=0.002, help="Yield rotation (rad)")
    ap.add_argument("--hinge_My", type=float, default=200e3, help="Yield moment (N·m)")
    ap.add_argument("--hinge_post_ratio", type=float, default=0.02, help="Post-yield stiffness ratio")

    ap.add_argument("--energy", action="store_true", help="Include ENER output")

    args=ap.parse_args()

    tv, is_g = read_eq_csv(Path(args.eq))
    if not tv:
        raise SystemExit("Earthquake CSV is empty")

    t_end = tv[-1][0]
    t_end_eff = t_end + 1.0e-9
    tv_mps2=[(t, val*G if is_g else val) for t,val in tv]

    if args.dt is None:
        dts=[tv[i+1][0]-tv[i][0] for i in range(len(tv)-1)]
        dt=min(dts) if dts else 0.0025
    else:
        dt=args.dt

    dt = round(dt, 12)

    # Geometry
    L, H = args.L, args.H
    col_b, col_h, col_t = args.col_b, args.col_h, args.col_t
    beam_b, beam_h, beam_t = args.beam_b, args.beam_h, args.beam_t

    nx_col, ny_col, nz = args.nx_col, args.ny_col, args.nz
    nx_beam, ny_beam = args.nx_beam, args.ny_beam

    nodes={}
    nid=1
    eid=1

    # Left column: x ∈ [0, col_b], y ∈ [0, H], z ∈ [0, col_t]
    n_colL, grid_colL, nid = structured_block_nodes(
        0.0, col_b, 0.0, H, 0.0, col_t, nx_col, ny_col, nz, nid)
    e_colL, eid = structured_block_hex_elems(grid_colL, nx_col, ny_col, nz, eid)

    # Right column: x ∈ [L-col_b, L], y ∈ [0, H], z ∈ [0, col_t]
    n_colR, grid_colR, nid = structured_block_nodes(
        L-col_b, L, 0.0, H, 0.0, col_t, nx_col, ny_col, nz, nid)
    e_colR, eid = structured_block_hex_elems(grid_colR, nx_col, ny_col, nz, eid)

    # Beam: x ∈ [col_b, L-col_b], y ∈ [H-beam_h, H], z ∈ [0, beam_t]
    n_beam, grid_beam, nid = structured_block_nodes(
        col_b, L-col_b, H-beam_h, H, 0.0, beam_t, nx_beam, ny_beam, nz, nid)
    e_beam, eid = structured_block_hex_elems(grid_beam, nx_beam, ny_beam, nz, eid)

    for d in (n_colL, n_colR, n_beam):
        nodes.update(d)

    # Joint faces (connection zones)
    # Left joint: column inner face at x=col_b, y ∈ [H-beam_h, H]
    face_L_col = select_face_nodes(n_colL, x=col_b, y_range=(H-beam_h, H))
    face_L_bem = select_face_nodes(n_beam, x=col_b, y_range=(H-beam_h, H))

    # Right joint: column inner face at x=L-col_b, y ∈ [H-beam_h, H]
    face_R_col = select_face_nodes(n_colR, x=L-col_b, y_range=(H-beam_h, H))
    face_R_bem = select_face_nodes(n_beam, x=L-col_b, y_range=(H-beam_h, H))

    # Base fixity
    base_L = select_face_nodes(n_colL, y=0.0)
    base_R = select_face_nodes(n_colR, y=0.0)
    base_all = sorted(set(base_L + base_R))

    all_solid_nodes = sorted(list(n_colL.keys()) + list(n_colR.keys()) + list(n_beam.keys()))

    # REF/ROT nodes for rigid bodies at joint centers
    zc = col_t/2.0
    yjc = H - beam_h/2.0

    refrot = {
        "REF_L_COL": (col_b, yjc, zc),
        "ROT_L_COL": (col_b, yjc, zc),
        "REF_L_BEM": (col_b, yjc, zc),
        "ROT_L_BEM": (col_b, yjc, zc),
        "REF_R_COL": (L-col_b, yjc, zc),
        "ROT_R_COL": (L-col_b, yjc, zc),
        "REF_R_BEM": (L-col_b, yjc, zc),
        "ROT_R_BEM": (L-col_b, yjc, zc),
    }
    ref_ids={}
    for name, xyz in refrot.items():
        ref_ids[name]=nid
        nodes[nid]=xyz
        nid += 1

    # SPRING2 elements (two hinges)
    spr_elems=[]
    for a,b in [("ROT_L_COL","ROT_L_BEM"), ("ROT_R_COL","ROT_R_BEM")]:
        spr_elems.append((eid, ref_ids[a], ref_ids[b]))
        eid += 1

    # Bilinear hinge curve M(theta)
    th_y=args.hinge_theta_y
    My=args.hinge_My
    k_el = My/th_y
    k_pl = args.hinge_post_ratio * k_el
    th2 = 10.0*th_y
    M2  = My + k_pl*(th2-th_y)
    hinge_curve=[(0.0,0.0),(th_y,My),(th2,M2)]

    out_path = Path(args.out)

    # Hinge map for postprocessing
    hinge_map={
        "hinges":[
            {"name":"H_L","rotA":ref_ids["ROT_L_COL"],"rotB":ref_ids["ROT_L_BEM"]},
            {"name":"H_R","rotA":ref_ids["ROT_R_COL"],"rotB":ref_ids["ROT_R_BEM"]},
        ],
        "curve_theta_M": hinge_curve,
        "units":{"theta":"rad","M":"N*m"},
        "geometry": {
            "H": H,
            "L": L,
            "column_section": {"b": col_b, "h": col_h, "t": col_t},
            "beam_section": {"b": beam_b, "h": beam_h, "t": beam_t}
        }
    }
    out_path.with_suffix(".hinge_map.json").write_text(json.dumps(hinge_map, indent=2))

    with open(out_path,"w") as f:
        f.write("*HEADING\n")
        f.write(f"Portal frame (Problema 4): H={H}m, L={L}m\n")
        f.write(f"** Columns: {col_b*100:.0f}x{col_h*100:.0f}cm, Beam: {beam_b*100:.0f}x{beam_h*100:.0f}cm\n")
        f.write("** Generated by generate_portal_eq_ccx.py\n")
        f.write("*NODE\n")
        for nid0,(x,y,z) in sorted(nodes.items()):
            f.write(f"{nid0}, {x:.8f}, {y:.8f}, {z:.8f}\n")

        f.write("** SOLIDS\n")
        f.write("*ELEMENT, TYPE=C3D8R, ELSET=E_COL_L\n")
        for e in e_colL:
            eid0,*nn=e
            f.write(f"{eid0}, " + ", ".join(str(n) for n in nn) + "\n")
        f.write("*ELEMENT, TYPE=C3D8R, ELSET=E_COL_R\n")
        for e in e_colR:
            eid0,*nn=e
            f.write(f"{eid0}, " + ", ".join(str(n) for n in nn) + "\n")
        f.write("*ELEMENT, TYPE=C3D8R, ELSET=E_BEAM\n")
        for e in e_beam:
            eid0,*nn=e
            f.write(f"{eid0}, " + ", ".join(str(n) for n in nn) + "\n")

        f.write("** HINGE SPRINGS\n")
        f.write("*ELEMENT, TYPE=SPRING2, ELSET=E_HINGE\n")
        for eid0,n1,n2 in spr_elems:
            f.write(f"{eid0}, {n1}, {n2}\n")

        # Sets
        write_nset(f,"BASE_ALL",base_all)
        write_nset(f,"ALL_SOLID",all_solid_nodes)
        write_nset(f,"FACE_L_COL",face_L_col)
        write_nset(f,"FACE_L_BEM",face_L_bem)
        write_nset(f,"FACE_R_COL",face_R_col)
        write_nset(f,"FACE_R_BEM",face_R_bem)
        for name, nid0 in ref_ids.items():
            write_nset(f,name,[nid0])
        write_nset(f,"ROT_ALL",[ref_ids["ROT_L_COL"],ref_ids["ROT_L_BEM"],ref_ids["ROT_R_COL"],ref_ids["ROT_R_BEM"]])

        # Material
        f.write("*MATERIAL, NAME=MAT\n")
        f.write("*ELASTIC\n")
        f.write(f"{args.E:.6e}, {args.nu:.6f}\n")
        f.write("*DENSITY\n")
        f.write(f"{args.rho:.6f}\n")

        f.write("*SOLID SECTION, ELSET=E_COL_L, MATERIAL=MAT\n")
        f.write("*SOLID SECTION, ELSET=E_COL_R, MATERIAL=MAT\n")
        f.write("*SOLID SECTION, ELSET=E_BEAM,  MATERIAL=MAT\n")

        # Rigid bodies
        f.write("** RIGID BODIES ON JOINT FACES\n")
        f.write(f"*RIGID BODY, NSET=FACE_L_COL, REF NODE={ref_ids['REF_L_COL']}, ROT NODE={ref_ids['ROT_L_COL']}\n")
        f.write(f"*RIGID BODY, NSET=FACE_L_BEM, REF NODE={ref_ids['REF_L_BEM']}, ROT NODE={ref_ids['ROT_L_BEM']}\n")
        f.write(f"*RIGID BODY, NSET=FACE_R_COL, REF NODE={ref_ids['REF_R_COL']}, ROT NODE={ref_ids['ROT_R_COL']}\n")
        f.write(f"*RIGID BODY, NSET=FACE_R_BEM, REF NODE={ref_ids['REF_R_BEM']}, ROT NODE={ref_ids['ROT_R_BEM']}\n")

        # Translation compatibility
        def write_eq(refA, refB, dof):
            f.write("*EQUATION\n2\n")
            f.write(f"{refA}, {dof}, 1.0, {refB}, {dof}, -1.0\n")
        for dof in (1,2):
            write_eq(ref_ids["REF_L_COL"], ref_ids["REF_L_BEM"], dof)
            write_eq(ref_ids["REF_R_COL"], ref_ids["REF_R_BEM"], dof)

        # Boundary conditions
        f.write("*BOUNDARY\n")
        f.write("BASE_ALL, 1, 3, 0.\n")
        for ref in ("REF_L_COL","REF_L_BEM","REF_R_COL","REF_R_BEM"):
            f.write(f"{ref}, 3, 3, 0.\n")
        for rot in ("ROT_L_COL","ROT_L_BEM","ROT_R_COL","ROT_R_BEM"):
            f.write(f"{rot}, 1, 2, 0.\n")

        # Hinge law
        f.write("*SPRING, ELSET=E_HINGE, NONLINEAR\n")
        f.write("3, 3\n")
        for th,M in hinge_curve:
            f.write(f"{th:.8e}, {M:.8e}\n")

        # Amplitude
        write_amplitude(f, "EQX", tv_mps2)

        # Step 1: gravity
        f.write("*STEP\n*STATIC\n")
        f.write("*DLOAD\n")
        f.write(f"E_COL_L, GRAV, {G:.8e}, 0., -1., 0.\n")
        f.write(f"E_COL_R, GRAV, {G:.8e}, 0., -1., 0.\n")
        f.write(f"E_BEAM,  GRAV, {G:.8e}, 0., -1., 0.\n")
        f.write("*NODE PRINT, NSET=ROT_ALL\nU\n")
        f.write("*NODE FILE\nU\n")
        f.write("*EL FILE\nS, E, ENER\n")
        f.write("*END STEP\n")

        # Step 2: dynamic
        ninc = int(math.ceil(t_end_eff / dt)) + 50
        f.write(f"*STEP, INC={ninc}\n*DYNAMIC\n")
        f.write(f"{dt:.12f}, {t_end_eff:.12f}, {dt:.12f}, {dt:.12f}\n")

        if args.damp_alpha != 0.0 or args.damp_beta != 0.0:
            f.write(f"*DAMPING,ALPHA={args.damp_alpha:.8e},BETA={args.damp_beta:.8e}\n")
        f.write("*DLOAD\n")
        f.write(f"E_COL_L, GRAV, {G:.8e}, 0., -1., 0.\n")
        f.write(f"E_COL_R, GRAV, {G:.8e}, 0., -1., 0.\n")
        f.write(f"E_BEAM,  GRAV, {G:.8e}, 0., -1., 0.\n")
        f.write("E_COL_L, GRAV, 1.0, 1., 0., 0., EQX\n")
        f.write("E_COL_R, GRAV, 1.0, 1., 0., 0., EQX\n")
        f.write("E_BEAM,  GRAV, 1.0, 1., 0., 0., EQX\n")

        f.write("*NODE PRINT, NSET=ROT_ALL, FREQUENCY=20\nU\n")
        f.write("*NODE FILE, FREQUENCY=20\nU\n")
        f.write("*EL FILE, FREQUENCY=20\n")
        if args.energy:
            f.write("S, E, ENER\n")
        else:
            f.write("S, E\n")
        f.write("*END STEP\n")

    print(f"Wrote {out_path}")
    print(f"  Geometry: H={H}m, L={L}m")
    print(f"  Columns:  {col_b*100:.0f} × {col_h*100:.0f} cm")
    print(f"  Beam:     {beam_b*100:.0f} × {beam_h*100:.0f} cm")
    print(f"  Nodes={len(nodes)} solids={len(e_colL)+len(e_colR)+len(e_beam)} springs={len(spr_elems)}")
    print(f"  EQ: t_end={t_end:.3f}s dt={dt:.6f}s units={'g' if is_g else 'm/s^2'}")

if __name__ == "__main__":
    main()
