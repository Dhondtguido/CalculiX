#!/usr/bin/env python3
"""
v2: Fixes nonpositive Jacobian by using CalculiX/CGX hexa node ordering.

Hexa connectivity used (CGX-style):
  1:(0,0,0) 2:(0,1,0) 3:(1,1,0) 4:(1,0,0)
  5:(0,0,1) 6:(0,1,1) 7:(1,1,1) 8:(1,0,1)

Usage:
  python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out portal_eq.inp --gravity
  ccx -i portal_eq
"""
import argparse, csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

Vec = Tuple[float, float, float]

def vadd(a: Vec, b: Vec) -> Vec:
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def vscale(a: Vec, s: float) -> Vec:
    return (a[0]*s, a[1]*s, a[2]*s)

@dataclass
class Block:
    nodes: Dict[Tuple[int,int,int], int]
    node_coords: Dict[int, Vec]
    elements: List[Tuple[int, Tuple[int,int,int,int,int,int,int,int]]]
    face_nodes: Dict[str, List[int]]

def build_block(node_start: int, elem_start: int, origin: Vec, ex: Vec, ey: Vec, ez: Vec,
                nx: int, ny: int, nz: int) -> Tuple[Block, int, int]:
    nodes: Dict[Tuple[int,int,int], int] = {}
    node_coords: Dict[int, Vec] = {}
    nid = node_start

    for k in range(nz):
        zk = k/(nz-1) if nz>1 else 0.0
        for j in range(ny):
            yj = j/(ny-1) if ny>1 else 0.0
            for i in range(nx):
                xi = i/(nx-1) if nx>1 else 0.0
                p = origin
                p = vadd(p, vscale(ex, xi))
                p = vadd(p, vscale(ey, yj))
                p = vadd(p, vscale(ez, zk))
                nodes[(i,j,k)] = nid
                node_coords[nid] = p
                nid += 1

    def N(i:int,j:int,k:int)->int:
        return nodes[(i,j,k)]

    elements: List[Tuple[int, Tuple[int,int,int,int,int,int,int,int]]] = []
    eid = elem_start

    for k in range(nz-1):
        for j in range(ny-1):
            for i in range(nx-1):
                conn = (
                    N(i,   j,   k),
                    N(i,   j+1, k),
                    N(i+1, j+1, k),
                    N(i+1, j,   k),
                    N(i,   j,   k+1),
                    N(i,   j+1, k+1),
                    N(i+1, j+1, k+1),
                    N(i+1, j,   k+1),
                )
                elements.append((eid, conn))
                eid += 1

    face_nodes: Dict[str, List[int]] = {}
    face_nodes["i0"] = [nodes[(0, j, k)] for k in range(nz) for j in range(ny)]
    face_nodes["i1"] = [nodes[(nx-1, j, k)] for k in range(nz) for j in range(ny)]
    face_nodes["j0"] = [nodes[(i, 0, k)] for k in range(nz) for i in range(nx)]
    face_nodes["j1"] = [nodes[(i, ny-1, k)] for k in range(nz) for i in range(nx)]
    face_nodes["k0"] = [nodes[(i, j, 0)] for j in range(ny) for i in range(nx)]
    face_nodes["k1"] = [nodes[(i, j, nz-1)] for j in range(ny) for i in range(nx)]

    return Block(nodes, node_coords, elements, face_nodes), nid, eid

def read_eq_csv(path: Path) -> Tuple[List[float], List[float]]:
    t: List[float] = []
    ax: List[float] = []
    with path.open("r", newline="") as f:
        r = csv.DictReader(f)
        if not r.fieldnames:
            raise ValueError("empty csv")
        cols = [c.strip() for c in r.fieldnames]
        if "t" not in cols or "ax_mps2" not in cols:
            raise ValueError(f"{path} must have columns t, ax_mps2. Found {cols}")
        for row in r:
            t.append(float(row["t"]))
            ax.append(float(row["ax_mps2"]))
    if len(t) < 2:
        raise ValueError("need >=2 samples")
    return t, ax

def centroid(nids: List[int], coords: Dict[int, Vec]) -> Vec:
    sx = sy = sz = 0.0
    for n in nids:
        x,y,z = coords[n]
        sx += x; sy += y; sz += z
    m = float(len(nids))
    return (sx/m, sy/m, sz/m)

def write_nset(f, name: str, ids: List[int]):
    f.write(f"*NSET, NSET={name}\n")
    for i in range(0, len(ids), 16):
        f.write(", ".join(str(x) for x in ids[i:i+16]) + "\n")

def write_elset(f, name: str, ids: List[int]):
    f.write(f"*ELSET, ELSET={name}\n")
    for i in range(0, len(ids), 16):
        f.write(", ".join(str(x) for x in ids[i:i+16]) + "\n")

def write_inp(out_inp: Path, t: List[float], ax: List[float], *,
              L: float, H: float, b: float, h: float, t_z: float,
              nL: int, nH: int, nx: int, nz: int,
              E: float, nu: float, rho: float,
              hinge_my: float, hinge_k_el: float, hinge_k_post: float,
              dt0: float | None, gravity: bool) -> None:
    nid, eid = 1, 1

    colL, nid, eid = build_block(nid, eid, (0.0, 0.0, 0.0), (b,0.0,0.0), (0.0,H,0.0), (0.0,0.0,t_z),
                                 nx=nx, ny=nH+1, nz=nz)
    colR, nid, eid = build_block(nid, eid, (L-b, 0.0, 0.0), (b,0.0,0.0), (0.0,H,0.0), (0.0,0.0,t_z),
                                 nx=nx, ny=nH+1, nz=nz)
    beam, nid, eid = build_block(nid, eid, (0.0, H-h, 0.0), (L,0.0,0.0), (0.0,h,0.0), (0.0,0.0,t_z),
                                 nx=nL+1, ny=nx, nz=nz)

    E_ALL = [e for e,_ in colL.elements] + [e for e,_ in colR.elements] + [e for e,_ in beam.elements]
    COL_L_BASE = colL.face_nodes["j0"]
    COL_R_BASE = colR.face_nodes["j0"]
    COL_L_TOP  = colL.face_nodes["j1"]
    COL_R_TOP  = colR.face_nodes["j1"]
    BEM_L = beam.face_nodes["i0"]
    BEM_R = beam.face_nodes["i1"]

    # ref/rot nodes
    REF_L_COL, ROT_L_COL = 1001, 1101
    REF_L_BEM, ROT_L_BEM = 1002, 1102
    REF_R_COL, ROT_R_COL = 1003, 1103
    REF_R_BEM, ROT_R_BEM = 1004, 1104

    coords_all: Dict[int, Vec] = {}
    for blk in (colL, colR, beam):
        coords_all.update(blk.node_coords)

    c_L_COL = centroid(COL_L_TOP, coords_all)
    c_L_BEM = centroid(BEM_L, coords_all)
    c_R_COL = centroid(COL_R_TOP, coords_all)
    c_R_BEM = centroid(BEM_R, coords_all)

    T = t[-1]
    if dt0 is None:
        dt0 = min(0.005, max(1e-4, t[1]-t[0]))

    th_y = hinge_my / hinge_k_el
    th_max = max(10*th_y, 0.02)
    m_max = hinge_my + hinge_k_post*(th_max - th_y)
    table = [(-th_max,-m_max), (-th_y,-hinge_my), (0.0,0.0), (th_y,hinge_my), (th_max,m_max)]

    with out_inp.open("w") as f:
        f.write("** Portal frame (solids) + M-theta hinges, EQ via *DLOAD GRAV AMPLITUDE\n")
        f.write("*HEADING\nPortal EQ test (v2)\n\n")
        f.write("*NODE\n")
        for blk in (colL, colR, beam):
            for nid0 in sorted(blk.node_coords):
                x,y,z = blk.node_coords[nid0]
                f.write(f"{nid0}, {x:.8e}, {y:.8e}, {z:.8e}\n")
        for nid_extra, p in [
            (REF_L_COL, c_L_COL),(ROT_L_COL,c_L_COL),
            (REF_L_BEM, c_L_BEM),(ROT_L_BEM,c_L_BEM),
            (REF_R_COL, c_R_COL),(ROT_R_COL,c_R_COL),
            (REF_R_BEM, c_R_BEM),(ROT_R_BEM,c_R_BEM),
        ]:
            f.write(f"{nid_extra}, {p[0]:.8e}, {p[1]:.8e}, {p[2]:.8e}\n")
        f.write("\n")

        f.write("*ELEMENT, TYPE=C3D8R, ELSET=E_COL_L\n")
        for eid0, conn in colL.elements:
            f.write(f"{eid0}, " + ", ".join(str(n) for n in conn) + "\n")
        f.write("*ELEMENT, TYPE=C3D8R, ELSET=E_COL_R\n")
        for eid0, conn in colR.elements:
            f.write(f"{eid0}, " + ", ".join(str(n) for n in conn) + "\n")
        f.write("*ELEMENT, TYPE=C3D8R, ELSET=E_BEAM\n")
        for eid0, conn in beam.elements:
            f.write(f"{eid0}, " + ", ".join(str(n) for n in conn) + "\n")
        f.write("\n")

        f.write("*ELEMENT, TYPE=SPRING2, ELSET=E_HINGE_L\n9001, 1101, 1102\n")
        f.write("*ELEMENT, TYPE=SPRING2, ELSET=E_HINGE_R\n9002, 1103, 1104\n\n")

        write_elset(f, "E_ALL", E_ALL)
        write_nset(f, "COL_L_BASE", COL_L_BASE)
        write_nset(f, "COL_R_BASE", COL_R_BASE)
        write_nset(f, "COL_L_TOP", COL_L_TOP)
        write_nset(f, "COL_R_TOP", COL_R_TOP)
        write_nset(f, "BEM_L", BEM_L)
        write_nset(f, "BEM_R", BEM_R)
        write_nset(f, "ROT_ALL", [1101,1102,1103,1104])
        write_nset(f, "REF_TOP", [1003])
        f.write("\n")

        f.write("*MATERIAL, NAME=STEEL\n*ELASTIC\n")
        f.write(f"{E:.8e}, {nu:.6f}\n*DENSITY\n{rho:.8e}\n\n")
        f.write("*SOLID SECTION, ELSET=E_ALL, MATERIAL=STEEL\n\n")

        f.write("*RIGID BODY, NSET=COL_L_TOP, REF NODE=1001, ROT NODE=1101\n")
        f.write("*RIGID BODY, NSET=BEM_L,     REF NODE=1002, ROT NODE=1102\n")
        f.write("*RIGID BODY, NSET=COL_R_TOP, REF NODE=1003, ROT NODE=1103\n")
        f.write("*RIGID BODY, NSET=BEM_R,     REF NODE=1004, ROT NODE=1104\n\n")

        for d in (1,2,3):
            f.write("*EQUATION\n2\n")
            f.write(f"1001, {d}, 1.0, 1002, {d}, -1.0\n")
            f.write("*EQUATION\n2\n")
            f.write(f"1003, {d}, 1.0, 1004, {d}, -1.0\n")
        f.write("\n")

        f.write("*SPRING, ELSET=E_HINGE_L, NONLINEAR\n3, 3\n")
        for th,m in table:
            f.write(f"{th:.8e}, {m:.8e}\n")
        f.write("*SPRING, ELSET=E_HINGE_R, NONLINEAR\n3, 3\n")
        for th,m in table:
            f.write(f"{th:.8e}, {m:.8e}\n")
        f.write("\n")

        f.write("*BOUNDARY\n")
        f.write("COL_L_BASE, 1, 3, 0.\n")
        f.write("COL_R_BASE, 1, 3, 0.\n\n")

        f.write("*AMPLITUDE, NAME=EQX\n")
        for ti, ai in zip(t, ax):
            f.write(f"{ti:.8e}, {ai:.8e}\n")
        f.write("\n")

        if gravity:
            f.write("*STEP\n*STATIC\n*DLOAD\n")
            f.write("E_ALL, GRAV, 9.81, 0., -1., 0.\n")
            f.write("*END STEP\n\n")

        f.write("*STEP\n*DYNAMIC\n")
        f.write(f"{dt0:.8e}, {T:.8e}\n")
        f.write("*DLOAD, AMPLITUDE=EQX\n")
        f.write("E_ALL, GRAV, 1.0, 1., 0., 0.\n\n")
        f.write("*NODE FILE\nU\n")
        f.write("*NODE PRINT, NSET=ROT_ALL\nU\n")
        f.write("*NODE PRINT, NSET=REF_TOP\nU\n")
        f.write("*END STEP\n")

    print(f"Wrote {out_inp} (T={T:.3f}s, dt0={dt0:g}s)")
    print(f"Elems solids={len(E_ALL)} + 2 springs")
    print(f"theta_y={th_y:.3e} rad, My={hinge_my:.3e} N*m")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eq", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("portal_eq.inp"))
    ap.add_argument("--L", type=float, default=4.0)
    ap.add_argument("--H", type=float, default=3.0)
    ap.add_argument("--b", type=float, default=0.30)
    ap.add_argument("--h", type=float, default=0.30)
    ap.add_argument("--t_z", type=float, default=0.20)
    ap.add_argument("--nL", type=int, default=30)
    ap.add_argument("--nH", type=int, default=30)
    ap.add_argument("--nx", type=int, default=4)
    ap.add_argument("--nz", type=int, default=3)
    ap.add_argument("--E", type=float, default=210e9)
    ap.add_argument("--nu", type=float, default=0.30)
    ap.add_argument("--rho", type=float, default=7850.0)
    ap.add_argument("--hinge_my", type=float, default=2.0e5)
    ap.add_argument("--hinge_k_el", type=float, default=1.0e8)
    ap.add_argument("--hinge_k_post", type=float, default=1.0e6)
    ap.add_argument("--dt0", type=float, default=None)
    ap.add_argument("--gravity", action="store_true")
    args = ap.parse_args()

    t, ax = read_eq_csv(args.eq)
    write_inp(args.out, t, ax,
              L=args.L, H=args.H, b=args.b, h=args.h, t_z=args.t_z,
              nL=args.nL, nH=args.nH, nx=args.nx, nz=args.nz,
              E=args.E, nu=args.nu, rho=args.rho,
              hinge_my=args.hinge_my, hinge_k_el=args.hinge_k_el, hinge_k_post=args.hinge_k_post,
              dt0=args.dt0, gravity=args.gravity)

if __name__ == "__main__":
    main()
