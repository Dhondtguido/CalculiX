"""Verify local BFOR shear/moment magnitudes and M1/M2 releases."""
import csv
import math
import re
import sys

NUMBER = re.compile(r"[+-]?\d+\.\d+E[+-]\d+", re.I)
def records(path, field):
    found=False
    rows={}
    for line in open(path,encoding="utf-8",errors="replace"):
        if line.startswith(" -4  "+field):
            found=True
            continue
        if found and line.startswith(" -3"):
            found=False
        if found and line.startswith(" -1"):
            m=re.match(r"\s*(\d+)(.*)",line[3:])
            if m:
                vals=[float(v) for v in NUMBER.findall(m.group(2))]
                if vals:
                    rows[int(m.group(1))]=vals
    return rows

def close(x, expected, rel=.02, abs_tol=1e-8):
    return math.isfinite(x) and math.isclose(
        x,expected,rel_tol=rel,abs_tol=abs_tol)

fail=False
for name in ("32_b31e","33_b31t"):
    csvpath=name+"_forces.csv"
    with open(csvpath,newline="",encoding="utf-8") as f:
        data=list(csv.DictReader(f))
    root=next((row for row in data if int(row["Station_Pct"])==0),None)
    tip=next((row for row in data if int(row["Station_Pct"])==100),None)
    if root is None or tip is None:
        sys.exit("Missing station 0/100 for "+name)
    csv_v=abs(float(root["Vz_Shear"]))
    csv_m=abs(float(root["My_Bending"]))
    csv_endm=abs(float(tip["My_Bending"]))
    frd=records(name+"_force_units.frd","BFOR")
    if 1 not in frd or 2 not in frd:
        sys.exit("Missing BFOR field at beam nodes for "+name)
    frd_v=abs(frd[1][2])
    frd_m=abs(frd[1][4])
    print("{}: CSV Vz={:.5g}N My={:.5g}Nm, FRD Vz={:.5g} My={:.5g}".format(
        name,csv_v,csv_m,frd_v,frd_m))
    if not all((close(csv_v,1000),close(csv_m,1000),
                close(frd_v,1000),close(frd_m,1000),
                close(csv_endm,0,abs_tol=1e-3))):
        fail=True

z0=records("34_b31e_fixed_fixed.frd","DISP")[2][2]
z1=records("35_b31e_release_m1.frd","DISP")[2][2]
z2=records("36_b31e_release_m2.frd","DISP")[2][2]
Iy=.1*.2**3/12
baseline=-1000*2**3/(192*2.1e11*Iy)
print("Releases: fixed={:.9g}, M1={:.9g}, M2={:.9g}, theory={:.9g}".format(
    z0,z1,z2,baseline))
if not close(z0,baseline,rel=.03) or not close(z2,z0,rel=.03):
    fail=True
if not (z1<0 and abs(z1)>abs(z0)*1.1):
    fail=True
if fail:
    sys.exit("FAIL: BFOR units / Abaqus M1/M2 axis regression")
print("PASS: two force-output and two bending-release comparisons")
