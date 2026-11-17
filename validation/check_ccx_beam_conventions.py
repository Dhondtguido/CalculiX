"""B31T/B31E analytical tip-displacement regression (SI units)."""
import math
import re
import sys

FLOAT = re.compile(r"[+-]?\d+\.\d+E[+-]\d+", re.I)

def tip(path):
    reading = False
    for line in open(path, encoding="utf8", errors="replace"):
        if line.startswith(" -4  DISP"):
            reading = True
            continue
        if reading and line.startswith(" -3"):
            reading = False
        if reading and line.startswith(" -1"):
            m = re.match(r"\s*(\d+)(.*)", line[3:])
            if m and int(m.group(1)) == 2:
                numbers = [float(v) for v in FLOAT.findall(m.group(2))]
                if len(numbers) >= 3:
                    return numbers[:3]
    raise RuntimeError("No node 2 DISP in "+path)

E, nu, L, F = 2.1e11, .3, 1., 1000.
A = .1*.2
G = E/(2*(1+nu))
k = 10*(1+nu)/(12+11*nu)
Iy = .1*.2**3/12
Iz = .2*.1**3/12
Iellipse = math.pi*.05*.1**3/4
Ibox = (.1*.2**3-(.1-2*.008)*(.2-2*.012)**3)/12
bending = lambda I: -F*L**3/(3*E*I)
cases = [
    ("22_b31e_rect_strong",2,bending(Iy)),
    ("23_b31e_rect_weak",1,bending(Iz)),
    ("24_b31e_default_normal",1,bending(Iy)),
    ("25_b31t_rect_strong",2,bending(Iy)-F*L/(k*G*A)),
    ("26_b31e_ellipse",2,bending(Iellipse)),
    ("27_b31e_native_box",2,bending(Ibox))
]
failed=False
for name, component, expected in cases:
    actual=tip(name+".frd")[component]
    ratio=actual/expected
    print("{}: result={:.9g}, theory={:.9g}, ratio={:.5f}".format(
        name,actual,expected,ratio))
    if not math.isfinite(actual) or not math.isclose(
        actual,expected,rel_tol=.03,abs_tol=1.e-9):
        failed=True
if failed:
    sys.exit("FAIL: B31T/B31E geometry or orientation regression")
print("PASS: six analytical beam cases")
