"""Check native I/L dimensions, section reference and CCX OFFSET1/2.

The four experiments load beam-node reference lines axially and check
the additional bending flexibility caused by centroid eccentricity.
"""
import math
import re
import sys

E=2.1e11
F=1000.
L=1.
NUM=re.compile(r"[+-]?\d+\.\d+E[+-]\d+",re.I)

def disp(path):
    active=False
    with open(path,encoding="utf-8",errors="replace") as f:
        for line in f:
            if line.startswith(" -4  DISP"):
                active=True
            elif active and line.startswith(" -3"):
                active=False
            elif active and line.startswith(" -1"):
                m=re.match(r"\s*(\d+)(.*)",line[3:])
                if m and int(m.group(1))==2:
                    values=[float(x) for x in NUM.findall(m.group(2))]
                    if len(values)>=3:
                        return values[0]
    raise RuntimeError("No node 2 DISP in "+path)

def props(rects):
    # Each rectangle is (width along section 1, depth along section 2,
    # section-1 centroid, section-2 centroid).
    area=sum(w*h for w,h,y,z in rects)
    yc=sum(w*h*y for w,h,y,z in rects)/area
    zc=sum(w*h*z for w,h,y,z in rects)/area
    iy=sum(w*h**3/12+w*h*(z-zc)**2 for w,h,y,z in rects)
    iz=sum(h*w**3/12+w*h*(y-yc)**2 for w,h,y,z in rects)
    iyz=sum(w*h*(y-yc)*(z-zc) for w,h,y,z in rects)
    theta=.5*math.atan2(-2*iyz,iy-iz) if abs(iyz)>1e-12 else 0.
    iy_p=.5*(iy+iz)+.5*(iy-iz)*math.cos(2*theta)-iyz*math.sin(2*theta)
    iz_p=.5*(iy+iz)-.5*(iy-iz)*math.cos(2*theta)+iyz*math.sin(2*theta)
    return area,yc,zc,iy_p,iz_p,theta

def axial_expected(pr, offy, offz):
    area,_,_,iy,iz,theta=pr
    oy=offy*math.cos(theta)+offz*math.sin(theta)
    oz=-offy*math.sin(theta)+offz*math.cos(theta)
    return F*L/E*(1/area+oy**2/iz+oz**2/iy)

# I input: l,h,b1,b2,t1,t2,t3
h,b1,b2,t1,t2,t3,l=.2,.12,.10,.012,.016,.008,.015
w=max(b1,b2)
i_rects=[(b1,t1,w/2,t1/2),
         (t3,h-t1-t2,w/2,t1+(h-t1-t2)/2),
         (b2,t2,w/2,h-t2/2)]
pi=props(i_rects)
i_result=axial_expected(pi,0.,pi[2]-l)

# L input: a,b,t1,t2
a,b,tf,tw=.12,.20,.012,.008
pl=props([(a,tf,a/2,tf/2),(tw,b-tf,tw/2,tf+(b-tf)/2)])
l_result=axial_expected(pl,pl[1],pl[2])

# Rect, CCX offsets point towards the reference line from section axis.
rect_props=props([(.1,.2,.05,.1)])
r_result=axial_expected(rect_props,-.25*.1,+.25*.2)

# Native BOX: a=.1 along n1, b=.2 along n2.
a,b,t=.1,.2,.01
box_props=props([
    (a,t,a/2,t/2),
    (a,t,a/2,b-t/2),
    (t,b-2*t,t/2,b/2),
    (t,b-2*t,a-t/2,b/2)])
b_result=axial_expected(box_props,-.25*a,+.25*b)

checks=[("28_b31e_abaqus_i_ref",i_result),
        ("29_b31e_abaqus_l_ref",l_result),
        ("30_b31e_ccx_offsets",r_result),
        ("31_b31e_ccx_box_offsets",b_result)]
failed=False
for name,expected in checks:
    actual=disp(name+".frd")
    ratio=actual/expected
    print("{}: axial={:.9g}, analytical={:.9g}, ratio={:.6f}".format(
        name,actual,expected,ratio))
    if not (math.isfinite(actual) and
            math.isclose(actual,expected,rel_tol=.03,abs_tol=1e-9)):
        failed=True
if failed:
    sys.exit("FAIL: native I/L reference or CCX offset conventions")
print("PASS: four native I/L and CCX offset benchmarks")
