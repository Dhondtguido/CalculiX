ccx_portal_eq_generator_v8
=========================

Goal
----
Generate a *beam* portal-frame model in CalculiX with **rotational hinges implemented via SPRING2 on DOF 6**
(using the same "piecewise M-θ" idea as your rot_spring2 tests), and drive it with a **single horizontal
earthquake acceleration history**.

Why this version works better
-----------------------------
- No solid bricks => no Jacobian sign headaches.
- No *RIGID BODY MPCs on beam nodes (ccx forbids that).
- We explicitly **activate the rotational DOF** at hinge nodes by adding an MPC-equation that references DOF 6
  (this is needed because ccx internally expands B31 beams into solids and only creates "knot" rotational DOFs
  if they appear in an SPC/MPC or in moment loading).

Files
-----
- generate_portal_eq_ccx.py   Generates portal_eq.inp (and hinge_map.json).
- earthquake_example.csv      Example a_x(t) in m/s^2 (t in s).
- post_hinge_theta.py         Extracts hinge Δθ(t) from the dummy nodes, and computes M(Δθ) from the spring table.
- run_all.sh                  Convenience runner (set CCX env var to your CalculiX executable).

Quick start
-----------
1) unzip -o ccx_portal_eq_generator_v8.zip -d ccx_portal_eq_generator_v8
2) cd ccx_portal_eq_generator_v8/ccx_portal_eq_generator_v8
3) python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out portal_eq.inp --gravity
4) CCX=../../../src/CalculiX bash run_all.sh
5) python3 post_hinge_theta.py portal_eq.dat portal_eq.inp

Notes
-----
- Units: m, N, kg, s (E in Pa, density in kg/m^3, moments in N*m).
- Base acceleration is modeled as an equivalent body-acceleration using *DLOAD, GRAV with AMPLITUDE=EQX.
