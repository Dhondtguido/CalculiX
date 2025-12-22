Portal-frame earthquake test generator for CalculiX (ccx)

Files:
  - generate_portal_eq_ccx.py        generator -> writes portal_eq.inp
  - earthquake_example.csv           synthetic accelerogram (t, ax_mps2)
  - post_hinge_theta.py              reads .dat and prints hinge rotations

Quick start:
  python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out portal_eq.inp --gravity
  ccx -i portal_eq
  python3 post_hinge_theta.py portal_eq.dat
