ccx_portal_eq_generator v2

Fix: use CGX/CalculiX hexa node order to avoid nonpositive Jacobian.

Run:
  python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out portal_eq.inp --gravity
  ccx -i portal_eq
  python3 post_hinge_theta.py portal_eq.dat
