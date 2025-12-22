ccx_portal_eq_generator v4

Main fix: REF/ROT node IDs are allocated above max solid node id (no collisions).
This was the reason you still had "nonpositive jacobian" even after fixing hexa ordering.

Run:
  unzip -o ccx_portal_eq_generator_v4.zip -d ccx_portal_eq_generator_v4
  cd ccx_portal_eq_generator_v4

  python3 generate_portal_eq_ccx.py --eq earthquake_example.csv --out portal_eq.inp --gravity
  ../../src/CalculiX -i portal_eq
  python3 post_hinge_theta.py portal_eq.dat portal_eq.inp

Tip:
  Open portal_eq.inp and check the header lines:
    ** max_solid_nid=...
    ** REF/ROT: ...
If REF/ROT IDs are > max_solid_nid, you're safe.
