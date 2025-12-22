ccx_portal_eq_generator v6

Files:
- generate_portal_eq_ccx.py
- earthquake_example.csv
- run_all.sh
- post_hinge_theta.py

From ~/CalculiX/tests:
  unzip -o ccx_portal_eq_generator_v6.zip -d ccx_portal_eq_generator_v6
  cd ccx_portal_eq_generator_v6
  bash run_all.sh

Sanity:
  grep -n "ENERGY PRINT" portal_eq.inp
