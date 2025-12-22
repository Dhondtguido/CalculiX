ccx_portal_eq_generator v7

Fix:
- Removes invalid "*ENERGY PRINT" (not supported in CCX)
- Uses "*EL PRINT, ELSET=E_ALL" + "ENER" in STEP 1 and STEP 2
  to satisfy CCX nonlinear implicit dynamics requirement.

Run (from ~/CalculiX/tests):
  unzip -o ccx_portal_eq_generator_v7.zip -d ccx_portal_eq_generator_v7
  cd ccx_portal_eq_generator_v7/ccx_portal_eq_generator_v7
  bash run_all.sh

If needed:
  CCX=../../../src/CalculiX bash run_all.sh
