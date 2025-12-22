ccx_portal_eq_generator v9 (SOLIDS + rigid-body hinges)

This is the "it should just work" path:
- Avoids rotational SPRING2 DOF 4..6 on beam nodes (often not active in CCX global DOF system).
- Uses rigid-body ROT nodes where U3 == theta_z [rad]. Hinges are SPRING2 on DOF 3 between ROT nodes.

Run:
  unzip -o ccx_portal_eq_generator_v9.zip -d ccx_portal_eq_generator_v9
  cd ccx_portal_eq_generator_v9/ccx_portal_eq_generator_v9

  # auto-detect CCX, or force it:
  CCX=../../../src/CalculiX bash run_all.sh

Outputs:
- portal_eq.inp / .dat / .frd / .sta
- hinge_map.json
- printed: theta_z at ROT nodes and dtheta for left/right hinges
