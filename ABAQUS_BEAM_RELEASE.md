# Abaqus-compatible beam end releases for B31T/B31E

## Scope

This development branch adds the **Abaqus/Standard `*RELEASE`**
model-data keyword for native direct beams (`B31T`, `B31E`) and
explicit UB31 user beams. Existing `*USER BEAM RELEASE` syntax remains
available; legacy CalculiX `B31` is **not** converted or altered.

The parser routes `*RELEASE` to the existing element end-release
bitmask and static-condensation machinery, not to new stiffness code.

## Abaqus form

```inp
*ELEMENT, TYPE=B31T, ELSET=FRAME
1, 1, 2
2, 2, 3
*RELEASE
FRAME, S1, M1-M2
2, S2, ALLM
```

- Field 1: existing element number or ELSET name.
- Field 2: `S1` = first connectivity node, `S2` = second.
- Field 3: `M1`, `M2`, `T`, `M1-M2`, `M1-T`, `M2-T`, `ALLM`.
- Rotations: `M1` = local bending about beam local y (rotation DOF 5);
  `M2` = local bending about beam local z (rotation DOF 6);
  `T` = torsion about beam local x (rotation DOF 4).
- The `*RELEASE` parser accepts **only** the seven Abaqus combinations
  above and only `S1`/`S2`. It rejects numerical bitmasks and
  translational releases. `INPUT=` is not implemented.

Place the keyword after element/ELSET definitions and before `*STEP`.
Unreleased DOFs remain connected. Do not fully release every member
rotation connected to an otherwise unrestrained node: a free rotational
DOF can create a singular system. Support boundary conditions must
suppress any unused nodal rotations.

## Validation

Use `validation/20_native_b31t_abaqus_release.inp`.
A simply supported, 2 m long beam carries 1000 N at midspan.
Expected |Rz| at each end = 500 N and zero end bending moments for
the released bending axes. Midspan deflection can be compared with
`PL^3/(48EI)` plus the Timoshenko shear term `PL/(4 k G A)`.
The deck exercises both an ELSET target and an individual element ID.

**Status:** source implementation and test input only. A rebuilt CCX
binary and a successful runtime regression are required before relying
on the functionality in engineering production.

## Beam offsets are a separate feature

Existing UB31 code also supports *nonstandard* section-header parameters:
`OFFSET=(y,z)` for a uniform local y/z shift and
`OFFSET1=(x,y,z)` / `OFFSET2=(x,y,z)` for independent local nodal
eccentricities. These are interpreted by `src/userbeamsections.f` and
`src/e_c3d_ub31.f` (using rigid offset kinematic coupling).

They are **not** standard Abaqus `*BEAM SECTION OFFSET` syntax.
In Abaqus, `*BEAM SECTION OFFSET` is associated with
`*BEAM GENERAL SECTION` and specifies a cross-section origin shift.
It is not interchangeable with independently varying member-end
eccentricities. No change to offset syntax is included in this branch.

## Optional Mecway release-axis compatibility (diagnostic)

Some Mecway-generated B31T input decks have been observed to export the
single-axis `rotv` and `rotw` releases using the opposite `M1`/`M2`
labels relative to the local basis supplied in the `*BEAM SECTION` card.
This is currently a **hypothesis awaiting a rerun against Mecway results**,
not a verified general Mecway bug.

To test without manually changing the exported `.inp`, enable this opt-in
environment setting **before starting the CCX process**:

```cmd
set CCX_MECWAY_SWAP_RELEASE_AXES=1
ccx_2.23.exe testb31t
```

For Mecway started through the Windows GUI, set a persistent **user**
environment variable once using PowerShell, then close and reopen Mecway:

```powershell
[Environment]::SetEnvironmentVariable('CCX_MECWAY_SWAP_RELEASE_AXES', '1', 'User')
```

The setting swaps `M1` and `M2` (including `M1-T` and `M2-T`)
**only** for standard `*RELEASE` lines. `T`, `M1-M2`, and `ALLM`
remain unchanged. The `*USER BEAM RELEASE` extension is unaffected.

**Important:** With this setting enabled, CCX intentionally interprets
`*RELEASE` differently from standard Abaqus. Do **not** enable it for
ordinary Abaqus decks or unrelated models. With the variable absent
or set to another value, the original Abaqus-compatible behaviour is retained.
To turn it off for the Windows GUI, remove the user environment variable:

```powershell
[Environment]::SetEnvironmentVariable('CCX_MECWAY_SWAP_RELEASE_AXES', $null, 'User')
```

Reopen Mecway after either environment change. This workaround does **not**
convert legacy `B31` to `B31T` or fix missing Mecway B31T export support.
Mecway's exporter remains the proper place for a permanent mapping fix.

Regression: `validation/21_mecway_release_axis_swap.inp` is run with the
setting both disabled and enabled in the Windows build workflow. Compare the
resulting displacement fields to confirm the opt-in changes the released
bending direction.
