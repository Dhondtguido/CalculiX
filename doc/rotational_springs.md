# Rotational SPRING2 support

CalculiX now supports **rotational** springs for `*SPRING` definitions used with
`SPRING2` elements. Degrees of freedom **4–6** (nodal rotations) are accepted
in addition to the existing translational DOFs **1–3**.

## Key points

- **Supported element:** `SPRING2` (two-node spring elements).
- **DOF range:** `1..6`. Rotational DOFs are `4..6`.
- **Units:** rotation is in **radians** and the spring force is a **moment** (N·m).
- **Active DOFs required:** rotational DOFs must be active at the involved nodes
  (typically introduced by beam-type elements). If not, CalculiX aborts with:
  `Rotational spring DOF requested but rotational DOF not active for these nodes.`
- **Nonlinear tables:** the M–θ table is handled exactly like translational spring
  tables (piecewise-linear interpolation and tangents).

## Minimal example (linear)

```inp
*ELEMENT, TYPE=SPRING2, ELSET=E_SPRING
1, 1, 2
*SPRING, ELSET=E_SPRING
6, 6
1000.
```
