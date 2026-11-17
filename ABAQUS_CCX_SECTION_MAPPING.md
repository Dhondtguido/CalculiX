# Native B31T/B31E section semantics — CCX and Abaqus 2016

The standard CalculiX beam section keywords retain **native CCX**
definitions. The additional I/L profiles use Abaqus 2016 cross-section
library field order. This is not a claim of full Abaqus compatibility.

## Section shapes

- `RECT=a,b`: section 1 width and section 2 depth (native CCX).
- `CIRC=d1,d2`: native CCX ellipse full diameters. **Not** the
  Abaqus `CIRC=r` convention; convert Abaqus radius to `d=2r`.
- `PIPE=R,t`: outer radius and thickness (both compatible).
- `BOX=a,b,t1,t2,t3,t4`: outside sections along 1,2 and right,
  top, left, bottom thicknesses (native CCX).
- `I=l,h,b1,b2,t1,t2,t3`: Abaqus reference height `l` above
  lower flange, total section depth `h`, lower flange width `b1`,
  upper flange width `b2`, lower/upper thicknesses `t1/t2`,
  web thickness `t3`. The line of nodes is at the specified `l`;
  the beam kernel uses the computed centroid offset.
- `L=a,b,t1,t2`: Abaqus horizontal leg length, vertical leg
  length, horizontal leg thickness, vertical leg thickness.
  The section reference lies at the lower-left outer corner;
  the offset to the centroid is accounted for.

Explicit `*USER BEAM SECTION` retains the former non-Abaqus I
(`h,b_top,t_top,b_bottom,t_bottom,t_web`) and L
(`a,b,t_equal`) layouts to avoid silently changing legacy input.

`T` and `CHAN` remain nonstandard extensions rather than
pretending to be documented Abaqus library section names.

## Orientation and offsets

The element tangent is connectivity node 1 -> node 2.
`n1` comes from the explicit *BEAM SECTION normal projected
perpendicular to the tangent, or native default `(0,0,-1)`.
`n2=t x n1`. Always export explicit normals from Mecway.

Native CalculiX `OFFSET1=x`, `OFFSET2=y` are dimensionless
fractions of the section's span in the respective local axis.
The solver translates these to physical centroid offsets,
including the I and L reference position, before rotating offsets
to internal principal axes.

Old vector offsets supported by `*USER BEAM SECTION` are
not native CCX OFFSET1/2: their values must not be interchanged.

## Release and force-output limitations

`*RELEASE` supports Abaqus S1/S2 and M1/M2/T codes for
symmetric sections. The solver currently rotates asymmetric
sections to principal axes before release condensation. Single-
axis M1/M2 releases on such sections therefore need validation
in geometric axes before claiming Abaqus equivalence.

`*USER BEAM OUTPUT` and FRD beam-force results have not
been calibrated against Mecway for units, orientations and signs.
A successful static-equilibrium test does not prove their full
convention parity. Units follow the consistent input system.

This implementation is elastic; it does not add Abaqus
nonlinear cross-section integration, thin-wall warping,
`GENERAL`, `HEX`, `ARBITRARY` or `TRAPEZOID`.

## Verification

`validation/check_abaqus_section_offsets.py` checks four
finite-element axial tip displacements under eccentric reference
loading for I, L, RECT+OFFSET1/2 and BOX+OFFSET1/2 against
independent section geometry and Euler-Bernoulli beam theory.
