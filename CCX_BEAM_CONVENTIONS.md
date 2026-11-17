# CCX-compatible section geometry and orientation (B31T/B31E)

The original CalculiX `B31` expanded beam is unchanged.
Native `B31T` and `B31E` use `*BEAM SECTION`. The explicit
`UB31`/`*USER BEAM SECTION` interface retains its legacy fields.

## Section axes

The beam tangent `t` runs from connectivity node 1 to 2.
The third line of `*BEAM SECTION` supplies section local axis
`n1` as a global vector. The solver projects `n1` perpendicular
to `t`, normalises it, and calculates `n2 = t x n1`.
Internally, `e2=n1` and `e3=n2`. Omitted n1 defaults to
`(0,0,-1)` (CalculiX convention). An explicit zero vector or an
orientation parallel to the beam tangent is rejected; specify a
valid normal instead. The direct beam does not implement all knot
and `*NORMAL` semantics of native expanded `B31`.

## Native beam section geometric data

| Shape | Fields | Definition |
|---|---|---|
| RECT | `a,b` | Dimensions along n1 and n2, b defaults to a |
| CIRC | `d1,d2` | Ellipse diameters along n1/n2, d2 defaults to d1 |
| PIPE | `R,t` | Outer radius and wall thickness |
| BOX | `a,b,t1,t2,t3,t4` | Outer width along n1, height along n2; right, top, left, bottom thicknesses |

Native BOX geometry is converted internally to the prior
`(h,width,bottom,left,top,right)` shape routine. Open BOX with
a zero-thickness wall is retained as BOX. Older Mecway-style BOX
cards with the bespoke legacy field order must be converted
explicitly before analysis; the solver no longer silently infers
channel geometry from a BOX card.

Native CIRC includes ellipses and uses exact elastic section A,
second moments and Saint-Venant J; the elliptical Timoshenko
shear-correction factor remains an approximation.

`I`, `L`, `T` and `CHAN` remain custom UB31 section
extensions and do **not** claim Abaqus's I/L parameter ordering.

## Further verification required

- Principal-axis rotation in asymmetric sections means the internal
  solver axes need explicit validation against end-release conventions.
- Original CCX `SECTION FORCES`, internal force/moment units, sign
  conventions and CGX export are not yet certified as equivalent.
- Cross-section integration with nonlinear material behavior,
  CCX knot behavior, thin-walled warping and full Abaqus beam
  library support remain outside this change.

## Verification

Run validation decks `22`–`27` with the revised executable,
then invoke `python validation/check_ccx_beam_conventions.py`
from the output directory. Six analytical tip displacements are
checked with a 3% tolerance. Full runtime validation is required
before use on safety-critical engineering projects.
