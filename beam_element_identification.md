# Beam element identification (current codebase)

## Element types and mapping
- Input labels `B31/B31R` and `B32/B32R` are accepted by the *ELEMENT reader (`src/elements.f` lines 130–193), including `B21 -> B31` and `T2D2 -> T3D2` remapping for 2D beams/trusses.
- Allocation records beam element sizes and expansion info (`src/allocation.f` lines 1036–1058), noting that `B31` expands to a C3D8I-style topology (11-node expansion) and `B32` expands to C3D20-style topology.
- During `gen3delem`, beam labels are converted into 3D solid labels with a beam flag in the 7th character (`src/gen3delem.f` lines 378–390):
  - `B31`/`T3D2` → `C3D8I B`
  - `B31R` → `C3D8R B`
  - `B32`/`T3D3` → `C3D20 B`
  - `B32R` → `C3D20RB`

## DOFs
- After remapping, the beam elements are treated as 3D elements with 3 DOFs per node (see `C3D8I`/`C3D20` handling in `src/mastructcs.c` lines 205–213).

## Stiffness and mass assembly
- Element assembly for non-`U*` elements (which includes the remapped beam elements) calls the standard 3D element routine `e_c3d` (`src/mafillsmcs.f` lines 204–218).
- `e_c3d` handles `C3D8I` and `C3D20` element topology selection and uses the standard 3D element path (`src/e_c3d.f` lines 134–176). This is where the stiffness/mass for the remapped beam elements is currently computed.

## Section forces / output pipeline
- Section forces for beam elements are assembled in `map3dto1d2d` (see `src/map3dto1d2d.f` lines 313–489). This routine integrates stresses over the beam cross section to compute forces/moments.
- The *SECTION PRINT (`SOF`/`SOM`) pipeline runs through `resultsprint.f`, which calls `extrapolate` (for nodal stresses) and then `printoutface` for the actual output (`src/resultsprint.f` lines 125–153).
- `printoutface` emits the section output headers and iterates over the requested section set (`src/printoutface.f` lines 72–140).

## Beam section property handling
- `*BEAM SECTION` parsing assigns materials, orientations, and section shape to beam elements (`src/beamsections.f` lines 33–217), while `*BEAM SECTION, SECTION=GENERAL` is restricted to `U1` elements (`src/beamgeneralsections.f` lines 234–276).
