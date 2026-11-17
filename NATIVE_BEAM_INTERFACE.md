# Native B31T / B31E Interface

This development branch keeps the historical CalculiX `B31` element unchanged and adds two native-style direct beam names backed by the UB31 kernel:

- `B31T` — 2-node 3D Timoshenko beam (shear deformation enabled using section shear factors).
- `B31E` — 2-node 3D Euler-Bernoulli beam (transverse shear deformation suppressed).
- `UB31` — the existing explicit user-element interface remains available for development and backward compatibility.

## Input syntax

```text
*ELEMENT, TYPE=B31T, ELSET=EBEAM
1, 1, 2

*MATERIAL, NAME=STEEL
*ELASTIC
2.1E11, 0.3

*BEAM SECTION, ELSET=EBEAM, MATERIAL=STEEL, SECTION=RECT
0.1, 0.2
0.0, 0.0, 1.0
```

Change `TYPE=B31T` to `TYPE=B31E` to select Euler-Bernoulli behavior. No `*USER ELEMENT` or `*USER BEAM SECTION` card is required for these two native-style names.

## Backward compatibility

`TYPE=B31` is deliberately **not** redirected. It remains the original CalculiX expanded beam element so existing input files retain their historical behavior.

Internally, the parser maps:

```text
B31T -> UB31T -> UB31 kernel (Timoshenko)
B31E -> UB31E -> UB31 kernel (Euler-Bernoulli)
```

The internal `T`/`E` marker is retained in `lakon(5:5)` and the UB31 kernel uses it to select the formulation. Plain `UB31` continues to use its existing section-card setting.

## Current scope

The native interface reuses the current UB31 section reader and therefore supports the section families already implemented there (`RECT`, `CIRC`, `PIPE`, `I`, `T`, `CHAN`, `L`, `BOX`). Native `I` and `L` geometry follows Abaqus 2016; `*USER BEAM SECTION` keeps the older I/L parameters. Abaqus T sections can be represented by degenerate native `I` profiles. `SECTION=GENERAL` is not yet claimed as supported by this interface. See [CCX_BEAM_CONVENTIONS.md](CCX_BEAM_CONVENTIONS.md) for native geometric data, axis rules, and limitations.

### Mecway open-section compatibility

Native `*BEAM SECTION, SECTION=BOX` now uses standard CalculiX field
ordering. A BOX section with one missing wall remains an open BOX:
it is no longer silently remapped to `CHAN`. Older Mecway input
that uses the custom `h,b,t_bottom,t_left,t_top,t_right` ordering
must be converted explicitly by the exporter or bridge before
the solver is run. Do not reuse legacy BOX data without mapping.

The interface still requires rebuilt CCX numerical regression
before production use.
