# Beam Element Validation Tests

**Created**: 2025-12-21
**Purpose**: Validate beam element formulation, deflections, and check for shear locking

---

## Test Suite

### Test 1: Cantilever Beam (`test1_cantilever.inp`)

**Description**: Classic cantilever beam with tip load

**Geometry**:
- Length: L = 1.0 m
- Section: Rectangular 0.01 m × 0.1 m (b × h)
- Elements: 5 × B31 (2-node beam elements)

**Material**:
- Steel: E = 210 GPa, ν = 0.3

**Loading**:
- Tip load: P = 1000 N (downward, -z direction)

**Boundary Conditions**:
- Fixed support at x = 0 (all 6 DOFs constrained)

**Analytical Solution** (Euler-Bernoulli):
```
I = b*h³/12 = 0.01 * 0.1³ / 12 = 8.333×10⁻⁷ m⁴
δ_tip = P*L³/(3*E*I) = 1000 * 1.0³ / (3 * 210×10⁹ * 8.333×10⁻⁷)
δ_tip = 0.001905 m = 1.905 mm
```

**Expected CalculiX Result**:
- Tip displacement (node 6, U3): ≈ 1.90-1.92 mm
- Error: < 2% (due to Timoshenko vs EB difference)

**How to Run**:
```bash
../../src/CalculiX -i test1_cantilever
grep "node.*6" test1_cantilever.dat
```

**Pass Criteria**:
- U3 at node 6: 1.85 mm < U3 < 1.95 mm

---

### Test 2: Deep Beam (`test2_deepbeam.inp`)

**Description**: Short, deep beam to check shear deformation and locking

**Geometry**:
- Length: L = 0.5 m
- Section: Rectangular 0.05 m × 0.2 m (b × h)
- **L/h ratio**: 0.5/0.2 = 2.5 (deep beam!)
- Elements: 4 × B31

**Material**:
- Steel: E = 210 GPa, ν = 0.3, G = 80.77 GPa

**Loading**:
- Midspan load: P = 10000 N (downward)

**Boundary Conditions**:
- Simple supports at ends (vertical displacement constrained)

**Analytical Solution** (Timoshenko):
```
I = 0.05 * 0.2³ / 12 = 3.333×10⁻⁵ m⁴
A = 0.05 * 0.2 = 0.01 m²
κ = 5/6 (shear coefficient for rectangular section)

δ_flexure = P*L³/(48*E*I) = 10000 * 0.5³ / (48 * 210×10⁹ * 3.333×10⁻⁵)
          = 3.7×10⁻⁵ m = 0.037 mm

δ_shear = P*L/(4*κ*G*A) = 10000 * 0.5 / (4 * 5/6 * 80.77×10⁹ * 0.01)
        = 1.86×10⁻⁵ m = 0.0186 mm

δ_total = δ_flexure + δ_shear = 0.0556 mm
```

**Expected CalculiX Result**:
- Midspan displacement (node 3, U3): ≈ 0.055-0.057 mm

**How to Run**:
```bash
../../src/CalculiX -i test2_deepbeam
grep "node.*3" test2_deepbeam.dat
```

**Pass Criteria**:
- U3 at node 3: 0.050 mm < U3 < 0.060 mm
- **Shear locking indicator**: If U3 << 0.050 mm, beam is locking!

**Why This Test Matters**:
- With L/h = 2.5, shear deformation is ≈ 33% of total deflection
- If formulation locks, shear will be underestimated → deflection too small
- A good Timoshenko formulation should capture both flexure and shear

---

## Validation Procedure

### Manual Validation

1. **Compile CalculiX** (if not already done):
   ```bash
   cd ../../src
   make
   ```

2. **Run Test 1**:
   ```bash
   cd ../tests/beams
   ../../src/CalculiX -i test1_cantilever
   ```

3. **Check Results**:
   ```bash
   # Look for node 6 displacement in output
   grep -A5 "displacements" test1_cantilever.dat
   # or
   tail -100 test1_cantilever.dat | grep -A2 "node.*6"
   ```

4. **Compare with Analytical**:
   - Test 1: δ_analytical = 1.905 mm
   - Test 2: δ_analytical = 0.0556 mm

5. **Repeat for Test 2**:
   ```bash
   ../../src/CalculiX -i test2_deepbeam
   grep -A5 "displacements" test2_deepbeam.dat
   ```

### Automated Validation (TODO)

Create a Python/Bash script `validate_tests.sh`:

```bash
#!/bin/bash

echo "Running Beam Element Validation Tests..."
echo "=========================================="

# Test 1: Cantilever
echo ""
echo "Test 1: Cantilever Beam"
../../src/CalculiX -i test1_cantilever > /dev/null 2>&1
DISP1=$(grep "node.*6" test1_cantilever.dat | awk '{print $4}')
echo "  Tip displacement: $DISP1 mm (expected: 1.905 mm)"

# Test 2: Deep Beam
echo ""
echo "Test 2: Deep Beam"
../../src/CalculiX -i test2_deepbeam > /dev/null 2>&1
DISP2=$(grep "node.*3" test2_deepbeam.dat | awk '{print $4}')
echo "  Midspan displacement: $DISP2 mm (expected: 0.0556 mm)"

# Check pass/fail
# (Add numerical comparison here)
```

---

## Expected Output Example

### Beam Section Forces

With the improved `resultsmech_u1.f`, you should see output like:

```
Beam element    1  section forces (local coords):
  Node 1: N, Vy, Vz, T, My, Mz =  0.00000E+00  5.00000E+02  0.00000E+00  0.00000E+00  2.50000E+02  0.00000E+00
  Node 2: N, Vy, Vz, T, My, Mz =  0.00000E+00  5.00000E+02  0.00000E+00  0.00000E+00  1.50000E+02  0.00000E+00
```

**Interpretation**:
- `N`: Axial force
- `Vy`, `Vz`: Shear forces
- `T`: Torsion
- `My`, `Mz`: Bending moments

---

## Known Issues / Limitations

1. **No Compiler Available**:
   - This test suite requires compiling CalculiX
   - If gfortran not available, tests cannot run

2. **Shear Locking (To Be Addressed)**:
   - Current formulation may exhibit shear locking for very deep beams
   - Test 2 will reveal if this is an issue
   - If locking detected, implement SRI (Selective Reduced Integration)

3. **Output Format**:
   - CalculiX .dat format may vary by version
   - Grep patterns may need adjustment

---

## Interpreting Results

### Good Results

✅ Test 1: δ_tip ≈ 1.90 mm (within 2%)
✅ Test 2: δ_mid ≈ 0.055 mm (within 10%)

### Potential Issues

❌ **Test 1 fails**: Check:
- Material properties (E, I)
- Boundary conditions (fixed = all 6 DOFs)
- Load application

❌ **Test 2 significantly underestimates** (e.g., δ < 0.04 mm):
- **Shear locking detected!**
- Formulation needs SRI or assumed shear strain

❌ **Test 2 significantly overestimates** (e.g., δ > 0.07 mm):
- Check shear coefficient (κ should be 5/6 for rectangular)
- Check material G (should be E/(2*(1+ν)))

---

## Next Steps

1. Run tests once CalculiX is compiled
2. Document results
3. If shear locking detected in Test 2:
   - Implement SRI in `e_c3d_u1.f`
   - Retest
4. If results pass:
   - Add more comprehensive tests (modal, dynamic)
   - Consider P-Delta test for NLGEOM

---

## References

- Timoshenko beam theory: Timoshenko & Gere, "Mechanics of Materials"
- Shear locking: T.J.R. Hughes, "The Finite Element Method"
- CalculiX documentation: http://www.dhondt.de/
