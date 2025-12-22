# CalculiX Beam Element - Análisis e Identificación

**Fecha**: 2025-12-21
**Objetivo**: Mejorar elementos beam tipo Abaqus (B31/B32) en este fork de CalculiX

---

## 1. IDENTIFICACIÓN DEL ELEMENTO BEAM ACTUAL

### 1.1 Tipo de Elemento
- **Nombre interno**: User Element U1
- **Archivo principal**: `src/e_c3d_u1.f` (líneas 19-650)
- **Archivo de resultados**: `src/resultsmech_u1.f` (líneas 19-648)
- **Referencia**: Yunhua Luo, "An Efficient 3D Timoshenko Beam Element with Consistent Shape Functions", Adv. Theor. Appl. Mech., Vol. 1, 2008, no. 3, 95-106

### 1.2 Tipos de Beam Soportados
Encontrados en `src/elements.f` y `src/allocation.f`:
- **B31**: 2 nodos, lineal
- **B32**: 3 nodos, cuadrático
- **B31R**: 2 nodos, reduced integration
- **B32R**: 3 nodos, reduced integration

### 1.3 Características del Elemento U1

#### Grados de Libertad
- **2 nodos** (nope=2, línea 279 en e_c3d_u1.f)
- **6 DOF por nodo** (ndof=6, línea 280): 3 traslaciones (u,v,w) + 3 rotaciones (φ,ψ,θ)
- **Total**: 12 DOF por elemento

#### Propiedades de Sección (líneas 268-277 en e_c3d_u1.f)
Almacenadas en `prop(index+1..10)`:
1. `a`: Área de la sección transversal
2. `xi11`: Momento de inercia sobre eje local 1 (Iy)
3. `xi12`: Momento de inercia mixto (debe ser 0)
4. `xi22`: Momento de inercia sobre eje local 2 (Iz)
5. `xk`: **Coeficiente de corte Timoshenko** (kappa)
6-8. `e2(1:3)`: Vector dirección local 2 (dado por usuario como "e1" en input)
9. `offset1`: Offset en dirección 1 (debe ser 0 para U1)
10. `offset2`: Offset en dirección 2 (debe ser 0 para U1)

#### Sistema de Coordenadas Local
- **e1**: Paralelo al eje de la viga (xl(j,2) - xl(j,1))
- **e2**: Dado por usuario en *BEAM SECTION (vector normal)
- **e3**: e1 × e2 (producto cruz)
- **Transformación**: Matriz tm(3,3) de global a local (líneas 333-337)

---

## 2. FORMULACIÓN ACTUAL

### 2.1 Stiffness Matrix (K)

**Ubicación**: `src/e_c3d_u1.f`, líneas 479-529

**Tipo**: Timoshenko beam element

**Términos intermedios** (líneas 479-488):
```fortran
y1 = xk*um*a*e*xi11*(12.d0*e*xi11 + xk*um*a*dl*dl)
y2 = (12.d0*e*xi11 - xk*um*a*dl*dl)**2
y3 = 4.d0*e*xi11*((xk*um*a)**2*dl**4 + 3.d0*xk*um*a*dl*dl*e*xi11 + 36.d0*(e*xi11)**2)

z1 = xk*um*a*e*xi22*(12.d0*e*xi22 + xk*um*a*dl*dl)
z2 = (12.d0*e*xi22 - xk*um*a*dl*dl)**2
z3 = 4.d0*e*xi22*((xk*um*a)**2*dl**4 + 3.d0*xk*um*a*dl*dl*e*xi22 + 36.d0*(e*xi22)**2)
```

Donde:
- `e`: Módulo de Young
- `um`: Módulo de corte G = E/(2(1+ν))
- `xk`: Coeficiente de corte Timoshenko (kappa)
- `a`: Área
- `dl`: Longitud del elemento
- `xi11`, `xi22`: Momentos de inercia

**Matriz de rigidez local S'** (12×12, líneas 499-529):
- Términos axiales (s(1,1), s(7,7)): E*A/L
- Términos de flexión en plano e2 (con corte): Usa y1, y2, y3
- Términos de flexión en plano e3 (con corte): Usa z1, z2, z3
- Término torsional (s(4,4)): um*(xi11+xi22)/dl

**Transformación a global**: S = T^T * S' * T (líneas 625-643)

### 2.2 Mass Matrix (M)

**Ubicación**: `src/e_c3d_u1.f`, líneas 533-574

**Tipo**: **MASA CONSISTENTE** (ya implementada!)

**Matriz de masa local SM'** (12×12, líneas 543-563):
```fortran
c1 = (1/3) * rho * detJ
c2 = 2 * c1

Diagonal:
  sm(1,1) = sm(2,2) = sm(3,3) = a*c2  (traslación)
  sm(4,4) = (xi11+xi22)*c2              (rotación axial)
  sm(5,5) = xi22*c2                     (rotación)
  sm(6,6) = xi11*c2                     (rotación)
  sm(7,7) = sm(8,8) = sm(9,9) = a*c2  (traslación nodo 2)
  sm(10,10) = (xi11+xi22)*c2            (rotación nodo 2)
  sm(11,11) = xi22*c2
  sm(12,12) = xi11*c2

Off-diagonal (acoplamiento entre nodos):
  sm(1,7) = a*c1
  sm(2,8) = a*c1
  sm(3,9) = a*c1
  sm(4,10) = (xi11+xi22)*c1
  sm(5,11) = xi22*c1
  sm(6,12) = xi11*c1
```

**Transformación a global**: SM = T^T * SM' * T (líneas 602-620)

### 2.3 Recuperación de Fuerzas Internas

**Ubicación**: `src/resultsmech_u1.f`, líneas 497-507

**Section forces en sistema local** (2 integration points = 2 nodos):
```fortran
stre(1) = e*a*emec(1)              ! N:  Fuerza axial
stre(2) = e*xi11*emec(2)           ! M2: Momento sobre eje 2 (My)
stre(3) = e*xi22*emec(3)           ! M1: Momento sobre eje 1 (Mz)
stre(4) = xk*um*a*emec(4)          ! Q1: Corte en dirección 1 (Vy)
stre(5) = xk*um*a*emec(5)          ! Q2: Corte en dirección 2 (Vz)
stre(6) = xk*um*(xi11+xi22)*emec(6)! T:  Torsión
```

Almacenadas en `stx(1:6, jj, i)` donde jj=1,2 (nodos 1 y 2).

---

## 3. PROBLEMAS IDENTIFICADOS

### 3.1 Shear Locking ⚠️
**Problema**: La formulación actual usa integración completa para los términos de corte.

**Evidencia**:
- Los términos y1, y2, y3, z1, z2, z3 incluyen xk*um*a (deformación de corte)
- No hay integración selectiva reducida (SRI) visible
- Para vigas esbeltas (L/h grande), el término de corte domina incorrectamente

**Solución necesaria**: Implementar Selective Reduced Integration (SRI):
- Integración completa (2 puntos Gauss) para flexión (E*I)
- Integración reducida (1 punto Gauss) para corte (kappa*G*A)

### 3.2 Output de Fuerzas Internas 📊
**Problema**: No hay output estructurado para end forces.

**Situación actual**:
- `stx(1:6, jj, i)` contiene N, M2, M1, Q1, Q2, T en jj=1,2
- NO se imprimen de forma clara para el usuario
- `sectionprints.f` maneja *SECTION PRINT pero no está conectado a beams

**Solución necesaria**:
- Extender `sectionprints.f` o crear output específico para beams
- Formato: Element_ID, Node1: N Vy Vz T My Mz, Node2: N Vy Vz T My Mz

### 3.3 P-Delta / Geometric Stiffness 🔄
**Problema**: No implementado para beams.

**Búsqueda**:
- `nonlingeo.c` existe pero beams no participan en NLGEOM
- e_c3d_u1.f rechaza iperturb (líneas 351-356): "no second order calculation"
- e_c3d_u1.f rechaza buckling (líneas 462-466)

**Solución necesaria**:
- Implementar matriz de rigidez geométrica K_geo
- K_geo basada en N axial (fuerza de compresión/tensión)
- Añadir a K material: K_total = K_mat + K_geo

### 3.4 Masa Consistente ✅
**Estado**: **YA IMPLEMENTADA CORRECTAMENTE**

**Verificación**:
- Masa consistente en líneas 543-563 de e_c3d_u1.f
- Incluye translación y rotación
- Transformada a global correctamente

**Acción**: Solo verificar que se usa (no hay lumping forzado)

---

## 4. ARCHIVOS RELACIONADOS

### 4.1 Elemento Beam
- `e_c3d_u1.f`: Stiffness y mass matrix del elemento U1
- `resultsmech_u1.f`: Cálculo de strains, stresses, section forces
- `e_c3d_u.f`, `e_c3d_us3.f`: Otros user elements (pueden tener Timoshenko también)

### 4.2 Secciones Beam
- `beamsections.f`: Lee *BEAM SECTION del input deck
- `beamgeneralsections.f`: Secciones generales (PIPE, BOX, etc.)
- `beamintscheme.f`: Esquema de integración para secciones no-rectangulares
- `beamextscheme.f`: Extrapolación

### 4.3 Infraestructura
- `beammpc.f`: MPCs para beams
- `calcmass.f`: Cálculo de masa (llama beamintscheme si es beam)
- `masss.f`: Lee *MASS del input
- `lump.f`: Lumping (evitar para beams)
- `sectionprints.f`: Output de *SECTION PRINT (extender para beams)
- `resultsforc*.f`: Output de fuerzas

### 4.4 Expansión 3D
- `gen3delem.f`: Expande beams 1D a elementos 3D para visualización
- `gen3dmpc.f`, `gen3dboun.f`: MPCs y boundary conditions para beams
- `map3dto1d2d.f`: Mapeo de resultados (incluye section forces para beams, líneas 313, 458, 661)

---

## 5. INPUTS Y COMPATIBILIDAD

### 5.1 Input Deck Típico
```
*BEAM SECTION, ELSET=beams, MATERIAL=steel, SECTION=RECT
height, width

*ORIENTATION, NAME=beam_ori
direction_x, direction_y, direction_z
```

### 5.2 Propiedades Calculadas
- En `gen3delem.f` se calculan A, Iy, Iz, J basadas en SECTION (CIRC, RECT, PIPE, BOX)
- Para RECT: A=b*h, Iy=b*h³/12, Iz=h*b³/12, J=función de b y h
- Para CIRC: A=π*r², I=π*r⁴/4, J=π*r⁴/2
- **xk** (kappa): Debe calcularse según sección (típicamente 5/6 para RECT, ~0.9 para CIRC)

### 5.3 Compatibilidad a Mantener
- *BEAM SECTION con MATERIAL, ORIENTATION, SECTION
- Secciones: CIRC, RECT, PIPE, BOX, general
- Offsets (actualmente no soportados en U1)
- Nodal thickness (thicke con valor -1.0)

---

## 6. CÁLCULO DE STIFFNESS Y MASS

### 6.1 Dónde se Calcula Ke
**Rutina**: `e_c3d_u1` (src/e_c3d_u1.f)
**Cuándo**: Cuando `stiffness=1` (línea 471)
**Llamada desde**:
- `mafillsm.f` (para static/dynamic stiffness assembly)
- `results.c` y derivados (para recalcular)

### 6.2 Dónde se Calcula Me
**Rutina**: `e_c3d_u1` (src/e_c3d_u1.f)
**Cuándo**: Cuando `mass=1` (línea 533)
**Llamada desde**:
- `mafillsm.f` (mass matrix assembly)

### 6.3 Dónde se Imprimen Outputs
**Stresses/Forces**:
- `resultsmech_u1` calcula stx(1:6, 1:2, nelem)
- Almacenados en arrays globales
- Output en `.dat` via `results.c`

**Section forces**:
- `map3dto1d2d.f` tiene lógica para section forces de beams (líneas 313, 458, 661)
- NO está completamente implementado el output estructurado

---

## 7. PLAN DE MEJORA

### 7.1 Prioridad Alta
1. **Shear locking**: Implementar SRI en e_c3d_u1.f
2. **Output forces**: Extender sectionprints.f para beams (N, Vy, Vz, T, My, Mz en ambos extremos)

### 7.2 Prioridad Media
3. **Masa consistente**: Verificar que no se fuerza lumping (calcmass.f, lump.f)

### 7.3 Prioridad Baja (Opcional)
4. **P-Delta**: Implementar K_geo si NLGEOM aplica

### 7.4 Tests Necesarios
- Cantilever EB (P en punta): δ = PL³/(3EI), θ = PL²/(2EI)
- Simply supported (carga distribuida w): δ_max = 5wL⁴/(384EI)
- Deep beam (L/h pequeño): Verificar NO locking
- Modal analysis: 1er modo cantilever

---

## 8. NOTAS TÉCNICAS

### 8.1 Correcciones en Paper de Luo
El código tiene comentarios (líneas 492-497 en e_c3d_u1.f) indicando errores en el paper original:
1. Línea 4 del Apéndice A: k(3,3)=k(3,9) debería ser k(3,3)=-k(3,9)
2. Línea 6 del Apéndice A: k(4,4)=k(7,7) debería ser k(4,4)=k(10,10)

### 8.2 Limitaciones Actuales del U1
- No soporta xi12 ≠ 0 (momento de inercia mixto)
- No soporta offset1, offset2 ≠ 0
- No soporta orientación distinta de la dada (iorien debe ser 0 o usar props)
- No soporta thermal analysis (ithermal ≥ 2)
- No soporta body forces
- No soporta coriolis
- No soporta buckling
- No soporta second order (iperturb)
- No soporta initial strains (iprestr)

---

**Fin del análisis de identificación**
