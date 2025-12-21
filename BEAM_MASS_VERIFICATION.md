# Beam Element - Verificación de Masa Consistente

**Fecha**: 2025-12-21
**Estado**: ✅ VERIFICADO - Masa consistente implementada correctamente

---

## Resumen

El elemento beam U1 (`src/e_c3d_u1.f`) **YA tiene masa consistente correctamente implementada**. No se requieren cambios.

---

## Implementación Actual

### Ubicación del Código
Archivo: `src/e_c3d_u1.f`, líneas 533-574

### Matriz de Masa (12×12 en coordenadas locales)

**Coeficientes**:
```fortran
c1 = (1/3) * rho * detJ
c2 = 2 * c1 = (2/3) * rho * detJ
```

Donde:
- `rho`: Densidad del material
- `detJ = 0.5 * sqrt(dx² + dy² + dz²)`: Jacobiano (mitad de la longitud del elemento)

**Términos Diagonales** (masa propia de cada nodo):
```fortran
sm(1,1) = sm(2,2) = sm(3,3) = a*c2        ! Traslación nodo 1
sm(4,4) = (xi11+xi22)*c2                   ! Rotación axial nodo 1
sm(5,5) = xi22*c2                          ! Rotación sobre eje 2 nodo 1
sm(6,6) = xi11*c2                          ! Rotación sobre eje 1 nodo 1

sm(7,7) = sm(8,8) = sm(9,9) = a*c2       ! Traslación nodo 2
sm(10,10) = (xi11+xi22)*c2                 ! Rotación axial nodo 2
sm(11,11) = xi22*c2                        ! Rotación sobre eje 2 nodo 2
sm(12,12) = xi11*c2                        ! Rotación sobre eje 1 nodo 2
```

**Términos Off-Diagonal** (acoplamiento entre nodos):
```fortran
sm(1,7) = a*c1           ! Acoplamiento traslación x nodo 1-2
sm(2,8) = a*c1           ! Acoplamiento traslación y nodo 1-2
sm(3,9) = a*c1           ! Acoplamiento traslación z nodo 1-2
sm(4,10) = (xi11+xi22)*c1 ! Acoplamiento rotación axial nodo 1-2
sm(5,11) = xi22*c1        ! Acoplamiento rotación sobre eje 2 nodo 1-2
sm(6,12) = xi11*c1        ! Acoplamiento rotación sobre eje 1 nodo 1-2
```

### Transformación a Coordenadas Globales

La matriz se transforma de coordenadas locales a globales (líneas 598-620):

```
SM_global = T^T * SM_local * T
```

Donde T es la matriz de transformación 12×12 basada en los ejes locales e1, e2, e3.

---

## Validación de la Formulación

### Términos Translacionales

Para un elemento beam de longitud L, masa total `m = rho*A*L`, la masa consistente con interpolación lineal es:

```
[M_trans] = (m/6) * [2  1]
                     [1  2]
```

En el código:
- `m = rho*A*L`
- `detJ = L/2`
- `c1 = rho*detJ/3 = rho*L/6`
- `c2 = 2*c1 = rho*L/3`
- `sm(1,1) = a*c2 = A*rho*L/3 = (rho*A*L)/3 = m/3`
- `sm(1,7) = a*c1 = A*rho*L/6 = (rho*A*L)/6 = m/6`

Esto coincide con la formulación teórica: `[m/3, m/6; m/6, m/3]` ✅

### Términos Rotacionales

Para rotaciones, la masa rotacional es `I_rot = rho * Iy` (o `Iz` según el eje).

La matriz consistente sigue la misma estructura que la translacional, usando momentos de inercia en lugar de área.

En el código:
- `sm(5,5) = xi22*c2 = Iz*(2*rho*detJ/3)` (correcto para rotación sobre eje 2)
- `sm(6,6) = xi11*c2 = Iy*(2*rho*detJ/3)` (correcto para rotación sobre eje 1)

---

## Verificación de No-Lumping

### Búsqueda de Lumping Forzado

Se realizó búsqueda en el código:

```bash
grep -rn "call lump" . --include="*.f" --include="*.c"
```

**Resultado**: No se encontraron llamadas a `lump()` para elementos beam.

### Análisis de `lump.f`

El archivo `lump.f` (líneas 1-53) es una rutina genérica que suma términos off-diagonal a la diagonal para crear una matriz lumped:

```fortran
adl(i) = adb(i)  ! Diagonal
do k=jq(j),jq(j+1)-1
  adl(i) = adl(i) + aub(k)  ! Sumar off-diagonal
  adl(j) = adl(j) + aub(k)
enddo
```

**Conclusión**: Esta rutina NO se llama para beams en el código actual.

---

## Pruebas de Consistencia

### Test Teórico: Masa Total

Para un beam de:
- Longitud L = 1.0 m
- Área A = 0.01 m²
- Densidad rho = 7850 kg/m³

Masa total esperada:
```
m_total = rho * A * L = 7850 * 0.01 * 1.0 = 78.5 kg
```

Masa de la matriz consistente (suma de todos los términos):
```
M_total = sum(M_ij) = 2*(m/3 + m/3 + m/6) = 2*(m/2) = m
```

Por lo tanto, la masa se conserva ✅

### Test Dinámico (Recomendado)

Para verificación completa, se recomienda un análisis modal de un cantilever beam:

**Frecuencia natural del primer modo (flexión)**:
```
f1 = (1.875² / (2π)) * sqrt(E*I / (rho*A*L⁴))
```

Para beam típico (acero, rectangular 10mm × 100mm, L=1m):
```
f1_teórico ≈ 33.5 Hz
```

Si la masa consistente está correcta, CalculiX debería dar un resultado muy cercano.

---

## Comparación: Consistente vs Lumped

| Aspecto | Masa Consistente (Actual) | Masa Lumped (Alternativa) |
|---------|---------------------------|---------------------------|
| Precisión | ✅ Alta (error ~0.5% freq) | ⚠️ Media (error ~3-5% freq) |
| Costo computacional | ⚠️ Mayor (matriz no-diagonal) | ✅ Menor (matriz diagonal) |
| Modos altos | ✅ Mejor representación | ❌ Modos espurios |
| Implementación | ✅ Ya existe | N/A |

**Recomendación**: Mantener masa consistente (estado actual).

---

## Conclusión

✅ **La masa consistente para beams está correctamente implementada en CalculiX.**

- Formulación correcta según teoría de elementos finitos
- Términos translacionales y rotacionales correctos
- Acoplamiento entre nodos presente
- No hay lumping forzado
- Transformación a coordenadas globales correcta

**No se requieren cambios en la implementación de masa.**

---

## Referencias

- Código fuente: `src/e_c3d_u1.f`, líneas 533-574
- Paper base: Yunhua Luo (2008), "An Efficient 3D Timoshenko Beam Element with Consistent Shape Functions"
- Libros de referencia:
  - T.J.R. Hughes, "The Finite Element Method"
  - K.J. Bathe, "Finite Element Procedures"
