# CalculiX Benchmark - Problema 4 (Tarea3_DC_2015)

Implementación del **Problema 4** de la Tarea 3 (DC 2015): Análisis cíclico de pórtico de 1 piso con rótulas plásticas bajo excitación sísmica senoidal.

## Especificaciones del Problema

### Geometría
- **Altura columnas**: H = 3.0 m
- **Luz del pórtico**: L = 5.0 m

### Secciones transversales (del Problema 2)

**Columnas (S1)**:
- b × h = 40 cm × 60 cm
- 4 barras φ=20mm en ambas caras
- Hormigón: σy = 4.2 tonf/cm² (≈42 MPa)

**Vigas (S2)**:
- b × h = 25 cm × 50 cm
- 3 barras φ=20mm cara superior
- 2 barras φ=20mm cara inferior

### Excitación sísmica
Señal senoidal:
```
a_g(t) = A · cos(0.2π·t) · sin(4π·t)
```

- **Duración**: t ∈ [0, 10] segundos
- **Amplitudes**: A = 0.1g, 0.2g, 0.3g, ... (hasta colapso)

### Modelo
- **Rótulas plásticas**: Base de columnas + extremos de vigas
- **Columnas**: Interacción M vs N, comportamiento elasto-plástico ideal (sin degradación)
- **Vigas**: Solo flexión, modelo SHM con degradación de rigidez/resistencia
- **Masa**: Ajustada para período fundamental T₀ = 0.5 seg

## Quick Start

### 1. Generar señal sísmica

```bash
cd /work/tests/portal4_compare

# Generar earthquake para A=0.1g
python3 generate_earthquake_senoidal.py --A 0.1 --out earthquake_A0.1g.csv
```

**Parámetros**:
- `--A`: Amplitud en unidades de g (e.g., 0.1, 0.2, 0.3)
- `--dt`: Paso de tiempo (default: 0.0025 s)
- `--duration`: Duración (default: 10.0 s)
- `--units`: Salida en 'g' o 'mps2' (default: g)

### 2. Generar input de CalculiX

```bash
python3 generate_portal_eq_ccx.py \
    --eq earthquake_A0.1g.csv \
    --out portal_A0.1g.inp \
    --H 3.0 \
    --L 5.0 \
    --col_b 0.40 \
    --col_h 0.60 \
    --beam_b 0.25 \
    --beam_h 0.50
```

**Parámetros clave**:
- `--H`, `--L`: Geometría del pórtico (m)
- `--col_b`, `--col_h`: Sección columna (m)
- `--beam_b`, `--beam_h`: Sección viga (m)
- `--hinge_theta_y`, `--hinge_My`: Curva M-θ de rótulas

### 3. Ejecutar CalculiX

```bash
# Compilar CalculiX si no existe
cd /work/src
make -j16 CalculiX

# Ejecutar
cd /work/tests/portal4_compare/new
/work/src/CalculiX -i portal_A0.1g
```

### 4. Post-procesar hinges

```bash
cd /work/tests/portal4_compare/new
python3 ../post_hinge_theta.py portal_A0.1g.dat portal_A0.1g.hinge_map.json
```

Genera: `portal_A0.1g.H_L.mtheta.csv`, `portal_A0.1g.H_R.mtheta.csv`

### 5. Comparar baseline vs new

```bash
cd /work/tests/portal4_compare
python3 compare_mtheta.py
python3 plot_mtheta.py
```

## Benchmark Automático (Múltiples Amplitudes)

Ejecutar análisis paramétrico para varias amplitudes:

```bash
bash run_benchmark_problema4.sh 0.1 0.2 0.3 0.4
```

Esto:
1. Genera señales sísmicas para A=0.1g, 0.2g, 0.3g, 0.4g
2. Crea inputs de CalculiX para cada amplitud
3. Ejecuta simulaciones (si CCX disponible)
4. Post-procesa hinges
5. Compara con baseline (si existe)

## Estructura de Archivos

```
/work/tests/portal4_compare/
├── generate_earthquake_senoidal.py     # Generador de señal sísmica
├── generate_portal_eq_ccx.py          # Generador de .inp para Problema 4
├── post_hinge_theta.py                # Post-procesador de hinges
├── compare_mtheta.py                  # Comparación baseline vs new
├── plot_mtheta.py                     # Generador de plots M-θ
├── run_benchmark_problema4.sh         # Script maestro (múltiples A)
│
├── earthquake_A0.1g.csv               # Señales sísmicas generadas
├── earthquake_A0.2g.csv
├── ...
│
├── portal_A0.1g.inp                   # Inputs de CalculiX
├── portal_A0.1g.hinge_map.json
├── portal_A0.2g.inp
├── ...
│
├── baseline/                          # Outputs de referencia
│   ├── portal_A0.1g.dat
│   ├── portal_A0.1g.H_L.mtheta.csv
│   └── ...
│
├── new/                               # Outputs de new run
│   ├── portal_A0.1g.dat
│   ├── portal_A0.1g.H_L.mtheta.csv
│   └── ...
│
└── report/                            # Reportes de comparación
    ├── hinge_compare.txt
    ├── H_L.mtheta.png
    ├── H_R.mtheta.png
    └── benchmark.log
```

## Modelo de Elementos Finitos

### Geometría 3D con sólidos
- **Columnas**: Bloques de elementos C3D8R (40×60×20 cm)
- **Viga**: Bloque de elementos C3D8R (25×50×20 cm)
- **Malla**: 2×16×1 elems (columnas), 20×2×1 elems (viga)

### Rótulas plásticas (hinges)
Implementadas con:
- **RIGID BODY** en zonas de conexión
- **ROT nodes**: U3 = θ_z (rotación relativa)
- **SPRING2** entre pares de ROT nodes:
  - H_L: Conexión columna izq. - viga
  - H_R: Conexión columna der. - viga

### Curva M-θ de rótulas
Por defecto (simplificada, **actualizar con análisis de sección RC**):
```
θ = 0.000 rad  →  M = 0 N·m
θ = 0.002 rad  →  M = 200 kN·m  (yield)
θ = 0.020 rad  →  M = 236 kN·m  (post-yield)
```

**TODO**: Calcular curva real basándose en:
- Diagrama de interacción M-N para columnas (S1)
- Capacidad flexural de vigas (S2)
- Acero φ=20mm, σy_steel ≈ 420 MPa
- Hormigón σy ≈ 42 MPa

### Análisis dinámico
- **Step 1 (STATIC)**: Gravedad
- **Step 2 (DYNAMIC)**:
  - t = 0 → 10 s
  - Δt = 0.0025 s (fijo)
  - Damping: β = 2×10⁻⁴
  - Excitación: Amplitud EQX (dirección horizontal +X)

## Métricas de Comparación

Para cada hinge, se calculan:

| Métrica | Descripción |
|---------|-------------|
| **max\|Δt\|** | Máxima diferencia de tiempo |
| **max\|Δθ\|** | Máxima diferencia de rotación (rad) |
| **RMS(Δθ)** | RMS de diferencias de rotación |
| **max\|ΔM\|** | Máxima diferencia de momento (N·m) |
| **RMS(ΔM)** | RMS de diferencias de momento |

## Parámetros Avanzados

### Material (hormigón)
```bash
--E 25e9        # Módulo de Young (Pa), default: 25 GPa
--nu 0.20       # Razón de Poisson
--rho 2400      # Densidad (kg/m³)
```

### Damping
```bash
--damp_alpha 0.0      # Damping proporcional a masa
--damp_beta 2.0e-4    # Damping proporcional a rigidez (default)
```

### Hinge (curva M-θ)
```bash
--hinge_theta_y 0.002       # Rotación de fluencia (rad)
--hinge_My 200e3            # Momento de fluencia (N·m)
--hinge_post_ratio 0.02     # Ratio rigidez post-yield / elástica
```

### Malla
```bash
--nx_col 2      # Elementos a lo ancho de columna
--ny_col 16     # Elementos a lo alto de columna
--nx_beam 20    # Elementos a lo largo de viga
--ny_beam 2     # Elementos a lo alto de viga
--nz 1          # Elementos en dirección fuera del plano
```

## Análisis de Resultados

### Señal sísmica
```python
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 10, 4001)
A = 0.1
a_g = A * np.cos(0.2*np.pi*t) * np.sin(4*np.pi*t)

plt.plot(t, a_g)
plt.xlabel('Time (s)')
plt.ylabel('Acceleration (g)')
plt.title(f'Senoidal Earthquake: A={A}g')
plt.grid()
plt.show()
```

### Curvas M-θ
```python
import pandas as pd
import matplotlib.pyplot as plt

# Leer datos
base = pd.read_csv('baseline/portal_A0.1g.H_L.mtheta.csv')
new = pd.read_csv('new/portal_A0.1g.H_L.mtheta.csv')

# Plotear
plt.figure(figsize=(10, 6))
plt.plot(base['theta_rel_rad'], base['M_Nm']/1e3, 'b-', label='Baseline', linewidth=2)
plt.plot(new['theta_rel_rad'], new['M_Nm']/1e3, 'r--', label='New', linewidth=2)
plt.xlabel('Rotation θ (rad)')
plt.ylabel('Moment M (kN·m)')
plt.title('M-θ Hinge Response (H_L)')
plt.legend()
plt.grid(alpha=0.3)
plt.show()
```

## Próximos Pasos

### 1. Calcular curva M-θ real para secciones RC
- Diagrama de interacción M-N (columnas)
- Capacidad a flexión (vigas)
- Considerar endurecimiento por deformación
- Modelo con degradación (vigas)

### 2. Calibrar masa para T₀=0.5s
- Cálculo del período fundamental
- Ajustar densidad o masa concentrada

### 3. Análisis paramétrico completo
- Ejecutar para A = 0.1g, 0.2g, ..., hasta colapso
- Identificar amplitud de colapso
- Curvas IDA (Incremental Dynamic Analysis)

### 4. Implementar modelo SHM para vigas
- Degradación de rigidez
- Degradación de resistencia
- Efecto pinching

### 5. Validación
- Comparar con solución analítica (si existe)
- Comparar con otros códigos (SAP2000, OpenSees, etc.)

## Referencias

- **Problema 2**: Secciones de hormigón armado (S1, S2)
- **Problema 4**: Especificaciones del benchmark (Tarea3_DC_2015)
- CalculiX documentation: http://www.calculix.de
- Sivaselvan & Reinhorn (2000): Modelo SHM con degradación

## Contacto / Issues

Para problemas con este benchmark:
- Revisa los logs en `report/`
- Verifica que CalculiX esté compilado correctamente
- Asegúrate de que SPOOLES esté instalado

