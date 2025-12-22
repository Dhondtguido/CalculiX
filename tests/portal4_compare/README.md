# CalculiX Portal Benchmark - Comparison Suite

Automated benchmark comparison for CalculiX "Pórtico Problema 4" (portal with solid elements, ROT hinges, and SPRING2).

## Overview

This suite compares a **baseline run** (reference outputs from a zip file) against a **new run** (using the current CalculiX executable).

**Outputs:**
- `report/hinge_compare.txt` - Detailed metrics per hinge
- `report/H_*.mtheta.png` - M-θ plots (baseline vs new)
- `report/run.log` - Execution log

## Quick Start

### Option 1: Automated (recommended)

```bash
cd /work/tests/portal4_compare

# Set paths (if different from defaults)
export ZIP_FILE=/path/to/ccx_portal_solid_hinge_eq_v11_incfix.zip
export CCX_EXEC=/work/src/CalculiX

# Run everything
bash run_benchmark.sh
```

### Option 2: Manual step-by-step

```bash
cd /work/tests/portal4_compare

# 1. Extract baseline
unzip -o /path/to/baseline.zip -d baseline/

# 2. Generate baseline CSVs
cd baseline
python3 ../post_hinge_theta.py portal_eq.dat portal_eq.hinge_map.json
cd ..

# 3. Prepare new run
rsync -a baseline/ new/ --exclude="*.dat" --exclude="*.csv"

# 4. Run CalculiX
cd new
/work/src/CalculiX -i portal_eq
cd ..

# 5. Generate new CSVs
cd new
python3 ../post_hinge_theta.py portal_eq.dat
cd ..

# 6. Compare and plot
python3 compare_mtheta.py
python3 plot_mtheta.py
```

## Scripts

### `run_benchmark.sh`
Master orchestration script. Runs all steps automatically.

**Environment variables:**
- `ZIP_FILE` - Path to baseline zip (default: `/mnt/data/ccx_portal_solid_hinge_eq_v11_incfix.zip`)
- `CCX_EXEC` - CalculiX executable (default: `/work/src/CalculiX`)

### `post_hinge_theta.py`
Extracts hinge moment-rotation data from `.dat` files.

**Usage:**
```bash
python3 post_hinge_theta.py <file.dat> [hinge_map.json]
```

**Output:**
- `<basename>.H_<name>.mtheta.csv` for each hinge
- Columns: `t`, `theta_rel_rad`, `M_Nm`

### `compare_mtheta.py`
Compares baseline vs new CSV files and generates metrics.

**Output:**
- `report/hinge_compare.txt` with per-hinge metrics:
  - max|Δt|, max|Δθ|, RMS(Δθ), max|ΔM|, RMS(ΔM)
  - Summary of worst hinges

**Usage:**
```bash
python3 compare_mtheta.py
```

### `plot_mtheta.py`
Generates M-θ comparison plots.

**Usage:**
```bash
# Plot first 2 hinges (default)
python3 plot_mtheta.py

# Plot specific hinges
python3 plot_mtheta.py H_L H_R
```

**Output:**
- PNG files (if matplotlib available): `report/H_*.mtheta.png`
- CSV files (fallback): `report/H_*.mtheta_combined.csv`

## Directory Structure

```
/work/tests/portal4_compare/
├── baseline/              # Extracted baseline run
│   ├── portal_eq.inp
│   ├── portal_eq.dat      # Reference outputs
│   └── *.mtheta.csv       # Generated CSV files
├── new/                   # New simulation run
│   ├── portal_eq.inp
│   ├── portal_eq.dat      # New outputs
│   └── *.mtheta.csv
├── report/                # Comparison results
│   ├── hinge_compare.txt
│   ├── H_*.mtheta.png
│   └── run.log
├── run_benchmark.sh       # Master script
├── post_hinge_theta.py
├── compare_mtheta.py
└── plot_mtheta.py
```

## Prerequisites

### Required:
- Python 3.6+
- CalculiX executable (compiled from `/work/src/`)
- Baseline zip file with:
  - `.inp` input file(s)
  - `.dat` reference outputs
  - Optional: `hinge_map.json`, run scripts

### Optional:
- `matplotlib` for PNG plots (gracefully degrades to CSV export)
- `rsync` (uses fallback if not available)

## Metrics Explained

For each hinge, we compute:

| Metric | Description |
|--------|-------------|
| **max\|Δt\|** | Maximum time difference between datasets |
| **max\|Δθ\|** | Maximum absolute rotation difference (rad) |
| **RMS(Δθ)** | Root-mean-square rotation difference (rad) |
| **max\|ΔM\|** | Maximum absolute moment difference (N·m) |
| **RMS(ΔM)** | Root-mean-square moment difference (N·m) |

Differences are computed at baseline time points using linear interpolation on the new data.

## Troubleshooting

### Error: Baseline zip not found
```bash
export ZIP_FILE=/correct/path/to/baseline.zip
bash run_benchmark.sh
```

### Error: CalculiX executable not found
Compile CalculiX first:
```bash
cd /work/src
make -j$(nproc) CalculiX
```

Or set custom path:
```bash
export CCX_EXEC=/path/to/CalculiX
```

### Warning: No hinge data found in .dat file
Check that:
1. `.dat` file contains SPRING2/hinge element outputs
2. Element output is enabled in `.inp` file
3. File format matches expectations in `post_hinge_theta.py`

You may need to adapt the parser for your specific output format.

### No matplotlib available
Script will automatically generate CSV files instead of PNG plots.
You can install matplotlib (if you have sudo):
```bash
sudo apt-get install python3-matplotlib
```

## Customization

### Adding custom hinge names
Create `baseline/portal_eq.hinge_map.json`:
```json
{
  "SPRING2_123": "H_L",
  "SPRING2_456": "H_R"
}
```

### Changing DAT parser
Edit `post_hinge_theta.py` function `parse_dat_file()` to match your output format.

### Modifying metrics
Edit `compare_mtheta.py` function `compare_hinges()` to add custom metrics.

## Example Output

```
============================================
SUMMARY
============================================
Total hinges compared: 4

Worst by max|Δθ|:
  H_L                  : 1.234567e-05 rad

Worst by max|ΔM|:
  H_R                  : 2.345678e+02 N·m
============================================
```

## License

Part of CalculiX test suite. See main repository LICENSE.

## Contact

For issues with this benchmark suite, check:
- Main CalculiX repository
- CalculiX discourse: https://calculix.discourse.group
