#!/bin/bash
#
# Master script to run CalculiX portal benchmark comparison
# Usage: bash run_benchmark.sh
#

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LOG_FILE="report/run.log"
mkdir -p report

echo "CalculiX Portal Benchmark - Comparison Script" | tee "$LOG_FILE"
echo "=============================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Configuration
ZIP_FILE="${ZIP_FILE:-/mnt/data/ccx_portal_solid_hinge_eq_v11_incfix.zip}"
CCX_EXEC="${CCX_EXEC:-/work/src/CalculiX}"

echo "Configuration:" | tee -a "$LOG_FILE"
echo "  Baseline zip: $ZIP_FILE" | tee -a "$LOG_FILE"
echo "  CalculiX exe: $CCX_EXEC" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# ============================================
# Step 1: Extract baseline
# ============================================
echo "[1/7] Extracting baseline..." | tee -a "$LOG_FILE"

if [ ! -f "$ZIP_FILE" ]; then
    echo "ERROR: Baseline zip not found: $ZIP_FILE" | tee -a "$LOG_FILE"
    echo "Please set ZIP_FILE environment variable or place zip at:" | tee -a "$LOG_FILE"
    echo "  $ZIP_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

rm -rf baseline/*
unzip -q -o "$ZIP_FILE" -d baseline/
echo "  ✓ Baseline extracted" | tee -a "$LOG_FILE"

# ============================================
# Step 2: Generate baseline CSVs
# ============================================
echo "[2/7] Generating baseline hinge CSVs..." | tee -a "$LOG_FILE"

cd baseline

# Auto-detect .dat file
DAT_FILE=$(ls *.dat 2>/dev/null | head -1)
if [ -z "$DAT_FILE" ]; then
    echo "ERROR: No .dat file found in baseline" | tee -a "../$LOG_FILE"
    exit 1
fi

echo "  Using: $DAT_FILE" | tee -a "../$LOG_FILE"

# Auto-detect hinge map (optional)
HINGE_MAP=""
if [ -f "*.hinge_map.json" ]; then
    HINGE_MAP=$(ls *.hinge_map.json | head -1)
    echo "  Hinge map: $HINGE_MAP" | tee -a "../$LOG_FILE"
fi

# Run post-processor
python3 ../post_hinge_theta.py "$DAT_FILE" $HINGE_MAP 2>&1 | tee -a "../$LOG_FILE"

# Count generated CSVs
CSV_COUNT=$(ls *.mtheta.csv 2>/dev/null | wc -l)
echo "  ✓ Generated $CSV_COUNT baseline CSV file(s)" | tee -a "../$LOG_FILE"

cd ..

# ============================================
# Step 3: Verify CalculiX executable
# ============================================
echo "[3/7] Verifying CalculiX executable..." | tee -a "$LOG_FILE"

if [ ! -x "$CCX_EXEC" ]; then
    echo "ERROR: CalculiX executable not found or not executable: $CCX_EXEC" | tee -a "$LOG_FILE"
    echo "Please compile CalculiX or set CCX_EXEC environment variable" | tee -a "$LOG_FILE"
    exit 1
fi

echo "  ✓ Found: $CCX_EXEC" | tee -a "$LOG_FILE"
"$CCX_EXEC" -v 2>&1 | head -3 | tee -a "$LOG_FILE" || true

# ============================================
# Step 4: Prepare new run directory
# ============================================
echo "[4/7] Preparing new run directory..." | tee -a "$LOG_FILE"

rm -rf new/*

# Copy only inputs (exclude outputs)
rsync -a baseline/ new/ \
    --exclude="*.dat" \
    --exclude="*.sta" \
    --exclude="*.frd" \
    --exclude="*.cvg" \
    --exclude="*.12d" \
    --exclude="spooles.out" \
    --exclude="*.mtheta.csv" \
    --exclude="*.log"

echo "  ✓ Inputs copied to new/" | tee -a "$LOG_FILE"

# ============================================
# Step 5: Run new simulation
# ============================================
echo "[5/7] Running new simulation..." | tee -a "$LOG_FILE"

cd new

# Detect input file (without .inp extension)
INP_FILE=$(ls *.inp 2>/dev/null | head -1)
if [ -z "$INP_FILE" ]; then
    echo "ERROR: No .inp file found in new/" | tee -a "../$LOG_FILE"
    exit 1
fi

INPUT_NAME="${INP_FILE%.inp}"
echo "  Input: $INPUT_NAME" | tee -a "../$LOG_FILE"

# Check for run script
if [ -f "run_all.sh" ]; then
    echo "  Using run_all.sh" | tee -a "../$LOG_FILE"
    CCX="$CCX_EXEC" bash run_all.sh 2>&1 | tee -a "../$LOG_FILE"
elif [ -f "RUNME" ] || [ -f "RUNME.sh" ]; then
    RUN_SCRIPT=$(ls RUNME* | head -1)
    echo "  Using $RUN_SCRIPT" | tee -a "../$LOG_FILE"
    CCX="$CCX_EXEC" bash "$RUN_SCRIPT" 2>&1 | tee -a "../$LOG_FILE"
else
    echo "  Running directly: $CCX_EXEC -i $INPUT_NAME" | tee -a "../$LOG_FILE"
    "$CCX_EXEC" -i "$INPUT_NAME" 2>&1 | tee -a "../$LOG_FILE"
fi

# Verify output
NEW_DAT=$(ls *.dat 2>/dev/null | head -1)
if [ ! -f "$NEW_DAT" ]; then
    echo "ERROR: Simulation failed, no .dat file generated" | tee -a "../$LOG_FILE"
    exit 1
fi

echo "  ✓ Simulation completed: $NEW_DAT" | tee -a "../$LOG_FILE"

cd ..

# ============================================
# Step 6: Generate new run CSVs
# ============================================
echo "[6/7] Generating new run hinge CSVs..." | tee -a "$LOG_FILE"

cd new

python3 ../post_hinge_theta.py "$NEW_DAT" $HINGE_MAP 2>&1 | tee -a "../$LOG_FILE"

CSV_COUNT_NEW=$(ls *.mtheta.csv 2>/dev/null | wc -l)
echo "  ✓ Generated $CSV_COUNT_NEW new CSV file(s)" | tee -a "../$LOG_FILE"

cd ..

# ============================================
# Step 7: Compare and generate reports
# ============================================
echo "[7/7] Generating comparison reports..." | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "Running comparison..." | tee -a "$LOG_FILE"
python3 compare_mtheta.py 2>&1 | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "Generating plots..." | tee -a "$LOG_FILE"
python3 plot_mtheta.py 2>&1 | tee -a "$LOG_FILE"

# ============================================
# Final summary
# ============================================
echo "" | tee -a "$LOG_FILE"
echo "=============================================" | tee -a "$LOG_FILE"
echo "BENCHMARK COMPLETE" | tee -a "$LOG_FILE"
echo "=============================================" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Output files:" | tee -a "$LOG_FILE"
echo "  - report/hinge_compare.txt  (detailed metrics)" | tee -a "$LOG_FILE"
echo "  - report/*.png              (M-theta plots)" | tee -a "$LOG_FILE"
echo "  - report/run.log            (execution log)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Display summary from comparison report if exists
if [ -f "report/hinge_compare.txt" ]; then
    echo "Quick summary from comparison:" | tee -a "$LOG_FILE"
    tail -20 report/hinge_compare.txt | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "Done!" | tee -a "$LOG_FILE"
