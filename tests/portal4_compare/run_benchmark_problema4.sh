#!/bin/bash
#
# Master script for Problema 4 benchmark comparison
# Generates inputs for multiple amplitudes and compares baseline vs new
#
# Usage:
#   bash run_benchmark_problema4.sh [amplitude1] [amplitude2] ...
#
# Example:
#   bash run_benchmark_problema4.sh 0.1 0.2 0.3
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LOG_FILE="report/benchmark.log"
mkdir -p report

echo "========================================" | tee "$LOG_FILE"
echo "Problema 4 - Portal Benchmark" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Default amplitudes if none specified
AMPLITUDES=("${@:-0.1 0.2}")

echo "Amplitudes to test: ${AMPLITUDES[@]}" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# ============================================
# Step 1: Generate earthquake signals
# ============================================
echo "[1/5] Generating earthquake signals..." | tee -a "$LOG_FILE"

for A in "${AMPLITUDES[@]}"; do
    echo "  A = ${A}g" | tee -a "$LOG_FILE"
    python3 generate_earthquake_senoidal.py \
        --A "$A" \
        --dt 0.0025 \
        --duration 10.0 \
        --units g \
        --out "earthquake_A${A}g.csv" 2>&1 | tee -a "$LOG_FILE"
done

echo "" | tee -a "$LOG_FILE"

# ============================================
# Step 2: Generate CalculiX inputs
# ============================================
echo "[2/5] Generating CalculiX inputs..." | tee -a "$LOG_FILE"

for A in "${AMPLITUDES[@]}"; do
    echo "  portal_A${A}g.inp" | tee -a "$LOG_FILE"
    python3 generate_portal_eq_ccx.py \
        --eq "earthquake_A${A}g.csv" \
        --out "portal_A${A}g.inp" \
        --H 3.0 \
        --L 5.0 \
        --col_b 0.40 \
        --col_h 0.60 \
        --beam_b 0.25 \
        --beam_h 0.50 2>&1 | tee -a "$LOG_FILE"
done

echo "" | tee -a "$LOG_FILE"

# ============================================
# Step 3: Check for baseline data
# ============================================
echo "[3/5] Checking for baseline data..." | tee -a "$LOG_FILE"

BASELINE_DIR="baseline"
if [ ! -d "$BASELINE_DIR" ] || [ -z "$(ls -A $BASELINE_DIR 2>/dev/null)" ]; then
    echo "WARNING: No baseline data found in $BASELINE_DIR/" | tee -a "$LOG_FILE"
    echo "Please provide baseline outputs or run with CCX to generate them." | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "To generate baseline, you can:" | tee -a "$LOG_FILE"
    echo "  1. Extract baseline zip to baseline/" | tee -a "$LOG_FILE"
    echo "  2. Or run: CCX=/path/to/CalculiX bash run_all.sh" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"

    # Ask if user wants to continue without baseline
    read -p "Continue without baseline comparison? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Exiting." | tee -a "$LOG_FILE"
        exit 1
    fi
    SKIP_COMPARISON=true
else
    echo "  ✓ Baseline data found" | tee -a "$LOG_FILE"
    SKIP_COMPARISON=false
fi

echo "" | tee -a "$LOG_FILE"

# ============================================
# Step 4: Run new simulations (if CCX available)
# ============================================
echo "[4/5] Running new simulations..." | tee -a "$LOG_FILE"

CCX_EXEC="${CCX:-/work/src/CalculiX}"

if [ ! -x "$CCX_EXEC" ]; then
    echo "WARNING: CalculiX executable not found: $CCX_EXEC" | tee -a "$LOG_FILE"
    echo "Set CCX environment variable or compile CalculiX first." | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    SKIP_RUN=true
else
    echo "  Using: $CCX_EXEC" | tee -a "$LOG_FILE"
    SKIP_RUN=false

    mkdir -p new

    for A in "${AMPLITUDES[@]}"; do
        echo "" | tee -a "$LOG_FILE"
        echo "  Running A=${A}g..." | tee -a "$LOG_FILE"

        # Copy input to new/
        cp "portal_A${A}g.inp" "new/"
        cp "portal_A${A}g.hinge_map.json" "new/"

        cd new

        BASE="portal_A${A}g"

        # Clean old outputs
        rm -f "${BASE}".{dat,frd,sta,cvg,12d,log,msg,prt} 2>/dev/null || true

        # Run CalculiX
        echo "    Running CCX..." | tee -a "../$LOG_FILE"
        set +e
        "$CCX_EXEC" -i "$BASE" 2>&1 | tee "../report/${BASE}.ccx.log"
        ccx_rc=${PIPESTATUS[0]}
        set -e

        if [[ $ccx_rc -ne 0 ]]; then
            echo "    WARNING: CCX exited with code $ccx_rc" | tee -a "../$LOG_FILE"
        else
            echo "    ✓ CCX completed" | tee -a "../$LOG_FILE"
        fi

        # Post-process hinges
        if [ -f "${BASE}.dat" ]; then
            echo "    Post-processing hinges..." | tee -a "../$LOG_FILE"
            python3 ../post_hinge_theta.py "${BASE}.dat" "${BASE}.hinge_map.json" 2>&1 | tee -a "../$LOG_FILE"
        else
            echo "    ERROR: No .dat file generated" | tee -a "../$LOG_FILE"
        fi

        cd ..
    done
fi

echo "" | tee -a "$LOG_FILE"

# ============================================
# Step 5: Compare baseline vs new
# ============================================
if [ "$SKIP_COMPARISON" = false ] && [ "$SKIP_RUN" = false ]; then
    echo "[5/5] Comparing baseline vs new..." | tee -a "$LOG_FILE"

    for A in "${AMPLITUDES[@]}"; do
        echo "" | tee -a "$LOG_FILE"
        echo "  A=${A}g:" | tee -a "$LOG_FILE"

        # Check if CSV files exist
        BASE_CSV=$(ls baseline/portal_A${A}g.*.mtheta.csv 2>/dev/null | head -1)
        NEW_CSV=$(ls new/portal_A${A}g.*.mtheta.csv 2>/dev/null | head -1)

        if [ -z "$BASE_CSV" ] || [ -z "$NEW_CSV" ]; then
            echo "    WARNING: Missing CSV files, skipping comparison" | tee -a "$LOG_FILE"
            continue
        fi

        # Run comparison (would need to adapt compare_mtheta.py for multiple runs)
        echo "    Comparing hinge outputs..." | tee -a "$LOG_FILE"
        # python3 compare_mtheta.py --baseline baseline/ --new new/ --amplitude "$A" 2>&1 | tee -a "$LOG_FILE"
    done

    echo "" | tee -a "$LOG_FILE"
    echo "  Note: Detailed comparison requires adapted compare_mtheta.py" | tee -a "$LOG_FILE"
else
    echo "[5/5] Skipping comparison (baseline or new run not available)" | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "BENCHMARK COMPLETE" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "End: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Generated files:" | tee -a "$LOG_FILE"
echo "  Earthquakes: earthquake_A*.csv" | tee -a "$LOG_FILE"
echo "  Inputs:      portal_A*.inp" | tee -a "$LOG_FILE"
if [ "$SKIP_RUN" = false ]; then
    echo "  Outputs:     new/portal_A*.dat" | tee -a "$LOG_FILE"
    echo "  CSV data:    new/portal_A*.*.mtheta.csv" | tee -a "$LOG_FILE"
fi
echo "  Logs:        report/*.log" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

echo "Done!" | tee -a "$LOG_FILE"
