#!/usr/bin/env python3
"""
Compare baseline vs new hinge M-theta curves.
Generates metrics and reports differences.

Usage:
    python3 compare_mtheta.py

Expects:
    - baseline/*.mtheta.csv files
    - new/*.mtheta.csv files (matching names)

Output:
    - report/hinge_compare.txt with detailed metrics
"""

import os
import csv
import glob
import sys
from pathlib import Path
import math


def load_mtheta_csv(csv_path):
    """
    Load M-theta CSV file.
    Returns: [(t, theta, M), ...]
    """
    data = []

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = float(row['t'])
            theta = float(row['theta_rel_rad'])
            M = float(row['M_Nm'])
            data.append((t, theta, M))

    return data


def interpolate(data, t_target):
    """
    Linear interpolation to find theta and M at t_target.
    Returns: (theta, M) or None if t_target out of range
    """
    if not data:
        return None

    # Sort by time
    data = sorted(data, key=lambda x: x[0])

    t_min = data[0][0]
    t_max = data[-1][0]

    if t_target < t_min or t_target > t_max:
        return None

    # Find bracketing points
    for i in range(len(data) - 1):
        t1, theta1, M1 = data[i]
        t2, theta2, M2 = data[i + 1]

        if t1 <= t_target <= t2:
            if t2 == t1:
                return (theta1, M1)

            # Linear interpolation
            frac = (t_target - t1) / (t2 - t1)
            theta = theta1 + frac * (theta2 - theta1)
            M = M1 + frac * (M2 - M1)
            return (theta, M)

    return None


def compare_hinges(baseline_data, new_data):
    """
    Compare two hinge datasets.
    Returns: dict with metrics
    """
    metrics = {
        'n_baseline': len(baseline_data),
        'n_new': len(new_data),
        'max_dt': 0.0,
        'max_dtheta': 0.0,
        'rms_dtheta': 0.0,
        'max_dM': 0.0,
        'rms_dM': 0.0,
        'comparable_points': 0
    }

    if not baseline_data or not new_data:
        return metrics

    # Get common time range
    t_base = [t for t, _, _ in baseline_data]
    t_new = [t for t, _, _ in new_data]

    t_min = max(min(t_base), min(t_new))
    t_max = min(max(t_base), max(t_new))

    if t_min >= t_max:
        return metrics

    # Sample at baseline time points within common range
    dtheta_sq_sum = 0.0
    dM_sq_sum = 0.0
    count = 0

    for t_base, theta_base, M_base in baseline_data:
        if t_base < t_min or t_base > t_max:
            continue

        # Interpolate new data at this time point
        new_interp = interpolate(new_data, t_base)
        if new_interp is None:
            continue

        theta_new, M_new = new_interp

        # Compute differences
        dtheta = abs(theta_new - theta_base)
        dM = abs(M_new - M_base)

        metrics['max_dtheta'] = max(metrics['max_dtheta'], dtheta)
        metrics['max_dM'] = max(metrics['max_dM'], dM)

        dtheta_sq_sum += dtheta ** 2
        dM_sq_sum += dM ** 2
        count += 1

    metrics['comparable_points'] = count

    if count > 0:
        metrics['rms_dtheta'] = math.sqrt(dtheta_sq_sum / count)
        metrics['rms_dM'] = math.sqrt(dM_sq_sum / count)

    # Max time difference (just check endpoints)
    metrics['max_dt'] = abs(max(t_base) - max(t_new))

    return metrics


def format_metric(value, precision=6):
    """Format metric value with scientific notation if needed."""
    if abs(value) < 1e-3 or abs(value) > 1e6:
        return f"{value:.{precision}e}"
    else:
        return f"{value:.{precision}f}"


def main():
    baseline_dir = Path("baseline")
    new_dir = Path("new")
    report_dir = Path("report")

    report_dir.mkdir(exist_ok=True)

    # Find all baseline CSV files
    baseline_csvs = list(baseline_dir.glob("*.mtheta.csv"))

    if not baseline_csvs:
        print("ERROR: No baseline *.mtheta.csv files found")
        sys.exit(1)

    print(f"Found {len(baseline_csvs)} baseline CSV file(s)")

    # Compare each hinge
    results = []

    for baseline_csv in sorted(baseline_csvs):
        hinge_name = baseline_csv.stem.replace('.mtheta', '')
        new_csv = new_dir / baseline_csv.name

        print(f"\nComparing: {hinge_name}")

        if not new_csv.exists():
            print(f"  WARNING: No matching new file: {new_csv.name}")
            continue

        # Load data
        baseline_data = load_mtheta_csv(baseline_csv)
        new_data = load_mtheta_csv(new_csv)

        print(f"  Baseline: {len(baseline_data)} points")
        print(f"  New:      {len(new_data)} points")

        # Compare
        metrics = compare_hinges(baseline_data, new_data)

        results.append({
            'hinge': hinge_name,
            'metrics': metrics,
            'baseline_csv': baseline_csv,
            'new_csv': new_csv
        })

        print(f"  Comparable points: {metrics['comparable_points']}")
        print(f"  max|Δθ|: {format_metric(metrics['max_dtheta'])} rad")
        print(f"  max|ΔM|: {format_metric(metrics['max_dM'])} N·m")

    # Write detailed report
    report_path = report_dir / "hinge_compare.txt"

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CALCULIX HINGE COMPARISON REPORT\n")
        f.write("Baseline vs New Run\n")
        f.write("=" * 80 + "\n\n")

        for result in results:
            hinge = result['hinge']
            m = result['metrics']

            f.write(f"\n{'=' * 80}\n")
            f.write(f"HINGE: {hinge}\n")
            f.write(f"{'=' * 80}\n")
            f.write(f"Baseline file:      {result['baseline_csv'].name}\n")
            f.write(f"New file:           {result['new_csv'].name}\n")
            f.write(f"Baseline points:    {m['n_baseline']}\n")
            f.write(f"New points:         {m['n_new']}\n")
            f.write(f"Comparable points:  {m['comparable_points']}\n")
            f.write(f"\nMETRICS:\n")
            f.write(f"  max|Δt|         = {format_metric(m['max_dt'])} s\n")
            f.write(f"  max|Δθ|         = {format_metric(m['max_dtheta'])} rad\n")
            f.write(f"  RMS(Δθ)         = {format_metric(m['rms_dtheta'])} rad\n")
            f.write(f"  max|ΔM|         = {format_metric(m['max_dM'])} N·m\n")
            f.write(f"  RMS(ΔM)         = {format_metric(m['rms_dM'])} N·m\n")

        # Summary: worst hinges
        f.write(f"\n\n{'=' * 80}\n")
        f.write("SUMMARY: TOP WORST HINGES\n")
        f.write(f"{'=' * 80}\n\n")

        # Sort by max|Δθ|
        sorted_by_theta = sorted(results, key=lambda r: r['metrics']['max_dtheta'], reverse=True)

        f.write("By max|Δθ|:\n")
        for i, result in enumerate(sorted_by_theta[:5], 1):
            m = result['metrics']
            f.write(f"  {i}. {result['hinge']:20s} : {format_metric(m['max_dtheta'])} rad\n")

        f.write("\n")

        # Sort by max|ΔM|
        sorted_by_M = sorted(results, key=lambda r: r['metrics']['max_dM'], reverse=True)

        f.write("By max|ΔM|:\n")
        for i, result in enumerate(sorted_by_M[:5], 1):
            m = result['metrics']
            f.write(f"  {i}. {result['hinge']:20s} : {format_metric(m['max_dM'])} N·m\n")

        f.write("\n" + "=" * 80 + "\n")

    print(f"\n\nReport written to: {report_path}")

    # Return worst hinge info for console summary
    if results:
        worst_theta = sorted_by_theta[0]
        worst_M = sorted_by_M[0]

        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        print(f"Total hinges compared: {len(results)}")
        print(f"\nWorst by max|Δθ|:")
        print(f"  {worst_theta['hinge']:20s} : {format_metric(worst_theta['metrics']['max_dtheta'])} rad")
        print(f"\nWorst by max|ΔM|:")
        print(f"  {worst_M['hinge']:20s} : {format_metric(worst_M['metrics']['max_dM'])} N·m")
        print("=" * 80)


if __name__ == '__main__':
    main()
