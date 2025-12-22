#!/usr/bin/env python3
"""
Generate M-theta plots for baseline vs new comparison.
Requires matplotlib, falls back to CSV export if not available.

Usage:
    python3 plot_mtheta.py [hinge1] [hinge2] ...

If no hinges specified, plots first 2 found.
"""

import sys
import csv
import glob
from pathlib import Path

# Try to import matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available, will generate CSV files instead")


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

    return sorted(data)


def plot_comparison_matplotlib(hinge_name, baseline_data, new_data, output_path):
    """
    Generate M-theta plot using matplotlib.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Extract theta and M
    theta_base = [theta for _, theta, _ in baseline_data]
    M_base = [M for _, _, M in baseline_data]

    theta_new = [theta for _, theta, _ in new_data]
    M_new = [M for _, _, M in new_data]

    # Plot
    ax.plot(theta_base, M_base, 'b-', linewidth=2, label='Baseline', alpha=0.7)
    ax.plot(theta_new, M_new, 'r--', linewidth=2, label='New', alpha=0.7)

    ax.set_xlabel('Rotation θ (rad)', fontsize=12)
    ax.set_ylabel('Moment M (N·m)', fontsize=12)
    ax.set_title(f'M-θ Comparison: {hinge_name}', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"  Plot saved: {output_path}")


def export_comparison_csv(hinge_name, baseline_data, new_data, output_path):
    """
    Export combined CSV for manual plotting.
    """
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['dataset', 't', 'theta_rad', 'M_Nm'])

        for t, theta, M in baseline_data:
            writer.writerow(['baseline', f'{t:.6e}', f'{theta:.6e}', f'{M:.6e}'])

        for t, theta, M in new_data:
            writer.writerow(['new', f'{t:.6e}', f'{theta:.6e}', f'{M:.6e}'])

    print(f"  CSV exported: {output_path}")


def main():
    baseline_dir = Path("baseline")
    new_dir = Path("new")
    report_dir = Path("report")

    report_dir.mkdir(exist_ok=True)

    # Find baseline CSV files
    baseline_csvs = sorted(baseline_dir.glob("*.mtheta.csv"))

    if not baseline_csvs:
        print("ERROR: No baseline *.mtheta.csv files found")
        sys.exit(1)

    # Determine which hinges to plot
    if len(sys.argv) > 1:
        # User specified hinges
        hinge_names = sys.argv[1:]
        selected_csvs = []
        for hinge in hinge_names:
            found = False
            for csv_path in baseline_csvs:
                if hinge in csv_path.stem:
                    selected_csvs.append(csv_path)
                    found = True
                    break
            if not found:
                print(f"WARNING: Hinge '{hinge}' not found, skipping")
    else:
        # Default: first 2 hinges
        selected_csvs = baseline_csvs[:2]

    if not selected_csvs:
        print("ERROR: No hinges to plot")
        sys.exit(1)

    print(f"Generating {'plots' if HAS_MATPLOTLIB else 'CSV exports'} for {len(selected_csvs)} hinge(s):")

    for baseline_csv in selected_csvs:
        hinge_name = baseline_csv.stem.replace('.mtheta', '')
        new_csv = new_dir / baseline_csv.name

        print(f"\n{hinge_name}:")

        if not new_csv.exists():
            print(f"  WARNING: No matching new file, skipping")
            continue

        # Load data
        baseline_data = load_mtheta_csv(baseline_csv)
        new_data = load_mtheta_csv(new_csv)

        print(f"  Baseline: {len(baseline_data)} points")
        print(f"  New:      {len(new_data)} points")

        # Generate output
        if HAS_MATPLOTLIB:
            output_path = report_dir / f"{hinge_name}.mtheta.png"
            plot_comparison_matplotlib(hinge_name, baseline_data, new_data, output_path)
        else:
            output_path = report_dir / f"{hinge_name}.mtheta_combined.csv"
            export_comparison_csv(hinge_name, baseline_data, new_data, output_path)

    print("\nDone!")


if __name__ == '__main__':
    main()
