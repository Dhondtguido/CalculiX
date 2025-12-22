#!/usr/bin/env python3
"""
Generate senoidal earthquake signal for Problema 4 (Tarea3_DC_2015).

Signal: a_g(t) = A * cos(0.2π*t) * sin(4π*t)
Duration: t ∈ [0, 10] seconds
Amplitude: A in units of g (0.1, 0.2, ..., until collapse)

Usage:
    python3 generate_earthquake_senoidal.py --A 0.1 --dt 0.0025 --out earthquake_A0.1g.csv
"""

import argparse
import csv
import math
from pathlib import Path

G = 9.81  # m/s²

def senoidal_signal(t, A_g):
    """
    Generate senoidal earthquake signal.

    Parameters:
        t: time value (seconds)
        A_g: amplitude in units of g

    Returns:
        a_g: acceleration in units of g
    """
    return A_g * math.cos(0.2 * math.pi * t) * math.sin(4 * math.pi * t)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--A", type=float, required=True,
                    help="Amplitude in units of g (e.g., 0.1, 0.2, 0.3, ...)")
    ap.add_argument("--dt", type=float, default=0.0025,
                    help="Time step in seconds (default: 0.0025)")
    ap.add_argument("--duration", type=float, default=10.0,
                    help="Duration in seconds (default: 10.0)")
    ap.add_argument("--out", type=str, default=None,
                    help="Output CSV filename (default: earthquake_A{A}g.csv)")
    ap.add_argument("--units", type=str, choices=['g', 'mps2'], default='g',
                    help="Output units: 'g' or 'mps2' (default: g)")

    args = ap.parse_args()

    # Generate output filename if not specified
    if args.out is None:
        args.out = f"earthquake_A{args.A}g.csv"

    out_path = Path(args.out)

    # Generate time array
    n_points = int(args.duration / args.dt) + 1
    t_values = [i * args.dt for i in range(n_points)]

    # Generate acceleration signal (in g)
    a_g_values = [senoidal_signal(t, args.A) for t in t_values]

    # Convert to desired units
    if args.units == 'mps2':
        a_out_values = [a * G for a in a_g_values]
        unit_label = "m/s^2"
    else:
        a_out_values = a_g_values
        unit_label = "g"

    # Write CSV
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['t', 'acc'])
        for ti, ai in zip(t_values, a_out_values):
            writer.writerow([f'{ti:.6f}', f'{ai:.8e}'])

    max_abs_a = max(abs(a) for a in a_out_values)

    print(f"Generated: {out_path}")
    print(f"  Points:   {len(t_values)}")
    print(f"  Duration: {t_values[-1]:.3f} s")
    print(f"  dt:       {args.dt:.6f} s")
    print(f"  A:        {args.A:.3f} g")
    print(f"  Units:    {unit_label}")
    print(f"  max|a|:   {max_abs_a:.6f} {unit_label}")

    # Plot if matplotlib and numpy available
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Convert to numpy arrays
        t_np = np.array(t_values)
        a_g_np = np.array(a_g_values)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

        # Time history
        ax1.plot(t_np, a_g_np, 'b-', linewidth=1)
        ax1.set_xlabel('Time (s)', fontsize=11)
        ax1.set_ylabel('Acceleration (g)', fontsize=11)
        ax1.set_title(f'Senoidal Earthquake Signal: A = {args.A}g', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.axhline(0, color='k', linewidth=0.5)

        # Spectrum (FFT)
        n = len(t_np)
        freq = np.fft.fftfreq(n, args.dt)
        fft_vals = np.fft.fft(a_g_np)

        # Only positive frequencies
        pos_mask = freq > 0
        freq_pos = freq[pos_mask]
        fft_mag = np.abs(fft_vals[pos_mask]) / n

        ax2.plot(freq_pos, fft_mag, 'r-', linewidth=1)
        ax2.set_xlabel('Frequency (Hz)', fontsize=11)
        ax2.set_ylabel('FFT Magnitude', fontsize=11)
        ax2.set_title('Frequency Content', fontsize=12)
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim([0, 5])  # Show up to 5 Hz

        plt.tight_layout()
        plot_path = out_path.with_suffix('.png')
        plt.savefig(plot_path, dpi=150)
        plt.close()

        print(f"  Plot:     {plot_path}")
    except (ImportError, ModuleNotFoundError):
        pass


if __name__ == '__main__':
    main()
