#!/usr/bin/env python3
"""
Post-processor for CalculiX hinge outputs.
Extracts moment-rotation curves from .dat files and generates CSV files.

Usage:
    python3 post_hinge_theta.py <dat_file> <hinge_map.json>

Output:
    Generates <basename>.H_<hinge_name>.mtheta.csv for each hinge
    with columns: t, theta_rel_rad, M_Nm
"""

import sys
import os
import json
import re
import csv
from pathlib import Path


def parse_dat_file(dat_path):
    """
    Parse CalculiX .dat file to extract hinge data.
    Returns dict: {hinge_name: [(t, theta, M), ...]}
    """
    hinges = {}

    if not os.path.exists(dat_path):
        print(f"WARNING: {dat_path} not found")
        return hinges

    with open(dat_path, 'r') as f:
        content = f.read()

    # Look for hinge sections - these may vary by version
    # Common patterns:
    # - "HINGE H_L" or "Element HINGE_L" sections
    # - Lines with "rotation" and "moment" data

    # Try to find element-based output first
    # Pattern: look for sections like:
    # element ... hinge ...
    # time = ...
    # rotation = ...
    # moment = ...

    current_hinge = None
    current_time = None
    current_theta = None
    current_moment = None

    for line in content.split('\n'):
        line = line.strip()

        # Match element/hinge identifier
        # Examples: "element 123, type SPRING2", "HINGE H_L", etc.
        if 'hinge' in line.lower() or 'spring2' in line.lower():
            # Try to extract hinge name
            match = re.search(r'H_\w+|HINGE[_\s]+(\w+)', line, re.IGNORECASE)
            if match:
                hinge_name = match.group(0).replace(' ', '_')
                if hinge_name not in hinges:
                    hinges[hinge_name] = []
                current_hinge = hinge_name

        # Match time
        if 'time' in line.lower() and '=' in line:
            match = re.search(r'time\s*=\s*([\d.eE+-]+)', line, re.IGNORECASE)
            if match:
                current_time = float(match.group(1))

        # Match rotation (in radians)
        if 'rotation' in line.lower() and '=' in line:
            match = re.search(r'rotation\s*=\s*([\d.eE+-]+)', line, re.IGNORECASE)
            if match:
                current_theta = float(match.group(1))

        # Match moment (in N·m)
        if 'moment' in line.lower() and '=' in line:
            match = re.search(r'moment\s*=\s*([\d.eE+-]+)', line, re.IGNORECASE)
            if match:
                current_moment = float(match.group(1))

        # When we have a complete data point, save it
        if (current_hinge and current_time is not None and
            current_theta is not None and current_moment is not None):
            hinges[current_hinge].append((current_time, current_theta, current_moment))
            # Reset for next data point
            current_theta = None
            current_moment = None

    return hinges


def load_hinge_map(map_path):
    """
    Load hinge mapping from JSON file (if exists).
    Format: {"H_L": "Hinge_Left", ...}
    Returns dict or None if file doesn't exist.
    """
    if not os.path.exists(map_path):
        print(f"INFO: No hinge map found at {map_path}, using auto-detected names")
        return None

    with open(map_path, 'r') as f:
        return json.load(f)


def write_hinge_csv(basename, hinge_name, data):
    """
    Write hinge data to CSV file.
    Filename: <basename>.H_<hinge_name>.mtheta.csv
    Columns: t, theta_rel_rad, M_Nm
    """
    csv_path = f"{basename}.{hinge_name}.mtheta.csv"

    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['t', 'theta_rel_rad', 'M_Nm'])

        for t, theta, moment in sorted(data):
            writer.writerow([f'{t:.6e}', f'{theta:.6e}', f'{moment:.6e}'])

    print(f"Written: {csv_path} ({len(data)} points)")
    return csv_path


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    dat_file = sys.argv[1]
    hinge_map_file = sys.argv[2] if len(sys.argv) > 2 else None

    # Get basename (without extension)
    basename = os.path.splitext(dat_file)[0]

    print(f"Processing: {dat_file}")

    # Parse DAT file
    hinges = parse_dat_file(dat_file)

    if not hinges:
        print("WARNING: No hinge data found in .dat file")
        print("This could mean:")
        print("  1. The file format is different than expected")
        print("  2. No SPRING2/hinge elements in the model")
        print("  3. Need to check element output format")
        return

    # Load hinge mapping (optional)
    hinge_map = load_hinge_map(hinge_map_file) if hinge_map_file else None

    # Write CSV files
    print(f"\nFound {len(hinges)} hinge(s):")
    for hinge_name, data in hinges.items():
        mapped_name = hinge_map.get(hinge_name, hinge_name) if hinge_map else hinge_name
        write_hinge_csv(basename, mapped_name, data)

    print("\nDone!")


if __name__ == '__main__':
    main()
