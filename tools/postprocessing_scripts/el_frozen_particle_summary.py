#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize el_frozen_particles_step*.csv diagnostics."
    )
    parser.add_argument(
        "csv_files",
        nargs="+",
        help="One or more frozen-field particle CSV files.",
    )
    return parser.parse_args()


def summarize_csv(path: Path) -> None:
    rows = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(row)

    if not rows:
        print(f"{path}: no particle rows")
        return

    re_values = [float(row["re_p"]) for row in rows]
    found_values = [float(row["found_count"]) for row in rows]
    force_mag = [
        math.sqrt(
            float(row["fx"]) ** 2 + float(row["fy"]) ** 2 + float(row["fz"]) ** 2
        )
        for row in rows
    ]

    print(path)
    print(f"  particle count      = {len(rows)}")
    print(
        "  Re_p avg/min/max   = "
        f"{sum(re_values) / len(re_values):.10e} "
        f"{min(re_values):.10e} "
        f"{max(re_values):.10e}"
    )
    print(
        "  |F| avg/min/max    = "
        f"{sum(force_mag) / len(force_mag):.10e} "
        f"{min(force_mag):.10e} "
        f"{max(force_mag):.10e}"
    )
    print(
        "  found_count avg    = "
        f"{sum(found_values) / len(found_values):.10e}"
    )


def main() -> None:
    args = parse_args()
    for name in args.csv_files:
        summarize_csv(Path(name))


if __name__ == "__main__":
    main()
