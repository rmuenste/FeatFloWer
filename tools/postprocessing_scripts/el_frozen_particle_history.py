#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract history for selected particle ids across el_frozen_particles_step*.csv files."
    )
    parser.add_argument(
        "csv_glob",
        help="Quoted glob for frozen-field CSV files, e.g. 'build.../el_frozen_particles_step*.csv'",
    )
    parser.add_argument(
        "particle_ids",
        nargs="+",
        type=int,
        help="Particle ids to extract.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional CSV file to write the combined history table.",
    )
    return parser.parse_args()


def step_key(path: Path) -> int:
    stem = path.stem
    suffix = stem.split("step")[-1]
    return int(suffix)


def read_rows(path: Path, selected_ids: set[int]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            particle_id = int(row["id"])
            if particle_id in selected_ids:
                row_copy = dict(row)
                row_copy["step"] = str(step_key(path))
                row_copy["source_file"] = str(path)
                rows.append(row_copy)
    return rows


def print_summary(rows: list[dict[str, str]], particle_id: int) -> None:
    particle_rows = [row for row in rows if int(row["id"]) == particle_id]
    if not particle_rows:
        print(f"id {particle_id}: no rows found")
        return

    first = particle_rows[0]
    last = particle_rows[-1]
    print(
        f"id {particle_id}: steps {first['step']} -> {last['step']} | "
        f"rows = {len(particle_rows)}"
    )
    print(
        "  first: "
        f"x={float(first['x']):.6e} y={float(first['y']):.6e} z={float(first['z']):.6e} | "
        f"pv=({float(first['pvx']):.6e}, {float(first['pvy']):.6e}, {float(first['pvz']):.6e}) | "
        f"u=({float(first['ux']):.6e}, {float(first['uy']):.6e}, {float(first['uz']):.6e}) | "
        f"Re_p={float(first['re_p']):.6e}"
    )
    print(
        "  last : "
        f"x={float(last['x']):.6e} y={float(last['y']):.6e} z={float(last['z']):.6e} | "
        f"pv=({float(last['pvx']):.6e}, {float(last['pvy']):.6e}, {float(last['pvz']):.6e}) | "
        f"u=({float(last['ux']):.6e}, {float(last['uy']):.6e}, {float(last['uz']):.6e}) | "
        f"Re_p={float(last['re_p']):.6e}"
    )


def write_output(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "step",
        "id",
        "x",
        "y",
        "z",
        "radius",
        "pvx",
        "pvy",
        "pvz",
        "ux",
        "uy",
        "uz",
        "slipx",
        "slipy",
        "slipz",
        "fx",
        "fy",
        "fz",
        "re_p",
        "found_count",
        "owner_rank_sum",
        "elem_sum",
        "source_file",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def main() -> None:
    args = parse_args()
    selected_ids = set(args.particle_ids)
    csv_paths = sorted(Path().glob(args.csv_glob), key=step_key)

    if not csv_paths:
        raise SystemExit(f"No files matched glob: {args.csv_glob}")

    rows: list[dict[str, str]] = []
    for path in csv_paths:
        rows.extend(read_rows(path, selected_ids))

    rows.sort(key=lambda row: (int(row["id"]), int(row["step"])))

    print(f"matched files = {len(csv_paths)}")
    for particle_id in sorted(selected_ids):
        print_summary(rows, particle_id)

    if args.output is not None:
        write_output(args.output, rows)
        print(f"wrote history table to {args.output}")


if __name__ == "__main__":
    main()
