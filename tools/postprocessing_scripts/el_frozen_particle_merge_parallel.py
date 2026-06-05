#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


FILENAME_RE = re.compile(r"el_frozen_particles_rank\d+_step(\d{6})\.csv$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge rank-local el_frozen_particles_rankXXXX_stepXXXXXX.csv files "
            "into one global CSV per timestep."
        )
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        default=".",
        help="Directory containing rank-local frozen-field particle CSV files.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="el_frozen_merged",
        help="Directory to write merged el_frozen_particles_stepXXXXXX.csv files into.",
    )
    return parser.parse_args()


def discover_step_files(input_dir: Path) -> dict[str, list[Path]]:
    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in sorted(input_dir.glob("el_frozen_particles_rank*_step*.csv")):
        match = FILENAME_RE.fullmatch(path.name)
        if match is None:
            continue
        grouped[match.group(1)].append(path)
    return dict(grouped)


def merge_step(step: str, files: list[Path], output_dir: Path) -> None:
    header: list[str] | None = None
    rows_by_id: dict[int, list[str]] = {}

    for path in files:
        with path.open("r", newline="") as handle:
            reader = csv.reader(handle)
            file_header = next(reader, None)
            if file_header is None:
                continue

            if header is None:
                header = file_header
            elif file_header != header:
                raise ValueError(f"Header mismatch in {path}")

            for row in reader:
                if not row:
                    continue
                particle_id = int(row[0].strip())
                if particle_id in rows_by_id:
                    raise ValueError(
                        f"Duplicate particle id {particle_id} encountered while merging step {step}"
                    )
                rows_by_id[particle_id] = row

    if header is None:
        return

    output_path = output_dir / f"el_frozen_particles_step{step}.csv"
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for particle_id in sorted(rows_by_id):
            writer.writerow(rows_by_id[particle_id])

    print(
        f"{output_path} <- {len(files)} rank files, {len(rows_by_id)} particles"
    )


def main() -> int:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.is_dir():
        print(f"Input directory does not exist: {input_dir}", file=sys.stderr)
        return 1

    grouped = discover_step_files(input_dir)
    if not grouped:
        print(
            f"No rank-local frozen-field CSV files found in {input_dir}",
            file=sys.stderr,
        )
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    for step, files in sorted(grouped.items()):
        merge_step(step, files, output_dir)

    print(f"Merged files written to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
