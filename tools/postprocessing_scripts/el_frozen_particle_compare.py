#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare two el_frozen_particles_step*.csv diagnostics files."
    )
    parser.add_argument("csv_a", help="Earlier frozen-field particle CSV file.")
    parser.add_argument("csv_b", help="Later frozen-field particle CSV file.")
    parser.add_argument(
        "--show-removed",
        type=int,
        default=16,
        help="Maximum number of removed particle ids to print (default: 16).",
    )
    return parser.parse_args()


def read_particles(path: Path) -> dict[int, dict[str, float]]:
    particles: dict[int, dict[str, float]] = {}
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            particle_id = int(row["id"])
            particles[particle_id] = {
                "x": float(row["x"]),
                "y": float(row["y"]),
                "z": float(row["z"]),
                "pvx": float(row["pvx"]),
                "pvy": float(row["pvy"]),
                "pvz": float(row["pvz"]),
                "found_count": float(row["found_count"]),
            }
    return particles


def vec_norm3(ax: float, ay: float, az: float) -> float:
    return math.sqrt(ax * ax + ay * ay + az * az)


def avg_min_max(values: list[float]) -> tuple[float, float, float]:
    if not values:
        return 0.0, 0.0, 0.0
    return sum(values) / len(values), min(values), max(values)


def main() -> None:
    args = parse_args()
    path_a = Path(args.csv_a)
    path_b = Path(args.csv_b)

    particles_a = read_particles(path_a)
    particles_b = read_particles(path_b)

    ids_a = set(particles_a)
    ids_b = set(particles_b)
    common_ids = sorted(ids_a & ids_b)
    removed_ids = sorted(ids_a - ids_b)
    added_ids = sorted(ids_b - ids_a)

    displacements: list[float] = []
    velocity_changes: list[float] = []
    found_count_a: list[float] = []
    found_count_b: list[float] = []

    for particle_id in common_ids:
        row_a = particles_a[particle_id]
        row_b = particles_b[particle_id]

        displacements.append(
            vec_norm3(
                row_b["x"] - row_a["x"],
                row_b["y"] - row_a["y"],
                row_b["z"] - row_a["z"],
            )
        )
        velocity_changes.append(
            vec_norm3(
                row_b["pvx"] - row_a["pvx"],
                row_b["pvy"] - row_a["pvy"],
                row_b["pvz"] - row_a["pvz"],
            )
        )
        found_count_a.append(row_a["found_count"])
        found_count_b.append(row_b["found_count"])

    disp_avg, disp_min, disp_max = avg_min_max(displacements)
    dvel_avg, dvel_min, dvel_max = avg_min_max(velocity_changes)
    found_a_avg, _, _ = avg_min_max(found_count_a)
    found_b_avg, _, _ = avg_min_max(found_count_b)

    print(f"{path_a} -> {path_b}")
    print(f"  particle count a/b        = {len(particles_a)} / {len(particles_b)}")
    print(f"  common particle ids       = {len(common_ids)}")
    print(f"  removed particle ids      = {len(removed_ids)}")
    print(f"  added particle ids        = {len(added_ids)}")
    print(
        "  displacement avg/min/max = "
        f"{disp_avg:.10e} {disp_min:.10e} {disp_max:.10e}"
    )
    print(
        "  |delta v| avg/min/max    = "
        f"{dvel_avg:.10e} {dvel_min:.10e} {dvel_max:.10e}"
    )
    print(
        "  found_count avg a/b      = "
        f"{found_a_avg:.10e} / {found_b_avg:.10e}"
    )

    if removed_ids:
        shown_removed = removed_ids[: max(0, args.show_removed)]
        print(
            "  removed ids shown        = "
            + ", ".join(str(particle_id) for particle_id in shown_removed)
        )
        if len(removed_ids) > len(shown_removed):
            print(f"  removed ids omitted      = {len(removed_ids) - len(shown_removed)}")

    if added_ids:
        shown_added = added_ids[: max(0, args.show_removed)]
        print(
            "  added ids shown          = "
            + ", ".join(str(particle_id) for particle_id in shown_added)
        )
        if len(added_ids) > len(shown_added):
            print(f"  added ids omitted        = {len(added_ids) - len(shown_added)}")


if __name__ == "__main__":
    main()
