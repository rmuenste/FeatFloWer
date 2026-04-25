#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


VECTOR_FIELDS = {
    "position": ("x", "y", "z"),
    "particle_velocity": ("pvx", "pvy", "pvz"),
    "carrier_velocity": ("ux", "uy", "uz"),
    "slip_velocity": ("slipx", "slipy", "slipz"),
    "force": ("fx", "fy", "fz"),
}

SCALAR_FIELDS = ("re_p", "found_count", "owner_rank_sum", "elem_sum", "radius")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare a serial frozen-field particle CSV with a merged parallel "
            "frozen-field particle CSV."
        )
    )
    parser.add_argument("serial_csv", help="Serial global particle CSV file.")
    parser.add_argument("parallel_csv", help="Merged parallel global particle CSV file.")
    parser.add_argument(
        "--top",
        type=int,
        default=12,
        help="Number of largest-difference particle ids to print per metric.",
    )
    return parser.parse_args()


def vec_norm3(ax: float, ay: float, az: float) -> float:
    return math.sqrt(ax * ax + ay * ay + az * az)


def read_particles(path: Path) -> dict[int, dict[str, float]]:
    particles: dict[int, dict[str, float]] = {}
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            particle_id = int(row["id"])
            particles[particle_id] = {key: float(value) for key, value in row.items() if key}
    return particles


def avg_min_max(values: list[float]) -> tuple[float, float, float]:
    if not values:
        return 0.0, 0.0, 0.0
    return sum(values) / len(values), min(values), max(values)


def top_entries(entries: list[tuple[float, int]], count: int) -> list[tuple[float, int]]:
    return sorted(entries, reverse=True)[: max(0, count)]


def main() -> None:
    args = parse_args()
    serial_path = Path(args.serial_csv)
    parallel_path = Path(args.parallel_csv)

    serial_particles = read_particles(serial_path)
    parallel_particles = read_particles(parallel_path)

    serial_ids = set(serial_particles)
    parallel_ids = set(parallel_particles)
    common_ids = sorted(serial_ids & parallel_ids)
    only_serial = sorted(serial_ids - parallel_ids)
    only_parallel = sorted(parallel_ids - serial_ids)

    print(f"serial   = {serial_path}")
    print(f"parallel = {parallel_path}")
    print(f"  particle count serial/parallel = {len(serial_particles)} / {len(parallel_particles)}")
    print(f"  common particle ids            = {len(common_ids)}")
    print(f"  only serial particle ids       = {len(only_serial)}")
    print(f"  only parallel particle ids     = {len(only_parallel)}")

    if not common_ids:
        return

    vector_stats: dict[str, list[float]] = {name: [] for name in VECTOR_FIELDS}
    scalar_stats: dict[str, list[float]] = {name: [] for name in SCALAR_FIELDS}
    vector_top: dict[str, list[tuple[float, int]]] = {name: [] for name in VECTOR_FIELDS}
    scalar_top: dict[str, list[tuple[float, int]]] = {name: [] for name in SCALAR_FIELDS}

    for particle_id in common_ids:
        serial_row = serial_particles[particle_id]
        parallel_row = parallel_particles[particle_id]

        for metric, components in VECTOR_FIELDS.items():
            delta = vec_norm3(
                parallel_row[components[0]] - serial_row[components[0]],
                parallel_row[components[1]] - serial_row[components[1]],
                parallel_row[components[2]] - serial_row[components[2]],
            )
            vector_stats[metric].append(delta)
            vector_top[metric].append((delta, particle_id))

        for metric in SCALAR_FIELDS:
            delta = abs(parallel_row[metric] - serial_row[metric])
            scalar_stats[metric].append(delta)
            scalar_top[metric].append((delta, particle_id))

    for metric in ("position", "particle_velocity", "carrier_velocity", "slip_velocity", "force"):
        avg_v, min_v, max_v = avg_min_max(vector_stats[metric])
        print(
            f"  delta {metric:17s} avg/min/max = "
            f"{avg_v:.10e} {min_v:.10e} {max_v:.10e}"
        )

    for metric in ("re_p", "found_count", "owner_rank_sum", "elem_sum", "radius"):
        avg_v, min_v, max_v = avg_min_max(scalar_stats[metric])
        print(
            f"  delta {metric:17s} avg/min/max = "
            f"{avg_v:.10e} {min_v:.10e} {max_v:.10e}"
        )

    for metric in ("carrier_velocity", "particle_velocity", "force", "position", "found_count", "elem_sum"):
        if metric in vector_top:
            entries = top_entries(vector_top[metric], args.top)
        else:
            entries = top_entries(scalar_top[metric], args.top)

        print(f"  top {len(entries)} ids by delta {metric}:")
        for delta, particle_id in entries:
            print(f"    id={particle_id} delta={delta:.10e}")


if __name__ == "__main__":
    main()
