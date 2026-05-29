#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections
import math
import pathlib
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read particle_force.log from a running q2p1_dns_drag case and report "
            "the last fully available timesteps."
        )
    )
    parser.add_argument(
        "logfile",
        nargs="?",
        default="particle_force.log",
        help="Path to particle_force.log",
    )
    parser.add_argument(
        "--last",
        type=int,
        default=3,
        help="Number of complete timesteps to report (default: 3)",
    )
    return parser.parse_args()


def load_records(logfile: pathlib.Path):
    records = collections.OrderedDict()

    with logfile.open("r", encoding="ascii", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 14:
                # Likely a partially written line while the simulation is running.
                continue

            try:
                time_key = parts[0]
                ip = int(parts[1])
                fx = float(parts[2])
                fy = float(parts[3])
                fz = float(parts[4])
            except ValueError:
                continue

            if time_key not in records:
                records[time_key] = []
            records[time_key].append((ip, fx, fy, fz, line_no))

    return records


def expected_particle_count(records) -> int:
    counts = [len(entries) for entries in records.values() if entries]
    if not counts:
        return 0

    freq = collections.Counter(counts)
    best_count = max(freq.items(), key=lambda item: (item[1], item[0]))[0]
    return best_count


def summarize_complete_steps(records, nlast: int):
    nparticles = expected_particle_count(records)
    if nparticles == 0:
        return nparticles, []

    complete = []
    for time_key, entries in records.items():
        if len(entries) != nparticles:
            continue

        ip_values = [entry[0] for entry in entries]
        if len(set(ip_values)) != nparticles:
            continue

        sum_fx = math.fsum(entry[1] for entry in entries)
        sum_fy = math.fsum(entry[2] for entry in entries)
        sum_fz = math.fsum(entry[3] for entry in entries)
        complete.append((time_key, sum_fx, sum_fy, sum_fz))

    return nparticles, complete[-nlast:]


def main() -> int:
    args = parse_args()
    logfile = pathlib.Path(args.logfile)

    if not logfile.exists():
        print(f"error: file not found: {logfile}", file=sys.stderr)
        return 1

    if args.last <= 0:
        print("error: --last must be positive", file=sys.stderr)
        return 1

    records = load_records(logfile)
    nparticles, last_steps = summarize_complete_steps(records, args.last)

    if nparticles == 0:
        print("No valid particle force records found.")
        return 1

    print(f"Log file: {logfile}")
    print(f"Expected particles per complete timestep: {nparticles}")

    if not last_steps:
        print("No complete timesteps available yet.")
        return 1

    print("time          sumFx           sumFy           sumFz")
    for time_key, sum_fx, sum_fy, sum_fz in last_steps:
        print(f"{time_key:>10}  {sum_fx: .8e}  {sum_fy: .8e}  {sum_fz: .8e}")

    if len(last_steps) >= 2:
        first = last_steps[0][3]
        last = last_steps[-1][3]
        rel_change = 0.0 if first == 0.0 else (last - first) / first
        abs_change = 0.0 if first == 0.0 else (last - first)
        print(f"relative change in sumFz over window: {rel_change: .6e}")
        print(f"absolute change in sumFz over window: {abs_change: .6e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
