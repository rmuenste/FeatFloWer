#!/usr/bin/env python3

from __future__ import annotations

import argparse
import pathlib
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read bulk_flow.log from a running q2p1_dns_drag case and report "
            "the last fully available timesteps."
        )
    )
    parser.add_argument(
        "logfile",
        nargs="?",
        default="bulk_flow.log",
        help="Path to bulk_flow.log",
    )
    parser.add_argument(
        "--last",
        type=int,
        default=3,
        help="Number of complete timesteps to report (default: 3)",
    )
    return parser.parse_args()


def parse_bool_flag(value: str) -> bool:
    value = value.strip().upper()
    if value == "T":
        return True
    if value == "F":
        return False
    raise ValueError(value)


def load_records(logfile: pathlib.Path):
    records = []

    with logfile.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 5:
                continue

            try:
                time_value = float(parts[0])
                u_sup = float(parts[1])
                u_fluid = float(parts[2])
                fluid_fraction = float(parts[3])
                if len(parts) >= 6:
                    valid = parse_bool_flag(parts[4])
                    enabled = parse_bool_flag(parts[5])
                else:
                    flags = parts[4].strip().upper()
                    if len(flags) != 2:
                        raise ValueError(flags)
                    valid = parse_bool_flag(flags[0])
                    enabled = parse_bool_flag(flags[1])
            except ValueError:
                continue

            records.append(
                (time_value, u_sup, u_fluid, fluid_fraction, valid, enabled)
            )

    return records


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
    if not records:
        print("No valid bulk flow records found.")
        return 1

    last_steps = records[-args.last :]

    print(f"Log file: {logfile}")
    print("time          U_sup           U_fluid         fluid_frac   valid enabled")
    for time_value, u_sup, u_fluid, fluid_fraction, valid, enabled in last_steps:
        print(
            f"{time_value:10.6f}  {u_sup: .8e}  {u_fluid: .8e}  "
            f"{fluid_fraction: .8e}   {str(valid)[0]}      {str(enabled)[0]}"
        )

    if len(last_steps) >= 2:
        first_u_sup = last_steps[0][1]
        last_u_sup = last_steps[-1][1]
        first_u_fluid = last_steps[0][2]
        last_u_fluid = last_steps[-1][2]

        rel_u_sup = 0.0 if first_u_sup == 0.0 else (last_u_sup - first_u_sup) / first_u_sup
        rel_u_fluid = (
            0.0
            if first_u_fluid == 0.0
            else (last_u_fluid - first_u_fluid) / first_u_fluid
        )

        print(f"relative change in U_sup over window:   {rel_u_sup: .6e}")
        print(f"relative change in U_fluid over window: {rel_u_fluid: .6e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
