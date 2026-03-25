#!/usr/bin/env python3
"""
Summarize one steady DNS drag-calibration case from q2p1_dns_drag output.

This script combines two time-history files:

- ``particle_force.log``:
  Written from the FBM hydrodynamic force computation path.
  Each complete timestep contains one row per particle with force, torque,
  position, and velocity. This script uses only the hydrodynamic force
  components and forms the timestep totals

  - ``sumFx``: total hydrodynamic force in x-direction over all particles
  - ``sumFy``: total hydrodynamic force in y-direction over all particles
  - ``sumFz``: total hydrodynamic force in z-direction over all particles

  In the current cube test case, the imposed forcing is along +z, so
  ``sumFz`` is the main drag/resistance signal. ``sumFx`` and ``sumFy``
  are symmetry checks and should stay small compared to ``sumFz``.

- ``bulk_flow.log``:
  Written from a cubature-based bulk-flow monitor in the CFD field.
  The velocity is projected onto the imposed forcing direction, so the
  logged values are not tied to z in principle, even though the current
  case is z-forced.

  Per timestep the file stores:

  - ``U_sup``:
    Superficial velocity in forcing direction,

    ``U_sup = (1 / V_total) * integral_domain (u · e_force) dV``

    where ``e_force`` is the unit vector in forcing direction. This is the
    bulk transport velocity based on the full box volume.

  - ``U_fluid``:
    Fluid-phase average velocity in forcing direction,

    ``U_fluid = (1 / V_fluid) * integral_domain phi_f * (u · e_force) dV``

    where ``phi_f`` is the local fluid volume fraction reconstructed from
    the FBM occupancy field. Because particles occupy part of the box,
    ``U_fluid`` is typically larger than ``U_sup``.

  - ``fluid_fraction``:
    Global fluid volume fraction,

    ``fluid_fraction = V_fluid / V_total``

    The corresponding solids fraction is

    ``solids_fraction = 1 - fluid_fraction``.

  - ``valid``:
    Indicates that the bulk monitor completed a valid FE/cubature pass for
    that timestep.

  - ``enabled``:
    Indicates that constant forcing was enabled and the monitor was active.

What the script does:

- identifies complete timesteps in ``particle_force.log``
- keeps only valid+enabled timesteps from ``bulk_flow.log``
- matches both files by timestep
- takes the last ``N`` matched timesteps
- computes tail averages and simple drift indicators

Reported summary quantities:

- ``mean_sumFx``, ``mean_sumFy``, ``mean_sumFz``:
  Tail-averaged total hydrodynamic force components on the particle array.

- ``mean_U_sup``:
  Tail-averaged superficial velocity in forcing direction.

- ``mean_U_fluid``:
  Tail-averaged fluid-phase mean velocity in forcing direction.

- ``mean_fluid_fraction``:
  Tail-averaged global fluid fraction. In a fixed-particle geometry this
  should normally be constant in time.

- ``mean_solids_fraction``:
  Complement of the fluid fraction. This is the global particle hold-up /
  volume fraction of solids in the domain.

- ``drag_over_U_sup``:
  ``mean_sumFz / mean_U_sup``

  This is a simple effective resistance measure based on superficial
  velocity. It is often a useful first calibration target because
  superficial velocity is a common bulk quantity in drag closures.

- ``drag_over_U_fluid``:
  ``mean_sumFz / mean_U_fluid``

  This is the corresponding effective resistance measure based on
  fluid-phase mean velocity.

Tail drift checks:

- ``rel_drift_sumFz``
- ``rel_drift_U_sup``
- ``rel_drift_U_fluid``

Each is computed as

``(last_value - first_value) / first_value``

over the selected tail window. These numbers are intended as practical
steady-state indicators. Small values mean the run is asymptotically flat
over the averaging window and the resulting summary row is suitable for use
as a DNS calibration sample.

CSV line:

The final ``CSV:...`` line is a compact machine-readable dataset row that can
be copied into a spreadsheet, a fitting script, or a growing drag-calibration
database.
"""

from __future__ import annotations

import argparse
import collections
import math
import pathlib
import statistics
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine particle_force.log and bulk_flow.log from q2p1_dns_drag "
            "and emit one tail-averaged steady-state summary row."
        )
    )
    parser.add_argument(
        "--force-log",
        default="particle_force.log",
        help="Path to particle_force.log",
    )
    parser.add_argument(
        "--bulk-log",
        default="bulk_flow.log",
        help="Path to bulk_flow.log",
    )
    parser.add_argument(
        "--last",
        type=int,
        default=20,
        help="Number of matching tail timesteps to average (default: 20)",
    )
    parser.add_argument(
        "--label",
        default="dns_drag_case",
        help="Case label written into the summary output",
    )
    parser.add_argument(
        "--rho",
        type=float,
        default=1.0,
        help="Fluid density used for Reynolds-number evaluation (default: 1.0)",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=1.0,
        help="Dynamic viscosity used for Reynolds-number evaluation (default: 1.0)",
    )
    parser.add_argument(
        "--box-volume",
        type=float,
        default=1.0,
        help="Total domain volume used to infer particle volume from solids fraction (default: 1.0)",
    )
    return parser.parse_args()


def parse_bool_flag(value: str) -> bool:
    value = value.strip().upper()
    if value == "T":
        return True
    if value == "F":
        return False
    raise ValueError(value)


def load_force_steps(logfile: pathlib.Path):
    records = collections.OrderedDict()

    with logfile.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 14:
                continue

            try:
                time_key = parts[0]
                ip = int(parts[1])
                fx = float(parts[2])
                fy = float(parts[3])
                fz = float(parts[4])
            except ValueError:
                continue

            records.setdefault(time_key, []).append((ip, fx, fy, fz))

    counts = [len(entries) for entries in records.values() if entries]
    if not counts:
        return 0, {}

    expected_particles = max(
        collections.Counter(counts).items(), key=lambda item: (item[1], item[0])
    )[0]

    complete = {}
    for time_key, entries in records.items():
        if len(entries) != expected_particles:
            continue

        ip_values = [entry[0] for entry in entries]
        if len(set(ip_values)) != expected_particles:
            continue

        complete[time_key] = {
            "sumFx": math.fsum(entry[1] for entry in entries),
            "sumFy": math.fsum(entry[2] for entry in entries),
            "sumFz": math.fsum(entry[3] for entry in entries),
        }

    return expected_particles, complete


def load_bulk_steps(logfile: pathlib.Path):
    steps = {}

    with logfile.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            parts = stripped.split()
            if len(parts) < 5:
                continue

            try:
                time_key = f"{float(parts[0]):.6f}"
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

            if not (valid and enabled):
                continue

            steps[time_key] = {
                "U_sup": u_sup,
                "U_fluid": u_fluid,
                "fluid_fraction": fluid_fraction,
            }

    return steps


def rel_drift(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    first = values[0]
    last = values[-1]
    return 0.0 if first == 0.0 else (last - first) / first


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else 0.0


def main() -> int:
    args = parse_args()
    force_log = pathlib.Path(args.force_log)
    bulk_log = pathlib.Path(args.bulk_log)

    if not force_log.exists():
        print(f"error: file not found: {force_log}", file=sys.stderr)
        return 1
    if not bulk_log.exists():
        print(f"error: file not found: {bulk_log}", file=sys.stderr)
        return 1
    if args.last <= 0:
        print("error: --last must be positive", file=sys.stderr)
        return 1

    expected_particles, force_steps = load_force_steps(force_log)
    bulk_steps = load_bulk_steps(bulk_log)

    if expected_particles == 0 or not force_steps:
        print("error: no complete force timesteps found", file=sys.stderr)
        return 1
    if not bulk_steps:
        print("error: no valid bulk-flow timesteps found", file=sys.stderr)
        return 1

    common_times = sorted(set(force_steps).intersection(bulk_steps), key=float)
    if not common_times:
        print("error: no matching timesteps between force and bulk logs", file=sys.stderr)
        return 1

    tail_times = common_times[-args.last :]

    sum_fx = [force_steps[t]["sumFx"] for t in tail_times]
    sum_fy = [force_steps[t]["sumFy"] for t in tail_times]
    sum_fz = [force_steps[t]["sumFz"] for t in tail_times]
    u_sup = [bulk_steps[t]["U_sup"] for t in tail_times]
    u_fluid = [bulk_steps[t]["U_fluid"] for t in tail_times]
    fluid_fraction = [bulk_steps[t]["fluid_fraction"] for t in tail_times]

    mean_sum_fx = mean(sum_fx)
    mean_sum_fy = mean(sum_fy)
    mean_sum_fz = mean(sum_fz)
    mean_u_sup = mean(u_sup)
    mean_u_fluid = mean(u_fluid)
    mean_fluid_fraction = mean(fluid_fraction)

    solids_fraction = 1.0 - mean_fluid_fraction
    drag_over_u_sup = 0.0 if mean_u_sup == 0.0 else mean_sum_fz / mean_u_sup
    drag_over_u_fluid = 0.0 if mean_u_fluid == 0.0 else mean_sum_fz / mean_u_fluid
    mean_drag_per_particle = mean_sum_fz / expected_particles

    particle_volume = 0.0
    sphere_diameter = 0.0
    sphere_radius = 0.0
    if expected_particles > 0 and args.box_volume > 0.0:
        particle_volume = solids_fraction * args.box_volume / expected_particles
        if particle_volume > 0.0:
            sphere_diameter = (6.0 * particle_volume / math.pi) ** (1.0 / 3.0)
            sphere_radius = 0.5 * sphere_diameter

    re_sup = 0.0
    re_fluid = 0.0
    if args.mu != 0.0 and sphere_diameter > 0.0:
        re_sup = args.rho * mean_u_sup * sphere_diameter / args.mu
        re_fluid = args.rho * mean_u_fluid * sphere_diameter / args.mu

    beta_sup = 0.0
    beta_fluid = 0.0
    if particle_volume > 0.0 and mean_u_sup != 0.0:
        beta_sup = mean_drag_per_particle / (particle_volume * mean_u_sup)
    if particle_volume > 0.0 and mean_u_fluid != 0.0:
        beta_fluid = mean_drag_per_particle / (particle_volume * mean_u_fluid)

    print(f"Case label: {args.label}")
    print(f"Force log: {force_log}")
    print(f"Bulk log: {bulk_log}")
    print(f"Matched tail timesteps used: {len(tail_times)}")
    print(f"Expected particles per timestep: {expected_particles}")
    print(f"Time window: {tail_times[0]} -> {tail_times[-1]}")
    print(f"rho                 = {args.rho: .8e}")
    print(f"mu                  = {args.mu: .8e}")
    print(f"box_volume          = {args.box_volume: .8e}")
    print("")
    print("Tail-averaged summary")
    print(f"mean_sumFx           = {mean_sum_fx: .8e}")
    print(f"mean_sumFy           = {mean_sum_fy: .8e}")
    print(f"mean_sumFz           = {mean_sum_fz: .8e}")
    print(f"mean_drag_per_particle = {mean_drag_per_particle: .8e}")
    print(f"mean_U_sup           = {mean_u_sup: .8e}")
    print(f"mean_U_fluid         = {mean_u_fluid: .8e}")
    print(f"mean_fluid_fraction  = {mean_fluid_fraction: .8e}")
    print(f"mean_solids_fraction = {solids_fraction: .8e}")
    print(f"particle_volume      = {particle_volume: .8e}")
    print(f"sphere_diameter      = {sphere_diameter: .8e}")
    print(f"sphere_radius        = {sphere_radius: .8e}")
    print(f"Re_sup               = {re_sup: .8e}")
    print(f"Re_fluid             = {re_fluid: .8e}")
    print(f"beta_sup             = {beta_sup: .8e}")
    print(f"beta_fluid           = {beta_fluid: .8e}")
    print(f"drag_over_U_sup      = {drag_over_u_sup: .8e}")
    print(f"drag_over_U_fluid    = {drag_over_u_fluid: .8e}")
    print("")
    print("Tail drift check")
    print(f"rel_drift_sumFz      = {rel_drift(sum_fz): .8e}")
    print(f"rel_drift_U_sup      = {rel_drift(u_sup): .8e}")
    print(f"rel_drift_U_fluid    = {rel_drift(u_fluid): .8e}")
    print("")
    print(
        "CSV:"
        f"{args.label},"
        f"{tail_times[0]},{tail_times[-1]},{len(tail_times)},{expected_particles},"
        f"{mean_fluid_fraction:.10e},{solids_fraction:.10e},"
        f"{mean_sum_fx:.10e},{mean_sum_fy:.10e},{mean_sum_fz:.10e},"
        f"{mean_drag_per_particle:.10e},{particle_volume:.10e},"
        f"{sphere_diameter:.10e},{sphere_radius:.10e},"
        f"{mean_u_sup:.10e},{mean_u_fluid:.10e},"
        f"{re_sup:.10e},{re_fluid:.10e},"
        f"{beta_sup:.10e},{beta_fluid:.10e},"
        f"{drag_over_u_sup:.10e},{drag_over_u_fluid:.10e},"
        f"{rel_drift(sum_fz):.10e},{rel_drift(u_sup):.10e},{rel_drift(u_fluid):.10e}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
