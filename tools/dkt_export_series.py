#!/usr/bin/env python3
"""Export DKT (drafting-kissing-tumbling) trajectory series from DNS rundirs.

Reads the per-step ``particle_force.log`` written by ``applications/q2p1_dkt``
and emits small two-column ASCII series suitable for publication on the
FeatFloWer benchmark website. Companion to ``tools/d21_plot_brenner.py`` and
``tools/compare_tencate.py``.

``particle_force.log`` columns (one line per particle per step)::

    time ip fx fy fz tx ty tz px py pz vx vy vz

The DKT case has exactly two particles: ip=1 is the leading (lower) sphere and
ip=2 is the trailing (upper) sphere.

Derived quantities, all nondimensional with sphere diameter d = 1:

``tilt``
    Angle of the pair axis away from vertical,
    ``theta = arccos(dz / |dx|)`` with ``dx = x_trailer - x_leader``.
    0 deg is the drafting/kissing column, 90 deg a horizontal pair, and
    values above 90 deg mean the spheres have exchanged roles.
``separation``
    Centre-to-centre distance in diameters. Exactly 1.0 during contact.
``vz_leader`` / ``vz_trailer``
    Vertical velocity of each sphere.
``xz_leader`` / ``xz_trailer``
    Trajectory in the x-z plane (the motion is planar to ~1e-9 in y).

Contact threshold: separation <= 1.005 d. This is the convention recorded in
``docs/md_docs/dns_validation_datasheet.md`` (row ``dkt_tkiss_correction``);
changing it changes the reported kissing time, so it is a named constant.

Usage::

    python3 tools/dkt_export_series.py --out-dir /path/to/ff-redesign/scripts/source-data/dkt
"""

from __future__ import annotations

import argparse
import math
import os
import sys

# Contact detection threshold in diameters (see module docstring).
CONTACT_THRESHOLD = 1.005
# Separation is considered re-established above this value.
SEPARATION_THRESHOLD = 1.02

# Website run id -> rundir, relative to the repository root.
RUNS = {
    "fric-dh8": "q2p1_dns_rundir_dkt_offset",
    "fric-dh16": "q2p1_dns_rundir_dkt16_offset",
    "nofric": "q2p1_dns_rundir_dkt_nofric_long",
    "axisym": "q2p1_dns_rundir_dkt_smoke",
}

# Series name -> index into the per-step record produced by `derive`.
SERIES = {
    "tilt": 1,
    "separation": 2,
    "vz_leader": 7,
    "vz_trailer": 8,
}
TRAJECTORIES = {"xz_leader": (3, 4), "xz_trailer": (5, 6)}


def read_force_log(path):
    """Return [(time, leader_cols, trailer_cols), ...] sorted by time."""
    steps = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            time = float(parts[0])
            ip = int(parts[1])
            steps.setdefault(time, {})[ip] = [float(v) for v in parts[2:]]
    rows = []
    for time in sorted(steps):
        step = steps[time]
        if 1 in step and 2 in step:
            rows.append((time, step[1], step[2]))
    return rows


def derive(rows):
    """Map raw records onto (t, tilt, separation, x1, z1, x2, z2, vz1, vz2)."""
    series = []
    for time, leader, trailer in rows:
        # Columns after `time ip`: fx fy fz tx ty tz px py pz vx vy vz
        p1, p2 = leader[6:9], trailer[6:9]
        v1, v2 = leader[9:12], trailer[9:12]
        dx, dy, dz = p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]
        distance = math.sqrt(dx * dx + dy * dy + dz * dz)
        if distance > 0.0:
            tilt = math.degrees(math.acos(max(-1.0, min(1.0, dz / distance))))
        else:
            tilt = 0.0
        series.append((time, tilt, distance, p1[0], p1[2], p2[0], p2[2], v1[2], v2[2]))
    return series


def phase_events(series):
    """Locate the kissing, separation-onset and horizontal-pair events."""
    events = {}
    kiss = next((s for s in series if s[2] <= CONTACT_THRESHOLD), None)
    if kiss:
        events["t_kiss"] = kiss[0]
        events["tilt_at_kiss"] = kiss[1]
        onset = next(
            (s for s in series if s[0] > kiss[0] and s[2] > SEPARATION_THRESHOLD), None
        )
        if onset:
            events["t_separation_onset"] = onset[0]
    horizontal = next((s for s in series if s[1] >= 90.0), None)
    if horizontal:
        events["t_horizontal"] = horizontal[0]
    final = series[-1]
    events["t_end"] = final[0]
    events["tilt_end"] = final[1]
    events["separation_end"] = final[2]
    return events


def write_series(path, pairs):
    with open(path, "w") as handle:
        for x, y in pairs:
            handle.write(f"{x:.6f} {y:.8g}\n")


def export_run(run_id, rundir, out_dir, stride):
    log = os.path.join(rundir, "particle_force.log")
    if not os.path.isfile(log):
        raise SystemExit(f"missing force log: {log}")
    series = derive(read_force_log(log))
    kept = series[::stride]

    for name, index in SERIES.items():
        write_series(
            os.path.join(out_dir, f"{run_id}_{name}.dat"),
            [(s[0], s[index]) for s in kept],
        )
    for name, (xi, zi) in TRAJECTORIES.items():
        write_series(
            os.path.join(out_dir, f"{run_id}_{name}.dat"),
            [(s[xi], s[zi]) for s in kept],
        )

    events = phase_events(series)
    print(f"{run_id:<12} {len(series):>6} steps -> {len(kept):>5} points")
    summary = "  ".join(f"{k}={v:.4g}" for k, v in events.items())
    print(f"{'':<12} {summary}")
    return events


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out-dir",
        required=True,
        help="directory to write the two-column .dat series into",
    )
    parser.add_argument(
        "--repo-root",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="FeatFloWer repository root holding the DKT rundirs",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=5,
        help="keep every Nth timestep (default 5, i.e. dt=0.025)",
    )
    args = parser.parse_args(argv)

    os.makedirs(args.out_dir, exist_ok=True)
    for run_id, rundir in RUNS.items():
        export_run(
            run_id,
            os.path.join(args.repo_root, rundir),
            args.out_dir,
            args.stride,
        )
    print(f"\nwrote {len(RUNS) * (len(SERIES) + len(TRAJECTORIES))} series to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
