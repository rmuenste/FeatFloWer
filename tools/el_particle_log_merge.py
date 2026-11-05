#!/usr/bin/env python3
"""Merge the per-rank EL trajectory logs into one DNS-style particle_force.log.

An EL run with SimPar@ELParticleLogFreq > 0 writes el_particles_rank<NNNN>.log
on every worker (locally owned particles only, no collectives). This script
concatenates them, sorts by (time, ip) and writes particle_force.log with the
DNS header, so tools/d32_ladder_analysis.py reads an EL rundir unchanged.

Usage:
  python3 tools/el_particle_log_merge.py <rundir> [--out particle_force.log]
"""
import argparse
import glob
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--out", default="particle_force.log")
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.rundir, "el_particles_rank*.log")))
    if not files:
        sys.exit(f"no el_particles_rank*.log in {a.rundir}")
    rows = []
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 14:
                    continue
                rows.append((float(parts[0]), int(parts[1]), line.rstrip("\n")))
    rows.sort(key=lambda r: (r[0], r[1]))
    out = os.path.join(a.rundir, a.out)
    with open(out, "w") as fh:
        fh.write("# time ip fx fy fz tx ty tz px py pz vx vy vz\n")
        for _, _, line in rows:
            fh.write(line + "\n")
    times = sorted({r[0] for r in rows})
    ids = {r[1] for r in rows}
    per_t = len(rows) / max(len(times), 1)
    print(f"{out}: {len(rows)} rows, {len(times)} snapshots, {len(ids)} particle ids, "
          f"{per_t:.1f} rows/snapshot from {len(files)} rank files")
    if len(ids) != round(per_t):
        print("WARNING: rows per snapshot != distinct ids (a particle logged twice or missing in a snapshot)")


if __name__ == "__main__":
    main()
