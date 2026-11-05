#!/usr/bin/env python3
"""Export the ten Cate E1-E4 settling traces of the certified campaign runs for the website.

Reads ``particle_force.log`` of each ``q2p1_dns_rundir_e<case>_l<level>`` rundir (one line
per step: ``time ip fx fy fz tx ty tz px py pz vx vy vz``) and writes, per case and level,
two headerless two-column text files in the format the ff-redesign sedimentation
converter consumes:

    vel<CASE>_L<n>.txt   time [s]   v_z [m/s]      (negative = settling)
    pos<CASE>_L<n>.txt   time [s]   z_centre [m]   (raw; the site normalises (z-0.0075)/0.015)

plus ``campaign_runs.csv`` with the provenance of every exported series (rundir, Slurm job,
binary, step count, peak settling velocity) so the site can cite the run it plots.

Usage:
    python3 tools/tencate_export_series.py [--root <FeatFloWer worktree>] --out <dir>

Only rundirs that exist and have finished (MaxNumStep lines in the log) are exported; a
missing level is reported, not invented.
"""
import argparse
import csv
import glob
import os
import re
import sys

CASES = ["E1", "E2", "E3", "E4"]
LEVELS = [2, 3, 4]
COL_T, COL_PZ, COL_VZ = 0, 10, 13


def read_log(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            rows.append((float(f[COL_T]), float(f[COL_PZ]), float(f[COL_VZ])))
    return rows


def deck_value(rundir, key):
    with open(os.path.join(rundir, "_data", "q2p1_param.dat")) as fh:
        for line in fh:
            m = re.match(rf"\s*{re.escape(key)}\s*=\s*(\S+)", line)
            if m:
                return m.group(1)
    return None


def provenance(rundir):
    job = sorted(glob.glob(os.path.join(rundir, "slurm-*.out")))
    jobid = re.search(r"slurm-(\d+)\.out", job[-1]).group(1) if job else ""
    binary = ""
    sb = os.path.join(rundir, "job.sbatch")
    if os.path.exists(sb):
        with open(sb) as fh:
            for line in fh:
                if line.startswith("mpirun"):
                    binary = line.split()[3]
    return jobid, binary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    prov = []
    for case in CASES:
        for lvl in LEVELS:
            rundir = os.path.join(args.root, f"q2p1_dns_rundir_{case.lower()}_l{lvl}")
            log = os.path.join(rundir, "particle_force.log")
            if not os.path.exists(log):
                print(f"{case} L{lvl}: no rundir/log at {rundir} - skipped", file=sys.stderr)
                continue
            rows = read_log(log)
            nsteps = int(deck_value(rundir, "SimPar@MaxNumStep"))
            if len(rows) < nsteps:
                print(f"{case} L{lvl}: {len(rows)} of {nsteps} steps in the log - unfinished, skipped",
                      file=sys.stderr)
                continue
            rows = rows[:nsteps]
            for prefix, col in (("vel", 2), ("pos", 1)):
                with open(os.path.join(args.out, f"{prefix}{case}_L{lvl}.txt"), "w") as fh:
                    for r in rows:
                        fh.write(f"{r[0]:.6f} {r[col]:.9g}\n")
            tpk, _, vpk = min(rows, key=lambda r: r[2])
            jobid, binary = provenance(rundir)
            prov.append(dict(case=case, level=f"L{lvl}", rundir=os.path.basename(rundir),
                             slurm_job=jobid, binary=binary, steps=len(rows),
                             dt_s=deck_value(rundir, "SimPar@TimeStep"),
                             t_end_s=f"{rows[-1][0]:.3f}", peak_vz_m_s=f"{vpk:.6f}",
                             t_peak_s=f"{tpk:.3f}"))
            print(f"{case} L{lvl}: {len(rows)} steps, peak {vpk:.5f} m/s at t={tpk:.3f} s")

    with open(os.path.join(args.out, "campaign_runs.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prov[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(prov)
    print(f"{len(prov)} series written to {args.out}")


if __name__ == "__main__":
    main()
