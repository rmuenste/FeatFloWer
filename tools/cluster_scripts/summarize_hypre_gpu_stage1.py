#!/usr/bin/env python3
"""Summarize HYPRE_TIMING and FF_TIMING records from solver benchmark runs.

Each case directory below <result_root>/runs is named
    [<app>-]<backend>-np<ranks>-type<solver_type>-rep<repetition>
and contains a run.log. One TSV row is emitted per (case, solver label):
solver labels are the values of the solver= field in the timing records,
e.g. gmres-amg / amg (Hypre), crs-p-t2 / crs-u-t5 (coarse dispatch wrappers),
mg-pressure / mg-velocity (whole MG stage solves).

Residual-based pass/fail gates apply only to the Hypre solver labels, whose
solver tolerances are known; the other records carry defect ratios or 0 and
are gated on application exit and finiteness alone.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import sys
from pathlib import Path


CASE_RE = re.compile(
    r"(?:(?P<app>[A-Za-z0-9_]+)-)?(?P<backend>cpu|gpu)"
    r"-np(?P<ranks>\d+)-type(?P<type>\d+)"
    r"(?:-vt(?P<vtype>\d+))?-rep(?P<rep>\d+)"
)
TIMING_RE = re.compile(
    r"(?:HYPRE_TIMING|FF_TIMING) solver=(?P<solver>\S+) call=(?P<call>\d+)"
    r" setup_s=\s*(?P<setup>[0-9.Ee+\-]+)"
    r" solve_s=\s*(?P<solve>[0-9.Ee+\-]+)"
    r" total_s=\s*(?P<total>[0-9.Ee+\-]+)"
    r" iterations=(?P<iterations>\d+)"
    r" rel_res=\s*(?P<residual>\S+)"
)
FORCE_RE = re.compile(r"BenchForce:\s*(?P<values>.*)")

# Residual gates for solvers with a known convergence tolerance.
RESIDUAL_LIMITS = {
    "gmres-amg": 1.0e-5,
    "amg": 1.0e-7,
    "pcg-amg": 1.0e-7,
}

FIELDNAMES = [
    "app",
    "backend",
    "ranks",
    "solver_type",
    "velo_type",
    "repetition",
    "solver",
    "samples",
    "setup_median_s",
    "solve_median_s",
    "total_median_s",
    "iterations_median",
    "max_rel_res",
    "run_exit_code",
    "status",
    "last_bench_force",
    "log",
]


def read_case_env(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key.strip()] = value.strip()
    return values


def parse_log(path: Path, warmup_calls: int) -> list[dict[str, object]]:
    match = CASE_RE.fullmatch(path.parent.name)
    if match is None:
        raise ValueError(f"Unexpected case directory: {path.parent.name}")

    records_by_solver: dict[str, list[dict[str, str]]] = {}
    last_force = ""
    for line in path.read_text(errors="replace").splitlines():
        timing = TIMING_RE.search(line)
        if timing:
            record = timing.groupdict()
            records_by_solver.setdefault(record["solver"], []).append(record)
        force = FORCE_RE.search(line)
        if force:
            last_force = force.group("values").strip()

    case_env = read_case_env(path.parent / "case.env")
    run_exit_code = case_env.get("run_exit_code", "")
    app_failed = bool(run_exit_code) and run_exit_code != "0"

    base_row: dict[str, object] = {
        "app": match.group("app") or "",
        "backend": match.group("backend"),
        "ranks": int(match.group("ranks")),
        "solver_type": int(match.group("type")),
        "velo_type": int(match.group("vtype") or 1),
        "repetition": int(match.group("rep")),
        "run_exit_code": run_exit_code,
        "last_bench_force": last_force,
        "log": str(path),
    }

    if not records_by_solver:
        row = dict(base_row)
        row.update(
            {
                "solver": "",
                "samples": 0,
                "setup_median_s": "",
                "solve_median_s": "",
                "total_median_s": "",
                "iterations_median": "",
                "max_rel_res": "",
                "status": "app_failed" if app_failed else "no_samples",
            }
        )
        return [row]

    rows: list[dict[str, object]] = []
    for solver, records in sorted(records_by_solver.items()):
        timed = [r for r in records if int(r["call"]) > warmup_calls]
        residuals = [float(r["residual"]) for r in records]

        if not all(math.isfinite(value) for value in residuals):
            status = "nonfinite"
            max_residual: object = "NaN"
        else:
            max_residual = max(residuals)
            limit = RESIDUAL_LIMITS.get(solver)
            if limit is not None and max_residual > limit:
                status = "residual_failed"
            else:
                status = "ok"
        if app_failed:
            status = "app_failed"

        row = dict(base_row)
        row.update(
            {
                "solver": solver,
                "samples": len(timed),
                "setup_median_s": "",
                "solve_median_s": "",
                "total_median_s": "",
                "iterations_median": "",
                "max_rel_res": max_residual,
                "status": status,
            }
        )
        if timed:
            row.update(
                {
                    "setup_median_s": statistics.median(
                        float(r["setup"]) for r in timed
                    ),
                    "solve_median_s": statistics.median(
                        float(r["solve"]) for r in timed
                    ),
                    "total_median_s": statistics.median(
                        float(r["total"]) for r in timed
                    ),
                    "iterations_median": statistics.median(
                        int(r["iterations"]) for r in timed
                    ),
                }
            )
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_root", type=Path)
    parser.add_argument("--warmup-calls", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    for path in sorted((args.result_root / "runs").glob("*/run.log")):
        rows.extend(parse_log(path, args.warmup_calls))

    output = args.output.open("w", newline="") if args.output else None
    stream = output if output is not None else sys.stdout
    try:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if output is not None:
            output.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
