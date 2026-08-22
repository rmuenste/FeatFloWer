#!/usr/bin/env python3
"""Partition and launch the q2p1_cc cylinder benchmark."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run q2p1_cc with one master rank and N-1 worker ranks."
    )
    parser.add_argument("-n", "--num-processors", type=int, default=4)
    parser.add_argument("--project", default="_adc/2D_FAC/2Dbench.prj")
    parser.add_argument("--mesh-folder", default="NEWFAC")
    parser.add_argument("--submeshes", type=int, default=1)
    parser.add_argument("--skip-partition", action="store_true")
    parser.add_argument("--log", default="run_q2p1_cc.log")
    return parser.parse_args()


def find_partitioner(script_path: Path, run_dir: Path) -> Path:
    candidates = [
        run_dir / "PyPartitioner.py",
        script_path.parents[2] / "tools" / "PyPartitioner.py",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("PyPartitioner.py was not found in the run or source tree")


def final_force(log_path: Path) -> str | None:
    match = None
    pattern = re.compile(r"^BenchForce:\s+(.+)$")
    with log_path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            current = pattern.match(line)
            if current:
                match = current.group(1).strip()
    return match


def main() -> int:
    args = parse_args()
    if args.num_processors < 2:
        raise SystemExit("q2p1_cc needs one master rank and at least one worker rank")

    run_dir = Path.cwd()
    executable = run_dir / "q2p1_cc"
    parameter_file = run_dir / "_data" / "q2p1_param.dat"
    project = run_dir / args.project
    for required in (executable, parameter_file, project):
        if not required.exists():
            raise FileNotFoundError(f"required runtime file is missing: {required}")

    environment = os.environ.copy()
    old_library_path = environment.get("LD_LIBRARY_PATH", "")
    environment["LD_LIBRARY_PATH"] = f"{run_dir}:{old_library_path}"

    if not args.skip_partition:
        partitioner = find_partitioner(Path(__file__).resolve(), run_dir)
        subprocess.run(
            [
                sys.executable,
                str(partitioner),
                str(args.num_processors - 1),
                "1",
                str(args.submeshes),
                args.mesh_folder,
                args.project,
            ],
            cwd=run_dir,
            env=environment,
            check=True,
        )

    log_path = run_dir / args.log
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            ["mpirun", "-np", str(args.num_processors), str(executable)],
            cwd=run_dir,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )

    force = final_force(log_path)
    if force:
        print(f"Final BenchForce: {force}")
    print(f"Log: {log_path}")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
