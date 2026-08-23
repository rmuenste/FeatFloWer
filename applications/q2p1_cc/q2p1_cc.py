#!/usr/bin/env python3
"""Partition and launch the q2p1_cc cylinder benchmark.

The parameter file _data/q2p1_param.dat is the single source of truth:
the partitioning inputs (project file, mesh folder, submesh count) are
read from its SimPar@ keys, so the mesh that gets partitioned is always
the mesh the solver will load.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run q2p1_cc with one master rank and N-1 worker ranks. "
        "Project file, mesh folder, and submesh count are taken from "
        "_data/q2p1_param.dat."
    )
    parser.add_argument("-n", "--num-processors", type=int, default=4)
    parser.add_argument("--skip-partition", action="store_true")
    parser.add_argument("--log", default="run_q2p1_cc.log")
    return parser.parse_args()


def read_deck_values(parameter_file: Path) -> dict[str, str]:
    """Read the SimPar@ keys that define the mesh from the parameter file."""
    wanted = {"ProjectFile", "MeshFolder", "SubMeshNumber"}
    pattern = re.compile(r"^\s*SimPar@(\w+)\s*=\s*([^!]*)")
    values: dict[str, str] = {}
    with parameter_file.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = pattern.match(line)
            if match and match.group(1) in wanted:
                values[match.group(1)] = match.group(2).strip().strip('"').strip("'")
    missing = sorted(wanted - values.keys())
    if missing:
        raise SystemExit(
            f"{parameter_file} does not define SimPar@{', SimPar@'.join(missing)}"
        )
    return values


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
    for required in (executable, parameter_file):
        if not required.exists():
            raise FileNotFoundError(f"required runtime file is missing: {required}")

    deck = read_deck_values(parameter_file)
    project = run_dir / deck["ProjectFile"]
    if not project.exists():
        raise FileNotFoundError(
            f"SimPar@ProjectFile from {parameter_file} is missing: {project}"
        )
    print(
        f"Partitioning per parameter file: project={deck['ProjectFile']}, "
        f"mesh folder={deck['MeshFolder']}, submeshes={deck['SubMeshNumber']}",
        flush=True,
    )

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
                deck["SubMeshNumber"],
                deck["MeshFolder"],
                deck["ProjectFile"],
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
