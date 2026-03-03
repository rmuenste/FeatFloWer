#!/usr/bin/env python3
"""
Helper launcher for q2p1_fc_ext.

It updates the local parameter file, partitions the mesh according to the
requested rank count and finally executes the solver via mpirun.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def _import_partitioner():
    search_roots = [SCRIPT_DIR]
    # When running from the source tree the partitioner still lives in tools/
    if len(SCRIPT_DIR.parents) >= 2:
        repo_root = SCRIPT_DIR.parents[2]
        search_roots.append(repo_root / "tools")

    for root in search_roots:
        candidate = root / "partitioner" / "__init__.py"
        if candidate.exists():
            sys.path.insert(0, str(root))
            import partitioner

            return partitioner

    sys.exit(
        "Unable to locate the 'partitioner' module. Make sure the partitioner "
        "directory is available next to this script or in the repository tools folder."
    )


partitioner = _import_partitioner()


def detect_node_count():
    """Return number of compute nodes by inspecting scheduler env vars/host files."""
    for env_var in ("SLURM_STEP_NUM_NODES", "SLURM_JOB_NUM_NODES"):
        val = os.environ.get(env_var)
        if val:
            try:
                nodes = int(val)
                if nodes > 0:
                    return nodes
            except ValueError:
                pass

    nodelist = os.environ.get("SLURM_JOB_NODELIST")
    if nodelist:
        try:
            output = subprocess.check_output(
                ["scontrol", "show", "hostnames", nodelist],
                universal_newlines=True,
            )
            hosts = [line.strip() for line in output.splitlines() if line.strip()]
            if hosts:
                return len(set(hosts))
        except (OSError, subprocess.SubprocessError):
            pass

    hostfile = os.environ.get("HOSTFILE") or os.environ.get("PBS_NODEFILE")
    if hostfile:
        try:
            hosts = set()
            with open(hostfile, "r") as handle:
                for raw in handle:
                    line = raw.split("#", 1)[0].strip()
                    if not line:
                        continue
                    token = line.split()[0]
                    hosts.add(token.split(":", 1)[0])
            if hosts:
                return len(hosts)
        except OSError:
            pass

    rankfile = os.environ.get("OMPI_MCA_rankfile")
    if rankfile:
        pattern = re.compile(r"rank\s+\d+\s*=\s*([^\s]+)")
        try:
            hosts = set()
            with open(rankfile, "r") as handle:
                for raw in handle:
                    match = pattern.search(raw)
                    if match:
                        hosts.add(match.group(1))
            if hosts:
                return len(hosts)
        except OSError:
            pass

    return 0


def parse_partition_format(param_file: Path) -> str:
    fmt = "legacy"
    try:
        with param_file.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.split("!")[0].strip()
                if not line:
                    continue
                lower = line.lower()
                if lower.startswith("simpar@partitionformat"):
                    parts = line.split("=", 1)
                    if len(parts) == 2:
                        candidate = parts[1].strip().strip('"').strip("'").lower()
                        if candidate in ("legacy", "json"):
                            fmt = candidate
                    break
    except OSError:
        pass
    return fmt


def parse_recursive_partitioning(param_file: Path) -> bool:
    default = True
    try:
        with param_file.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.split("!")[0].strip()
                if not line:
                    continue
                lower = line.lower()
                if lower.startswith("simpar@recursivepartitioning"):
                    parts = line.split("=", 1)
                    if len(parts) == 2:
                        token = parts[1].strip().strip('"').strip("'").lower()
                        if token in ("yes", "true", "on"):
                            return True
                        if token in ("no", "false", "off"):
                            return False
                    break
    except OSError:
        pass
    return default


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Partition a project file and launch q2p1_fc_ext."
    )
    parser.add_argument(
        "-n",
        "--number-of-ranks",
        dest="ranks",
        type=int,
        required=True,
        help="Total number of MPI ranks (master + workers).",
    )
    parser.add_argument(
        "-f",
        "--input-file",
        dest="project_file",
        required=True,
        help="Path to the *.prj project file that should be partitioned.",
    )
    parser.add_argument(
        "--node-groups",
        dest="node_groups",
        type=int,
        default=None,
        help="Override the automatically detected number of compute nodes "
        "used for recursive (two-level) partitioning.",
    )
    args = parser.parse_args()

    if args.ranks < 2:
        parser.error("number of ranks must be at least 2 (one master + workers).")

    return args


def update_simulation_parameters(project_entry: str, mesh_folder: str) -> None:
    param_file = SCRIPT_DIR / "_data" / "q2p1_param.dat"
    if not param_file.exists():
        sys.exit(f"Parameter file '{param_file}' not found.")

    replacements = {
        "SimPar@ProjectFile": f'SimPar@ProjectFile = "{project_entry}"\n',
        "SimPar@MeshFolder": f'SimPar@MeshFolder = "{mesh_folder}"\n',
    }
    seen_keys = {key: False for key in replacements}
    updated_lines = []

    with param_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.lstrip()
            replaced_line = False
            for key, new_line in replacements.items():
                if stripped.startswith(key):
                    updated_lines.append(new_line)
                    seen_keys[key] = True
                    replaced_line = True
                    break
            if not replaced_line:
                updated_lines.append(line)

    for key, new_line in replacements.items():
        if not seen_keys[key]:
            updated_lines.append(new_line)

    with param_file.open("w", encoding="utf-8") as handle:
        handle.writelines(updated_lines)


def main():
    args = parse_arguments()
    mesh_folder = "NEWFAC"
    orig_cwd = Path.cwd()

    os.chdir(SCRIPT_DIR)

    try:
        project_arg = Path(args.project_file)
        if project_arg.is_absolute():
            project_path = project_arg
            project_entry = str(project_arg)
        else:
            project_entry = args.project_file
            project_path = (SCRIPT_DIR / project_arg).resolve()

        if not project_path.exists():
            sys.exit(f"Project file '{args.project_file}' not found (looked in {project_path}).")

        update_simulation_parameters(project_entry, mesh_folder)
        param_file = SCRIPT_DIR / "_data" / "q2p1_param.dat"
        partition_format = parse_partition_format(param_file)
        recursive_mode = parse_recursive_partitioning(param_file)

        node_groups = 1
        if recursive_mode:
            if args.node_groups is not None:
                node_groups = max(1, args.node_groups)
            else:
                node_groups = detect_node_count()
                if node_groups <= 0:
                    node_groups = 1

        worker_partitions = args.ranks - 1
        print(
            f"[fc_ext_start] Partitioning '{project_entry}' into "
            f"{worker_partitions} parts using {node_groups} subdomains per node "
            f"(mesh folder '{mesh_folder}', format '{partition_format}', "
            f"recursive {'enabled' if recursive_mode else 'disabled'})."
        )
        partitioner.partition(
            worker_partitions,
            1,
            node_groups,
            mesh_folder,
            project_entry,
            partition_format=partition_format,
        )

        cmd = ["mpirun", "-np", str(args.ranks), "./q2p1_fc_ext"]
        print(f"[fc_ext_start] Launching solver: {' '.join(cmd)}")
        subprocess.check_call(cmd)
    finally:
        os.chdir(orig_cwd)


if __name__ == "__main__":
    main()
