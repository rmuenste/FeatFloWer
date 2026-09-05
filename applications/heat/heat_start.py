#!/usr/bin/env python
# vim: set filetype=python
"""
A python launcher script for a FeatFloWer application
"""
import os
import shutil
import sys
from pathlib import Path
try:
    sys.path.append(os.environ['FF_PY_HOME'])
except:
    pass
try:
    import partitioner
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
    import partitioner
import getopt
import subprocess
import re
try:
    from e3d_layout import RunLayout, absolute_from_invocation, CASE_SEED_FILES
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'e3d_scripts'))
    from e3d_layout import RunLayout, absolute_from_invocation, CASE_SEED_FILES

HEAT_DIRS = ('_data', '_mesh', '_vtk', '_dump', 'start')
HEAT_SEEDS = CASE_SEED_FILES + (Path('_data/q2p1_param.dat'),)

#===============================================================================
#                          setup the case folder 
#===============================================================================
def checked_case_path(case, relative):
    """Reject runtime symlinks before writing or removing case data."""
    path = case
    for part in Path(relative).parts:
        path = path / part
        if path.is_symlink():
            raise ValueError(f'Runtime path must not be a symlink: {path}')
    return path


def paths_overlap(first, second):
    first, second = first.resolve(), second.resolve()
    return first == second or first in second.parents or second in first.parents

#===============================================================================
#                      Function: Usage
#===============================================================================
def usage():
    print("Usage: heat_start.py [options]")
    print("Where options can be:")
    print("[-h, --help]: prints this message")
    print("[-n, --num-processors]: defines the number of parallel jobs to be used")
    print("[-f, --in-folder]: input folder containing heat.s3d and optional meshDir")
    print("[-C, --case]: case directory (default: current directory)")
    print("[-u, --use-srun]: launch with srun instead of mpirun")
    print("Geometry paths are case-relative or absolute; geometry is not copied.")

def parse_partition_format(param_file):
    """
    Returns 'legacy' or 'json' based on SimPar@PartitionFormat entry.
    Defaults to legacy on missing/invalid values.
    """
    fmt = "legacy"
    try:
        with open(param_file, "r") as f:
            for raw_line in f:
                line = raw_line.split("!")[0].strip()
                if not line:
                    continue
                lower_line = line.lower()
                if lower_line.startswith("simpar@partitionformat"):
                    parts = line.split("=", 1)
                    if len(parts) == 2:
                        candidate = parts[1].strip().strip('"').strip("'").lower()
                        if candidate in ("legacy", "json"):
                            fmt = candidate
                    break
    except OSError:
        pass
    return fmt

def detect_node_count():
    """
    Returns the number of compute nodes assigned to this run.
    Tries scheduler env vars first, then host/rank files.
    """
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
            with open(hostfile, "r") as f:
                for raw in f:
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
            with open(rankfile, "r") as f:
                for raw in f:
                    match = pattern.search(raw)
                    if match:
                        hosts.add(match.group(1))
            if hosts:
                return len(hosts)
        except OSError:
            pass

    return 0

def main(argv):
    inputFile = '_data/meshDir/file.prj'
    inputCaseFolder = ''
    useSrun = False
    caseDir = ''

    numProcessors = -1

    try:
        opts, args = getopt.getopt(argv, "huf:n:C:",
                                ["help", "use-srun", "in-folder=", "num-processors=", "case="])
    except getopt.GetoptError as error:
        raise ValueError(str(error)) from error
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ("-n", "--num-processors"):
            numProcessors = int(arg)
        elif opt in ("-f", "--in-folder"):
            inputCaseFolder = arg
        elif opt in ('-u', '--use-srun'):
            useSrun = True
        elif opt in ('-C', '--case'):
            caseDir = arg
        else:
            usage()
            sys.exit(2)

    if args or not inputCaseFolder or numProcessors < 2:
        raise ValueError('Specify -f INPUT and -n RANKS (at least 2); no positional arguments')
    invocation = Path.cwd()
    layout = RunLayout(
        case_dir=absolute_from_invocation(caseDir, invocation) or invocation,
        install_dir_env='FF_HEAT_HOME', runtime_dirs=HEAT_DIRS, seed_files=HEAT_SEEDS)
    direct = Path(absolute_from_invocation(inputCaseFolder, invocation))
    project = direct if direct.exists() else layout.resolve_input(inputCaseFolder)
    project = project.resolve()
    inputCaseFolder = str(project)
    inputCaseFile = project / 'heat.s3d'
    if not inputCaseFile.is_file():
        raise ValueError(f'Missing input configuration: {inputCaseFile}')
    mesher, solver = layout.exe('s3d_mesher'), layout.exe('heat')
    for executable in (mesher, solver):
        if not os.access(executable, os.X_OK):
            raise ValueError(f'Not executable: {executable}')
    for relative in HEAT_DIRS + tuple(str(p) for p in HEAT_SEEDS) + ('_data/heat.s3d',):
        checked_case_path(layout.case_dir, relative)
    mesh = checked_case_path(layout.case_dir, '_data/meshDir')
    if ((project / 'meshDir').exists() and paths_overlap(project / 'meshDir', mesh)
            or paths_overlap(inputCaseFile, mesh)):
        raise ValueError('Input mesh/configuration overlaps runtime _data/meshDir')
    for seed in HEAT_SEEDS:
        if not (layout.case_dir / seed).is_file() and not (layout.install_dir / seed).is_file():
            raise ValueError(f'Missing runtime default: {seed}')
    for name in ('HOSTFILE', 'PBS_NODEFILE', 'OMPI_MCA_rankfile'):
        if os.environ.get(name):
            os.environ[name] = absolute_from_invocation(os.environ[name], invocation)
    print(layout.describe())
    if (layout.case_dir / '_data/q2p1_param.dat').exists():
        print('[layout] retaining case _data/q2p1_param.dat')
    layout.prepare_case()
    layout.enter_case()
    destination = layout.case_dir / '_data/heat.s3d'
    if not destination.exists() or not os.path.samefile(inputCaseFile, destination):
        shutil.copyfile(inputCaseFile, destination)
    if mesh.exists():
        shutil.rmtree(mesh)
    mesher_status = subprocess.call([mesher, '-a', 'heat'])
    if mesher_status:
        print(f'Mesher exited with status {mesher_status}; checking mesh/fallback')

    # Check if mesh was generated, if not check if there is one provided in the case folder. If not EXIT!
    if not os.path.exists("_data/meshDir"):
      if os.path.exists(inputCaseFolder + "/meshDir"):
        shutil.copytree(inputCaseFolder + "/meshDir","_data/meshDir")
      else:
        print("Error: No mesh automatically generated and no <meshDir> " + 
              "folder present the case folder " + inputCaseFolder)
        sys.exit(2)
      
    if not Path(inputFile).is_file():
        raise ValueError(f'Missing mesh project: {inputFile}')

    # Call the partitioner
    param_file = Path("_data") / "q2p1_param.dat"
    partition_format = parse_partition_format(str(param_file))
    node_groups = detect_node_count()
    if node_groups <= 0:
        node_groups = 1
    print("Partitioner configuration: %d worker parts, %d node group(s), format '%s'" %
          (numProcessors - 1, node_groups, partition_format))
    try:
        partitioner.partition(
            numProcessors - 1,
            1,
            node_groups,
            "NEWFAC",
            str(inputFile),
            partition_format=partition_format
        )
    except TypeError:
        # Legacy partitioner versions do not expose the partition_format keyword.
        partitioner.partition(numProcessors - 1, 1, node_groups, "NEWFAC", str(inputFile))

    # Configure the launch command and start a simulation
    launchCommand = ""
    if useSrun:
        launchCommand = ['srun', solver]
    else:
        launchCommand = ['mpirun', '-np', str(numProcessors), solver]

    # Start the simulation as a subprocess
    exitCode = subprocess.call(launchCommand)
    if exitCode != 0:
      sys.exit(exitCode)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except (ValueError, OSError) as error:
        print(f'Error: {error}', file=sys.stderr)
        sys.exit(2)
