"""Run stage: parameter patching and simulation launch delegation."""

import os
import re
import time
from typing import Optional

from featflower_test.errors import RunError
from featflower_test.models import (
    LaunchConfig,
    LevelConfig,
    RunConfig,
    SlurmConfig,
    StageResult,
    StageStatus,
)


def patch_parameter_file(
    file_in: str,
    file_out: str,
    level: int,
) -> None:
    """Copy parameter file, replacing SimPar@MaxMeshLevel.

    Mirrors legacy ``moveAndSetLevel`` from q2p1_ctest_start.py:139-145.
    """
    max_level_str = f"SimPar@MaxMeshLevel = {level}"
    pattern = re.compile(r"^\s*SimPar@MaxMeshLevel\s*=.*", re.MULTILINE)

    with open(file_in, "r") as fh:
        content = fh.read()

    new_content = pattern.sub(max_level_str, content)

    os.makedirs(os.path.dirname(file_out), exist_ok=True)
    with open(file_out, "w") as fh:
        fh.write(new_content)


def run_simulation(
    config: RunConfig,
    repo_root: str,
    log_dir: str,
    runner: Optional[object] = None,
    dry_run: bool = False,
) -> list:
    """Execute the simulation for each configured level.

    Parameters
    ----------
    config : RunConfig
        Run section from the test definition.
    repo_root : str
        Repository root directory.
    log_dir : str
        Directory for stage log files.
    runner : Runner, optional
        Runner instance (local or SLURM). If None, uses local runner.
    dry_run : bool
        If True, print commands but do not execute.

    Returns
    -------
    list of StageResult
        One result per level.
    """
    results = []
    levels = config.levels or LevelConfig()
    launch = config.launch or LaunchConfig()
    workdir = os.path.join(repo_root, launch.workdir) if launch.workdir else repo_root

    for i in range(levels.count):
        level = levels.start + i
        log_path = os.path.join(log_dir, f"run-level-{level}.log")
        os.makedirs(os.path.dirname(log_path), exist_ok=True)

        t0 = time.time()

        try:
            # Patch parameter file for this level
            if levels.parameter_patch:
                file_in = os.path.join(
                    repo_root, levels.parameter_patch.get("file_in", "")
                )
                file_out = os.path.join(
                    workdir, levels.parameter_patch.get("file_out", "_data/q2p1_param.dat")
                )
                if dry_run and not os.path.isfile(file_in):
                    pass  # skip patching in dry-run when source is absent
                else:
                    patch_parameter_file(file_in, file_out, level)

            if dry_run:
                with open(log_path, "w") as fh:
                    fh.write(f"[DRY RUN] level={level}\n")
                    fh.write(f"workdir: {workdir}\n")
                    fh.write(f"command: mpirun -np {launch.mpi_ranks} {launch.command}\n")
                results.append(StageResult(
                    stage=f"run-level-{level}",
                    status=StageStatus.PASSED,
                    duration_s=time.time() - t0,
                    log_path=log_path,
                    metadata={"level": level, "dry_run": True},
                ))
                continue

            # Delegate to runner
            if runner is not None:
                stage_result = runner.execute(
                    launch=launch,
                    slurm=config.slurm,
                    workdir=workdir,
                    log_path=log_path,
                    level=level,
                )
                results.append(stage_result)
                if stage_result.status == StageStatus.FAILED:
                    break
            else:
                # Fallback: import local runner
                from featflower_test.runners.local import LocalRunner
                local = LocalRunner()
                stage_result = local.execute(
                    launch=launch,
                    slurm=config.slurm,
                    workdir=workdir,
                    log_path=log_path,
                    level=level,
                )
                results.append(stage_result)
                if stage_result.status == StageStatus.FAILED:
                    break

        except (OSError, FileNotFoundError) as exc:
            results.append(StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=RunError(
                    str(exc),
                    stage=f"run-level-{level}",
                ).to_test_error(),
                metadata={"level": level},
            ))
            break

        except RunError as exc:
            results.append(StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=exc.to_test_error(),
                metadata={"level": level},
            ))
            break

    return results
