"""Partition stage: mesh partitioning via PyPartitioner."""

import glob
import os
import subprocess
import time

from featflower_test.errors import PartitionError
from featflower_test.models import PartitionConfig, StageResult, StageStatus


_SUCCESS_MARKER = "The partitioning was successful!"


def run_partition(
    config: PartitionConfig,
    repo_root: str,
    workdir: str,
    log_path: str,
    dry_run: bool = False,
) -> StageResult:
    """Run mesh partitioning and verify output.

    Parameters
    ----------
    config : PartitionConfig
        Partition section from the test definition.
    repo_root : str
        Repository root directory.
    workdir : str
        Absolute working directory for partitioning (application build dir).
    log_path : str
        Path to write partition log.
    dry_run : bool
        If True, print command but do not execute.

    Returns
    -------
    StageResult
    """
    t0 = time.time()
    metadata = {}

    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    # Resolve command template with {repo_root} and {workdir}
    cmd_str = config.command.format(repo_root=repo_root, workdir=workdir)

    # Set LD_LIBRARY_PATH to include workdir (for libmetis.so)
    env = os.environ.copy()
    ld_path = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = workdir + (":" + ld_path if ld_path else "")

    metadata["command"] = cmd_str
    metadata["workdir"] = workdir

    try:
        with open(log_path, "w") as log_fh:
            log_fh.write(f"=== partition ===\n$ {cmd_str}\nworkdir: {workdir}\n\n")
            log_fh.flush()

            if dry_run:
                log_fh.write("[DRY RUN] skipping execution\n")
                return StageResult(
                    stage="partition",
                    status=StageStatus.PASSED,
                    duration_s=time.time() - t0,
                    log_path=log_path,
                    metadata=metadata,
                )

            result = subprocess.run(
                cmd_str.split(),
                cwd=workdir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
            )

            log_fh.write(result.stdout)
            log_fh.write(f"\n--- exit code: {result.returncode} ---\n")
            log_fh.flush()

            if result.returncode != 0:
                raise PartitionError(
                    f"Partitioner failed (rc={result.returncode})",
                    stage="partition",
                    exit_code=result.returncode,
                    log_path=log_path,
                )

            # Check for success marker (from part_main.py:172)
            if _SUCCESS_MARKER not in result.stdout:
                raise PartitionError(
                    f"Partitioner did not report success marker: '{_SUCCESS_MARKER}'",
                    stage="partition",
                    log_path=log_path,
                )

            # Verify partition file count
            if config.expected_partition_count > 0:
                _verify_partition_count(workdir, config.expected_partition_count, log_path)

    except PartitionError as exc:
        return StageResult(
            stage="partition",
            status=StageStatus.FAILED,
            duration_s=time.time() - t0,
            log_path=log_path,
            error=exc.to_test_error(),
            metadata=metadata,
        )

    return StageResult(
        stage="partition",
        status=StageStatus.PASSED,
        duration_s=time.time() - t0,
        log_path=log_path,
        metadata=metadata,
    )


def _verify_partition_count(workdir: str, expected: int, log_path: str) -> None:
    """Check that the expected number of GRID*.tri files were created."""
    # Look in _mesh/*/sub0001/ for GRID*.tri files
    pattern = os.path.join(workdir, "_mesh", "*", "sub0001", "GRID*.tri")
    grid_files = sorted(glob.glob(pattern))

    # Filter out the base GRID.tri (keep only numbered GRID0001.tri etc.)
    numbered = [f for f in grid_files if os.path.basename(f) != "GRID.tri"]

    if len(numbered) != expected:
        raise PartitionError(
            f"Expected {expected} partition files (GRID*.tri), found {len(numbered)}: "
            f"{[os.path.basename(f) for f in numbered]}",
            stage="partition",
            log_path=log_path,
        )
