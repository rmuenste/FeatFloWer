"""SLURM runner: sbatch submission + sacct polling.

Mirrors the legacy ``submitAndObserveSync`` from q2p1_ctest_start.py:93-134.
"""

import os
import subprocess
import tempfile
import time
from typing import List, Optional

from featflower_test.errors import RunError, SlurmError, TimeoutError
from featflower_test.models import LaunchConfig, SlurmConfig, StageResult, StageStatus
from featflower_test.runners import Runner

_POLL_INTERVAL = 10  # seconds, matching legacy

# SLURM states that indicate completion
_TERMINAL_STATES = {
    "COMPLETED", "FAILED", "CANCELLED", "TIMEOUT",
    "NODE_FAIL", "PREEMPTED", "OUT_OF_MEMORY",
}


class SlurmRunner(Runner):
    """Execute simulations via SLURM batch submission."""

    def __init__(self, poll_interval: int = _POLL_INTERVAL):
        self.poll_interval = poll_interval

    def execute(
        self,
        launch: LaunchConfig,
        slurm: Optional[SlurmConfig],
        workdir: str,
        log_path: str,
        level: int,
    ) -> StageResult:
        if slurm is None:
            return StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                error=SlurmError(
                    "SLURM configuration missing in test definition",
                    stage=f"run-level-{level}",
                ).to_test_error(),
            )

        t0 = time.time()
        metadata = {"level": level}

        os.makedirs(os.path.dirname(log_path), exist_ok=True)

        try:
            # Generate job script
            script_content = _generate_job_script(launch, slurm)
            script_path = os.path.join(workdir, f"job_level_{level}.sh")
            with open(script_path, "w") as fh:
                fh.write(script_content)
            metadata["script_path"] = script_path

            # Submit via sbatch
            job_id = _submit_job(script_path, workdir)
            metadata["slurm_job_id"] = job_id

            with open(log_path, "w") as log_fh:
                log_fh.write(f"=== SLURM run level={level} ===\n")
                log_fh.write(f"Job ID: {job_id}\n")
                log_fh.write(f"Script: {script_path}\n\n")
                log_fh.flush()

                # Poll sacct until terminal state
                final_state = _poll_job(job_id, log_fh, self.poll_interval)
                metadata["slurm_state"] = final_state

                log_fh.write(f"\n--- final state: {final_state} ---\n")

        except (SlurmError, TimeoutError) as exc:
            return StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=exc.to_test_error(),
                metadata=metadata,
            )

        if final_state != "COMPLETED":
            error = _map_slurm_state_to_error(final_state, level, log_path)
            return StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=error.to_test_error(),
                metadata=metadata,
            )

        return StageResult(
            stage=f"run-level-{level}",
            status=StageStatus.PASSED,
            duration_s=time.time() - t0,
            log_path=log_path,
            metadata=metadata,
        )


def _generate_job_script(
    launch: LaunchConfig,
    slurm: SlurmConfig,
) -> str:
    """Generate a SLURM batch script."""
    lines = ["#!/bin/bash"]
    lines.append(f"#SBATCH --partition={slurm.partition}")
    lines.append(f"#SBATCH --nodes={slurm.nodes}")
    lines.append(f"#SBATCH --ntasks={slurm.ntasks}")
    lines.append(f"#SBATCH --cpus-per-task={slurm.cpus_per_task}")
    lines.append(f"#SBATCH --time={slurm.time}")
    if slurm.constraint:
        lines.append(f"#SBATCH --constraint={slurm.constraint}")
    if slurm.mem_per_cpu:
        lines.append(f"#SBATCH --mem-per-cpu={slurm.mem_per_cpu}")
    lines.append(f"#SBATCH --job-name=FF-test")
    lines.append("")
    lines.append("set -euo pipefail")
    lines.append("")
    lines.append(f"mpirun -np {launch.mpi_ranks} {launch.command}")
    lines.append("")
    return "\n".join(lines)


def _submit_job(script_path: str, workdir: str) -> str:
    """Submit a job script and return the SLURM job ID."""
    try:
        output = subprocess.check_output(
            ["sbatch", script_path],
            cwd=workdir,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise SlurmError(
            f"sbatch failed (rc={exc.returncode}): {exc.output}",
            stage="run",
            exit_code=exc.returncode,
        )

    # Output format: "Submitted batch job 12345"
    parts = output.strip().split()
    if len(parts) < 4:
        raise SlurmError(
            f"Unexpected sbatch output: {output.strip()}",
            stage="run",
        )
    return parts[-1]


def _poll_job(job_id: str, log_fh, poll_interval: int) -> str:
    """Poll sacct until the job reaches a terminal state.

    Matches legacy ``submitAndObserveSync`` polling pattern.
    """
    sacct_cmd = [
        "sacct", "--format=State",
        f"--jobs={job_id}",
        "--noheader", "--parsable2",
    ]

    while True:
        try:
            output = subprocess.check_output(
                sacct_cmd,
                stderr=subprocess.STDOUT,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            raise SlurmError(
                f"sacct query failed: {exc.output}",
                stage="run",
            )

        states = [s.strip() for s in output.strip().splitlines() if s.strip()]

        if not states:
            state = "PENDING"
        else:
            # Take the first state line (batch job state)
            state = states[0]

        log_fh.write(f"Job {job_id} state: {state}\n")
        log_fh.flush()

        if state in _TERMINAL_STATES:
            return state

        time.sleep(poll_interval)


def _map_slurm_state_to_error(state: str, level: int, log_path: str):
    """Map SLURM terminal state to appropriate error type."""
    stage = f"run-level-{level}"
    if state == "TIMEOUT":
        return TimeoutError(
            f"SLURM job timed out (state={state})",
            stage=stage,
            log_path=log_path,
        )
    elif state in ("FAILED", "NODE_FAIL", "OUT_OF_MEMORY"):
        return RunError(
            f"SLURM job failed (state={state})",
            stage=stage,
            log_path=log_path,
        )
    else:
        return SlurmError(
            f"SLURM job ended with state={state}",
            stage=stage,
            log_path=log_path,
        )
