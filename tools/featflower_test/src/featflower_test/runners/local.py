"""Local runner: mpirun via subprocess."""

import os
import subprocess
import time
from typing import Optional

from featflower_test.errors import RunError, TimeoutError
from featflower_test.models import LaunchConfig, SlurmConfig, StageResult, StageStatus
from featflower_test.runners import Runner


class LocalRunner(Runner):
    """Execute simulations locally via mpirun."""

    def __init__(self, timeout: int = 3600):
        self.timeout = timeout

    def execute(
        self,
        launch: LaunchConfig,
        slurm: Optional[SlurmConfig],
        workdir: str,
        log_path: str,
        level: int,
    ) -> StageResult:
        t0 = time.time()
        cmd = ["mpirun", "-np", str(launch.mpi_ranks)] + launch.command.split()
        output_file = f"simulation_output_level_{level}.log"
        output_path = os.path.join(workdir, output_file)

        os.makedirs(os.path.dirname(log_path), exist_ok=True)

        try:
            with open(log_path, "w") as log_fh:
                log_fh.write(f"=== run level={level} ===\n")
                log_fh.write(f"$ {' '.join(cmd)}\nworkdir: {workdir}\n\n")
                log_fh.write(f"capturing stdout: {output_path}\n\n")
                log_fh.flush()

            with open(output_path, "w") as output_fh:
                result = subprocess.run(
                    cmd,
                    cwd=workdir,
                    stdout=output_fh,
                    stderr=subprocess.STDOUT,
                    timeout=self.timeout,
                )

            with open(log_path, "a") as log_fh:
                _append_captured_output(log_fh, output_path)
                log_fh.write(f"\n--- exit code: {result.returncode} ---\n")
                log_fh.flush()

        except subprocess.TimeoutExpired:
            with open(log_path, "a") as log_fh:
                _append_captured_output(log_fh, output_path)
                log_fh.write(f"\n--- timeout after {self.timeout}s ---\n")
            return StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=TimeoutError(
                    f"Simulation timed out after {self.timeout}s",
                    stage=f"run-level-{level}",
                ).to_test_error(),
                metadata={"level": level, "output_file": output_file},
            )

        if result.returncode != 0:
            return StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=RunError(
                    f"Simulation failed (rc={result.returncode})",
                    stage=f"run-level-{level}",
                    exit_code=result.returncode,
                    log_path=log_path,
                ).to_test_error(),
                metadata={"level": level, "output_file": output_file},
            )

        return StageResult(
            stage=f"run-level-{level}",
            status=StageStatus.PASSED,
            duration_s=time.time() - t0,
            log_path=log_path,
            metadata={"level": level, "output_file": output_file},
        )


def _append_captured_output(log_fh, output_path: str) -> None:
    """Mirror solver stdout into the stage log after preserving metric output."""
    log_fh.write("--- captured stdout ---\n")
    try:
        with open(output_path, "r") as output_fh:
            for line in output_fh:
                log_fh.write(line)
    except OSError as exc:
        log_fh.write(f"<could not read captured stdout: {exc}>\n")
