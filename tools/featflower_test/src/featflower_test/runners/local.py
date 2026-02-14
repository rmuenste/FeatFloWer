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

        os.makedirs(os.path.dirname(log_path), exist_ok=True)

        try:
            with open(log_path, "w") as log_fh:
                log_fh.write(f"=== run level={level} ===\n")
                log_fh.write(f"$ {' '.join(cmd)}\nworkdir: {workdir}\n\n")
                log_fh.flush()

                result = subprocess.run(
                    cmd,
                    cwd=workdir,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    timeout=self.timeout,
                )

                log_fh.write(f"\n--- exit code: {result.returncode} ---\n")
                log_fh.flush()

        except subprocess.TimeoutExpired:
            return StageResult(
                stage=f"run-level-{level}",
                status=StageStatus.FAILED,
                duration_s=time.time() - t0,
                log_path=log_path,
                error=TimeoutError(
                    f"Simulation timed out after {self.timeout}s",
                    stage=f"run-level-{level}",
                ).to_test_error(),
                metadata={"level": level},
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
                metadata={"level": level},
            )

        return StageResult(
            stage=f"run-level-{level}",
            status=StageStatus.PASSED,
            duration_s=time.time() - t0,
            log_path=log_path,
            metadata={"level": level},
        )
