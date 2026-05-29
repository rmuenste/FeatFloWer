"""Execution runners (local and SLURM)."""

import abc
from typing import Optional

from featflower_test.models import LaunchConfig, SlurmConfig, StageResult


class Runner(abc.ABC):
    """Abstract base class for simulation runners."""

    @abc.abstractmethod
    def execute(
        self,
        launch: LaunchConfig,
        slurm: Optional[SlurmConfig],
        workdir: str,
        log_path: str,
        level: int,
    ) -> StageResult:
        """Execute a simulation and return the stage result."""
        ...
