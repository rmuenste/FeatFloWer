"""Typed error hierarchy for the FeatFloWer test system."""

from typing import Optional

from featflower_test.models import ErrorCategory, TestError


class FeatFlowerTestError(Exception):
    """Base exception for all FeatFloWer test errors."""

    category = ErrorCategory.INFRA_ERROR

    def __init__(
        self,
        message: str,
        stage: str = "",
        exit_code: Optional[int] = None,
        log_path: Optional[str] = None,
    ):
        super().__init__(message)
        self.stage = stage
        self.exit_code = exit_code
        self.log_path = log_path

    def to_test_error(self) -> TestError:
        return TestError(
            category=self.category,
            stage=self.stage,
            message=str(self),
            exit_code=self.exit_code,
            log_path=self.log_path,
        )


class ConfigurationError(FeatFlowerTestError):
    """YAML schema, missing env vars, or invalid configuration."""
    category = ErrorCategory.CONFIGURATION_ERROR


class SubmoduleError(FeatFlowerTestError):
    """Git submodule sync/update failure."""
    category = ErrorCategory.SUBMODULE_ERROR


class BuildError(FeatFlowerTestError):
    """CMake configure or build failure."""
    category = ErrorCategory.BUILD_ERROR


class PartitionError(FeatFlowerTestError):
    """Mesh partitioning failure."""
    category = ErrorCategory.PARTITION_ERROR


class RunError(FeatFlowerTestError):
    """Simulation execution failure."""
    category = ErrorCategory.RUN_ERROR


class MetricsError(FeatFlowerTestError):
    """Metric extraction or parsing failure."""
    category = ErrorCategory.METRICS_ERROR


class SlurmError(FeatFlowerTestError):
    """SLURM submission or job failure."""
    category = ErrorCategory.SLURM_ERROR


class TimeoutError(FeatFlowerTestError):
    """Execution exceeded time limit."""
    category = ErrorCategory.TIMEOUT
