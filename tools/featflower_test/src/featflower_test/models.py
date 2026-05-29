"""Core data models for the FeatFloWer test system."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class StageStatus(Enum):
    """Status of a pipeline stage."""
    PENDING = "pending"
    RUNNING = "running"
    PASSED = "passed"
    FAILED = "failed"
    SKIPPED = "skipped"


class ErrorCategory(Enum):
    """Categorized error types matching the design document."""
    SUCCESS = "SUCCESS"
    CONFIGURATION_ERROR = "CONFIGURATION_ERROR"
    SUBMODULE_ERROR = "SUBMODULE_ERROR"
    BUILD_ERROR = "BUILD_ERROR"
    PARTITION_ERROR = "PARTITION_ERROR"
    RUN_ERROR = "RUN_ERROR"
    METRICS_ERROR = "METRICS_ERROR"
    REFERENCE_MISMATCH = "REFERENCE_MISMATCH"
    SLURM_ERROR = "SLURM_ERROR"
    TIMEOUT = "TIMEOUT"
    INFRA_ERROR = "INFRA_ERROR"


class ComparisonResult(Enum):
    """Outcome of a metric comparison."""
    PASS = "pass"
    FAIL = "fail"
    SKIP = "skip"


# ---------------------------------------------------------------------------
# Error / stage result
# ---------------------------------------------------------------------------

@dataclass
class TestError:
    """Structured error information attached to a stage failure."""
    category: ErrorCategory
    stage: str
    message: str
    exit_code: Optional[int] = None
    log_path: Optional[str] = None
    first_error_line: Optional[str] = None
    timestamp: Optional[str] = None


@dataclass
class StageResult:
    """Outcome of executing a single pipeline stage."""
    stage: str
    status: StageStatus
    duration_s: Optional[float] = None
    log_path: Optional[str] = None
    error: Optional[TestError] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Metric extraction
# ---------------------------------------------------------------------------

@dataclass
class MetricValue:
    """A single named numeric value extracted from output."""
    name: str
    value: float
    raw: str = ""


@dataclass
class MetricExtractionResult:
    """Result of extracting metrics from a single metric definition."""
    metric_id: str
    values: List[MetricValue] = field(default_factory=list)
    source_file: Optional[str] = None
    source_line: Optional[str] = None
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

@dataclass
class ComparisonDetail:
    """Per-value comparison outcome."""
    metric_id: str
    value_name: str
    measured: float
    reference: float
    tolerance: float
    difference: float
    result: ComparisonResult


@dataclass
class ComparisonReport:
    """Aggregated comparison results for a test run."""
    overall: ComparisonResult
    details: List[ComparisonDetail] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Configuration sub-dataclasses (parsed from YAML)
# ---------------------------------------------------------------------------

@dataclass
class ModuleConfig:
    """Environment modules to load."""
    modules: List[str] = field(default_factory=list)


@dataclass
class GitConfig:
    """Git/submodule configuration."""
    submodule_mode: str = "strict"
    update_recursive: bool = True


@dataclass
class SetupConfig:
    """Setup stage configuration."""
    modules: List[str] = field(default_factory=list)
    git: Optional[GitConfig] = None


@dataclass
class BuildStep:
    """A single build step (cmake_configure or cmake_build)."""
    kind: str
    build_dir: str = "build-release"
    source_dir: str = "."
    options: List[str] = field(default_factory=list)
    targets: List[str] = field(default_factory=list)
    jobs: int = 8


@dataclass
class BuildConfig:
    """Build stage configuration."""
    workspace: str = "build-release"
    cache_policy: str = "reuse_if_compatible"
    steps: List[BuildStep] = field(default_factory=list)


@dataclass
class LevelConfig:
    """Mesh level iteration configuration."""
    start: int = 2
    count: int = 1
    parameter_patch: Optional[Dict[str, str]] = None


@dataclass
class PartitionConfig:
    """Mesh partitioning configuration."""
    tool: str = "python"
    command: str = ""
    expected_partition_count: int = 0


@dataclass
class LaunchConfig:
    """Simulation launch configuration."""
    workdir: str = ""
    mpi_ranks: int = 4
    command: str = ""


@dataclass
class SlurmConfig:
    """SLURM job submission parameters."""
    partition: str = "med"
    nodes: int = 1
    ntasks: int = 4
    cpus_per_task: int = 1
    time: str = "00:30:00"
    constraint: str = ""
    mem_per_cpu: str = ""


@dataclass
class RunConfig:
    """Run stage configuration."""
    levels: Optional[LevelConfig] = None
    partition: Optional[PartitionConfig] = None
    launch: Optional[LaunchConfig] = None
    slurm: Optional[SlurmConfig] = None


@dataclass
class MetricCompareConfig:
    """Comparison settings for a metric."""
    type: str = "tolerance"
    tolerance: Dict[str, float] = field(default_factory=dict)


@dataclass
class MetricConfig:
    """A single metric extraction definition."""
    id: str = ""
    parser: str = "keyword_columns"
    file: str = ""
    keyword: str = ""
    columns: Dict[str, int] = field(default_factory=dict)
    occurrence: str = "last"
    numeric_format: str = "fortran_d_or_e"
    compare: Optional[MetricCompareConfig] = None
    filter: Optional[Dict[str, Any]] = None


@dataclass
class ReferencesConfig:
    """Reference/baseline file pointers."""
    baseline: str = ""


@dataclass
class OutputsConfig:
    """Output artifact configuration."""
    logs: List[str] = field(default_factory=list)
    artifacts: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Top-level test definition
# ---------------------------------------------------------------------------

@dataclass
class TestDefinition:
    """Complete parsed YAML test definition."""
    schema_version: str
    id: str
    name: str
    setup: SetupConfig
    build: BuildConfig
    run: RunConfig
    metrics: List[MetricConfig]
    references: ReferencesConfig
    suite: str = "smoke"
    priority: str = "normal"
    enabled: bool = True
    outputs: Optional[OutputsConfig] = None


# ---------------------------------------------------------------------------
# Run metadata
# ---------------------------------------------------------------------------

@dataclass
class RunMetadata:
    """Metadata captured during a test run."""
    run_id: str
    test_id: str
    commit_sha: str = ""
    branch: str = ""
    dirty: bool = False
    submodule_shas: Dict[str, str] = field(default_factory=dict)
    loaded_modules: List[str] = field(default_factory=list)
    cmake_options: Dict[str, str] = field(default_factory=dict)
    slurm_job_ids: List[str] = field(default_factory=list)
    stages: List[StageResult] = field(default_factory=list)
    metrics: List[MetricExtractionResult] = field(default_factory=list)
    comparison: Optional[ComparisonReport] = None
    errors: List[TestError] = field(default_factory=list)
