"""YAML loading, schema validation, and path resolution."""

import os
import warnings
from typing import Dict, List, Optional

import yaml

from featflower_test.errors import ConfigurationError
from featflower_test.models import (
    BuildConfig,
    BuildStep,
    GitConfig,
    LaunchConfig,
    LevelConfig,
    MetricCompareConfig,
    MetricConfig,
    OutputsConfig,
    PartitionConfig,
    ReferencesConfig,
    RunConfig,
    SetupConfig,
    SlurmConfig,
    TestDefinition,
)

_SUPPORTED_VERSIONS = {"0.2"}
_REQUIRED_KEYS = {
    "schema_version", "id", "name", "setup", "build", "run", "metrics", "references",
}
_VALID_BUILD_KINDS = {"cmake_configure", "cmake_build"}


def load_test_definition(yaml_path: str, repo_root: str) -> TestDefinition:
    """Load and validate a YAML test definition.

    Parameters
    ----------
    yaml_path : str
        Path to the YAML test definition file.
    repo_root : str
        Absolute path to the FeatFloWer repository root.

    Returns
    -------
    TestDefinition

    Raises
    ------
    ConfigurationError
        On schema violations or missing required keys.
    """
    try:
        with open(yaml_path, "r") as fh:
            raw = yaml.safe_load(fh)
    except OSError as exc:
        raise ConfigurationError(
            f"Cannot read test definition: {exc}",
            stage="validate",
        )
    except yaml.YAMLError as exc:
        raise ConfigurationError(
            f"Invalid YAML syntax: {exc}",
            stage="validate",
        )

    if not isinstance(raw, dict):
        raise ConfigurationError(
            "Test definition must be a YAML mapping",
            stage="validate",
        )

    # -- schema_version --
    sv = str(raw.get("schema_version", ""))
    if sv not in _SUPPORTED_VERSIONS:
        raise ConfigurationError(
            f"Unsupported schema_version '{sv}'; expected one of {_SUPPORTED_VERSIONS}",
            stage="validate",
        )

    # -- required keys --
    missing = _REQUIRED_KEYS - set(raw.keys())
    if missing:
        raise ConfigurationError(
            f"Missing required keys: {sorted(missing)}",
            stage="validate",
        )

    # -- parse sections --
    setup = _parse_setup(raw.get("setup", {}))
    build = _parse_build(raw.get("build", {}))
    run = _parse_run(raw.get("run", {}))
    metrics = _parse_metrics(raw.get("metrics", []))
    references = _parse_references(raw.get("references", {}))
    outputs = _parse_outputs(raw.get("outputs")) if "outputs" in raw else None

    # -- warn on missing baseline --
    baseline_abs = os.path.join(repo_root, references.baseline)
    if references.baseline and not os.path.isfile(baseline_abs):
        warnings.warn(
            f"Baseline file not found: {baseline_abs}",
            stacklevel=2,
        )

    return TestDefinition(
        schema_version=sv,
        id=raw["id"],
        name=raw["name"],
        setup=setup,
        build=build,
        run=run,
        metrics=metrics,
        references=references,
        suite=raw.get("suite", "smoke"),
        priority=raw.get("priority", "normal"),
        enabled=raw.get("enabled", True),
        outputs=outputs,
    )


# ---------------------------------------------------------------------------
# Section parsers
# ---------------------------------------------------------------------------

def _parse_setup(data: dict) -> SetupConfig:
    git_data = data.get("git")
    git = None
    if git_data:
        subs = git_data.get("submodules", {})
        git = GitConfig(
            submodule_mode=subs.get("mode", "strict"),
            update_recursive=subs.get("update_recursive", True),
        )
    return SetupConfig(
        modules=data.get("modules", []),
        git=git,
    )


def _parse_build(data: dict) -> BuildConfig:
    steps = []
    for step_data in data.get("steps", []):
        kind = step_data.get("kind", "")
        if kind not in _VALID_BUILD_KINDS:
            raise ConfigurationError(
                f"Invalid build step kind '{kind}'; expected one of {_VALID_BUILD_KINDS}",
                stage="validate",
            )
        steps.append(BuildStep(
            kind=kind,
            build_dir=step_data.get("build_dir", "build-release"),
            source_dir=step_data.get("source_dir", "."),
            options=step_data.get("options", []),
            targets=step_data.get("targets", []),
            jobs=step_data.get("jobs", 8),
        ))
    return BuildConfig(
        workspace=data.get("workspace", "build-release"),
        cache_policy=data.get("cache_policy", "reuse_if_compatible"),
        steps=steps,
    )


def _parse_run(data: dict) -> RunConfig:
    levels = None
    levels_data = data.get("levels")
    if levels_data:
        levels = LevelConfig(
            start=levels_data.get("start", 2),
            count=levels_data.get("count", 1),
            parameter_patch=levels_data.get("parameter_patch"),
        )

    partition = None
    part_data = data.get("partition")
    if part_data:
        partition = PartitionConfig(
            tool=part_data.get("tool", "python"),
            command=part_data.get("command", ""),
            expected_partition_count=part_data.get("expected_partition_count", 0),
        )

    launch = None
    launch_data = data.get("launch")
    if launch_data:
        launch = LaunchConfig(
            workdir=launch_data.get("workdir", ""),
            mpi_ranks=launch_data.get("mpi_ranks", 4),
            command=launch_data.get("command", ""),
        )

    slurm = None
    slurm_data = data.get("slurm")
    if slurm_data:
        slurm = SlurmConfig(
            partition=slurm_data.get("partition", "med"),
            nodes=slurm_data.get("nodes", 1),
            ntasks=slurm_data.get("ntasks", 4),
            cpus_per_task=slurm_data.get("cpus_per_task", 1),
            time=slurm_data.get("time", "00:30:00"),
            constraint=slurm_data.get("constraint", ""),
            mem_per_cpu=slurm_data.get("mem_per_cpu", ""),
        )

    return RunConfig(
        levels=levels,
        partition=partition,
        launch=launch,
        slurm=slurm,
    )


def _parse_metrics(data: list) -> List[MetricConfig]:
    metrics = []
    for item in data:
        compare = None
        cmp_data = item.get("compare")
        if cmp_data:
            compare = MetricCompareConfig(
                type=cmp_data.get("type", "tolerance"),
                tolerance=cmp_data.get("tolerance", {}),
            )
        metrics.append(MetricConfig(
            id=item.get("id", ""),
            parser=item.get("parser", "keyword_columns"),
            file=item.get("file", ""),
            keyword=item.get("keyword", ""),
            columns=item.get("columns", {}),
            occurrence=item.get("occurrence", "last"),
            numeric_format=item.get("numeric_format", "fortran_d_or_e"),
            compare=compare,
            filter=item.get("filter"),
        ))
    return metrics


def _parse_references(data: dict) -> ReferencesConfig:
    return ReferencesConfig(
        baseline=data.get("baseline", ""),
    )


def _parse_outputs(data: Optional[dict]) -> Optional[OutputsConfig]:
    if data is None:
        return None
    return OutputsConfig(
        logs=data.get("logs", []),
        artifacts=data.get("artifacts", []),
    )
