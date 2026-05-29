"""Build stage: ordered CMake configure + build steps."""

import os
import subprocess
import time

from featflower_test.errors import BuildError
from featflower_test.models import BuildConfig, StageResult, StageStatus


def run_build(
    config: BuildConfig,
    repo_root: str,
    log_path: str,
    dry_run: bool = False,
) -> StageResult:
    """Execute ordered CMake build steps.

    Parameters
    ----------
    config : BuildConfig
        Build section from the test definition.
    repo_root : str
        Repository root directory.
    log_path : str
        Path to write build log.
    dry_run : bool
        If True, print commands but do not execute.

    Returns
    -------
    StageResult
    """
    t0 = time.time()
    metadata = {"cmake_options": {}, "commands": []}

    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    try:
        with open(log_path, "w") as log_fh:
            for step in config.steps:
                if step.kind == "cmake_configure":
                    _run_cmake_configure(step, repo_root, log_fh, metadata, dry_run)
                elif step.kind == "cmake_build":
                    _run_cmake_build(step, repo_root, log_fh, metadata, dry_run)
    except BuildError as exc:
        return StageResult(
            stage="build",
            status=StageStatus.FAILED,
            duration_s=time.time() - t0,
            log_path=log_path,
            error=exc.to_test_error(),
            metadata=metadata,
        )

    return StageResult(
        stage="build",
        status=StageStatus.PASSED,
        duration_s=time.time() - t0,
        log_path=log_path,
        metadata=metadata,
    )


def _run_cmake_configure(step, repo_root, log_fh, metadata, dry_run):
    """Run a cmake configure step."""
    build_dir = os.path.join(repo_root, step.build_dir)
    source_dir = os.path.join(repo_root, step.source_dir)

    cmd = ["cmake"]
    for opt in step.options:
        cmd.append(f"-D{opt}")
        # Record option in metadata
        if "=" in opt:
            key, val = opt.split("=", 1)
            metadata["cmake_options"][key] = val
    cmd.extend(["-S", source_dir, "-B", build_dir])

    metadata["commands"].append(" ".join(cmd))
    log_fh.write(f"=== cmake configure ===\n$ {' '.join(cmd)}\n\n")
    log_fh.flush()

    if dry_run:
        log_fh.write("[DRY RUN] skipping execution\n\n")
        return

    os.makedirs(build_dir, exist_ok=True)
    result = subprocess.run(
        cmd,
        cwd=repo_root,
        stdout=log_fh,
        stderr=subprocess.STDOUT,
    )
    log_fh.write(f"\n--- exit code: {result.returncode} ---\n\n")
    log_fh.flush()

    if result.returncode != 0:
        raise BuildError(
            f"cmake configure failed (rc={result.returncode})",
            stage="build",
            exit_code=result.returncode,
            log_path=log_fh.name,
        )


def _run_cmake_build(step, repo_root, log_fh, metadata, dry_run):
    """Run a cmake build step for each target."""
    build_dir = os.path.join(repo_root, step.build_dir)

    for target in step.targets:
        cmd = [
            "cmake", "--build", build_dir,
            "--target", target,
            "--", f"-j{step.jobs}",
        ]

        metadata["commands"].append(" ".join(cmd))
        log_fh.write(f"=== cmake build: {target} ===\n$ {' '.join(cmd)}\n\n")
        log_fh.flush()

        if dry_run:
            log_fh.write("[DRY RUN] skipping execution\n\n")
            continue

        result = subprocess.run(
            cmd,
            cwd=repo_root,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
        )
        log_fh.write(f"\n--- exit code: {result.returncode} ---\n\n")
        log_fh.flush()

        if result.returncode != 0:
            raise BuildError(
                f"cmake build target '{target}' failed (rc={result.returncode})",
                stage="build",
                exit_code=result.returncode,
                log_path=log_fh.name,
            )
