"""Setup stage: environment validation and submodule handling."""

import os
import subprocess
import time

from featflower_test.errors import ConfigurationError, SubmoduleError
from featflower_test.models import (
    SetupConfig,
    StageResult,
    StageStatus,
)


def run_setup(
    config: SetupConfig,
    repo_root: str,
    log_path: str,
) -> StageResult:
    """Validate environment and ensure submodules are current.

    Parameters
    ----------
    config : SetupConfig
        Setup section from the test definition.
    repo_root : str
        Repository root directory.
    log_path : str
        Path to write stage log output.

    Returns
    -------
    StageResult
    """
    t0 = time.time()
    metadata = {}

    try:
        _check_env_vars(repo_root)
        submodule_shas = _handle_submodules(config, repo_root, log_path)
        metadata["submodule_shas"] = submodule_shas
    except (ConfigurationError, SubmoduleError) as exc:
        return StageResult(
            stage="setup",
            status=StageStatus.FAILED,
            duration_s=time.time() - t0,
            log_path=log_path,
            error=exc.to_test_error(),
        )

    return StageResult(
        stage="setup",
        status=StageStatus.PASSED,
        duration_s=time.time() - t0,
        log_path=log_path,
        metadata=metadata,
    )


def _check_env_vars(repo_root: str) -> None:
    """Check required environment variables."""
    mesh_dir = os.environ.get("Q2P1_MESH_DIR")
    if not mesh_dir:
        raise ConfigurationError(
            "Q2P1_MESH_DIR environment variable is not set",
            stage="setup",
        )
    if not os.path.isdir(mesh_dir):
        raise ConfigurationError(
            f"Q2P1_MESH_DIR points to non-existent directory: {mesh_dir}",
            stage="setup",
        )


def _handle_submodules(
    config: SetupConfig,
    repo_root: str,
    log_path: str,
) -> dict:
    """Sync and update submodules, return SHA mapping."""
    log_lines = []

    # git submodule sync --recursive
    try:
        result = subprocess.run(
            ["git", "submodule", "sync", "--recursive"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=120,
        )
        log_lines.append(f"=== git submodule sync --recursive (rc={result.returncode}) ===\n")
        log_lines.append(result.stdout)
        if result.stderr:
            log_lines.append(result.stderr)
        if result.returncode != 0:
            raise SubmoduleError(
                f"git submodule sync failed (rc={result.returncode}): {result.stderr}",
                stage="setup",
                exit_code=result.returncode,
            )
    except subprocess.TimeoutExpired:
        raise SubmoduleError(
            "git submodule sync timed out after 120s",
            stage="setup",
        )

    # git submodule update --init --recursive
    try:
        result = subprocess.run(
            ["git", "submodule", "update", "--init", "--recursive"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=300,
        )
        log_lines.append(f"\n=== git submodule update --init --recursive (rc={result.returncode}) ===\n")
        log_lines.append(result.stdout)
        if result.stderr:
            log_lines.append(result.stderr)
        if result.returncode != 0:
            raise SubmoduleError(
                f"git submodule update failed (rc={result.returncode}): {result.stderr}",
                stage="setup",
                exit_code=result.returncode,
            )
    except subprocess.TimeoutExpired:
        raise SubmoduleError(
            "git submodule update timed out after 300s",
            stage="setup",
        )

    # Record submodule SHAs
    try:
        result = subprocess.run(
            ["git", "submodule", "status", "--recursive"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=30,
        )
        log_lines.append(f"\n=== git submodule status --recursive ===\n")
        log_lines.append(result.stdout)
    except subprocess.TimeoutExpired:
        pass

    # Write log
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as fh:
        fh.writelines(log_lines)

    # Parse SHAs from status output
    shas = {}
    if result.returncode == 0:
        for line in result.stdout.strip().splitlines():
            # Format: " <sha> <path> (<desc>)" or "+<sha> <path>"
            parts = line.strip().split()
            if len(parts) >= 2:
                sha = parts[0].lstrip("+-")
                path = parts[1]
                shas[path] = sha

    return shas
