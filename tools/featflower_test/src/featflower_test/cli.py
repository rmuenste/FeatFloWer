"""CLI entry point for featflower-test.

Commands:
  validate <yaml>                      Load + validate, exit 0/1
  run <yaml> [--slurm] [--dry-run]     Full pipeline
  status <run-id>                      Print stage statuses
  compare <run-id>                     Re-run comparison from stored metrics
  report <run-id> --format json|text   Generate report
"""

import argparse
import os
import subprocess
import sys
import time
import uuid
from datetime import datetime

from featflower_test import __version__
from featflower_test.comparison import compare_metrics
from featflower_test.config import load_test_definition
from featflower_test.errors import FeatFlowerTestError
from featflower_test.models import (
    ComparisonResult,
    MetricExtractionResult,
    RunMetadata,
    StageResult,
    StageStatus,
)
from featflower_test.parsers.keyword_columns import extract_keyword_columns
from featflower_test.report import generate_report
from featflower_test.results import ResultStore
from featflower_test.stages.build import run_build
from featflower_test.stages.partition import run_partition
from featflower_test.stages.run import run_simulation
from featflower_test.stages.setup import run_setup


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="featflower-test",
        description="FeatFloWer automated test runner",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command")

    # validate
    p_val = sub.add_parser("validate", help="Validate a YAML test definition")
    p_val.add_argument("yaml", help="Path to YAML test definition")
    p_val.add_argument("--repo-root", default=None, help="Repository root (auto-detected)")

    # run
    p_run = sub.add_parser("run", help="Run a test through the full pipeline")
    p_run.add_argument("yaml", help="Path to YAML test definition")
    p_run.add_argument("--slurm", action="store_true", help="Use SLURM runner")
    p_run.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    p_run.add_argument("--skip-build", action="store_true", help="Skip the build stage")
    p_run.add_argument("--repo-root", default=None, help="Repository root (auto-detected)")
    p_run.add_argument("--results-dir", default="results", help="Results base directory")

    # status
    p_status = sub.add_parser("status", help="Show status of a run")
    p_status.add_argument("run_id", help="Run ID")
    p_status.add_argument("--test-id", default=None, help="Test ID (auto-detected)")
    p_status.add_argument("--results-dir", default="results", help="Results base directory")

    # compare
    p_cmp = sub.add_parser("compare", help="Re-run comparison from stored metrics")
    p_cmp.add_argument("run_id", help="Run ID")
    p_cmp.add_argument("--test-id", default=None, help="Test ID")
    p_cmp.add_argument("--results-dir", default="results", help="Results base directory")

    # report
    p_rep = sub.add_parser("report", help="Generate a report")
    p_rep.add_argument("run_id", help="Run ID")
    p_rep.add_argument("--test-id", default=None, help="Test ID")
    p_rep.add_argument("--format", choices=["text", "json"], default="text")
    p_rep.add_argument("--results-dir", default="results", help="Results base directory")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "validate":
        return cmd_validate(args)
    elif args.command == "run":
        return cmd_run(args)
    elif args.command == "status":
        return cmd_status(args)
    elif args.command == "compare":
        return cmd_compare(args)
    elif args.command == "report":
        return cmd_report(args)

    return 0


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_validate(args) -> int:
    """Load and validate a YAML test definition."""
    repo_root = _resolve_repo_root(args.repo_root)
    try:
        td = load_test_definition(args.yaml, repo_root)
        print(f"OK: {td.id} ({td.name})")
        print(f"  schema_version: {td.schema_version}")
        print(f"  suite: {td.suite}")
        print(f"  build steps: {len(td.build.steps)}")
        print(f"  metrics: {len(td.metrics)}")
        print(f"  baseline: {td.references.baseline}")
        return 0
    except FeatFlowerTestError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1


def cmd_run(args) -> int:
    """Full pipeline: validate -> setup -> build -> partition -> run -> extract -> compare."""
    repo_root = _resolve_repo_root(args.repo_root)
    results_dir = os.path.join(repo_root, args.results_dir)

    # Validate
    try:
        td = load_test_definition(args.yaml, repo_root)
    except FeatFlowerTestError as exc:
        print(f"Validation failed: {exc}", file=sys.stderr)
        return 1

    # Create run
    run_id = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:8]
    store = ResultStore(results_dir, run_id, td.id)
    store.ensure_dirs()

    metadata = RunMetadata(run_id=run_id, test_id=td.id)
    errors = []
    all_stages = []

    print(f"Run ID: {run_id}")
    print(f"Test:   {td.id} ({td.name})")
    print()

    # -- Setup --
    print("Stage: setup ... ", end="", flush=True)
    setup_result = run_setup(td.setup, repo_root, store.stage_log_path("setup"))
    all_stages.append(setup_result)
    print(setup_result.status.value)
    if setup_result.status == StageStatus.FAILED:
        errors.append(setup_result.error)
        metadata.stages = [_stage_to_dict(s) for s in all_stages]
        metadata.errors = errors
        store.save_metadata(metadata)
        store.save_errors(errors)
        return 1

    metadata.submodule_shas = setup_result.metadata.get("submodule_shas", {})

    # -- Build --
    if args.skip_build:
        print("Stage: build ... skipped")
        all_stages.append(StageResult(stage="build", status=StageStatus.SKIPPED))
    else:
        print("Stage: build ... ", end="", flush=True)
        build_result = run_build(
            td.build, repo_root, store.stage_log_path("build"), dry_run=args.dry_run,
        )
        all_stages.append(build_result)
        print(build_result.status.value)
        if build_result.status == StageStatus.FAILED:
            errors.append(build_result.error)
            metadata.stages = [_stage_to_dict(s) for s in all_stages]
            metadata.errors = errors
            store.save_metadata(metadata)
            store.save_errors(errors)
            return 1
        metadata.cmake_options = build_result.metadata.get("cmake_options", {})

    # -- Partition --
    if td.run.partition:
        print("Stage: partition ... ", end="", flush=True)
        launch = td.run.launch
        workdir = os.path.join(repo_root, launch.workdir) if launch and launch.workdir else repo_root
        part_result = run_partition(
            td.run.partition, repo_root, workdir,
            store.stage_log_path("partition"),
            dry_run=args.dry_run,
        )
        all_stages.append(part_result)
        print(part_result.status.value)
        if part_result.status == StageStatus.FAILED:
            errors.append(part_result.error)
            metadata.stages = [_stage_to_dict(s) for s in all_stages]
            metadata.errors = errors
            store.save_metadata(metadata)
            store.save_errors(errors)
            return 1

    # -- Run (per level) --
    runner = None
    if args.slurm:
        from featflower_test.runners.slurm import SlurmRunner
        runner = SlurmRunner()

    print("Stage: run ... ", end="", flush=True)
    run_results = run_simulation(
        td.run, repo_root, store.stages_dir,
        runner=runner, dry_run=args.dry_run,
    )
    all_stages.extend(run_results)

    run_failed = any(r.status == StageStatus.FAILED for r in run_results)
    print("failed" if run_failed else "passed")

    if run_failed:
        for r in run_results:
            if r.error:
                errors.append(r.error)
        metadata.stages = [_stage_to_dict(s) for s in all_stages]
        metadata.errors = errors
        store.save_metadata(metadata)
        store.save_errors(errors)
        return 1

    # Collect SLURM job IDs from stage metadata
    for r in run_results:
        jid = r.metadata.get("slurm_job_id")
        if jid:
            metadata.slurm_job_ids.append(jid)

    # -- Extract metrics --
    print("Stage: metrics ... ", end="", flush=True)
    extraction_results = []
    launch = td.run.launch
    workdir = os.path.join(repo_root, launch.workdir) if launch and launch.workdir else repo_root

    for metric in td.metrics:
        file_path = os.path.join(workdir, metric.file)
        if metric.parser == "keyword_columns":
            result = extract_keyword_columns(file_path, metric)
        else:
            result = MetricExtractionResult(
                metric_id=metric.id,
                error=f"Unsupported parser: {metric.parser}",
            )
        extraction_results.append(result)

    metric_errors = [r for r in extraction_results if r.error]
    print(f"{'failed' if metric_errors else 'passed'} ({len(extraction_results)} metrics)")
    metadata.metrics = extraction_results
    store.save_metrics(extraction_results)

    # -- Compare --
    baseline_path = os.path.join(repo_root, td.references.baseline) if td.references.baseline else ""

    if baseline_path and os.path.isfile(baseline_path):
        print("Stage: compare ... ", end="", flush=True)
        comparison = compare_metrics(extraction_results, td.metrics, baseline_path)
        metadata.comparison = comparison
        store.save_comparison(comparison)
        print(comparison.overall.value)
    else:
        print("Stage: compare ... skipped (no baseline)")

    # -- Record git info --
    metadata.commit_sha = _git_rev_parse(repo_root)
    metadata.branch = _git_branch(repo_root)
    metadata.dirty = _git_dirty(repo_root)

    # -- Save metadata --
    metadata.stages = [_stage_to_dict(s) for s in all_stages]
    metadata.errors = errors
    store.save_metadata(metadata)
    if errors:
        store.save_errors(errors)

    # -- Print report --
    print()
    print(generate_report(store, fmt="text"))

    # Exit code based on comparison
    if metadata.comparison and metadata.comparison.overall == ComparisonResult.FAIL:
        return 1
    return 0


def cmd_status(args) -> int:
    """Print stage statuses from a stored run."""
    store = _find_store(args)
    if store is None:
        return 1
    metadata = store.load_metadata()
    if metadata is None:
        print(f"No metadata found for run {args.run_id}", file=sys.stderr)
        return 1

    print(f"Run:  {metadata.get('run_id', '?')}")
    print(f"Test: {metadata.get('test_id', '?')}")
    for s in metadata.get("stages", []):
        print(f"  {s['stage']}: {s['status']}")
    return 0


def cmd_compare(args) -> int:
    """Re-run comparison from stored metrics (not implemented in Phase 1a)."""
    print("compare command: re-run comparison from stored metrics")
    print("(requires stored metrics and test definition)")
    return 0


def cmd_report(args) -> int:
    """Generate and print a report."""
    store = _find_store(args)
    if store is None:
        return 1
    print(generate_report(store, fmt=args.format))
    return 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_repo_root(explicit: str = None) -> str:
    """Resolve repository root from explicit arg or git."""
    if explicit:
        return os.path.abspath(explicit)
    try:
        output = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return output.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return os.getcwd()


def _git_rev_parse(repo_root: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root, stderr=subprocess.DEVNULL, text=True,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def _git_branch(repo_root: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=repo_root, stderr=subprocess.DEVNULL, text=True,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def _git_dirty(repo_root: str) -> bool:
    try:
        result = subprocess.run(
            ["git", "diff", "--quiet"],
            cwd=repo_root, stderr=subprocess.DEVNULL,
        )
        return result.returncode != 0
    except FileNotFoundError:
        return False


def _stage_to_dict(s: StageResult) -> dict:
    """Convert StageResult to a plain dict for JSON serialization."""
    d = {
        "stage": s.stage,
        "status": s.status.value,
    }
    if s.duration_s is not None:
        d["duration_s"] = s.duration_s
    if s.log_path:
        d["log_path"] = s.log_path
    return d


def _find_store(args) -> ResultStore:
    """Locate a ResultStore from CLI args.

    For status/compare/report, we need to find the test_id if not given.
    """
    results_dir = getattr(args, "results_dir", "results")
    run_dir = os.path.join(results_dir, "runs", args.run_id, "tests")

    test_id = getattr(args, "test_id", None)
    if test_id is None:
        if not os.path.isdir(run_dir):
            print(f"Run directory not found: {run_dir}", file=sys.stderr)
            return None
        entries = os.listdir(run_dir)
        if len(entries) == 0:
            print(f"No tests found in {run_dir}", file=sys.stderr)
            return None
        test_id = entries[0]

    return ResultStore(results_dir, args.run_id, test_id)


if __name__ == "__main__":
    sys.exit(main())
