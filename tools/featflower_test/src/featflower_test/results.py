"""Result directory layout and JSON persistence.

Layout:
  results/runs/<run-id>/tests/<test-id>/
    stages/          stage log files
    artifacts/       collected output artifacts
    metadata.json    run metadata
    metrics.json     extracted metrics
    compare.json     comparison results
    errors.json      collected errors
"""

import json
import os
from dataclasses import asdict
from typing import List, Optional

from featflower_test.models import (
    ComparisonReport,
    MetricExtractionResult,
    RunMetadata,
    StageResult,
    TestError,
)


def _enum_serializer(obj):
    """JSON serializer that handles Enum values."""
    if hasattr(obj, "value"):
        return obj.value
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


class ResultStore:
    """Manages the result directory layout and persistence."""

    def __init__(self, base_dir: str, run_id: str, test_id: str):
        self.base_dir = base_dir
        self.run_id = run_id
        self.test_id = test_id
        self.test_dir = os.path.join(
            base_dir, "runs", run_id, "tests", test_id,
        )
        self.stages_dir = os.path.join(self.test_dir, "stages")
        self.artifacts_dir = os.path.join(self.test_dir, "artifacts")

    def ensure_dirs(self) -> None:
        """Create the result directory structure."""
        os.makedirs(self.stages_dir, exist_ok=True)
        os.makedirs(self.artifacts_dir, exist_ok=True)

    def stage_log_path(self, stage_name: str) -> str:
        """Return the log path for a given stage."""
        return os.path.join(self.stages_dir, f"{stage_name}.log")

    def save_metadata(self, metadata: RunMetadata) -> str:
        """Save run metadata to JSON."""
        path = os.path.join(self.test_dir, "metadata.json")
        self._write_json(path, asdict(metadata))
        return path

    def save_metrics(self, metrics: List[MetricExtractionResult]) -> str:
        """Save extracted metrics to JSON."""
        path = os.path.join(self.test_dir, "metrics.json")
        self._write_json(path, [asdict(m) for m in metrics])
        return path

    def save_comparison(self, report: ComparisonReport) -> str:
        """Save comparison report to JSON."""
        path = os.path.join(self.test_dir, "compare.json")
        self._write_json(path, asdict(report))
        return path

    def save_errors(self, errors: List[TestError]) -> str:
        """Save error list to JSON."""
        path = os.path.join(self.test_dir, "errors.json")
        self._write_json(path, [asdict(e) for e in errors])
        return path

    def load_metadata(self) -> Optional[dict]:
        """Load metadata.json if it exists."""
        path = os.path.join(self.test_dir, "metadata.json")
        return self._read_json(path)

    def load_metrics(self) -> Optional[list]:
        """Load metrics.json if it exists."""
        path = os.path.join(self.test_dir, "metrics.json")
        return self._read_json(path)

    def load_comparison(self) -> Optional[dict]:
        """Load compare.json if it exists."""
        path = os.path.join(self.test_dir, "compare.json")
        return self._read_json(path)

    def _write_json(self, path: str, data) -> None:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            json.dump(data, fh, indent=2, default=_enum_serializer)
            fh.write("\n")

    def _read_json(self, path: str) -> Optional[dict]:
        if not os.path.isfile(path):
            return None
        with open(path, "r") as fh:
            return json.load(fh)
