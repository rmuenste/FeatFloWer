"""Tests for the tolerance comparison engine."""

import os

import pytest

from featflower_test.comparison import compare_metrics
from featflower_test.models import (
    ComparisonResult,
    MetricCompareConfig,
    MetricConfig,
    MetricExtractionResult,
    MetricValue,
)


def _metric_config(tolerance_drag=1e-4, tolerance_lift=1e-4):
    return MetricConfig(
        id="drag_lift",
        parser="keyword_columns",
        file="_data/prot.txt",
        keyword="BenchForce:",
        columns={"drag": 2, "lift": 3},
        compare=MetricCompareConfig(
            type="tolerance",
            tolerance={"drag": tolerance_drag, "lift": tolerance_lift},
        ),
    )


def _extraction_result(drag=5.579535, lift=0.010618):
    return MetricExtractionResult(
        metric_id="drag_lift",
        values=[
            MetricValue(name="drag", value=drag, raw=str(drag)),
            MetricValue(name="lift", value=lift, raw=str(lift)),
        ],
    )


class TestCompareMetrics:
    def test_pass_within_tolerance(self, baseline_yaml):
        result = compare_metrics(
            [_extraction_result(5.579535, 0.010618)],
            [_metric_config()],
            baseline_yaml,
        )
        assert result.overall == ComparisonResult.PASS
        assert len(result.details) == 2
        for d in result.details:
            assert d.result == ComparisonResult.PASS

    def test_fail_outside_tolerance(self, baseline_yaml):
        result = compare_metrics(
            [_extraction_result(5.580, 0.010618)],
            [_metric_config()],
            baseline_yaml,
        )
        assert result.overall == ComparisonResult.FAIL
        drag_detail = next(d for d in result.details if d.value_name == "drag")
        assert drag_detail.result == ComparisonResult.FAIL
        assert drag_detail.difference == pytest.approx(abs(5.580 - 5.579535))

    def test_pass_at_tolerance_boundary(self, baseline_yaml):
        # Exactly at tolerance boundary
        result = compare_metrics(
            [_extraction_result(5.579535 + 1e-4, 0.010618)],
            [_metric_config()],
            baseline_yaml,
        )
        assert result.overall == ComparisonResult.PASS

    def test_extraction_error_causes_fail(self, baseline_yaml):
        err_result = MetricExtractionResult(
            metric_id="drag_lift",
            error="File not found",
        )
        result = compare_metrics(
            [err_result],
            [_metric_config()],
            baseline_yaml,
        )
        assert result.overall == ComparisonResult.FAIL

    def test_missing_reference_value_skipped(self, tmp_path):
        # Baseline has no lift value
        baseline = tmp_path / "partial.yaml"
        baseline.write_text("metrics:\n  drag_lift:\n    drag: 5.579535\n")

        result = compare_metrics(
            [_extraction_result()],
            [_metric_config()],
            str(baseline),
        )
        # drag passes, lift skipped
        drag_d = next(d for d in result.details if d.value_name == "drag")
        lift_d = next(d for d in result.details if d.value_name == "lift")
        assert drag_d.result == ComparisonResult.PASS
        assert lift_d.result == ComparisonResult.SKIP
