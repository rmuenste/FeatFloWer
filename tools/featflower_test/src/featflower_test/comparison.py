"""Tolerance-based comparison engine.

Loads a baseline YAML file, extracts reference values by metric_id,
and compares each named value:  abs(measured - reference) <= tolerance.
"""

from typing import Dict, List, Optional

import yaml

from featflower_test.models import (
    ComparisonDetail,
    ComparisonReport,
    ComparisonResult,
    MetricConfig,
    MetricExtractionResult,
)


def load_baseline(baseline_path: str) -> Dict:
    """Load a baseline YAML file and return its contents as a dict."""
    with open(baseline_path, "r") as fh:
        return yaml.safe_load(fh)


def compare_metrics(
    extraction_results: List[MetricExtractionResult],
    metric_configs: List[MetricConfig],
    baseline_path: str,
) -> ComparisonReport:
    """Compare extracted metrics against baseline reference values.

    Parameters
    ----------
    extraction_results : list of MetricExtractionResult
        Metrics extracted from the current run.
    metric_configs : list of MetricConfig
        Metric definitions (carry tolerance info).
    baseline_path : str
        Path to baseline YAML file with reference values.

    Returns
    -------
    ComparisonReport
    """
    baseline = load_baseline(baseline_path)
    baseline_metrics = baseline.get("metrics", {})

    details: List[ComparisonDetail] = []
    overall = ComparisonResult.PASS

    # Build lookup: metric_id → MetricConfig
    config_by_id = {m.id: m for m in metric_configs}

    for result in extraction_results:
        if result.error:
            details.append(ComparisonDetail(
                metric_id=result.metric_id,
                value_name="(extraction failed)",
                measured=0.0,
                reference=0.0,
                tolerance=0.0,
                difference=0.0,
                result=ComparisonResult.FAIL,
            ))
            overall = ComparisonResult.FAIL
            continue

        cfg = config_by_id.get(result.metric_id)
        if cfg is None or cfg.compare is None:
            continue

        ref_block = baseline_metrics.get(result.metric_id, {})
        tolerances = cfg.compare.tolerance

        for mv in result.values:
            ref_val = ref_block.get(mv.name)
            tol = tolerances.get(mv.name, 0.0)

            if ref_val is None:
                details.append(ComparisonDetail(
                    metric_id=result.metric_id,
                    value_name=mv.name,
                    measured=mv.value,
                    reference=0.0,
                    tolerance=tol,
                    difference=0.0,
                    result=ComparisonResult.SKIP,
                ))
                continue

            ref_val = float(ref_val)
            diff = abs(mv.value - ref_val)
            passed = diff <= tol
            result_status = ComparisonResult.PASS if passed else ComparisonResult.FAIL

            if result_status == ComparisonResult.FAIL:
                overall = ComparisonResult.FAIL

            details.append(ComparisonDetail(
                metric_id=result.metric_id,
                value_name=mv.name,
                measured=mv.value,
                reference=ref_val,
                tolerance=tol,
                difference=diff,
                result=result_status,
            ))

    return ComparisonReport(overall=overall, details=details)
