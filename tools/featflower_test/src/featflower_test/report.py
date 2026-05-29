"""Text and JSON report generation."""

import json
from typing import Optional

from featflower_test.results import ResultStore, _enum_serializer


def generate_report(store: ResultStore, fmt: str = "text") -> str:
    """Generate a report from stored results.

    Parameters
    ----------
    store : ResultStore
        Initialised result store pointing to a completed run.
    fmt : str
        Output format: "text" or "json".

    Returns
    -------
    str
        The formatted report.
    """
    metadata = store.load_metadata()
    metrics = store.load_metrics()
    comparison = store.load_comparison()

    if fmt == "json":
        return _json_report(metadata, metrics, comparison)
    return _text_report(metadata, metrics, comparison)


def _json_report(metadata, metrics, comparison) -> str:
    report = {
        "metadata": metadata,
        "metrics": metrics,
        "comparison": comparison,
    }
    return json.dumps(report, indent=2, default=_enum_serializer)


def _text_report(metadata, metrics, comparison) -> str:
    lines = []

    if metadata:
        lines.append(f"Run:    {metadata.get('run_id', '?')}")
        lines.append(f"Test:   {metadata.get('test_id', '?')}")
        lines.append(f"Commit: {metadata.get('commit_sha', '?')}")
        lines.append(f"Branch: {metadata.get('branch', '?')}")
        lines.append("")

        stages = metadata.get("stages", [])
        if stages:
            lines.append("Stages:")
            for s in stages:
                dur = s.get("duration_s")
                dur_str = f" ({dur:.1f}s)" if dur is not None else ""
                lines.append(f"  {s['stage']}: {s['status']}{dur_str}")
            lines.append("")

        errors = metadata.get("errors", [])
        if errors:
            lines.append("Errors:")
            for e in errors:
                lines.append(f"  [{e.get('category', '?')}] {e.get('stage', '?')}: {e.get('message', '?')}")
            lines.append("")

    if metrics:
        lines.append("Metrics:")
        for m in metrics:
            mid = m.get("metric_id", "?")
            if m.get("error"):
                lines.append(f"  {mid}: ERROR - {m['error']}")
            else:
                vals = m.get("values", [])
                parts = [f"{v['name']}={v['value']}" for v in vals]
                lines.append(f"  {mid}: {', '.join(parts)}")
        lines.append("")

    if comparison:
        overall = comparison.get("overall", "?")
        lines.append(f"Comparison: {overall}")
        for d in comparison.get("details", []):
            status = d.get("result", "?")
            symbol = "PASS" if status == "pass" else "FAIL"
            lines.append(
                f"  [{symbol}] {d.get('metric_id', '?')}.{d.get('value_name', '?')}: "
                f"measured={d.get('measured', '?')} ref={d.get('reference', '?')} "
                f"tol={d.get('tolerance', '?')} diff={d.get('difference', '?')}"
            )
        lines.append("")

    return "\n".join(lines)
