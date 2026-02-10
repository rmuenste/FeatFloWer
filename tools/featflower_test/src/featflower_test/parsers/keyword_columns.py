"""Keyword + column-index metric extractor.

Replaces the legacy ``get_col_data`` and ``getLogEntry`` functions from
``q2p1_ctest_start.py`` (lines 200-266).

Behaviour (matching legacy):
  1. Read file, find all lines containing ``keyword``.
  2. Select by ``occurrence`` ("last" → ``matching_lines[-1]``,
     matching legacy ``reversed()`` search).
  3. ``line.split()`` then extract by 0-indexed column indices
     (matching legacy ``val[int(idx)]``).
  4. Parse values through ``parse_fortran_float``.
  5. Return ``MetricExtractionResult`` with typed ``MetricValue`` objects.
"""

from typing import Dict, Optional

from featflower_test.models import MetricConfig, MetricExtractionResult, MetricValue
from featflower_test.parsers.fortran_numbers import parse_fortran_float


def extract_keyword_columns(
    file_path: str,
    metric: MetricConfig,
) -> MetricExtractionResult:
    """Extract named metric values from *file_path* using keyword + column indices.

    Parameters
    ----------
    file_path : str
        Absolute path to the log/protocol file.
    metric : MetricConfig
        Metric definition with ``keyword``, ``columns`` (name→index),
        and ``occurrence``.

    Returns
    -------
    MetricExtractionResult
    """
    try:
        with open(file_path, "r") as fh:
            lines = fh.readlines()
    except OSError as exc:
        return MetricExtractionResult(
            metric_id=metric.id,
            source_file=file_path,
            error=f"Cannot read file: {exc}",
        )

    matching = [line for line in lines if metric.keyword in line]

    if not matching:
        return MetricExtractionResult(
            metric_id=metric.id,
            source_file=file_path,
            error=f"Keyword '{metric.keyword}' not found in {file_path}",
        )

    # Select occurrence (default: last, matching legacy reversed() search)
    if metric.occurrence == "first":
        selected = matching[0]
    else:
        selected = matching[-1]

    tokens = selected.split()

    values = []
    for name, col_idx in metric.columns.items():
        if col_idx >= len(tokens):
            return MetricExtractionResult(
                metric_id=metric.id,
                source_file=file_path,
                source_line=selected.strip(),
                error=(
                    f"Column index {col_idx} out of range "
                    f"(line has {len(tokens)} tokens)"
                ),
            )
        raw = tokens[col_idx]
        try:
            val = parse_fortran_float(raw)
        except ValueError:
            return MetricExtractionResult(
                metric_id=metric.id,
                source_file=file_path,
                source_line=selected.strip(),
                error=f"Cannot parse '{raw}' as float for metric '{name}'",
            )
        values.append(MetricValue(name=name, value=val, raw=raw))

    return MetricExtractionResult(
        metric_id=metric.id,
        values=values,
        source_file=file_path,
        source_line=selected.strip(),
    )
