"""Tests for keyword_columns metric extractor."""

import os

import pytest

from featflower_test.models import MetricConfig, MetricCompareConfig
from featflower_test.parsers.keyword_columns import extract_keyword_columns


def _make_metric(keyword="BenchForce:", columns=None, occurrence="last"):
    if columns is None:
        columns = {"drag": 2, "lift": 3}
    return MetricConfig(
        id="drag_lift",
        parser="keyword_columns",
        file="_data/prot.txt",
        keyword=keyword,
        columns=columns,
        occurrence=occurrence,
    )


class TestExtractKeywordColumns:
    def test_last_occurrence(self, prot_file):
        metric = _make_metric()
        result = extract_keyword_columns(prot_file, metric)
        assert result.error is None
        assert len(result.values) == 2
        # Last BenchForce line: BenchForce: 0.3000E+01 5.579535 0.010618 ...
        drag = next(v for v in result.values if v.name == "drag")
        lift = next(v for v in result.values if v.name == "lift")
        assert drag.value == pytest.approx(5.579535)
        assert lift.value == pytest.approx(0.010618)

    def test_first_occurrence(self, prot_file):
        metric = _make_metric(occurrence="first")
        result = extract_keyword_columns(prot_file, metric)
        assert result.error is None
        drag = next(v for v in result.values if v.name == "drag")
        # First BenchForce line: ... 5.579535 ...
        assert drag.value == pytest.approx(5.579535)

    def test_keyword_not_found(self, prot_file):
        metric = _make_metric(keyword="NONEXISTENT_KEYWORD:")
        result = extract_keyword_columns(prot_file, metric)
        assert result.error is not None
        assert "not found" in result.error

    def test_column_out_of_range(self, prot_file):
        metric = _make_metric(columns={"too_far": 99})
        result = extract_keyword_columns(prot_file, metric)
        assert result.error is not None
        assert "out of range" in result.error

    def test_missing_file(self, tmp_path):
        metric = _make_metric()
        result = extract_keyword_columns(
            str(tmp_path / "nonexistent.txt"), metric,
        )
        assert result.error is not None
        assert "Cannot read" in result.error

    def test_force_acting_keyword(self, tmp_path):
        """Test with 'Force acting' keyword used in design doc example."""
        prot = tmp_path / "prot.txt"
        prot.write_text(
            "some header\n"
            "Force acting on body: 0.100E+01 0.200E+01 0.300E+01 0.400E+01 "
            "0.500E+01 0.600E+01 5.5795D+00 1.0618D-02\n"
        )
        metric = MetricConfig(
            id="force_acting",
            keyword="Force acting",
            columns={"drag": 10, "lift": 11},
            occurrence="last",
        )
        result = extract_keyword_columns(str(prot), metric)
        # Columns are 0-indexed after split
        # "Force acting on body: 0.100E+01 0.200E+01 0.300E+01 0.400E+01
        #  0.500E+01 0.600E+01 5.5795D+00 1.0618D-02"
        # splits to 12 tokens: [0]=Force [1]=acting ... [10]=5.5795D+00 [11]=1.0618D-02
        assert result.error is None
        drag = next(v for v in result.values if v.name == "drag")
        lift = next(v for v in result.values if v.name == "lift")
        assert drag.value == pytest.approx(5.5795)
        assert lift.value == pytest.approx(0.010618)
