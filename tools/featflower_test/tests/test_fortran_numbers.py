"""Tests for Fortran number parsing."""

import pytest

from featflower_test.parsers.fortran_numbers import (
    extract_fortran_floats,
    parse_fortran_float,
)


class TestParseFortranFloat:
    def test_standard_e_notation(self):
        assert parse_fortran_float("1.23E+04") == pytest.approx(1.23e4)

    def test_fortran_d_notation(self):
        assert parse_fortran_float("0.7350D-09") == pytest.approx(7.35e-10)

    def test_lowercase_d(self):
        assert parse_fortran_float("1.5d+02") == pytest.approx(150.0)

    def test_plain_float(self):
        assert parse_fortran_float("5.579535") == pytest.approx(5.579535)

    def test_integer_string(self):
        assert parse_fortran_float("42") == pytest.approx(42.0)

    def test_negative(self):
        assert parse_fortran_float("-3.14") == pytest.approx(-3.14)

    def test_negative_d_exponent(self):
        assert parse_fortran_float("-0.1940D-09") == pytest.approx(-1.94e-10)

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            parse_fortran_float("not_a_number")


class TestExtractFortranFloats:
    def test_bench_force_line(self):
        line = "BenchForce: 0.1000E+01 5.579535 0.010618 0.2100E+01 0.0062E-01"
        result = extract_fortran_floats(line)
        assert len(result) == 5
        assert result[0] == pytest.approx(1.0)
        assert result[1] == pytest.approx(5.579535)
        assert result[2] == pytest.approx(0.010618)

    def test_fortran_d_line(self):
        line = "   0  0.1281E-08   1.000"
        result = extract_fortran_floats(line)
        assert len(result) == 3
        assert result[0] == pytest.approx(0.0)
        assert result[1] == pytest.approx(1.281e-9)

    def test_empty_line(self):
        assert extract_fortran_floats("no numbers here!") == []

    def test_mixed_d_and_e(self):
        line = "0.5D+02 1.0E-03"
        result = extract_fortran_floats(line)
        assert len(result) == 2
        assert result[0] == pytest.approx(50.0)
        assert result[1] == pytest.approx(0.001)
