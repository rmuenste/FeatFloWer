"""Fortran numeric format parsing.

FeatFloWer logs use Fortran-style exponents where 'D' or 'd' denotes
powers of ten (e.g. 0.7350D-09 = 7.35e-10).  This module normalises
those tokens so Python's float() can handle them.
"""

import re
from typing import List

# Matches signed/unsigned floats with optional D/d/E/e exponent.
_FORTRAN_FLOAT_RE = re.compile(
    r"[+-]?"
    r"(?:\d+\.?\d*|\.\d+)"
    r"(?:[DdEe][+-]?\d+)?"
)


def parse_fortran_float(s: str) -> float:
    """Parse a single Fortran-style float string to Python float.

    Replaces 'D'/'d' exponent markers with 'E' before conversion.
    """
    normalised = s.replace("D", "E").replace("d", "e")
    return float(normalised)


def extract_fortran_floats(line: str) -> List[float]:
    """Extract all Fortran-style floats from a text line."""
    tokens = _FORTRAN_FLOAT_RE.findall(line)
    results = []
    for tok in tokens:
        try:
            results.append(parse_fortran_float(tok))
        except ValueError:
            continue
    return results
