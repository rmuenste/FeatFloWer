"""Shared test fixtures."""

import os

import pytest


@pytest.fixture
def fixtures_dir():
    """Return path to the test fixtures directory."""
    return os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def prot_file(fixtures_dir):
    """Return path to the sample prot_fac2d.txt fixture."""
    return os.path.join(fixtures_dir, "prot_fac2d.txt")


@pytest.fixture
def param_file(fixtures_dir):
    """Return path to the sample parameter file."""
    return os.path.join(fixtures_dir, "q2p1_param_sample.dat")


@pytest.fixture
def guide01_yaml(fixtures_dir):
    """Return path to the guide01 YAML fixture."""
    return os.path.join(fixtures_dir, "guide01.yaml")


@pytest.fixture
def baseline_yaml(fixtures_dir):
    """Return path to the baseline YAML fixture."""
    return os.path.join(fixtures_dir, "baseline_guide01.yaml")
