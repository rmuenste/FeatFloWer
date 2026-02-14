"""Tests for YAML config loading and validation."""

import os

import pytest

from featflower_test.config import load_test_definition
from featflower_test.errors import ConfigurationError


class TestLoadTestDefinition:
    def test_load_valid_guide01(self, guide01_yaml):
        repo_root = os.path.dirname(os.path.dirname(guide01_yaml))
        td = load_test_definition(guide01_yaml, repo_root)
        assert td.id == "test_guide01"
        assert td.name == "Test Guide 01"
        assert td.schema_version == "0.2"
        assert len(td.build.steps) == 2
        assert td.build.steps[0].kind == "cmake_configure"
        assert td.build.steps[1].kind == "cmake_build"
        assert len(td.metrics) == 1
        assert td.metrics[0].id == "drag_lift"
        assert td.metrics[0].keyword == "BenchForce:"
        assert td.metrics[0].columns == {"drag": 2, "lift": 3}

    def test_load_missing_file(self, tmp_path):
        with pytest.raises(ConfigurationError, match="Cannot read"):
            load_test_definition(
                str(tmp_path / "nonexistent.yaml"),
                str(tmp_path),
            )

    def test_invalid_schema_version(self, tmp_path):
        yaml_file = tmp_path / "bad_version.yaml"
        yaml_file.write_text(
            "schema_version: '9.9'\nid: x\nname: x\n"
            "setup: {}\nbuild: {steps: []}\n"
            "run: {}\nmetrics: []\nreferences: {}\n"
        )
        with pytest.raises(ConfigurationError, match="schema_version"):
            load_test_definition(str(yaml_file), str(tmp_path))

    def test_missing_required_keys(self, tmp_path):
        yaml_file = tmp_path / "missing_keys.yaml"
        yaml_file.write_text(
            "schema_version: '0.2'\nid: x\nname: x\n"
        )
        with pytest.raises(ConfigurationError, match="Missing required"):
            load_test_definition(str(yaml_file), str(tmp_path))

    def test_invalid_build_step_kind(self, tmp_path):
        yaml_file = tmp_path / "bad_step.yaml"
        yaml_file.write_text(
            "schema_version: '0.2'\nid: x\nname: x\n"
            "setup: {}\n"
            "build:\n  steps:\n    - kind: invalid_kind\n"
            "run: {}\nmetrics: []\nreferences: {}\n"
        )
        with pytest.raises(ConfigurationError, match="Invalid build step kind"):
            load_test_definition(str(yaml_file), str(tmp_path))

    def test_invalid_yaml_syntax(self, tmp_path):
        yaml_file = tmp_path / "bad_syntax.yaml"
        yaml_file.write_text("{{invalid yaml::")
        with pytest.raises(ConfigurationError, match="Invalid YAML"):
            load_test_definition(str(yaml_file), str(tmp_path))

    def test_partition_config_parsed(self, guide01_yaml):
        repo_root = os.path.dirname(os.path.dirname(guide01_yaml))
        td = load_test_definition(guide01_yaml, repo_root)
        assert td.run.partition is not None
        assert td.run.partition.expected_partition_count == 3
        assert "PyPartitioner" in td.run.partition.command

    def test_slurm_config_optional(self, guide01_yaml):
        repo_root = os.path.dirname(os.path.dirname(guide01_yaml))
        td = load_test_definition(guide01_yaml, repo_root)
        # guide01.yaml fixture doesn't have slurm section
        assert td.run.slurm is None

    def test_level_config_parsed(self, guide01_yaml):
        repo_root = os.path.dirname(os.path.dirname(guide01_yaml))
        td = load_test_definition(guide01_yaml, repo_root)
        assert td.run.levels is not None
        assert td.run.levels.start == 2
        assert td.run.levels.count == 1
        assert td.run.levels.parameter_patch is not None
