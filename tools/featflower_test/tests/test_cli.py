"""Tests for the CLI entry point."""

import os

import pytest

from featflower_test.cli import main


class TestCLIValidate:
    def test_validate_valid_yaml(self, guide01_yaml):
        rc = main(["validate", guide01_yaml, "--repo-root",
                    os.path.dirname(os.path.dirname(guide01_yaml))])
        assert rc == 0

    def test_validate_missing_file(self, tmp_path):
        rc = main(["validate", str(tmp_path / "nope.yaml"),
                    "--repo-root", str(tmp_path)])
        assert rc == 1

    def test_validate_invalid_yaml(self, tmp_path):
        bad = tmp_path / "bad.yaml"
        bad.write_text("schema_version: '9.9'\nid: x\nname: x\n"
                       "setup: {}\nbuild: {steps: []}\n"
                       "run: {}\nmetrics: []\nreferences: {}\n")
        rc = main(["validate", str(bad), "--repo-root", str(tmp_path)])
        assert rc == 1


class TestCLINoCommand:
    def test_no_command_returns_1(self):
        rc = main([])
        assert rc == 1


class TestCLIVersion:
    def test_version_flag(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(["--version"])
        assert exc_info.value.code == 0
