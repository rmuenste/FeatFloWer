"""Tests for parameter file patching (SimPar@MaxMeshLevel)."""

import os

import pytest

from featflower_test.stages.run import patch_parameter_file


class TestPatchParameterFile:
    def test_replaces_level(self, param_file, tmp_path):
        out_path = str(tmp_path / "q2p1_param.dat")
        patch_parameter_file(param_file, out_path, 5)

        with open(out_path) as fh:
            content = fh.read()

        assert "SimPar@MaxMeshLevel = 5" in content
        # Other lines preserved
        assert 'SimPar@MeshFolder = "NEWFAC"' in content
        assert "SimPar@TimeSteps = 100" in content

    def test_replaces_level_from_different_value(self, param_file, tmp_path):
        out_path = str(tmp_path / "q2p1_param.dat")
        patch_parameter_file(param_file, out_path, 2)

        with open(out_path) as fh:
            content = fh.read()

        assert "SimPar@MaxMeshLevel = 2" in content
        # Original had level 3, should be gone
        assert "SimPar@MaxMeshLevel = 3" not in content

    def test_output_dir_created(self, param_file, tmp_path):
        out_path = str(tmp_path / "subdir" / "_data" / "q2p1_param.dat")
        patch_parameter_file(param_file, out_path, 4)
        assert os.path.isfile(out_path)

    def test_preserves_other_content(self, param_file, tmp_path):
        out_path = str(tmp_path / "q2p1_param.dat")
        patch_parameter_file(param_file, out_path, 7)

        with open(out_path) as fh:
            lines = fh.readlines()

        # Should have same number of lines
        with open(param_file) as fh:
            orig_lines = fh.readlines()

        assert len(lines) == len(orig_lines)
