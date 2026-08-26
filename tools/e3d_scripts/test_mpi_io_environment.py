#!/usr/bin/env python3

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))

from mpi_io_environment import (  # noqa: E402
    configure_mpi_io_environment,
    detect_openmpi_io_components,
)


class DetectComponentsTests(unittest.TestCase):
    def test_detects_all_io_components(self):
        output = """
                  MCA io: ompio (MCA v2.1.0, API v2.0.0)
                  MCA io: romio321 (MCA v2.1.0, API v2.0.0)
                 """
        status, components = detect_openmpi_io_components(
            which=lambda unused: "/opt/openmpi/bin/ompi_info",
            check_output=lambda *args, **kwargs: output,
        )
        self.assertEqual("available", status)
        self.assertEqual(("ompio", "romio321"), components)

    def test_reports_missing_openmpi(self):
        status, components = detect_openmpi_io_components(
            which=lambda unused: None)
        self.assertEqual("no_openmpi", status)
        self.assertEqual((), components)


class ConfigureEnvironmentTests(unittest.TestCase):
    @staticmethod
    def detector():
        return "available", ("ompio", "romio321")

    def test_selects_romio_and_sets_defaults(self):
        environment = {}
        status = configure_mpi_io_environment(
            environment, emit=lambda unused: None, detector=self.detector)
        self.assertEqual("configured_romio", status)
        self.assertEqual("romio321", environment["OMPI_MCA_io"])
        self.assertEqual("16777216", environment["ROMIO_CB_BUFFER_SIZE"])
        self.assertEqual("enable", environment["ROMIO_DS_WRITE"])

    def test_preserves_explicit_romio_tuning(self):
        environment = {
            "OMPI_MCA_io": "romio321",
            "ROMIO_CB_BUFFER_SIZE": "33554432",
            "ROMIO_DS_WRITE": "disable",
        }
        configure_mpi_io_environment(
            environment, emit=lambda unused: None, detector=self.detector)
        self.assertEqual("33554432", environment["ROMIO_CB_BUFFER_SIZE"])
        self.assertEqual("disable", environment["ROMIO_DS_WRITE"])

    def test_non_romio_selection_does_not_add_romio_tuning(self):
        environment = {"OMPI_MCA_io": "ompio"}
        status = configure_mpi_io_environment(
            environment, emit=lambda unused: None, detector=self.detector)
        self.assertEqual("configured_non_romio", status)
        self.assertNotIn("ROMIO_CB_BUFFER_SIZE", environment)
        self.assertNotIn("ROMIO_DS_WRITE", environment)

    def test_rejects_stale_component_override(self):
        environment = {"OMPI_MCA_io": "romio341"}
        with self.assertRaisesRegex(RuntimeError, "romio341"):
            configure_mpi_io_environment(
                environment, emit=lambda unused: None, detector=self.detector)


if __name__ == "__main__":
    unittest.main()
