"""
Custom build hook: copies libmetis.so from $FEATFLOWER_BUILD_DIR into the
package at install time so the library travels with the Python package and no
runtime environment variable is needed.
"""

import os
import shutil
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py

_PKG_DIR = Path("src") / "featflower_partitioner"
_LIB_NAME = "libmetis.so"
_LIB_IN_BUILD = Path("extern") / "libraries" / "metis-5.1.0" / "libmetis" / _LIB_NAME


class BuildPyWithMetis(build_py):
    def run(self):
        build_dir = os.environ.get("FEATFLOWER_BUILD_DIR")
        if build_dir:
            src = Path(build_dir) / _LIB_IN_BUILD
            dst = _PKG_DIR / _LIB_NAME
            if src.exists():
                shutil.copy2(src, dst)
                print(f"featflower-partitioner: bundled {src} -> {dst}")
            else:
                print(
                    f"featflower-partitioner: WARNING — FEATFLOWER_BUILD_DIR is set "
                    f"but {src} was not found.  libmetis.so will not be bundled."
                )
        else:
            print(
                "featflower-partitioner: FEATFLOWER_BUILD_DIR not set — "
                "skipping libmetis.so bundling.  Set the variable and reinstall "
                "to bundle the library, or ensure libmetis.so is on the system path."
            )
        super().run()


setup(cmdclass={"build_py": BuildPyWithMetis})
