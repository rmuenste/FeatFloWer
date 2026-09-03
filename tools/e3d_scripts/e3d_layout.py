"""
Run layout for the e3d launchers: installation directory versus case directory.

The FeatFloWer solvers address every file relative to their working
directory with fixed prefixes (``_data/``, ``_mesh/``, ``_vtk/`` ...).  The
launchers therefore have to run the solvers with the working directory set to
the case folder.  Historically the case folder *was* the installation folder,
so every simulation carried its own copy of the binaries.

``RunLayout`` separates the two locations:

``install_dir``
    Where the executables, the ``partitioner`` package, the default parameter
    templates (``_data_BU``), the default YAML plans and ``_data/MG.dat`` live.
    Defaults to the directory containing the launcher script; can be
    overridden with the environment variable ``FF_GENDIE_HOME``.

``case_dir``
    Where the simulation inputs and outputs live.  Defaults to the current
    working directory, which keeps the classic "run inside the install
    folder" workflow working unchanged.

Input files (parameter templates, YAML plans) are resolved case-first: a file
with the same relative name in the case folder overrides the shipped copy.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

INSTALL_DIR_ENV = "FF_GENDIE_HOME"

# Output skeleton the SSE/gendie solvers expect to exist in the case folder.
CASE_RUNTIME_DIRS = (
    "_data",
    "_mesh",
    "_vtk",
    "_dump",
    "_1D",
    "_hist",
    "_prot0",
    "_prot1",
    "_RTD",
    "start",
)

# Static files every case needs; seeded from the install when missing.
#  _data/MG.dat              lookup table for CreateDumpStructures (umbrella/dump)
#  start/sampleRigidBody.xml rigid-body world read by init_fc_rigid_body
#  start/data.TXT            FullC0ntact parameters read by the same init
CASE_SEED_FILES = (
    Path("_data") / "MG.dat",
    Path("start") / "sampleRigidBody.xml",
    Path("start") / "data.TXT",
)


class RunLayout:
    """Knows where the installation and the case folder are."""

    def __init__(self, install_dir=None, case_dir=None, log=print):
        if install_dir is None:
            install_dir = os.environ.get(INSTALL_DIR_ENV)
        if install_dir is None:
            install_dir = Path(sys.argv[0]).resolve().parent
        self.install_dir = Path(install_dir).resolve()
        self.case_dir = Path(case_dir if case_dir is not None else Path.cwd()).resolve()
        self._log = log
        self._reported = set()

    # ------------------------------------------------------------------ #
    # Queries
    # ------------------------------------------------------------------ #
    def exe(self, name):
        """Absolute path of an executable shipped in the installation."""
        exe_name = name + ".exe" if sys.platform == "win32" else name
        path = self.install_dir / exe_name
        if not path.exists():
            raise FileNotFoundError(
                f"Executable '{exe_name}' not found in installation directory "
                f"{self.install_dir}"
            )
        return str(path)

    def resolve_input(self, relative, must_exist=True):
        """
        Resolve an input file or directory.

        Absolute paths are returned as they are.  Relative paths are looked
        up in the case folder first and then in the installation.  The
        location used is reported once per file so that a stale case-local
        override never goes unnoticed.
        """
        candidate = Path(relative)
        if candidate.is_absolute():
            if must_exist and not candidate.exists():
                raise FileNotFoundError(f"Input '{relative}' not found")
            return candidate

        case_copy = self.case_dir / candidate
        install_copy = self.install_dir / candidate
        if case_copy.exists():
            self._report(relative, "case", case_copy)
            return case_copy
        if install_copy.exists():
            self._report(relative, "install", install_copy)
            return install_copy
        if must_exist:
            raise FileNotFoundError(
                f"Input '{relative}' not found. Looked in case folder "
                f"({case_copy}) and installation ({install_copy})."
            )
        return case_copy

    def _report(self, relative, origin, path):
        key = str(relative)
        if key in self._reported:
            return
        self._reported.add(key)
        if self.case_dir == self.install_dir:
            return
        self._log(f"[layout] {relative}: using {origin} copy ({path})")

    # ------------------------------------------------------------------ #
    # Case folder handling
    # ------------------------------------------------------------------ #
    def prepare_case(self):
        """Create the output skeleton and seed static files. Idempotent."""
        self.case_dir.mkdir(parents=True, exist_ok=True)
        for name in CASE_RUNTIME_DIRS:
            (self.case_dir / name).mkdir(exist_ok=True)
        for rel in CASE_SEED_FILES:
            target = self.case_dir / rel
            if target.exists():
                continue
            source = self.install_dir / rel
            if not source.exists():
                raise FileNotFoundError(
                    f"Cannot seed '{rel}': not found in installation {self.install_dir}"
                )
            shutil.copyfile(str(source), str(target))
            self._log(f"[layout] seeded {rel} from installation")

    def enter_case(self):
        """Change into the case folder. Call once, after paths given on the
        command line have been made absolute."""
        os.chdir(str(self.case_dir))

    def describe(self):
        return (
            f"[layout] installation: {self.install_dir}\n"
            f"[layout] case folder:  {self.case_dir}"
        )


def absolute_from_invocation(path_str, invocation_dir):
    """
    Turn a path given on the command line into an absolute path relative to
    the directory the launcher was invoked from.  Empty strings stay empty.
    """
    if not path_str:
        return path_str
    path = Path(path_str)
    if path.is_absolute():
        return str(path)
    return str((Path(invocation_dir) / path).resolve())
