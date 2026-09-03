"""
Tests for the installation/case folder layout of the e3d launchers.

Two levels:

* unit tests of ``e3d_layout.RunLayout``;
* an end-to-end run of ``e3d_start_yaml.py`` against a fake installation
  whose "executables" are tiny shell scripts.  This exercises the real
  launcher code (option parsing, layout setup, folder setup, partitioner,
  YAML plan execution) without a compiled solver and finishes in seconds.

Run with:  python3 -m pytest tools/e3d_scripts/tests -q
"""

import os
import shutil
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPTS_DIR.parents[1]
sys.path.insert(0, str(SCRIPTS_DIR))

from e3d_layout import RunLayout, absolute_from_invocation, CASE_RUNTIME_DIRS  # noqa: E402


# --------------------------------------------------------------------------- #
# Unit tests
# --------------------------------------------------------------------------- #
def make_install(root):
    install = root / "install"
    (install / "_data").mkdir(parents=True)
    (install / "_data" / "MG.dat").write_text("1 1 1\n")
    (install / "start").mkdir()
    (install / "start" / "sampleRigidBody.xml").write_text("<world/>\n")
    (install / "start" / "data.TXT").write_text("params\n")
    (install / "_data_BU").mkdir()
    (install / "_data_BU" / "tpl.dat").write_text("install template\n")
    return install


def test_exe_returns_absolute_install_path(tmp_path):
    install = make_install(tmp_path)
    exe = install / "q2p1_sse"
    exe.write_text("")
    layout = RunLayout(install_dir=install, case_dir=tmp_path / "case")
    assert layout.exe("q2p1_sse") == str(exe)
    with pytest.raises(FileNotFoundError):
        layout.exe("does_not_exist")


def test_prepare_case_creates_skeleton_and_seeds_mg(tmp_path):
    install = make_install(tmp_path)
    case = tmp_path / "case"
    layout = RunLayout(install_dir=install, case_dir=case, log=lambda *_: None)
    layout.prepare_case()
    for name in CASE_RUNTIME_DIRS:
        assert (case / name).is_dir()
    assert (case / "_data" / "MG.dat").read_text() == "1 1 1\n"
    assert (case / "start" / "sampleRigidBody.xml").read_text() == "<world/>\n"
    assert (case / "start" / "data.TXT").read_text() == "params\n"

    # idempotent and never overwrites an existing MG.dat
    (case / "_data" / "MG.dat").write_text("case copy\n")
    layout.prepare_case()
    assert (case / "_data" / "MG.dat").read_text() == "case copy\n"


def test_resolve_input_prefers_case_copy(tmp_path):
    install = make_install(tmp_path)
    case = tmp_path / "case"
    messages = []
    layout = RunLayout(install_dir=install, case_dir=case, log=messages.append)
    layout.prepare_case()

    resolved = layout.resolve_input("_data_BU/tpl.dat")
    assert resolved == install / "_data_BU" / "tpl.dat"
    assert any("install copy" in m for m in messages)

    (case / "_data_BU").mkdir()
    (case / "_data_BU" / "tpl.dat").write_text("case template\n")
    layout = RunLayout(install_dir=install, case_dir=case, log=messages.append)
    resolved = layout.resolve_input(Path("_data_BU") / "tpl.dat")
    assert resolved == case / "_data_BU" / "tpl.dat"
    assert any("case copy" in m for m in messages)

    with pytest.raises(FileNotFoundError):
        layout.resolve_input("_data_BU/missing.dat")
    assert layout.resolve_input("_data_BU/missing.dat", must_exist=False) == \
        case / "_data_BU" / "missing.dat"

    absolute = install / "_data_BU" / "tpl.dat"
    assert layout.resolve_input(absolute) == absolute


def test_absolute_from_invocation(tmp_path):
    assert absolute_from_invocation("", tmp_path) == ""
    assert absolute_from_invocation("/abs/x", tmp_path) == "/abs/x"
    assert absolute_from_invocation("rel/y", tmp_path) == str((tmp_path / "rel" / "y").resolve())


# --------------------------------------------------------------------------- #
# End-to-end run of the launcher against a fake installation
# --------------------------------------------------------------------------- #
FAKE_MPIRUN = textwrap.dedent(
    """\
    #!/bin/sh
    # fake mpirun: drop "-np N" and run the program in the current directory
    while [ $# -gt 0 ]; do
      case "$1" in
        -np|-n) shift 2 ;;
        *) break ;;
      esac
    done
    exec "$@"
    """
)

# The fake solver records where it ran, which parameter file it saw and
# writes into the output directories exactly like the real solver would.
FAKE_SOLVER = textwrap.dedent(
    """\
    #!/bin/sh
    echo "$(pwd) $*" >> solver_calls.txt
    test -f _data/q2p1_param.dat || { echo "no param file" >&2; exit 3; }
    test -f _data/Extrud3D.dat || { echo "no Extrud3D.dat" >&2; exit 3; }
    test -f _data/MG.dat || { echo "no MG.dat" >&2; exit 3; }
    test -f start/sampleRigidBody.xml || { echo "no rigid body file" >&2; exit 3; }
    test -d _mesh/NEWFAC || { echo "no partitioned mesh" >&2; exit 3; }
    head -1 _data/q2p1_param.dat >> solver_calls.txt
    echo "proto" > _data/prot.txt
    echo "0 0" > _1D/extrud3d_0000.res
    exit 0
    """
)

FAKE_MESHER = textwrap.dedent(
    """\
    #!/bin/sh
    # real s3d_mesher would create _data/meshDir; leave it to the launcher's
    # fallback that copies meshDir from the project folder
    echo "$(pwd)" >> mesher_calls.txt
    exit 0
    """
)


def _write_exe(path, content):
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def fake_install(tmp_path):
    """A minimal bin/q2p1_gendie look-alike with stub executables."""
    install = tmp_path / "install"
    install.mkdir()
    for name in ("e3d_start_yaml.py", "e3d_layout.py", "e3d_die.yaml"):
        shutil.copy(SCRIPTS_DIR / name, install / name)
    shutil.copytree(REPO_ROOT / "tools" / "partitioner", install / "partitioner")
    shutil.copytree(REPO_ROOT / "applications" / "q2p1_sse" / "_data_BU", install / "_data_BU")
    (install / "_data").mkdir()
    shutil.copy(REPO_ROOT / "_data" / "MG.dat", install / "_data" / "MG.dat")
    shutil.copytree(REPO_ROOT / "applications" / "q2p1_sse" / "start", install / "start")
    _write_exe(install / "q2p1_sse", FAKE_SOLVER)
    _write_exe(install / "s3d_mesher", FAKE_MESHER)
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    _write_exe(bin_dir / "mpirun", FAKE_MPIRUN)
    return install, bin_dir


def _metis_library():
    """The partitioner needs libmetis.so; take it from any configured build."""
    for build in sorted(REPO_ROOT.glob("build*")):
        candidate = build / "extern" / "libraries" / "metis-5.1.0" / "libmetis" / "libmetis.so"
        if candidate.exists():
            return candidate
    return None


def _run_launcher(install, bin_dir, cwd, args, env_extra=None):
    env = dict(os.environ)
    env["PATH"] = str(bin_dir) + os.pathsep + env["PATH"]
    env["PYTHONUNBUFFERED"] = "1"
    if env_extra:
        env.update(env_extra)
    cmd = [sys.executable, str(install / "e3d_start_yaml.py")] + args
    return subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=True, text=True, timeout=600)


@pytest.mark.skipif(_metis_library() is None, reason="needs a built libmetis.so for the partitioner")
def test_launcher_runs_case_outside_install(fake_install, tmp_path):
    install, bin_dir = fake_install
    metis = _metis_library()
    shutil.copy(metis, install / "libmetis.so")

    # case folder holds nothing but the e3d input folder
    case = tmp_path / "cases" / "die_01"
    case.mkdir(parents=True)
    shutil.copytree(REPO_ROOT / "applications" / "q2p1_sse" / "_ianus" / "DIE" / "RH_LOC",
                    case / "e3d_input")
    # a case-local override of one template must win over the shipped copy
    (case / "_data_BU").mkdir()
    override = case / "_data_BU" / "q2p1_paramV_DIE_1.dat"
    shutil.copy(install / "_data_BU" / "q2p1_paramV_DIE_1.dat", override)
    override.write_text("! CASE-LOCAL OVERRIDE\n" + override.read_text())

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    result = _run_launcher(
        install, bin_dir, elsewhere,
        ["-C", str(case), "-f", str(case / "e3d_input"), "-y", "e3d_die.yaml", "-n", "4",
         "-a", "0"],
    )
    assert result.returncode == 0, result.stdout + result.stderr

    # nothing was written into the installation or the invocation directory
    assert not (install / "solver_calls.txt").exists()
    assert not (install / "e3d.log").exists()
    assert not (elsewhere / "e3d.log").exists()
    assert list(elsewhere.iterdir()) == []

    # everything landed in the case folder
    assert (case / "e3d.log").exists()
    assert (case / "_data" / "MG.dat").exists()
    assert (case / "_data" / "q2p1_param.dat").exists()
    assert (case / "_mesh" / "NEWFAC").is_dir()
    assert (case / "_1D" / "extrud3d_0000.res").exists()
    mesher_calls = (case / "mesher_calls.txt").read_text().strip().splitlines()
    assert mesher_calls == [str(case)]

    calls = (case / "solver_calls.txt").read_text()
    # every solver invocation ran inside the case folder with the install binary
    for line in calls.splitlines():
        if line.startswith("/"):
            assert line.startswith(str(case) + " "), line
    # plan: init uses DIE_0 (shipped), main and final use DIE_1 (case override)
    assert "CASE-LOCAL OVERRIDE" in calls
    assert "using case copy" in result.stdout
    assert "using install copy" in result.stdout


@pytest.mark.skipif(_metis_library() is None, reason="needs a built libmetis.so for the partitioner")
def test_launcher_classic_invocation_inside_install(fake_install, tmp_path):
    """Running from inside the installation must behave as before."""
    install, bin_dir = fake_install
    shutil.copy(_metis_library(), install / "libmetis.so")
    (install / "_ianus").mkdir()
    shutil.copytree(REPO_ROOT / "applications" / "q2p1_sse" / "_ianus" / "DIE" / "RH_LOC",
                    install / "_ianus" / "DIE" / "RH_LOC")

    result = _run_launcher(install, bin_dir, install,
                           ["-f", "_ianus/DIE/RH_LOC", "-n", "4", "-y", "e3d_die.yaml", "-a", "0"])
    assert result.returncode == 0, result.stdout + result.stderr
    assert (install / "e3d.log").exists()
    assert (install / "solver_calls.txt").exists()
    assert (install / "_mesh" / "NEWFAC").is_dir()
    # no layout chatter when install and case coincide
    assert "using install copy" not in result.stdout
