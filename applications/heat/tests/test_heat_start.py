"""Launcher integration tests: no MPI installation or compiled solver required."""
import json
import hashlib
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'tools/e3d_scripts'))
from e3d_layout import RunLayout


def test_repository_heat_case_is_complete():
    example = ROOT / 'applications/heat/_ianus/HEAT'
    assert (example / 'heat.s3d').is_file()
    assert (example / 'meshDir/file.prj').is_file()
    xml_text = (example / 'sampleRigidBody.xml').read_text().lstrip()
    root = ET.fromstring(xml_text)
    description = root.find('BoundaryDescription')
    assert description is not None
    boundaries = description.findall('BoundaryShape')
    assert int(description.attrib['ncomponents']) == len(boundaries)
    assert boundaries
    for boundary in boundaries:
        assert (example / boundary.attrib['meshFile']).is_file()


def executable(path, body):
    path.write_text('#!' + sys.executable + '\n' + body)
    path.chmod(0o755)


@pytest.fixture
def setup(tmp_path):
    install = tmp_path / 'install space'
    install.mkdir()
    shutil.copy(ROOT / 'applications/heat/heat_start.py', install)
    shutil.copy(ROOT / 'tools/e3d_scripts/e3d_layout.py', install)
    for rel in ('_data/MG.dat', '_data/q2p1_param.dat',
                'start/sampleRigidBody.xml', 'start/data.TXT'):
        dest = install / rel
        dest.parent.mkdir(exist_ok=True)
        dest.write_text('SimPar@PartitionFormat = json\n')
    (install / 'partitioner.py').write_text(
        'import json\nfrom pathlib import Path\n'
        'def partition(*args, **kwargs):\n'
        '    Path("partition.json").write_text(json.dumps([args, kwargs]))\n')
    executable(install / 's3d_mesher',
               'import os, sys\nfrom pathlib import Path\n'
               'assert sys.argv[1:] == ["-a", "heat"]\n'
               'if os.environ.get("GENERATE"):\n'
               '    Path("_data/meshDir").mkdir()\n'
               '    Path("_data/meshDir/file.prj").write_text("generated")\n'
               'sys.exit(int(os.environ.get("MESHER_EXIT", "0")))\n')
    executable(install / 'heat',
               'import os, sys\nfrom pathlib import Path\n'
               'import xml.etree.ElementTree as ET\n'
               'for ref in Path("_data/heat.s3d").read_text().splitlines():\n'
               '    if ref.startswith("/"): assert Path(ref).is_file(), ref\n'
               'for node in ET.parse("start/sampleRigidBody.xml").findall("./BoundaryDescription/BoundaryShape"):\n'
               '    assert Path(node.attrib["meshFile"]).is_file()\n'
               'Path("_vtk/result").write_text(str(Path.cwd()))\n'
               'sys.exit(int(os.environ.get("SOLVER_EXIT", "0")))\n')
    bins = tmp_path / 'bins'
    bins.mkdir()
    for name in ('mpirun', 'srun'):
        executable(bins / name,
                   'import os, sys, json\nfrom pathlib import Path\n'
                   'Path("launch.json").write_text(json.dumps(sys.argv))\n'
                   'args = sys.argv[1:]\n'
                   'if args[0] == "-np": args = args[2:]\n'
                   'os.execv(args[0], args)\n')
    case = tmp_path / 'case'
    case.mkdir()
    project = case / 'input'
    project.mkdir()
    (project / 'meshDir').mkdir()
    (project / 'meshDir/file.prj').write_text('fallback')
    (project / 'geometry').mkdir()
    (project / 'geometry/steel.off').write_text('geometry')
    (case / 'keep.OFF').write_text('original')
    shared = tmp_path / 'shared.off'
    shared.write_text('shared')
    (project / 'heat.s3d').write_text(
        '[E3DGeometryData/Machine/Element_1]\nScrewOFF(2) =\n'
        '# geometry list\ngeometry/steel.off\n' + str(shared) + '\n'
        'SensorOFF(1) =\ngeometry/steel.off\nTemperature = 320\n')
    (project / 'sampleRigidBody.xml').write_text(
        '\n<?xml version="1.0"?><root><ApplicationSettings cgalConfigFile="mesh_names.offs"/>'
        '<RigidBodyList><RigidBody meshFile="unused.off"/></RigidBodyList>'
        '<BoundaryDescription ncomponents="1">'
        '<BoundaryShape type="1" meshFile="wall_1.off"/>'
        '</BoundaryDescription></root>')
    (project / 'wall_1.off').write_text('input boundary geometry\n')
    env = dict(os.environ, PATH=str(bins) + os.pathsep + os.environ['PATH'])
    env.pop('FF_HEAT_HOME', None)
    for name in ('SLURM_STEP_NUM_NODES', 'SLURM_JOB_NUM_NODES', 'SLURM_JOB_NODELIST',
                 'HOSTFILE', 'PBS_NODEFILE', 'OMPI_MCA_rankfile'):
        env.pop(name, None)
    return install, case, project, env


def run(setup, *args, cwd=None, extra=None):
    install, case, project, env = setup
    return subprocess.run(
        [sys.executable, str(install / 'heat_start.py'), *args],
        cwd=cwd or case.parent, env=dict(env, **(extra or {})),
        text=True, capture_output=True, timeout=30)


@pytest.mark.parametrize('srun', [False, True])
@pytest.mark.parametrize('generate', [False, True])
def test_separate_case(setup, srun, generate):
    install, case, project, env = setup
    (case / '_data').mkdir()
    params = case / '_data/q2p1_param.dat'
    params.write_text('SimPar@PartitionFormat = legacy\n')
    (case / 'start').mkdir()
    (case / 'start/sampleRigidBody.xml').write_text('<stale-case-copy/>\n')
    (case / 'wall_1.off').write_text('stale boundary geometry\n')
    result = run(setup, '-C', str(case), '-f', 'input', '-n', '4',
                 *(['-u'] if srun else []),
                 extra={'GENERATE': '1'} if generate else {'MESHER_EXIT': '7'})
    assert result.returncode == 0, result.stdout + result.stderr
    assert params.read_text().endswith('legacy\n')
    assert str(project / 'geometry/steel.off') in (case / '_data/heat.s3d').read_text()
    assert '\ngeometry/steel.off\n' in (project / 'heat.s3d').read_text()
    xml = ET.parse(case / 'start/sampleRigidBody.xml')
    assert xml.find('./BoundaryDescription/BoundaryShape').get('meshFile') == str(project / 'wall_1.off')
    assert xml.find('./ApplicationSettings').get('cgalConfigFile') == 'mesh_names.offs'
    assert xml.find('./RigidBodyList/RigidBody').get('meshFile') == 'unused.off'
    assert (case / 'wall_1.off').read_text() == 'stale boundary geometry\n'
    assert (project / 'geometry/steel.off').read_text() == 'geometry'
    assert (case / 'keep.OFF').read_text() == 'original'
    assert not (case / 'steel.off').exists()
    assert not (install / '_vtk').exists()
    assert (project / 'meshDir/file.prj').read_text() == 'fallback'
    args, kwargs = json.loads((case / 'partition.json').read_text())
    assert args[:3] == [3, 1, 1]
    assert kwargs == {'partition_format': 'legacy'}
    launch = json.loads((case / 'launch.json').read_text())
    assert launch[-1] == str(install / 'heat')
    assert Path(launch[0]).name == ('srun' if srun else 'mpirun')


def test_solver_failure_preserves_geometry(setup):
    _, case, project, _ = setup
    result = run(setup, '-f', 'input', '-n', '2', cwd=case, extra={'SOLVER_EXIT': '9'})
    assert result.returncode == 9
    assert (case / 'keep.OFF').exists()
    assert (case / '_data/q2p1_param.dat').read_text().endswith('json\n')


@pytest.mark.parametrize('ranks', ['1', '-1', 'abc'])
def test_invalid_ranks_do_not_prepare_case(setup, ranks):
    _, case, project, _ = setup
    result = run(setup, '-C', str(case), '-f', str(project), '-n', ranks)
    assert result.returncode == 2
    assert not (case / '_data').exists()


def test_overlap_preserves_mesh(setup):
    _, case, project, _ = setup
    shutil.copytree(project, case / '_data')
    result = run(setup, '-C', str(case), '-f', str(case / '_data'), '-n', '2')
    assert result.returncode == 2
    assert 'separate from generated' in result.stderr
    assert (case / '_data/meshDir/file.prj').read_text() == 'fallback'


def test_symlink_preserves_external_data(setup):
    _, case, project, _ = setup
    (case / '_data').symlink_to(project, target_is_directory=True)
    result = run(setup, '-C', str(case), '-f', str(project), '-n', '2')
    assert result.returncode == 2
    assert 'symlink' in result.stderr
    assert (project / 'meshDir/file.prj').exists()


def test_classic_and_install_lookup(setup):
    install, case, project, _ = setup
    (install / 'example').symlink_to(project, target_is_directory=True)
    result = run(setup, '-f', 'example', '-n', '2', cwd=install)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (install / '_vtk/result').read_text() == str(install)
    result = run(setup, '-C', str(case), '-f', 'example', '-n', '2')
    assert result.returncode == 0, result.stdout + result.stderr


def test_same_configuration_rejected_without_modifying_source(setup):
    _, case, project, _ = setup
    (case / '_data').mkdir()
    shutil.copy(project / 'heat.s3d', case / '_data/heat.s3d')
    shutil.copy(project / 'sampleRigidBody.xml', case / '_data/sampleRigidBody.xml')
    shutil.copytree(project / 'geometry', case / '_data/geometry')
    shutil.copy(project / 'wall_1.off', case / '_data/wall_1.off')
    result = run(setup, '-C', str(case), '-f', str(case / '_data'), '-n', '2',
                 extra={'GENERATE': '1'})
    assert result.returncode == 2
    assert 'separate from generated' in result.stderr
    assert (case / '_data/heat.s3d').read_bytes() == (project / 'heat.s3d').read_bytes()


@pytest.mark.parametrize('missing', ['input', 'rigid-body', 'default', 'mesh'])
def test_missing_files(setup, missing):
    install, case, project, _ = setup
    target = {'input': project / 'heat.s3d',
              'rigid-body': project / 'sampleRigidBody.xml',
              'default': install / '_data/MG.dat',
              'mesh': project / 'meshDir/file.prj'}[missing]
    target.unlink()
    result = run(setup, '-C', str(case), '-f', str(project), '-n', '2')
    assert result.returncode == 2
    assert not (case / 'launch.json').exists()


def test_custom_layout_environment_and_seeds(setup, monkeypatch):
    install, case, _, _ = setup
    monkeypatch.setenv('FF_HEAT_HOME', str(install))
    monkeypatch.setenv('FF_GENDIE_HOME', '/unrelated')
    layout = RunLayout(case_dir=case, install_dir_env='FF_HEAT_HOME',
                       runtime_dirs=('_data',), seed_files=('_data/q2p1_param.dat',))
    assert layout.install_dir == install
    layout.prepare_case()
    params = case / '_data/q2p1_param.dat'
    params.write_text('override')
    layout.prepare_case()
    assert params.read_text() == 'override'
    assert not (case / 'start').exists()


def test_relative_hostfile_and_invocation_input(setup):
    install, case, project, _ = setup
    invocation = case.parent
    (invocation / 'hosts').write_text('node-a\nnode-b\n')
    shutil.copytree(project, invocation / 'input')
    (project / 'heat.s3d').write_text('missing geometry if case input wins')
    result = run(setup, '-C', case.name, '-f', 'input', '-n', '4',
                 extra={'HOSTFILE': 'hosts', 'FF_HEAT_HOME': str(install)})
    assert result.returncode == 0, result.stdout + result.stderr
    args, _ = json.loads((case / 'partition.json').read_text())
    assert args[2] == 2


def test_legacy_partitioner(setup):
    install, case, project, _ = setup
    (install / 'partitioner.py').write_text(
        'from pathlib import Path\n'
        'def partition(a, b, c, d, e):\n'
        '    Path("legacy-used").touch()\n')
    result = run(setup, '-f', 'input', '-n', '2', cwd=case)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (case / 'legacy-used').exists()


def test_staging_preserves_parameters(tmp_path):
    if not shutil.which('cmake'):
        pytest.skip('requires cmake')
    runtime = tmp_path / 'stage'
    cmd = ['cmake', '-DHEAT_SOURCE=' + str(ROOT / 'applications/heat'),
           '-DHEAT_RUNTIME=' + str(runtime), '-DFF_SOURCE=' + str(ROOT),
           '-P', str(ROOT / 'applications/heat/stage_defaults.cmake')]
    subprocess.run(cmd, check=True, capture_output=True)
    params = runtime / '_data/q2p1_param.dat'
    params.write_text('user settings')
    subprocess.run(cmd, check=True, capture_output=True)
    assert params.read_text() == 'user settings'
    assert (runtime / '_data/MG.dat').is_file()
    assert (runtime / 'start/data.TXT').is_file()


@pytest.mark.parametrize('problem', ['missing', 'space', 'long', 'incomplete', 'xml'])
def test_geometry_preflight_preserves_case(setup, problem):
    _, case, project, _ = setup
    if problem == 'xml':
        (project / 'sampleRigidBody.xml').write_text('<broken')
    else:
        ref = 'absent.off'
        if problem == 'space':
            ref = 'has space.off'
            (project / ref).touch()
        elif problem == 'long':
            ref = 'x' * 201 + '.off'
            (project / ref).touch()
        value = '' if problem == 'incomplete' else ref + '\n'
        (project / 'heat.s3d').write_text(
            '[E3DGeometryData/Machine/Element_1]\nScrewOFF(1) =\n' + value)
    result = run(setup, '-C', str(case), '-f', str(project), '-n', '2')
    assert result.returncode == 2
    assert not (case / '_data').exists()


def test_absolute_xml_paths_and_two_runtime_cases(setup):
    _, case, project, _ = setup
    boundary = project / 'boundary & shared.off'
    boundary.write_text('boundary')
    xml = ET.ElementTree(ET.fromstring((project / 'sampleRigidBody.xml').read_text().lstrip()))
    xml.find('./BoundaryDescription/BoundaryShape').set('meshFile', str(boundary))
    xml.write(project / 'sampleRigidBody.xml')
    originals = {p: p.read_bytes() for p in (project / 'heat.s3d', project / 'sampleRigidBody.xml')}
    for target in (case, case.parent / 'second run'):
        result = run(setup, '-C', str(target), '-f', str(project), '-n', '2')
        assert result.returncode == 0, result.stdout + result.stderr
        staged = ET.parse(target / 'start/sampleRigidBody.xml')
        assert staged.find('./BoundaryDescription/BoundaryShape').get('meshFile') == str(boundary)
        assert not (target / 'wall_1.off').exists()
    assert {p: p.read_bytes() for p in originals} == originals


@pytest.mark.skipif(not os.environ.get('HEAT_SMOKE_INSTALL'),
                    reason='set HEAT_SMOKE_INSTALL to an installed bin/heat directory')
def test_real_heat_four_ranks(tmp_path):
    """Two EWIKON steps; heat suppresses visualization output on step one."""
    install = Path(os.environ['HEAT_SMOKE_INSTALL']).resolve()
    example = install / '_ianus/HEAT'
    if not example.is_dir():
        pytest.fail('Real heat smoke test requires the installed _ianus/HEAT case')
    case = tmp_path / 'case'
    project = example
    source_files = [p for p in example.rglob('*') if p.is_file()]
    before = {p: hashlib.sha256(p.read_bytes()).digest() for p in source_files}
    (case / '_data').mkdir(parents=True)
    parameters = (install / '_data/q2p1_param.dat').read_text()
    for key, value in {'SimPar@MaxMeshLevel': '1', 'SimPar@MaxNumStep': '2',
                       'SimPar@OutputFreq': '0.5d0', 'SimPar@MaxSimTime': '1.0d0',
                       'Velo@MGCrsSolverType': '1', 'Pres@MGCrsSolverType': '1'}.items():
        parameters = re.sub(r'^' + re.escape(key) + r'\s*=.*$', key + ' = ' + value,
                            parameters, flags=re.MULTILINE)
    (case / '_data/q2p1_param.dat').write_text(parameters)
    env = dict(os.environ, PYTHONUNBUFFERED='1', OMP_NUM_THREADS='1')
    env.pop('FF_HEAT_HOME', None)
    env['LD_LIBRARY_PATH'] = '/usr/lib/x86_64-linux-gnu:' + env.get('LD_LIBRARY_PATH', '')
    with (tmp_path / 'heat.log').open('w') as log:
        result = subprocess.run(
            [sys.executable, str(install / 'heat_start.py'), '-C', str(case),
             '-f', str(project), '-n', '4'], cwd=tmp_path, env=env,
            stdout=log, stderr=subprocess.STDOUT, timeout=300)
    assert result.returncode == 0, str(tmp_path / 'heat.log')
    assert (case / '_data/prot.txt').is_file()
    assert any(p.suffix in ('.vtk', '.vtu') for p in (case / '_vtk').rglob('*'))
    assert not any((install / '_vtk').iterdir())
    assert not list(case.rglob('*.off'))
    assert not list(case.rglob('*.OFF'))
    assert {p: hashlib.sha256(p.read_bytes()).digest() for p in source_files} == before
    xml = ET.parse(case / 'start/sampleRigidBody.xml')
    boundaries = xml.findall('./BoundaryDescription/BoundaryShape')
    assert len(boundaries) == 30
    assert all(Path(b.get('meshFile')).is_absolute() for b in boundaries)
