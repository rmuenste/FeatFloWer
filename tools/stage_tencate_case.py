#!/usr/bin/env python3
"""Stage a ten Cate E1-E4 sedimentation rundir from an existing template.

Clones a working rundir (default: q2p1_dns_rundir_e4_l3), excludes run
outputs, and applies the per-case fluid properties from ten Cate (2002)
Table I to BOTH sides that must agree (the d0_visc lesson):

  _data/q2p1_param.dat : Prop@Density, Prop@Viscosity, SimPar@MaxNumStep,
                         SimPar@MaxMeshLevel, SimPar@TimeStep
  example.json         : fluidDensity_, fluidViscosity_

Also writes a job.sbatch. Step counts cover the PIV time window of each
case (tc-ref/) plus margin at dt = 1 ms; --dt rescales them.

Examples:
  stage_tencate_case.py --case E1
  stage_tencate_case.py --case E4 --level 4 --dt 0.0005   # dt-ladder run
  for c in E1 E2 E3; do tools/stage_tencate_case.py --case $c; done
"""
import argparse
import os
import re
import shutil
import sys

# rho_f [kg/m^3], mu [Pa s], steps at dt=1ms (PIV window end + margin)
CASES = {
    'E1': dict(rho=970.0, mu=0.373, steps=4300),   # PIV ends 4.11 s
    'E2': dict(rho=965.0, mu=0.212, steps=2700),   # PIV ends 2.50 s
    'E3': dict(rho=962.0, mu=0.113, steps=1800),   # PIV ends 1.63 s
    'E4': dict(rho=960.0, mu=0.058, steps=1300),   # PIV ends 1.20 s
}

EXCLUDE = {'_dump', '_gmv', '_vtk', 'testresults', 'paraview', '_sol',
           'solution'}
# run-output dirs the solver expects to exist in a fresh rundir
RECREATE = ('_dump', '_gmv', '_vtk', 'testresults')


def edit_file(path, subs):
    """Apply (pattern, replacement, expected_count) regex subs to a file."""
    text = open(path).read()
    for pat, rep, n in subs:
        text, cnt = re.subn(pat, rep, text, flags=re.M)
        if cnt != n:
            sys.exit('%s: pattern %r matched %d times, expected %d'
                     % (path, pat, cnt, n))
    open(path, 'w').write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True, choices=sorted(CASES))
    ap.add_argument('--level', type=int, default=3)
    ap.add_argument('--dt', type=float, default=0.0010)
    ap.add_argument('--steps', type=int, default=None,
                    help='override step count (default: case table scaled by dt)')
    ap.add_argument('--template', default='q2p1_dns_rundir_e4_l3')
    ap.add_argument('--dest', default=None)
    ap.add_argument('--backupfreq', type=int, default=100000,
                    help='SimPar@BackUpFreq: counts OUTPUT EVENTS (a dump is '
                         'written each time a VTK output happens) and '
                         'SimPar@BackUpNum rotates a bounded slot set; high '
                         'value = effectively no restart dumps, correct for '
                         'gate runs that are never restarted (default '
                         '100000)')
    ap.add_argument('--outputlevel', default='MAX-1',
                    help='SimPar@OutputLevel; MAX-1 suffices for '
                         'numerics/physics runs, MAX+1 (once-refined mesh, '
                         '~8x file size) only for figure-quality viz '
                         '(default MAX-1)')
    ap.add_argument('--viz', action='store_true',
                    help='keep the template SimPar@OutputFreq (VTK output '
                         'on); default OFF sets OutputFreq = 1d6, i.e. '
                         'effectively no VTK')
    ap.add_argument('--np', type=int, default=32)
    ap.add_argument('--partition', default='short')
    ap.add_argument('--walltime', default='02:00:00')
    ap.add_argument('--mem', default='25G',
                    help='sbatch --mem; L4 needs 60G')
    ap.add_argument('--binary', default='build-dns-pe-serial/applications/'
                    'q2p1_bench_sedimentation/q2p1_bench_sedimentation')
    a = ap.parse_args()

    c = CASES[a.case]
    steps = a.steps if a.steps else int(round(c['steps'] * 0.0010 / a.dt))
    tag = '%s_l%d' % (a.case.lower(), a.level)
    if a.dt != 0.0010:
        tag += '_dt%s' % ('%g' % (a.dt * 1e3)).replace('.', 'p')
    dest = a.dest or 'q2p1_dns_rundir_' + tag
    if any(p.startswith('build') for p in
           os.path.normpath(os.path.abspath(dest)).split(os.sep)):
        sys.exit('refusing to stage into a build tree: %s '
                 '(buried rundirs under build-*/ are how 60 GB of VTK '
                 'gets lost)' % dest)
    if os.path.exists(dest):
        sys.exit('%s already exists; remove it or pass --dest' % dest)

    shutil.copytree(a.template, dest, symlinks=True,
                    ignore=lambda d, names: [n for n in names
                                             if n in EXCLUDE
                                             or n.endswith(('.log', '.png'))
                                             or n.startswith('slurm-')])
    for d in RECREATE:
        os.makedirs(os.path.join(dest, d), exist_ok=True)

    # mu in the deck as integer milli-Pa-s (matches 58d-3 style); rho as NNNd0
    mu_milli = round(c['mu'] * 1000)
    if abs(mu_milli - c['mu'] * 1000) > 1e-9:
        sys.exit('mu %g not representable as integer milli-Pa-s' % c['mu'])
    deck = os.path.join(dest, '_data', 'q2p1_param.dat')
    # NOTE: numeric Prop@ lines must NOT carry inline comments (list-directed read)
    deck_subs = [
        (r'^SimPar@MaxMeshLevel = .*$',
         'SimPar@MaxMeshLevel = %d' % a.level, 1),
        (r'^SimPar@TimeStep = .*$',
         'SimPar@TimeStep = %sd0' % ('%.8g' % a.dt), 1),
        (r'^SimPar@MaxNumStep = .*$',
         'SimPar@MaxNumStep = %d' % steps, 1),
        # BackUpFreq counts output events (dump written per VTK output);
        # BackUpNum rotates bounded slots, so dumps are size-bounded.
        # High BackUpFreq = effectively no restart dumps (right for gate runs).
        (r'^SimPar@BackUpFreq = .*$',
         'SimPar@BackUpFreq = %d' % a.backupfreq, 1),
        # OutputLevel MAX+1 writes on a once-refined mesh (~8x file size);
        # MAX-1 suffices for numerics/physics runs.
        (r'^SimPar@OutputLevel = .*$',
         'SimPar@OutputLevel = %s' % a.outputlevel, 1),
        (r'^Prop@Density = .*$',
         'Prop@Density = %dd0,30d0' % round(c['rho']), 1),
        (r'^Prop@Viscosity = .*$',
         'Prop@Viscosity = %dd-3,1d0' % mu_milli, 1),
    ]
    if not a.viz:
        # viz off by default: effectively no VTK output
        deck_subs.append((r'^SimPar@OutputFreq = .*$',
                          'SimPar@OutputFreq = 1d6', 1))
    edit_file(deck, deck_subs)
    edit_file(os.path.join(dest, 'example.json'), [
        (r'"fluidDensity_":\s*[0-9.eE+-]+',
         '"fluidDensity_": %.1f' % c['rho'], 1),
        (r'"fluidViscosity_":\s*[0-9.eE+-]+',
         '"fluidViscosity_": %.3f' % c['mu'], 1),
    ])

    root = os.path.abspath('.')
    dest_abs = os.path.join(root, dest)
    with open(os.path.join(dest, 'job.sbatch'), 'w') as f:
        f.write('''#!/bin/bash
#SBATCH --job-name=dns_%s
#SBATCH --partition=%s
#SBATCH --nodes=1
#SBATCH --ntasks=%d
#SBATCH --cpus-per-task=1
#SBATCH --mem=%s
#SBATCH --prefer=nx
#SBATCH --time=%s
#SBATCH --output=%s/slurm-%%j.out
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
cd %s
FREE_GB=$(df -BG --output=avail /data/warehouse17 | tail -1 | tr -dc 0-9); [ "$FREE_GB" -lt "${FF_MIN_FREE_GB:-100}" ] && { echo "REFUSING to run: only ${FREE_GB}G free on warehouse17"; exit 1; }
mpirun -np %d %s > run_slurm.log 2>&1
''' % (tag, a.partition, a.np, a.mem, a.walltime, dest_abs, dest_abs,
       a.np, os.path.join(root, a.binary)))

    print('staged %s: rho_f=%g mu=%g level=%d dt=%g steps=%d (t_end=%g s)'
          % (dest, c['rho'], c['mu'], a.level, a.dt, steps, steps * a.dt))


if __name__ == '__main__':
    main()
