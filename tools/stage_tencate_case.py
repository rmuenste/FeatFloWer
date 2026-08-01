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

EXCLUDE = {'_dump', '_gmv', '_vtk', 'testresults'}


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
    if os.path.exists(dest):
        sys.exit('%s already exists; remove it or pass --dest' % dest)

    shutil.copytree(a.template, dest, symlinks=True,
                    ignore=lambda d, names: [n for n in names
                                             if n in EXCLUDE
                                             or n.endswith(('.log', '.png'))
                                             or n.startswith('slurm-')])
    for d in EXCLUDE:
        os.makedirs(os.path.join(dest, d), exist_ok=True)

    # mu in the deck as integer milli-Pa-s (matches 58d-3 style); rho as NNNd0
    mu_milli = round(c['mu'] * 1000)
    if abs(mu_milli - c['mu'] * 1000) > 1e-9:
        sys.exit('mu %g not representable as integer milli-Pa-s' % c['mu'])
    deck = os.path.join(dest, '_data', 'q2p1_param.dat')
    # NOTE: numeric Prop@ lines must NOT carry inline comments (list-directed read)
    edit_file(deck, [
        (r'^SimPar@MaxMeshLevel = .*$',
         'SimPar@MaxMeshLevel = %d' % a.level, 1),
        (r'^SimPar@TimeStep = .*$',
         'SimPar@TimeStep = %sd0' % ('%.8g' % a.dt), 1),
        (r'^SimPar@MaxNumStep = .*$',
         'SimPar@MaxNumStep = %d' % steps, 1),
        (r'^Prop@Density = .*$',
         'Prop@Density = %dd0,30d0' % round(c['rho']), 1),
        (r'^Prop@Viscosity = .*$',
         'Prop@Viscosity = %dd-3,1d0' % mu_milli, 1),
    ])
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
mpirun -np %d %s > run_slurm.log 2>&1
''' % (tag, a.partition, a.np, a.mem, a.walltime, dest_abs, dest_abs,
       a.np, os.path.join(root, a.binary)))

    print('staged %s: rho_f=%g mu=%g level=%d dt=%g steps=%d (t_end=%g s)'
          % (dest, c['rho'], c['mu'], a.level, a.dt, steps, steps * a.dt))


if __name__ == '__main__':
    main()
