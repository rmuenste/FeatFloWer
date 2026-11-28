#!/usr/bin/env python3
"""Profile per cell layer of a case-11 run: gas fraction, pressure gradient and the Ergun gradient.

Reads the ASCII vtu files of the last output step of RUNDIR (_vtk/res_node_*.<step>.vtu, one per
rank; arrays Points, connectivity, ELEpsilonF per cell, Pressure_V per point) of a run on a mesh of
horizontal cell layers, and prints for the lowest layers

  eps_f     mean gas fraction of the cells of the layer
  G         pressure gradient over the layer, from the plane means of the pressure at its two faces
  Gcl       Ergun gradient at the gas fraction of the layer and the superficial velocity U
  G/Gcl     their ratio
  share     part of the whole pressure drop that arises in the layer

and three sums in units of the bed weight per area W/A: the whole pressure drop; the pressure drop a
bed of height --bed-height would have with the mean gradient of the core layers; the Ergun gradient
summed over all layers that hold particles. The first two layers above the inlet show what the
inlet rows and the wall layering of the packing do to a fixed-bed hold (RUNBOOK section 17).

usage: fl11_layer_profile.py RUNDIR U [--layers 12] [--core 3 8] [--bed-height 0.0238]
                             [--mu 1.8e-5 --rho 1.0 --d 0.87e-3 --wa 126.408]
"""
import argparse
import glob
import os
import re

import numpy as np

EPS_EMPTY = 0.99      # a layer with a larger gas fraction holds no bed


def ergun_gradient(eps, u, mu, rho, d):
    """Ergun pressure gradient (Pa/m) of a fixed bed of gas fraction eps at superficial velocity u"""
    return 150*(1 - eps)**2/eps**3*mu*u/d**2 + 1.75*(1 - eps)/eps**3*rho*u*u/d


def last_step_files(rundir):
    files = sorted(glob.glob(os.path.join(rundir, '_vtk', 'res_node_*.vtu')))
    if not files:
        raise SystemExit(f'no _vtk/res_node_*.vtu in {rundir}')
    last = max(re.search(r'\.(\d+)\.vtu$', f).group(1) for f in files)
    return last, [f for f in files if f.endswith('.' + last + '.vtu')]


def data_arrays(path):
    """the named DataArray blocks of an ASCII vtu file as flat float arrays"""
    out, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            m = re.search(r'<DataArray[^>]*Name="([^"]+)"', line)
            if m:
                name, buf = m.group(1), []
            elif '</DataArray>' in line:
                if name:
                    out[name] = np.array(' '.join(buf).split(), dtype=float)
                name = None
            elif name:
                buf.append(line)
    return out


def layer_means(files):
    """(z of the pressure planes, plane means of the pressure, gas fraction per layer), over all ranks"""
    eps, pz = {}, {}
    for f in files:
        a = data_arrays(f)
        pts = a['Points'].reshape(-1, 3)
        con = a['connectivity'].astype(int).reshape(-1, 8)
        for z, e in zip(np.round(pts[con, 2].mean(axis=1), 7), a['ELEpsilonF']):
            eps.setdefault(z, []).append(e)
        for z, p in zip(np.round(pts[:, 2], 7), a['Pressure_V']):
            pz.setdefault(z, []).append(p)
    zp = np.array(sorted(pz))
    return zp, np.array([np.mean(pz[z]) for z in zp]), np.array([np.mean(eps[z]) for z in sorted(eps)])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rundir')
    ap.add_argument('U', type=float, help='superficial gas velocity of the hold (m/s)')
    ap.add_argument('--layers', type=int, default=12, help='number of layers printed')
    ap.add_argument('--core', type=int, nargs=2, default=(3, 8), metavar=('FIRST', 'LAST'),
                    help='layers (1-based, inclusive) whose mean is the core gradient')
    ap.add_argument('--bed-height', type=float, default=0.0238, help='bed height for the core extrapolation (m)')
    ap.add_argument('--mu', type=float, default=1.8e-5)
    ap.add_argument('--rho', type=float, default=1.0)
    ap.add_argument('--d', type=float, default=0.87e-3)
    ap.add_argument('--wa', type=float, default=126.408, help='bed weight per area W/A (Pa)')
    a = ap.parse_args()

    step, files = last_step_files(a.rundir)
    zp, p, eps = layer_means(files)
    h = zp[1] - zp[0]
    grad = -(p[1:] - p[:-1])/h
    gcl = ergun_gradient(eps, a.U, a.mu, a.rho, a.d)
    dp = p[0] - p[-1]
    print('%s: step %s, %d ranks, %d layers of %.2f mm, U %.2f; p_in %.3f Pa'
          % (os.path.basename(os.path.normpath(a.rundir)), step, len(files), len(eps), h*1e3, a.U, p[0]))
    print(' layer z_lo-z_hi mm  eps_f    G Pa/m   Gcl(eps)  G/Gcl   dp share %')
    for i in range(min(a.layers, len(eps))):
        print(' %2d  %5.2f-%5.2f  %.4f  %7.1f  %7.1f  %.3f   %5.2f'
              % (i + 1, zp[i]*1e3, zp[i + 1]*1e3, eps[i], grad[i], gcl[i], grad[i]/gcl[i], 100*grad[i]*h/dp))
    core = slice(a.core[0] - 1, a.core[1])
    print(' core layers %d-%d: eps %.4f G %.1f ; whole dp %.3f Pa = %.4f of W/A ; if every bed layer had the core '
          'gradient over H=%.2f mm: %.4f' % (a.core[0], a.core[1], eps[core].mean(), grad[core].mean(), dp, dp/a.wa,
                                              a.bed_height*1e3, grad[core].mean()*a.bed_height/a.wa))
    print(' sum over layers of Gcl(eps_layer)*h (layers with eps<%.2f): %.4f of W/A'
          % (EPS_EMPTY, (gcl[eps < EPS_EMPTY]*h).sum()/a.wa))


if __name__ == '__main__':
    main()
