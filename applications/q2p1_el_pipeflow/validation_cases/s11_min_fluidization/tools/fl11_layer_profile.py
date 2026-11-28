#!/usr/bin/env python3
"""Profile per cell layer of a case-11 run: gas fraction, pressure gradient and the Ergun gradient.

Reads the ASCII vtu files of the last output step of RUNDIR (_vtk/res_node_*.<step>.vtu, one per
rank; arrays Points, connectivity, ELEpsilonF per cell, Pressure_V per point) of a run on a mesh of
horizontal cell layers, and prints for the lowest layers

  eps_f     mean gas fraction of the cells of the layer (weighted with the cell base area)
  G         pressure gradient over the layer, from the area means of the pressure on its two planes
  Gcl       Ergun gradient at the gas fraction of the layer and the superficial velocity U
  G/Gcl     their ratio
  share     part of the whole pressure drop that arises in the layer

and three sums in units of the bed weight per area W/A: the whole pressure drop; the pressure drop a
bed of height --bed-height would have with the mean gradient of the core layers; the Ergun gradient
summed over all layers that hold particles. The first two layers above the inlet show what the
inlet rows and the wall layering of the packing do to a fixed-bed hold (RUNBOOK section 17). The mesh
must consist of horizontal cell layers (any layer heights, any partition); other meshes are refused.

usage: fl11_layer_profile.py RUNDIR U [--layers 12] [--core 3 8] [--bed-height 0.0238]
                             [--mu 1.8e-5 --rho 1.0 --d 0.87e-3 --wa 126.408]
       fl11_layer_profile.py --self-test
"""
import argparse
import glob
import os
import re
import sys
import tempfile

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


def quad_areas(xy):
    """areas of convex quadrilaterals given by four corners (n, 4, 2) in any order (shoelace after sorting the
    corners by their angle about the centroid)"""
    rel = xy - xy.mean(axis=1, keepdims=True)
    order = np.argsort(np.arctan2(rel[:, :, 1], rel[:, :, 0]), axis=1)
    x = np.take_along_axis(rel[:, :, 0], order, axis=1)
    y = np.take_along_axis(rel[:, :, 1], order, axis=1)
    return 0.5*np.abs((x*np.roll(y, -1, axis=1) - y*np.roll(x, -1, axis=1)).sum(axis=1))


def layer_means(files, ztol=1e-7):
    """(z of the planes, area mean of the pressure on every plane, gas fraction per layer), over all ranks.

    Every cell belongs to one rank and contributes its lower and its upper face, so the result does not
    depend on the partition: plane mean = sum over the faces of area x mean of the four corner pressures (the
    area mean of the bilinear interpolant on parallelogram faces), layer gas fraction = face-area weighted
    mean over the cells. Refuses a mesh whose cell faces are not horizontal or whose layers do not stack."""
    plane, layer = {}, {}
    for f in files:
        a = data_arrays(f)
        pts = a['Points'].reshape(-1, 3)
        con = a['connectivity'].astype(int).reshape(-1, 8)
        z = pts[con, 2]
        order = np.argsort(z, axis=1, kind='stable')
        lower, upper = np.take_along_axis(con, order[:, :4], axis=1), np.take_along_axis(con, order[:, 4:], axis=1)
        zl, zu = pts[lower, 2], pts[upper, 2]
        if np.ptp(zl, axis=1).max() > ztol or np.ptp(zu, axis=1).max() > ztol:
            raise SystemExit(f'{f}: cells without horizontal lower / upper faces (not a mesh of horizontal layers)')
        area = quad_areas(pts[lower, :2])
        for face, zf in ((lower, zl), (upper, zu)):
            pmean = a['Pressure_V'][face].mean(axis=1)
            for key, w, v in zip(np.round(zf.mean(axis=1), 7), area, pmean):
                acc = plane.setdefault(key, [0.0, 0.0])
                acc[0] += w
                acc[1] += w*v
        for key, w, e in zip(zip(np.round(zl.mean(axis=1), 7), np.round(zu.mean(axis=1), 7)), area, a['ELEpsilonF']):
            acc = layer.setdefault(key, [0.0, 0.0])
            acc[0] += w
            acc[1] += w*e
    layers = sorted(layer)
    zp = np.array([layers[0][0]] + [hi for _, hi in layers])
    if any(layers[i][1] != layers[i + 1][0] for i in range(len(layers) - 1)) or sorted(plane) != list(zp):
        raise SystemExit('the cell layers do not stack between common horizontal planes')
    p = np.array([plane[z][1]/plane[z][0] for z in zp])
    eps = np.array([layer[k][1]/layer[k][0] for k in layers])
    return zp, p, eps


def write_vtu(path, pts, cells, eps, pressure):
    """a minimal ASCII vtu file with the arrays this tool reads (self-test)"""
    def block(name, values):
        return f'<DataArray Name="{name}">\n' + ' '.join(repr(float(v)) for v in values) + '\n</DataArray>\n'
    with open(path, 'w') as fh:
        fh.write(block('Pressure_V', pressure) + block('ELEpsilonF', eps)
                 + block('Points', np.ravel(pts)) + block('connectivity', np.ravel(cells)))


def self_test():
    """two columns of cells (x 0..1, 1..2; y 0..1) in two layers of unequal height (planes z = 0, 1, 3):
    the same mesh in one file and split into two files that share the points at x = 1"""
    zs, xs = (0.0, 1.0, 3.0), (0.0, 1.0, 2.0)
    pts = np.array([(x, y, z) for z in zs for y in (0.0, 1.0) for x in xs])
    node = lambda ix, iy, iz: iz*6 + iy*3 + ix
    cells = np.array([[node(ix + dx, dy, iz + dz) for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)]
                      for iz in (0, 1) for ix in (0, 1)])
    eps = np.array([0.4, 0.5, 0.6, 0.7])
    ok = True
    with tempfile.TemporaryDirectory() as d:
        for name, pressure in (('lin', 10.0 - pts[:, 2]), ('quad', pts[:, 0]**2*(1.0 - pts[:, 2]))):
            one = os.path.join(d, name + '_1.vtu')
            write_vtu(one, pts, cells, eps, pressure)
            halves = []
            for k, sel in enumerate(([0, 2], [1, 3])):
                used = np.unique(cells[sel])
                local = {g: i for i, g in enumerate(used)}
                halves.append(os.path.join(d, f'{name}_2{k}.vtu'))
                write_vtu(halves[-1], pts[used], [[local[g] for g in c] for c in cells[sel]], eps[sel], pressure[used])
            zp, p, e = layer_means([one])
            zq, q, f = layer_means(halves)
            ok &= np.allclose(zp, zs) and np.allclose(zq, zs) and np.allclose(p, q, atol=1e-14) and np.allclose(e, f)
            ok &= np.allclose(e, [0.45, 0.65])
            if name == 'lin':
                ok &= np.allclose(-np.diff(p)/np.diff(zp), [1.0, 1.0])
            else:
                ok &= np.allclose(p, [1.5, 0.0, -3.0])
        # a cell with a tilted upper face is refused
        tilted = pts.copy()
        tilted[node(0, 0, 2), 2] += 0.1
        bad = os.path.join(d, 'tilted.vtu')
        write_vtu(bad, tilted, cells, eps, 10.0 - tilted[:, 2])
        try:
            layer_means([bad])
            ok = False
        except SystemExit:
            pass
        ok &= np.allclose(quad_areas(np.array([[(0, 0), (2, 3), (2, 0), (0, 3)]], dtype=float)), [6.0])
    print('fl11_layer_profile self-test', 'PASSED' if ok else 'FAILED')
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rundir', nargs='?')
    ap.add_argument('U', type=float, nargs='?', help='superficial gas velocity of the hold (m/s)')
    ap.add_argument('--self-test', action='store_true', help='run the self-test on a synthetic mesh')
    ap.add_argument('--layers', type=int, default=12, help='number of layers printed')
    ap.add_argument('--core', type=int, nargs=2, default=(3, 8), metavar=('FIRST', 'LAST'),
                    help='layers (1-based, inclusive) whose mean is the core gradient')
    ap.add_argument('--bed-height', type=float, default=0.0238, help='bed height for the core extrapolation (m)')
    ap.add_argument('--mu', type=float, default=1.8e-5)
    ap.add_argument('--rho', type=float, default=1.0)
    ap.add_argument('--d', type=float, default=0.87e-3)
    ap.add_argument('--wa', type=float, default=126.408, help='bed weight per area W/A (Pa)')
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if a.rundir is None or a.U is None:
        ap.error('RUNDIR and U are required')

    step, files = last_step_files(a.rundir)
    zp, p, eps = layer_means(files)
    h = np.diff(zp)
    grad = -np.diff(p)/h
    gcl = ergun_gradient(eps, a.U, a.mu, a.rho, a.d)
    dp = p[0] - p[-1]
    height = '%.2f mm' % (h[0]*1e3) if np.ptp(h) <= 1e-9*h.max() else '%.2f to %.2f mm' % (h.min()*1e3, h.max()*1e3)
    print('%s: step %s, %d ranks, %d layers of %s, U %.2f; p_in %.3f Pa'
          % (os.path.basename(os.path.normpath(a.rundir)), step, len(files), len(eps), height, a.U, p[0]))
    print(' layer z_lo-z_hi mm  eps_f    G Pa/m   Gcl(eps)  G/Gcl   dp share %')
    for i in range(min(a.layers, len(eps))):
        print(' %2d  %5.2f-%5.2f  %.4f  %7.1f  %7.1f  %.3f   %5.2f'
              % (i + 1, zp[i]*1e3, zp[i + 1]*1e3, eps[i], grad[i], gcl[i], grad[i]/gcl[i], 100*grad[i]*h[i]/dp))
    core = slice(a.core[0] - 1, a.core[1])
    print(' core layers %d-%d: eps %.4f G %.1f ; whole dp %.3f Pa = %.4f of W/A ; if every bed layer had the core '
          'gradient over H=%.2f mm: %.4f' % (a.core[0], a.core[1], eps[core].mean(), grad[core].mean(), dp, dp/a.wa,
                                              a.bed_height*1e3, grad[core].mean()*a.bed_height/a.wa))
    print(' sum over layers of Gcl(eps_layer)*h (layers with eps<%.2f): %.4f of W/A'
          % (EPS_EMPTY, (gcl*h)[eps < EPS_EMPTY].sum()/a.wa))


if __name__ == '__main__':
    sys.exit(main())
