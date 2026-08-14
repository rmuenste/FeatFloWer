#!/usr/bin/env python3
"""Generate a random non-overlapping fixed-sphere configuration for D3.1.

Random sequential addition in the periodic unit cell [0,1]^3 with
minimum-image distance checks. Surface-to-surface minimum gap is
`--mingap` diameters (Beetstra-style hard-sphere configs use ~0.05d).

Writes one "x y z" line per sphere (readVectorsFromFile format) plus a
'# meta' header line (ignored by the reader? NO - the PE reader parses
every line, so metadata goes to a sidecar .meta file instead).

Usage:
  gen_random_array.py --n 16 --phi 0.05 --seed 1 --out particles.xyz
  (radius is derived from n and phi; printed and stored in the sidecar)
"""
import argparse, math, random, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--phi', type=float, required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--mingap', type=float, default=0.05, help='min surface gap in diameters')
    ap.add_argument('--out', required=True)
    ap.add_argument('--maxtries', type=int, default=2000000)
    a = ap.parse_args()

    r = (a.phi * 3.0 / (4.0 * math.pi * a.n)) ** (1.0 / 3.0)
    dmin = 2.0 * r * (1.0 + a.mingap)     # center-to-center minimum
    rng = random.Random(a.seed)
    pts = []
    tries = 0
    while len(pts) < a.n:
        tries += 1
        if tries > a.maxtries:
            sys.exit('FAILED: placed %d/%d after %d tries (phi too high for RSA?)'
                     % (len(pts), a.n, tries))
        p = (rng.random(), rng.random(), rng.random())
        ok = True
        for q in pts:
            d2 = 0.0
            for c in range(3):
                d = abs(p[c] - q[c])
                d = min(d, 1.0 - d)       # minimum image
                d2 += d * d
            if d2 < dmin * dmin:
                ok = False
                break
        if ok:
            pts.append(p)

    # Periodic image copies: the FBM indicator is NOT periodic (datasheet
    # d31_periodic_indicator) - a sphere crossing a face loses the part
    # outside [0,1]^3. For FIXED arrays explicit image spheres restore the
    # geometry exactly. One line per image; .meta records parent indices.
    images = []
    for ip, p in enumerate(pts):
        shifts = [[0.0], [0.0], [0.0]]
        for c in range(3):
            if p[c] < r:
                shifts[c].append(1.0)
            elif p[c] > 1.0 - r:
                shifts[c].append(-1.0)
        for sx in shifts[0]:
            for sy in shifts[1]:
                for sz in shifts[2]:
                    if sx == sy == sz == 0.0:
                        continue
                    images.append((ip, (p[0] + sx, p[1] + sy, p[2] + sz)))

    with open(a.out, 'w') as f:
        for p in pts:
            f.write('%.12f %.12f %.12f\n' % p)
        for _, q in images:
            f.write('%.12f %.12f %.12f\n' % q)
    with open(a.out + '.meta', 'w') as f:
        f.write('n=%d phi=%.6f radius=%.9f mingap=%.3fd seed=%d tries=%d nimages=%d\n'
                % (a.n, a.phi, r, a.mingap, a.seed, tries, len(images)))
        for k, (ip, _) in enumerate(images):
            f.write('image %d parent %d\n' % (a.n + k, ip))
    print('n=%d phi=%.4f -> radius=%.6f (D=%.4f), mingap=%.3fd, seed=%d, %d tries, %d image spheres'
          % (a.n, a.phi, r, 2 * r, a.mingap, a.seed, tries, len(images)))


if __name__ == '__main__':
    main()
