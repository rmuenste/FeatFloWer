#!/usr/bin/env python3
"""Interpolated particle movie of a batch-settling run (spec 10) from FeatFloWer particle logs.
usage: batch_movie_ff.py RUNDIR OUTDIR [SUBSTEPS [FRAME0 [FRAME1]]]
Reads RUNDIR/particles.xyz (t = 0, id = line number) and RUNDIR/el_particles_rank*.log
(time ip f(3) t(3) p(3) v(3), every 0.1 s), interpolates each particle with a cubic Hermite spline
(positions and velocities of the two neighbouring log frames) at SUBSTEPS (default 8) sub-frames
per log interval and writes OUTDIR/frame_%05d.png (2560 x 1440).
Panels: column cut at the wall (y < 2.5 mm) and at the mid-plane (2.5 mm slab), both coloured by
the vertical particle velocity (grey band = the top-down slab); solid-fraction profile phi(z) (2 mm bins, centre counts); top-down
view of a 2.5 mm slab; mean settling speed of the wall ring and of the interior over time.
Encode: ffmpeg -framerate 30 -i OUTDIR/frame_%05d.png -c:v libx264 -crf 18 -pix_fmt yuv420p out.mp4"""
import os, sys, glob, math, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

run, out = sys.argv[1], sys.argv[2]
nsub = int(sys.argv[3]) if len(sys.argv) > 3 else 8
f0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
f1 = int(sys.argv[5]) if len(sys.argv) > 5 else 10**9
os.makedirs(out, exist_ok=True)

# ---- case constants (spec 10, phi0 0.15) and view choices
L, H, D = 20.0, 200.0, 1.0                 # mm
VREF = 8.172                               # mm/s, reference settling / upper-front speed
H0 = 160.0                                 # mm, initial suspension height
SLAB = 2.5                                 # mm, slab thickness of every cut (one fluid element)
ZS0 = 80.0                                 # mm, lower edge of the top-down slab
CORE = (40.0, 100.0)                       # mm, height range of the wall/interior statistics
VLIM = 20.0                                # mm/s, colour range
SURF, INK, INK2, GRID = '#1a1a19', '#f0efec', '#c3c2b7', '#3a3a37'
C_WALL, C_IN = '#199e70', '#c98500'        # categorical: wall ring, interior
cmap = LinearSegmentedColormap.from_list('vz', ['#4d9bff', '#2f62a8', '#77776f', '#a8523d', '#ff7a5c'])
norm = Normalize(-VLIM, VLIM)

# ---- data: frames[k] = (t, pos(N,3) mm, vel(N,3) mm/s), particles ordered by id
xyz = np.loadtxt(os.path.join(run, 'particles.xyz')) * 1e3
N = len(xyz)
raw = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log'))):
    a = np.loadtxt(f, comments='#', usecols=(0, 1, 8, 9, 10, 11, 12, 13), ndmin=2)
    for t in np.unique(np.round(a[:, 0], 6)):
        raw[t].append(a[np.round(a[:, 0], 6) == t])
frames = [(0.0, xyz, np.zeros_like(xyz))]
for t in sorted(raw):
    a = np.vstack(raw[t])
    if len(a) != N: continue                       # incomplete last frame of a killed run
    a = a[np.argsort(a[:, 1])]
    assert int(a[0, 1]) == 1 and int(a[-1, 1]) == N
    frames.append((t, a[:, 2:5] * 1e3, a[:, 5:8] * 1e3))
T = np.array([f[0] for f in frames])

def ring_stats(p, v):
    core = (p[:, 2] >= CORE[0]) & (p[:, 2] < CORE[1])
    wall = core & ((p[:, 0] < SLAB) | (p[:, 0] > L - SLAB) | (p[:, 1] < SLAB) | (p[:, 1] > L - SLAB))
    inner = core & ~wall
    return -v[wall, 2].mean(), -v[inner, 2].mean(), -v[core, 2].mean()
stats = np.array([ring_stats(p, v) for _, p, v in frames])

def hermite(k, s):
    t0, p0, v0 = frames[k]; t1, p1, v1 = frames[k + 1]; h = t1 - t0
    h00 = 2 * s**3 - 3 * s**2 + 1; h10 = s**3 - 2 * s**2 + s; h01 = -2 * s**3 + 3 * s**2; h11 = s**3 - s**2
    return t0 + s * h, h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1, (1 - s) * v0 + s * v1

# ---- figure
W_PX, H_PX, DPI = 2560, 1440, 100
fig = plt.figure(figsize=(W_PX / DPI, H_PX / DPI), dpi=DPI, facecolor=SURF)
pxmm = 6.5                                       # column scale, px per mm
ch, cw = H * pxmm / H_PX, L * pxmm / W_PX        # column axes size in figure fractions
yb = 0.055
ax_w = fig.add_axes([0.045, yb, cw, ch]); ax_m = fig.add_axes([0.045 + cw + 0.022, yb, cw, ch])
ax_p = fig.add_axes([0.045 + 2 * cw + 0.050, yb, 0.085, ch])
sq = 0.50                                        # top-down slab, fraction of the figure height
ax_t = fig.add_axes([0.335, 0.455, sq * H_PX / W_PX, sq])
ax_c = fig.add_axes([0.335, 0.405, sq * H_PX / W_PX, 0.014])
ax_s = fig.add_axes([0.675, 0.455, 0.295, sq])
for ax in (ax_w, ax_m, ax_t, ax_p, ax_s):
    ax.set_facecolor(SURF)
    for s in ax.spines.values(): s.set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=13, length=4, color=GRID)

def msize(ax, d_mm):                             # scatter size (points^2) of a circle of d_mm data units
    p0 = ax.transData.transform((0, 0)); p1 = ax.transData.transform((d_mm, 0))
    return ((p1[0] - p0[0]) * 72.0 / DPI) ** 2

for ax, ttl in ((ax_w, 'at the wall'), (ax_m, 'mid-plane')):
    ax.set_xlim(0, L); ax.set_ylim(0, H); ax.set_xticks([]); ax.set_aspect('equal', adjustable='box')
    ax.set_title(ttl, color=INK, fontsize=15, pad=10)
ax_w.set_ylabel('height (mm)', color=INK2, fontsize=14); ax_m.set_yticks([])
ax_p.set_xlim(0, 0.65); ax_p.set_ylim(0, H); ax_p.set_yticks([]); ax_p.set_xticks([0, 0.15, 0.6])
ax_p.set_xticklabels(['0', '0.15', '0.6']); ax_p.set_title('solid fraction', color=INK, fontsize=15, pad=10)
ax_p.axvline(0.15, color=GRID, lw=1)
ax_t.set_xlim(0, L); ax_t.set_ylim(0, L); ax_t.set_aspect('equal'); ax_t.set_xticks([0, 5, 10, 15, 20]); ax_t.set_yticks([0, 5, 10, 15, 20])
ax_t.set_xlabel('x (mm)', color=INK2, fontsize=14); ax_t.set_ylabel('y (mm)', color=INK2, fontsize=14)
ax_t.set_title('seen from above: slab at %.0f to %.1f mm height' % (ZS0, ZS0 + SLAB), color=INK, fontsize=15, pad=10)
ax_t.plot([SLAB, L - SLAB, L - SLAB, SLAB, SLAB], [SLAB, SLAB, L - SLAB, L - SLAB, SLAB], color=INK2, lw=1, ls=(0, (4, 4)), zorder=5)
ax_t.text(SLAB + 0.3, L - SLAB - 0.3, 'interior', color=INK, fontsize=13, va='top', zorder=6, bbox=dict(facecolor=SURF, edgecolor='none', pad=3, alpha=0.85))
ax_t.text(0.3, L - 0.3, 'wall ring', color=INK, fontsize=13, va='top', zorder=6, bbox=dict(facecolor=SURF, edgecolor='none', pad=3, alpha=0.85))
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_c, orientation='horizontal', ticks=[-20, -10, 0, 10, 20])
cb.outline.set_edgecolor(GRID); ax_c.tick_params(colors=INK2, labelsize=13, color=GRID)
ax_c.set_xlabel('vertical particle velocity (mm/s): blue sinks, red rises', color=INK2, fontsize=14)
ax_s.set_xlim(0, T[-1]); ax_s.set_ylim(-2, 20); ax_s.grid(color=GRID, lw=0.6)
ax_s.set_xlabel('time (s)', color=INK2, fontsize=14)
ax_s.set_title('mean settling speed between %.0f and %.0f mm height (mm/s)' % CORE, color=INK, fontsize=15, pad=10)
ax_s.axhline(VREF, color=INK2, lw=1.5, ls=(0, (6, 4)))
ax_s.text(T[-1] * 0.985, VREF + 0.35, 'one-dimensional reference %.2f' % VREF, color=INK2, fontsize=13, ha='right')
ln_in, = ax_s.plot([], [], color=C_IN, lw=2.5); ln_wall, = ax_s.plot([], [], color=C_WALL, lw=2.5)
ln_all, = ax_s.plot([], [], color=INK, lw=2.5)
tx_in = ax_s.text(0, 0, 'interior', color=INK, fontsize=14, va='bottom'); tx_wall = ax_s.text(0, 0, 'wall ring', color=INK, fontsize=14, va='bottom')
tx_all = ax_s.text(0, 0, 'all particles', color=INK, fontsize=14, va='bottom')
sc_w = ax_w.scatter([], [], s=msize(ax_w, D), lw=0); sc_m = ax_m.scatter([], [], s=msize(ax_m, D), lw=0)
sc_t = ax_t.scatter([], [], s=msize(ax_t, D), lw=1.5, edgecolors=SURF)
ln_p, = ax_p.plot([], [], color=INK, lw=2, drawstyle='steps-mid')
refs = [ax.axhline(H0, color=INK, lw=1.5, ls=(0, (6, 4))) for ax in (ax_w, ax_m, ax_p)]
tx_ref = ax_p.text(0.63, H0 + 2, 'reference front', color=INK2, fontsize=12, ha='right')
for ax in (ax_w, ax_m):
    ax.axhspan(ZS0, ZS0 + SLAB, color=INK, alpha=0.18, lw=0)
fig.text(0.335, 0.335, 'Batch settling, 18 335 spheres of 1 mm in a 20 x 20 x 200 mm column, solid fraction 0.15', color=INK, fontsize=19)
fig.text(0.335, 0.300, 'FeatFloWer Euler-Lagrange, wall-bounded bed, fluid mesh 8 x 8 x 80. Particle paths interpolated between logged frames (0.1 s).', color=INK2, fontsize=14)
fig.text(0.335, 0.262, 'Fewer particles sit next to the walls. The denser interior sinks at about twice the reference speed;', color=INK2, fontsize=14)
fig.text(0.335, 0.236, 'the liquid returns upward along the walls and carries particles with it, so the wall-ring mean is near zero.', color=INK2, fontsize=14)
tx_time = fig.text(0.335, 0.12, '', color=INK, fontsize=44, family='monospace')
zb = np.arange(0, H + 1e-9, 2.0); zc = 0.5 * (zb[1:] + zb[:-1]); vp = math.pi / 6 * D**3

nfr = (len(frames) - 1) * nsub + 1
for i in range(max(f0, 0), min(f1, nfr - 1) + 1):
    k, j = divmod(i, nsub)
    if k == len(frames) - 1: k, j = k - 1, nsub
    t, p, v = hermite(k, j / nsub)
    col = cmap(norm(v[:, 2]))
    for sc, m in ((sc_w, p[:, 1] < SLAB), (sc_m, np.abs(p[:, 1] - L / 2) < SLAB / 2)):
        o = np.argsort(-p[m, 1]); sc.set_offsets(p[m][o][:, [0, 2]]); sc.set_facecolors(col[m][o])
    m = (p[:, 2] >= ZS0) & (p[:, 2] < ZS0 + SLAB); o = np.argsort(p[m, 2])
    sc_t.set_offsets(p[m][o][:, :2]); sc_t.set_facecolors(col[m][o])
    ln_p.set_data(np.histogram(p[:, 2], zb)[0] * vp / (L * L * 2.0), zc)
    zr = H0 - VREF * t
    for r in refs: r.set_ydata([zr, zr])
    tx_ref.set_y(zr + 2)
    n = np.searchsorted(T, t + 1e-9); ts = np.append(T[:n], t)
    for q, (ln, tx) in enumerate(((ln_wall, tx_wall), (ln_in, tx_in), (ln_all, tx_all))):
        ys = np.append(stats[:n, q], np.interp(t, T, stats[:, q])); ln.set_data(ts, ys)
        tx.set_position((min(t + 0.04, T[-1] * 0.93), ys[-1] + 0.3))
    tx_time.set_text('t = %4.2f s' % t)
    fig.savefig(os.path.join(out, 'frame_%05d.png' % i), dpi=DPI, facecolor=SURF)
print('frames', nfr, 'log frames', len(frames), 'N', N)
