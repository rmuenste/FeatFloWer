#!/usr/bin/env python3
"""Particle movie of a case 11 run (minimum fluidisation, laterally periodic bed R) from the FeatFloWer
particle logs and the bed audit. Sibling of runs/batch10/tools/batch_movie_periodic_ff.py (same surface,
inks, colour map, typography, encode recipe); one image per logged particle frame, NO interpolation.

usage: fl11_movie.py RUNDIR OUTDIR [--frames F0 F1] [--stride k] [--stats] [--check] [--no-peer]
RUNDIR   particles.xyz (t = 0, id = line number), el_particles_rank*.log (time ip f(3) t(3) p(3) v(3);
         a particle changes its rank file when it crosses a partition: all files are read and indexed
         by (time, ip)), run_slurm.log (EL_BED_AUDIT, EL_NEWTON_PAIR, EL_INLET_TABLE), _data/q2p1_param.dat.
         Read only. A run in progress can be read: a partial last line and incomplete frames are dropped.
         A rank log whose time decreases holds the frames of an earlier run (a second run in the same
         directory appends to the logs): refused (fl11_submit.sh refuses to submit into such a directory).
OUTDIR   frame_%05d.png (2560 x 1440) and the cache fl11_movie_cache.npz (positions, velocities, audit
         series; rebuilt when a log of RUNDIR changed in size or time, reused by parallel workers).
--stride k     use every k-th log frame (default 1); the images are numbered consecutively.
--frames F0 F1 render the images F0 .. F1 (inclusive, index of the strided sequence) for parallel workers.
--stats  print the frame inventory, the per-step values of the post tool and the velocity percentiles
         that the colour range was chosen from; no rendering.
--check  no rendering: the check line (frames, first / last frame time, cadence, stride, the solver end
         t0 + K dt of the evaluation's clock and the protocol end), then the displacement of the particles between consecutive images in units of the
         diameter (minimum image in x and y), per velocity step, next to max |v| dt / d of the logged
         velocities (the minimum image is only a lower bound once |v| dt exceeds L / 2).
--no-peer  leave the peer-reported per-step means (MFiX, CFDEM) out of the trace.

Definitions. Everything that is a number of the evaluation comes from fl11_post.py (imported: evaluate,
u_of_t, div, finite), nothing is re-derived:
  dp*            dp_dyn A / W of every EL_BED_AUDIT line (2.5 ms), W = -weight_z, A = area_in
  floor force    particle momentum balance: floor_bal_z of the audit line, or (older logs) the
                 EL_NEWTON_PAIR reconstruction of fl11_post.evaluate
  bed top        z99 + r of the audit line;   bed porosity  eps_bed of the audit line (bed window)
  U              fl11_post.u_of_t of the echoed inlet table; 'pending' until the echo is complete (no
                 other source: the file _data/inlet_table.txt is not read)
  step mean      fl11_post.evaluate: the last 0.3 s of each intended hold of the echoed protocol,
                 validated audit rows; shown for measured steps only
  state          fl11_post.classify of every step: fluidised (dp* > 0.95), fixed (dp* < 0.8 of the
                 measured plateau), intermediate (neither), not classified (no measured plateau),
                 unmeasurable; only the fixed set is called 'fixed bed' and shaded in the trace
  trace end      the protocol end of the evaluation (fl11_post: the last table row, or the end of a
                 finished run when later)
  running mean   centred mean over RUNMEAN = 0.1 s of the audit series (41 lines); display only
Display only (not evaluation numbers):
  cuts           slab of SLAB mm; a sphere is drawn as a disc at true scale with the radius of its largest
                 section inside the slab: r for a centre inside the slab, the chord radius
                 sqrt(r^2 - s^2) for a centre at distance s < r outside it. (The case 10 tools draw full
                 discs for centres inside the slab only; with SLAB = 0 this is a plane cut.) Discs that
                 cross a periodic plane are drawn on both sides.
  phi(z)         exact sphere-slab overlap in DZ = 1 mm bins (as the case 10 tool, peer's batch_eval.py)
  colour         vertical particle velocity, diverging map centred at 0, fixed range +-VLIM (see VLIM)
Encode: ffmpeg -framerate 30 -i OUTDIR/frame_%05d.png -vf tpad=stop_mode=clone:stop_duration=2
        -c:v libx264 -crf 18 -profile:v high -pix_fmt yuv420p -movflags +faststart out.mp4
Render job (check, 4 workers, encode; rendering runs in Slurm, not on the gateway):
        sbatch /data/warehouse17/rmuenste/runs/fluid11/movie/render_R_proto_movie.sbatch [RUN [T_END [STRIDE]]]"""
import os, sys, glob, io, json, math, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fl11_post as post

ap = argparse.ArgumentParser()
ap.add_argument('run'); ap.add_argument('out')
ap.add_argument('--frames', type=int, nargs=2, default=(0, 10**9))
ap.add_argument('--stride', type=int, default=1)
ap.add_argument('--stats', action='store_true'); ap.add_argument('--check', action='store_true')
ap.add_argument('--no-peer', action='store_true')
A = ap.parse_args()
run, out = A.run, A.out

# ---- case constants (spec 11, bed R) and view choices; style constants as in the case 10 tools
L, HCOL, D = 15.0, 100.0, 0.87             # mm: periodic box edge, column height, particle diameter
R = D / 2
ZVIEW = 40.0                               # mm, height shown (bed top <= 31.3 mm after the first 0.2 s)
SLAB = D                                   # mm, slab thickness of every cut (one diameter: single particles
                                           # and voids stay visible in the packed bed; 2 d is a closed wall)
CUT = L / 2                                # mm, lateral centre of the two column cuts
ZS0 = 12.0                                 # mm, lower edge of the top-down slab (mid height of the bed)
DZ = 1.0                                   # mm, bins of phi(z)
RUNMEAN = 0.1                              # s, running mean of the traces
# colour range, chosen from the vz percentiles of R_proto_hrz (0.1 s frames; fl11_movie.py --stats):
#   U 0.30 m/s: p1 -0.20, p5 -0.11, p99 +0.051, p99.9 +0.055 m/s;  U 0.24 m/s: -0.13, -0.07, +0.030, +0.034
# The bed rises slowly as a plug (never faster than 0.06 m/s) and falls fast. +-0.10 keeps the rise
# (up to half of the red side) apart from rest and saturates only the fastest falling particles:
# |vz| > 0.10 for 6 % of the particles at 0.30 m/s, 5 % at 0.28-0.26, 3 % at 0.24, 1 % at 0.22, none below.
VLIM = 0.10                                # m/s
SURF, INK, INK2, GRID = '#1a1a19', '#f0efec', '#c3c2b7', '#3a3a37'
C_FF, C_TOP = '#199e70', '#c98500'         # categorical: FeatFloWer pressure drop, bed top
LABEL = 'FeatFloWer + pe, Gidaspow, 9 789 particles, bed 15 x 15 mm periodic'
# peer-reported per-step means of dp* (not recomputed here), U 0.30 ... 0.04 m/s. MFiX: message of
# 2026-11-20. CFDEM: run fluid_C_gid_inlet2 with the uniform inlet (peer RUNBOOK, 2026-11-20/21); the
# earlier top-driven CFDEM run (0.948 ... 0.159) is void (gas entered as jets) and is not shown.
PEER_U = [0.30 - 0.02 * k for k in range(14)]
PEER = (('MFiX-DEM (run F_none)', 'o', [0.989, 0.999, 0.993, 0.988, 1.005, 0.997, 0.993,
                                        0.917, 0.790, 0.667, 0.547, 0.430, 0.318, 0.208]),
        ('CFDEM (run fluid_C_gid_inlet2)', 's', [0.964, 0.959, 0.950, 0.956, 0.948, 0.953, 0.936,
                                                 0.821, 0.709, 0.599, 0.492, 0.388, 0.286, 0.188]))
CACHE = os.path.join(out, 'fl11_movie_cache.npz')
CACHE_FORMAT = 'fl11_movie cache 3'          # with fl11_post.EVAL_VERSION and the deck part of the signature


# ---- data: read once, cache in OUTDIR
def signature():
    fs = sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log')))
    fs += [os.path.join(run, q) for q in ('particles.xyz', 'run_slurm.log', '_data/q2p1_param.dat', 'example.json')]
    sig = []
    for f in fs:
        if os.path.exists(f):
            st = os.stat(f); sig.append([os.path.basename(f), st.st_size, st.st_mtime_ns])
    # the cache format and the evaluation version: a cache of an older tool is never reused
    return json.dumps({'format': CACHE_FORMAT, 'evaluation': post.EVAL_VERSION, 'files': sig})


def read_log(f):
    """rows of one particle log: tick (microseconds), ip, position (mm), velocity (m/s); refuses a log whose
    time decreases (the frames of an earlier run in the same directory, followed by those of a new one)"""
    with open(f, 'rb') as fh:
        data = fh.read()
    data = data[:data.rfind(b'\n') + 1]                # a running job may leave a partial last line
    if data.count(b'\n') < 2:
        return None
    a = np.loadtxt(io.BytesIO(data), comments='#', usecols=(0, 1, 8, 9, 10, 11, 12, 13), ndmin=2)
    back = np.flatnonzero(np.diff(a[:, 0]) < 0)
    if len(back):
        j = back[0]
        sys.exit('refused: %s: the time decreases from %.6g to %.6g s (data row %d): the file holds particle frames '
                 'of an earlier run (a second run in the same directory appends to the rank logs). Render a run '
                 'directory with the logs of one run only (stage a fresh run directory; fl11_submit.sh refuses to '
                 'submit into a directory with rank logs).' % (f, a[j, 0], a[j + 1, 0], j + 2))
    return (np.round(a[:, 0] * 1e6).astype(np.int64), a[:, 1].astype(np.int64),
            (a[:, 2:5] * 1e3).astype(np.float32), a[:, 5:8].astype(np.float32))


STATES = ('fluidised', 'fixed', 'intermediate', 'unclassified', 'unmeasurable')   # fl11_post.classify
LABELS = {'fluidised': 'fluidised', 'fixed': 'fixed bed', 'intermediate': 'intermediate',
          'unclassified': 'not classified', 'unmeasurable': 'unmeasurable'}


def build_cache(sig):
    xyz = np.loadtxt(os.path.join(run, 'particles.xyz')) * 1e3
    n = len(xyz)
    parts = [q for q in (read_log(f) for f in sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log')))) if q]
    if parts:
        tick = np.concatenate([q[0] for q in parts]); ip = np.concatenate([q[1] for q in parts])
        p = np.vstack([q[2] for q in parts]); v = np.vstack([q[3] for q in parts])
        del parts
        assert ip.min() >= 1 and ip.max() <= n, 'particle id outside 1 .. %d' % n
        ut, fi = np.unique(tick, return_inverse=True)
        cnt = np.bincount(fi, minlength=len(ut))
        seen = np.zeros((len(ut), n), bool); seen[fi, ip - 1] = True
        ok = (cnt == n) & seen.all(axis=1)             # every particle exactly once
        del seen
        bad = np.flatnonzero(~ok)
        if len(bad):
            tail = len(bad) and bad[0] == len(ut) - len(bad)
            print('note: %d incomplete log frame(s) dropped, t = %s s%s' % (
                len(bad), ' '.join('%.4f' % (ut[b] * 1e-6) for b in bad[:6]),
                ' (end of a run in progress)' if tail else '  WARNING: not only at the end'), file=sys.stderr)
        new = np.cumsum(ok) - 1; keep = ok[fi]
        nf = int(ok.sum())
        pos = np.empty((nf + 1, n, 3), np.float32); vel = np.zeros((nf + 1, n, 3), np.float32)
        pos[new[fi[keep]] + 1, ip[keep] - 1] = p[keep]; vel[new[fi[keep]] + 1, ip[keep] - 1] = v[keep]
        T = np.concatenate([[0.0], ut[ok] * 1e-6])
    else:
        pos = np.empty((1, n, 3), np.float32); vel = np.zeros((1, n, 3), np.float32); T = np.array([0.0])
    pos[0] = xyz

    # audit series and step table: fl11_post.evaluate (validated audit rows, the holds of the echoed protocol,
    # the per-step evaluation with its guarded ratios and coverage; nothing re-derived here)
    W = area = float('nan'); aud = np.zeros((0, 8)); srow = np.zeros((0, 10)); table = []; echo = False
    horizon = solver_end = float('nan')
    if os.path.exists(os.path.join(run, 'run_slurm.log')):
        ev = post.evaluate(run, D * 1e-3)
        echo = ev['source'] == 'echo' and not ev['protocol_errors']
        table = list(ev['knots']) if echo else []       # the inlet display waits for the complete echo
        horizon, solver_end = (ev['horizon'] if echo else float('nan')), ev['clock'].t_last
        rows = ev['clock'].sorted
        area, W = ev['area'], ev['W']
        if rows:
            # display columns: t dp* floor/W bed top (absolute z, mm; the frame coordinates) eps zmax n U, each
            # row normalised with its own weight and area
            aud = np.array([(r['t'], post.div(r['dp_dyn'], post.div(-r['wz'], r['a_in'])), post.div(r['fbal'], -r['wz']),
                             (r['z99'] + R * 1e-3) * 1e3, r['eps'], r['zmax'] * 1e3, r['n'],
                             r['U'] if echo else float('nan')) for r in rows])
        srow = np.array([(x['t0'], x['t1'], x['U'], x['dps'], x['stds'], x['H'] * 1e3, x['eps'], x['fl'],
                          float(x['measured']), STATES.index(c)) for x, c in zip(ev['steps'], ev['classes'])]).reshape(-1, 10)
    os.makedirs(out, exist_ok=True)
    tmp = os.path.join(out, '.fl11_movie_cache.%d.npz' % os.getpid())
    np.savez(tmp, sig=sig, T=T, pos=pos, vel=vel, aud=aud, srow=srow, table=np.array(table).reshape(-1, 2),
             W=W, area=area, echo=echo, horizon=horizon, solver_end=solver_end)
    os.replace(tmp, CACHE)


sig = signature()
fresh = False
if os.path.exists(CACHE):
    with np.load(CACHE) as z:
        fresh = str(z['sig']) == sig
    if not fresh:
        print('note: cache %s not reused (older cache format or evaluation version, or changed inputs): rebuilt' % CACHE,
              file=sys.stderr)
if not fresh:
    build_cache(sig)
with np.load(CACHE) as z:
    T, POS, VEL, AUD, SROW, TABLE = z['T'], z['pos'], z['vel'], z['aud'], z['srow'], z['table']
    W, AREA, ECHO = float(z['W']), float(z['area']), bool(z['echo'])
    HORIZON, SOLVER_END = float(z['horizon']), float(z['solver_end'])
N = POS.shape[1]
table = [tuple(r) for r in TABLE]
TA = AUD[:, 0]                                         # audit columns: t dp* floor/W H eps zmax n U
sel = np.arange(0, len(T), A.stride)                   # log frames of the rendered sequence
nfr = len(sel)
TEND = HORIZON if post.finite(HORIZON) else T[-1]     # protocol end of the evaluation (frames only without the echo)
STATE = [STATES[int(r[9])] for r in SROW]              # display state of every step (fl11_post.classify)
FIXSET = [r for r, c in zip(SROW, STATE) if c == 'fixed']
FLUSET = [r for r, c in zip(SROW, STATE) if c == 'fluidised']


def audit_at(t):
    """audit line of the frame time (the logs share the step clock), None when there is none"""
    k = np.searchsorted(TA, t - 1e-7)
    return AUD[k] if k < len(TA) and abs(TA[k] - t) < 1e-6 else None


def step_at(t):
    for r in SROW:
        if r[0] - 1e-9 <= t <= r[1] + 1e-9: return r
    return None


def runmean(y):
    """centred mean over RUNMEAN; returns the times and the means of the full windows"""
    if len(TA) < 3: return TA[:0], TA[:0]
    n = int(round(RUNMEAN / (TA[1] - TA[0]))) + 1
    if len(TA) < n: return TA[:0], TA[:0]
    return TA[n // 2:len(TA) - (n - 1 - n // 2)], np.convolve(y, np.ones(n) / n, 'valid')


def groups():
    """protocol intervals for the tables: preparation, then the velocity steps"""
    g = []
    if len(SROW):
        if SROW[0, 0] > 0: g.append((0.0, SROW[0, 0], SROW[0, 2], 'prep'))
        g += [(r[0], r[1], r[2], '%4d' % k) for k, r in enumerate(SROW)]
    else:
        g.append((0.0, T[-1], float('nan'), ' all'))
    return g


if A.stats:
    print('%s: N %d, %d frames (t = 0 and %d log frames), t %.4f .. %.4f s, cadence %s s' % (
        os.path.basename(os.path.abspath(run)), N, len(T), len(T) - 1, T[0], T[-1],
        ('%.4f' % np.median(np.diff(T))) if len(T) > 1 else '-'))
    print('audit lines %d%s, W %.6e N, A %.6e m2, W/A %.3f Pa; inlet table %s; solver end %.4f s, protocol end %s' % (
        len(TA), (' (t %.4f .. %.4f s)' % (TA[0], TA[-1])) if len(TA) else '', W, AREA, post.div(W, AREA),
        'from the log echo' if ECHO else 'pending (echo absent or incomplete)', SOLVER_END,
        ('%.4f s' % HORIZON) if post.finite(HORIZON) else 'unknown'))
    print('\nper step (fl11_post.evaluate over the last 0.3 s, state of fl11_post.classify):\n step  U [m/s]  window      dp*      std   H [mm]  eps_bed  floor/W  dp*+floor  state')
    for k, r in enumerate(SROW):
        print(' %3d  %6.3f  %5.2f-%5.2f  %7.4f  %7.4f  %6.2f  %7.4f  %7.4f  %8.4f   %s' % (
            k, r[2], r[1] - 0.3, r[1], r[3], r[4], r[5], r[6], r[7], r[3] + r[7],
            LABELS[STATE[k]]))
    print('\nvertical particle velocity of the log frames (m/s), colour range +-%.2f:\n step  U [m/s]  frames     p0.1      p1      p5     p50     p95     p99   p99.9  max|v|  share |vz| > range  highest centre [mm]' % VLIM)
    for t0, t1, u, name in groups():
        m = (T > t0 + 1e-9) & (T <= t1 + 1e-9)
        if not m.any(): continue
        vz = VEL[m][:, :, 2].ravel()
        print(' %s  %6.3f  %6d  ' % (name, u, m.sum()) + ' '.join('%7.3f' % q for q in np.percentile(vz, [0.1, 1, 5, 50, 95, 99, 99.9])) +
              '  %6.3f  %8.4f  %18.2f' % (np.sqrt((VEL[m].astype(float)**2).sum(-1)).max(), (np.abs(vz) > VLIM).mean(), POS[m][:, :, 2].max()))
    sys.exit(0)

if A.check:
    print('check: frames %d t_first %.4f t_last %.4f cadence %.4f stride %d solver_end %.4f protocol_end %.4f' % (
        nfr, T[sel[0]], T[sel[-1]], np.median(np.diff(T[sel])) if nfr > 1 else 0.0, A.stride, SOLVER_END, HORIZON))
    if nfr < 2: sys.exit(0)
    dmax = np.zeros(nfr - 1); d99 = np.zeros(nfr - 1); d50 = np.zeros(nfr - 1); vdt = np.zeros(nfr - 1); who = np.zeros(nfr - 1, int)
    for j in range(nfr - 1):
        a, b = sel[j], sel[j + 1]
        dlt = (POS[b] - POS[a]).astype(float); dlt[:, :2] -= L * np.round(dlt[:, :2] / L)   # minimum image
        s = np.sqrt((dlt**2).sum(1)) / D
        dmax[j], d99[j], d50[j], who[j] = s.max(), np.percentile(s, 99), np.median(s), s.argmax() + 1
        sp = np.maximum(np.sqrt((VEL[a].astype(float)**2).sum(1)), np.sqrt((VEL[b].astype(float)**2).sum(1)))
        vdt[j] = sp.max() * 1e3 * (T[b] - T[a]) / D
    te = T[sel[1:]]
    print('check: displacement between consecutive images in units of d = %.2f mm (3-D, minimum image in x and y; p99 and median over the particles; each column: worst image pair of the interval)' % D)
    print('  step  U [m/s]   t [s]        pairs   max [d]   p99 [d]  median [d]  max |v| dt [d]')
    for t0, t1, u, name in groups():
        m = (te > t0 + 1e-9) & (te <= t1 + 1e-9)
        if m.any():
            print('  %s  %6.3f  %5.2f-%5.2f  %6d  %8.3f  %8.3f  %10.4f  %14.3f' % (name, u, t0, t1, m.sum(), dmax[m].max(), d99[m].max(), d50[m].max(), vdt[m].max()))
    j = dmax.argmax()
    print('check: largest displacement %.3f d (%.3f mm) between t = %.4f and %.4f s, particle %d' % (dmax[j], dmax[j] * D, T[sel[j]], te[j], who[j]))
    k = np.flatnonzero(te > 0.2 + 1e-9)
    if len(k):
        j = k[dmax[k].argmax()]
        print('check: after the collapse of the start bed (t > 0.2 s): largest %.3f d between t = %.4f and %.4f s' % (dmax[j], T[sel[j]], te[j]))
    if vdt.max() * D > L / 2:
        print('check: max |v| dt reaches %.2f d > L / 2 = %.2f d: the minimum-image displacement is only a lower bound there' % (vdt.max(), L / 2 / D))
    sys.exit(0)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
os.makedirs(out, exist_ok=True)
cmap = LinearSegmentedColormap.from_list('vz', ['#4d9bff', '#2f62a8', '#77776f', '#a8523d', '#ff7a5c'])
norm = Normalize(-VLIM, VLIM)

# ---- figure (pixel geometry; surface and typography of the case 10 movies)
W_PX, H_PX, DPI = 2560, 1440, 100
fig = plt.figure(figsize=(W_PX / DPI, H_PX / DPI), dpi=DPI, facecolor=SURF)
def axpx(x0, y0, w, h): return fig.add_axes([x0 / W_PX, y0 / H_PX, w / W_PX, h / H_PX])
def fx(x): return x / W_PX
def fy(y): return y / H_PX
PXMM = 30.0                                      # scale of the column cuts, px per mm
TD = 450.0                                       # edge of the top-down slab, px
PXT = TD / L                                      # scale of the top-down slab
YB, CW, CH = 125, L * PXMM, ZVIEW * PXMM
ax_w = axpx(105, YB, CW, CH); ax_m = axpx(105 + CW + 40, YB, CW, CH)
ax_p = axpx(105 + 2 * CW + 85, YB, 180, CH)
XR = 1385                                        # left edge of the right block
ax_t = axpx(XR, YB + CH - TD, TD, TD)
XT = XR + TD + 60                                # text column
ax_c = axpx(XT, 815, W_PX - 50 - XT, 14)
TW = 1055
ax_s = axpx(XR, YB, TW, 480); ax_u = axpx(XR, YB + 500, TW, 85)
ax_h = ax_s.twinx()
for ax in (ax_w, ax_m, ax_t, ax_p, ax_s, ax_u, ax_h):
    ax.set_facecolor(SURF)
    for s in ax.spines.values(): s.set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=13, length=4, color=GRID)
ax_h.patch.set_visible(False)

def msize(pxmm, d_mm):                           # scatter size (points^2) of a disc of diameter d_mm
    return (np.asarray(d_mm) * pxmm * 72.0 / DPI) ** 2
DASH = (0, (6, 4)); DOT = (0, (1.5, 2.5))
box = dict(facecolor=SURF, edgecolor='none', pad=3, alpha=0.85)

# column cuts: periodic sides dashed, floor solid
for ax, ttl, lab in ((ax_w, 'x-z cut, %.2f mm slab' % SLAB, 'x (mm)'), (ax_m, 'y-z cut, %.2f mm slab' % SLAB, 'y (mm)')):
    ax.set_xlim(0, L); ax.set_ylim(0, ZVIEW); ax.set_xticks([0, 5, 10, 15]); ax.set_yticks(np.arange(0, ZVIEW + 1, 5))
    ax.set_title(ttl, color=INK, fontsize=15, pad=10); ax.set_xlabel(lab, color=INK2, fontsize=14)
    for side in ('left', 'right'):
        ax.spines[side].set_linestyle((0, (4, 4))); ax.spines[side].set_color(INK2); ax.spines[side].set_linewidth(1.0)
    ax.spines['bottom'].set_color(INK2); ax.spines['bottom'].set_linewidth(1.6); ax.spines['top'].set_visible(False)
    for q in (ZS0, ZS0 + SLAB): ax.axhline(q, color=INK2, lw=1, ls=DOT, zorder=5)
ax_w.set_ylabel('height (mm)', color=INK2, fontsize=14); ax_m.set_yticklabels([])
ax_p.set_xlim(0, 0.72); ax_p.set_ylim(0, ZVIEW); ax_p.set_yticks(np.arange(0, ZVIEW + 1, 5)); ax_p.set_yticklabels([])
ax_p.set_xticks([0, 0.3, 0.6]); ax_p.set_xticklabels(['0', '0.3', '0.6'])
ax_p.set_title('solid fraction', color=INK, fontsize=15, pad=10); ax_p.set_xlabel('1 mm bins', color=INK2, fontsize=14)
for q in (0.3, 0.6): ax_p.axvline(q, color=GRID, lw=1)

# top-down slab
ax_t.set_xlim(0, L); ax_t.set_ylim(0, L); ax_t.set_xticks([0, 5, 10, 15]); ax_t.set_yticks([0, 5, 10, 15])
ax_t.set_xlabel('x (mm)', color=INK2, fontsize=14); ax_t.set_ylabel('y (mm)', color=INK2, fontsize=14)
ax_t.set_title('seen from above: slab at %.1f to %.1f mm height' % (ZS0, ZS0 + SLAB), color=INK, fontsize=15, pad=10)
for s in ax_t.spines.values(): s.set_linestyle((0, (4, 4))); s.set_color(INK2); s.set_linewidth(1.0)
c0, c1 = CUT - SLAB / 2, CUT + SLAB / 2          # the two column cuts, drawn in the top-down view
for q in (c0, c1):
    ax_t.axhline(q, color=INK2, lw=1, ls=DOT, zorder=5); ax_t.axvline(q, color=INK2, lw=1, ls=DOT, zorder=5)
ax_t.text(0.25, c1 + 0.2, 'x-z cut', color=INK, fontsize=13, va='bottom', zorder=6, bbox=box)
ax_t.text(c1 + 0.2, L - 0.25, 'y-z cut', color=INK, fontsize=13, va='top', zorder=6, bbox=box)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_c, orientation='horizontal',
                  ticks=[-VLIM, -VLIM / 2, 0, VLIM / 2, VLIM])
cb.outline.set_edgecolor(GRID); ax_c.tick_params(colors=INK2, labelsize=13, color=GRID)
ax_c.set_xlabel('vertical particle velocity (m/s): blue falls, red rises', color=INK2, fontsize=14)

# trace over the protocol: inlet velocity (strip), dp* and dp* + floor (left axis), bed top (right axis)
YS, YH = (0.0, 1.6), (-3.0, 33.0)                # bed top 24 .. 32 mm in the band above the pressure traces
for ax in (ax_s, ax_u):
    ax.set_xlim(0, TEND); ax.set_xticks(np.arange(0, TEND + 1e-9, 0.5)); ax.grid(color=GRID, lw=0.6)
ax_u.set_xticklabels([]); ax_u.set_ylim(0, 0.34); ax_u.set_yticks([0, 0.1, 0.2, 0.3])
ax_u.set_ylabel('U (m/s)', color=INK2, fontsize=14)
ax_u.set_title('inlet velocity U, pressure drop dp* = dp A / W and bed top over the protocol', color=INK, fontsize=15, pad=10)
if table: ax_u.plot([q[0] for q in table], [q[1] for q in table], color=INK, lw=2)
ax_s.set_ylim(*YS); ax_s.set_yticks(np.arange(0, 1.5 + 1e-9, 0.25))
ax_s.set_xlabel('time (s)', color=INK2, fontsize=14); ax_s.set_ylabel('dp* = dp A / W', color=INK2, fontsize=14)
ax_s.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '%g' % v))
ax_h.set_xlim(0, TEND); ax_h.set_ylim(*YH); ax_h.set_yticks([24, 26, 28, 30, 32])
ax_h.set_ylabel('bed top (mm)', color=C_TOP, fontsize=14); ax_h.yaxis.set_label_coords(1.045, (28 - YH[0]) / (YH[1] - YH[0]))
ax_h.tick_params(axis='y', colors=C_TOP)
ax_s.axhline(1.0, color=INK2, lw=1.5, ls=DASH, zorder=2)
for r in FIXSET:                                 # the fixed set of fl11_post.classify, step by step
    for ax in (ax_s, ax_u):
        ax.axvspan(r[0], r[1], color=INK, alpha=0.05, lw=0, zorder=0)
if FIXSET and FLUSET and FLUSET[-1][1] <= FIXSET[0][0] + 1e-9:   # one transition: label both sides
    t_fix = FIXSET[0][0]
    for ax in (ax_s, ax_u): ax.axvline(t_fix, color=INK2, lw=1, zorder=1)
    ax_u.text(t_fix - 0.08, 0.035, 'fluidised', color=INK, fontsize=14, ha='right', va='bottom')
    ax_u.text(t_fix + 0.08, 0.305, 'fixed bed', color=INK, fontsize=14, ha='left', va='top')
if not table:
    ax_u.text(0.02, 0.5, 'inlet velocity: table echo pending', color=INK2, fontsize=14, transform=ax_u.transAxes, va='center')
handles = []
if len(TA):
    handles.append(ax_s.plot(TA, AUD[:, 1], color=C_FF, lw=1, alpha=0.40, zorder=3, label='FeatFloWer: pressure drop dp*, every 2.5 ms')[0])
    tm, ym = runmean(AUD[:, 1])
    handles.append(ax_s.plot(tm, ym, color=C_FF, lw=2.5, zorder=5, label='FeatFloWer: dp* averaged over %g s' % RUNMEAN)[0])
    tm, ym = runmean(AUD[:, 1] + AUD[:, 2])
    handles.append(ax_s.plot(tm, ym, color=INK, lw=1.6, zorder=4, label='FeatFloWer: dp* + force on the floor / W, averaged over %g s' % RUNMEAN)[0])
    handles.append(ax_h.plot(TA, AUD[:, 3], color=C_TOP, lw=2, zorder=4, label='FeatFloWer: bed top (right axis)')[0])
if not A.no_peer and len(SROW):
    for name, mk, val in PEER:
        xs, ys = [], []
        for r in SROW:
            k = [i for i, u in enumerate(PEER_U) if abs(u - r[2]) < 1e-6]
            if k: xs.append(r[1] - 0.15); ys.append(val[k[0]])
        handles.append(ax_s.plot(xs, ys, mk, ms=7, mfc=SURF, mec=INK, mew=1.3, zorder=6, ls='none',
                                 label=name + ': peer-reported step mean')[0])
leg = ax_s.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.012, 0.025), fontsize=13, frameon=True,
                  facecolor=SURF, edgecolor='none', framealpha=1.0, labelcolor=INK, handlelength=2.2, borderpad=0.4, labelspacing=0.42)
leg.set_zorder(9)                                # above the time cursor
ax_s.text(TEND * 0.988, 1.0 + 0.025, 'dp* = 1: buoyant bed weight per area, W / A = %.2f Pa' % (W / AREA) if len(TA) else '', color=INK2, fontsize=13,
          va='bottom', ha='right', zorder=6)
cur = [ax.axvline(0, color=INK, lw=1.3, zorder=8) for ax in (ax_s, ax_u)]

# text block
YT = YB + CH                                     # top of the panels
fig.text(fx(XT), fy(YT - 29), 'Minimum fluidisation of a gas-solid bed', color=INK, fontsize=19)
cad = np.median(np.diff(T[sel])) if nfr > 1 else 0.0
for q, s in enumerate((LABEL[:LABEL.index(' bed 15')], LABEL[LABEL.index('bed 15'):],
                       'Inlet velocity lowered by 0.02 m/s every 0.5 s.',
                       'One image per logged particle frame (%g ms), no interpolation.' % (cad * 1e3),
                       'Disc size: largest section of the sphere inside the cut slab.')):
    fig.text(fx(XT), fy(YT - 63 - 26 * q), s, color=INK2, fontsize=13)
tx_time = fig.text(fx(XT), fy(YT - 245), '', color=INK, fontsize=44, family='monospace')
tx_reg = fig.text(fx(XT + 425), fy(YT - 238), '', color=INK2, fontsize=15)
XV1, XV2, XUN = XT + 340, XT + 470, XT + 482
fig.text(fx(XV1), fy(YT - 285), 'now', color=INK2, fontsize=13, ha='right'); fig.text(fx(XV2), fy(YT - 285), 'step mean', color=INK2, fontsize=13, ha='right')
ROWS = (('inlet velocity U', 'm/s'), ('pressure drop dp*', ''), ('bed top z99 + r', 'mm'), ('bed porosity', ''))
tx_now, tx_mean = [], []
for q, (name, unit) in enumerate(ROWS):
    y = fy(YT - 317 - 31 * q)
    fig.text(fx(XT), y, name, color=INK, fontsize=15); fig.text(fx(XUN), y, unit, color=INK2, fontsize=13)
    tx_now.append(fig.text(fx(XV1), y, '', color=INK, fontsize=15, ha='right', family='monospace'))
    tx_mean.append(fig.text(fx(XV2), y, '', color=INK2, fontsize=15, ha='right', family='monospace'))
fig.text(fx(XT), fy(YT - 317 - 31 * 4), 'now: bed audit line of the frame; step mean: last 0.3 s of the step', color=INK2, fontsize=12)

sc_w = ax_w.scatter([], [], s=1, lw=0.6, edgecolors=SURF); sc_m = ax_m.scatter([], [], s=1, lw=0.6, edgecolors=SURF)
sc_t = ax_t.scatter([], [], s=1, lw=0.8, edgecolors=SURF)
ln_p, = ax_p.plot([], [], color=INK, lw=2, drawstyle='steps-mid')
tops = [ax.axhline(0, color=C_TOP, lw=1.8, ls=DASH, zorder=6) for ax in (ax_w, ax_m, ax_p)]
tx_top = ax_p.text(0.70, 0, 'bed top', color=C_TOP, fontsize=13, ha='right', va='bottom')
tx_clip = ax_w.text(0.3, ZVIEW - 0.4, '', color=INK, fontsize=13, va='top', zorder=7, bbox=box)

nb = int(round(HCOL / DZ)); zc = (np.arange(nb) + 0.5) * DZ
assert D <= DZ
def cap(h):                                        # volume of a sphere below height h above its bottom
    h = np.clip(h, 0.0, D); return math.pi * h * h * (3 * R - h) / 3
def profile(z):                                    # phi(z), exact sphere-slab overlap (D <= DZ: two bins at most)
    zb = z - R; b0 = np.floor(zb / DZ).astype(int)
    v0 = cap((b0 + 1) * DZ - zb) - cap(b0 * DZ - zb); v1 = math.pi * D**3 / 6 - cap((b0 + 1) * DZ - zb)
    v = np.zeros(nb + 2)
    np.add.at(v, np.clip(b0, -1, nb) + 1, v0); np.add.at(v, np.clip(b0 + 1, -1, nb) + 1, v1)
    return v[1:-1] / (L * L * DZ)

def cut(p, col, dep, lo, hi, lat, order, periodic):
    """discs of the slab lo <= p[dep] <= hi: centres (lat columns), radius of the largest section inside the
    slab, colours; drawn in the given depth order; periodic images across the listed lateral columns"""
    s = p[:, dep]; dist = np.maximum(np.maximum(lo - s, s - hi), 0.0)
    m = dist < R
    o = np.argsort(order * s[m])
    xy = p[m][o][:, lat].astype(float); rc = np.sqrt(R * R - dist[m][o]**2); cc = col[m][o]
    for a in periodic:
        add = []                                       # images of this axis include those of the axes before
        for sgn, g in ((L, xy[:, a] < rc), (-L, xy[:, a] > L - rc)):
            q = xy[g].copy(); q[:, a] += sgn; add.append((q, rc[g], cc[g]))
        xy = np.vstack([xy] + [q[0] for q in add]); rc = np.concatenate([rc] + [q[1] for q in add])
        cc = np.vstack([cc] + [q[2] for q in add])
    return xy, rc, cc

for j in range(max(A.frames[0], 0), min(A.frames[1], nfr - 1) + 1):
    k = sel[j]; t = T[k]; p = POS[k].astype(float); vz = VEL[k][:, 2]
    col = cmap(norm(vz))
    # x-z cut: slab in y around CUT, seen along +y (nearest drawn last); y-z cut: slab in x, seen along +x
    for sc, dep, lat in ((sc_w, 1, [0, 2]), (sc_m, 0, [1, 2])):
        xy, rc, cc = cut(p, col, dep, c0, c1, lat, -1.0, (0,))
        sc.set_offsets(xy); sc.set_sizes(msize(PXMM, 2 * rc)); sc.set_facecolors(cc)
    xy, rc, cc = cut(p, col, 2, ZS0, ZS0 + SLAB, [0, 1], 1.0, (0, 1))
    sc_t.set_offsets(xy); sc_t.set_sizes(msize(PXT, 2 * rc)); sc_t.set_facecolors(cc)
    ln_p.set_data(profile(p[:, 2]), zc)
    a = audit_at(t); r = step_at(t)
    u = post.u_of_t(table, t) if table else float('nan')
    tx_time.set_text('t = %5.3f s' % t)
    tx_now[0].set_text('%.3f' % u if table else 'pending')
    for q, (i, fmt) in enumerate(((1, '%.3f'), (3, '%.2f'), (4, '%.3f'))):
        tx_now[q + 1].set_text(fmt % a[i] if a is not None else '')
    full = r is not None and r[8] > 0
    for q, (i, fmt) in enumerate(((2, '%.3f'), (3, '%.3f'), (5, '%.2f'), (6, '%.3f'))):
        tx_mean[q].set_text(fmt % r[i] if full else '')
    tx_reg.set_text(LABELS[STATES[int(r[9])]] if r is not None else ('preparation' if len(SROW) and t < SROW[0, 0] else ''))
    top = a[3] if a is not None else None
    for q in tops:
        q.set_visible(top is not None)
        if top is not None: q.set_ydata([top, top])
    tx_top.set_visible(top is not None and top < ZVIEW - 1.5)
    if top is not None: tx_top.set_y(top + 0.35)
    zhi = p[:, 2].max() + R
    tx_clip.set_text('start bed (loose, top at %.0f mm) collapsing' % zhi if zhi > ZVIEW else '')
    for q in cur: q.set_xdata([t, t])
    fig.savefig(os.path.join(out, 'frame_%05d.png' % j), dpi=DPI, facecolor=SURF)
print('frames', nfr, 'log frames', len(T), 'stride', A.stride, 'N', N, 'rendered', max(A.frames[0], 0), min(A.frames[1], nfr - 1))
