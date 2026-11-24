"""z-profile per cell layer from the ASCII vtu output of a case-11 hold: gas fraction, pressure gradient,
closure gradient at the layer gas fraction. usage: layer_profile.py RUN_DIR U [nlayers_print]"""
import glob, re, sys
import numpy as np
run, U = sys.argv[1], float(sys.argv[2]); nprint = int(sys.argv[3]) if len(sys.argv) > 3 else 12
mu, rho, d = 1.8e-5, 1.0, 0.87e-3
gcl = lambda e: 150*(1-e)**2/e**3*mu*U/d**2 + 1.75*(1-e)/e**3*rho*U*U/d
files = sorted(glob.glob(run + '/_vtk/res_node_*.vtu'))
last = max(re.search(r'\.(\d+)\.vtu$', f).group(1) for f in files)
files = [f for f in files if f.endswith('.' + last + '.vtu')]
def arrays(f):
    out = {}; name = None; buf = []
    for l in open(f):
        m = re.search(r'<DataArray[^>]*Name="([^"]+)"', l)
        if m: name = m.group(1); buf = []; continue
        if '</DataArray>' in l:
            if name: out[name] = np.array(' '.join(buf).split(), dtype=float)
            name = None; continue
        if name: buf.append(l)
    return out
eps = {}; pz = {}
for f in files:
    a = arrays(f); pts = a['Points'].reshape(-1, 3); con = a['connectivity'].astype(int).reshape(-1, 8)
    zc = pts[con, 2].mean(axis=1)
    for z, e in zip(np.round(zc, 7), a['ELEpsilonF']): eps.setdefault(z, []).append(e)
    for z, p in zip(np.round(pts[:, 2], 7), a['Pressure_V']): pz.setdefault(z, []).append(p)
zl = np.array(sorted(eps)); e = np.array([np.mean(eps[z]) for z in zl])
zp = np.array(sorted(pz)); p = np.array([np.mean(pz[z]) for z in zp])
h = zp[1] - zp[0]; G = -(p[1:] - p[:-1]) / h
print('%s: step %s, %d ranks, %d layers of %.2f mm, U %.2f; p_in %.3f Pa' % (run.rstrip('/').split('/')[-1], last, len(files), len(zl), h*1e3, U, p[0]))
print(' layer z_lo-z_hi mm  eps_f    G Pa/m   Gcl(eps)  G/Gcl   dp share %')
dp = p[0] - p[-1]
for i in range(nprint):
    print(' %2d  %5.2f-%5.2f  %.4f  %7.1f  %7.1f  %.3f   %5.2f' % (i+1, zp[i]*1e3, zp[i+1]*1e3, e[i], G[i], gcl(e[i]), G[i]/gcl(e[i]), 100*G[i]*h/dp))
core = slice(2, 8)
print(' core layers 3-8: eps %.4f G %.1f ; whole dp %.3f Pa = %.4f of W/A ; if every bed layer had the core gradient over H=%.2f mm: %.4f' % (e[core].mean(), G[core].mean(), dp, dp/126.408, 23.80, G[core].mean()*0.0238/126.408))
print(' sum over layers of Gcl(eps_layer)*h (layers with eps<0.99): %.4f of W/A' % (sum(gcl(x)*h for x in e if x < 0.99)/126.408))
