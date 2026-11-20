#!/usr/bin/env python3
"""Decompose the fixed-bed core-gradient deficit of a case-11 run at one output frame.
usage: core_deficit.py RUNDIR FRAMESTEP TIME U   (example: R_proto_hrz 29999 7.5 0.04)
Reads _vtk/res_node_*.<FRAMESTEP>.vtu (ELEpsilonF per element, Velocity at the vertices), the particle
logs at TIME (position, PE force = particle_total), the EL_BED_AUDIT line nearest to TIME (window,
eps_bed, gradz_win). Replicates the kernel sample of the code: deen_poly kernel of radius 2.5 d, 27-point
Gauss quadrature per element (the code's CB3H(9)), periodic images in x and y, eps_f piecewise constant
per element (exact), velocity trilinear from the vertices (the code uses the Q2 field: approximation).
Model A, particles at rest: F_tot = F_d + F_p, grad_meas = sum(frac F_tot)/V_win."""
import sys, glob, math, re
import xml.etree.ElementTree as ET
R, FR, T, U = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
d = 0.87e-3; r = d / 2; mu = 1.8e-5; rho = 1.0; rhop = 860.0; g = 9.81; L = 0.015
vp = math.pi * d**3 / 6; width = 2.5 * d

def arr(piece, name):
    for a in piece.iter('DataArray'):
        if a.get('Name') == name:
            return [float(x) for x in a.text.split()]
    raise KeyError(name)

elems = {}; vel = {}
h = 2.5e-3
for f in glob.glob(f'{R}/_vtk/res_node_*.{FR}.vtu'):
    piece = ET.parse(f).getroot().find('UnstructuredGrid').find('Piece')
    pts = arr(piece, 'Points'); v = arr(piece, 'Velocity'); con = [int(x) for x in arr(piece, 'connectivity')]
    eps = arr(piece, 'ELEpsilonF')
    P = [tuple(pts[3*i:3*i+3]) for i in range(len(pts)//3)]
    for i, p in enumerate(P):
        vel[(round(p[0]/h), round(p[1]/h), round(p[2]/h))] = v[3*i:3*i+3]
    for e in range(len(eps)):
        c = [P[j] for j in con[8*e:8*e+8]]
        lo = tuple(min(q[k] for q in c) for k in range(3))
        elems[(round(lo[0]/h), round(lo[1]/h), round(lo[2]/h))] = eps[e]
nx = round(L/h); nz = max(k[2] for k in elems) + 1
print(f'{len(elems)} elements ({nx} x {nx} x {nz}), {len(vel)} vertices')

def uvert(ix, iy, iz):
    # periodic vertex index in x and y (the vertex at x = L carries the same value as x = 0)
    return vel.get((ix % nx, iy % nx, iz)) or vel.get((nx if ix % nx == 0 else ix % nx, nx if iy % nx == 0 else iy % nx, iz))


def u_at(x, y, z):
    ix = math.floor(x/h); iy = math.floor(y/h); iz = min(max(math.floor(z/h), 0), nz-1)
    fx = x/h - ix; fy = y/h - iy; fz = z/h - iz
    out = [0.0, 0.0, 0.0]
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = (fx if a else 1-fx)*(fy if b else 1-fy)*(fz if c else 1-fz)
                vv = uvert(ix+a, iy+b, iz+c)
                for k in range(3): out[k] += w*vv[k]
    return out

gp = [-math.sqrt(0.6), 0.0, math.sqrt(0.6)]; gw = [5/9, 8/9, 5/9]
def kernel(dist):
    q = dist/width
    return 0.0 if q >= 1 else (1+3*q)*(1-q)**3

# particles
parts = []
for f in glob.glob(f'{R}/el_particles_rank*.log'):
    for l in open(f):
        if l[0] == '#': continue
        c = l.split()
        if abs(float(c[0]) - T) < 1e-9:
            parts.append(dict(x=float(c[8]), y=float(c[9]), z=float(c[10]), fz=float(c[4])))
# audit line
best = None
for l in open(f'{R}/run_slurm.log', errors='replace'):
    if l.startswith('EL_BED_AUDIT'):
        c = l.split(); t = float(c[2])
        if best is None or abs(t-T) < abs(best[0]-T):
            best = (t, float(c[26]), float(c[28]), float(c[29]), float(c[31]), float(c[14]), float(c[10]))
ta, eps_bed, zlo, zhi, gw_meas, dpdyn, area = best
print(f'{len(parts)} particles at t {T}; audit t {ta}: eps_bed {eps_bed:.5f} window {zlo*1e3:.3f}-{zhi*1e3:.3f} mm gradz_win {gw_meas:.1f} Pa/m')

def slab(zc, a, b):
    lo = max(a-zc, -r); hi = min(b-zc, r)
    return 0.0 if hi <= lo else math.pi*(r*r*(hi-lo) - (hi**3-lo**3)/3)

def drag_B(eps, s):
    if eps <= 0.8:
        return vp*(150*(1-eps)*mu/(eps*d*d) + 1.75*rho*s/d)
    re_ = rho*eps*s*d/mu
    cds = 24*mu/(rho*eps*d)*(1+0.15*re_**0.687) if re_ < 1000 else 0.44*s
    return vp*0.75*cds*rho*eps*eps**-2.65/d

def grad_cl(eps, Usup):
    s = Usup/eps
    return (1-eps)/vp*drag_B(eps, s)*s/eps

wbg = (rhop-rho)*vp*g
rows = []
for p in parts:
    fr = slab(p['z'], zlo, zhi)/vp
    num = 0.0; den = 0.0; ux = [0.0, 0.0, 0.0]
    ic = math.floor(p['x']/h); jc = math.floor(p['y']/h); kc = math.floor(p['z']/h)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                k = kc + dk
                if k < 0 or k >= nz: continue
                i = ic + di; j = jc + dj
                e = elems[(i % nx, j % nx, k)]
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            qx = (i + 0.5 + 0.5*gp[a])*h; qy = (j + 0.5 + 0.5*gp[b])*h; qz = (k + 0.5 + 0.5*gp[c])*h
                            w = gw[a]*gw[b]*gw[c]*kernel(math.dist((qx, qy, qz), (p['x'], p['y'], p['z'])))
                            if w <= 0: continue
                            num += w*e; den += w
                            uu = u_at(qx % L, qy % L, qz)
                            for m in range(3): ux[m] += w*uu[m]
    eps_p = min(1.0, max(0.2, num/den)); u_p = [x/den for x in ux]
    s = math.sqrt(sum(x*x for x in u_p))
    fd = drag_B(eps_p, s)*u_p[2]
    rows.append(dict(fr=fr, eps=eps_p, uz=u_p[2], s=s, fd=fd, ftot=p['fz'] + wbg))
budget = None
step_t = None
for l in open(f'{R}/run_slurm.log', errors='replace'):
    if l.startswith('EL_FORCE_BUDGET'):
        c = l.split(); st = int(c[2])
        if abs(st*2.5e-4 - T) < 1e-9: budget = (float(c[6]), float(c[10]))
print(f'replication check over all {len(rows)} particles: sum of replicated F_d,z {sum(x["fd"] for x in rows):.5e} N against EL_FORCE_BUDGET drag z {budget[0]:.5e} N (ratio {sum(x["fd"] for x in rows)/budget[0]:.4f}); pressure force z {budget[1]:.5e} N = V_p sum g_i')
rows = [x for x in rows if x['fr'] > 0]
W = sum(x['fr'] for x in rows)
mean = lambda k: sum(x['fr']*x[k] for x in rows)/W
Vwin = area*(zhi - zlo)
n_eff = W/Vwin
grad_meas = sum(x['fr']*x['ftot'] for x in rows)/Vwin
print(f'window particles (slab-weighted) {W:.1f}, solid fraction from them {W*vp/Vwin:.5f} (1 - eps_bed = {1-eps_bed:.5f})')
gcl = grad_cl(eps_bed, U)
print(f'closure gradient at eps_bed {eps_bed:.4f}, U {U}: {gcl:.1f} Pa/m; gradz_win / closure = {gw_meas/gcl:.4f}')
e_mean = mean('eps')
es = sorted(x['eps'] for x in rows)
print(f'(a) kernel-sampled eps at the window particles: mean {e_mean:.5f} (median {es[len(es)//2]:.4f}, 5-95 % {es[len(es)//20]:.4f}-{es[-len(es)//20]:.4f}); minus eps_bed {e_mean-eps_bed:+.5f}')
fcl = lambda e: drag_B(e, U/eps_bed)*U/eps_bed
fa = fcl(e_mean)/fcl(eps_bed)
fb = mean_cl = sum(x['fr']*fcl(x['eps']) for x in rows)/W/fcl(e_mean)
uz = mean('uz')
print(f'(c) sampled slip u_z: mean {uz:.5f} m/s against U/eps_bed {U/eps_bed:.5f} (ratio {uz/(U/eps_bed):.4f}); against U/<eps_p> {U/e_mean:.5f}')
fdm = mean('fd')
fc = fdm/(sum(x['fr']*fcl(x['eps']) for x in rows)/W)
print(f'drag per particle: closure at eps_bed {fcl(eps_bed):.4e} N; (a) at <eps_p> x{fa:.4f}; (b) mean over eps_p x{fb:.4f}; (c) with the sampled slip x{fc:.4f}; product {fa*fb*fc:.4f}')
# model A: grad = n <F_d> / eps if F_p = V_p grad; the measured F_tot gives F_d = F_tot - V_p grad_meas
# model A in a fixed bed: -grad p = n <F_d> / eps if the pressure force is V_p grad p
g_pred = n_eff*fdm/eps_bed
print(f'(d) model-A identity: n<F_d>/eps_bed from the replicated drag {g_pred:.1f} Pa/m against gradz_win {gw_meas:.1f} (ratio {gw_meas/g_pred:.4f}); closure {gcl:.1f}')
print(f'split of gradz_win/closure = {gw_meas/gcl:.4f}: (a) porosity {fa:.4f} x (b) nonlinearity {fb:.4f} x (c) slip {fc:.4f} x (d+e) rest {gw_meas/gcl/(fa*fb*fc):.4f}')
# deposited (element) porosity over the window, z-overlap weighted
num = 0.0; den = 0.0
for (i, j, k), e in elems.items():
    z0, z1 = k*h, (k+1)*h
    ov = min(z1, zhi) - max(z0, zlo)
    if ov > 0: num += ov*e; den += ov
eps_dep = num/den
print(f'deposited element porosity over the window {eps_dep:.5f} (exact slab {eps_bed:.5f}); identity with it: n<F_d>/eps_dep {n_eff*fdm/eps_dep:.1f} Pa/m, gradz_win / that {gw_meas/(n_eff*fdm/eps_dep):.4f}')
# porosity bias in the window interior (particles more than one kernel radius from the window ends)
inner = [x for x, p in zip(rows, [q for q in parts if slab(q['z'], zlo, zhi) > 0]) if p['z'] > zlo + width and p['z'] < zhi - width]
wi = sum(x['fr'] for x in inner)
zin0, zin1 = zlo + width, zhi - width
solid_in = sum(slab(q['z'], zin0, zin1) for q in parts)
eps_in = 1 - solid_in/(area*(zin1 - zin0))
print(f'interior ({zin0*1e3:.2f}-{zin1*1e3:.2f} mm): exact slab porosity {eps_in:.5f}, kernel-sampled mean at its particles {sum(x["fr"]*x["eps"] for x in inner)/wi:.5f}, sampled slip / (U/eps) {sum(x["fr"]*x["uz"] for x in inner)/wi/(U/eps_in):.4f}')
# interior gradient from the vertex pressure (Pressure_V, the projected P1 field) between vertex planes
pv = {}
for f in glob.glob(f'{R}/_vtk/res_node_*.{FR}.vtu'):
    piece = ET.parse(f).getroot().find('UnstructuredGrid').find('Piece')
    pts = arr(piece, 'Points'); pr = arr(piece, 'Pressure_V')
    for i in range(len(pr)):
        key = (round(pts[3*i]/h) % nx, round(pts[3*i+1]/h) % nx, round(pts[3*i+2]/h))
        pv.setdefault(key, []).append(pr[i])
def plane(k):
    v = [sum(pv[(i, j, k)])/len(pv[(i, j, k)]) for i in range(nx) for j in range(nx)]
    return sum(v)/len(v)
for k0, k1 in ((1, 9), (2, 8)):
    z0, z1 = k0*h, k1*h
    gi = (plane(k0) - plane(k1))/(z1 - z0)
    e01 = 1 - sum(slab(q['z'], z0, z1) for q in parts)/(area*(z1 - z0))
    print(f'vertex-pressure gradient {z0*1e3:.1f}-{z1*1e3:.1f} mm: {gi:.1f} Pa/m, exact porosity {e01:.5f}, closure {grad_cl(e01, U):.1f}, ratio {gi/grad_cl(e01, U):.4f}')
