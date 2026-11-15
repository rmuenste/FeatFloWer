#!/usr/bin/env python3
"""E2 column post-processing from run_slurm.log and the particle logs.
usage: e2_post.py <rundir> [--u0 0.002 --u1 0.02 --tramp 0.4] [--fit 0.003 0.0075]
Per step: inlet U(t), sum of hydrodynamic particle forces (drag_z, pressure_z from EL_FORCE_BUDGET),
force-based dp = F_z / A (A = pi 0.01385^2, CFDEM forces.txt observable), domain superficial
velocity uf_super_z (EL_MEAN_SLIP; equals U_in for a conserved flux), mean particle velocity.
Fixed-bed slope d(dp)/dU by least squares over the fit window; bed weight plateau; bed height
from the particle logs (top = 99.5 percentile of z + r)."""
import argparse, glob, math, re, os
import numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('rundir'); ap.add_argument('--u0',type=float,default=0.002); ap.add_argument('--u1',type=float,default=0.02)
ap.add_argument('--tramp',type=float,default=0.4); ap.add_argument('--fit',type=float,nargs=2,default=[0.003,0.0075]); ap.add_argument('--r',type=float,default=4.969479e-4)
a=ap.parse_args(); A=math.pi*0.01385**2
num=r'[-+0-9.Ee]+'
dt=None; fb={}; ms={}
for line in open(os.path.join(a.rundir,'run_slurm.log'),errors='replace'):
    if line.startswith('TimeStep ='): dt=float(line.split('=')[1].replace('D','E'))
    elif line.startswith('EL_FORCE_BUDGET'):
        m=re.search(r'step= (\d+) drag=\s*(%s)\s+(%s)\s+(%s) pressure=\s*(%s)\s+(%s)\s+(%s).*grav_buoy=\s*(%s)\s+(%s)\s+(%s)'%((num,)*9),line)
        fb[int(m.group(1))]=(float(m.group(4)),float(m.group(7)),float(m.group(10)))
    elif line.startswith('EL_MEAN_SLIP t='):
        m=re.search(r't=\s*(%s) up=\s*%s\s+%s\s+(%s) uf_intr=.*uf_super=\s*%s\s+%s\s+(%s).*eps=\s*(%s)'%((num,)*8),line)
        ms[round(float(m.group(1))/dt)]=(float(m.group(2)),float(m.group(3)),float(m.group(4)))
steps=sorted(set(fb)&set(ms)); t=np.array(steps)*dt
U=a.u0+(a.u1-a.u0)*np.minimum(1,t/a.tramp)
drag=np.array([fb[s][0] for s in steps]); pres=np.array([fb[s][1] for s in steps]); grav=np.array([fb[s][2] for s in steps])
upz=np.array([ms[s][0] for s in steps]); ufz=np.array([ms[s][1] for s in steps])
dpd=drag/A; dpt=(drag+pres)/A; W=-grav.mean()/A
print(f"{os.path.basename(a.rundir.rstrip('/'))}: steps {len(steps)}, dt {dt}, bed weight/A {W:.2f} Pa (grav_buoy {grav.mean():.5e} N)")
m=(U>=a.fit[0])&(U<=a.fit[1])
for name,y in (('drag only',dpd),('drag + pressure force',dpt)):
    s,i=np.polyfit(U[m],y[m],1); print(f"  slope d(dp)/dU [{name}] over U in ({a.fit[0]}, {a.fit[1]}): {s:.1f} Pa s/m (intercept {i:.2f} Pa), dp at U=0.005: {np.interp(0.005,U,y):.2f}, at 0.01: {np.interp(0.01,U,y):.2f}, at 0.02: {np.interp(0.02,U,y):.2f} Pa")
print(f"  flux check uf_super_z / U_in: at U=0.005 {np.interp(0.005,U,ufz/U):.4f}, 0.01 {np.interp(0.01,U,ufz/U):.4f}, 0.02 {np.interp(0.02,U,ufz/U):.4f}")
print(f"  total hydro force / bed weight at end: {dpt[-1]/W:.4f}; max {np.max(dpt/W):.4f} at U {U[np.argmax(dpt/W)]:.5f}; mean particle w at end {upz[-1]:.3e}, max |w| {np.abs(upz).max():.3e} m/s")
above=np.where(dpt>=W)[0]
print("  first U with hydro force >= bed weight:", f"{U[above[0]]:.5f} m/s (t {t[above[0]]:.4f})" if len(above) else "not reached")
fr={}
for f in glob.glob(os.path.join(a.rundir,'el_particles_rank*.log')):
    for line in open(f):
        if line[0]=='#': continue
        c=line.split(); fr.setdefault(float(c[0]),{})[int(c[1])]=float(c[10])
ts=sorted(fr)
if ts:
    print("  bed height (t, U_in, H = p99.5(z)+r, mean z):")
    for tt in ts[::max(1,len(ts)//8)]+[ts[-1]]:
        z=np.array(list(fr[tt].values())); print(f"    {tt:.3f}  {a.u0+(a.u1-a.u0)*min(1,tt/a.tramp):.5f}  {np.percentile(z,99.5)+a.r:.6f}  {z.mean():.6f}")
