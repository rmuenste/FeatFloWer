#!/usr/bin/env python3
"""Mean particle settling velocity v_z(t) of the start-up case from the per-rank particle logs.
usage: startup_post.py <rundir> [--box 0.02 0.02 0.06] [--margin 0.004] [--ref "t63_ms v_inf"]
Particles are selected once, from their first logged position: middle third in z and at least
`margin` from the side walls. Prints N, v_inf (mean over the last 20 % of the samples), t63 (first
time |v| reaches 63.2 % of |v_inf|, linear interpolation) and a coarse table."""
import argparse, glob, numpy as np

def fit_tau(t,v,tmax):
    """least-squares fit of v = A (1 - exp(-t/tau)) on t <= tmax: scan tau, A linear"""
    import numpy as np
    t=np.asarray(t); v=np.asarray(v); m=t<=tmax; t=t[m]; v=v[m]; best=None
    for tau in np.geomspace(1e-4,1e-1,3000):
        b=1-np.exp(-t/tau); A=(b@v)/(b@b); r=((A*b-v)**2).sum()
        if best is None or r<best[0]: best=(r,tau,A)
    return best[1],best[2]
ap=argparse.ArgumentParser(); ap.add_argument('rundir'); ap.add_argument('--tfit',type=float,default=0.02); ap.add_argument('--box',nargs=3,type=float,default=[0.02,0.02,0.06]); ap.add_argument('--margin',type=float,default=0.004); a=ap.parse_args()
rows=[]
for f in sorted(glob.glob(a.rundir+'/el_particles_rank*.log')):
    d=np.loadtxt(f,comments='#',usecols=(0,1,8,9,10,13,4))  # t ip px py pz vz fz
    if d.size: rows.append(d.reshape(-1,7))
D=np.vstack(rows); t0=D[:,0].min(); first=D[D[:,0]==t0]
L=a.box; m=(first[:,4]>L[2]/3)&(first[:,4]<2*L[2]/3)&(first[:,2]>a.margin)&(first[:,2]<L[0]-a.margin)&(first[:,3]>a.margin)&(first[:,3]<L[1]-a.margin)
ids=set(first[m,1].astype(int)); sel=np.isin(D[:,1].astype(int),list(ids)); D=D[sel]
ts=np.unique(D[:,0]); v=np.array([D[D[:,0]==t,5].mean() for t in ts]); nlast=max(1,len(ts)//5); vinf=v[-nlast:].mean()
x=np.abs(v)/abs(vinf); k=int(np.argmax(x>=0.632120558)); t63=ts[k] if k==0 else ts[k-1]+(0.632120558-x[k-1])/(x[k]-x[k-1])*(ts[k]-ts[k-1])
print(f"{a.rundir}: selected {len(ids)} particles, {len(ts)} samples, t = {ts[0]:.5f} .. {ts[-1]:.5f}")
tau,A=fit_tau(ts,v,a.tfit); print(f"  fit on t <= {1e3*a.tfit:.0f} ms: tau {1e3*tau:.4f} ms, amplitude {A:.6e} m/s")
print(f"  v_inf {vinf:.6e} m/s (sd over the last {nlast} samples {v[-nlast:].std():.2e}); t63 {1e3*t63:.4f} ms")
for i in range(0,len(ts),max(1,len(ts)//12)): print(f"    t {1e3*ts[i]:7.3f} ms  v {v[i]: .6e}  v/v_inf {v[i]/vinf:.4f}")
