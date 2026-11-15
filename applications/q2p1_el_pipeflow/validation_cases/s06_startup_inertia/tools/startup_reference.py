#!/usr/bin/env python3
"""1-D start-up of a homogeneous suspension from rest, zero net flux (eps u + phi v = 0), for the
momentum formulations of FeatFloWer / CFDEM. Slip s = u - v, v = -eps s. Di Felice chi(Re) as
coded (el_forces.f90, difelice_std): per-particle drag F(s) on slip s at porosity eps.
  A, Bfull : (rho_f phi + rho_p eps)       ds/dt = drho g - n F(s) / (eps phi)
  none+pf+fb (unweighted inertia, -eps grad p): (rho_f phi/eps + rho_p eps) ds/dt = same
  B        : eps phi rho_p                 ds/dt = phi drho g - n F(s) / eps
usage: startup_reference.py --phi 0.4 --d 1e-3 --rhof 1000 --rhop 1100 --mu 1e-3 [--g 9.81] [--tend 0.05] [--table]
Prints per formulation: s_inf, v_inf = -eps s_inf, t63 (time to 63.2 % of s_inf), linear tau at s -> 0."""
import argparse, math

def fit_tau(t,v,tmax):
    """least-squares fit of v = A (1 - exp(-t/tau)) on t <= tmax: scan tau, A linear"""
    import numpy as np
    t=np.asarray(t); v=np.asarray(v); m=t<=tmax; t=t[m]; v=v[m]; best=None
    for tau in np.geomspace(1e-4,1e-1,3000):
        b=1-np.exp(-t/tau); A=(b@v)/(b@b); r=((A*b-v)**2).sum()
        if best is None or r<best[0]: best=(r,tau,A)
    return best[1],best[2]
ap=argparse.ArgumentParser(); ap.add_argument('--tfit',type=float,default=0.02); ap.add_argument('--phi',type=float,required=True); ap.add_argument('--d',type=float,required=True)
ap.add_argument('--rhof',type=float,required=True); ap.add_argument('--rhop',type=float,required=True); ap.add_argument('--mu',type=float,required=True)
ap.add_argument('--g',type=float,default=9.81); ap.add_argument('--tend',type=float,default=0.05); ap.add_argument('--table',action='store_true'); a=ap.parse_args()
phi=a.phi; eps=1-phi; d=a.d; vp=math.pi*d**3/6; n=phi/vp; drho=a.rhop-a.rhof
def F(s):
    s=abs(s)
    if s<1e-14: return 0.0
    re=a.rhof*eps*d*s/a.mu; cd=(0.63+4.8/math.sqrt(re))**2; chi=3.7-0.65*math.exp(-0.5*(1.5-math.log10(re))**2)
    return 0.5*a.rhof*cd*math.pi*d*d/4*s*s*eps**(2-chi)
forms={'A/Bfull':(a.rhof*phi+a.rhop*eps, lambda s: drho*a.g-n*F(s)/(eps*phi)),
       'none+pf+fb':(a.rhof*phi/eps+a.rhop*eps, lambda s: drho*a.g-n*F(s)/(eps*phi)),
       'B':(eps*phi*a.rhop, lambda s: phi*drho*a.g-n*F(s)/eps)}
def run(m,rhs,dt=1e-6):
    s=0.0; t=0.0; out=[(0.0,0.0)]
    while t<a.tend:
        k1=rhs(s)/m; k2=rhs(s+0.5*dt*k1)/m; k3=rhs(s+0.5*dt*k2)/m; k4=rhs(s+dt*k3)/m
        s+=dt*(k1+2*k2+2*k3+k4)/6; t+=dt; out.append((t,s))
    return out
res={}
for name,(m,rhs) in forms.items():
    lo,hi=0.0,1.0
    for _ in range(200):
        mid=0.5*(lo+hi)
        if rhs(mid)>0: lo=mid
        else: hi=mid
    sinf=lo; tr=run(m,rhs); t63=next(t for t,s in tr if s>=0.632120558*sinf)
    k0=n*F(1e-9)/1e-9   # beta at s -> 0
    tau0=(m*eps*phi/k0) if name!='B' else (m*eps/k0)
    res[name]=(sinf,t63,tr); print(f"{name:12s} s_inf {sinf:.6e} m/s  v_inf {-eps*sinf:.6e}  Re_inf {a.rhof*eps*d*sinf/a.mu:.3f}  t63 {t63*1e3:.4f} ms  tau(s->0) {tau0*1e3:.4f} ms")
fits={}
for k,(sinf,t63,tr) in res.items():
    tt=[x[0] for x in tr[::250]]; vv=[-eps*x[1] for x in tr[::250]]; fits[k]=fit_tau(tt,vv,a.tfit); print(f"{k:12s} fit on t <= {1e3*a.tfit:.0f} ms: tau {1e3*fits[k][0]:.4f} ms, amplitude {fits[k][1]:.6e} m/s")
print("fitted tau ratios to A/Bfull: "+", ".join(f"{k} {v[0]/fits['A/Bfull'][0]:.4f}" for k,v in fits.items()))
ref=res['A/Bfull'][1]; print("t63 ratios to A/Bfull: "+", ".join(f"{k} {v[1]/ref:.4f}" for k,v in res.items()))
if a.table:
    print("# t  v_A/Bfull  v_none+pf+fb  v_B")
    N=len(res['B'][2])
    for i in range(0,N,max(1,N//200)): print(f"{res['B'][2][i][0]:.6e} "+" ".join(f"{-eps*res[k][2][i][1]:.6e}" for k in ('A/Bfull','none+pf+fb','B')))
