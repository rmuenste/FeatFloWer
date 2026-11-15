#!/usr/bin/env python3
"""E1 reference built on the actual eps field: series in z (the flux is uniform along z by
continuity), parallel in r. For a set of elements find U_sup with  <n F(U/eps, eps)>_V = w rho f_b,
n = (1-eps)/V_p, coded difelice_std. Variants: uniform at eps_mean; z-series over the whole
cross-section (eps averaged per z layer); rings in parallel, each a z-series of its own elements.
usage: e1_series.py <rundir> <frame> --d D --fb 200 500 1000 [--share 0.977]"""
import argparse, glob, math, re
import numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('rundir'); ap.add_argument('frame'); ap.add_argument('--d',type=float,required=True)
ap.add_argument('--fb',type=float,nargs='+',default=[200,500,1000]); ap.add_argument('--share',type=float,default=1.0,help='fraction of rho f_b carried by drag (rest: wall friction)')
ap.add_argument('--rho',type=float,default=10.0); ap.add_argument('--mu',type=float,default=1.5e-3); ap.add_argument('--nr',type=int,default=10); a=ap.parse_args()
def arr(t,n,dt=float):
    m=re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>'%n,t,re.S); return np.array(m.group(1).split(),dtype=dt)
E=[];V=[];R=[];Z=[]
for f in sorted(glob.glob(f"{a.rundir}/_vtk/res_node_*.{a.frame}.vtu")):
    t=open(f,errors='replace').read(); P=arr(t,'Points').reshape(-1,3); c=arr(t,'connectivity',int).reshape(-1,8); X=P[c]; vol=np.zeros(len(c))
    for i,j,k,l in [(0,1,3,4),(1,2,3,6),(1,3,4,6),(1,4,5,6),(3,4,6,7)]:
        vol+=np.abs(np.einsum('ij,ij->i',np.cross(X[:,j]-X[:,i],X[:,k]-X[:,i]),X[:,l]-X[:,i]))/6
    cc=X.mean(axis=1); E.append(arr(t,'ELEpsilonF')); V.append(vol); R.append(np.hypot(cc[:,0],cc[:,1])); Z.append(cc[:,2])
e=np.concatenate(E); v=np.concatenate(V); r=np.concatenate(R); z=np.concatenate(Z); Vp=math.pi*a.d**3/6
def F(u,eps):
    re_=np.maximum(1e-12,a.rho*eps*a.d*u/a.mu); cd=(0.63+4.8/np.sqrt(re_))**2; chi=3.7-0.65*np.exp(-0.5*(1.5-np.log10(re_))**2)
    return 0.5*a.rho*cd*math.pi*a.d**2/4*u*u*eps**(2-chi)
def series(eps,w,G):          # one U for all entries, volume-weighted mean resistance = G
    lo,hi=1e-9,100.0
    for _ in range(100):
        m=math.sqrt(lo*hi)
        if np.sum(w*(1-eps)/Vp*F(m/eps,eps))/w.sum()<G: lo=m
        else: hi=m
    return lo
zl=np.unique(np.round(z,7)); em=np.sum(e*v)/v.sum()
ez=np.array([np.sum(e[np.isclose(z,q,atol=2e-7)]*v[np.isclose(z,q,atol=2e-7)])/v[np.isclose(z,q,atol=2e-7)].sum() for q in zl]); vz=np.array([v[np.isclose(z,q,atol=2e-7)].sum() for q in zl])
print(f"z layers {len(zl)}  eps(z): min {ez.min():.4f} max {ez.max():.4f} median {np.median(ez):.4f}; layers with eps > 0.5: {np.sum(ez>0.5)} (z = {', '.join('%.1f'%(1e3*q) for q in zl[ez>0.5])} mm)")
edges=np.linspace(0,r.max()*1.0001,a.nr+1)
for fb in a.fb:
    G=a.share*a.rho*fb
    Uu=series(np.array([em]),np.array([1.0]),G); Uz=series(ez,vz,G)
    Ur=0.0
    for i in range(a.nr):
        m=(r>=edges[i])&(r<edges[i+1]); Ur+=series(e[m],v[m],G)*v[m].sum()
    Ur/=v.sum()
    Ue=series(e,v,G)
    print(f"f_b {fb:6.0f} share {a.share}: uniform eps_mean {Uu:.6f} | z-series {Uz:.6f} | all elements in series {Ue:.6f} | rings parallel, series within ring {Ur:.6f}")
