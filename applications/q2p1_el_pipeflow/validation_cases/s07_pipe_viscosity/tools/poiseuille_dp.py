#!/usr/bin/env python3
"""Pressure gradient of developed pipe flow from a VTK frame of q2p1_el_pipeflow: layer-mean
pressure against z, least-squares slope over a z window, compared with Hagen-Poiseuille
dp/dz = -8 mu U / R^2. Also the maximum axial velocity over U (2 when developed).
usage: poiseuille_dp.py <rundir> <frame> --uin U --mu MU --radius R [--zwin z0 z1]"""
import argparse, glob, re, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('rundir'); ap.add_argument('frame')
ap.add_argument('--uin',type=float,required=True); ap.add_argument('--mu',type=float,required=True)
ap.add_argument('--radius',type=float,required=True); ap.add_argument('--zwin',type=float,nargs=2,default=None)
a=ap.parse_args()
def arr(t,n,dt=float):
    m=re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>'%n,t,re.S); return np.array(m.group(1).split(),dtype=dt)
Z=[];P=[];V=[];W=[]
for f in sorted(glob.glob(f"{a.rundir}/_vtk/res_node_*.{a.frame}.vtu")):
    t=open(f,errors='replace').read(); X0=arr(t,'Points').reshape(-1,3); c=arr(t,'connectivity',int).reshape(-1,8); X=X0[c]; vol=np.zeros(len(c))
    for i,j,k,l in [(0,1,3,4),(1,2,3,6),(1,3,4,6),(1,4,5,6),(3,4,6,7)]:
        vol+=np.abs(np.einsum('ij,ij->i',np.cross(X[:,j]-X[:,i],X[:,k]-X[:,i]),X[:,l]-X[:,i]))/6
    Z.append(X.mean(axis=1)[:,2]); P.append(arr(t,'Pressure_V')[c].mean(axis=1)); V.append(vol); W.append(arr(t,'Velocity').reshape(-1,3)[:,2])
z=np.concatenate(Z);p=np.concatenate(P);v=np.concatenate(V);w=np.concatenate(W)
zl=np.unique(np.round(z,7)); pl=np.array([np.sum(p[np.isclose(z,q,atol=2e-7)]*v[np.isclose(z,q,atol=2e-7)])/v[np.isclose(z,q,atol=2e-7)].sum() for q in zl])
z0,z1=a.zwin if a.zwin else (zl[0]+0.5*(zl[-1]-zl[0]), zl[0]+0.9*(zl[-1]-zl[0]))
m=(zl>=z0)&(zl<=z1); s,i=np.polyfit(zl[m],pl[m],1); ref=-8*a.mu*a.uin/a.radius**2
print(f"{a.rundir} frame {a.frame}: dp/dz over z in ({1e3*z0:.1f}, {1e3*z1:.1f}) mm = {s:.6e} Pa/m; Hagen-Poiseuille with mu {a.mu:g}: {ref:.6e} Pa/m; ratio {s/ref:.4f}; max u_z / U_in {w.max()/a.uin:.4f}; p(first layer) - p(last layer) {pl[0]-pl[-1]:.6e} Pa")
