#!/usr/bin/env python3
"""z profile of a frame: per element layer the area means of eps, u_z (vertex mean) and eps*u_z,
and the layer-mean pressure. usage: e2_zprofile.py <rundir> <frame> --uin U [--every 4]"""
import argparse, glob, re, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('rundir'); ap.add_argument('frame'); ap.add_argument('--uin',type=float,required=True); ap.add_argument('--every',type=int,default=4); a=ap.parse_args()
def arr(t,n,dt=float):
    m=re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>'%n,t,re.S); return np.array(m.group(1).split(),dtype=dt)
E=[];U=[];V=[];Z=[];P=[]
for f in sorted(glob.glob(f"{a.rundir}/_vtk/res_node_*.{a.frame}.vtu")):
    t=open(f,errors='replace').read(); X0=arr(t,'Points').reshape(-1,3); c=arr(t,'connectivity',int).reshape(-1,8); X=X0[c]; vol=np.zeros(len(c))
    for i,j,k,l in [(0,1,3,4),(1,2,3,6),(1,3,4,6),(1,4,5,6),(3,4,6,7)]:
        vol+=np.abs(np.einsum('ij,ij->i',np.cross(X[:,j]-X[:,i],X[:,k]-X[:,i]),X[:,l]-X[:,i]))/6
    E.append(arr(t,'ELEpsilonF')); U.append(arr(t,'Velocity').reshape(-1,3)[c][:,:,2].mean(axis=1)); V.append(vol); Z.append(X.mean(axis=1)[:,2]); P.append(arr(t,'Pressure_V')[c].mean(axis=1))
e=np.concatenate(E);u=np.concatenate(U);v=np.concatenate(V);z=np.concatenate(Z);p=np.concatenate(P)
zl=np.unique(np.round(z,7)); print(f"U_in {a.uin}; columns: z mm, eps, <u_z>/U_in, <eps u_z>/U_in, <p> Pa")
for q in zl[::a.every]:
    m=np.isclose(z,q,atol=2e-7); w=v[m]; print(f"  {1e3*q:6.2f}  {np.sum(e[m]*w)/w.sum():.4f}  {np.sum(u[m]*w)/w.sum()/a.uin:.4f}  {np.sum(e[m]*u[m]*w)/w.sum()/a.uin:.4f}  {np.sum(p[m]*w)/w.sum():9.3f}")
print(f"domain mean eps*u_z / U_in {np.sum(e*u*v)/v.sum()/a.uin:.4f}; p(z_min layer) - p(z_max layer) = {np.sum(p[np.isclose(z,zl[0],atol=2e-7)]*v[np.isclose(z,zl[0],atol=2e-7)])/v[np.isclose(z,zl[0],atol=2e-7)].sum()-np.sum(p[np.isclose(z,zl[-1],atol=2e-7)]*v[np.isclose(z,zl[-1],atol=2e-7)])/v[np.isclose(z,zl[-1],atol=2e-7)].sum():.3f} Pa")
