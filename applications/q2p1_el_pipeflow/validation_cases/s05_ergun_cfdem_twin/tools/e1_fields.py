#!/usr/bin/env python3
"""E1 field post-processing from the solve-level VTU pieces of one frame.
usage: e1_fields.py <rundir> <frame, e.g. 00099> --fb 500 --d 9.862961e-4 [--rsplit 0.012]
Per element: eps (ELEpsilonF, cell), u_z = mean of the 8 vertex values, volume from the
hexahedron (6-tetra split), centroid radius. Reports the bed integrals, the core / wall-ring
split, and the closure applied element by element to the actual eps field (n_e F(u_e) = rho f_b
with n_e = (1-eps_e)/V_p, no wall friction), i.e. what the coded drag predicts for this field."""
import argparse, glob, math, re, sys
import numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('rundir'); ap.add_argument('frame'); ap.add_argument('--fb',type=float,required=True)
ap.add_argument('--d',type=float,required=True); ap.add_argument('--rsplit',type=float,default=0.012)
ap.add_argument('--rho',type=float,default=10.0); ap.add_argument('--mu',type=float,default=1.5e-3); ap.add_argument('--L',type=float,default=0.0554); a=ap.parse_args()
def arr(txt,name,dtype=float):
    m=re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>'%name,txt,re.S); return np.array(m.group(1).split(),dtype=dtype)
E=[];UZ=[];VOL=[];RC=[]
for f in sorted(glob.glob(f"{a.rundir}/_vtk/res_node_*.{a.frame}.vtu")):
    t=open(f,errors='replace').read()
    P=arr(t,'Points').reshape(-1,3); con=arr(t,'connectivity',int).reshape(-1,8); u=arr(t,'Velocity').reshape(-1,3); e=arr(t,'ELEpsilonF')
    X=P[con]                                   # (ne,8,3)
    c=X.mean(axis=1)
    tets=[(0,1,3,4),(1,2,3,6),(1,3,4,6),(1,4,5,6),(3,4,6,7)]   # 5-tet split of a hex (VTK ordering)
    vol=np.zeros(len(con))
    for i,j,k,l in tets:
        vol+=np.abs(np.einsum('ij,ij->i',np.cross(X[:,j]-X[:,i],X[:,k]-X[:,i]),X[:,l]-X[:,i]))/6
    E.append(e);UZ.append(u[con][:,:,2].mean(axis=1));VOL.append(vol);RC.append(np.hypot(c[:,0],c[:,1]))
e=np.concatenate(E);uz=np.concatenate(UZ);v=np.concatenate(VOL);r=np.concatenate(RC)
V=v.sum(); Vp=math.pi*a.d**3/6
print(f"elements {len(e)}  V_mesh {V:.6e} m^3  cross-section {V/a.L:.6e} m^2 = {V/a.L/(math.pi*0.01385**2):.5f} of pi 0.01385^2")
print(f"eps_mean {np.sum(e*v)/V:.6f}  eps min/max {e.min():.4f} / {e.max():.4f}  int(1-eps)dV {np.sum((1-e)*v):.6e}  = {np.sum((1-e)*v)/Vp:.1f} particles of d {a.d:.6e}")
U=np.sum(e*uz*v)/V
print(f"uf_super_z {U:.6e}")
def F(u,eps):
    re_=np.maximum(1e-12,a.rho*eps*a.d*u/a.mu); cd=(0.63+4.8/np.sqrt(re_))**2
    chi=3.7-0.65*np.exp(-0.5*(1.5-np.log10(re_))**2)
    return 0.5*a.rho*cd*math.pi*a.d**2/4*u*u*eps**(2-chi)
def solve(eps,G):
    lo=np.full_like(eps,1e-9); hi=np.full_like(eps,100.0); n=(1-eps)/Vp
    for _ in range(120):
        m=np.sqrt(lo*hi); low=n*F(m,eps)<G; lo=np.where(low,m,lo); hi=np.where(low,hi,m)
    return eps*lo
Uloc=solve(e,a.rho*a.fb); Upred=np.sum(Uloc*v)/V
Uuni=solve(np.array([np.sum(e*v)/V]),a.rho*a.fb)[0]
print(f"closure at eps_mean (uniform bed): {Uuni:.6e}  measured/that {U/Uuni:.4f}")
print(f"closure element by element on the eps field (no wall friction): {Upred:.6e}  measured/that {U/Upred:.4f}")
for name,m in (('core r<%.1f mm'%(a.rsplit*1e3),r<a.rsplit),('ring r>%.1f mm'%(a.rsplit*1e3),r>=a.rsplit)):
    vm=v[m].sum(); print(f"{name}: volume share {vm/V:.4f}  eps {np.sum(e[m]*v[m])/vm:.4f}  eps*u_z {np.sum(e[m]*uz[m]*v[m])/vm:.6e} ({np.sum(e[m]*uz[m]*v[m])/vm/U:.3f} of mean)  local-closure prediction {np.sum(Uloc[m]*v[m])/vm:.6e}")
edges=np.linspace(0,0.01385,11)
print("radial profile (r_mid mm, eps, eps*u_z/mean, local-closure/mean):")
for i in range(10):
    m=(r>=edges[i])&(r<edges[i+1])
    if m.any(): vm=v[m].sum(); print(f"  {0.5e3*(edges[i]+edges[i+1]):6.2f}  {np.sum(e[m]*v[m])/vm:.4f}  {np.sum(e[m]*uz[m]*v[m])/vm/U:.3f}  {np.sum(Uloc[m]*v[m])/vm/U:.3f}")
