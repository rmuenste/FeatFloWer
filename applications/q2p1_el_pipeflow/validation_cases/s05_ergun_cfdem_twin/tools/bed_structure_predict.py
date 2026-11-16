#!/usr/bin/env python3
"""Fixed-bed pressure-drop slope of a bed AS IT IS against "closure at the core porosity x
equivalent height". usage: bed_structure_predict.py <particles.xyz> --d D [--R 0.01385]
  [--dz 1.065e-3] [--rwall 0.012] [--eps-core E --H H]

Each sphere is split into equal-volume sample points and binned into z layers and two rings
(core r < rwall, wall ring). Di Felice chi(Re) as coded in el_forces.f90 / DiFeliceDrag.C at a small
U (Stokes-like, slope per unit U). Per bin: g = n F(U/eps, eps) / U [unweighted, G = n F] and
g/eps [consistent model A, G = n F / eps]. Models: uniform reference, z series (layer eps), z series
with the two rings of each layer in parallel (equal pressure gradient)."""
import argparse, math, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('xyz'); ap.add_argument('--d',type=float,required=True)
ap.add_argument('--R',type=float,default=0.01385); ap.add_argument('--dz',type=float,default=1.065e-3)
ap.add_argument('--rwall',type=float,default=0.012); ap.add_argument('--eps-core',type=float); ap.add_argument('--H',type=float)
ap.add_argument('--U',type=float,default=5e-3); ap.add_argument('--rho',type=float,default=10.0); ap.add_argument('--mu',type=float,default=1.5e-3); a=ap.parse_args()
U=a.U; Vp=math.pi*a.d**3/6
def F(u,eps):
    re=a.rho*eps*a.d*u/a.mu; cd=(0.63+4.8/math.sqrt(re))**2; chi=3.7-0.65*math.exp(-0.5*(1.5-math.log10(re))**2)
    return 0.5*a.rho*cd*math.pi*a.d**2/4*u*u*eps**(2-chi)
def g(eps,cons):          # pressure gradient per unit superficial velocity of a bin
    eps=min(eps,1-1e-9); return (1-eps)/Vp*F(U/eps,eps)/U/(eps if cons else 1.0)
X=np.loadtxt(a.xyz)[:,:3]; N=len(X)
rng=np.random.default_rng(1); m=400; q=rng.normal(size=(m,3)); q/=np.linalg.norm(q,axis=1)[:,None]; q*=(a.d/2)*rng.random(m)[:,None]**(1/3)
P=(X[:,None,:]+q[None,:,:]).reshape(-1,3); w=Vp/m
nz=int(math.ceil((X[:,2].max()+a.d/2)/a.dz)); iz=np.clip((P[:,2]/a.dz).astype(int),0,nz-1); r=np.hypot(P[:,0],P[:,1]); ir=(r>=a.rwall).astype(int)
A=[math.pi*a.rwall**2, math.pi*(a.R**2-a.rwall**2)]; At=sum(A)
vs=np.zeros((nz,2)); np.add.at(vs,(iz,ir),w)
epsb=1-vs/(np.array(A)[None,:]*a.dz); epsl=1-vs.sum(axis=1)/(At*a.dz)
core=slice(3,nz-4); ec=a.eps_core if a.eps_core else float(np.mean(epsl[core])); H=a.H if a.H else N*Vp/(At*(1-ec))
print(f"N {N}, d {a.d}, layers {nz} x dz {a.dz*1e3:.3f} mm; core eps (layers 3..{nz-5}) {np.mean(epsl[core]):.4f}; reference eps {ec:.4f}, H {H*1e3:.2f} mm")
print("layer eps (all | core ring | wall ring):"); 
for k in range(nz): print(f"  {k:2d} z {1e3*(k+.5)*a.dz:5.2f}  {epsl[k]:.3f}  {epsb[k,0]:.3f}  {epsb[k,1]:.3f}")
for cons in (False,True):
    ref=g(ec,cons)*H; ser=sum(g(e,cons) for e in epsl)*a.dz
    par=sum(1.0/sum(A[i]/At/g(epsb[k,i],cons) for i in range(2)) for k in range(nz))*a.dz
    mid=np.mean([g(e,cons) for e in epsl[core]])/g(ec,cons)
    print(f"{'consistent G = nF/eps' if cons else 'unweighted G = nF     '}: uniform {ref:8.0f} Pa s/m | z series {ser:8.0f} ({ser/ref:.3f}) | + rings parallel {par:8.0f} ({par/ref:.3f}) | core layers local/uniform {mid:.3f}")
