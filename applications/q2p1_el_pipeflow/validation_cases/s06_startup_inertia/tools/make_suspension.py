#!/usr/bin/env python3
"""Uniform suspension for the start-up case: jittered FCC lattice at a target solid fraction in a
closed box. usage: make_suspension.py --size LX LY LZ --d D --phi 0.4 [--jitter 0.35] [--seed 1] -o particles.xyz
The lattice constant follows from phi (4 spheres per cell); each centre is displaced uniformly by
up to jitter * (nearest-neighbour gap)/2 per coordinate direction norm-bounded, so no two spheres
can touch; centres keep a distance of 0.55 d from the walls. Prints N and the bulk solid fraction."""
import argparse, math, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('--size',nargs=3,type=float,required=True); ap.add_argument('--d',type=float,required=True)
ap.add_argument('--phi',type=float,default=0.4); ap.add_argument('--jitter',type=float,default=0.35); ap.add_argument('--seed',type=int,default=1)
ap.add_argument('-o',default='particles.xyz'); a=ap.parse_args()
d=a.d; vp=math.pi*d**3/6; al=(4*vp/a.phi)**(1/3); nn=al/math.sqrt(2); gap=nn-d
basis=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])*al
n=[int(math.ceil(L/al))+1 for L in a.size]
I=np.array(np.meshgrid(*[np.arange(k) for k in n],indexing='ij')).reshape(3,-1).T*al
P=(I[:,None,:]+basis[None,:,:]).reshape(-1,3)
# centre the lattice in the box
P+= (np.array(a.size)-(P.max(axis=0)+P.min(axis=0)))/2
rng=np.random.default_rng(a.seed); v=rng.normal(size=P.shape); v/=np.linalg.norm(v,axis=1)[:,None]
P+=v*(rng.random(len(P))[:,None]**(1/3))*a.jitter*gap/2
m=np.all((P>0.55*d)&(P<np.array(a.size)-0.55*d),axis=1); P=P[m]
np.savetxt(a.o,P,fmt='%.10e')
V=np.prod(a.size); print(f"N {len(P)}  lattice constant {al/d:.4f} d  nearest neighbour {nn/d:.4f} d  min gap after jitter >= {(1-a.jitter)*gap/d:.4f} d")
print(f"box solid fraction {len(P)*vp/V:.5f}; nominal lattice {a.phi}")
z=P[:,2]; lo,hi=a.size[2]/3,2*a.size[2]/3; mm=(z>lo)&(z<hi); print(f"middle third: N {mm.sum()}  solid fraction {mm.sum()*vp/(V/3):.5f}")
