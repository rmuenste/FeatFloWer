#!/usr/bin/env python3
"""Uniform suspension for the start-up case, wall to wall: jittered FCC lattice at solid fraction
phi whose outermost lattice planes sit 0.55 d from the side walls, so the cross-section has no
loose wall region. The box width follows from the number of lattice planes.
usage: make_suspension_tight.py [--d 1e-3] [--phi 0.4] [--planes 24] [--lz 0.06] [--seed 3] -o particles.xyz
Prints the box width L (use it for the mesh: make_box_mesh.py --size L L LZ), N and the solid
fraction of the box and of its middle third."""
import argparse, math, numpy as np
ap=argparse.ArgumentParser(); ap.add_argument('--d',type=float,default=1e-3); ap.add_argument('--phi',type=float,default=0.4)
ap.add_argument('--planes',type=int,default=24); ap.add_argument('--lz',type=float,default=0.06); ap.add_argument('--seed',type=int,default=3)
ap.add_argument('--jitter',type=float,default=0.35); ap.add_argument('-o',default='particles.xyz'); a=ap.parse_args()
d=a.d; vp=math.pi*d**3/6; al=(4*vp/a.phi)**(1/3); hp=al/2; off=0.55*d; L=2*off+(a.planes-1)*hp; LZ=a.lz
nzp=int((LZ-2*off)/hp)+1
I=np.array(np.meshgrid(np.arange(a.planes),np.arange(a.planes),np.arange(nzp),indexing='ij')).reshape(3,-1).T
I=I[(I.sum(axis=1)%2)==0]; P=off+I*hp
P[:,2]+=(LZ-(P[:,2].max()+P[:,2].min()))/2
gap=al/math.sqrt(2)-d; rng=np.random.default_rng(a.seed); v=rng.normal(size=P.shape); v/=np.linalg.norm(v,axis=1)[:,None]
P+=v*(rng.random(len(P))[:,None]**(1/3))*a.jitter*gap/2
assert P[:,:2].min()>0.5*d and P[:,:2].max()<L-0.5*d
np.savetxt(a.o,P,fmt='%.10e'); z=P[:,2]; m=(z>LZ/3)&(z<2*LZ/3)
print(f"L {L!r}  N {len(P)}  box solid fraction {len(P)*vp/(L*L*LZ):.5f}  middle third {m.sum()*vp/(L*L*LZ/3):.5f}  min wall gap {(P[:,:2].min()-0.5*d)/d:.3f} d")
