#!/usr/bin/env python3
"""E1 like-for-like prediction: uniform bed, FeatFloWer difelice_std as coded (el_forces.f90),
n F(u) = w rho f_b with w = 1 (unweighted momentum, 'FF now') or w = eps (consistent).
usage: e1_predict.py --eps E --d D --n N_per_m3 --fb 200 500 1000"""
import argparse, math
ap=argparse.ArgumentParser(); ap.add_argument('--eps',type=float,required=True); ap.add_argument('--d',type=float,required=True)
ap.add_argument('--n',type=float,required=True); ap.add_argument('--fb',type=float,nargs='+',default=[200,500,1000])
ap.add_argument('--rho',type=float,default=10.0); ap.add_argument('--mu',type=float,default=1.5e-3); a=ap.parse_args()
def F(u,eps,d):
    re=max(1e-12,a.rho*eps*d*u/a.mu); cd=(0.63+4.8/math.sqrt(re))**2
    chi=3.7-0.65*math.exp(-0.5*(1.5-math.log10(re))**2)
    return 0.5*a.rho*cd*math.pi*d*d/4*u*u*eps**(2-chi)
def solve(G):
    lo,hi=1e-9,10.0
    for _ in range(200):
        m=math.sqrt(lo*hi)
        if a.n*F(m,a.eps,a.d)<G: lo=m
        else: hi=m
    return a.eps*lo
for fb in a.fb:
    print(f"f_b {fb:6.0f}: U_sup FF-now {solve(a.rho*fb):.6f}  consistent {solve(a.eps*a.rho*fb):.6f}")
