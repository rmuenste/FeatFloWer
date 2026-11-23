#!/usr/bin/env python3
"""Kinematic-wave (Kynch) reference for the batch-settling case, with CFDEM's Di Felice drag.

Steady homogeneous settling, zero net volume flux, model A/B/Bfull share the same balance:
    F_drag(u_r, eps) = V_p (rho_p - rho_f) g eps
    F_drag = 1/8 Cd rho_f pi d^2 eps^2 |u_r| u_r eps^(-chi),  Re = eps |u_r| d / nu,
    Cd = (0.63 + 4.8/sqrt(Re))^2,  chi = 3.7 - 0.65 exp(-(1.5 - log10 Re)^2 / 2)   (DiFeliceDrag.C)
u_r is the interstitial slip. Particle velocity v_p = -eps u_r; flux f(phi) = phi eps u_r.
Upper front (suspension / clear fluid) speed: eps0 u_r0.
Lower front (sediment / suspension) speed: (f0 - f*) / (phi* - phi0), f* = 0 for a packed bed at rest.
usage: batch_reference.py [phi0 ...]   (env: D RHOP RHOF MU PHIPACK H0 G)
"""
import os, sys, math
def env(k, v): return float(os.environ.get(k, v))
d, rhop, rhof, mu = env('D', 1e-3), env('RHOP', 2700), env('RHOF', 1000), env('MU', 0.05)
g, phip, H0 = env('G', 9.81), env('PHIPACK', 0.60), env('H0', 0.16)
nu = mu / rhof; V = math.pi * d**3 / 6

def drag(ur, eps):
    Re = max(eps * ur * d / nu, 1e-30)
    Cd = (0.63 + 4.8 / math.sqrt(Re))**2
    chi = 3.7 - 0.65 * math.exp(-0.5 * (1.5 - math.log10(Re))**2)
    return 0.125 * Cd * rhof * math.pi * d**2 * eps**2 * ur * ur * eps**(-chi), Re, chi

def slip(eps):
    rhs = V * (rhop - rhof) * g * eps
    lo, hi = 1e-12, 100.0
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if drag(mid, eps)[0] < rhs: lo = mid
        else: hi = mid
    return mid

def flux(phi): return phi * (1 - phi) * slip(1 - phi)

if __name__ == '__main__':
    uSt = (rhop - rhof) * g * d**2 / (18 * mu)
    u1 = slip(1.0)
    print('d %g  rho_p %g  rho_f %g  mu %g  Stokes velocity %.6g m/s' % (d, rhop, rhof, mu, uSt))
    print('single particle (Di Felice form): u_t %.6g m/s = %.4f Stokes, Re_t %.4g' % (u1, u1 / uSt, u1 * d / nu))
    print('tau_p = rho_p d^2 / (18 mu) = %.4g s' % (rhop * d**2 / (18 * mu)))
    print('%6s %10s %7s %7s %11s %11s %9s %9s %8s' % ('phi0', 'u_r', 'Re', 'chi', 'v_upper', 'v_lower', 't_meet', 'z_meet', 'mode'))
    for a in (sys.argv[1:] or ['0.10', '0.15', '0.20']):
        phi0 = float(a); eps = 1 - phi0; ur = slip(eps); _, Re, chi = drag(ur, eps)
        f0 = phi0 * eps * ur; vu = eps * ur; vl = f0 / (phip - phi0)
        tm = H0 / (vu + vl); zm = vl * tm
        # Kynch mode I check: chord from (phi0,f0) to (phip,0) must stay below the flux curve
        ok = all(flux(phi0 + (phip - phi0) * i / 200) >= f0 * (1 - i / 200) - 1e-15 for i in range(1, 200))
        print('%6.3f %10.6f %7.4f %7.4f %11.6f %11.6f %9.3f %9.5f %8s' % (phi0, ur, Re, chi, vu, vl, tm, zm, 'I' if ok else 'not I'))
