"""Drag closures of case 11 (minimum fluidisation) for the post-processing tools.

The one Python copy of the per-particle drag coefficient of EL_FORCES (model-A form, the
coefficient B of F_d = B (u_f - v)) and of the fixed-bed pressure gradient and minimum
fluidisation velocity derived from it. Used by fl11_post.py (gates F2, F3) and core_deficit.py
(deficit decomposition); the arithmetic is the one both tools carried before (review r2 P3:
one shared module instead of two copies).

  drag_b(eps, slip, d, rhof, mu, closure)       B of one particle at porosity eps and slip speed
                                                 slip (m/s): ergun, wenyu, or gidaspow (ergun for
                                                 eps <= 0.8, wenyu above, the switch of EL_FORCES)
  closure_grad(u, eps, d, rhof, mu, closure)    |dp/dz| of a fixed bed at superficial velocity u:
                                                 n F_d(u/eps) / eps (n = (1-eps)/V_p)
  closure_umf(eps, d, rhop, rhof, mu, g, closure)  u where closure_grad meets the buoyant weight
                                                 per volume (1-eps)(rhop-rhof) g (bisection)

Outside 0 < eps <= 1 (and for a non-finite input) the gradient and U_mf are NaN, never an
exception: the callers report such a value as unavailable."""
import math

CLOSURES = ('gidaspow', 'ergun', 'wenyu')


def drag_b(eps, slip, d, rhof, mu, closure='gidaspow'):
    vp = math.pi * d ** 3 / 6

    def ergun():
        return vp * (150 * (1 - eps) * mu / (eps * d * d) + 1.75 * rhof * slip / d)

    def wenyu():
        re_ = rhof * eps * slip * d / mu
        cds = 24 * mu / (rhof * eps * d) * (1 + 0.15 * re_ ** 0.687) if re_ < 1000 else 0.44 * slip
        return vp * 0.75 * cds * rhof * eps * eps ** -2.65 / d
    if closure == 'ergun':
        return ergun()
    if closure == 'wenyu':
        return wenyu()
    if closure == 'gidaspow':
        return ergun() if eps <= 0.8 else wenyu()
    raise ValueError('unknown closure %r (one of %s)' % (closure, ', '.join(CLOSURES)))


def _valid(eps, *more):
    return all(isinstance(x, (int, float)) and math.isfinite(x) for x in (eps,) + more) and 0 < eps <= 1


def closure_grad(u, eps, d, rhof, mu, closure='gidaspow'):
    """|dp/dz| of a fixed bed: n F_d(U/eps)/eps with the drag coefficient of EL_FORCES (model-A form)"""
    if not _valid(eps, u):
        return float('nan')
    vp = math.pi * d ** 3 / 6
    ur = u / eps
    return (1 - eps) / vp * drag_b(eps, ur, d, rhof, mu, closure) * ur / eps


def closure_umf(eps, d, rhop, rhof, mu, g, closure='gidaspow'):
    if not _valid(eps):
        return float('nan')
    target = (1 - eps) * (rhop - rhof) * g
    a, b = 1e-6, 10.0
    for _ in range(200):
        m = 0.5 * (a + b)
        if closure_grad(m, eps, d, rhof, mu, closure) < target:
            a = m
        else:
            b = m
    return 0.5 * (a + b)
