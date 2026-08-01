#!/usr/bin/env python3
"""Spatial/temporal error decomposition for the ten Cate D1.3 ladder runs.

Fits the additive model  err(L, dt) = S(L) + c * dt**p  to the measured
peak-velocity errors (vs ten Cate Table II printed experimental peaks),
per case, for temporal orders p = 1 and p = 2 (the observed order is not
yet pinned; two dt points per level cannot distinguish them).

Peaks are the argmin of SED_BENCH_VEL v_z from the run logs listed below
(datasheet rows e4_ladder_*, e4_l*_dt*, e1_l*). Rerun after new ladder
points land and update MEASURED.

Also prints Richardson extrapolations h->0 at dt->0 (orders 1 and 2 in
h) as an estimate of the converged solution -- for E1 compare against
ten Cate's own LBM (S1 ratio 0.894), not the experiment.
"""
import numpy as np

# case: reference peak [m/s] (Table II ratio x Table I u_inf)
REF = {'E4': -0.955 * 0.128, 'E1': -0.947 * 0.038}

# (level, D/h, dt [ms], measured peak [m/s], rundir)
MEASURED = {
    'E4': [
        (2, 11.4, 1.0, -0.12763, 'q2p1_dns_rundir_e4_l2'),
        (3, 23.9, 1.0, -0.12323, 'q2p1_dns_rundir_e4_l3'),
        (3, 23.9, 0.5, -0.12381, 'q2p1_dns_rundir_e4_l3_dt0p5'),
        (4, 49.1, 1.0, -0.12033, 'q2p1_dns_rundir_e4_l4'),
        (4, 49.1, 0.5, -0.12120, 'q2p1_dns_rundir_e4_l4_dt0p5'),
    ],
    'E1': [
        (3, 23.9, 1.0, -0.03484, 'q2p1_dns_rundir_e1_l3'),
        (3, 23.9, 0.5, -0.03499, 'q2p1_dns_rundir_e1_l3_dt0p5'),
        (4, 49.1, 1.0, -0.03410, 'q2p1_dns_rundir_e1_l4'),
    ],
}


def fit(case, p):
    """LSQ fit err = S(level) + c*dt^p; returns (S dict, c, residuals)."""
    data = MEASURED[case]
    levels = sorted({d[0] for d in data})
    err = np.array([(d[3] / REF[case] - 1) * 100 for d in data])  # in %
    A = np.zeros((len(data), len(levels) + 1))
    for i, d in enumerate(data):
        A[i, levels.index(d[0])] = 1.0
        A[i, -1] = d[2] ** p
    x, *_ = np.linalg.lstsq(A, err, rcond=None)
    res = A @ x - err
    return {L: x[k] for k, L in enumerate(levels)}, x[-1], res


def main():
    for case in ('E4', 'E1'):
        print('=== %s (reference peak %.5f m/s, Table II printed)' % (case, REF[case]))
        for _, _, dt, v, run in MEASURED[case]:
            print('    %-34s dt=%.2gms  err=%+.2f%%' % (run, dt, (v / REF[case] - 1) * 100))
        for p in (1, 2):
            S, c, res = fit(case, p)
            terms = '  '.join('S(L%d)=%+.2f' % (L, S[L]) for L in sorted(S))
            print('  p=%d fit: %s  T(1ms)=%+.2f pp  |res|max=%.2f pp'
                  % (p, terms, c, np.abs(res).max()))
        # Richardson h->0 using dt-corrected (dt->0) level values, p=1 temporal
        S, c, _ = fit(case, 1)
        Ls = sorted(S)
        if len(Ls) >= 2:
            s_fine, s_coarse = S[Ls[-1]], S[Ls[-2]]  # h halves per level
            for q, name in ((1, 'O(h) '), (2, 'O(h^2)')):
                s_inf = s_fine + (s_fine - s_coarse) / (2 ** q - 1)
                v_inf = REF[case] * (1 + s_inf / 100)
                print('  Richardson %s: S_inf=%+.2f pp -> v_inf=%.5f m/s (ratio %.3f)'
                      % (name, s_inf, v_inf, abs(v_inf) / (abs(REF[case]) /
                         {'E4': 0.955, 'E1': 0.947}[case])))
        print()


if __name__ == '__main__':
    main()
