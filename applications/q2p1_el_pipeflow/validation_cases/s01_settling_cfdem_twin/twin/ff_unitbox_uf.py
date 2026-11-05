"""u_f at the particle for the unit-box FF settling runs (sc.FF cases), using twin_post's inversion
with the FF case parameters instead of the CFDEM ones.  usage: ff_unitbox_uf.py <case> <run> [<oneway run>]"""
import sys, os
sys.path.insert(0, "/data/warehouse17/rmuenste/code/liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/01_settling/twin")
sys.path.insert(0, "/data/warehouse17/rmuenste/code/liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/01_settling")
import settling_cases as sc
import twin_post as tp
case, run = sys.argv[1], sys.argv[2]
one = sys.argv[3] if len(sys.argv) > 3 else None
c = sc.FF[case]
tp.params = lambda _case: dict(d=c["d"], rho_p=c["rho_p"], rho_f=sc.FF_RHO_F, nu=c["nu"], g=sc.FF_G)
rows = tp.read_ff(run, case, one, "run_slurm.log")
tp.summarise(rows, case, "FeatFloWer(unit box)", "ffstd")
if one is None:
    eps = [r[3] for r in rows if r[3] == r[3]]
    print(f"  eps inverted: first {eps[0]:.6f}, last {eps[-1]:.6f}, n {len(eps)}")
