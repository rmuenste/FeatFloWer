"""G5 and continuity-source means for a 04_voidfraction VB run from run_slurm.log.
usage: vft_g5.py <rundir> [--tmin 0.002]
G5 = uf_super_z / (N Vp |w| / V_box), N Vp |w| / V_box = 4.24115e-4 m/s (SPEC 04; V_box 0.01 m^3,
MOM = 4.24115e-6 m^4/s). Source lines (EL_VA_SOURCE/MOMENT/TSRC) may carry a trailing 'corr= k'
when ELVASourceCorrectors > 0: k = 0 predictor, 1..n corrector solves, -1 final assembly (no
solve). Means use the last SOLVE pass per step (largest k >= 0); the final assembly is reported
separately together with the EL_VA_RESIDUAL history."""
import re, sys, math, os
run = sys.argv[1]; tmin = float(sys.argv[sys.argv.index("--tmin")+1]) if "--tmin" in sys.argv else 0.002
REF, MOM, V = 4.24115e-4, 4.24115e-6, 0.01
num = r"([-+0-9.Ee]+)"
ms = re.compile(r"^EL_MEAN_SLIP t=\s*"+num+r".*uf_super=\s*"+num+r"\s+"+num+r"\s+"+num+r"(?:.*uz_unw=\s*"+num+r")?")
src = re.compile(r"^EL_VA_SOURCE t=\s*"+num+r" step= (\d+) int_g=\s*"+num+r" gbar=\s*"+num+r" L1_dedt=\s*"+num+r" L1_adv=\s*"+num+r"(?:\s+corr= (-?\d+))?")
mom = re.compile(r"^EL_VA_MOMENT t=\s*"+num+r" step= (\d+) M1raw=\s*"+num+r" Q2=\s*"+num+r" adv_net=\s*"+num+r" S2=\s*"+num+r"(?:\s+corr= (-?\d+))?")
tsrc = re.compile(r"^EL_VA_TSRC t=\s*"+num+r" step= (\d+) dedt_net=\s*"+num+r" L1_dcont=\s*"+num+r" scheme= (\S+)(?:\s+corr= (-?\d+))?")
cmp = re.compile(r"^EL_VA_COMPAT t=\s*"+num+r" step= (\d+) int_g_raw=\s*"+num+r" int_g_corr=\s*"+num+r" W=\s*"+num+r" z_alpha=\s*"+num+r" mshift=\s*"+num+r" mode= (\S+)(?:\s+corr= (-?\d+))?")
res = re.compile(r"^EL_VA_RESIDUAL t=\s*"+num+r" step= (\d+) corr= (-?\d+) L1=\s*"+num+r" max=\s*"+num+r" int_g=\s*"+num)
g5, uz = [], []
S, M, T, Sfin, Mfin, RS, C = {}, {}, {}, {}, {}, {}, {}
def keep(d, k, corr, val):
    c = int(corr) if corr is not None else 0
    if c < 0: return True
    if k not in d or c >= d[k][0]: d[k] = (c, val)
    return False
for line in open(os.path.join(run, "run_slurm.log"), errors="replace"):
    m = ms.match(line)
    if m:
        t = float(m.group(1))
        if t >= tmin:
            g5.append(float(m.group(4))/REF)
            if m.group(5): uz.append((t, float(m.group(5))))
        continue
    m = src.match(line)
    if m:
        v = (float(m.group(1)), float(m.group(3)), float(m.group(5)), float(m.group(6)))
        if keep(S, int(m.group(2)), m.group(7), v): Sfin[int(m.group(2))] = v
        continue
    m = mom.match(line)
    if m:
        v = (float(m.group(1)), float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6)))
        if keep(M, int(m.group(2)), m.group(7), v): Mfin[int(m.group(2))] = v
        continue
    m = tsrc.match(line)
    if m:
        keep(T, int(m.group(2)), m.group(6), (float(m.group(1)), float(m.group(3)), float(m.group(4)), m.group(5)))
        continue
    m = cmp.match(line)
    if m:
        keep(C, int(m.group(2)), m.group(9), (float(m.group(1)), float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6)), float(m.group(7))))
        continue
    m = res.match(line)
    if m:
        RS.setdefault(int(m.group(2)), {})[int(m.group(3))] = (float(m.group(1)), float(m.group(4)), float(m.group(5)), float(m.group(6)))
npass = max([v[0] for v in S.values()] + [0])
S = {k: v[1] for k, v in S.items()}; M = {k: v[1] for k, v in M.items()}; T = {k: v[1] for k, v in T.items()}; C = {k: v[1] for k, v in C.items()}
mean = lambda x: sum(x)/len(x) if x else float("nan")
sd = lambda x: math.sqrt(sum((v-mean(x))**2 for v in x)/(len(x)-1)) if len(x) > 1 else float("nan")
print(f"{os.path.basename(os.path.normpath(run))}: {len(g5)} EL_MEAN_SLIP samples t >= {tmin}" + (f"; source means from solve pass corr={npass}" if npass else ""))
print(f"  G5 (uf_super_z / {REF:g})  mean {mean(g5):.4f}  sd {sd(g5):.4f}  last {g5[-1] if g5 else float('nan'):.4f}")
ss = [S[k] for k in sorted(S) if S[k][0] >= tmin]
print(f"  EL_VA_SOURCE ({len(ss)} steps): int_g mean {mean([s[1] for s in ss]):.4e}  L1_dedt mean {mean([s[2] for s in ss]):.4e}  L1_adv mean {mean([s[3] for s in ss]):.4e}")
if M:
    mm = [M[k] for k in sorted(M) if M[k][0] >= tmin]
    print(f"  EL_VA_MOMENT ({len(mm)} steps): M1raw mean {mean([x[1] for x in mm]):.4e}  Q2 mean {mean([x[2] for x in mm]):.4e}  adv_net mean {mean([x[3] for x in mm]):.4e}  S2 mean {mean([x[4] for x in mm]):.4e}")
    ig = mean([s[1] for s in ss]); an = mean([x[3] for x in mm])
    print(f"  int_g + adv_net mean {ig+an:.4e}")
if T:
    tt = [T[k] for k in sorted(T) if T[k][0] >= tmin]
    print(f"  EL_VA_TSRC ({len(tt)} steps, scheme {tt[-1][3] if tt else '?'}): dedt_net mean {mean([x[1] for x in tt]):.4e}  L1_dcont mean {mean([x[2] for x in tt]):.4e}")
if C:
    cc = [C[k] for k in sorted(C) if C[k][0] >= tmin]
    print(f"  EL_VA_COMPAT ({len(cc)} steps): int_g_raw mean {mean([x[1] for x in cc]):.4e}  int_g_corr mean {mean([x[2] for x in cc]):.3e} (max |.| {max(abs(x[2]) for x in cc):.3e})  W mean {mean([x[3] for x in cc]):.4e}  z_alpha mean {mean([x[4] for x in cc]):.4f}  mshift mean {mean([x[5] for x in cc]):.4e}  mshift/MOM {mean([x[5] for x in cc])/MOM:.4f}")
if RS:
    for c in sorted({c for k in RS for c in RS[k]}):
        rr = [RS[k][c] for k in sorted(RS) if c in RS[k] and RS[k][c][0] >= tmin]
        print(f"  EL_VA_RESIDUAL corr {c:2d} ({len(rr)} steps): L1 mean {mean([x[1] for x in rr]):.4e}  max mean {mean([x[2] for x in rr]):.4e}  int_g mean {mean([x[3] for x in rr]):.4e}")
    if Sfin:
        sf = [Sfin[k] for k in sorted(Sfin) if Sfin[k][0] >= tmin]
        print(f"  final-assembly source (corr -1): int_g mean {mean([x[1] for x in sf]):.4e}")
if uz:
    uzd = {round(t, 9): v for t, v in uz}
    ref = Sfin if Sfin else S; refM = Mfin if Mfin else M
    resid, rel = [], []
    for k in sorted(refM):
        t = refM[k][0]
        if t < tmin or round(t, 9) not in uzd or k not in ref: continue
        zs = C[k][4] if k in C else 0.5   # sink position: z_alpha (alpha removal) or the box centre (uniform)
        lhs = uzd[round(t, 9)]*V; rhs = -refM[k][1] + zs*ref[k][1]
        resid.append(lhs-rhs); rel.append((lhs-rhs)/abs(rhs) if rhs else float("nan"))
    if resid:
        print(f"  moment identity ({'final assembly' if Sfin else 'predictor'}, n={len(resid)}): |lhs-rhs| max {max(abs(r) for r in resid):.3e}, rel rms {math.sqrt(mean([r*r for r in rel])):.3e}")
    print(f"  integral u_z / MOM mean {mean([v*V/MOM for t, v in uz]):.4f} (uz_unw*V/MOM)")
