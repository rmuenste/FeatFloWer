"""G5 and continuity-source means for a 04_voidfraction VB run from run_slurm.log.
usage: vft_g5.py <rundir> [--tmin 0.002]   G5 = uf_super_z / (N Vp |w| / V_box) with N Vp |w|/V = 4.24115e-4 m/s"""
import re, sys, math, os
run = sys.argv[1]; tmin = float(sys.argv[sys.argv.index("--tmin")+1]) if "--tmin" in sys.argv else 0.002
REF = 4.24115e-4   # N Vp |w| / V_box  (SPEC 04; V_box = 0.01 m^3, MOM = 4.24115e-6 m^4/s)
MOM = 4.24115e-6
num = r"([-+0-9.Ee]+)"
ms = re.compile(r"^EL_MEAN_SLIP t=\s*"+num+r".*uf_super=\s*"+num+r"\s+"+num+r"\s+"+num+r"(?:.*uz_unw=\s*"+num+r")?")
src = re.compile(r"^EL_VA_SOURCE t=\s*"+num+r" step= (\d+) int_g=\s*"+num+r" gbar=\s*"+num+r" L1_dedt=\s*"+num+r" L1_adv=\s*"+num)
mom = re.compile(r"^EL_VA_MOMENT t=\s*"+num+r" step= (\d+) M1raw=\s*"+num+r" Q2=\s*"+num+r" adv_net=\s*"+num+r" S2=\s*"+num)
g5, uz, S, M = [], [], {}, {}
for line in open(os.path.join(run, "run_slurm.log"), errors="replace"):
    m = ms.match(line)
    if m:
        t = float(m.group(1))
        if t >= tmin:
            g5.append(float(m.group(4))/REF)
            if m.group(5): uz.append((t, float(m.group(5))))
        continue
    m = src.match(line)
    if m: S[int(m.group(2))] = (float(m.group(1)), float(m.group(3)), float(m.group(5)), float(m.group(6))); continue
    m = mom.match(line)
    if m: M[int(m.group(2))] = (float(m.group(1)), float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6)))
mean = lambda x: sum(x)/len(x) if x else float("nan")
sd = lambda x: math.sqrt(sum((v-mean(x))**2 for v in x)/(len(x)-1)) if len(x) > 1 else float("nan")
print(f"{os.path.basename(os.path.normpath(run))}: {len(g5)} EL_MEAN_SLIP samples t >= {tmin}")
print(f"  G5 (uf_super_z / {REF:g})  mean {mean(g5):.4f}  sd {sd(g5):.4f}  last {g5[-1] if g5 else float('nan'):.4f}")
ss = [S[k] for k in sorted(S) if S[k][0] >= tmin]
print(f"  EL_VA_SOURCE ({len(ss)} steps): int_g mean {mean([s[1] for s in ss]):.4e}  L1_dedt mean {mean([s[2] for s in ss]):.4e}  L1_adv mean {mean([s[3] for s in ss]):.4e}")
if M:
    mm = [M[k] for k in sorted(M) if M[k][0] >= tmin]
    print(f"  EL_VA_MOMENT ({len(mm)} steps): M1raw mean {mean([x[1] for x in mm]):.4e}  Q2 mean {mean([x[2] for x in mm]):.4e}  adv_net mean {mean([x[3] for x in mm]):.4e}  S2 mean {mean([x[4] for x in mm]):.4e}")
    ig = mean([s[1] for s in ss]); an = mean([x[3] for x in mm])
    print(f"  int_g + adv_net mean {ig+an:.4e}")
if uz:
    V = 0.01
    uzd = {round(t, 9): v for t, v in uz}
    # identity per step: integral u_z dV = -M1raw + 0.5 int_g   (u_z integral = uz_unw * V)
    res, rel = [], []
    for k in sorted(M):
        t = M[k][0]
        if t < tmin or round(t, 9) not in uzd or k not in S: continue
        lhs = uzd[round(t, 9)]*V; rhs = -M[k][1] + 0.5*S[k][1]
        res.append(lhs-rhs); rel.append((lhs-rhs)/abs(rhs) if rhs else float("nan"))
    print(f"  moment identity (n={len(res)}): |lhs-rhs| max {max(abs(r) for r in res):.3e}, rel rms {math.sqrt(mean([r*r for r in rel])):.3e}")
    print(f"  integral u_z / MOM mean {mean([v*V/MOM for t, v in uz]):.4f} (uz_unw*V/MOM)")
