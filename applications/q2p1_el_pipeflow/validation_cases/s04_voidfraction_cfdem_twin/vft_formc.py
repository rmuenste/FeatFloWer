#!/usr/bin/env python3
# usage: vft_formc.py <run_slurm.log> [...]  -- form-C readings of a VB run from the EL_VA_FLUX / EL_VA_FLUXMOM
# lines (final pass of each step) and the campaign G5 from EL_MEAN_SLIP; means over t >= 0.002 and over the settled
# window t >= 0.1 (the kernel deposition ramps up over the first ~50 steps). MOM = N V_p |w| / V_box (SPEC 04),
# V_box = 0.01 m^3; m1rate is the storage first-moment rate (m^4/s) and should equal N V_p |w| = MOM V_box once
# settled; 'ident' is the per-step ratio uf_eps_z V_box / m1rate, the identity form C enforces (1 to solver tolerance).
import re, sys, statistics as st
MOM = 4.24115e-4
VBOX = 0.01
def num(s): return float(s.replace('D','E'))
def parse(path):
    slip=[]; flux={}; mom={}
    for line in open(path):
        if line.startswith('EL_MEAN_SLIP t='):
            t=num(re.search(r't=\s*(\S+)',line).group(1)); m=re.search(r'uf_super=\s*(\S+)\s+(\S+)\s+(\S+)',line)
            slip.append((t,num(m.group(3))))
        elif line.startswith('EL_VA_FLUX t=') or line.startswith('EL_VA_FLUXMOM t='):
            d=dict(re.findall(r'(\w+)=\s*([-+0-9.EeDdNa]+)',line))
            key='flux' if line.startswith('EL_VA_FLUX t=') else 'mom'
            t=num(d['t']); corr=d.get('corr','0')
            if corr!='-1': continue
            (flux if key=='flux' else mom).setdefault(t,{}).update(d)
    return slip, flux, mom
def mean(xs): return (st.mean(xs), st.pstdev(xs)) if xs else (float('nan'),0)
for path in sys.argv[1:]:
    slip,flux,mom=parse(path)
    g5=[u/MOM for t,u in slip if t>=0.002]; g5l=[u/MOM for t,u in slip if t>=0.1]
    ts=sorted(flux)
    ue=[num(mom[t]['uf_eps_z'])/MOM for t in ts if t>=0.002 and t in mom]
    up=[num(mom[t]['uf_p0_z'])/MOM for t in ts if t>=0.002 and t in mom]
    res=[num(flux[t]['res1_L1']) for t in ts if t>=0.002]; resL=[num(flux[t]['resL_L1']) for t in ts if t>=0.002]
    sk=[abs(num(flux[t]['sink_eps'])) for t in ts if t>=0.002]; phi=[abs(num(flux[t]['Phi'])) for t in ts]; S=[abs(num(flux[t]['S'])) for t in ts]
    idz=[abs(num(mom[t]['id_z'])) for t in ts if t in mom]
    uel=[num(mom[t]['uf_eps_z'])/MOM for t in ts if t>=0.1 and t in mom]
    m1l=[num(mom[t]['m1rate'])/(MOM*VBOX) for t in ts if t>=0.1 and t in mom]
    ident=[num(mom[t]['uf_eps_z'])*VBOX/num(mom[t]['m1rate']) for t in ts if t>=0.002 and t in mom]
    name=path.rstrip('/').split('/')[-2] if '/' in path else path
    print(f"{name}: steps={len(slip)} G5(P0,EL_MEAN_SLIP)={mean(g5)[0]:.4f}+-{mean(g5)[1]:.4f} (t>=0.1: {mean(g5l)[0]:.4f}) | "
          f"G5(eps_h)={mean(ue)[0]:.5f}+-{mean(ue)[1]:.5f} (t>=0.1: {mean(uel)[0]:.5f}+-{mean(uel)[1]:.5f}) G5(P0,Q2rule)={mean(up)[0]:.5f} | "
          f"m1rate/(N V_p|w|) t>=0.1: {mean(m1l)[0]:.6f} | ident={mean(ident)[0]:.8f}+-{mean(ident)[1]:.1e} | "
          f"res1_L1 max={max(res) if res else 0:.2e} resL max={max(resL) if resL else 0:.2e} | |sink_eps| max={max(sk) if sk else 0:.2e} "
          f"|Phi| max={max(phi) if phi else 0:.1e} |S| max={max(S) if S else 0:.1e} |id_z| max={max(idz) if idz else 0:.1e}")
