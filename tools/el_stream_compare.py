"""Bitwise comparison of the EL_* diagnostic streams of two q2p1_el_pipeflow logs.
usage: el_stream_compare.py <reference.log> <candidate.log>
Per tag: line counts, worst relative difference over the numeric tokens the reference line has
(extra trailing tokens in the candidate are reported, not compared)."""
import re, sys
num = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
TAGS = ("EL_MOMENTUM_ELEMINT", "EL_MEAN_SLIP t=", "EL_MEAN_SLIP_TAVG", "EL_NEWTON_PAIR", "EL_VOLUME_CONSERVATION",
        "EL_FORCE_BUDGET", "EL_FEEDBACK_CONSERVATION", "EL_TERMINAL_VEL", "EL_CONTACT_STATS", "EL_VA_SOURCE",
        "EL_SUSP_STRESS", "EL_CONT_STRESS", "EL_HALO_RECORDS", "EL_LUB_IMPULSE", "EL_VA_MOMENT")
def load(p):
    d = {t: [] for t in TAGS}
    for line in open(p, errors="replace"):
        for t in TAGS:
            if line.startswith(t):
                d[t].append([float(x.replace("D", "E").replace("d", "e")) for x in num.findall(line[len(t):])])
                break
    return d
a, b = load(sys.argv[1]), load(sys.argv[2])
allbit = True
for t in TAGS:
    la, lb = a[t], b[t]
    if not la and not lb: continue
    worst, ex, extra = 0.0, "", 0
    for ra, rb in zip(la, lb):
        if len(rb) > len(ra): extra = max(extra, len(rb)-len(ra))
        for x, y in zip(ra, rb):
            if x == y: continue
            r = abs(x-y)/max(abs(x), abs(y), 1e-300)
            if r > worst: worst, ex = r, f"{x:.6E}->{y:.6E}"
    status = "bitwise" if worst == 0.0 and len(la) == len(lb) else "DIFF"
    if status == "DIFF": allbit = False
    print(f"  {t:24s} lines {len(la)}/{len(lb)}  worst rel diff {worst:.2e} {('('+ex+')') if ex else ''} {('+'+str(extra)+' extra tokens') if extra else ''}  {status}")
print("ALL COMMON STREAMS BITWISE" if allbit else "NOT BITWISE")
