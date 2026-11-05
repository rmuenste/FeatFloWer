import re, sys
D="/anvme/workspace/k115ce12-featflower/rundirs/q2p1_dns_rundir_d62_v1b/"
segs=[("run_slurm_seg1_t0.0-41.3.log",0.0,40.01),("run_slurm_seg2_t40.0-81.2.log",40.01,80.01),("run_slurm.log",80.01,1e9)]
pat=re.compile(r"DNS_PART_AXIS time=\s*(\S+)")
n=0
for f,t0,t1 in segs:
    for line in open(D+f):
        if "DNS_PART_AXIS" in line:
            t=float(pat.search(line).group(1))
            if t0<=t<t1: sys.stdout.write(line); n+=1
sys.stderr.write("lines %d\n"%n)
