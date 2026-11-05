# Tier-2 momentum_conservation: re-baseline for branch feature/el-volume-averaged

Binary: worktree FeatFloWer-va, build-va, FF 5d8c2734 (feature/dns-validation
fa6c242d + merged feature/el-rz-postfix) with pe feature/el-box-walls-on-master
c02cf28 (= pe master 1f60291 + boxWallsEnabled_). Job 149704, 14 s, 702 PERIODIC_COMM.

Against the certified reference (job 137542, pe 1fef6b8) the EL streams are NO
LONGER bitwise: pe master (PRs #35-#38) now parses particles.xyz as double, so
the seed position is the literal 0.4137/0.5273/0.6191 instead of the former
float-rounded 0.41370001435/0.52730000019/0.61909997463 (1.4e-8 shift). All
differences are at that level (see numbers below); EL_CONTACT_STATS identical.
The same shift is the "NEW DKT reference head" of the DNS campaign.

reference/run_g0_twin_va_2026-10-14.log is therefore the bitwise reference for
every subsequent change on this branch (volume-averaged keys OFF must reproduce
it exactly).
    EL_MOMENTUM_ELEMINT      lines 60/60  worst relative difference 6.87e-06  (-1.019193E-12->-1.019200E-12)
    EL_MEAN_SLIP             lines 60/60  worst relative difference 6.78e-06  (1.031871E-08->1.031878E-08)
    EL_NEWTON_PAIR           lines 60/60  worst relative difference 0.00e+00  ()
    EL_VOLUME_CONSERVATION   lines 60/60  worst relative difference 0.00e+00  ()
