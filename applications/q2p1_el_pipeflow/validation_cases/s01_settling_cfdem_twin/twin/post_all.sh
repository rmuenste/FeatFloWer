#!/bin/sh
# Post-process every finished twin run with the peer's twin_post.py (one-way first, then two-way
# with its one-way partner). Prints the summary block per run; twin_series.csv lands in each dir.
R=$(cd "$(dirname "$0")" && pwd)
P=/data/warehouse17/rmuenste/code/liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/01_settling/twin/twin_post.py
for c in G2 L1 L2 L3; do for l in lev1 lev2 lev3 lev2_wf10p7; do
  one=$R/${c}_oneway_$l; two=$R/${c}_twoway_$l
  for d in $one $two; do
    [ -e "$d/run_slurm.log" ] || continue
    tail -1 "$d"/slurm-*.out 2>/dev/null | grep -q "exit=0" || { echo "== $(basename $d): NOT FINISHED / exit!=0"; continue; }
    echo "== $(basename $d)"
    case $d in *_oneway_*) python3 -I $P ff "$d" --case $c --log run_slurm.log ;;
                *) python3 -I $P ff "$d" --case $c --oneway "$one" --log run_slurm.log ;; esac
  done
done; done
