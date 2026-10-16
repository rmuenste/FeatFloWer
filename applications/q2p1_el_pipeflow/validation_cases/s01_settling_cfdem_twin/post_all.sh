#!/bin/sh
# Post-process every finished 01_settling run with the spec's settling_post.py.
S=/data/warehouse17/rmuenste/code/liggghts-dem/LIGGGHTS-PUBLIC/md_docs/case_specs/01_settling
R=$(cd "$(dirname "$0")" && pwd)
for d in "$R"/*_*way_*; do
  n=$(basename "$d"); [ -s "$d/run_slurm.log" ] || continue
  case=${n%%_*}; tag=${n##*_}
  case $tag in std) drag=difelice_std;; df) drag=difelice;; sn) drag=schiller_naumann;; esac
  echo "=== $n"; python3 -I "$S/settling_post.py" ff "$d/run_slurm.log" --case "$case" --drag "$drag" 2>&1
done
