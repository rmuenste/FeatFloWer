#!/bin/sh
# Clone a finished EL run directory as a fresh, unstarted one (inputs only; copies, no hardlinks).
# usage: pk_clone_run.sh <source rundir> <new rundir>
# Left out: logs, VTK frames, particle logs, provenance and job files of the source run, and the
# derived files of its evaluation. Empty output directories are recreated.
set -eu
SRC=$(CDPATH= cd -- "$1" && pwd); DST=$2
[ -e "$DST" ] && { echo "exists: $DST" >&2; exit 1; }
mkdir -p "$DST"
rsync -a --exclude='run_slurm.log' --exclude='slurm-*.out' --exclude='versions.txt' \
  --exclude='job.sbatch' --exclude='el_particles_rank*.log' --exclude='twin_series.csv' \
  --exclude='/_vtk/*' --exclude='/_gmv/*' --exclude='*.png' --exclude='*.csv' --exclude='eval_*' \
  "$SRC"/ "$DST"/
mkdir -p "$DST/_vtk" "$DST/_gmv"
echo "cloned $SRC -> $DST ($(du -sh "$DST" | cut -f1))"
