#!/bin/sh
# Submit every staged 01_settling run (4 ranks, 1 node, short partition) that has no run_slurm.log yet.
# Usage: sh submit_all.sh [probe]   -> 'probe' submits only M-L2_oneway_std
set -eu
VA=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va
export EL_BIN=$VA/build-va/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
R=$(cd "$(dirname "$0")" && pwd)
if [ "${1:-}" = probe ]; then DIRS="$R/M-L2_oneway_std"; else DIRS=$(ls -d "$R"/*_*way_*); fi
for d in $DIRS; do
  [ -e "$d/run_slurm.log" ] && { echo "skip (has log) $(basename "$d")"; continue; }
  sh "$VA/tools/el_slurm_submit_np.sh" "$d" "s01_$(basename "$d")" 00:30:00 short 4 1 8G
done
