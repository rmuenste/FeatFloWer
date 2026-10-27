#!/bin/sh
# Submit every staged 01_settling same-mesh twin run (5 ranks = master + 4 z-slabs) and the
# three unit-box follow-ups (4 ranks) that have no run_slurm.log yet.
# Usage: sh submit_all.sh [probe]   -> 'probe' submits only L1_oneway_lev1 and L1_oneway_lev3
set -eu
VA=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer-va
export EL_BIN=$VA/build-va/applications/q2p1_el_pipeflow/q2p1_el_pipeflow
R=$(cd "$(dirname "$0")" && pwd)
if [ "${1:-}" = probe ]; then DIRS="$R/L1_oneway_lev1 $R/L1_oneway_lev3"; else DIRS=$(ls -d "$R"/*_*way_lev*/ "$R"/fu_*/ | sed "s:/\$::"); fi
for d in $DIRS; do
  [ -e "$d/job.sbatch" ] && { echo "skip (already submitted) $(basename "$d")"; continue; }
  n=$(basename "$d")
  case $n in
    fu_c_*) sh "$VA/tools/el_slurm_submit_np.sh" "$d" "s01t_$n" 02:00:00 short 4 1 100G ;;
    fu_*)   sh "$VA/tools/el_slurm_submit_np.sh" "$d" "s01t_$n" 02:00:00 short 4 1 8G ;;
    *)      sh "$VA/tools/el_slurm_submit_np.sh" "$d" "s01t_$n" 01:00:00 short 5 1 8G ;;
  esac
done
