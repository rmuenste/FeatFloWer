#!/bin/bash
# cleanup_scratch.sh -- remove the GPU pressure-Poisson scratch data on tardis's node-local /scratch:
#   the instance .npz files, the raw .bin dumps (bin/), the slab tmp dir, the hypre source/build tree and tarball.
# KEEPS the installed CUDA hypre (/scratch/rmuenste/gpu-poisson/hypre-cuda, small; the driver binary is kept in its bin/).
# Lists everything with sizes BEFORE deleting. Must run inside a Slurm job on tardis (the data is node-local), after all
# results have been copied to gpu-poisson-a100/results/. Usage: bash cleanup_scratch.sh [--dry-run]
set -u
SCR=/scratch/rmuenste/gpu-poisson
DRY=0; [ "${1:-}" = "--dry-run" ] && DRY=1
echo "cleanup of $SCR on $(hostname), $(date -Is), dry-run=$DRY, job ${SLURM_JOB_ID:-none}"
[ -d "$SCR" ] || { echo "$SCR does not exist here (run this on tardis)"; exit 1; }
echo "## /scratch before"; df -h /scratch | tail -1
echo "## contents before"; du -sh "$SCR"/* 2>/dev/null; du -sh "$SCR"
TARGETS=()
for t in "$SCR"/*.npz "$SCR"/bin "$SCR"/tmp "$SCR"/hypre-src "$SCR"/hypre-v3.1.0-src.tar.gz "$SCR"/hypre_boomer_pcg; do [ -e "$t" ] && TARGETS+=("$t"); done
echo "## to delete (${#TARGETS[@]} entries)"
for t in "${TARGETS[@]}"; do du -sh "$t"; done
echo "## kept"; du -sh "$SCR"/hypre-cuda 2>/dev/null || echo "(no hypre-cuda install present)"
if [ $DRY = 0 ]; then
  for t in "${TARGETS[@]}"; do rm -rf "$t" && echo "deleted $t"; done
else
  echo "(dry run: nothing deleted)"
fi
echo "## contents after"; ls -la "$SCR"; du -sh "$SCR"/* 2>/dev/null; du -sh "$SCR"
echo "## /scratch after"; df -h /scratch | tail -1
echo "cleanup done $(date -Is)"
