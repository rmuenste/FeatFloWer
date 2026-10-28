#!/bin/sh
# Write <rundir>/versions.txt from the git state of the source tree the binary was built in.
# usage: el_write_versions.sh <rundir> <binary> [note]
# The source worktree is the parent of the build directory (.../<worktree>/build-*/applications/...).
# Called by el_slurm_submit_np.sh at every submission, so a run directory cloned from another
# never keeps the clone's provenance.
set -e
RUNDIR=${1:?rundir}; BIN=${2:?binary}; NOTE=${3:-}
SRC=$(echo "$BIN" | sed 's#/build[^/]*/.*##')
{
  printf 'FeatFloWer HEAD %s (%s)%s\n' "$(git -C "$SRC" rev-parse HEAD)" "$(git -C "$SRC" rev-parse --abbrev-ref HEAD)" "${NOTE:+ $NOTE}"
  printf 'libs/pe HEAD %s\n' "$(git -C "$SRC/libs/pe" rev-parse HEAD)"
  printf 'dirty source files %s\n' "$(git -C "$SRC" status --short | grep -vc '^??' || true)"
  printf 'binary %s (mtime %s; %s)\n' "$BIN" "$(stat -c %y "$BIN")" "$(basename "$(dirname "$(dirname "$(dirname "$BIN")")")")"
  printf 'written %s by el_write_versions.sh\n' "$(date '+%F %T')"
} > "$RUNDIR/versions.txt"
