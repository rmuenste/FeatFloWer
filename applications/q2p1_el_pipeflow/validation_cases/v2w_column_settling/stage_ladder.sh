#!/bin/sh
# stage_ladder.sh -- stage the EL column-settling twin of DNS D3.2 as
# independent run directories under <target-root>:
#   q2p1_el_rundir_v2w_ut            single sphere at (0,0,20): u_t reference
#   q2p1_el_rundir_v2w_n{20,40,80,120}_s{1,2,3}   the twelve DNS clouds, byte-identical
# Decks/json are the committed case files; only particles.xyz changes per run.
# Usage: stage_ladder.sh <target-root>
set -eu
[ $# -eq 1 ] || { echo "usage: $0 <target-root>" >&2; exit 2; }
ROOT=$(CDPATH= cd -- "$1" && pwd)
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$HERE/../../../.." && pwd)
STAGE=$REPO/tools/el_stage_rundir.sh
"$STAGE" -f v2w_column_settling "$ROOT/q2p1_el_rundir_v2w_ut" >/dev/null
echo "staged $ROOT/q2p1_el_rundir_v2w_ut (u_t reference)"
for n in n20 n40 n80 n120; do
  for s in s1 s2 s3; do
    d=$ROOT/q2p1_el_rundir_v2w_${n}_${s}
    "$STAGE" -f v2w_column_settling "$d" >/dev/null
    cp "$HERE/clouds/${n}_${s}.xyz" "$d/particles.xyz"
    echo "staged $d ($(wc -l < "$d/particles.xyz") particles)"
  done
done
