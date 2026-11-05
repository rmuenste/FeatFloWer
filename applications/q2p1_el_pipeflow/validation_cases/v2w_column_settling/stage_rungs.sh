#!/bin/sh
# stage_rungs.sh -- discriminator rungs for the column twin (RUNBOOK section 7):
#   k15  : kernel width delta = 1.5 d (ELKernelWidthFactor 1.5), level 2; u_t + all 12 clouds
#   l3   : MaxMeshLevel 3 (h = d/8, the DNS cell size), delta = 2.5 d; u_t + n40 x3 + n120 x3
# Usage: stage_rungs.sh <target-root>
set -eu
[ $# -eq 1 ] || { echo "usage: $0 <target-root>" >&2; exit 2; }
ROOT=$(CDPATH= cd -- "$1" && pwd)
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$HERE/../../../.." && pwd)
STAGE=$REPO/tools/el_stage_rundir.sh
stage() { # tag cloud-or-ut sedexpr
  d=$ROOT/q2p1_el_rundir_v2w${1}_$2
  "$STAGE" -f v2w_column_settling "$d" >/dev/null
  [ "$2" = ut ] || cp "$HERE/clouds/$2.xyz" "$d/particles.xyz"
  sed -i -e "$3" "$d/_data/q2p1_param.dat"
  echo "staged $d"
}
K='s/^SimPar@ELKernelWidthFactor = .*/SimPar@ELKernelWidthFactor = 1.5d0/'
L='s/^SimPar@MaxMeshLevel = .*/SimPar@MaxMeshLevel = 3/'
stage k15 ut "$K"
for n in n20 n40 n80 n120; do for s in s1 s2 s3; do stage k15 ${n}_${s} "$K"; done; done
stage l3 ut "$L"
for n in n40 n120; do for s in s1 s2 s3; do stage l3 ${n}_${s} "$L"; done; done
