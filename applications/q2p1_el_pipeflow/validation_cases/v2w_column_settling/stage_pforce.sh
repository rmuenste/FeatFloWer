#!/bin/sh
# stage_pforce.sh -- Model-A pressure-force rung (deck-only, no code change):
#   ELFluidGravity = No  (fluid gravity absorbed into the pressure; in a closed
#                         box a uniform body force is a pure gradient, so the
#                         dynamics are unchanged) and
#   ELPressureForce = Yes (particle gets -V grad p_dyn, which inside a settled
#                         cloud carries the suspension-weight gradient phi*drho*g;
#                         analytic buoyancy (rho_p-rho_f) V g stays).
# u_t + N=40 x3 + N=120 x3 at level 2, delta = 2.5 d.  Usage: stage_pforce.sh <target-root>
set -eu
[ $# -eq 1 ] || { echo "usage: $0 <target-root>" >&2; exit 2; }
ROOT=$(CDPATH= cd -- "$1" && pwd)
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$HERE/../../../.." && pwd)
STAGE=$REPO/tools/el_stage_rundir.sh
P='s/^SimPar@ELFluidGravity = .*/SimPar@ELFluidGravity = No/; s/^SimPar@ELPressureForce = .*/SimPar@ELPressureForce = Yes/'
stage() { d=$ROOT/q2p1_el_rundir_v2wpf_$1; "$STAGE" -f v2w_column_settling "$d" >/dev/null; [ "$1" = ut ] || cp "$HERE/clouds/$1.xyz" "$d/particles.xyz"; sed -i -e "$P" "$d/_data/q2p1_param.dat"; echo "staged $d"; }
stage ut
for n in n40 n120; do for s in s1 s2 s3; do stage ${n}_${s}; done; done
