#!/bin/sh
# stage_postfix_sweep.sh -- stage the v2_rz_settling POST-FIX re-run
# (periodic Q2 coupling fixed on 2026-08-03, see RUNBOOK section
# "Post-fix re-run") as independent run directories.
#
# Produces, under <target-root>:
#   q2p1_el_rundir_v2pf_smoke   phi=0.20, 200 steps  (gate before the sweep)
#   q2p1_el_rundir_v2pf_phi05   phi=0.05, 50000 steps to t=100
#   q2p1_el_rundir_v2pf_phi10   phi=0.10
#   q2p1_el_rundir_v2pf_phi15   phi=0.15
#   q2p1_el_rundir_v2pf_phi20   phi=0.20
#   q2p1_el_rundir_v2pf_u0      single particle U0 reference, 5000 steps to t=10
#
# The decks are DERIVED from the committed case deck exactly as sweep 3
# (jobs 135866-135870) derived them: production horizon, per-phi counter
# force phi*(rho_p-rho_f)*g, ELMomentumFix=Yes (legacy convection),
# OutputFreq 10, PE vtk every 250 steps. Nothing else changes, so the only
# difference to sweep 3 is the binary (periodic coupling present).
#
# Usage:  stage_postfix_sweep.sh <target-root>      (root must exist)

set -eu
[ $# -eq 1 ] || { echo "usage: $0 <target-root>" >&2; exit 2; }
ROOT=$(CDPATH= cd -- "$1" && pwd)
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$HERE/../../../.." && pwd)
STAGE=$REPO/tools/el_stage_rundir.sh

# phi -> counter force z = phi * (rho_p - rho_f) * g = phi * 1 * 0.981
force_for() {
    case $1 in
        0.05) echo 0.04905d0 ;;
        0.10) echo 0.0981d0 ;;
        0.15) echo 0.14715d0 ;;
        0.20) echo 0.1962d0 ;;
        *) echo "bad phi $1" >&2; exit 1 ;;
    esac
}

stage_one() {   # name phi steps tmax vf_json
    name=$1; phi=$2; steps=$3; tmax=$4; vf=$5
    dir=$ROOT/q2p1_el_rundir_$name
    "$STAGE" -f v2_rz_settling "$dir" >/dev/null
    deck=$dir/_data/q2p1_param.dat
    json=$dir/example.json
    sed -i \
        -e "s/^SimPar@MaxSimTime = .*/SimPar@MaxSimTime = ${tmax}d0/" \
        -e "s/^SimPar@MaxNumStep = .*/SimPar@MaxNumStep = ${steps}/" \
        -e "s/^SimPar@ConstantForcing = .*/SimPar@ConstantForcing = 0d0,0d0,$(force_for "$phi")/" \
        -e "s/^SimPar@OutputFreq = .*/SimPar@OutputFreq = 10.0d0/" \
        "$deck"
    printf 'SimPar@ELMomentumFix = Yes\n' >> "$deck"
    sed -i \
        -e "s/\"timesteps_\": [0-9]*/\"timesteps_\": ${steps}/" \
        -e "s/\"vtk_\": false/\"vtk_\": true/" \
        -e "s/\"visspacing_\": [0-9]*/\"visspacing_\": 250/" \
        -e "s/\"volumeFraction_\": [0-9.]*/\"volumeFraction_\": ${vf}/" \
        "$json"
    echo "staged $dir  (phi=$phi steps=$steps)"
}

stage_one v2pf_smoke 0.20   200   0.4001 0.2
stage_one v2pf_phi05 0.05 50000 100.001 0.05
stage_one v2pf_phi10 0.10 50000 100.001 0.1
stage_one v2pf_phi15 0.15 50000 100.001 0.15
stage_one v2pf_phi20 0.20 50000 100.001 0.2

# U0 reference: one particle from particles.xyz, no counter force, t=10.
u0=$ROOT/q2p1_el_rundir_v2pf_u0
"$STAGE" -f v2_rz_settling "$u0" >/dev/null
sed -i \
    -e "s/^SimPar@MaxSimTime = .*/SimPar@MaxSimTime = 10.001d0/" \
    -e "s/^SimPar@MaxNumStep = .*/SimPar@MaxNumStep = 5000/" \
    -e "s/^SimPar@UseConstantForcing = .*/SimPar@UseConstantForcing = No/" \
    -e "s/^SimPar@OutputFreq = .*/SimPar@OutputFreq = 10.0d0/" \
    "$u0/_data/q2p1_param.dat"
printf 'SimPar@ELMomentumFix = Yes\n' >> "$u0/_data/q2p1_param.dat"
sed -i \
    -e "s/\"timesteps_\": [0-9]*/\"timesteps_\": 5000/" \
    -e "s/\"vtk_\": false/\"vtk_\": true/" \
    -e "s/\"visspacing_\": [0-9]*/\"visspacing_\": 250/" \
    -e "s/\"seedMode_\": \"random\",/\"seedMode_\": \"file\",/" \
    -e "/\"volumeFraction_\":/d" \
    "$u0/example.json"
echo "staged $u0  (U0 reference)"
