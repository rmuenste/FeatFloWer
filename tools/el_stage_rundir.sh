#!/bin/sh
# el_stage_rundir.sh -- stage an Euler-Lagrange validation case into an
# independent run directory so that concurrent runs never share (and
# clobber) a single staged _data/q2p1_param.dat, example.json, etc.
#
# Usage:
#   tools/el_stage_rundir.sh [-f] <case-name> <target-dir>
#
#   <case-name>   subdirectory of applications/q2p1_el_pipeflow/validation_cases/
#                 (e.g. pipe_hp_check, v3_ss_frozen, v1b_tencate_settling)
#   <target-dir>  directory to stage into (created if missing; an existing
#                 NON-empty directory is refused unless -f is given)
#
# Environment overrides:
#   EL_BIN      absolute path of the q2p1_el_pipeflow binary printed in the
#               final mpirun hint (never copied)
#   EL_REF_DIR  reference staged directory whose _adc symlink (if any) is
#               replicated into the target dir
#
# All case inputs are taken from the repository SOURCE tree, never from a
# shared staged build directory.

set -eu

usage() {
    echo "usage: $0 [-f] <case-name> <target-dir>" >&2
    exit 2
}

FORCE=0
if [ "${1:-}" = "-f" ]; then
    FORCE=1
    shift
fi
[ $# -eq 2 ] || usage
CASE=$1
TARGET=$2

# Resolve repo root from this script's location (tools/..)
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)

CASES_ROOT=$REPO/applications/q2p1_el_pipeflow/validation_cases
CASE_DIR=$CASES_ROOT/$CASE

EL_BIN=${EL_BIN:-/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/q2p1_el_pipeflow}
EL_REF_DIR=${EL_REF_DIR:-/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow}

if [ ! -d "$CASE_DIR" ]; then
    echo "error: unknown case '$CASE' (no directory $CASE_DIR)" >&2
    echo "available cases:" >&2
    ls "$CASES_ROOT" >&2
    exit 1
fi
for f in q2p1_param.dat example.json particles.xyz; do
    if [ ! -f "$CASE_DIR/$f" ]; then
        echo "error: case '$CASE' is missing $f" >&2
        exit 1
    fi
done
if [ ! -d "$CASE_DIR/mesh" ]; then
    echo "error: case '$CASE' has no mesh/ directory" >&2
    exit 1
fi
if [ ! -f "$REPO/_data/MG.dat" ]; then
    echo "error: repo-level _data/MG.dat not found under $REPO" >&2
    exit 1
fi

# Refuse to stage over an existing non-empty directory unless -f was given.
if [ -e "$TARGET" ]; then
    if [ ! -d "$TARGET" ]; then
        echo "error: target '$TARGET' exists and is not a directory" >&2
        exit 1
    fi
    if [ -n "$(ls -A "$TARGET" 2>/dev/null)" ] && [ "$FORCE" -ne 1 ]; then
        echo "error: target '$TARGET' is not empty (use -f to stage anyway)" >&2
        exit 1
    fi
fi

mkdir -p "$TARGET"
TARGET=$(CDPATH= cd -- "$TARGET" && pwd)

echo "Staging case '$CASE'"
echo "  from: $CASE_DIR"
echo "  into: $TARGET"

# --- configuration files ---------------------------------------------------
mkdir -p "$TARGET/_data"
cp "$CASE_DIR/q2p1_param.dat" "$TARGET/_data/q2p1_param.dat"
cp "$REPO/_data/MG.dat"       "$TARGET/_data/MG.dat"
cp "$CASE_DIR/example.json"   "$TARGET/example.json"
cp "$CASE_DIR/particles.xyz"  "$TARGET/particles.xyz"

# --- mesh: coarse project folders + partitions -------------------------------
# mesh/ may itself be a symlink into another case (e.g. v3_ss_frozen ->
# ../pipe_hp_check/mesh), hence cp -RL to dereference everything.
for entry in "$CASE_DIR"/mesh/*; do
    [ -e "$entry" ] || continue
    name=$(basename -- "$entry")
    [ "$name" = "partitions" ] && continue
    rm -rf "$TARGET/${name:?}"
    cp -RL "$entry" "$TARGET/$name"
done

mkdir -p "$TARGET/_mesh"
if [ -d "$CASE_DIR/mesh/partitions" ]; then
    for entry in "$CASE_DIR"/mesh/partitions/*; do
        [ -e "$entry" ] || continue
        name=$(basename -- "$entry")
        rm -rf "$TARGET/_mesh/${name:?}"
        cp -RL "$entry" "$TARGET/_mesh/$name"
    done
fi

# --- output directories the solver expects at startup ------------------------
for d in _vtk _gmv _ns _dump _sol _sol/1 solution testresults protocols; do
    mkdir -p "$TARGET/$d"
done

# --- _adc symlink (mesh repository), mirrored from the reference staged dir --
if [ -L "$EL_REF_DIR/_adc" ]; then
    adc_target=$(readlink "$EL_REF_DIR/_adc")
    rm -f "$TARGET/_adc"
    ln -s "$adc_target" "$TARGET/_adc"
    echo "  _adc -> $adc_target"
fi

echo ""
echo "Staged OK. Run with (binary is referenced in place, not copied):"
echo ""
echo "  BIN=$EL_BIN"
echo "  mpirun --oversubscribe --wdir $TARGET -np 28 \"\$BIN\""
