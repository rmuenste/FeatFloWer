#!/bin/bash
# Build the Layer-1 FBM query stress test against an existing pe static library.
#   bash tools/fbm_query_stress/fbm_query_stress.build.sh [<pe build tree>] [<output>]
set -e
cd "$(dirname "$0")/../.."   # repo root
TREE=${1:-build-twin-fbmaccel}
OUT=${2:-/tmp/rmuenste/fbm_query_stress}
# two tree layouts: an FF build tree (libs/pe/...) or a standalone pe tree (./config.h)
if [ -f "$TREE/libs/pe/config.h" ]; then CFG=$TREE/libs/pe/config.h; LIB=$TREE/libs/pe/lib/libpe.a; DEPS=$TREE/_deps
else CFG=$TREE/config.h; LIB=$TREE/lib/libpe.a; DEPS=$TREE/_deps; fi
CGALINC=""; CGALLIB=""
if grep -q "PE_USE_CGAL 1\|define PE_USE_CGAL" "$CFG" 2>/dev/null || [ -n "${WITH_CGAL:-}" ]; then
  # match pe's own CGAL flags exactly (flags.make of the pe_static target):
  # PE_USE_CGAL changes the TriangleMesh layout, so the test TU must agree with the library,
  # and CGAL_DISABLE_GMP is what keeps CGAL's CORE headers compiling with gcc 13.
  CGALINC="-isystem /sfw/cgal/6.1.2/gcc13.2.0/include -DPE_USE_CGAL -DCGAL_DISABLE_GMP=1"; CGALLIB=""
fi
/sfw/openmpi/gcc13.2.x/4.1.6/ucx-threaded-noverbs/bin/mpicxx \
  -std=c++17 -O2 -w -DNDEBUG -DBOOST_BIND_GLOBAL_PLACEHOLDERS \
  -DPE_SERIAL_MODE -DPE_USE_EIGEN -DMOBILE_INFINITE=1 \
  -Dpe_CONSTRAINT_SOLVER=pe::response::HardContactAndFluid \
  -include $CFG $CGALINC \
  -Ilibs/pe -I$DEPS/eigen-src -I$DEPS/json-src/include \
  tools/fbm_query_stress/fbm_query_stress.cpp $LIB \
  -lboost_thread -lboost_filesystem -lboost_program_options -lboost_random -lboost_system $CGALLIB \
  -o "$OUT"
echo "built $OUT against $TREE"
